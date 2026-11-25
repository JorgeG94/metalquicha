# The SCF: guess, RIJCOSX, stability, second-order and direct minimisation

How a reference is started and finished, and what each route refuses. Read from `AGENTS.md`, which says when.

## RIJCOSX

RIJCOSX is not a method name either. `keywords.scf.exchange: "cosx"` beside
`keywords.scf.density_fitting` builds the Coulomb matrix from the fit and the
exchange matrix by chain-of-spheres quadrature (`mqc_czt_cosx`), for a
restricted Hartree-Fock or global-hybrid energy on the CPU backend; everything
else is refused by name. Three things to know before touching it:

- **The reference for it is a fitted Coulomb matrix with exact exchange**
  (PySCF's `only_dfj`), not Hartree-Fock and not another code's seminumerical
  exchange. The fit is then common to both numbers and only the quadrature
  separates them. Its validation tolerances are 5e-5 to 1.5e-4 hartree because
  they bound an approximation, not an SCF.
- **It needs a `-jfit` or `-jkfit` auxiliary basis.** With an `-rifit` set it
  is millihartrees off, because RI-JK's Coulomb and exchange errors cancel and
  RIJCOSX has no exchange fit to cancel against.
- **Two thirds of a run is one integral**, `<mu| 1/|r - r_g| |nu>` over
  batches of grid points, and it does not come from libfint: `mqc_czt_esp_batch`
  is a McMurchie-Davidson kernel written for it, about three and a half times
  faster than `int1e_grids` and checked against it to 1e-11 for every shell
  type through g. Profile before touching it (`perf` with `-e cpu-clock`
  works here once `linux-tools-generic` is installed): the things that
  mattered were a logical reduction and a `minval` that do not vectorise, a
  double-precision `1/sqrt`, and adding Hermite terms one to a pass. Its
  Coulomb matrix comes from `df_jfit`, the same Coulomb-only fit a pure
  functional uses. Against RI-JK it is slower on small systems, level by
  twenty waters in def2-SVP, and runs where the fitted tensor does not fit.
  See `mqc_docs/source/rijcosx.rst`.

## Initial guess

Initial guess is `keywords.scf.guess`: `core`, `gwh`, `sac`, `sad`, or `auto`.
`auto` resolves per backend - `sad` on the CPU path, `gwh` on cuEST - so the two
can differ without either knowing what the other chose.

## Wavefunction stability, and the second-order SCF

`keywords.scf.stability` asks whether a converged SCF is a *minimum* rather than
merely a stationary point, which is not the same question and is not answered by
any convergence test: a saddle point converges, reports an ordinary energy, and
says nothing. `keywords.scf.second_order` is the other half: it converges the
closed-shell reference by trust-region Newton once DIIS has got it close, so it
can *escape* a saddle rather than only detect one.

**Both are core functionality and work in a default build**, with nothing
optional fetched. If you are adding to the stability analysis, add to
`mqc_czt_native_stability`.

Two modules, in `backends/cenzontle/stability/`:

- `mqc_czt_ov_hessian` is `(A+B)` -- the real singlet orbital-rotation Hessian --
  as a matrix-free operator over the non-redundant occupied-virtual rotations.
  It needs nothing optional and is always compiled. It is **recovered from the
  coupled-perturbed operator rather than built**: `nuclear_apply` computes
  `Delta^-1 K`, so the Hessian is `Delta * (x - apply(x))`, and that module's
  header derives it. The redundant occupied-occupied rows the response operator
  carries are excluded, because an eigensolver handed them finds `n_occ^2` exact
  zeros and calls every reference marginally stable.
- `mqc_czt_native_stability` diagonalises it with `mqc_davidson` -- the CI
  solver, unmodified except for two optional labels so a verbose run does not
  claim to be a CI. `davidson_flat` already works against an abstract
  `sigma_operator_t` over a flat vector, its preconditioner is exactly the
  orbital-energy diagonal this wants, and its no-guess start is the unit vector
  on the smallest diagonal element, which over this space is the HOMO-LUMO
  rotation. The adapter is one type with one `apply`. Because
  `sigma_operator_t%apply` carries an `error_t`, this path needs **no module
  state** and no seed, and it **always reports the eigenvalue**. It is checked
  against a dense diagonalisation of the same operator, built column by column
  through it, which is a reference that depends on no eigensolver being right.
  (A second implementation of this one solve, borrowed from a library rather
  than written here, exists on a separate branch purely as a cross-check; it is
  not in this one and nothing here needs it.)

The second-order SCF is two more modules plus a phase:

- `src/methods/scf/mqc_orbital_rotation` holds `rotation_matrix` (the exact
  matrix exponential by scaling and squaring) and `level_shifted_step` (the
  level-shifted trust-region Newton step with the saddle escape), plus
  `MIN_CURVATURE`, `SADDLE_CURVATURE`, `MAX_ROTATION`, `MIN_ROTATION`,
  `ENERGY_RESOLUTION` and `TRUST_GROWTH`. **Lifted out of `mqc_czt_mcscf`, not
  copied**: CASSCF's `newton_step` is now just the mapping between its
  `(rows, cols)` pair list and a flat vector, and its numbers are unchanged.
  `rotation_matrix` is still public from `mqc_czt_mcscf`, use-associated, so
  existing callers and tests are untouched.
- `backends/cenzontle/scf/mqc_czt_soscf` holds the three things specific to an
  SCF: the gradient out of the Fock matrix, the semicanonicalisation, and the
  Krylov solve of the Newton equations against the matrix-free Hessian.
- `run_soscf_phase`, private in `mqc_czt_rhf`, is the loop. It lives there
  because it needs `assemble_fock`, `commutator` and the `rhf_operators_t` /
  `rhf_state_t` views, and because convergence must be the *same*
  `scf_convergence_t` rule the DIIS phase used or the two energies are not
  comparable.

Things to know before touching any of it:

- **The factor is four, and it cancels out of the step.** The gradient is
  `4 F_ai` and the Hessian is `4 (A+B)`, so a Newton step `-H^-1 g` is right
  even with the factor dropped -- and the *curvature* is then a quarter of the
  truth, which silently retunes `MIN_CURVATURE` and `SADDLE_CURVATURE`.
  `test_mqc_czt_soscf.f90` differences the energy numerically to fix both the
  factor and the sign; do not change either on the strength of a comment.
- **`(A+B)` is the Hessian only at semicanonical orbitals.** Its first term is
  the orbital-energy gaps, which is `delta_ij F_ab - delta_ab F_ij` only when
  the Fock matrix is diagonal within each block. A rotation destroys that, so
  every Newton iteration re-diagonalises the occupied and virtual blocks first.
  That moves no density and costs no Fock build.
- **The curvature the SCF reports is an upper bound.** It is a Rayleigh-Ritz
  value from a Krylov subspace, so a negative one proves a saddle and a positive
  one proves nothing. The convergence test uses it one-sidedly, to *refuse*
  convergence at a saddle. `keywords.scf.stability` is what proves a minimum.
- **DIIS first, always.** A Newton step from an initial guess diverges. The
  handover is at `keywords.scf.soscf_start` on the commutator, defaulting to
  1e-2, and the iteration it happened at is reported.
- **Refused, not approximated**: a continuum solvent (no surface-charge
  response in the Hessian), and a Fock projector (forbidden rotations are not
  excluded from the parameter space). Restricted Kohn-Sham *is* supported,
  through the same kernel the analytic Hessian's coupled-perturbed solve uses.
- **Unrestricted is `run_soscf_phase_uhf`**, over one vector of alpha's
  rotations and then beta's. Its Hessian is `uhf_rotation_hessian_t` in
  `mqc_czt_soscf`: twice `response_product_uhf`'s spin-blocked `(A+B)`, with
  the gradient `2 F_ai` per spin, both fixed by finite differences in
  `test_mqc_czt_soscf`. Direct builds only: a fitted or stored-integral
  unrestricted SCF is refused. It searches the spin-symmetry-breaking
  rotations, so a singlet run unrestricted can leave the restricted solution:
  the stretched N2 below goes to -108.7504, 0.179 below the restricted
  minimum. `keywords.scf.stability` is still closed-shell only.

## Geometric direct minimisation

`keywords.scf.gdm` (or `accelerator: "gdm"`) is the third way to finish an
SCF, restricted or unrestricted: DIIS opens, and below `gdm_start` (1e-1) an L-BFGS
minimisation over the orbital rotations takes over, one Fock build per step, a
step kept only if the energy did not rise. Van Voorhis and Head-Gordon, 2002.

- `src/methods/scf/mqc_scf_gdm` is the numerics -- history, preconditioner,
  accept rule, step cap -- and is **terco's `trc_scf_gdm`, name for name**. A
  fix to one belongs in the other; their tests are the same fourteen.
- `run_gdm_phase`, private in `mqc_czt_rhf` beside `run_soscf_phase` and
  shaped like it, is the loop. Every Fock build in it is a full one, because
  the accept rule compares energies to 1e-10.
- **The orbitals are not semicanonicalised between steps**, unlike the
  second-order phase. The L-BFGS history compares gradients across successive
  orbital frames with no transport, which holds only while the frames differ
  by one small rotation; mixing within the occupied or virtual block would
  make the history meaningless. Once, at the end.
- **It is a minimiser, not a saddle escape.** On the stretched N2 below it
  lands on DIIS's saddle, and `gdm/n2_stretched_gdm` pins that on purpose.
- `run_gdm_phase_uhf` is the same loop over two spins: one vector, alpha's
  rotations then beta's, one history, gradient `2 F_ai` per spin with the
  preconditioner's factor two to match. **A change to one phase's accept or
  restart logic belongs in the other.**
- Refused by name: with `second_order`, COSX (the energy is not
  variational, so the Fock matrix is not its gradient), PCM, a Fock projector
  and cuEST. On the terco backend the request is forwarded to terco's own GDM
  (`trc_set_gdm`).

**A real case where DIIS is not enough**, worth keeping in mind as the thing
these two features are for: N2 at 1.6 A in 6-31G, RHF. DIIS converges from the
core, GWH and SAD guesses alike, to a stationary point the stability analysis
reports as a *saddle* -- curvature -1.29e-1 from core, -6.77e-2 from the other
two. With `second_order: true` the same run reaches -108.5718 hartree, 0.020 to
0.242 hartree lower, with a curvature of -1e-12: a minimum. It costs 225 Fock
builds against DIIS's 12, and DIIS's answer is wrong.
