# Methods: what each one is, and what it refuses

Semi-empirical through coupled cluster, the terco backend, effective core
potentials, analytic derivatives, and continuum solvation. Read from `AGENTS.md`, which says when.

## Semi-empirical

Semi-empirical, via tblite:

- `GFN1` - GFN1-xTB (faster, older parametrization)
- `GFN2` - GFN2-xTB (recommended, more accurate)

## Cenzontle, the CPU ab initio backend

Ab initio on the CPU is **cenzontle** (`backends/cenzontle/`, modules `mqc_czt_*`,
built by `-DMQC_ENABLE_CZT=ON`) -- named for the bird of four hundred voices,
because one backend speaks HF, DFT, MP2, CC, MCSCF, EFP and SAPT.

It is **not** a wrapper around an integral library, which is why it stopped
being called `libcint`. The integrals come from **libfint**, an all-Fortran port
of libcint, so there is no C compiler and no second integral backend. libfint
exports the `libcint_fortran` module and the `libcint_*`/`LIBCINT_*` symbols
under libcint's own spelling, and sources under `backends/cenzontle/` `use` it
directly. Those spellings are libfint's export names and are left alone; only our
own names are `czt`. libfint's own test suite checks it against libcint value for
value, so this repository does not build against libcint. ECP integrals exist
only in libfint.

## Effective core potentials

**An ECP has an analytic gradient and no analytic Hessian.** Its first
derivative enters every gradient through `mqc_czt_gradient` and nowhere
else: `one_electron_deriv` with `DERIV_NUC` adds `<nabla mu|U|nu>` over every
centre, and `nuclear_attraction_deriv_at` adds one atom's potential to that
atom's `-Z/r` (`hellmann_feynman_gradient` does the same for the SCF, without
storing the matrix). A gradient that builds its
core-Hamiltonian derivative from those two carries the potential with nothing
further; one that calls `iprinv_deriv_at` directly does not. The Hessian's
first-order pieces go through them too (`add_hcore_deriv` in
`mqc_czt_hessian`), and `nuclear_repulsion_hessian` takes the molecule's
charges, not atomic numbers. Second derivatives are not in libfint, so every analytic
Hessian is skipped when `mol%necpbas > 0` and the semi-numerical path
differences the gradient instead; `ecp_refuses_derivatives` is the backstop
inside each Hessian routine. MCSCF, SAPT, EFP, MakeFP, FMO and EFMO refuse a
potential in the schema (`check_ecp_supported`): only `run_czt_hf` passes one
to the molecule builder, and a new entry point that builds its own molecules
belongs on that list until it does too.
**Test with the heavy atom late in the atom list**: `LIBCINT_ATOM_OF` indexes a
`bas` row directly while the `env` slots want `+ 1`, and confusing the two made
`ecp_on_atom` right for the first three atoms only, which every molecule with
its heavy atom first passes.

## Method by method

The methods:

- `HF` - restricted and unrestricted; direct, in-core, or density-fitted
  The unrestricted fitted build is `build_fock_df_uhf`: the restricted
  `build_fock_df` called once for J on the total density and once per spin
  for K, with each spin density's eigenvectors standing in for orbitals, so
  it is right for a guess or a mixed density too. Always the dense tensor.
  Coupled cluster and a double hybrid over a fitted unrestricted reference
  are refused by name; a range-separated open-shell gradient is refused
  fitted or not.
  **ROHF and ROKS are `model.restricted_open`, and run inside `run_czt_uhf`**
  (`restricted_open=`): the two spin Fock matrices are built as for UHF, then
  replaced by Roothaan's effective matrix (`roothaan_fock`), diagonalised
  once, with the DIIS error taken against the total density. After
  convergence each spin's orbitals are rotated to diagonalise that spin's own
  Fock matrix within its occupied and virtual blocks. **That rotation is why
  nothing downstream knows about ROHF**: the result is shaped as an
  unrestricted one, and `C_occ eps C_occ^T` is then `D F D` per spin, whose
  sum is the ROHF energy-weighted density, so the UHF gradient is the ROHF
  gradient unchanged. It also means **alpha's first `n_beta` orbitals are not
  the doubly occupied ones** -- alpha's occupied block mixes closed and open
  -- so project with `D_b` and `D_a - D_b`, never by orbital index. Anything
  that assumes a diagonal Fock matrix across occupied and virtual (MP2, CC,
  double hybrids, TDDFT) is refused in the bridge; GDM and the second-order
  SCF are refused in the driver. OH under a pure functional does not
  converge, here or in PySCF: its open orbital is one of a degenerate pair.
- `DFT` - restricted and unrestricted Kohn-Sham over libxc, through
  `model.functional`. LDA through double hybrid, including range-separated
  hybrids, whose long-range exchange needs the direct build
- `MP2`, `SCS-MP2`, `SOS-MP2`, `RI-MP2` - restricted, or unrestricted as an
  energy only (`run_czt_ump2`, `run_czt_uri_mp2`, exact or fitted reference;
  frozen orbitals are counted per spin). Gradients and Hessians are RHF only,
  with one exception: **RI-MP2 over an exact UHF reference has an analytic
  gradient**, `czt_uri_mp2_gradient` in `mqc_czt_uri_mp2_gradient`, and its
  Hessian is that gradient differenced. It is written from the Hylleraas
  functional per spin, not translated from the closed-shell routine, and its
  header states the functional. Three things to know: the coupled-perturbed
  operator is `response_product_uhf`, exact and integral-direct, so a fitted
  reference is refused rather than given the wrong operator; the Z-vector
  solve is a Galerkin Krylov one and not conjugate gradients, because **an
  unrestricted reference is often a saddle** (triplet water in 6-31G is) and
  its gradient is still defined; and the reference's two-electron term and
  the cross term come from `two_electron_gradient` called at
  reference-plus-relaxed and at relaxed alone, the energy being quadratic.
  PySCF has no analytic gradient for its fitted unrestricted MP2, so the
  validation references are its energy differenced and extrapolated
  (`pyscf_uri_mp2_gradient` in the generator), good to 1e-9. **Central
  differences at 1e-3 bohr are 3e-7 off on a stiff bond**: extrapolate
  before concluding a gradient is wrong.
  Analytic Hessians: MP2 and RI-MP2 over an exact RHF reference (frozen core
  allowed), and a double hybrid's perturbative term, fitted when
  `model.aux_basis` is named (`mqc_czt_ri_mp2_hessian`, no four-index array).
  RI-MP2 also over a density-fitted reference (`ref_aux` to
  `ri_mp2_correlation_hessian`, added to `df_rhf_hessian`): the block reaches
  the reference's operator in five places -- `t1_batch`, `t1x_stack`,
  `skeleton_matrices`, `solve_response` and `pass1_reference_share` -- and
  each has a fitted branch through `df_operator_t`. **A sixth use of the
  reference's integrals needs a fitted branch too**, or the Hessian is that of
  no energy. A double hybrid over a fitted Kohn-Sham reference takes the
  same branches with `xc`. SCS/SOS, conventional MP2 over a fitted reference,
  and a frozen-core double hybrid take the semi-numerical path. See
  `mqc_docs/source/analytic_hessians.rst`.
  **A double hybrid's Hessian spends its time on the functional, not on
  MP2**, and three of its kernel routines ran one perturbation at a time or
  on one thread: `ks_dl_pair_kernel` made two passes over the grid per
  Cartesian coordinate per density, `xc_kernel_deriv` and
  `xc_potential_hessian` were serial. On 1,3,5-trifluorobenzene, B2PLYP in
  6-31G, ladder A was 195 s of a 343 s run on four threads, and 75 s on
  fifty-two. `ks_dl_pair_kernel_many` takes every perturbation in two passes
  (`xc_kernel_apply_many`, `xc_kernel2_apply_many`), and the other two are
  threaded over grid blocks as their siblings `xc_potential_deriv` and
  `xc_hessian` already were: 47 s and 164 s, the Hessian the same to 1e-11.
  The conventional (not RI) double-hybrid Hessian still calls
  `ks_dl_pair_kernel` in a loop. **The nuclear responses are solved once**: the
  reference's Hessian returns its first-order orbitals (`mo1=`) and
  `ri_mp2_correlation_hessian` takes them (`mo1_in=`) in place of its own
  solve, for RI-MP2 and for a fitted double hybrid alike. That is right only
  because both solves had the same equations, operator and tolerance; the
  conventional MP2 block solves at 1e-13 where the reference solves at the
  deck's tolerance, so it does not share; a range-separated double hybrid
  never reaches the analytic path. `test_mqc_ri_mp2_hessian_dh` holds the
  shared block to the block's own solve, and `test_mqc_xc_hessian` holds the
  threaded kernel routines to one thread -- **ctest and the validation runner
  both run one thread**, so nothing else exercises a parallel region. **Read the phase lines before assuming
  which method a slow Hessian is**: plain RI-MP2 on the same molecule has a
  ladder A of 3 s.
  **Density-fitted Hartree-Fock and Kohn-Sham have their own analytic
  Hessian** (`df_rhf_hessian`, which takes an optional `xc`, with the fitted
  terms in `mqc_czt_df_hessian`): it differentiates the fitted energy as
  written, shares `rhf_hessian`'s one-electron terms and assembly, and reaches
  the response solver through one optional argument, `fit`. Under a functional
  the fitted terms take its exchange fraction and the kernel is added by
  `response_xc_kernel`, the same routine the unfitted operator calls; a
  range-separated hybrid is declined. **Do not check a Kohn-Sham Hessian
  against differenced gradients**: the gradient carries the grid's response
  and the Hessian holds the grid fixed, which is 4e-4 on water and looks like
  a bug. The reference is PySCF's
  `df.hessian` with `auxbasis_response = 2`; at 1 the fit's second-order
  response is missing and it is 1e-4 off. The fitted terms are written
  twice: `df_two_electron_derivatives` in the module is what runs, and
  `dense_two_electron_derivatives` in `test_mqc_czt_df_hessian` is the
  formulas as written, one full array per perturbation. **Change the dense
  one first and make the other agree** -- the test holds them to 1e-9, which
  is every sign and index convention at once. The solver's operator is `df_response_mo`, in the MO
  basis, and only the virtual-occupied block of what it returns is read by
  the solve and only the occupied-occupied block by the assembly.
- `NEO-HF` - nuclear-electronic orbital Hartree-Fock: `keywords.neo` with
  `quantum_nuclei` (0-based indices or symbols) and `nuclear_basis` (default
  `pb4-d`, from `basis_sets/neo/`) gives the named protons orbitals of their
  own; closed-shell HF or DFT, with `epc` (`17-1`/`17-2`) adding the
  electron-proton correlation functional. Validated against PySCF-NEO
  (github.com/theorychemyang/pyscf); see `mqc_docs/source/neo.rst`
- `CCSD`, `CCSD(T)`, `RI-CCSD`, `RI-CCSD(T)` - RHF reference, spin-adapted
  (spatial orbitals) by default. `keywords.cc.spin_adapted: false` selects the
  spin-orbital formulation instead; the two are exact for a closed shell and
  agree to machine precision, so that flag chooses how a number is computed and
  not which number. Spatial is the default because it is roughly sixteen times
  smaller and several times faster.
  **Analytic gradient: RI-CCSD and RI-CCSD(T) only**, closed shell over an
  exact reference, frozen core allowed (`mqc_czt_ri_ccsd_gradient`). The Lambda equations
  (`mqc_czt_rcc_lambda`) and the densities (`mqc_czt_rcc_density`) are PySCF's
  term for term, written with `contract` from `mqc_tensor_contract` so each
  line carries the reference's index string. The two places that dominated
  have been rewritten by hand and the rest left readable: the all-virtual
  block of the two-particle density, and the triples' terms. The gradient
  prints its phases; **read them before optimising anything else** -- a Lambda
  iteration already costs what an amplitude iteration does.
  The relaxation and assembly are `relaxed_fitted_gradient`, shared with
  the RI-MP2 gradient. PySCF has no fitted CCSD gradient, so the reference is
  `tools/cpu_validation/ri_ccsd_gradient.py`, itself held to finite
  differences. Conventional integrals, open shells and a fitted reference are
  refused by name. (T) enters as three sets of terms from
  `mqc_czt_rcc_triples_terms` -- sources in the Lambda equations, a diagonal
  on the one-particle density, additions to the two-particle density -- built
  one set of three virtuals at a time, as the energy is, at two and a half
  to four and a half times the energy's cost (a water tetramer and adenine,
  both in cc-pVDZ). **(T) is not invariant to a
  rotation among the occupied orbitals or among the virtuals**, so the
  gradient takes the off-diagonal relaxed density in those blocks from the
  Lagrangian, divided by the orbital-energy gap, and skips degenerate pairs;
  using the amplitudes' value there is right for CCSD and wrong for (T).

## The terco backend

**The terco backend runs RI-CCSD and RI-CCSD(T) too** (`trc_run_rccsd`,
through `run_coupled_cluster` in `mqc_terco_driver`): closed shell over
Hartree-Fock, energies only. The correlation energy comes back as the CPU
backend's two parts, `t1 t1` and `t2`, through `trc_rccsd_energy_parts`
(terco ABI 1.29), so `energy%cc` reads the same from either backend.
Conventional coupled cluster on that backend is
refused by name. **The terco bridge is not built in CI or in a default
build**: compile it with `-DMQC_ENABLE_TERCO=ON -DTERCO_ROOT=<checkout>`; with
`-DTERCO_LIBRARY` pointed at any file the library target compiles without a
terco build, which checks the interface and nothing else.

## Method names and prefixes

An `RI-`/`DF-` prefix parses to the same method type as the bare name, so the
intent is recovered in the reader by `method_wants_density_fitting` while the
spelling still exists. `ccsd` and `ccsd(t)` are separate method types, so the
triples survive the parse.

## Solvation Models

Through tblite, for xTB only:

- `ALPB` - Analytical Linearized Poisson-Boltzmann
- `GBSA` - Generalized Born with Solvent-Accessible Surface Area

Continuum solvation on the CPU ab initio path (`mqc_czt_pcm.f90`):

- `CPCM` - conductor-like, the scaled-conductor approximation
- `IEFPCM` - the integral-equation formalism

Both work with Hartree-Fock and DFT, restricted and unrestricted; see
`mqc_docs/source/continuum_solvation.rst`. 40+ solvents.

**MP2 in a continuum needs `keywords.pcm.correlation`**, and without it is
refused: `"pte"` solvates the reference and runs any MP2 on its orbitals;
`"ptd"` runs the reference and closed-shell `ri-mp2` in the gas phase, then
solves the continuum once against the relaxed MP2 density and adds that as
`energy%solvation` (`solvation_energy` in the output). The context is built
and then switched off for the SCF (`pcm_ctx%enabled = .false.`), so anything
else that reads it under PTD sees the gas phase. The relaxed density is
`czt_ri_mp2_gradient(..., relaxed_density=)`, which returns before any
derivative integral, with `fit_response` following the reference: **the
gradient's default fitted response is 2.4e-7 hartree off here** and is why an
exact-reference case would miss PySCF. `test_mqc_czt_mp2` holds the density to
a finite-field derivative. PTED needs the surface-charge response in the
Z-vector and is not implemented.
