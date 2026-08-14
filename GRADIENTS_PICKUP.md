# Analytic gradients — pickup notes

Written 2026-08-14. For starting a fresh session without re-deriving what was
already settled. Delete when the ladder is finished.

## Where things stand

**Pushed, rebased onto current main, all green:**

| branch | tip | contents |
|---|---|---|
| `feat/cpu-gradients` | `50d750e6d` | libcint pin bump, RHF + UHF gradients, threading/symmetry |
| `feat/cpu-gradients-df` | `b7eb16303` | DF-RHF gradient |
| `feat/cpu-gradients-dft` | `5a76bd7a2` | RKS/UKS **LDA** |
| `feat/cpu-gradients-gga` | `2deae8294` | AO 2nd derivatives, then GGA + hybrid-GGA |

Linear stack, each contains the previous. Ready for stacked PRs.

**Current branch: `feat/cpu-gradients-mp2`**, branched off `-gga`.
Working tree has **uncommitted** harness edits to `validation/check_scf_gradient.f90`
(see "Fix first" below). Nothing else outstanding.

`mqc-geoopt` is a separate, already-pushed stack (the DL-FIND geometry
optimizer). Jorge rebased it onto main and added HDLC residues — refetch
before touching it.

Also pushed separately: `JorgeG94/libcint` `develop` at `3c78069`, which added
the Fortran bindings for `int1e_ipkin/ipnuc/iprinv`, `int3c2e_ip1/ip2`,
`int2c2e_ip1` and their optimizers. mqc pins that commit.

## Done since (branch `feat/cpu-gradients-validation`, tip `7ffa9d80c`)

Branched off `-gga`, was briefly named `-mp2`. Not pushed yet.

- **The gradient is reachable from a deck.** `run_libcint_hf` refused
  `driver: Gradient`; it now calls `libcint_scf_gradient` after the SCF, for
  every combination the SCF supports. Post-HF alongside a gradient is refused.
- **The output JSON carries the components**, Hartree/Bohr, one `[x, y, z]` per
  atom, at all four writer sites. `gradient_norm` stays.
- **`run_validation.py` compares element by element** against
  `expected_gradient` (`gradient_tolerance`, default 1e-7) and reports the
  translational residual for anything with components -- a failure only where
  the deck sets `check_translation`.
- **Ten generated CPU gradient cases**, one per code path, in
  `GRADIENT_CASES` / `GRADIENT_FRAGMENTED_CASES` in `tools/cpu_validation`.
  204 CPU + 20 xTB + 59 unit tests pass.
- **Compiler**: GCC 13 works (`build-g13`), configured with the source dirs
  reused from `/home/jorge/dev/metalquicha/build/_deps`. The `build` and
  `build-xc` trees still refuse to configure -- the shared dftd4 checkout moved
  to v4.2.0, which rejects gfortran 15.1.
- **The check_scf_gradient filter was never broken** -- the binary under test
  was stale. It now also covers the unrestricted DF branch case.

Two things learned that are worth keeping:

- **PySCF's KS gradient omits the grid-weight derivative by default.** Set
  `grid_response = True` or the two codes differ by ~1e-4 with nothing wrong.
- **The SCF stopping point enters a gradient linearly.** At the suite default
  of 1e-8 the cases sat ~1e-7 from PySCF; at 1e-11, 4e-12 to 7e-9.

## Done: MP2 (branch `feat/cpu-gradients-mp2`, tip `8d356e82b`, pushed)

`backends/libcint/mqc_libcint_mp2_gradient.f90`, wired into the bridge, three
generated validation cases, and `validation/check_mp2_gradient.f90` for finite
differences. Agrees with PySCF to 3e-11..2e-9 and with finite difference to
9e-7..5e-6.

Two things from it that generalise:

- **PySCF's MP2 gradient is itself ~4e-8 wrong** on water/cc-pVDZ, because its
  Z-vector goes through `pyscf.scf.cphf.solve`, whose Krylov solver carries
  error `tol` does not control. `dense_zvector` in `tools/cpu_validation` is the
  monkeypatch. Without it the disagreement reads as a missing term on our side.
- **The t2 = 0 limit is the cheapest test there is.** It must reproduce the
  Hartree-Fock gradient exactly, and it isolated two bugs in the shared
  skeleton that the correlation algebra would have masked.

## In progress: RI-MP2 (branch `feat/cpu-gradients-ri-mp2`, off `-mp2`)

**The reference does not exist.** `pyscf.mp.dfmp2.DFMP2.nuc_grad_method` and
`dfmp2_native`'s both `raise NotImplementedError`, and there is no
`pyscf/grad/dfmp2.py` or `pyscf/df/grad/mp2.py`. Every gradient in this ladder
so far was pinned component-wise against PySCF to ~1e-9; this one cannot be.
That changes the discipline rather than lowering it:

1. **Finite difference of mqc's own RI-MP2 energy is now the primary check**,
   not the corroborating one. It is good to ~1e-6 with the step scaling
   verified (halve the step, the deviation must fall by four).
2. **The `t2 = 0` limit must reproduce the Hartree-Fock gradient exactly.** The
   reference here is an exact-ERI RHF -- RI-MP2 decks fit only the correlation
   (`aux_only` in the generator) -- so the whole one-particle skeleton is the
   already-validated one, and this check is free.
3. **Translational invariance**, as always, and worth as little as always.
4. **The conventional MP2 gradient is the limit of this one** as the auxiliary
   basis saturates. Not a tolerance, but the difference should sit at the size
   of the fitting error in the energy (a few tenths of a millihartree on a
   matched RIFIT set), and a wrong dJ term will not.

Write the finite-difference harness first and check it against the
*conventional* MP2 gradient, where the answer is known. A reference that is
itself unvalidated is worse than none.

**What changes and what does not.** The reference is exact, so `doo`, `dvv`,
the Lagrangian's `veff` term, the Z-vector, the relaxed density, the
energy-weighted density, and every one-electron and Hartree-Fock two-electron
term are structurally identical to `mqc_libcint_mp2_gradient`. Three things
change:

- the amplitudes come from `(ia|jb) = sum_P B^P_ia B^P_jb`, which
  `run_libcint_ri_mp2` already builds through `build_df_tensor`;
- the two-particle density never becomes a four-index AO object. It stays as
  `Gamma3^P_ia = sum_jb Gamma_iajb B^P_jb`, which is `n_o n_v n_aux`;
- that object contracts against `three_centre_deriv` (both `which=1` and
  `which=2`) and `two_centre_deriv`, all three already on
  `feat/cpu-gradients-df` and PySCF-validated in isolation.

The derivative to be careful with is the metric's:
`d(ia|jb) = 2 sum [d(ia|P)] c^P_jb - sum c^P_ia dJ_PQ c^Q_jb`, with
`c = J^-1 (·|·)`. The factor of two is the `(ia)<->(jb)` symmetry of `Gamma`,
which holds because `t2_ijab = t2_jiba`. Getting the metric term's sign or
factor wrong is invisible to translational invariance and to a matching energy,
and finite difference is what catches it -- which is the whole reason for
point 1 above.

Memory: `Gamma_iajb` is `n_o^2 n_v^2`, which is the ceiling here and is not the
`n_ao^4` the conventional path had. Block over `i` when it matters.

## Validation discipline — do not skip

Three checks, because they fail differently. This is not ceremony; each one
caught something the others missed.

1. **Finite difference of mqc's own energy.** The only check that cannot be
   fooled by a misunderstanding shared with PySCF. Consider step-scaling: if the
   deviation falls by 4x per halving of the step, the residual is the difference
   formula and not the derivative.
2. **PySCF component-wise**, fed *this repo's* basis JSON via
   `sys.path.insert(0, "tools/cpu_validation"); from gen_cpu_validation import bse_to_pyscf`.
   PySCF's own basis tables differ in the 8th decimal on Pople sets and it looks
   exactly like a gradient bug. Python: `/home/jorge/dev/metalquicha/.venv/bin/python`.
   Note `mf.Gradients()` raises NotImplementedError unless you import the grad
   module explicitly (`from pyscf.grad import rhf as rhf_grad`).
3. **Translational invariance**, `sum_A dE/dR_A = 0`. Cheap, and it caught a
   missing Hellmann-Feynman term in the HF gradient. **But it is structurally
   blind to exchange-correlation terms** — it passed at 2e-15 over an LDA
   weight-derivative error worth 0.77 Hartree/Bohr, and again over a GGA index
   pairing wrong by 6%. Never treat it passing as evidence a new XC term is right.

**A matching total energy proves nothing about a gradient.** In both XC failures
the converged energy agreed with PySCF to all ten printed digits while the
gradient was badly wrong.

**Validate a new integral layer on its own before building on it.** The
three-centre derivatives and the AO second derivatives were each checked
against PySCF in isolation first. An ordering or normalisation error there does
not announce itself in an assembled gradient — it produces a plausible number
that disagrees with finite differences for reasons that could be anywhere.

**Compare error floors, not absolute tolerances.** The AO second derivatives
agree with PySCF to 3.53e-7 relative, and the pre-existing AO *values* agree to
3.52e-7 — the documented basis-normalisation floor. Landing *on* the existing
floor rather than above it is what says a new quantity adds no error of its own.

## Known gaps, deliberately left

- **DF-UHF does not exist.** `run_libcint_uhf` has no `aux` argument and says so
  ("No density fitting here yet"). The unrestricted DF *gradient* branch is
  written and its channel weights are pinned to 9e-16 by driving a closed shell
  through it as two half-filled channels, but genuinely open-shell behaviour is
  unexercised and will stay that way until the SCF exists. Not a basis-set issue.
- **DF gradient memory.** `df_two_electron_gradient` builds full
  `(nao, nao, naux)` intermediates. That is a memory ceiling before it is a speed
  problem — block over the auxiliary index. Belongs as a follow-up commit on
  `feat/cpu-gradients-df`, which means rebasing the branches above it.
- **No prescreening on the exact-ERI gradient.** Every quartet is computed
  however negligible. mqc is at parity with PySCF single-threaded and ~6x faster
  on 40 threads at (H2O)x3/cc-pVDZ, but that advantage narrows and then inverts
  on a large sparse system where PySCF's Schwarz screening changes the scaling.
  The ordinary Schwarz bound is **not** rigorous for a differentiated integral —
  the bound needs rederiving, which is why it was not bolted on.
- **CCSD gradients need Λ equations**, which mqc does not have at all (zero
  occurrences of "lambda" in the CC module). That is its own project, comparable
  in size to the T solver, before a gradient is even started.

## Workflow (Jorge's, follow it)

- One branch per gradient, branched from the branch carrying the machinery it
  builds on — not from main.
- **Commit only after PySCF validation.** Speed work — threading, permutational
  symmetry, screening — goes in a **separate follow-up commit**, never mixed into
  the correctness commit.
- Run `pre-commit` before pushing; cmake-format and fprettify rewrite files, and
  the fixes must be staged *into* the commit (re-check with
  `pre-commit run --from-ref HEAD~1 --to-ref HEAD`).
- No Claude session links in commit messages. `Co-Authored-By` stays.
