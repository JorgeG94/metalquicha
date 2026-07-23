# cuEST backend — logbook

Working notes for the `feat/cuest` branch: what is done, what is trusted, what
is next. Read this first when picking the work back up. The reference document
for *how the backend works* is `CUEST.md`; this file is about *state and
history*.

## One-line status

A GPU quantum chemistry backend (NVIDIA cuEST) is wired into metalquicha behind
the `qc_method_t` interface. HF and DFT, energies / gradients / Hessians,
closed and open shell, serial and fragmented, all validated on an H200. The
only major thing implemented-but-never-run is multi-GPU (more than one rank).

Branch: `feat/cuest`, based on `main` at `67cae8f2`. Never pushed as of the
last entry, so history has been rewritten freely (session links stripped);
force-push not yet needed.

## Environment (NCI Gadi)

```sh
module load gcc/13.2.0 cuda      # or gcc/15.1.0; see the gfortran note below
```

- **Compile/link** works on the login node. **Running cuEST** needs sm_80+, so
  it fails on the login node's V100 with `CUEST_STATUS_UNSUPPORTED_ARCHITECTURE`
  (status 11) at `cuestCreate`. That is expected, not a bug.
- **Run on** `gpuhopper` (H100) / an H200 node / `dgxa100` (A100).
- The cuEST package is a sibling: `../libcuest-linux-x86_64-0.2.0.4_cuda12-archive`.

Build:

```sh
FC=mpifort CC=mpicc cmake -B build \
    -DMQC_ENABLE_TBLITE=OFF -DMQC_ENABLE_CUEST=ON \
    -DCUEST_ROOT=/scratch/bm55/jlv900/dev/libcuest-linux-x86_64-0.2.0.4_cuda12-archive
cmake --build build -j
ctest --test-dir build -R mqc          # 29 tests, all host-side, no GPU needed
```

`-DMQC_CUEST_BINDINGS=fetch` clones the bindings from
[mod_cuest](https://github.com/JorgeG94/mod_cuest) instead of using the vendored
copy. `-DMQC_ENABLE_TBLITE=ON` also works but is slower to build.

## How to restart: run these first

Every claim below is backed by one of these. On a fresh GPU node, from the repo
root:

```sh
./validation/run_cuest_validation.sh     # 9 energies: HF/DFT, closed+open shell
./build/check_gradient                   # analytic vs finite-diff gradient (HF)
./build/check_gradient pbe0              #   "" for a functional
./validation/run_fragment_test.sh        # water hexamer MBE(2), 21 fragments
./validation/run_guess_comparison.sh     # core vs GWH vs SAC
./build/mqc validation/inputs/hess_water_hf.mqc   # Hessian, freqs, IR, dipole
```

Expected numbers are in `CUEST.md` (Validation section) and reproduced in the
commit messages. If any energy moves by more than ~1e-9 Ha from the committed
value, something regressed — chase that before doing anything new.

## What is done and TRUSTED (validated on hardware)

| capability | evidence |
|---|---|
| RHF / RKS energy, 20 functionals | 9-row validation sweep, refs in CUEST.md |
| Analytic gradients (RHF/RKS) | 3e-8 Ha/Bohr vs finite difference |
| Hessians, frequencies, thermochem | water freqs 1758/3939/4052, textbook HF |
| Dipole + IR intensities | dipole 2.09 D; sum rule satisfied to 3e-6 a.u. |
| Open shell UHF / UKS | OH doublet, O2 triplet; <S^2> correct; UKS gradient |
| Fragmented energies + gradients | water hexamer MBE(2), 20-mer (157 frags) |
| Three basis formats (json/gbs/gamess) | json vs gbs cross-checked to 1e-10 |
| Initial guesses: core / GWH / SAC | guess comparison table; GWH is default |
| fpm build excludes cuEST cleanly | fpm compiles all of src/, only BLAS link left |

Key validated facts worth not re-deriving:
- `def2-SVP` sits ~0.066 Ha ABOVE `cc-pVDZ` for every method (different O
  primitive sets, core-dominated). Not a bug; cost me an hour once.
- OH with a **core guess** converges cleanly to the WRONG state (A2Sigma+,
  ~4.3 eV up, <S^2> a respectable 0.758). GWH fixes it. This is why the guess
  matters and why GWH is the default.
- O2 at HF gives ~1 eV bond energy vs 5.2 experimental — HF's famous
  multireference failure, reproduced correctly.
- SAC (superposed atomic coefficients) is implemented but was NOT better than
  GWH in the one comparison run — worse on O2. Measured, not assumed. GWH stays
  default.

## What is implemented but UNTESTED on hardware

- **Multi-GPU (more than one rank).** Each rank binds to
  `mod(node_local_rank, device_count)` and logs it
  (`cuEST: node-local rank N bound to GPU N of M`). Never run past one rank.
  First test to do with >1 GPU:
  ```sh
  mpirun -np 4 ./build/mqc validation/inputs/frag_water20_hf.mqc \
      2>&1 | grep -iE "bound to GPU|Final|ERROR"
  ```
  The parallel energy MUST equal a serial run of the same input.
- **SAC on DFT.** The empty-beta-channel fix (needed for the H atom) went in
  after the last SAC run, so `sac_water_pbe` has never completed. Quick check.
- **H-cap redistribution with cuEST.** Water clusters have no broken bonds.
  Fragmenting something covalent (the glycine inputs) would exercise it.

## What is NOT implemented

- Empirical dispersion (requesting it errors, deliberately).
- ECPs, PCM. cuEST supports both; neither wired. BSE JSON already carries ECP
  data, so the parser half is partly there.
- Raman. Design is written up in CUEST.md ("Towards Raman"). Groundwork —
  the dipole via `cuestMultipoleCompute` — is done and validated. Two routes:
  semi-numerical (alpha = dmu/dF by finite field, works for all functionals)
  and analytic CPHF (`cuestDFMOIntegralsCompute` gives the MO-basis 3-index DF
  tensor, so the orbital Hessian is small enough to build and solve directly;
  HF-only, since cuEST exposes the XC potential but never the kernel f_xc).

## Architecture, in one breath

- `src/methods/mqc_method_hf.F90` / `mqc_method_dft.F90` — the `qc_method_t`
  implementations. NO preprocessor conditionals; they `use mqc_cuest_bridge`.
- `src/methods/mqc_cuest_iface.f90` — the settings type, always compiled.
- `src/methods/mqc_cuest_bridge_stub.f90` — no-op bridge (fpm / no-cuEST).
- `cuest_backend/backend/mqc_cuest_bridge.f90` — real bridge (CMake + cuEST).
- `cuest_backend/backend/` — the actual backend (context, integrals, scf,
  gradient, grid, functionals, atomic_guess, driver).
- `cuest_backend/bindings/` — generated cuEST + CUDA Fortran bindings.

The whole `cuest_backend/` tree is OUTSIDE `src/` on purpose: fpm globs `src/`
and cannot link cuEST, and its dependency scanner is preprocessor-blind (it
reads `use` lines through `#ifdef`). The stub/real bridge split is what lets
fpm build the CPU/xTB path with zero cuEST references. Do not move the backend
back under `src/` without understanding this.

## Gotchas that cost real time

- **cuEST has no SCF, no ERI path.** It gives S/T/V, J[D], K[C_occ], Vxc,
  multipoles, and all the `*DerivativeCompute` entry points. J and K are ALWAYS
  density-fitted, so every run needs a JKFIT auxiliary basis. The SCF, DIIS,
  orthogonalization, guess, and energy expression are all ours (host, pic-blas).
- **Host vs device pointers are both `type(c_ptr)`.** Shell exponents/coeffs and
  pair-list coordinates are HOST; every matrix in/out is DEVICE. Wrong = segfault.
- **Workspace query and compute must get identical arguments.** A mismatch (the
  beta-exchange bug, `374bade2`) throws `CUEST_STATUS_EXCEPTION`. There is an
  audit loop in that commit's history if it recurs.
- **DF derivative wants the HALF density** (`densityScale=2`) for RHF, and the
  TOTAL density (`densityScale=0.5`, `coeffScale=-0.5`, two coeff matrices) for
  UKS. Getting it wrong is a clean factor of two the FD gradient catches.
- **Exact-exchange fraction is QUERIED from the XC plan**, never tabulated, so
  hybrids/range-separated hybrids can't drift. Do not hardcode a table.
- **gfortran 13.2.0 miscompiles** intrinsic assignment from a
  `class(..), allocatable` function result (segfault) and silently drops
  settings passed to a `class(..), intent(inout)` dummy. Fixed in 15.1.0.
  The factory uses concrete locals + `allocate(..,source=)` to dodge it. If
  weird polymorphic behaviour reappears, suspect the compiler first.
- **FP64 is pinned** (`CUEST_NATIVE_FP64_MATH_MODE`) so MBE energy differences
  are reproducible. Runs are bit-identical. Don't relax it.
- **Scripted string-replace edits can silently no-op.** Two bugs this session
  came from a `patch()` that didn't match and did nothing. Assert the match.

## Config plumbing — the weak spot

Settings are hand-copied across three hops: input file -> `mqc_config_t` ->
`method_config_t` -> concrete `*_options_t`. Nothing catches a missed copy, and
it silently dropped `basis`, `aux_basis` and `functional` three separate times.
`test/test_mqc_config_roundtrip.f90` guards against it now (verified to fail on
a reintroduced bug). The proper fix — methods holding `method_config_t` directly
and deleting the duplicated `*_options_t` — is still not done; it is the highest-
value refactor if touching this area again.

## Suggested next steps, in order

1. **Multi-GPU** — the one untested capability. Start with the 4-rank command
   above. Likely just works; the binding is in place.
2. **SAC-on-DFT smoke test** — `./build/mqc validation/inputs/sac_water_pbe.mqc`,
   confirm it completes after the empty-beta fix.
3. **Raman** — semi-numerical first (finite field -> alpha -> finite-diff
   dalpha/dR), validate, then analytic CPHF for HF. Design in CUEST.md.
4. **Config refactor** — collapse `*_options_t` into `method_config_t`. Pure
   host, the round-trip test covers it.
5. Cosmetic: suppress IR intensities on trans/rot modes (they're provably
   projection artifacts — sum rule holds to 3e-6); tighten the near-zero
   frequency classifier that flip-flops "2 imaginary" vs "0".

## Commit map (feat/cuest, newest first)

```
b03c1e48 dipole print + dipole-derivative sum-rule check
d8b25a15 keep cuEST out of the fpm build (bridge stub/real split)
9cb81a4f dipole moment + IR intensities
2e3c5c27 empty beta channel in UKS XC (H atom); measured guess table
203e1adc guess comparison script
851b81dd wire SAC guess
6e8a8c19 O2 triplet result
d4c7295b open-shell result + initial-guess lesson
019ce281 BSE JSON basis reader, preferred
ff5ac018 GWH + SAC scaffolding
459c0d6e `unrestricted` keyword (also the UKS-reduces-to-RHF check)
374bade2 fix beta-spin exchange query/compute mismatch
0e1c9229 open-shell UHF/UKS
0cd3adef restructure docs; advertise backend in README
f83c8c09 delete libcint stub; config round-trip test
b2e4b6ce fetch bindings from mod_cuest; semi-numerical Hessians
db17ef59 document the backend (CUEST.md)
9a59b494 log GPU binding; fragmented DFT + gradient inputs
45f15d84 per-rank GPU binding; fragmented tests
d4e897f3 analytic gradients
3e1cda43 cuEST backend: density-fitted HF and DFT
```
