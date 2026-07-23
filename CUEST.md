# cuEST backend — GPU-accelerated HF and DFT

GPU quantum chemistry in metalquicha, built on
[NVIDIA cuEST](https://developer.nvidia.com/cuda/cuda-x-libraries/cuest) (v0.2.0)
with native Fortran bindings, analytic gradients, and — unusually — the whole
thing usable *inside* the fragmentation machinery, so every MBE/GMBE subsystem is
evaluated on the GPU.

> **Note:** density fitting is not optional. cuEST exposes no conventional
> four-index ERI path, so J and K are always fitted and every calculation needs a
> JKFIT auxiliary basis alongside the orbital basis. Energies therefore carry a
> small RI error (~1e-4 Ha, visible in the cc-pVDZ row below). This is a property
> of the engine, not a switch.

## Features

- **Full SCF**: RHF/RKS closed-shell and UHF/UKS open-shell, with density fitting (RI-JK)
- **DIIS**: commutator (FDS−SDF) DIIS in the orthogonal basis, configurable subspace
- **Gradients**: analytic, assembled from all five cuEST derivative APIs, restricted and unrestricted
- **Hessians**: central differences of the analytic gradients, plus frequencies and thermochemistry
- **Fragmented**: energies and gradients through MBE/GMBE, one GPU context per rank
- **Functionals**: 20 built-ins — SVWN5, BLYP, PBE, M06-L, r2SCAN, B97, B3LYP, B3LYP1, PBE0, M06, M06-2X, HSE06, CAM-B3LYP, LC-ωPBE, LC-ωPBEh, ωB97X, ωB97X-V, B97M-V, ωB97M-V, HF
- **Hybrids without hardcoding**: exact-exchange fractions are queried from cuEST's XC plan, never tabulated here
- **Validated**: reference energies across three basis sets and four functionals; gradients to 3e-8 Ha/Bohr against finite differences; water frequencies to textbook HF accuracy

## Requirements

- **cuEST 0.2.0** shared library (headers not needed — the bindings are pre-generated)
- **NVIDIA GPU**, compute capability 8.0+
- **CUDA Toolkit 12.x**
- **CMake 3.22+**, a Fortran compiler, MPI, BLAS/LAPACK

cuEST ships GPU code for sm_80 and newer only. It compiles and links anywhere,
but on an older card — including the Gadi login node's V100 — `cuestCreate`
returns `CUEST_STATUS_UNSUPPORTED_ARCHITECTURE`. That is expected, not a defect.

## Quick Start

### 1. Get the cuEST SDK

Download from [NVIDIA cuEST Downloads](https://developer.nvidia.com/cuest-downloads)
and unpack it. Only `lib/` is needed to build.

### 2. Build

```bash
module load gcc/13.2.0 cuda
FC=mpifort CC=mpicc cmake -B build \
    -DMQC_ENABLE_TBLITE=OFF \
    -DMQC_ENABLE_CUEST=ON \
    -DCUEST_ROOT=/path/to/libcuest-linux-x86_64-<ver>_cuda12-archive
cmake --build build -j
```

The Fortran bindings come from the vendored copy under `src/cuest_bindings/`.
To pull them from upstream instead:

```bash
cmake -B build -DMQC_CUEST_BINDINGS=fetch -DMQC_CUEST_BINDINGS_TAG=main ...
```

which clones [mod_cuest](https://github.com/JorgeG94/mod_cuest) at configure time.
The sources are identical either way. Fetching needs network access when CMake
runs, so configure on a login node.

### 3. Basis sets

Three formats are read, and `find_basis_file` takes the first that exists:
`.json` (Basis Set Exchange), then `.gbs` (Gaussian94), then `.txt` (GAMESS
`$DATA`). Provenance for the bundled sets is in `basis_sets/PROVENANCE.md`.

**JSON is the one to use.** It keeps a combined shell's coefficient sets
separate, so an `SP` shell splits cleanly into an S and a P sharing exponents:

```json
"angular_momentum": [0, 1],
"exponents":    ["...", "...", "..."],
"coefficients": [[s...], [p...]]
```

Gaussian94 cannot express that, so `.gbs` files must be downloaded pre-split
with `uncontract_spdf=1` or the reader rejects them. JSON also carries ECP data
in the same file.

```bash
curl "https://www.basissetexchange.org/api/basis/<name>/format/json/?elements=H,C,N,O" \
    -o basis_sets/<name>.json
```

The two readers are cross-checked against each other in the test suite: parsing
def2-SVP from `.json` and from `.gbs` must give identical exponents and
coefficients to 1e-10.

Orbital ↔ auxiliary pairings used by the validation inputs:

| Orbital | Auxiliary |
|---------|-----------|
| STO-3G, cc-pVDZ, def2-SVP | def2-universal-jkfit |

### 4. Run

```bash
python3 mqc_prep.py validation/inputs/dft_water_pbe0.json
./build/mqc validation/inputs/dft_water_pbe0.mqc
```

## Usage

The `%model` section selects the backend; everything else is ordinary metalquicha
input.

```
%model
method     = dft                    ! or hf
basis      = def2-svp
aux_basis  = def2-universal-jkfit   ! required
functional = pbe0                   ! dft only
end

%driver
type = Gradient                     ! Energy | Gradient | Hessian
end

%scf
maxiter   = 100
tolerance = 1e-08
end
```

Add a `%fragmentation` section as usual and every subsystem runs on the GPU — no
separate configuration.

**Notes**

- An unknown functional name errors with the full list of accepted names rather
  than falling back to a default.
- Open shell is automatic: `multiplicity /= 1`, or an odd electron count,
  selects UHF/UKS. `<S^2>` is reported so spin contamination is visible.
- `Hessian` costs 6N gradient evaluations, each a full SCF.
- Requesting empirical dispersion errors; it is not implemented.
- Each MPI rank binds to `mod(node_local_rank, device_count)` and logs which
  device it took.

## Validation

### Energies

Water, r(OH) = 0.959 Å, angle 107.3°, def2-universal-JKFIT throughout.

| Calculation | cuEST (Ha) | Reference |
|---|---|---|
| HF/STO-3G | -74.962092 | -74.96 |
| HF/cc-pVDZ | -76.026404 | -76.027 |
| HF/def2-SVP | -75.960665 | |
| SVWN5/def2-SVP | -75.795009 | |
| PBE/def2-SVP | -76.271689 | |
| B3LYP/def2-SVP | -76.320709 | |
| PBE0/def2-SVP | -76.275940 | |

The cc-pVDZ row differs from the un-fitted reference by 6e-4 Ha, which is the
expected RI fitting error.

def2-SVP sits ~0.066 Ha *above* cc-pVDZ for every method. Both contract to
`[3s2p1d]`, but oxygen's primitive sets differ — (7s4p1d) vs (9s4p1d) — and the
total energy is dominated by the 1s core. Applying that single offset reproduces
all four DFT references to a few mHa.

```bash
./validation/run_cuest_validation.sh
```

### Gradients

Against central finite differences, water/def2-SVP:

| | max deviation | gradient norm |
|---|---|---|
| HF | 3.5e-08 Ha/Bohr | 5.9e-02 |
| PBE0 | 3.4e-08 Ha/Bohr | 3.1e-02 |

Both sit at the finite-difference noise floor.

```bash
./build/check_gradient          # Hartree-Fock
./build/check_gradient pbe0     # any functional
```

### Frequencies

Water, HF/def2-SVP, Hessian from central differences of the analytic gradients:

| Mode | calc (cm⁻¹) | expt | ×0.90 |
|---|---|---|---|
| bend | 1758.3 | 1595 | 1582 |
| symmetric stretch | 3939.3 | 3657 | 3545 |
| asymmetric stretch | 4052.5 | 3756 | 3647 |

Six near-zero modes (translations and rotations) plus three real vibrations. The
8–10% overestimate is what Hartree-Fock does — no correlation, harmonic
approximation — and the standard ~0.90 scaling factor recovers experiment.

### Open shell

UHF/def2-SVP:

| system | multiplicity | energy (Ha) | `<S^2>` (exact) |
|---|---|---|---|
| OH radical | 2 | -75.325060 | 0.7548 (0.75) |
| O2 | 3 | -149.490204 | 2.0338 (2.00) |

OH implies D_e(H-OH) = 86.5 kcal/mol against an experimental ~125, the right
size for Hartree-Fock underbinding. O2 implies ~1 eV against an experimental
5.2 eV, which is Hartree-Fock's best-known failure: O2 has multireference
character a single determinant cannot capture. Reproducing a known failure
correctly is as much a check as reproducing a known success.

This one is also a cautionary tale about initial guesses. With a core guess the
same calculation converges tidily, with a respectable `<S^2>` of 0.7584, to
-75.167485 -- which is 4.29 eV higher, against OH's experimental A2Sigma+ <-
X2Pi gap of 4.05 eV. A wrong occupation can be perfectly self-consistent, so
the SCF settles onto the excited state and reports success. GWH fixes the
sigma/pi ordering and lands on the ground state.

### Fragmented

Water hexamer, MBE(2), 21 subsystems on one rank:

| | 1-body | MBE(2) | 2-body | time |
|---|---|---|---|---|
| HF/def2-SVP | -455.758908 | -455.834051 | -47.2 kcal/mol | 1.0 s |
| PBE/def2-SVP | -457.633788 | -457.756970 | -77.3 kcal/mol | 9.5 s |

Both two-body corrections are the right magnitude for six waters, and PBE
overbinding the hydrogen bonds relative to HF is expected GGA behaviour. A water
20-mer MBE(2) screens 210 dimers to 157 at a 6 Å cutoff and runs in seconds.

```bash
./validation/run_fragment_test.sh
```

## Architecture

```
src/cuest_bindings/            # Generated bindings -- do not hand-edit
├── cuest.f90                  # 129 functions, 289 enums, 2 workspace types
├── cuest_helpers.f90          # Typed wrappers over the void*+size_t API
├── cuda_runtime.f90           # CUDA runtime, compiler-agnostic
└── cuda_helpers.f90           # Error decoding, typed array copies

src/methods/cuest/
├── mqc_cuest_runtime.f90      # Status -> error_t, workspaces, device memory
├── mqc_cuest_context.f90      # Per-rank handle, GPU binding, scratch pools
├── mqc_cuest_basis.f90        # molecular_basis_type -> cuEST AO shells
├── mqc_cuest_grid.f90         # Treutler-Ahlrichs M4 radial x Lebedev angular
├── mqc_cuest_functionals.f90  # Functional name -> cuEST identifier
├── mqc_cuest_integrals.f90    # Per-molecule objects; S/T/V, J, K, Vxc, derivatives
├── mqc_cuest_scf.f90          # Closed-shell SCF: orthogonalization, DIIS, energy
├── mqc_cuest_gradient.f90     # Gradient assembly from the derivative APIs
└── mqc_cuest_driver.f90       # physical_fragment_t -> converged result

src/basis/
├── mqc_gbs_reader.f90         # Gaussian94 (.gbs) reader
└── mqc_basis_normalization.f90   # Shell contraction normalization

src/methods/
├── mqc_method_hf.F90          # Hartree-Fock, cuEST-backed
├── mqc_method_dft.F90         # Kohn-Sham, cuEST-backed
└── mqc_semi_numerical_hessian.f90  # Method-agnostic finite-difference Hessian

validation/
├── run_cuest_validation.sh    # Energy references, HF + DFT
├── run_fragment_test.sh       # Water hexamer MBE(2)
└── check_gradient.f90         # Analytic vs finite-difference gradients
```

### Design notes

**One context per MPI rank.** The cuEST handle and all plain device scratch live
in a single per-rank context, created lazily and reused by every fragment that
rank evaluates. Re-creating the handle per fragment would dominate a fragmented
run — 21 hexamer fragments complete in 1.0 s, which 21 `cuestCreate` calls could
not fit inside. Device scratch grows to a high-water mark and is never shrunk;
because fragments are size-sorted largest-first, in practice it is allocated once
per rank at the largest subsystem.

**Exact exchange is queried, never tabulated.** The XC plan is built first, then
interrogated for `EXCHANGE_SCALE`, `LRC_EXCHANGE_SCALE` and `LRC_OMEGA`, and those
feed the DF plan. Hybrids and range-separated hybrids therefore cannot end up with
mismatched Coulomb and XC operator definitions.

**FP64 is pinned.** The handle sets `CUEST_NATIVE_FP64_MATH_MODE` explicitly. MBE
differences fragment energies against each other, so reduced-precision or
nondeterministic integrals would surface directly as noise in the many-body
corrections. Repeat runs are bit-identical.

**Energy expression.** With K carrying the functional's exchange fraction and D
the total density,

```
F = H + J - K + Vxc
E = tr(D H) + 1/2 tr(D J) - 1/2 tr(D K) + Exc
```

which covers HF as the case `Vxc = Exc = 0`. Note `Exc` enters directly, not as
`1/2 tr(D Vxc)`: the familiar `1/2 tr(D (H+F))` form is correct for HF but double
counts XC.

**Gradients.**

```
dE/dR = dE_nuc/dR + D dH/dR - W dS/dR + dE_JK/dR + dE_xc/dR
```

with `W = 2 Σ_i ε_i C_ui C_vi` the energy-weighted density. No CPHF term is needed
— at convergence the orbitals are stationary and the Pulay term covers the
relaxation. cuEST's DF derivative differentiates `E_JK = s_D E_J[D] + s_C E_K[C]`
and for closed-shell RKS wants the **half** density with `s_D = 2`, `s_C = -1`,
not the total density the SCF carries.

### Open shell

Unrestricted follows cuEST's own spin conventions:

```
D^a = C_a C_a^T,  D^b = C_b C_b^T,  D^t = D^a + D^b   (no factor 2)
F^a = H + J[D^t] - K[C_a] + Vxc_a
E   = tr(D^t H) + 1/2 tr(D^t J) - 1/2 (tr(D^a K_a) + tr(D^b K_b)) + Exc
```

which reduce to the restricted expressions when `D^a = D^b = D/2` — the check
that the factors are right. DIIS stacks the two spin commutators into one error
vector. The DF derivative takes the two coefficient matrices concatenated with
`densityScale = 0.5` and `coefficientScale = -0.5`, against the restricted
case's half density with `densityScale = 2`.

The initial guess is the core Hamiltonian for both spins; the differing alpha
and beta occupations are what break the symmetry. A converged atomic guess
would both cut iterations and help difficult radicals.

## Not implemented

- Empirical dispersion. Requesting it errors rather than silently omitting it.
- Dipoles, and therefore IR intensities alongside a Hessian. cuEST has the
  multipole entry points; they are not called.
- ECPs, PCM. cuEST supports both; neither is wired up.

Multi-rank GPU binding is implemented and logged, but has not yet been exercised
beyond one rank.
