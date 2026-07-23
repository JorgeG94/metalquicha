# cuEST backend

GPU-accelerated Hartree-Fock and Kohn-Sham DFT in metalquicha, with integrals
from NVIDIA [cuEST](https://developer.nvidia.com/cuest), usable both for whole
molecules and inside the fragmentation (MBE/GMBE) machinery.

## Building

```sh
module load gcc/13.2.0 cuda
FC=mpifort CC=mpicc cmake -B build \
    -DMQC_ENABLE_TBLITE=OFF \
    -DMQC_ENABLE_CUEST=ON \
    -DCUEST_ROOT=/path/to/libcuest-linux-x86_64-<ver>_cuda12-archive
cmake --build build -j
```

Only the cuEST shared library is needed: the Fortran bindings are
pre-generated, so the headers are not required to build.

The bindings come from the vendored copy under `src/cuest_bindings/` by
default. To pull them from upstream instead:

```sh
cmake -B build -DMQC_CUEST_BINDINGS=fetch -DMQC_CUEST_BINDINGS_TAG=main ...
```

which clones [mod_cuest](https://github.com/JorgeG94/mod_cuest) at configure
time. The sources are identical either way -- the vendored copy is a straight
copy of them. Fetching needs network access when CMake runs, so configure on a
login node.

**cuEST ships GPU code for sm_80 and newer only.** It compiles and links
anywhere, but on an older card -- including the Gadi login node's V100 --
`cuestCreate` returns `CUEST_STATUS_UNSUPPORTED_ARCHITECTURE`. That is expected,
not a defect. Run on `gpuhopper` or `dgxa100`.

## Running

```
%model
method    = dft          ! or hf
basis     = def2-svp
aux_basis = def2-universal-jkfit
functional = pbe         ! dft only
end
```

`driver type = Gradient` gives analytic gradients, and `type = Hessian` gives a
Hessian from central differences of those gradients -- semi-numerical, so only
one derivative is taken numerically and the accuracy is far better than
differencing energies twice. It costs 6N gradient evaluations, i.e. 6N SCFs.

Everything works fragmented: set `%fragmentation` as usual and each subsystem
is evaluated on the GPU.

Twenty functionals are available (SVWN5, BLYP, PBE, M06-L, r2SCAN, B97, B3LYP,
PBE0, M06, M06-2X, HSE06, CAM-B3LYP, LC-wPBE, wB97X, wB97X-V, B97M-V, wB97M-V,
...); an unknown name errors with the full list rather than defaulting.

### Density fitting is mandatory

cuEST exposes no conventional four-index ERI path, so J and K are *always*
density fitted and an auxiliary JKFIT basis is required alongside the orbital
basis. This is a property of the engine, not a configurable approximation, and
it means every energy here carries a small RI fitting error (~1e-4 Ha, visible
in the cc-pVDZ row below).

### Basis sets

`.gbs` (Gaussian94) is preferred, `.txt` (GAMESS `$DATA`) still works. See
`basis_sets/PROVENANCE.md` -- in particular, download `.gbs` files from the
Basis Set Exchange with `uncontract_spdf=1`, since combined `SP` shells are not
supported.

## Validation

Water, r(OH) = 0.959 A, angle 107.3 deg, def2-universal-JKFIT throughout.

| Calculation | Total energy (Ha) | Reference |
|---|---|---|
| HF/STO-3G | -74.962092 | -74.96 |
| HF/cc-pVDZ | -76.026404 | -76.027 |
| HF/def2-SVP | -75.960665 | |
| SVWN5/def2-SVP | -75.795009 | |
| PBE/def2-SVP | -76.271689 | |
| B3LYP/def2-SVP | -76.320709 | |
| PBE0/def2-SVP | -76.275940 | |

def2-SVP sits ~0.066 Ha *above* cc-pVDZ for every method. Both contract to
`[3s2p1d]`, but oxygen's primitive sets differ -- (7s4p1d) vs (9s4p1d) -- and
total energy is dominated by the 1s core. Applying that single offset reproduces
all four DFT references to a few mHa.

**Gradients** against central finite differences (`build/check_gradient`):

| | max deviation | gradient norm |
|---|---|---|
| HF/def2-SVP | 3.5e-08 Ha/Bohr | 5.9e-02 |
| PBE0/def2-SVP | 3.4e-08 Ha/Bohr | 3.1e-02 |

Both sit at the finite-difference noise floor.

**Fragmented**, water hexamer MBE(2), 21 subsystems on one rank
(`validation/run_fragment_test.sh`):

| | 1-body | MBE(2) | 2-body | time |
|---|---|---|---|---|
| HF/def2-SVP | -455.758908 | -455.834051 | -47.2 kcal/mol | 1.0 s |
| PBE/def2-SVP | -457.633788 | -457.756970 | -77.3 kcal/mol | 9.5 s |

Both two-body corrections are the right magnitude for six waters, and PBE
overbinding the hydrogen bonds relative to HF is expected GGA behaviour. A water
20-mer MBE(2) screens 210 dimers to 157 at a 6 A cutoff and runs in seconds.

## Design notes

**One context per MPI rank.** The cuEST handle and all plain device scratch live
in a single per-rank context, created lazily and reused by every fragment that
rank evaluates. Re-creating the handle per fragment would dominate a fragmented
run -- 21 hexamer fragments complete in 1.0 s, which 21 `cuestCreate` calls could
not fit inside. Device scratch grows to a high-water mark and is never shrunk;
because fragments are size-sorted largest-first, in practice it is allocated once
per rank at the largest subsystem.

**Exact exchange is queried, never tabulated.** The XC plan is built first, then
interrogated for `EXCHANGE_SCALE`, `LRC_EXCHANGE_SCALE` and `LRC_OMEGA`, and
those feed the DF plan. Hybrids and range-separated hybrids therefore cannot end
up with mismatched Coulomb and XC operator definitions, and no table of
functional coefficients is maintained here.

**FP64 is pinned.** The handle sets `CUEST_NATIVE_FP64_MATH_MODE` explicitly.
MBE differences fragment energies against each other, so nondeterministic or
reduced-precision integrals would surface directly as noise in the many-body
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

**Gradients.** `dE_nuc + D dH/dR - W dS/dR + dE_JK/dR + dE_xc/dR`, with
`W = 2 sum_i eps_i C_ui C_vi`. No CPHF term is needed -- at convergence the
orbitals are stationary and the Pulay term covers the relaxation. cuEST's DF
derivative wants the *half* density with `densityScale = 2`, not the total
density the SCF carries; the finite-difference check exists largely to catch that.

## Not implemented

- Open shell (UKS). Odd electron counts error out rather than guessing.
- Empirical dispersion. Requesting it errors rather than silently omitting it.
- ECPs, PCM. cuEST supports both; neither is wired up.
- Dipoles, and therefore IR intensities alongside a Hessian. cuEST has the
  multipole entry points; they are not called.

Multi-rank GPU binding is implemented -- each rank binds to
`mod(node_local_rank, device_count)` and logs which device it took -- but has
not yet been exercised beyond one rank.
