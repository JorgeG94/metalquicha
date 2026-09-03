# Metalquicha - agent context

`CLAUDE.md` imports this file, so there is one copy and any agent that looks for
either name finds it.

## Project Overview

**Metalquicha** (Huastec/Tenek word for "sunflower") is a modern Fortran quantum chemistry framework specializing in fragment-based calculations using Many-Body Expansion (MBE) and Generalized Many-Body Expansion (GMBE) methods.

- **Language**: Fortran 2003+ (modern OOP features)
- **Build Systems**: CMake (primary), FPM (secondary)
- **Parallelization**: MPI + OpenMP
- **QC Engine**: tblite (GFN1-xTB, GFN2-xTB)
- **Code Documentation**: https://jorgeg94.github.io/metalquicha/
- **User Facing Documentation**: https://metalquicha.readthedocs.io/en/latest/

## Quick Start

```bash
# Build (a preset, which is the shortest way to a correct configuration)
cmake --preset default && cmake --build --preset default
ctest --preset mqc

# Or by hand
cmake -B build && cmake --build build -j

# Run (serial)
./mqc input.json

# Run (parallel)
mpirun -np 4 ./mqc input.json

# Run tests
ctest -R "mqc"
```

## Directory Structure

```
metalquicha/
├── app/main.f90                 # Entry point
├── src/
│   ├── core/                    # Data types, constants, elements
│   ├── io/                      # Config parsing, XYZ reading, JSON output
│   ├── fragmentation/           # MBE/GMBE engines, fragment generation
│   ├── methods/                 # QC methods (xTB, HF placeholder)
│   ├── vibrational/             # Frequency analysis, thermochemistry
│   ├── parallel/                # MPI utilities
│   ├── basis/                   # Basis set handling
│   ├── utils/                   # Error handling, finite differences
│   └── mqc_driver.f90           # Main calculation dispatcher
├── test/                        # Unit tests (test-drive framework)
├── validation/                  # Physics validation test cases
├── cmake/                       # CMake configuration
├── basis_sets/                  # Basis set data files
└── mqc_docs/                    # Sphinx documentation source
```

## Key Source Files

| File | Purpose |
|------|---------|
| `app/main.f90` | Entry point, MPI init, input parsing |
| `src/mqc_driver.f90` | Routes to fragmented/unfragmented workflows |
| `src/io/mqc_json_config_reader.f90` | Parses JSON input decks |
| `src/io/mqc_json_schema.f90` | Validates a deck before it is read |
| `src/io/mqc_config_types.f90` | `mqc_config_t` and the types it is built from |
| `src/fragmentation/mqc_mbe.f90` | Core MBE implementation (62KB) |
| `src/fragmentation/mqc_physical_fragment.f90` | Fragment representation, H-capping |
| `src/fragmentation/mqc_gmbe_utils.f90` | GMBE with PIE (overlapping fragments) |
| `src/methods/mqc_method_xtb.f90` | xTB via tblite library |
| `src/core/mqc_result_types.f90` | Result containers (energy, gradient, hessian) |
| `src/vibrational/mqc_vibrational_analysis.f90` | Frequency calculations |
| `src/vibrational/mqc_thermochemistry.f90` | Thermochemistry (RRHO) |

## Input Format (JSON)

Read by `mqc_json_config_reader.f90` straight into `mqc_config_t`:

```json
{
  "schema": {"name": "example", "version": "1.0"},
  "model": {"method": "gfn2"},
  "driver": "Energy",
  "molecules": [{
    "xyz": "water.xyz",
    "molecular_charge": 0,
    "molecular_multiplicity": 1
  }],
  "keywords": {
    "fragmentation": {
      "method": "MBE",
      "level": 2,
      "cutoff_method": "distance",
      "cutoffs": {"dimer": 5.0, "trimer": 4.0}
    }
  }
}
```

A molecule gives either `xyz` (a path, resolved relative to the deck) or
`symbols` plus a flat `geometry` list. Atom indices are 0-based. Bonds listed
in `connectivity` are marked broken automatically when their two atoms fall in
different fragments -- that is derived, not declared.

The `.mqc` format and its `mqc_prep.py` generator were removed in 0.2.0; see
`mqc_docs/source/input_files.rst` for the migration table.

## Core Concepts

### Many-Body Expansion (MBE)
- Breaks system into fragments (monomers, dimers, trimers...)
- E_MBE = ΣE_I - ΣΔE_IJ + ΣΔE_IJK - ...
- Truncate at level N for O(N^k) scaling

### Generalized MBE (GMBE)
- Handles overlapping fragments
- Uses Principle of Inclusion-Exclusion (PIE)
- Automatically computes intersection terms

### Cutting a covalent bond

`keywords.fragmentation.bond_breaking` is `"none"` by default and refuses a
partition that severs one. `"afo"` detaches it with an adjusted frozen orbital
instead: a model system around the bond is solved and localized, the orbital on
the bond is reduced to the detached atom's own basis functions, and that hybrid
is frozen -- empty in the fragment that gets nothing of the bond, occupied in
the one that gets all of it. Restricted to `embedding = "none"` for now.

The bookkeeping rule that matters: **a cut belongs to a group, not a fragment.**
A bond severed between two monomers is whole again inside the dimer holding both
ends, so that dimer carries no ghost, no frozen orbital and no electron shift
there. Deciding it per fragment and inheriting it is what made an earlier capped
version 11 Hartree wrong. See `mqc_czt_afo.f90` and `mqc_docs/source/fmo.rst`.

### Hydrogen Capping
- Broken covalent bonds get H-cap atoms
- Forces/Hessian elements redistributed to original heavy atoms
- Implemented in `mqc_physical_fragment.f90`

### Distance Screening
- Filter fragments by inter-monomer distance
- Per-level cutoffs (dimer, trimer, etc.)
- 5-10x speedup for large systems

## Key Types

```fortran
! Main result container (src/core/mqc_result_types.f90)
type :: calculation_result_t
  type(energy_t) :: energy           ! SCF, MP2, CC energies
  real(dp), allocatable :: gradient(:,:)
  real(dp), allocatable :: hessian(:,:)
  real(dp) :: dipole(3)
  logical :: has_energy, has_gradient, has_hessian, has_error
end type

! Fragment representation (src/fragmentation/mqc_physical_fragment.f90)
type :: physical_fragment_t
  integer :: n_atoms
  integer, allocatable :: element_numbers(:)    ! Atomic numbers (Z)
  real(dp), allocatable :: coordinates(:,:)     ! (3, n_atoms) in Bohr
  integer :: charge, multiplicity, nelec
  integer :: n_caps                             ! Number of H-caps
  integer, allocatable :: cap_replaces_atom(:)  ! Original atom each cap replaces
  integer, allocatable :: local_to_global(:)    ! Map to system indices
end type

! System geometry (src/fragmentation/mqc_physical_fragment.f90)
type :: system_geometry_t
  integer :: n_monomers, total_atoms
  integer, allocatable :: element_numbers(:)
  real(dp), allocatable :: coordinates(:,:)     ! (3, total_atoms) in Bohr
  integer :: charge, multiplicity
end type

! QC method base class (src/methods/mqc_method_base.f90)
type, abstract :: qc_method_t
contains
  procedure(calc_energy_interface), deferred :: calc_energy
  procedure(calc_gradient_interface), deferred :: calc_gradient
  procedure(calc_hessian_interface), deferred :: calc_hessian
end type
```

## Calculation Types

- `CALC_TYPE_ENERGY` - Single-point energy
- `CALC_TYPE_GRADIENT` - First derivatives
- `CALC_TYPE_HESSIAN` - Second derivatives + vibrational analysis

## Methods Available

Semi-empirical, via tblite:

- `GFN1` - GFN1-xTB (faster, older parametrization)
- `GFN2` - GFN2-xTB (recommended, more accurate)

Ab initio on the CPU is **cenzontle** (`backends/cenzontle/`, modules `mqc_czt_*`,
built by `-DMQC_ENABLE_CZT=ON`) -- named for the bird of four hundred voices,
because one backend speaks HF, DFT, MP2, CC, MCSCF, EFP and SAPT.

It is **not** a wrapper around either integral library, which is why it stopped
being called `libcint`. The integrals come from **libfint** by default -- an
all-Fortran port of libcint, so no C compiler -- and `-DMQC_USE_LIBFINT=OFF`
takes libcint itself. Both export the same `libcint_fortran` module and the same
`libcint_*`/`LIBCINT_*` symbols, so libfint drops in without an ABI change and no
source under `backends/cenzontle/` knows which one it linked. Those spellings are
therefore the real library's and are left alone; only our own names are `czt`.
The two are not interchangeable at compile time and only one CI job builds
libcint:

- `HF` - restricted and unrestricted; direct, in-core, or density-fitted
- `DFT` - restricted and unrestricted Kohn-Sham over libxc, through
  `model.functional`. LDA through double hybrid, including range-separated
  hybrids, whose long-range exchange needs the direct build
- `MP2`, `SCS-MP2`, `SOS-MP2`, `RI-MP2` - RHF reference only
- `CCSD`, `CCSD(T)`, `RI-CCSD`, `RI-CCSD(T)` - RHF reference, spin-adapted
  (spatial orbitals) by default. `keywords.cc.spin_adapted: false` selects the
  spin-orbital formulation instead; the two are exact for a closed shell and
  agree to machine precision, so that flag chooses how a number is computed and
  not which number. Spatial is the default because it is roughly sixteen times
  smaller and several times faster.

An `RI-`/`DF-` prefix parses to the same method type as the bare name, so the
intent is recovered in the reader by `method_wants_density_fitting` while the
spelling still exists. `ccsd` and `ccsd(t)` are separate method types, so the
triples survive the parse.

Initial guess is `keywords.scf.guess`: `core`, `gwh`, `sac`, `sad`, or `auto`.
`auto` resolves per backend - `sad` on the CPU path, `gwh` on cuEST - so the two
can differ without either knowing what the other chose.

## Solvation Models

Through tblite, for xTB only:

- `ALPB` - Analytical Linearized Poisson-Boltzmann
- `GBSA` - Generalized Born with Solvent-Accessible Surface Area

Continuum solvation on the CPU ab initio path (`mqc_czt_pcm.f90`):

- `CPCM` - conductor-like, the scaled-conductor approximation
- `IEFPCM` - the integral-equation formalism

Both work with Hartree-Fock and DFT, restricted and unrestricted; see
`mqc_docs/source/continuum_solvation.rst`. 40+ solvents.

## Dependencies

| Library | Purpose |
|---------|---------|
| tblite | xTB quantum chemistry engine |
| pic | Utility library (sorting, logging, timers) |
| pic-mpi | MPI wrapper layer |
| pic-blas | BLAS/LAPACK interface |
| test-drive | Unit testing framework |
| jsonfortran | I/O |

## Testing

```bash
# Run all tests
cd build && ctest

# Run specific test
./test_mqc_mbe

# Run the validation suite -- from inside validation/, every default is relative
cd validation && python3 run_validation.py
```

Test files follow pattern `test/test_mqc_*.f90`.

## Build System

`CMAKE_STYLE.md` is the guide. The short version: the top-level `CMakeLists.txt`
is an outline of four `include()`s, every dependency is declared in its own
`cmake/modules/Find<pkg>.cmake` through `mqc_fetch`, every option lives in
`cmake/MqcOptions.cmake`, and validation programs go through
`mqc_add_validation_program`. `CMakePresets.json` carries the configurations
worth remembering, `perlmutter` among them.

## Coding Conventions

See `FORTRAN_STYLE.md` for the complete style guide, and `CMAKE_STYLE.md`
for the build system. Key points:

- Module naming: `mqc_*.f90`
- Type naming: `*_t` suffix (e.g., `calculation_result_t`)
- Internal units: Bohr (length), Hartree (energy)
- Use `to_bohr()` / `to_angstrom()` for conversions
- Error handling via `error_t` type with `has_error()` / `get_message()`
- **Linter**: Uses `fortitude check` - all `use` statements must have `only` clause
- **CI/CD**: GitHub Actions workflows in `.github/workflows/` (build, test, coverage via Codecov)
- **Environment**: Conda environment defined in `environment.yml`

**Forbidden**: `GOTO`, arithmetic IF, implicit SAVE, COMMON blocks, EQUIVALENCE, fixed-form source, assumed-size arrays

### Comments

`!!` documents a Fortran entity, `!` is an ordinary comment. Nothing else --
`mqc_lint.py` enforces this and now runs **strict**, with no `--allow-predoc`,
so `!>` is refused rather than tolerated.

FORD's `!!` is *post-doc*: it documents what **precedes** it. A block therefore
goes on the lines after a `module`/`type`/`subroutine`/`function` statement, and
after a declaration -- never above one. Converting a block that sits above a
declaration by swapping `!>` for `!!` reattaches it to whatever is above; move
it instead. `ford.md` renders `src` and `app` only, so elsewhere the placement
is a convention for readers rather than something the published docs depend on.

A plain `!` line inside a `!!` block splits the docstring and FORD renders only
the first half. Put a marker after the complete block, never inside it.

A docstring says what the thing is and the facts a caller cannot see from the
signature -- shape, units, index base, when it is allocated, which way a sign
runs. Design rationale, counterfactuals, worked examples and history do not
belong there; if a note would only matter to whoever edits *this file*, demote
it to `!`.

`! TODO(mqc): <what is wrong, and the consequence>` is the only marker form.
Plain `!`, never `!!`, so it stays out of the rendered docs. Prefer recording a
defect this way over fixing it inside an unrelated change.

`test/` and `validation/` are deliberately left verbose; do not trim them.

## MPI Architecture

- **Global coordinator**: Distributes fragments across nodes
- **Node coordinator**: Manages workers within a node
- **Workers**: Compute individual fragments
- Tags defined in `src/parallel/mqc_mpi_tags.f90`
- **Hessian parallelization**: Finite difference displacements are distributed across MPI ranks, even for unfragmented systems

## Common Workflows

### Add a new QC method
1. Create `src/methods/mqc_method_newmethod.f90`
2. Extend `qc_method_t` base class
3. Implement `compute` procedure
4. Register in method factory

### Add a new input keyword
Check first whether it already exists. `keywords.scf.density_fitting`,
`keywords.scf.unrestricted` and most of `correlation_config_t` and `cc_config_t`
were already present and working when they were asked for.

Keyword objects: `keywords.scf` for the reference, `keywords.correlation` for what
every post-HF method shares, `keywords.cc` for what only an iterative one needs.

1. Add the field to the relevant type in `mqc_config_types.f90`, with its default
2. Add the key to its object's allow-list in `mqc_json_schema.f90` — the
   validator rejects anything it does not know, so a key added anywhere else
   first will be refused before the reader ever sees it
3. Read it in `mqc_json_config_reader.f90` (`optional_*` leaves the default alone
   when the key is absent, so there is no second table of defaults)
4. Handle in `mqc_config_adapter.f90`
5. Add a case to `test/test_mqc_json_reader.f90`

### Add a new test
1. Create `test/test_mqc_feature.f90`
2. Use test-drive framework
3. Add to CMakeLists.txt via `ADDTEST(mqc_feature)`

## Validation Tests

Decks live under `validation/inputs/<hardware>/<engine>/<method>/`, with the
geometries shared at `validation/inputs/sample_inputs/` and reached from a deck
with `../../../`. The CPU suite is **generated** - edit
`tools/cpu_validation/gen_cpu_validation.py` and rerun it; a deck added by hand
under `cpu/mqc/` is deleted by the next regeneration.

Validation categories under `cpu/mqc/` are `ccsd dft dh efp fmo fukui gradient
hessian mcscf mp2 pcm quao rhf sapt udft uhf` -- worth checking against before
assuming something is untested.

Examples, all under `validation/inputs/cpu/tblite/gfn1/`:
- `h3o.json` - Unfragmented hydronium
- `prism.json` - Water prism MBE(2)
- `overlapping_gly3.json` - Glycine tripeptide GMBE(1)
- `w20_isomer.json` - Water 20-mer MBE(2)

The CPU ab initio suite is the bulk of ~305 decks under `cpu/mqc/`: RHF across H-Ar
and eight basis sets, UHF, density fitting, MP2, RI-MP2, CCSD(T) and RI-CCSD(T),
and one fragmented case. References are PySCF fed this repository's own basis
JSON, not PySCF's internal tables - those differ in the eighth decimal on Pople
sets, which looks exactly like a bug in whichever code you are checking.

### Running it

`validation/run_validation.py`, and it has to be run **from inside
`validation/`** -- `--exe` defaults to `../build/mqc` and `--validation-dir` to
`.`, so it finds nothing from the repository root. There are two manifests and
the default is the smaller one:

```bash
cd validation
python3 run_validation.py                                       # validation_tests.json, 23 cases
python3 run_validation.py --manifest validation_tests_cpu.json  # 299 cases, the ab initio suite
python3 run_validation.py -t water_sto-3g                       # substring filter on the test name
python3 run_validation.py --mpi --np 4                          # fragmented on 4 ranks
```

**A build without tblite skips most of the default manifest, and skips are not
failures.** `-t h3o` on such a build runs zero tests and prints a green summary.
Check the "Total tests" line, not the colour. `validation_tests_cpu.json` is the
one that exercises a libfint/libcint-only build.

The runner deletes every `output_*.json` in `validation/` before it starts.
Those are working outputs rather than references, but do not leave anything you
want under that name.

## Output Format

JSON output (`output_<basename>.json`) contains:
- Total energy and per-fragment breakdown
- Gradients (if requested)
- Hessian (if requested)
- Vibrational frequencies and thermochemistry
- GMBE coefficients and intersection terms

## Performance

Four things worth knowing before optimising anything here, all of them measured.
A benchmark suite that checks them lives on `perf/benchmark-suite`.

* **The exchange-correlation quadrature dominates a DFT run** -- 89 per cent of
  an LDA one and 98 of a meta-GGA one, against a Fock build under two per cent.
  Threading it is worth 5-10x, and until that lands a pure GGA costs about 26x
  Hartree-Fock, which is backwards for a functional carrying no exact exchange.
* **Only xTB is pinned to one thread, and not for speed.** `needs_serial_execution`
  in `src/mqc_method_types.f90` names GFN1 and GFN2 and nothing else: threaded,
  tblite corrupts a result rather than failing. Every other method -- the whole
  cenzontle ab initio path -- keeps the threads the launcher gave it, on
  the unfragmented path and on the fragment workers alike, and a Fock build
  threads its own quartet loop. Clamping those was an old bug, not a policy:
  raising `OMP_NUM_THREADS` did nothing and the cause was nowhere near the
  command line. The measured "threaded MKL competes with the ranks" result --
  five times faster on one molecule, thirty-one per cent slower on four --
  was taken on that clamped xTB path; do not carry it over to ab initio.
* **Structure is a poor guide to where time goes.** Two separate optimisation
  attempts here were aimed at the wrong file on the strength of reading the
  code; one VTune hotspots collection settled it in minutes. `perf_event_paranoid`
  is 4 on the development box, so PMU counters are unavailable -- use
  `-collect hotspots -knob sampling-mode=sw`.
* **Repeat before believing a timing.** DFT energies repeat to under one per
  cent; gradients used to vary by sixteen, and a single-shot comparison across
  that produced a phantom twenty per cent regression during this work.

## Compiler Support

| Compiler | Status |
|----------|--------|
| gfortran 11+ | Full support |
| Intel ifort | Full support |
| Intel ifx | Full support |
| nvfortran | Partial (no tblite) |
| LLVM Flang | Partial (no tblite) |

## Useful Commands

```bash
# Rebuild after changes
cmake --build build -j

# Clean build
rm -rf build && cmake -B build && cmake --build build -j

# Generate coverage report (requires lcov, GCC only)
cmake -DCMAKE_BUILD_TYPE=Coverage-mqc -B build && cmake --build build && cmake --build build --target coverage

# Run with verbose output -- mqc has no such flag; it takes the deck and
# nothing else (a second argument is a hard error, app/main.f90). Verbosity is
# "system": {"logger": {"level": "verbose"}} in the deck. Levels: debug,
# verbose, info (default), performance, warning, error, knowledge.

# Run linter
fortitude check

# Run only metalquicha tests (skip dependency tests)
ctest --test-dir build -R "mqc"
```

## Common Issues / Gotchas

1. **Coverage build type**: Must use `Coverage-mqc` (not `Coverage`) - the flags are only applied for this specific build type in `cmake/CMakeLists.txt`

2. **Dependency tests running**: mostly fixed -- `cmake/MqcDependencies.cmake` turns off every dependency test suite that can be turned off, which is all but two. `mctc-lib` and `dftd4` `add_subdirectory("test")` unconditionally and cannot be configured away, so their ~38 cases remain. `ctest --preset mqc`, or `ctest -R "mqc"`, still filters to this project's own.

3. **Linter `only` clause**: The `fortitude` linter requires all `use` statements to have an `only` clause - code will fail CI without this

4. **Units mismatch**: Coordinates are stored in Bohr internally but input files often use Angstrom - always use `to_bohr()` / `to_angstrom()` for conversions

5. **MPI validation tests**: When running `run_validation.py --mpi`, unfragmented tests still use `np=1` while fragmented tests use multiple ranks. To test parallel Hessian on unfragmented systems, mark the test as `"type": "fragmented"` since displacements are distributed across ranks

6. **Coverage on CI**: GitHub Actions uses 2 cores, so MPI validation uses `--np 2`. For local testing with more cores, run validation manually with higher `--np`
