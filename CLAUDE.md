# Metalquicha - Claude Code Context

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
# Build
mkdir build && cd build
cmake ..
make -j

# Run (serial)
./mqc input.mqc

# Run (parallel)
mpirun -np 4 ./mqc input.mqc

# Run tests
ctest
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
| `src/io/mqc_config_parser.f90` | Parses `.mqc` input files |
| `src/fragmentation/mqc_mbe.f90` | Core MBE implementation (62KB) |
| `src/fragmentation/mqc_physical_fragment.f90` | Fragment representation, H-capping |
| `src/fragmentation/mqc_gmbe_utils.f90` | GMBE with PIE (overlapping fragments) |
| `src/methods/mqc_method_xtb.f90` | xTB via tblite library |
| `src/core/mqc_result_types.f90` | Result containers (energy, gradient, hessian) |
| `src/vibrational/mqc_vibrational_analysis.f90` | Frequency calculations |
| `src/vibrational/mqc_thermochemistry.f90` | Thermochemistry (RRHO) |

## Input Format (.mqc)

Section-based format parsed by `mqc_config_parser.f90`:

```
%schema
name = example
version = 1.0
index_base = 1
end

%model
method = gfn2
end

%driver
type = Energy
end

%structure
charge = 0
multiplicity = 1
end

%geometry
O  0.000  0.000  0.117
H  0.000  0.757 -0.469
H  0.000 -0.757 -0.469
end

%fragmentation
method = MBE
level = 2
cutoff_method = distance
end

%cutoffs
dimer = 5.0
trimer = 4.0
end
```

Use `python mqc_prep.py input.json` to convert JSON to `.mqc` format.

## Core Concepts

### Many-Body Expansion (MBE)
- Breaks system into fragments (monomers, dimers, trimers...)
- E_MBE = ΣE_I - ΣΔE_IJ + ΣΔE_IJK - ...
- Truncate at level N for O(N^k) scaling

### Generalized MBE (GMBE)
- Handles overlapping fragments
- Uses Principle of Inclusion-Exclusion (PIE)
- Automatically computes intersection terms

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

- `GFN1` - GFN1-xTB (faster, older parametrization)
- `GFN2` - GFN2-xTB (recommended, more accurate)

## Solvation Models

- `ALPB` - Analytical Linearized Poisson-Boltzmann
- `GBSA` - Generalized Born with Solvent-Accessible Surface Area
- `CPCM` - Conductor-like Polarizable Continuum Model

Supports 40+ solvents (water, methanol, acetonitrile, etc.)

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

# Run validation suite
python run_validation.py
```

Test files follow pattern `test/test_mqc_*.f90`.

## Coding Conventions

See `FORTRAN_STYLE.md` for the complete style guide. Key points:

- Module naming: `mqc_*.f90`
- Type naming: `*_t` suffix (e.g., `calculation_result_t`)
- Internal units: Bohr (length), Hartree (energy)
- Use `to_bohr()` / `to_angstrom()` for conversions
- Error handling via `error_t` type with `has_error()` / `get_message()`
- **Linter**: Uses `fortitude check` - all `use` statements must have `only` clause
- **CI/CD**: GitHub Actions workflows in `.github/workflows/` (build, test, coverage via Codecov)
- **Environment**: Conda environment defined in `environment.yml`

**Forbidden**: `GOTO`, arithmetic IF, implicit SAVE, COMMON blocks, EQUIVALENCE, fixed-form source, assumed-size arrays

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
1. Add to appropriate section in `mqc_config_parser.f90`
2. Add field to relevant config type
3. Handle in `mqc_config_adapter.f90`

### Add a new test
1. Create `test/test_mqc_feature.f90`
2. Use test-drive framework
3. Add to CMakeLists.txt via `ADDTEST(mqc_feature)`

## Validation Tests

Located in `validation/` directory:
- `h3o.mqc` - Unfragmented hydronium
- `prism.mqc` - Water prism MBE(2)
- `overlapping_gly3.mqc` - Glycine tripeptide GMBE(1)
- `w20_isomer.mqc` - Water 20-mer MBE(2)

## Output Format

JSON output (`output_<basename>.json`) contains:
- Total energy and per-fragment breakdown
- Gradients (if requested)
- Hessian (if requested)
- Vibrational frequencies and thermochemistry
- GMBE coefficients and intersection terms

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

# Run with verbose output
./build/mqc input.mqc --verbose

# Run linter
fortitude check

# Run only metalquicha tests (skip dependency tests)
ctest --test-dir build -R "mqc"
```

## Common Issues / Gotchas

1. **Coverage build type**: Must use `Coverage-mqc` (not `Coverage`) - the flags are only applied for this specific build type in `cmake/CMakeLists.txt`

2. **Dependency tests running**: Use `ctest -R "mqc"` to filter to only metalquicha tests; otherwise ctest runs tests from all dependencies (tblite, pic, etc.)

3. **Linter `only` clause**: The `fortitude` linter requires all `use` statements to have an `only` clause - code will fail CI without this

4. **Units mismatch**: Coordinates are stored in Bohr internally but input files often use Angstrom - always use `to_bohr()` / `to_angstrom()` for conversions

5. **MPI validation tests**: When running `run_validation.py --mpi`, unfragmented tests still use `np=1` while fragmented tests use multiple ranks. To test parallel Hessian on unfragmented systems, mark the test as `"type": "fragmented"` since displacements are distributed across ranks

6. **Coverage on CI**: GitHub Actions uses 2 cores, so MPI validation uses `--np 2`. For local testing with more cores, run validation manually with higher `--np`
