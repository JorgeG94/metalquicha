# The build system, the serial build and compilers

`CMAKE_STYLE.md` is the style guide; this is what is not obvious from it. Read from `AGENTS.md`, which says when.

## Build System

`CMAKE_STYLE.md` is the guide. The short version: the top-level `CMakeLists.txt`
is an outline of four `include()`s, every dependency is declared in its own
`cmake/modules/Find<pkg>.cmake` through `mqc_fetch`, every option lives in
`cmake/MqcOptions.cmake`, and validation programs go through
`mqc_add_validation_program`. `CMakePresets.json` carries the configurations
worth remembering, `perlmutter` among them.

### Building serially

`-DMQC_ENABLE_SERIAL=ON` turns off OpenMP, MPI and tblite together, leaving
BLAS/LAPACK as the only parallel-adjacent dependency. The three go together
because they are not independent: tblite is built with OpenMP and cannot be
configured without it, and a build with neither has no reason to carry MPI.

It exists for a compiler that cannot yet build the parallel forms, not for
speed -- nothing is threaded, nothing is distributed, and the fragment loop
runs one fragment at a time. Do not benchmark against it.

Two things it does that are not obvious. It sets
`CMAKE_DISABLE_FIND_PACKAGE_OpenMP`, because libfint runs its own
`find_package(OpenMP)` to make its shell-quartet workspace `!$omp
threadprivate` and would otherwise put the flag back. And it compiles
`compat/omp_lib_serial.f90`, which answers the OpenMP runtime API as one
thread: the `!$omp` directives vanish without the flag but the twenty-odd
unguarded `use omp_lib` statements do not, and gfortran leaves
`omp_get_max_threads_` in libgomp. Linking libgomp instead would work and would
be one line, but it restores the dependency the mode exists to remove and
reports the machine's core count in a build with no parallel region to use it.

The whole CPU validation suite -- 307 cases -- passes on a serial gfortran
build.

**LFortran** is why this exists. As of 0.66.0 every source file in this project
compiles and `libmetalquicha.a` links; the `mqc` executable does not, and
`.github/workflows/lfortran-build.yml` therefore builds the library target and
stops there. It is a pinned 0.66.0, deliberately -- LFortran moves fast enough
that a floating version turns the job red for reasons unrelated to the commit.

What stops the executable is two LFortran codegen bugs, both reached only when a
*program* forces a module's procedures to be instantiated, which is why the
library goes through and the binary does not:
`visit_StructConstructor() not implemented`, on a component default initialised
by a structure constructor -- pic-mpi's serial `comm_t`, `request_t` and `win_t`
-- and `Internal error: x is nullptr`, reached from `run_calculation`. Neither
is anything in this tree. Reproducers are kept outside this repository.

The widest-reaching gap is still that a continued `!$omp` directive is rejected
-- 359 of them across 23 files -- which is what makes serial the only shape
worth trying. LFortran also needs `--implicit-interface`,
`--mangle-underscore-external` and `--legacy-array-sections`, set in
`cmake/CMakeLists.txt`; each is load-bearing and each is commented there. The
second of those mangles `bind(c)` binding labels as well as Fortran externals,
which is what stops libxc's own Fortran examples from linking; the library
target does not build them.

## Compiler Support

| Compiler | Status |
|----------|--------|
| gfortran 11+ | Full support |
| Intel ifort | Full support |
| Intel ifx | Full support |
| nvfortran | Partial (no tblite) |
| LLVM Flang | Partial (no tblite) |
| LFortran 0.66 | Library only, serial build; see "Building serially" |

## Dependencies

| Library | Purpose |
|---------|---------|
| tblite | xTB quantum chemistry engine |
| pic | Utility library (sorting, logging, timers) |
| pic-mpi | MPI wrapper layer |
| pic-blas | BLAS/LAPACK interface |
| test-drive | Unit testing framework |
| jsonfortran | I/O |

## Dependency tests in ctest

`cmake/MqcDependencies.cmake` turns off every dependency test suite that can be turned off, which is all but two. `mctc-lib` and `dftd4` `add_subdirectory("test")` unconditionally and cannot be configured away, so their ~38 cases remain. `ctest --preset mqc`, or `ctest -R "mqc"`, filters to this project's own.

## Coverage

The build type is `Coverage-mqc`, not `Coverage`: the flags are applied for that type only, in `cmake/CMakeLists.txt`.

```bash
cmake -DCMAKE_BUILD_TYPE=Coverage-mqc -B build && cmake --build build && cmake --build build --target coverage   # needs lcov, GCC only
```
