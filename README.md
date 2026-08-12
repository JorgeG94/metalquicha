![License](https://img.shields.io/github/license/JorgeG94/metalquicha?style=for-the-badge)
![GitHub Repo stars](https://img.shields.io/github/stars/JorgeG94/metalquicha?style=for-the-badge&logo=github)
[![Issues](https://img.shields.io/github/issues/JorgeG94/metalquicha?style=for-the-badge)](https://github.com/JorgeG94/metalquicha/issues)
[![codecov](https://img.shields.io/codecov/c/github/JorgeG94/metalquicha?style=for-the-badge&logo=codecov)](https://codecov.io/gh/JorgeG94/metalquicha)
[![ReadTheDocs](https://img.shields.io/badge/docs-ReadTheDocs-8CA1AF?style=for-the-badge&logo=readthedocs&logoColor=white)](https://metalquicha.readthedocs.io/en/latest/)
[![FORD](https://img.shields.io/badge/docs-FORD-734F96?style=for-the-badge&logo=fortran&logoColor=white)](https://jorgeg94.github.io/metalquicha/)

# Met'al q'uicha (metalquicha)

<p align="center">
  <img src="images/sunflower.png" alt="Otter coding logo" title="Project logo" width="250">
</p>

Yes, this is AI generated (the image) if you know an artist, please let me know.



Met'al q'uicha (the Huastec (tenek) word for sunflower), which I'll just write as metalquicha, is a sample quantum chemistry backend
with focus on using the [pic](https://github.com/JorgeG94/pic) library and its derivatives:
[pic-mpi](https://github.com/JorgeG94/pic-mpi) and [pic-blas](https://github.com/JorgeG94/pic-blas)
which are Fortran based implementations of commonly used routines such as sorting algorithms,
array handling, strings, loggers, timers, etc.

The documentation is hosted at readthedocs, [here](https://metalquicha.readthedocs.io/en/latest/index.html).

Additionally, users can opt to try the [vapaa](https://github.com/jeffhammond/vapaa) backend for the `mpi_f08` module
to ensure cross compiler portability. Please report any issues associated here and in vapaa.

Metalquicha implements a naive backend for unfragmented and fragmented quantum chemistry
calculations. Three chemistry engines are available:

- [tblite](https://github.com/tblite/tblite) for semi-empirical xTB (GFN1, GFN2), on the CPU
- [libcint](https://github.com/sunqm/libcint) for Gaussian-basis ab initio on the
  CPU — Hartree-Fock and Kohn-Sham DFT, both restricted and unrestricted, plus
  MP2 and CCSD(T), each conventional or density-fitted. This one exists to be
  checked against as much as to be run: it gives the GPU path a second, independent implementation to
  disagree with, and every method in it is validated against PySCF.
- [NVIDIA cuEST](https://developer.nvidia.com/cuda/cuda-x-libraries/cuest) for
  Hartree-Fock and Kohn-Sham DFT on the GPU — energies, analytic gradients and
  Hessians, for whole molecules and for every fragment of an MBE/GMBE expansion.
  See **[CUEST.md](backends/cuest/CUEST.md)** for the full story: build instructions, the 20
  available functionals, and the validation numbers.

Both plug in behind the same `qc_method_t` interface, so fragmentation, screening
and many-body assembly are unchanged by the choice of engine.

If you are interested in contributing, please see [here](https://github.com/JorgeG94/pic/blob/main/contributing.md). Pic is the main project here and all the contributions fall downstream.

You can see [Project](https://github.com/users/JorgeG94/projects/4) for some information on development priorities and things being done!

## AI Disclaimer

The development of Metalquicha has been assisted by LLMs, such as ChatGPT, and Claude. The philosophy of "vibe coding" applied to this project is as follows:

- The programmer (Jorge), describes the overall architecture of a subroutine to be implemented and provides pseudocode
- The LLM produces an implementation that compiles
- The programmer writes a unit test for the function and validates the subroutine
- The LLM is asked to optimize the code while keeping the tests passing
- The programmer evaluates the code and evaluates if the routine needs to be redone or just upgraded by hand
- Either the programmer changes the code themselves or if they are lazy or cooking dinner while developing, they ask the LLM to try again

This was applied for routines such as the `mqc_finite_difference` module, which is pretty trivial to implement.

LLMs were also extensively used to add comments and basic documentation for the code. The idea is that
Metalquicha is a platform for development of fragmentation methods aimed to be suitable for everyone -
from students with no experience in Fortran and/or Quantum Chemistry to experienced researchers with
extensive expertise in both.

*Justification for LLM use*

I wanted to see to what extent LLMs can be used for Fortran code development. I can conclude that they are actually quite good.

## Can I use AI to study and work on this codebase?

Yes. But keep in mind that code reviews will still happen.

## Building

You will need an internet connection to download the dependencies. The main dependencies are:

- CMake
- A Fortran compiler
- An MPI installation
- A BLAS/LAPACK install
- TBLITE (will be downloaded automatically), for xTB
- NVIDIA cuEST and CUDA 12, for GPU Hartree-Fock/DFT (optional; see [CUEST.md](backends/cuest/CUEST.md))

`cmake -B build` with no options gives you tblite **and** the CPU ab initio path
— Hartree-Fock, DFT, MP2 and coupled cluster in a Gaussian basis, no GPU needed.
All three dependencies are fetched automatically, so **a default configure needs
network access.** On a machine without it, either point
`FETCHCONTENT_SOURCE_DIR_LIBCINT` and `FETCHCONTENT_SOURCE_DIR_LIBXC` at local
copies, or turn them off and build the xTB path alone.

| Option | Default | What it controls |
| --- | --- | --- |
| `-DMQC_ENABLE_TBLITE=` | `ON` | xTB (GFN1/GFN2) through tblite |
| `-DMQC_ENABLE_LIBCINT=` | `ON` | Gaussian integrals on the CPU, no GPU needed |
| `-DMQC_ENABLE_LIBXC=` | `ON` | Exchange-correlation functionals, so DFT |
| `-DMQC_ENABLE_HDF5=` | `OFF` | Binary checkpoints, to restart a gradient or Hessian |
| `-DMQC_ENABLE_CUEST=` | `OFF` | GPU Hartree-Fock/DFT; also needs `-DCUEST_ROOT=` |

Turning `LIBXC` off while leaving `LIBCINT` on is a supported build, but note
what it means: every deck naming a functional is refused, because there is
nothing to evaluate it with. That is the reason both default to on.

**tblite needs gfortran, ifort or ifx.** With nvfortran or LLVM Flang the
configure stops and says so; pass `-DMQC_ENABLE_TBLITE=OFF` and the rest of the
program — including the CPU ab initio path — builds normally.

The CPU backend exists mainly so results can be checked without a GPU; it is
validated against PySCF rather than tuned for speed.

You can then simply:

```
mkdir build
cd build
cmake ../
make -j
```

To build the GPU backend instead of (or alongside) xTB:

```
FC=mpifort CC=mpicc cmake -B build \
    -DMQC_ENABLE_TBLITE=OFF \
    -DMQC_ENABLE_CUEST=ON \
    -DCUEST_ROOT=/path/to/libcuest-linux-x86_64-<ver>_cuda12-archive
cmake --build build -j
```

Only the cuEST shared library is needed; the Fortran bindings are pre-generated
and vendored, and can optionally be fetched from
[mod_cuest](https://github.com/JorgeG94/mod_cuest) instead. Running needs a GPU
of compute capability 8.0 or newer. Full details in [CUEST.md](backends/cuest/CUEST.md).

### Notes on Fortran compiler compatibility

If you enable tblite (enabled by default at the moment) you are going to be blocked by which compilers does tblite
support. If you decide to not build tblite and just build the framework the code will work with most modern compilers.

Supported compilers:

Using TBlite: gcc, ifx, ifort

Without tblite, i.e. no quantum chemistry: gcc, nvfortran, flang(new), ifx, ifort

### Building with the Fortran Package Manager (FPM)

Before executing: The tblite package and some of its dependencies depend on `-lblas` which
is usually not installed, i.e. I use openblas or mkl. You will need to create a symlink
to `libblas.a`. You can do this by knowing where BLAS is installed and doing:

```
ln -s ${BLAS_ROOT}/lib/libopenblas.a ${LOCAL_BLAS_ROOT}/libblas.a
```

Then: `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LOCAL_BLAS_ROOT`

If you don't do this, then things will not work!

Simply then just do: `fpm install --prefix . --compiler mpifort --profile release`

#### Obtaining the FPM

Install the FPM following the [instructions](https://fpm.fortran-lang.org/install/index.html#install) and then simply: `fpm install`

## Running a calculation

Input is a JSON deck. Run it in serial:

```bash
./build/mqc validation/inputs/prism.json
```

or across ranks, which is how fragmented calculations are meant to be run:

```bash
mpirun -np 4 ./build/mqc validation/inputs/prism.json
```

A minimal deck:

```json
{
    "schema": { "name": "example", "version": "1.0" },
    "molecules": [
        { "xyz": "water.xyz", "molecular_charge": 0, "molecular_multiplicity": 1 }
    ],
    "model": { "method": "gfn2" },
    "driver": "Energy",
    "keywords": {
        "fragmentation": { "method": "MBE", "level": 2 }
    }
}
```

A molecule gives either `xyz` (a path, resolved relative to the deck) or
`symbols` plus a flat `geometry` list. Atom indices are 0-based. Bonds listed in
`connectivity` are marked broken automatically when their two atoms land in
different fragments -- that is derived, not declared.

`driver` is `Energy`, `Gradient` or `Hessian`.

For Hartree-Fock or DFT rather than xTB:

```json
"model": {
    "method": "hf",
    "basis": "def2-svp",
    "aux_basis": "def2-universal-jkfit",
    "functional": "pbe0"
}
```

`functional` applies to `dft` only. Which backend runs depends on the build:
cuEST when it is compiled in, otherwise libcint on the CPU. cuEST always
density-fits J and K, so `aux_basis` is required there. libcint has both paths
and uses exact integrals unless asked:

```json
"keywords": { "scf": { "density_fitting": true } }
```

The full keyword reference is in
[the documentation](https://metalquicha.readthedocs.io/en/latest/).

The `.mqc` text format and its `mqc_prep.py` generator were removed in 0.2.0.
See `mqc_docs/source/input_files.rst` for the migration table.

## Driving it from Python

The program can be driven from Python. Fortran still does the calculation and
still owns MPI; Python sets up the molecule, decides what to compute, and reads
the answers back.

Python runs on rank 0 and nowhere else. When it asks for a calculation, the
Fortran side spreads the work over the whole job and returns when it is done --
so a script that reads as single-threaded runs on as many nodes as it was
launched with, and `mpirun -np 64 python script.py` is a valid way to start one.

The interface loads `libmqc.so`, which is a separate target from the executable:

```bash
cmake -B build
cmake --build build --target mqc_shared
export PYTHONPATH=$PWD/python
```

The package looks for the library next to an in-tree `build/`; `MQC_LIBRARY`
overrides that with an explicit path.

```python
import mqc

with mqc.session():
    cluster = mqc.System.from_xyz("water20.xyz")
    cluster.auto_monomers()
    result = mqc.MBE(cluster, level=2, method="gfn2").run(label="w20")
    print(result.energy)
```

Everything happens inside `mqc.session()`, which starts MPI on entry and stops
it on exit.

Runnable examples are in `python/examples/`: `backends.py` covers standalone and
fragmented calculations for both xTB and Hartree-Fock and asserts against
reference energies, and `energy_screened_mbe.py` shows a two-pass calculation
that recomputes only the terms whose contribution exceeded a threshold.

Which methods are available depends on the build, though the default has all of
them: `gfn1`/`gfn2` need `MQC_ENABLE_TBLITE`, while `MQC_ENABLE_LIBCINT` brings
the CPU ab initio path — `hf`, `dft`, `mp2`, `ccsd`, `ccsd(t)`, and the `ri-`
spellings of the correlated ones. `dft` additionally needs `MQC_ENABLE_LIBXC`,
which is on by default too.

Gradients come from xTB and cuEST; the CPU backend refuses them rather than
returning something untested. Open shells run unrestricted for Hartree-Fock and
DFT, and are refused for MP2 and coupled cluster on the same principle — both
transforms want one set of orbitals, and an approximate answer is worse than
none.

Details, including density fitting and the current limitations, are in
[the Python interface documentation](https://metalquicha.readthedocs.io/en/latest/python_interface.html).
