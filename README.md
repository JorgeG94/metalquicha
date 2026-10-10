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
- [libfint](https://github.com/JorgeG94/libfint) for Gaussian-basis ab initio on
  the CPU — an all-Fortran port of libcint, and what a default build uses;
  `-DMQC_USE_LIBFINT=OFF` takes libcint itself instead — Hartree-Fock and Kohn-Sham DFT, both restricted and unrestricted, plus
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
| `-DMQC_ENABLE_CZT=` | `ON` | Gaussian integrals on the CPU, no GPU needed |
| `-DMQC_ENABLE_LIBXC=` | `ON` | Exchange-correlation functionals, so DFT |
| `-DMQC_ENABLE_HDF5=` | `OFF` | Binary checkpoints, to restart a gradient or Hessian |
| `-DMQC_ENABLE_CUEST=` | `OFF` | GPU Hartree-Fock/DFT; also needs `-DCUEST_ROOT=` |
| `-DMQC_ENABLE_MPI=` | `ON` | `OFF` builds against pic-mpi's single-rank backend, so no MPI install is needed |
| `-DMQC_USE_LIBFINT=` | `ON` | The all-Fortran port of libcint. `OFF` takes libcint itself, which needs C |
| `-DMQC_ENABLE_DLFIND=` | `OFF` | DL-FIND for geometry optimisation (LGPL-3, hence off) |

### The smallest useful build

```bash
cmake -DMQC_ENABLE_MPI=OFF -B build
cmake --build build -j
```

No MPI installation and no C compiler, and everything else is on by default:
tblite, libfint and libxc. That is xTB, Hartree-Fock, DFT, MP2 and coupled
cluster, with energies, gradients and Hessians — a complete shared-memory
quantum chemistry program that fragments as well, just without spreading the
fragments over ranks. Threads still work; it is one process.

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

## Measuring your machine

```bash
python3 benchmarks/run_benchmarks.py --exe build/mqc --record
```

Times this build on this machine and records a baseline; run it again after a
change and it reports what moved, judged against the spread it measured rather
than a fixed percentage. It ends by telling you what to run things with — the
thread count past which your machine stops improving, whether fragment work
wants threads or ranks here, and which of MP2 and RI-MP2 is actually faster on
your hardware. About twelve minutes, or three with `--quick`. See
[benchmarks/README.md](benchmarks/README.md).

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
them: `gfn1`/`gfn2` need `MQC_ENABLE_TBLITE`, while `MQC_ENABLE_CZT` brings
the CPU ab initio path — `hf`, `dft`, `mp2`, `ccsd`, `ccsd(t)`, and the `ri-`
spellings of the correlated ones. `dft` additionally needs `MQC_ENABLE_LIBXC`,
which is on by default too.

Gradients come from all three engines — xTB, the CPU path and cuEST — and are
validated case by case rather than assumed; MP2 gradients want
`keywords.correlation.freeze_core: false`, which the program will tell you.
Open shells run unrestricted for Hartree-Fock and DFT, and are refused for MP2
and coupled cluster on principle: both transforms want one set of orbitals, and
an approximate answer is worse than none.

Details, including density fitting and the current limitations, are in
[the Python interface documentation](https://metalquicha.readthedocs.io/en/latest/python_interface.html).

## Citing the software metalquicha builds on

Metalquicha is a layer over other people's work. If you publish results from it,
please cite the libraries your calculation used. Which ones that is depends on the
build and the method. A Kohn-Sham run also logs the papers for its functional,
read from libxc's own reference table.

### Electronic structure

- **libcint** (the integrals, `-DMQC_USE_LIBFINT=OFF`) and **libfint** (the
  default, an all-Fortran port of libcint): Q. Sun, *J. Comput. Chem.* **36**, 1664
  (2015), [doi:10.1002/jcc.23981](https://doi.org/10.1002/jcc.23981). Sources:
  [libcint](https://github.com/JorgeG94/libcint),
  [libfint](https://github.com/JorgeG94/libfint).
- **libxc** (exchange-correlation functionals): S. Lehtola, C. Steigemann,
  M. J. T. Oliveira and M. A. L. Marques, *SoftwareX* **7**, 1 (2018),
  [doi:10.1016/j.softx.2017.11.002](https://doi.org/10.1016/j.softx.2017.11.002).
  Cite the functional's own papers as well.
- **tblite** (GFN1-xTB and GFN2-xTB; [source](https://github.com/tblite/tblite)):
  - GFN2-xTB: C. Bannwarth, S. Ehlert and S. Grimme, *J. Chem. Theory Comput.*
    **15**, 1652 (2019),
    [doi:10.1021/acs.jctc.8b01176](https://doi.org/10.1021/acs.jctc.8b01176).
  - GFN1-xTB: S. Grimme, C. Bannwarth and P. Shushkov, *J. Chem. Theory Comput.*
    **13**, 1989 (2017),
    [doi:10.1021/acs.jctc.7b00118](https://doi.org/10.1021/acs.jctc.7b00118).
  - The xTB family: C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen,
    P. Pracht, J. Seibert, S. Spicher and S. Grimme, *WIREs Comput. Mol. Sci.*
    **11**, e1493 (2021), [doi:10.1002/wcms.1493](https://doi.org/10.1002/wcms.1493).
- **simple-dftd3** (D3 dispersion; [source](https://github.com/dftd3/simple-dftd3)):
  - D3: S. Grimme, J. Antony, S. Ehrlich and H. Krieg, *J. Chem. Phys.* **132**,
    154104 (2010), [doi:10.1063/1.3382344](https://doi.org/10.1063/1.3382344).
  - Becke-Johnson damping: S. Grimme, S. Ehrlich and L. Goerigk,
    *J. Comput. Chem.* **32**, 1456 (2011),
    [doi:10.1002/jcc.21759](https://doi.org/10.1002/jcc.21759).
- **dftd4** (D4 dispersion; [source](https://github.com/dftd4/dftd4)):
  - E. Caldeweyher, C. Bannwarth and S. Grimme, *J. Chem. Phys.* **147**, 034112
    (2017), [doi:10.1063/1.4993215](https://doi.org/10.1063/1.4993215).
  - E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth
    and S. Grimme, *J. Chem. Phys.* **150**, 154122 (2019),
    [doi:10.1063/1.5090222](https://doi.org/10.1063/1.5090222).
  - E. Caldeweyher, J.-M. Mewes, S. Ehlert and S. Grimme, *Phys. Chem. Chem.
    Phys.* **22**, 8499 (2020),
    [doi:10.1039/D0CP00502A](https://doi.org/10.1039/D0CP00502A).
- **cuEST** (the GPU backend): NVIDIA's library, linked from a local install.

### Geometry and conformers

- **DL-FIND** (geometry optimisation, through
  [libdlfind](https://github.com/JorgeG94/libdlfind)): J. Kästner, J. M. Carr,
  T. W. Keal, W. Thiel, A. Wander and P. Sherwood, *J. Phys. Chem. A* **113**,
  11856 (2009), [doi:10.1021/jp9028968](https://doi.org/10.1021/jp9028968).
- **CREST** (conformer search; [source](https://github.com/JorgeG94/crest)):
  - P. Pracht, F. Bohle and S. Grimme, *Phys. Chem. Chem. Phys.* **22**, 7169
    (2020), [doi:10.1039/C9CP06869D](https://doi.org/10.1039/C9CP06869D).
  - P. Pracht *et al.*, *J. Chem. Phys.* **160**, 114110 (2024),
    [doi:10.1063/5.0197592](https://doi.org/10.1063/5.0197592).

### Data

- **Basis Set Exchange** (every basis set extracted into `basis_sets/`, and the
  6-311++G(3df,2p) under `basis_sets/pople/`, which is generated from one of them):
  - B. P. Pritchard, D. Altarawy, B. Didier, T. D. Gibson and T. L. Windus,
    *J. Chem. Inf. Model.* **59**, 4814 (2019),
    [doi:10.1021/acs.jcim.9b00725](https://doi.org/10.1021/acs.jcim.9b00725).
  - D. Feller, *J. Comput. Chem.* **17**, 1571 (1996).
  - K. L. Schuchardt, B. T. Didier, T. Elsethagen, L. Sun, V. Gurumoorthi,
    J. Chase, J. Li and T. L. Windus, *J. Chem. Inf. Model.* **47**, 1045 (2007),
    [doi:10.1021/ci600510j](https://doi.org/10.1021/ci600510j).
- **AAMBS**, the minimal basis the QUAO analysis projects onto, transcribed from
  GAMESS: G. M. J. Barca *et al.*, *J. Chem. Phys.* **152**, 154102 (2020),
  [doi:10.1063/5.0005188](https://doi.org/10.1063/5.0005188). Its exponents' own
  sources are listed in `basis_sets/aambs/PROVENANCE.md`.
- **MINAO**, the minimal basis of the `minao` initial guess, taken from PySCF's
  copy of ANO-RCC: B. O. Roos, R. Lindh, P.-Å. Malmqvist, V. Veryazov and
  P.-O. Widmark, *J. Phys. Chem. A* **108**, 2851 (2004) and **109**, 6575 (2005).
  See `tools/minao/gen_minao_basis.py`.
- **Nuclear basis sets** for NEO (PB4-D to PB6-H): Q. Yu, F. Pavošević and
  S. Hammes-Schiffer, *J. Chem. Phys.* **152**, 244123 (2020). See
  `basis_sets/neo/PROVENANCE.md`.

### Validation

The reference numbers in `validation/` come from **PySCF**:
Q. Sun *et al.*, *J. Chem. Phys.* **153**, 024109 (2020),
[doi:10.1063/5.0006074](https://doi.org/10.1063/5.0006074), and Q. Sun *et al.*,
*WIREs Comput. Mol. Sci.* **8**, e1340 (2018),
[doi:10.1002/wcms.1340](https://doi.org/10.1002/wcms.1340). The NEO references come
from [Yang Yang's PySCF fork](https://github.com/theorychemyang/pyscf).

### Infrastructure

- **LAPACK** and a BLAS: E. Anderson *et al.*, *LAPACK Users' Guide*, 3rd ed.
  (SIAM, 1999), [doi:10.1137/1.9780898719604](https://doi.org/10.1137/1.9780898719604).
- **HDF5** (derivative checkpoints, optional): The HDF Group,
  [Hierarchical Data Format, version 5](https://www.hdfgroup.org/solutions/hdf5/).
- [json-fortran](https://github.com/jacobwilliams/json-fortran) (input and
  output), [test-drive](https://github.com/fortran-lang/test-drive) (the unit
  tests), and [pic](https://github.com/JorgeG94/pic),
  [pic-mpi](https://github.com/JorgeG94/pic-mpi) and
  [pic-blas](https://github.com/JorgeG94/pic-blas) (utilities, MPI and BLAS
  layers).
