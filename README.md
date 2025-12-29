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
calculations. Currently, metalquicha uses [tblite](https://github.com/tblite/tblite) as
its chemistry engine which performs energy calculations.

If you are interested in contributing, please see [here](https://github.com/JorgeG94/pic/blob/main/contributing.md). Pic is the main project here and all the contributions fall downstream.

## Building

You will need an internet connection to download the dependencies. The main dependencies are:

- CMake
- A Fortran compiler
- An MPI installation
- A BLAS/LAPACK install
- TBLITE (will be downloaded automatically)

You can then simply:

```
mkdir build
cd build
cmake ../
make -j
```
### Notes on Fortran compiler compatibility

If you enable tblite (enabled by default at the moment) you are going to be blocked by which compilers does tblite
support. If you decide to not build tblite and just build the framework the code will work with most modern compilers.

### Building with the Fortran Package Manager (FPM)

*FPM will only work if you are building with openblas, since the linking step is hardcoded.*

Simply then just do: `fpm install --prefix . --compiler mpifort --profile release`

#### Obtaining the FPM

Install the FPM following the [instructions](https://fpm.fortran-lang.org/install/index.html#install) and then simply: `fpm install`
