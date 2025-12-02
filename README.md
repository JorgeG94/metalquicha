# Met'al q'uicha (metalquicha)

<p align="center">
  <img src="images/sunflower.png" alt="Otter coding logo" title="Project logo" width="250">
</p>
Yes, this is AI generated (the image) if you know an artist, please let me know.


Met'al q'uicha (the Huastec (tenek) word for sunflower), which I'll just write as metalquicha is a sample quantum chemistry backend
with focus on using the [pic](https://github.com/JorgeG94/pic) library and its derivatives:
[pic-mpi](https://github.com/JorgeG94/pic-mpi) and [pic-blas](https://github.com/JorgeG94/pic-blas)
which are Fortran based implementations of commonly used routines such as sorting algorithms,
array handling, strings, loggers, timers, etc. While pic-mpi and pic-blas provide modern
Fortran interfaces to MPI and BLAS implementations in a portable way. Specifically, the MPI
library lets the user switch betwen the `mpi` and `mpi_f08` modules with ease.

Metalquicha implements a naive backend for unfragmented and fragmented quantum chemistry
calculations. Currently, metalquicha uses [tblite](https://github.com/tblite/tblite) as
its chemistry engine which performs energy calculations.

The main purpose of this package is to showcase the ability of being able to write
simple, powerful Fortran based programs that are able to access massively parallel
ecosystems with ease.

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

### Building with the Fortran Package Manager (FPM)

*FPM will only work if you are building with openblas, since the linking step is hardcoded.*

Simply then just do: `fpm install --prefix . --compiler mpifort --profile release`

#### Obtaining the FPM

Install the FPM following the [instructions](https://fpm.fortran-lang.org/install/index.html#install) and then simply: `fpm install`
