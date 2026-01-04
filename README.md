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

## Running a calculation

To run a calculation you need to process the JSON input into our `mqc` format. To do this, you can simply do:

```bash
python mqc_prep.py validation/inputs/prism.json
```

And this will generate a `prism.mqc`. Which can be simply run as `./build/mqc validation/inputs/prism.mqc` to be run
in serial mode. Or `mpirun -np 4 ./build/mqc validation/inputs/prism.mqc`.

A sample `mqc` file is shown below:

```
%schema
name = mqc-frag
version = 1.0
index_base = 0
units = angstrom
end  ! schema

%model
method = XTB-GFN1
basis = cc-pVDZ
aux_basis = cc-pVDZ-RIFIT
end  ! model

%driver
type = Energy
end  ! driver

%structure
charge = 0
multiplicity = 1
end  ! structure

%geometry
3

O 0 0 0.119262
H 0 0.763239 -0.477047
H 0 -0.763239 -0.477047
end  ! geometry

%scf
maxiter = 300
tolerance = 1e-06
end  ! scf
```

If you don't want to use the python script, you can modify this file by adding an xyz formatted geometry. Supported calculations are `Energy`, `Gradient`, and `Hessian`.
