.. _installation:

------------
Installation
------------

.. contents::


Obtaining Metalquicha
======================

Currently, the only way to obtain Metalquicha is to clone the git repository from the
`project web site <https://github.com/JorgeG94/metalquicha>`_. Simply:

.. code-block:: bash

   git clone git@github.com:JorgeG94/metalquicha.git

Building Metalquicha
=======================

As of now, Metalquicha uses tblite as its quanutm chemistry engine thus, it is the
most complex dependency to account for. Basically, you'll need one of the supported
Fortran compilers for tblite. I recommend using gfortran, ifort, of ifx. Flang, NVFortran,
and other compilers are not quite supported by tblite as of now.

You can build without tblite, but as of now there is no other quantum chemistry engine so
you'll just have an MPI distributed program. So for now, just build with tblite.

Metalquicha uses CMake as its build system, the minimum CMake version is 3.22. To build Metalquicha,
you will need the following dependencies:

- A Fortran compiler (supported by tblite)
- CMake
- An MPI implementation with support for the `mpi_f08` module
- BLAS and LAPACK libraries

To build Metalquicha, you can use the following commands:

.. code-block:: bash

   mkdir build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=$PATH_TO_YOUR_INSTALL ..
   make -j install

An example path to your INSTALL can be `$HOME/install/metalquicha`.

Building with the Fortran Package Manager (fpm)
===============================================

First, you can install the FPM by following the instructions at https://fpm.fortran-lang.org/install/index.html#install

Once you have fpm installed, you can build Metalquicha by running the following command:

.. code-block:: bash

   fpm install --prefix=$PATH_TO_YOUR_INSTALL --compiler mpifort

Running Metalquicha
=======================

Metalquicha has a two input file (ish) system. The first file is a JSON format that is compliant
with the QCSchema specificiation. These json files get converted to the mqc format which will
be then read by metalquicha as the main input file. The idea is that the user will only
interact with the JSON files and the mqc files will be generated automatically. However,
the mqc files are pretty human readable, so you can also edit them directly if you'd like.

In order to go from JSON to mqc, you need to use thte `mqc_prep.py` program. Simply run:

.. code-block:: bash

   mqc_prep.py input.json

You can see examples of the JSON input files in the `sample_json` directory of the
Metalquicha repository.

Once you have the mqc file, you can run Metalquicha simply by running:

.. code-block:: bash

   ./mqc input.mqc

This will run Metalquicha in serial mode (i.e. 1 MPI process). TBLite will use OpenMP
threads to enable parallelism within the single MPI process. To run Metalquicha in parallel
you can use `mpirun`, or `mpiexec`. For example, to run with 4 MPI processes:

.. code-block:: bash

   mpirun -n 4 ./mqc input.mqc

Metalquicha will run on multiple nodes, you can do that by using:

.. code-block:: bash

   mpirun -np 64 --map-by ppr:32:node ./mqc input.mqc

This will run Metalquicha with 64 MPI processes, with 32 processes per node.
