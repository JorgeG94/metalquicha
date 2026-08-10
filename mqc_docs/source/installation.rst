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

A simple way to install the dependencies is using conda, simply do:

.. code-block:: bash

   conda env create --name mqc --file environment.yml
   conda install -c conda-forge -n mqc gfortran=14.2.0
   conda activate mqc

Then, to build Metalquicha, you can use the following commands:

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

Metalquicha takes a single JSON input file, in a format compliant with the QCSchema
specification.

You can see examples of the JSON input files in the `sample_json` directory of the
Metalquicha repository.

.. note::

   **Changed in 0.2.0.** Earlier versions converted the JSON to an intermediate
   ``.mqc`` file with a ``mqc_prep.py`` helper, and read that. Both are gone;
   pass the ``.json`` file directly. See :ref:`migrating_from_mqc`.

Run Metalquicha by passing it your input:

.. code-block:: bash

   ./mqc input.json

This will run Metalquicha in serial mode (i.e. 1 MPI process). TBLite will use OpenMP
threads to enable parallelism within the single MPI process. To run Metalquicha in parallel
you can use `mpirun`, or `mpiexec`. For example, to run with 4 MPI processes:

.. code-block:: bash

   mpirun -n 4 ./mqc input.json

Metalquicha will run on multiple nodes, you can do that by using:

.. code-block:: bash

   mpirun -np 64 --map-by ppr:32:node ./mqc input.json

This will run Metalquicha with 64 MPI processes, with 32 processes per node.
