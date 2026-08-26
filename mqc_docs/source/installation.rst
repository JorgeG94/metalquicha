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

Optional components
===================

Some backends are chosen at configure time. To see what a binary actually ended up
with, run ``mqc --version``, which prints a ``features:`` line.

.. list-table::
   :header-rows: 1
   :widths: 28 12 60

   * - Option
     - Default
     - What it adds
   * - ``MQC_ENABLE_TBLITE``
     - ``ON``
     - xTB (GFN1, GFN2) through tblite
   * - ``MQC_ENABLE_LIBCINT``
     - ``ON``
     - The CPU ab initio path: Hartree-Fock, DFT, MP2, coupled cluster
   * - ``MQC_ENABLE_DLFIND``
     - ``OFF``
     - Geometry optimization, through DL-FIND
   * - ``MQC_ENABLE_CUEST``
     - ``OFF``
     - GPU integrals (needs ``CUEST_ROOT`` and an sm_80 card)
   * - ``MQC_ENABLE_HDF5``
     - ``OFF``
     - Binary checkpoints carrying gradients and Hessians
   * - ``MQC_ENABLE_CREST``
     - ``OFF``
     - Conformer and ensemble sampling, through CREST

``MQC_ENABLE_LIBCINT`` is named after the library the backend was first written
against, not the one it links now. A default build takes its integrals from
`libfint <https://github.com/JorgeG94/libfint>`_, an all-Fortran port of libcint
that this project also maintains; the flag switches the whole CPU ab initio path
on and off regardless of which of the two is underneath.

Two options are deliberately absent from the table above. Both default to ``ON``,
both are marked advanced in CMake, and neither is a knob to reach for:

``MQC_USE_LIBFINT``
   ``OFF`` takes libcint itself instead of the Fortran port. The two are
   bit-identical, so this changes what gets built and not what gets computed.
   libcint is being phased out of the default surface; the option stays because
   the two configurations compile different source -- libfint carries L (sp)
   shells and libcint cannot represent one -- and CI builds both.

``MQC_ENABLE_LIBXC``
   ``OFF`` removes the exchange-correlation functionals, so every deck naming
   one is refused. It survives for the coverage build, which disables libxc
   because libxc refuses that build type, rather than because a metalquicha
   without DFT is something to want.

``MQC_ENABLE_DLFIND`` is off for a licensing reason rather than a technical one:
DL-FIND is LGPL-3 and metalquicha is MIT. It is fetched and linked as a shared
library so the two stay separable, and turning it on is a choice the person
building makes. See :doc:`geometry_optimization`.

``MQC_ENABLE_CREST`` is off for the same reason -- CREST is LGPL-3.0 -- and for
a second one: it needs ``WITH_GFN0=ON`` to be useful, and CREST's own tblite pin
disagrees with the release this project builds against. See
:doc:`conformer_sampling`.

.. code-block:: bash

   cmake -B build -DMQC_ENABLE_DLFIND=ON

Each of the fetched dependencies can be pointed at a local clone instead of being
downloaded, which is how to configure on a cluster node with no network:

.. code-block:: bash

   cmake -B build -DMQC_ENABLE_DLFIND=ON \
         -DFETCHCONTENT_SOURCE_DIR_LIBDLFIND=/path/to/libdlfind \
         -DFETCHCONTENT_SOURCE_DIR_LIBFINT=/path/to/libfint

The variable is named after what is actually fetched, so a default build wants
``_LIBFINT``. ``_LIBCINT`` is the right one only alongside
``-DMQC_USE_LIBFINT=OFF``; setting it on a default configure does nothing and
the fetch still reaches for the network, which on a compute node with none is a
configure that fails rather than one that falls back.

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
