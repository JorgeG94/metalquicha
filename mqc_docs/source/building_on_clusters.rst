.. _building_on_clusters:

=====================
Building on a cluster
=====================

.. contents::

:doc:`installation` covers a laptop, where conda supplies the compiler and the
libraries and CMake finds everything by itself. A cluster is not that, and the
differences are not cosmetic: on more than one of them the default configure
produces a binary that builds, links, runs, and returns wrong numbers above one
thread.

This page collects what changes, and then a section per machine.

What a cluster changes
======================

**The compiler is a wrapper.** ``ftn``, ``mpifort``, ``mpiifx`` -- a script that
calls the real compiler with the site's MPI, BLAS and accelerator flags already
on the line. Point ``CMAKE_Fortran_COMPILER`` at the wrapper, never at the
compiler behind it, or the flags it was going to add are silently missing.

**Compute nodes usually have no network.** Every fetched dependency can be
pointed at a local clone instead, which is how to configure where FetchContent
cannot reach GitHub:

.. code-block:: bash

   cmake -B build -DFETCHCONTENT_SOURCE_DIR_TBLITE=/path/to/tblite \
                  -DFETCHCONTENT_SOURCE_DIR_LIBFINT=/path/to/libfint

Configuring on a login node and building on a compute node is usually simpler.

**There are often two BLAS libraries, and only one is thread-safe.** Vendor math
libraries ship a serial build and a threaded build. The serial one is not merely
"the one that does not thread" -- on Cray it is not re-entrant, and two OpenMP
threads calling into it at once corrupt each other's workspace. It does not
crash and it does not return an error: it returns different numbers. See the
Perlmutter section below for what that looks like from the outside, because it
looks like a chemistry problem rather than a linking one.

**Threads and ranks divide the node between them.** ``mqc`` is MPI over
fragments and OpenMP inside a calculation. A batch script that gives all cores
to the ranks and leaves ``OMP_NUM_THREADS`` unset, or the reverse, is the usual
reason a run is slower than the same work on a workstation.

.. _perlmutter:

Perlmutter (NERSC)
==================

Two builds are worth having and they are not the same build. The CPU one is
GNU, and gets everything: xTB, the ``cenzontle`` ab initio path, DFT, CREST.
The GPU one is nvfortran, gets cuEST, and gives up tblite to have it.

CPU build (PrgEnv-gnu)
----------------------

``PrgEnv-gnu`` is the default programming environment at NERSC, so on a fresh
login there is usually nothing to load at all:

.. code-block:: bash

   module load PrgEnv-gnu          # already loaded on a default login
   cmake --preset default
   cmake --build --preset default -j 16
   ctest --preset mqc

``ftn`` is the Fortran compiler and CMake reports it as ``GNU``, which is
correct -- it is a wrapper around ``gfortran`` with Cray MPICH and libsci
already attached. There is no need to name an MPI: the wrapper is the MPI.

To add CREST for :doc:`conformer_sampling`, which is off by default:

.. code-block:: bash

   cmake -B build -DMQC_ENABLE_CREST=ON -DMQC_ENABLE_TBLITE=ON \
                  -DWITH_TBLITE=ON -DWITH_GFN0=ON

GPU build (cuEST)
-----------------

There is a preset, and it exists because each of its settings works around
something that fails elsewhere with a message naming neither the machine nor
the flag that fixes it:

.. code-block:: bash

   module load PrgEnv-nvidia cudatoolkit craype-accel-nvidia80
   export CUEST_ROOT=/path/to/libcuest-linux-x86_64-<version>-archive
   cmake --preset perlmutter
   cmake --build --preset perlmutter -j 16

``MQC_ENABLE_TBLITE`` is ``OFF`` in that preset and cannot simply be turned back
on: toml-f, which is reached through tblite, has a backslash string literal
nvfortran rejects and no flag makes it accept. ``CMAKE_STYLE.md`` has the rest
of what the preset is working around.

Things that bite here
---------------------

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Symptom
     - Cause and fix
   * - SCFs that will not converge, or ``dsygvd`` returning ``info`` in the
       forties, only above one OpenMP thread
     - The **serial** libsci got linked ahead of the threaded one. Handled by
       the build now -- see below -- but do not undo it by pinning
       ``BLAS_LIBRARIES`` at ``libsci_gnu.so``
   * - Configure fails at *"Detecting C compiler ABI info"* with
       ``undefined reference to cudaMalloc@libcudart.so.13``
     - ``cray-mpich``'s CUDA GTL and ``cudatoolkit`` disagree about the CUDA
       major version, and ``craype-accel-nvidia80`` forces ``-lmpi_gtl_cuda``
       onto every link. ``ldd $CRAY_MPICH_DIR/../gtl/lib/libmpi_gtl_cuda.so``
       says which CUDA that GTL wants: ``cray-mpich/9.1.0`` needs CUDA 13,
       ``8.1.28`` through ``9.0.1`` are the CUDA 12 ones. Match the pair, or
       ``module unload craype-accel-nvidia80`` if GPU-aware MPI is not needed
       -- nothing in this program passes device pointers to MPI
   * - ``mpi_f08`` generics that do not resolve under nvfortran
     - ``PIC_USE_LEGACY_MPI=ON``, already in the preset
   * - ``MPI_Win_allocate``'s generic does not resolve
     - ``PIC_MPICH_PERLMUTTER=ON``, already in the preset. Nothing here uses
       one-sided communication, so dropping those wrappers is free

The libsci trap, in full
------------------------

Worth spelling out because it is invisible from every angle that normally
matters, and because it is the reason a threaded run on this machine could not
be trusted before it was found.

``find_package(MPI)`` learns what to link by asking ``ftn`` for its link line,
and it asks *without* ``-fopenmp``. What comes back names the serial libsci --
``libsci_gnu``, ``libsci_gnu_mpi`` -- while ``find_package(OpenMP)`` asks with
the flag and gets ``libsci_gnu_mp`` and ``libsci_gnu_mpi_mp``. Both end up on
the link line, and whichever is listed first is what every BLAS and LAPACK call
in the binary resolves to.

Cray's serial libsci is not re-entrant. A short probe -- the same generalized
eigenproblem solved 200 times from eight OpenMP threads, each result compared
against the serial answer:

.. code-block:: text

   threads=8  failures=0   / 200     threaded libsci (what ftn itself links)
   threads=8  failures=190 / 200     serial libsci linked first

Nothing about that surfaces as a threading problem. It reads as an SCF that
will not converge, or an eigensolver that failed on a perfectly ordinary
structure, and it never appears in a traceback naming the BLAS. On a conformer
search it meant every metadynamics trajectory died on its first gradient, while
the *same* tblite inside a binary without the serial libsci ran fine.

``cmake/MqcDependencies.cmake`` now puts the threaded pair on the link line
ahead of everything, taken from what ``find_package(OpenMP)`` already resolved.
A configure on this machine says so:

.. code-block:: text

   -- Cray libsci: threaded copy linked first (.../libsci_gnu_mpi_mp.so;.../libsci_gnu_mp.so);
      the serial libsci the MPI wrapper names is not thread-safe

If that line is absent on a Cray, the build is the old one. ``ldd`` on the
binary settles it: ``libsci_gnu_mp.so.6`` must come before ``libsci_gnu.so.6``.

Running
-------

Threads come from ``OMP_NUM_THREADS`` and are used by the ab initio path, by
the fragment workers, and by CREST's sampling. xTB is the one method pinned to
a single thread, and not for speed -- see ``AGENTS.md``.

.. code-block:: bash

   # interactive, one CPU node
   salloc -N 1 -C cpu -q interactive -t 60

   # single rank, threads across the socket
   export OMP_NUM_THREADS=32
   srun -n 1 -c 32 ./build/mqc input.json

   # fragments over ranks, threads inside each
   export OMP_NUM_THREADS=8
   srun -n 16 -c 8 --cpu-bind=cores ./build/mqc input.json

A conformer search is single-rank by construction -- CREST samples on one rank
with OpenMP underneath, and the gradient interface it calls expects every rank
to take part, so the two cannot both be satisfied. Give it one rank and all the
threads:

.. code-block:: bash

   export OMP_NUM_THREADS=32
   srun -n 1 -c 32 ./build/mqc conformers.json

Adding another machine
======================

Two things, in this order.

Add a configure preset to ``CMakePresets.json`` when a configuration has caught
someone twice -- that threshold is deliberate, and ``CMAKE_STYLE.md`` explains
it. A preset is a name for a set of flags, not documentation of why they are
needed.

Then add a section here, in the shape of the Perlmutter one: the modules, the
configure line, and a table of symptoms with their causes. The table is the
part that earns its keep. Everything in the Perlmutter one cost somebody a day,
and every entry describes a failure whose message names neither the machine nor
the flag that fixes it.
