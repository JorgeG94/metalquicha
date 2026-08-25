==================
Conformer sampling
==================

Metalquicha can link CREST_ and drive its conformer and ensemble sampling with
gradients served from mqc's own code, in the same process. No file is written
between the two and no subprocess is spawned.

.. _CREST: https://github.com/crest-lab/crest

.. note::

   There is no deck keyword for this yet. Everything below is reachable from
   Fortran that links ``libmetalquicha.a``, and the three ``check_crest_*``
   programs under ``validation/`` are working examples. A ``driver`` value that
   runs a search from a JSON input is the next piece of work, not a thing this
   page is describing.

What it is for
==============

A conformer search asks for an enormous number of gradients -- an iMTD-GC run
on a small molecule measures in the tens of thousands -- and almost none of
them need to be good. The arrangement this supports is the one that follows
from that: a semiempirical method does the sampling, and an ab initio method
sees only what survives.

On water, with an HF/STO-3G refinement, that split measured **1,432 mqc calls
against roughly 56,000 gradients in total**. The ratio is what makes ab initio
refinement affordable at all.

Building
========

.. code-block:: bash

   cmake -B build -DMQC_ENABLE_CREST=ON -DMQC_ENABLE_TBLITE=ON \
                  -DWITH_TBLITE=ON -DWITH_GFN0=ON

``MQC_ENABLE_CREST`` is ``OFF`` by default. Two things about the flags:

``WITH_GFN0=ON`` is not optional in practice. ``crest_search_imtdgc`` generates
its metadynamics length from a flexibility measure that calls GFN0-xTB for
Wiberg bond orders, and stops with *"Compiled without GFN0-xTB support!"*
without it. That is true regardless of who supplies the gradients.

The ``WITH_*`` options belong to CREST rather than to this project, and CREST
resolves its own dependencies -- toml-f, gfnff, gfn0, libpvol, lwoniom -- through
``config/modules/crest-utils.cmake``, which tries a local subproject, an
installed package, then a fetch. Nothing needs declaring here for any of them.

``tblite`` is the exception, because both projects use it and only one copy can
exist in a binary. CREST guards it with ``if(NOT TARGET "tblite::tblite")``, so
it reuses the one this project already provides.

.. warning::

   Upstream CREST pins tblite to a ``xtb_solvation`` branch whose last commit
   predates the release this project uses by about eighteen months, and the two
   disagree about ``eeq_guess``. The fork this integration builds against
   carries a one-line fix; a stock CREST will fail to compile against a modern
   tblite. ``MQC_CREST_REPOSITORY`` and ``MQC_CREST_TAG`` select which CREST is
   fetched.

Running on one rank
===================

CREST samples on a single rank with OpenMP threads underneath, while
``compute_energy_and_forces`` expects every rank to call it. A driver must
therefore refuse above one rank rather than deadlock.

This is checked at run time, not at build time, so one binary serves both and
``mpirun -n 1`` still works. ``MQC_ENABLE_MPI`` can stay ``ON``.

How the gradients cross
=======================

CREST dispatches every gradient through one ``select case`` over a job type.
The fork adds ``jobtype%callback``, which calls a procedure pointer the host
installs:

.. code-block:: fortran

   use crest_calculator, only: calculation_settings, jobtype

   type(calculation_settings) :: level

   level%id = jobtype%callback
   level%external_engrad => my_engrad

The pointer's interface takes plain arrays rather than any CREST derived type:

.. code-block:: fortran

   subroutine my_engrad(nat, at, xyz, energy, grad, iostatus)
      integer,  intent(in)  :: nat, at(nat)
      real(wp), intent(in)  :: xyz(3, nat)     !! bohr
      real(wp), intent(out) :: energy          !! Eh
      real(wp), intent(out) :: grad(3, nat)    !! Eh/bohr
      integer,  intent(out) :: iostatus        !! non-zero aborts the run
   end subroutine

Units are the trap worth naming. ``xyz`` arrives in **bohr** and ``grad`` must
be returned in **Eh/bohr**. A callback that gets this wrong returns a correct
energy with a gradient wrong by a constant factor, which nothing downstream
will notice.

Selecting ``jobtype%callback`` without installing a pointer reports a non-zero
status rather than returning zero, because a silent zero is indistinguishable
from a converged result and would propagate through the whole ensemble.

Configuring a run
=================

``systemdata`` carries on the order of a hundred fields and the algorithms read
more of them than their entry points name, so it is not built by hand. Instead,
hand CREST's own parser an argument vector built **in memory** -- no command
line and no configuration file -- exactly as ``crest_main`` does with the real
argv:

.. code-block:: fortran

   character(len=:), allocatable :: argv(:)

   call tim%init(20)
   allocate (character(len=64) :: argv(1))
   argv(1) = "start.xyz"
   call parseflags(env, argv, 1)

CREST then applies all of its own defaults, and only the calculator is
overridden afterwards. The starting structure is still a file; it is read once,
at setup, and every gradient after that crosses in memory.

Splitting sampling from refinement
==================================

Add a second calculation level tagged as a refinement stage, and register it.
CREST applies it to the ensemble the search leaves behind:

.. code-block:: fortran

   use crest_data, only: refine

   refine_level%id = jobtype%callback
   refine_level%external_engrad => mqc_refine_engrad
   refine_level%refine_lvl = refine%singlepoint
   call env%calc%add(refine_level)
   call env%addrefine(refine%singlepoint)

Level one stays whatever ``parseflags`` chose -- GFN2-xTB through tblite, in
process -- and level two is mqc. This is the same arrangement CREST uses
internally for two xTB levels in ``legacy_wrappers.f90``.

Checking a build
================

Three manual programs under ``validation/``, built only when
``MQC_ENABLE_CREST=ON``:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Program
     - What it establishes
   * - ``check_crest_callback``
     - A gradient crosses the boundary unaltered. Runs an HF/STO-3G gradient
       twice, directly and through CREST, and compares exactly.
   * - ``check_crest_search``
     - A whole iMTD-GC run can be configured from memory and driven through the
       hook.
   * - ``check_crest_refine``
     - The sampling/refinement split holds: xTB samples, mqc refines, and mqc's
       energies are the ones reported.

.. code-block:: bash

   cmake --build build --target check_crest_callback
   ./build/check_crest_callback

They are excluded from ``all`` and from ``ctest``, like the other
``validation/`` programs, because a full search takes minutes and writes an
ensemble into the working directory.

Licensing
=========

CREST is LGPL-3.0 and metalquicha is MIT, the same situation as
``MQC_ENABLE_DLFIND`` and the reason both are off by default. LGPL permits
linking from any licence, but a distributed binary carries the obligation that
a recipient can relink against a modified CREST -- in practice, dynamic linking
or shipping the objects. That is a decision about distribution rather than
about building, and it has not been made here.
