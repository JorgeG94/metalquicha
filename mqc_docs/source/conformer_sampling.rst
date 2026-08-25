==================
Conformer sampling
==================

Metalquicha can link CREST_ and drive its conformer and ensemble sampling with
gradients served from mqc's own code, in the same process. No file is written
between the two and no subprocess is spawned.

Without ``MQC_ENABLE_CREST``, asking for this driver is refused by name and
every other driver runs unchanged.

.. _CREST: https://github.com/crest-lab/crest

Running one
===========

.. code-block:: json

   {
     "schema": { "name": "mqc-frag", "version": "1.0" },
     "molecules": [
       { "xyz": "water.xyz", "molecular_charge": 0, "molecular_multiplicity": 1 }
     ],
     "model": { "method": "hf", "basis": "sto-3g" },
     "driver": "conformers"
   }

``"conformer"`` and ``"crest"`` are accepted as well. The ``model`` block is the
*refinement* method -- what the survivors are re-evaluated with. The sampling
method is CREST's own default, GFN2-xTB in process, and is not chosen here.

The search writes its ensemble into the working directory the way CREST always
does: ``crest_conformers.xyz``, ``crest_rotamers.xyz``, ``crest_best.xyz`` and
``crest.energies``, alongside the usual scratch. The absolute energies on the
comment line of each structure are the refinement method's -- for the deck
above, HF/STO-3G:

.. code-block:: text

     3
           -74.96203654
    O         -0.0000000000        0.0636623383        0.0000000001

``crest.energies`` holds energies relative to the lowest, in kcal/mol.

.. code-block:: text

   conformer sampling: CREST samples, this program refines
   ...
   number of unique conformers for further calc            1
   refinement calculations: 1432

Where it does not work
======================

``conformers`` runs from an input deck, for a single molecule, through the
``mqc`` executable. It is driven above ``run_calculation`` rather than inside
it -- the same position as ``optimize``, and for the same reason -- so a
multi-molecule deck, a session and the C API all refuse it by name rather than
letting it fall through.

It also refuses to run on more than one rank; see below.

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

How often refinement runs
=========================

.. important::

   The refinement level is applied after **every** optimization-and-sorting
   cycle, not once at the end. On the water deck above that is six passes, and
   it is why a molecule with a single conformer still costs 1,432 refinement
   calculations rather than one.

CREST's iMTD-GC workflow is a loop: metadynamics, optimize every snapshot,
CREGEN sorting, regular MDs on the survivors, genetic crossing, repeat until no
new lowest conformer appears. Anything registered through ``env%addrefine``
runs inside that loop, re-ranking the whole surviving set each time round.

At semiempirical or small-basis cost that is invisible. At RI-MP2 or DFT on a
real molecule it is the difference between an afternoon and a week, and it is
the first thing to look at if a run is taking longer than expected. There is
currently no deck keyword to move refinement to the end; changing it means
registering a different ``refine`` stage, which is a code change.

What a refinement call costs
============================

Each call is an **energy**, not a gradient. CREST re-ranks with the energy and
never reads the gradient at this stage, and an RI-MP2 gradient is several times
its energy -- so asking for one would be most of the cost of the thing you
actually want. The gradient handed back is zero.

Measured on the water deck, three runs of the same input:

.. list-table::
   :header-rows: 1
   :widths: 46 18 36

   * - Configuration
     - Wall time
     - Refinement
   * - sampling only, no refinement level
     - 18.9 s
     - --
   * - refinement asking for gradients, writing JSON per call
     - 26.7 s
     - 7.8 s, 5.4 ms/call
   * - refinement asking for energies, no per-call output
     - 20.6 s
     - 1.7 s, 1.2 ms/call

At this size refinement is **not** the bottleneck -- the sampling is 92% of the
run. That inverts completely with an ab initio method, which is the case the
split exists for.

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
``compute_energy_and_forces`` expects every rank to call it. The two cannot both
be satisfied, and the failure if it is not caught is a hang rather than an
error, so the driver refuses:

.. code-block:: text

   conformer sampling runs on a single rank. [...] Re-run without mpirun, or with -n 1.

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
