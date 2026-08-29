==========
Benchmarks
==========

Time this build on this machine, and say what changed since last time.

.. code-block:: bash

   cmake --build build -j
   python3 benchmarks/run_benchmarks.py --exe build/mqc --record   # first time
   python3 benchmarks/run_benchmarks.py --exe build/mqc            # after that

About twelve minutes on forty cores for a build in reasonable shape; ``--quick``
is three and skips the thread ladder and the fragmented cases. It is a
performance regression suite and a "does this build perform the way it should
here" report, because those want the same measurements.

Why each phase is there
=======================

Every phase exists because a plain wall-clock number was misleading during a
real investigation. None of it is generic benchmarking furniture.

**The build configuration, first.** The BLAS choice alone moves a density
functional run by a factor of five, is set by one CMake cache variable, and
cannot be recovered from inside the program. ``mqc --version`` reports it:

.. code-block:: text

   build: type=Release compiler=GNU 15.1.0 blas=Intel10_64lp_seq openmp=yes

A timing that does not say which BLAS was linked is not reproducible, so the
suite records it beside the numbers and warns when a baseline came from a
different build.

**Hartree-Fock beside every functional.** A pure functional carries no exact
exchange and should cost *less* than Hartree-Fock. It recently cost twenty-six
times as much, and no absolute number would have said so -- only the ratio to a
reference sharing the same integrals. The ratio is size dependent, since the
quadrature and the Fock build scale differently, so compare like with like.

**The same case at every thread count.** A stage that does not parallelise is
invisible at one thread count: it looks expensive rather than serial. The
exchange-correlation quadrature was serial for the whole life of the code, and
this ladder is the report that would have shown it.

**The correlated methods, each at its own size.** MP2 goes as the fifth power of
the basis and coupled cluster as the sixth or seventh, so no single system
serves them. MP2 and RI-MP2 run on ten waters in cc-pVDZ, the coupled-cluster
pair on five in 6-31G, each placed where its own work rather than the reference
SCF is what moves.

**A fragmented case, serial and under MPI.** Fragment work pins itself to one
OpenMP thread and spreads over ranks instead, so a change that helps a single
molecule can hurt the fragment path. One did: threaded BLAS is five times faster
on one molecule and thirty-one per cent slower on four ranks. Without this phase
that trade reads as a pure win.

What it recommends
==================

The run ends with settings rather than only numbers, each printing the
measurement behind it -- advice without its measurement is an assertion that
happens to be in a program.

* the thread count past which this machine stops improving, taken as the fewest
  threads within five per cent of the best time on the ladder
* a warning when parallel efficiency at that point is under sixty per cent,
  which means a serial stage and points at the ladder to find it
* whether fragment work wants threads or ranks *at the size measured*
* whether RI-MP2 is actually faster than conventional MP2 on this hardware, and
  it says so plainly when density fitting has stopped paying for itself

Judging a regression
====================

Against the measured spread, never a fixed percentage. A case repeating to 0.3
per cent that moves by 3 has regressed; one repeating to 16 per cent that moves
by 3 has not moved at all. One threshold for both is how a suite becomes noise
that people stop reading.

The band is twice the measured spread and never under five per cent, because the
spread understates the noise it stands in for: repeats inside one invocation run
back to back on a warm cache at a settled clock, while a baseline recorded days
earlier saw none of that. Cases repeating to under one per cent within a run
still drift about three between runs.

Stages are compared as well as totals. A correlated method's own work is a
minority of a run the reference SCF dominates, so a regression in the integral
transform moves its stage plainly and the total barely at all.

Baselines are per host and are not committed. Reference numbers in a repository
are meaningless on hardware unlike the machine that produced them.

Adding a case
=============

Add it to ``ENERGY_CASES``, ``CORRELATED_CASES`` or ``GRADIENT_CASES`` in
``benchmarks/run_benchmarks.py``. Each case carries its own geometry and basis,
because the methods do not share a cost scale. Decks are generated into a
temporary directory at run time rather than committed, since they differ only in
a couple of fields; the geometries are in ``benchmarks/geometries/``.

Keep the directory out of ``validation/``: the CPU validation suite is
generated, and its generator deletes decks under ``validation/inputs/cpu/mqc/``
that it did not write.
