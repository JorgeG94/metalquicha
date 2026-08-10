================
Python interface
================

Metalquicha can be driven from Python. The Fortran program still does the
calculation and still owns MPI; Python sets up the molecule, decides which
fragments to compute, and reads the answers back.

The division is deliberate. Python runs on rank 0 and nowhere else -- it never
sees the other ranks and never needs to. When it asks for a calculation, the
Fortran side distributes the work across the whole job and returns when the
work is done. So a script that looks single-threaded runs on as many nodes as
it was launched with, and ``mpirun -np 64 python script.py`` is a valid way to
start one.

Building it
===========

The interface loads ``libmqc.so``, which is a separate target from the
executable:

.. code-block:: bash

   cmake -B build -DMQC_ENABLE_TBLITE=ON -DMQC_ENABLE_LIBCINT=ON
   cmake --build build --target mqc_shared

Then point Python at it. The package looks for the library next to an in-tree
``build/`` directory, so this is usually enough:

.. code-block:: bash

   export PYTHONPATH=$PWD/python
   python3 -c "import mqc; print(mqc.__all__)"

and if the library is elsewhere, ``MQC_LIBRARY`` overrides the search:

.. code-block:: bash

   export MQC_LIBRARY=/path/to/libmqc.so

Which backends are available depends on how the library was built. ``gfn2`` and
``gfn1`` need ``MQC_ENABLE_TBLITE``; ``hf`` on the CPU needs
``MQC_ENABLE_LIBCINT``. A method the build does not have fails when it is run,
not when it is named.

A first calculation
===================

.. code-block:: python

   import mqc

   with mqc.session():
       water = mqc.System(
           symbols=["O", "H", "H"],
           coords=[[0.0, 0.0, 0.1008],
                   [0.0, 0.7725, -0.4678],
                   [0.0, -0.7725, -0.4678]],
       )
       water.set_monomers([[0, 1, 2]])
       result = mqc.MBE(water, level=0, method="gfn2").run(write_to_file=False)
       print(result.energy)

Everything happens inside ``mqc.session()``. It starts MPI on entry and stops
it on exit, and outside it there is nothing to talk to.

Note ``set_monomers`` even though nothing is being fragmented. A partition is
required whether or not the expansion uses one, and ``level=0`` means the whole
system is computed at once. For a molecule made of separate pieces,
``auto_monomers()`` finds them; it refuses a single covalently bonded molecule,
because where to cut such a thing is a chemical decision and not one to guess.

Fragmenting
===========

.. code-block:: python

   cluster = mqc.System.from_xyz("water20.xyz")
   cluster.auto_monomers()
   result = mqc.MBE(cluster, level=2, method="gfn2").run(label="w20")

With exactly two monomers an ``MBE(2)`` is not an approximation -- the
two-body term cancels the monomer terms and what is left is the supermolecular
energy. That identity is worth knowing because it is the cheapest available
test of the fragmentation machinery, and it is what
``python/examples/backends.py`` checks.

Examples
========

``python/examples/`` holds runnable scripts, not fragments:

``backends.py``
   Standalone and fragmented, xTB and Hartree-Fock -- the four paths through
   the program that the interface can reach. Each case asserts against a
   reference and the script exits non-zero if any is wrong, so it runs in CI
   rather than sitting there going quietly stale. The Hartree-Fock references
   are PySCF; the fragmented cases check the expansion against the same dimer
   computed whole.

``energy_screened_mbe.py``
   A two-pass calculation: the whole expansion at a cheap level, then the terms
   whose contribution exceeded a threshold recomputed with DFT. The criterion
   is energetic rather than geometric, which is the point -- a distance cutoff
   discards a strong interaction that happens to be far away.

What it cannot do
=================

Gradients are not available from the CPU Hartree-Fock backend; it refuses them
rather than returning something untested. Density fitting is likewise not
reachable from the Python side yet -- the auxiliary basis is accepted and
ignored, and the calculation runs with exact integrals.
