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

Density fitting
===============

Hartree-Fock on the CPU can fit J and K against an auxiliary basis instead of
computing exact four-index integrals:

.. code-block:: python

   mqc.MBE(water, level=0, method="hf", basis="sto-3g",
           aux_basis="def2-universal-jkfit", density_fitting=True)

Post-Hartree-Fock and DFT
^^^^^^^^^^^^^^^^^^^^^^^^^

``method`` is passed through to the deck unchanged, so every method the input
format accepts is reachable here by name:

.. code-block:: python

   mqc.MBE(water, level=0, method="mp2", basis="cc-pvdz").run()
   mqc.MBE(water, level=0, method="ccsd(t)", basis="cc-pvdz").run()
   mqc.MBE(water, level=0, method="ri-mp2", basis="cc-pvdz",
           aux_basis="cc-pvdz-rifit").run()

Kohn-Sham takes the functional as its own argument, because the functional is a
separate field from the method:

.. code-block:: python

   mqc.MBE(water, level=0, method="dft", functional="b3lyp", basis="cc-pvdz").run()
   mqc.MBE(water, level=0, method="dft", functional="b2plyp", basis="cc-pvdz",
           aux_basis="cc-pvdz-rifit").run()

See the functional table in :doc:`input_files` for the names, and note that any
libxc name works too.

Anything without a named argument
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``keywords`` merges straight into the deck's ``keywords`` block, so a setting does
not have to wait for a Python argument to be reachable:

.. code-block:: python

   mqc.MBE(water, level=0, method="ccsd(t)", basis="cc-pvdz",
           keywords={"cc": {"tolerance": 1e-10},
                     "correlation": {"n_frozen_core": 1},
                     "scf": {"guess": "sad"}}).run()

   mqc.MBE(water, level=0, method="dft", functional="tpss", basis="cc-pvdz",
           keywords={"dft": {"grid_level": 5}}).run()


It has to be asked for. Setting ``aux_basis`` alone does not turn it on,
because that name carries a default -- inferring the request from its presence
would mean every Hartree-Fock quietly fitted. The difference on water/sto-3g is
8.6e-5 Hartree: large enough to matter, small enough to be mistaken for
convergence noise.

cuEST is the other way round. It has no four-index path at all, so it fits
whether or not the flag is set, and the flag is ignored there.

In a deck the same switch is ``keywords.scf.density_fitting``.

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

Gradients are not available from the CPU Hartree-Fock backend. libcint has the
derivative entry points, so this is a gap rather than a wall, but an untested
gradient is worse than an absent one -- so it raises:

.. code-block:: python

   >>> mqc.MBE(water, level=0, method="hf", driver="gradient").run()
   MQCError: the CPU backend has no gradients yet; run an energy, or build with cuEST

xTB gradients work normally, and so do cuEST's. Only the CPU Hartree-Fock path
refuses.

``backend=`` picks the integral backend -- ``"cuest"``/``"gpu"``,
``"libcint"``/``"cpu"``, or omitted for the build's default. A request the build
cannot honour raises rather than falling back:

.. code-block:: python

   >>> mqc.MBE(water, level=0, method="hf", basis="sto-3g", backend="gpu").run()
   MQCError: This calculation needs the cuEST integral backend; build with CMake
   and -DMQC_ENABLE_CUEST=ON

Continuum solvation is a named argument, and cuEST-only -- the CPU backend has no
cavity and refuses it rather than returning a gas-phase energy:

.. code-block:: python

   >>> mqc.MBE(water, level=0, method="dft", functional="pbe0", basis="def2-svp",
   ...         pcm={"dielectric": 78.39}).run(write_to_file=False).energy

The dict takes the same keys as ``keywords.pcm`` in a deck, documented under
:ref:`input_files`. ``dielectric`` is required; everything else has a default.

Correlated methods on the CPU backend are closed-shell only. MP2, coupled
cluster and the double hybrids each need separate alpha and beta transforms for an
open-shell reference, and are refused in the same way rather than approximated.
Hartree-Fock and Kohn-Sham DFT both run unrestricted, and pick it up from the
multiplicity without being told:

.. code-block:: python

   >>> ch3 = mqc.System.from_xyz("ch3.xyz", multiplicity=2)
   >>> ch3.set_monomers([[0, 1, 2, 3]])
   >>> mqc.MBE(ch3, level=0, method="dft", functional="pbe",
   ...         basis="cc-pvdz").run(write_to_file=False).energy
   -39.7691396276

Nothing there asks for the unrestricted path: the multiplicity selects it. The
range-separated hybrids work on an open shell as well, with the attenuated
exchange pass run once per spin.
