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

Methods
=======

Every method the backends implement is reachable from here, and
``python/examples/methods.py`` runs each one and checks it against PySCF on the
same geometry -- so this section is a description of something that is tested
rather than of something that ought to work.

Post-Hartree-Fock and DFT
^^^^^^^^^^^^^^^^^^^^^^^^^

``method`` is passed through to the deck unchanged, so every method the input
format accepts is reachable here by name:

.. code-block:: python

   mqc.MBE(water, level=0, method="mp2", basis="cc-pvdz").run()
   mqc.MBE(water, level=0, method="ccsd(t)", basis="cc-pvdz").run()
   mqc.MBE(water, level=0, method="ri-mp2", basis="cc-pvdz",
           aux_basis="cc-pvdz-rifit").run()

A name the parser does not know, and one it knows but nothing implements, are
both refused when the settings are read -- with the accepted spellings listed
-- rather than reaching the method factory, which has no error to return and
can only stop the process:

.. code-block:: python

   >>> mqc.MBE(water, level=0, method="mp3", basis="sto-3g").run()
   MQCError: unknown model.method 'mp3'. Accepted: gfn1, gfn2, hf, dft, mp2
   (also scs-, sos-, ri-, df-), ccsd, ccsd(t) (also ri-, df-), casscf, casci,
   mcscf, sapt0, efp2

Kohn-Sham takes the functional as its own argument, because the functional is a
separate field from the method:

.. code-block:: python

   mqc.MBE(water, level=0, method="dft", functional="b3lyp", basis="cc-pvdz").run()
   mqc.MBE(water, level=0, method="dft", functional="b2plyp", basis="cc-pvdz",
           aux_basis="cc-pvdz-rifit").run()

See the functional table in :doc:`input_files` for the names, and note that any
libxc name works too.

Unrestricted is not asked for either: the multiplicity selects it, on the
system or on the monomer.

.. code-block:: python

   oh = mqc.System(symbols=["O", "H"], coords=[[0, 0, 0], [0, 0, 0.97]],
                   multiplicity=2)
   oh.set_monomers([[0, 1]], multiplicities=[2])
   mqc.MBE(oh, level=0, method="hf", basis="sto-3g").run()

The keyword blocks
^^^^^^^^^^^^^^^^^^

Each of the deck's ``keywords`` blocks is a named argument taking a dict of
that block's own keys: ``scf``, ``correlation``, ``cc``, ``mcscf``, ``dft``,
``guess``, ``xtb`` and ``pcm``. They are dicts rather than a scalar per
setting, so a key added to the schema is reachable without the signature
moving, and named rather than left to ``keywords`` because what a method needs
should be visible from ``help(mqc.MBE)`` -- an active space is not an obscure
setting for a CASSCF, it is the calculation.

.. code-block:: python

   # the frozen core, which is on by default
   mqc.MBE(water, level=0, method="mp2", basis="cc-pvdz",
           correlation={"freeze_core": False}).run()

   # the amplitude solver, and the spin-orbital formulation
   mqc.MBE(water, level=0, method="ccsd(t)", basis="cc-pvdz",
           cc={"tolerance": 1e-10, "spin_adapted": False}).run()

   # the active space, named directly or chosen by AVAS
   mqc.MBE(water, level=0, method="casscf", basis="cc-pvdz",
           mcscf={"n_active_electrons": 4, "n_active_orbitals": 4}).run()
   mqc.MBE(water, level=0, method="casscf", basis="cc-pvdz",
           mcscf={"avas": {"orbitals": ["O 2p"], "threshold": 0.2}}).run()

   # an active space cut into subspaces with occupation windows
   mqc.MBE(water, level=0, method="casci", basis="cc-pvdz",
           mcscf={"n_active_electrons": 4, "n_active_orbitals": 4,
                  "ormas": {"subspaces": [1, 3], "min_electrons": [2, 0],
                            "max_electrons": [4, 2]}}).run()

   # the quadrature, and the initial guess
   mqc.MBE(water, level=0, method="dft", functional="tpss", basis="cc-pvdz",
           dft={"grid_level": 5}, guess={"type": "sad"}).run()

Properties
^^^^^^^^^^

``properties`` is the deck's block of the same name, beside ``keywords``
rather than inside it, and the distinction is worth keeping: keywords say how
to compute the wave function and change the number that comes out, properties
ask for something further to be done with it and change nothing. That is why a
bonding analysis leaves the driver as an energy:

.. code-block:: python

   mqc.MBE(water, level=0, method="hf", basis="6-31g",
           properties={"bonding_analysis": {"type": "gms_quao"}}).run()

Interaction energies
^^^^^^^^^^^^^^^^^^^^

SAPT0 and EFP2 return the interaction between two systems rather than the
energy of one, so for both of them the partition *is* the physics -- the
monomers are the interacting pieces, and ``level=0`` is right because there is
nothing for an expansion to expand.

SAPT0 wants exactly two monomers, and its decomposition comes back on the
result. The terms travel through the output document rather than through the
energy, so the run has to write files:

.. code-block:: python

   >>> dimer.auto_monomers()
   >>> result = mqc.MBE(dimer, level=0, method="sapt0", basis="sto-3g").run(
   ...     label="dimer", write_to_file=True)
   >>> result.energy            # the total interaction
   0.0032483374
   >>> result.sapt["elst10"]    # and the term another code would quote
   -0.0033...

EFP2 solves no wave function of its own: each fragment carries one already,
computed when its potential was made. So it needs a potential per monomer,
named on the system, and ``driver="makefp"`` is what writes them -- it builds
one for the whole system it is given and names it after the label:

.. code-block:: python

   monomer.set_monomers([[0, 1, 2]])
   mqc.MBE(monomer, level=0, method="hf", basis="6-31g",
           driver="makefp").run(label="water", write_to_file=False)

   dimer.auto_monomers()
   dimer.set_fragment_potentials(["water.efp", "water.efp"])
   mqc.MBE(dimer, level=0, method="efp2").run(write_to_file=False).energy

Every monomer needs one. A mixed quantum/EFP system is a deck-only feature,
and a system with no potentials at all is refused rather than returning the
zero it used to.

Anything without a named argument
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``keywords`` merges straight into the deck's ``keywords`` block, and wins over
everything above, so a setting does not have to wait for a Python argument to
be reachable:

.. code-block:: python

   mqc.MBE(water, level=0, method="ccsd(t)", basis="cc-pvdz",
           keywords={"cc": {"tolerance": 1e-10},
                     "correlation": {"n_frozen_core": 1},
                     "scf": {"guess": "sad"}}).run()

Density fitting
^^^^^^^^^^^^^^^

Hartree-Fock on the CPU can fit J and K against an auxiliary basis instead of
computing exact four-index integrals:

.. code-block:: python

   mqc.MBE(water, level=0, method="hf", basis="sto-3g",
           aux_basis="def2-universal-jkfit", density_fitting=True)

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

``methods.py``
   Every method the backends implement, run through this interface and checked
   against PySCF fed this repository's own basis data -- Hartree-Fock,
   unrestricted, Kohn-Sham through the double hybrids, the MP2 family, coupled
   cluster with and without the triples and with and without fitting, CASCI,
   CASSCF, ORMAS, AVAS, SAPT0 and EFP2. Where an identity says more than a
   reference would, the identity is checked instead: RI against conventional,
   spin-adapted against spin-orbital, CASSCF against CASCI on the same space.
   The last cases are the refusals -- an unknown method, ``driver="optimize"``,
   EFP2 with no potentials -- which are there because each of them used to
   return a plausible number or take the interpreter down.

``energy_screened_mbe.py``
   A two-pass calculation: the whole expansion at a cheap level, then the terms
   whose contribution exceeded a threshold recomputed with DFT. The criterion
   is energetic rather than geometric, which is the point -- a distance cutoff
   discards a strong interaction that happens to be far away.

What it cannot do
=================

``driver="optimize"`` is the one driver a session cannot host. The optimizer
drives ``run_calculation`` rather than being driven by it, so it works from a
deck for a single molecule and nowhere else. It is refused before the workers
are told anything, so it raises rather than -- as it used to -- aborting the
communicator and taking the interpreter with it:

.. code-block:: python

   >>> mqc.MBE(water, level=0, method="hf", basis="sto-3g",
   ...         driver="optimize").run()
   MQCError: driver 'optimize' is not available through a session or the C
   API. ... Ask for 'energy' or 'gradient' here, or run the deck through the
   mqc executable.

Gradients are available from every backend, including the CPU one -- that used
to be the entry in this section and no longer is:

.. code-block:: python

   >>> mqc.MBE(water, level=0, method="hf", basis="sto-3g",
   ...         driver="gradient").run(label="grad", write_to_file=True).gradient_norm
   0.0861241914007507

``write_to_file=True`` is not incidental: ``gradient_norm`` is read back out of
``output_<label>.json``, so a run that wrote nothing returns ``None`` whether or
not it computed a gradient. The norm is all this layer exposes -- enough to tell
a run that produced a gradient from one that did not, not enough to tell a
correct gradient from one with a sign error. The components are in that same
JSON.

What is refused is a *particular combination*, not a backend, and the list is in
:doc:`capabilities` under gradient calculations -- an open-shell double hybrid,
a frozen-core MP2, spin-scaled MP2, coupled cluster, and a few others. Each
raises with a message naming the missing piece rather than returning a number
from a formula that does not apply.

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
