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

What this buys you over a deck is *decisions in between*. An input file states
one calculation; a script can run a cheap expansion, look at which terms
mattered, and spend the budget on those -- which is what
``python/examples/energy_screened_mbe.py`` does, and what no deck can express.

The whole interface, at a glance
================================

Four objects and a context manager:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Object
     - What it is
   * - ``mqc.session()``
     - The MPI lifetime. Everything happens inside it; outside it there is
       nothing to talk to.
   * - ``mqc.System``
     - A molecule: atoms, the partition into monomers, the bonds. Optionally
       bond orders, atomic charges, and effective fragment potentials.
   * - ``mqc.MBE``
     - A calculation over that system: the method, the settings, and the term
       list. ``level=0`` means "the whole thing, once".
   * - ``mqc.Result``
     - What came back -- the energy -- and the way to the rest of it: the
       per-term breakdown, the gap, the gradient norm, the SAPT decomposition.
   * - ``mqc.MQCError``
     - Anything Fortran refused to do, with the reason it gave.

.. code-block:: python

   import mqc

   with mqc.session():
       water = mqc.System.from_xyz("water.xyz")
       water.set_monomers([[0, 1, 2]])
       print(mqc.MBE(water, level=0, method="hf", basis="cc-pvdz").run().energy)

Building and importing
======================

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

Which methods are available depends on how the library was built. ``gfn2`` and
``gfn1`` need ``MQC_ENABLE_TBLITE``; everything ab initio needs
``MQC_ENABLE_LIBCINT``; the Kohn-Sham functionals need ``MQC_ENABLE_LIBXC``;
``backend="gpu"`` needs ``MQC_ENABLE_CUEST``. Asking for a method this build
does not carry raises, naming the CMake option -- it never quietly substitutes
another one.

Running under MPI
=================

**Every rank runs the script until ``session()`` is entered, and only rank 0
comes out of it.** Under ``mpirun -np 64 python script.py`` there are 64
interpreters; 63 of them enter ``session()`` and never return, spending the
rest of their lives inside Fortran waiting to be told what to compute.

.. code-block:: python

   print("before the session: every rank runs this")
   with mqc.session():
       print("inside: rank", mqc.rank(), "of", mqc.n_ranks())
       ...

Run that under ``-np 2`` and the first line appears twice, the second once.
The consequences are worth stating plainly:

* Anything above the session happens N times, on a shared filesystem. Read
  structures and parse arguments *inside*.
* Anything below it happens once, on rank 0.
* **The session must be closed.** An uncaught exception on rank 0 leaves the
  other 63 blocked on a broadcast that never arrives -- a hang rather than a
  traceback, and on a batch system a hang that burns the allocation.
  ``with mqc.session()`` closes it even when the body raises; the bare
  ``mqc.begin()`` / ``mqc.end()`` pair exists for a REPL.

A fragmented run is distributed automatically: the term list is handed out
across the job and the energies come back to rank 0. The number is the same
one a serial run produces -- a three-water MBE(2) gives
``-224.86311196000983`` at ``-np 1`` and at ``-np 2``, because the expansion is
a sum and the sum does not care who computed which term.

The molecule
============

Geometry
--------

Coordinates are Angstrom here and Bohr inside; the boundary converts. Atom
indices are 0-based, as Python expects.

.. code-block:: python

   water = mqc.System(
       symbols=["O", "H", "H"],
       coords=[[0.0, 0.0, 0.1008],
               [0.0, 0.7725, -0.4678],
               [0.0, -0.7725, -0.4678]],
       charge=0, multiplicity=1,
   )

   cluster = mqc.System.from_xyz("water20.xyz", charge=0, multiplicity=1)

The partition
-------------

A partition into monomers is required whether or not an expansion will use it:
the validator wants one, and for SAPT0 and EFP2 the monomers *are* the
physics. ``auto_monomers()`` finds the connected components under a
covalent-radius criterion -- exactly right for a cluster of intact molecules,
and it returns one monomer for anything held together by covalent bonds,
because where to cut such a thing is a chemical decision and not one to guess.

.. code-block:: python

   cluster.auto_monomers()                       # a cluster of molecules
   peptide.set_monomers([[0, 1, 2, 3], [4, 5, 6]])   # say it yourself
   peptide.set_bonds([(3, 4)])                       # and say what you cut

Monomer indices in a term list are 1-based -- that is the expansion's own
convention, and the breakdown CSV uses it too. Atom indices are 0-based
everywhere.

Bonds
-----

A bond crossing a monomer boundary is *derived* to be broken, from the
partition, and capped with hydrogen. What cannot be derived is a bond you never
mentioned, so a system that cuts undeclared bonds is refused rather than
silently fragmented into radicals:

.. code-block:: python

   >>> peptide.missing_bonds()      # look before you get there
   1
   >>> mqc.MBE(peptide, level=2, method="gfn2").run()
   MQCError: the monomers cut 1 bond(s) that were never declared, starting with
   atoms 3 and 4; declare them with mqc_system_set_bonds so the fragments are
   capped, or set unchecked_input if the partition is deliberate

``perceive_bonds()`` fills the list in from the geometry and declares it, which
is the quick answer when the connectivity is not in question. Passing no bonds
at all is a statement too -- ``set_bonds([])`` says "nothing to cut" -- and
that is different from never having said anything, which is what
``bonds_declared`` reports.

Bond orders and charges
-----------------------

Two questions a partition often depends on, answered on the system rather than
by a calculation you have to assemble yourself:

.. code-block:: python

   cluster.compute_bond_orders(variant="gfn2")   # one xTB single point
   cluster.bond_order(4, 7)                      # Wiberg-Mayer, 0-based

   cluster.compute_charges(scheme="chelpg", basis="6-31g")   # one RHF
   cluster.charges()

Both are computed once and read many times, because a caller trying twenty
trial fragmentations does not want twenty calculations. Bond orders separate a
real bond from a hydrogen bond in a way a distance rule cannot; charges come
from CHELPG or Mulliken, and the two disagree by design. See
:doc:`charges_and_bond_orders`.

Effective fragment potentials
-----------------------------

``set_fragment_potentials`` names one ``.efp`` file per monomer, which is what
makes ``method="efp2"`` mean anything -- see `Interaction energies`_ below.

Asking for a calculation
========================

``mqc.MBE`` holds the settings and the term list. Its arguments become the same
JSON document ``mqc`` reads from a deck, less the molecules, which is what
crosses to the other ranks:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Argument
     - Meaning
   * - ``level``
     - Expansion order. **0 is a single calculation on the whole system**, and
       is how every method that is not a fragmentation method is run.
   * - ``method``, ``basis``, ``aux_basis``, ``functional``
     - The model. Same spellings a deck uses.
   * - ``driver``
     - ``"energy"``, ``"gradient"``, ``"hessian"``, ``"makefp"``.
   * - ``density_fitting``
     - Fit J and K. Asked for, never inferred from ``aux_basis``.
   * - ``cutoffs``
     - Per-level distance screening, e.g. ``{"dimer": 5.0, "trimer": 4.0}``.
   * - ``backend``
     - ``"cuest"``/``"gpu"``, ``"libcint"``/``"cpu"``, or omitted for the
       build's default.
   * - ``scf``, ``correlation``, ``cc``, ``mcscf``, ``dft``, ``guess``,
       ``xtb``, ``pcm``
     - One dict per keyword block of the deck, taking that block's own keys.
   * - ``properties``
     - The deck's ``properties`` block -- things to report once the wave
       function exists.
   * - ``checkpoint``
     - Append each fragment as it finishes, and resume from what is there.
   * - ``allow_crap_scf``
     - Keep a non-converged SCF instead of failing. Off unless asked for.
   * - ``unchecked``
     - Skip the semantic checks on the system. For a partition you know is
       unusual and mean.
   * - ``verbosity``
     - ``debug``, ``verbose``, ``info`` (default), ``performance``,
       ``warning``, ``error``.
   * - ``keywords``, ``system_options``
     - Escape hatches for the deck's ``keywords`` and ``system`` blocks. They
       are merged last and win over everything above.

``settings()`` returns that document. Print it when a run surprises you:

.. code-block:: python

   >>> mqc.MBE(cluster, level=2, method="hf", basis="sto-3g").settings()["keywords"]
   {'fragmentation': {'method': 'MBE', 'level': 2}}

``run(label="mqc", write_to_file=True)`` computes across the whole job. The
label names the output -- ``output_<label>.json`` and
``output_<label>_fragments.csv`` -- and is rejected if it contains dots or
slashes, because Fortran treats it as a filename and would quietly eat them:
``"tau0.001"`` and ``"tau0.002"`` would otherwise write over each other.

Methods
=======

Every method the backends implement is reachable from here, and
``python/examples/methods.py`` runs each one and checks it against PySCF on the
same geometry -- so this section describes something that is tested rather than
something that ought to work.

Post-Hartree-Fock and DFT
-------------------------

``method`` is passed through to the deck unchanged, so every method the input
format accepts is reachable here by name:

.. code-block:: python

   mqc.MBE(water, level=0, method="mp2", basis="cc-pvdz").run()
   mqc.MBE(water, level=0, method="ccsd(t)", basis="cc-pvdz").run()
   mqc.MBE(water, level=0, method="ri-mp2", basis="cc-pvdz",
           aux_basis="cc-pvdz-rifit").run()

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - ``method``
     - Also needs
     - Notes
   * - ``gfn1``, ``gfn2``
     - --
     - tblite. Solvation through ``xtb={"solvent": "water"}``
   * - ``hf``
     - ``basis``
     - Restricted or unrestricted, from the multiplicity
   * - ``dft``
     - ``basis``, ``functional``
     - LDA through double hybrid; grid through ``dft={...}``
   * - ``mp2``, ``scs-mp2``, ``sos-mp2``
     - ``basis``
     - Frozen core by default
   * - ``ri-mp2`` / ``df-mp2``
     - ``basis``, ``aux_basis``
     - The prefix is the request; the type is the same
   * - ``ccsd``, ``ccsd(t)``
     - ``basis``
     - Spin-adapted by default; ``cc={"spin_adapted": False}`` for spin orbitals
   * - ``ri-ccsd``, ``ri-ccsd(t)``
     - ``basis``, ``aux_basis``
     - Fitted ladder
   * - ``casscf``, ``casci``
     - ``basis``, ``mcscf``
     - Space by counts, by AVAS, or cut into ORMAS subspaces
   * - ``sapt0``
     - ``basis``, two monomers
     - Interaction energy, decomposed
   * - ``efp2``
     - a potential per monomer
     - No SCF: the wave functions were solved when the potentials were made

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
system and on the monomer.

.. code-block:: python

   oh = mqc.System(symbols=["O", "H"], coords=[[0, 0, 0], [0, 0, 0.97]],
                   multiplicity=2)
   oh.set_monomers([[0, 1]], multiplicities=[2])
   mqc.MBE(oh, level=0, method="hf", basis="sto-3g").run()

The keyword blocks
------------------

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
----------

``properties`` is the deck's block of the same name, beside ``keywords``
rather than inside it, and the distinction is worth keeping: keywords say how
to compute the wave function and change the number that comes out, properties
ask for something further to be done with it and change nothing. That is why a
bonding analysis leaves the driver as an energy:

.. code-block:: python

   mqc.MBE(water, level=0, method="hf", basis="6-31g",
           properties={"bonding_analysis": {"type": "gms_quao"}}).run()

See :doc:`bonding_analysis` for what the tables mean.

The energy decomposition
^^^^^^^^^^^^^^^^^^^^^^^^

``energy_decomposition`` resolves the total energy onto atoms and atom pairs
in the quasi-atomic basis, and ``result.bonding`` hands the numbers back:

.. code-block:: python

   >>> result = mqc.MBE(water, level=0, method="hf", basis="6-31g",
   ...     properties={"bonding_analysis": {"type": "gms_quao",
   ...                                      "energy_decomposition": True}}
   ...     ).run(label="water", write_to_file=True)
   >>> result.bonding["energy_of_formation"]
   -0.20800323755355832
   >>> result.bonding["atoms"][0]
   {'index': 0, 'energy': -73.655036293688, 'free_atom': -74.78030989551682,
    'adaptation': 1.1252736018288232}
   >>> result.bonding["pairs"][0]
   {'atoms': [0, 1], 'energy': -1.186123764253729,
    'classical': -0.3292398999705206, 'interference': -0.8568838642832084}

Atoms are 0-based; each pair is written once, with ``i < j``, carrying the full
pair energy. Everything is Hartree -- the printed tables convert to
millihartree and kcal/mol, this does not.

**For a single determinant the pieces add up exactly.** Every term of the
decomposition belongs to one atom or to one pair, so

.. code-block:: python

   sum(a["energy"] for a in b["atoms"]) + sum(p["energy"] for p in b["pairs"])

is ``result.energy`` to rounding -- it is a regrouping of that sum, not a model
of it. ``python/examples/methods.py`` asserts exactly that, and it is a
stronger check than any external reference: another code's number would say the
analysis agrees with somebody, this says it did not lose or double-count a term
of its own energy.

For a correlated wave function the two differ, and that is a statement rather
than an error: the gap is the part of the density living outside the
quasi-atomic span being summed over. A CAS(4,4) on water in 6-31G leaves
23 millihartree there, and the difference is worth reading as "how much of this
wave function the decomposition describes".

Two costs to know about. ``energy_decomposition`` is off by default because it
needs the dense ``n_ao^4`` integral array -- eight hundred megabytes at a
hundred basis functions -- where the bonding tables alone need one-electron
integrals; ``result.bonding`` being ``None`` is how you can tell it did not
run. And ``no_sharing: True`` (which implies the decomposition) solves a full
valence CI over the quasi-atomic orbitals: ethane is eleven million
determinants and benzene is out of reach.

.. _Interaction energies:

Interaction energies
--------------------

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
   >>> sorted(result.sapt)
   ['delta_hf', 'disp20', 'e_int_hf_cp', 'elst10', 'exch10', 'exch10_s2',
    'exch_disp20', 'exch_ind20_r', 'exch_ind20_u', 'ind20_r', 'ind20_u', 'total']

The names are the ones another program quotes its terms in, which is the point
of keeping them: a term can only be checked against the same term. A cluster is
one SAPT calculation per pair -- see :doc:`sapt`.

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

Every monomer needs one. A mixed quantum/EFP system is a deck-only feature, and
a system with no potentials at all is refused rather than returning the zero it
used to. See :doc:`makefp`.

Density fitting
---------------

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

Anything without a named argument
---------------------------------

``keywords`` merges straight into the deck's ``keywords`` block, and wins over
everything above, so a setting does not have to wait for a Python argument to
be reachable -- including whole features, like the generalized expansion:

.. code-block:: python

   mqc.MBE(water, level=0, method="ccsd(t)", basis="cc-pvdz",
           keywords={"cc": {"tolerance": 1e-10},
                     "correlation": {"n_frozen_core": 1},
                     "scf": {"guess": "sad"}}).run()

   # GMBE over overlapping fragments: the fragmentation block is a block like
   # any other, so `method` and its settings can be replaced outright.
   mqc.MBE(cluster, level=1, method="hf", basis="sto-3g",
           keywords={"fragmentation": {"method": "GMBE",
                                       "allow_overlapping_fragments": True}}).run()

``system_options`` does the same for the deck's ``system`` block -- it is not
called ``system`` because that is the molecule, the first argument.

Choosing which terms to compute
===============================

An expansion's term list is a value you can look at and change before anything
is computed. This is the part of the interface a deck has no equivalent for.

.. code-block:: python

   >>> mbe = mqc.MBE(cluster, level=2, method="hf", basis="sto-3g")
   >>> mbe.terms()
   [(1,), (2,), (1, 2), (3,), (1, 3), (2, 3)]

   >>> mbe.keep(lambda t: len(t) == 1 or t == (1, 2))
   >>> mbe.terms()
   [(1,), (2,), (1, 2), (3,)]

The predicate sees one term at a time, as a tuple of 1-based monomer indices,
so ``len(t)`` is its level. **The list is then closed under subsets, and that
is not a convenience.** An n-body term's delta is its energy less the delta of
every proper subset, so a screen that keeps a trimer and drops one of its
dimers has not approximated the expansion -- it has made one that cannot be
assembled. Dropped subsets are put back, so the list that runs may be longer
than the one the predicate chose, and ``run`` refuses a list that is not
closed rather than assembling nonsense.

``set_terms(...)`` hands a list over outright, for a restart or a screen done
elsewhere; it closes the list the same way.

Distance screening does not need any of that -- ``cutoffs`` is a per-level
threshold applied before the terms are generated:

.. code-block:: python

   mqc.MBE(cluster, level=3, method="gfn2",
           cutoffs={"dimer": 5.0, "trimer": 4.0}).run(label="w20")

Reading the answer
==================

``run`` returns a ``Result``. The energy is on it; everything else is read back
out of the files the run wrote, so ``write_to_file=True`` is what makes the
rest available.

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Attribute
     - What it is
   * - ``energy``
     - The total, in Hartree. For SAPT0 and EFP2, the interaction.
   * - ``fingerprint``
     - Identity of the calculation: geometry, partition, bonds, method,
       thresholds -- and deliberately not the logging or the rank count. The
       thing to compare before reusing any of these energies.
   * - ``gap_ev``
     - HOMO-LUMO gap of the whole system, unfragmented runs only. A fragmented
       run returns ``None`` on purpose: no sum of fragment gaps is the
       molecule's gap.
   * - ``gradient_norm``
     - Norm of the nuclear gradient, from a ``driver="gradient"`` run.
   * - ``sapt``
     - The SAPT0 decomposition as a dict, or ``None``.
   * - ``bonding``
     - The intrinsic energy decomposition -- formation energy, per-atom and
       per-pair terms -- from a run that asked for one, or ``None``.
   * - ``breakdown()``
     - The per-term rows, as ``Term`` objects.
   * - ``unconverged()``
     - The terms whose SCF did not converge. Empty is the answer you want.

.. code-block:: python

   >>> result = mqc.MBE(cluster, level=2, method="hf", basis="sto-3g").run(label="w3")
   >>> result.energy
   -224.86311196000983
   >>> result.fingerprint
   '6bd1d2d41af0172c'
   >>> row = result.breakdown()[0]
   >>> row
   <Term (1, 2) delta=+3.444e-03>
   >>> row.level, row.distance, row.converged, round(row.gap_ev, 3)
   (2, 2.9, True, 26.101)

A ``Term`` carries ``monomers`` (1-based), ``energy``, ``delta`` -- the n-body
contribution, which is what a threshold is on -- ``distance``, ``converged``,
and ``homo``/``lumo``. ``converged`` is a tri-state: ``None`` means the method
did not report, which is not a claim that it converged.

Per-fragment gaps are worth more than they look. Gaps are not additive, so no
sum of them is a molecular property, but the tail of that column is where an
SCF struggles and where a single determinant is the weakest description --
which tends to be the same rows as the ``converged is False`` ones.

Checkpoints and restart
=======================

``checkpoint=`` appends each fragment as it finishes and resumes from what is
already there. A missing file is created, an existing one is resumed, and one
written by a *different* calculation is refused rather than reused -- that is
what the fingerprint is for.

.. code-block:: python

   mbe = mqc.MBE(cluster, level=3, method="hf", basis="cc-pvdz",
                 checkpoint="w20.json")
   mbe.run(label="w20")     # dies at term 4000 of 6000
   mbe.run(label="w20")     # picks up at 4001

The text format holds one energy per record and is crash-proof line by line. A
run that needs derivatives has to have HDF5, and naming the file ``.h5`` asks
for it outright -- which is what a large screening pass wants. A build without
``MQC_ENABLE_HDF5`` says so rather than silently writing text.

Recipes
=======

Screen cheap, refine expensive
------------------------------

The two-pass calculation the interface exists for: run the whole expansion with
something cheap, keep the terms that actually contributed, and recompute those
with something good.

.. code-block:: python

   with mqc.session():
       cluster = mqc.System.from_xyz("w20.xyz")
       cluster.auto_monomers()

       cheap = mqc.MBE(cluster, level=3, method="gfn2")
       screen = cheap.run(label="screen")

       big = {t.monomers for t in screen.breakdown() if abs(t.delta) > 1e-5}
       good = mqc.MBE(cluster, level=3, method="dft", functional="pbe0",
                      basis="def2-svp")
       good.keep(lambda t: len(t) == 1 or t in big)
       print(good.run(label="refined").energy)

The criterion is energetic rather than geometric, which is the point: a
distance cutoff discards a strong interaction that happens to be far away.
``python/examples/energy_screened_mbe.py`` is this, complete.

Let the chemistry choose the fragments
--------------------------------------

.. code-block:: python

   cluster.compute_bond_orders(variant="gfn2")
   cut_here = [(i, j) for i, j in candidate_bonds
               if cluster.bond_order(i, j) < 0.5]

One xTB single point, then as many trial partitions as you like off the same
matrix.

Rank the bonds by what actually holds them
------------------------------------------

The decomposition gives every pair an energy and splits it into the part an
electrostatic model could produce and the part it could not. Sorting on the
second is asking which contacts are covalent rather than which are close:

.. code-block:: python

   result = mqc.MBE(molecule, level=0, method="hf", basis="6-31g",
                    properties={"bonding_analysis": {"type": "gms_quao",
                                                     "energy_decomposition": True}}
                    ).run(label="bonds", write_to_file=True)

   pairs = sorted(result.bonding["pairs"], key=lambda p: p["interference"])
   for pair in pairs[:5]:
       i, j = pair["atoms"]
       print(symbols[i], symbols[j], pair["energy"], pair["interference"])

The atom entries answer the other half of the question -- ``adaptation`` is
what it cost to prepare each atom in the shape the molecule needs, and binding
is the pairs giving back more than the atoms paid.

SAPT0 over every pair of a cluster
----------------------------------

SAPT is a two-body theory, so a cluster is one calculation per pair rather than
one calculation:

.. code-block:: python

   import itertools

   monomers = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
   for a, b in itertools.combinations(range(len(monomers)), 2):
       pair = mqc.System.from_xyz("trimer.xyz")
       pair.set_monomers([monomers[a], monomers[b]])
       out = mqc.MBE(pair, level=0, method="sapt0", basis="jun-cc-pvdz").run(
           label=f"pair_{a}_{b}", write_to_file=True)
       print(a, b, out.energy, out.sapt["elst10"])

Note the labels: one file set per pair, and no dots in the name.

Examples
========

``python/examples/`` holds runnable scripts, not fragments. All of them assert
and exit non-zero on a wrong number, so they run in CI rather than sitting
there going quietly stale:

``backends.py``
   Standalone and fragmented, xTB and Hartree-Fock -- the four paths through
   the program that the interface can reach. The Hartree-Fock references are
   PySCF; the fragmented cases check the expansion against the same dimer
   computed whole.

``methods.py``
   Every method the backends implement, checked against PySCF fed this
   repository's own basis data -- Hartree-Fock, unrestricted, Kohn-Sham through
   the double hybrids, the MP2 family, coupled cluster with and without the
   triples and with and without fitting, CASCI, CASSCF, ORMAS, AVAS, SAPT0 and
   EFP2. Where an identity says more than a reference would, the identity is
   checked instead: RI against conventional, spin-adapted against spin-orbital,
   CASSCF against CASCI on the same space. The last cases are the refusals --
   an unknown method, ``driver="optimize"``, EFP2 with no potentials -- which
   are there because each of them used to return a plausible number or take the
   interpreter down. The energy decomposition is checked against itself: its
   atoms and pairs must add back up to the energy the run reported.

``energy_screened_mbe.py``
   The two-pass calculation above: GFN2 over every term, DFT over the ones
   whose contribution cleared a threshold, assembled as
   ``E_low(everything) + sum of (E_high - E_low)`` over the kept terms.

``screening_quality.py``, ``screening_quality_bonded.py``
   How good that screen actually is, measured against the answer it
   approximates: the low pass, the *full* high-level expansion as the
   reference, and then the screened total for a range of thresholds --
   arithmetic over the two breakdowns, so the sweep costs nothing once the
   reference exists. The second script is the covalent twin, where the
   partition cuts bonds and the term list means something different.

``charges.py``
   Mulliken against CHELPG, and where they disagree.

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
   >>> ch3.set_monomers([[0, 1, 2, 3]], multiplicities=[2])
   >>> mqc.MBE(ch3, level=0, method="dft", functional="pbe",
   ...         basis="cc-pvdz").run(write_to_file=False).energy
   -39.7691396276

Nothing there asks for the unrestricted path: the multiplicity selects it. The
range-separated hybrids work on an open shell as well, with the attenuated
exchange pass run once per spin.

When something goes wrong
=========================

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - What you see
     - What it is
   * - ``ImportError: could not load libmqc``
     - The shared library was not built (``--target mqc_shared``) or is
       somewhere the search does not look. ``MQC_LIBRARY`` overrides it.
   * - The job hangs and never exits
     - The session was not closed. Use ``with mqc.session()``; the other ranks
       are waiting on a broadcast that will now never come.
   * - ``the monomers cut N bond(s) that were never declared``
     - The partition cuts covalent bonds nobody mentioned, so those fragments
       would be uncapped radicals. Declare them with ``set_bonds``, or
       ``perceive_bonds()``, or pass ``unchecked=True`` if you mean it.
   * - ``MQCError: unknown model.method '...'``
     - A typo, or a method this program does not have. The message lists what
       it accepts.
   * - ``... needs the tblite library``, ``... needs the cuEST integral
       backend``
     - The method exists; this build does not carry its backend. Reconfigure
       with the CMake option named in the message.
   * - ``label '...' would write its output as '...'``
     - The label reached Fortran as a filename, which eats dots and slashes.
       Two runs would have written over each other.
   * - ``the screen kept nothing``
     - A ``keep`` predicate that matched no term. There is nothing to compute.
   * - ``SCF not converged in N cycles``
     - Raise ``scf={"maxiter": ...}``, change ``guess={"type": ...}``, or pass
       ``allow_crap_scf=True`` if a bad number beats no number -- and then read
       ``result.unconverged()`` to see which terms it applied to.

API reference
=============

Session
-------

.. list-table::
   :widths: 34 66

   * - ``mqc.session()``
     - Context manager. Starts MPI on entry, releases the workers and shuts it
       down on exit, including on an exception.
   * - ``mqc.begin()`` / ``mqc.end()``
     - The same pair, unmanaged, for a REPL.
   * - ``mqc.rank()`` / ``mqc.n_ranks()``
     - This rank, and how many there are. Rank 0 is the only one that returns
       from ``begin()``.
   * - ``mqc.ELEMENTS``
     - Element symbols, H through Og, indexed from zero.
   * - ``mqc.MQCError``
     - Raised for anything Fortran refused, with its message.

``mqc.System``
--------------

.. list-table::
   :widths: 46 54

   * - ``System(symbols=, coords=, numbers=, charge=0, multiplicity=1)``
     - Angstrom, 0-based atoms.
   * - ``System.from_xyz(path, charge=0, multiplicity=1)``
     - Read an .xyz file.
   * - ``set_geometry(...)``
     - Replace the atoms of an existing system.
   * - ``auto_monomers(tolerance=1.2)``
     - Connected components as monomers.
   * - ``set_monomers(monomers, charges=None, multiplicities=None)``
     - The partition, as lists of 0-based atom indices.
   * - ``set_bonds(bonds, orders=None)``
     - Connectivity as ``(i, j)`` pairs. Broken is derived from the partition.
   * - ``perceive_bonds(tolerance=1.2)``
     - Fill the bond list in from the geometry, and declare it.
   * - ``missing_bonds(tolerance=1.2)``
     - Cross-monomer bonds the geometry implies and you did not declare.
   * - ``compute_bond_orders(variant="gfn2", accuracy=0.0)``
     - One xTB single point, kept on the system.
   * - ``bond_orders()`` / ``bond_order(i, j)`` / ``has_bond_orders``
     - The Wiberg-Mayer matrix, one pair, whether it exists.
   * - ``compute_charges(scheme="chelpg", basis="6-31g")``
     - One RHF, kept on the system. Closed shell only.
   * - ``charges()`` / ``charge_on(i)`` / ``charge_scheme`` / ``has_charges``
     - The partial charges, one atom, which scheme produced them.
   * - ``set_fragment_potentials(paths)``
     - One ``.efp`` file per monomer, for ``method="efp2"``.
   * - ``n_atoms`` / ``n_monomers`` / ``n_bonds`` / ``bonds_declared``
     - Counts, and whether bonds were stated at all.

``mqc.MBE``
-----------

.. list-table::
   :widths: 46 54

   * - ``MBE(system, level=2, method="gfn2", ...)``
     - See the argument table above.
   * - ``settings()``
     - The JSON document that will be sent.
   * - ``terms()``
     - The terms this expansion will compute, as 1-based tuples.
   * - ``keep(predicate)``
     - Screen the term list, then close it under subsets.
   * - ``set_terms(terms)``
     - Hand a list over outright; also closed under subsets.
   * - ``run(label="mqc", write_to_file=True)``
     - Compute across the whole job, and return a ``Result``.

``mqc.Result`` and ``mqc.Term``
-------------------------------

``Result``: ``energy``, ``label``, ``wrote``, ``fingerprint``, ``gap_ev``,
``gradient_norm``, ``sapt``, ``bonding``, ``breakdown()``, ``breakdown_csv``,
``unconverged()``.

``Term``: ``monomers``, ``energy``, ``delta``, ``distance``, ``converged``,
``homo``, ``lumo``, ``gap_ev``, ``level``.
