Atomic charges and bond orders
==============================

Two properties exposed through the **Python interface**, and both also through a
deck. The Python surface is deliberate for both: neither is the output
of a calculation somebody asked for, they are inputs to deciding what
calculation to run. Working out where to cut a molecule means asking about the
same system many times over many trial partitions, which is a loop, not a JSON
key.

Both follow the same shape there -- compute once onto a system handle, read many
times.

Both have a second surface because they answer a second question. Once the
calculation you wanted has run, its charges and its bond orders are properties
*of that calculation*, and reading them off the density it already converged
costs no second SCF. Those are ``properties.charges`` and
``properties.bond_orders``, below, and they are the ones to reach for unless you
are in the trial-partition loop.

Three things called bond orders
-------------------------------

Before any of the numbers: this code computes three different quantities under
that name, they do not agree beyond a trend, and none of them is an
approximation to either of the others.

======================  ================================  ================================
Quantity                Where                             What it is
======================  ================================  ================================
**Wiberg--Mayer, xTB**  ``compute_bond_orders("gfn2")``   Mayer's definition inside a
                                                          semi-empirical Hamiltonian, in
                                                          its own minimal basis
**Mayer**               ``properties.bond_orders``, or    The same definition over the
                        ``compute_bond_orders("mayer")``  converged ab initio density and
                                                          the AO overlap
**QUAO kinetic**        ``properties.bonding_analysis``   Ruedenberg's definition in the
                                                          orthonormal quasi-atomic basis,
                                                          plus an energy-like partner
======================  ================================  ================================

Which to use. The xTB orders are for **deciding where to cut a molecule**: they
cost one semi-empirical single point, so a loop over trial partitions can afford
them. The Mayer orders are for **checking that decision**, and for reporting the
bonding of a calculation you were running anyway. The QUAO analysis is for
**taking a bonding picture apart** -- it comes with orbitals, an energy
decomposition and a table of which quasi-atomic orbitals are involved.

Wiberg's bond order, in the strict sense, is this same sum of squares in an
*orthonormal* basis, where ``S`` is the identity and the ``D S`` product
collapses. Over non-orthogonal atomic orbitals it is not Wiberg's quantity and
not anyone else's, so ``properties.bond_orders`` refuses ``"wiberg"`` by name
rather than quietly computing Mayer's. Over the quasi-atomic basis, which *is*
orthonormal, the QUAO analysis's population bond order already is it.

Mayer bond orders
-----------------

.. math::

   B_{AB} = \sum_{\mu \in A} \sum_{\nu \in B} (DS)_{\mu\nu} (DS)_{\nu\mu}

for a closed shell, with :math:`D` the total density. For an open shell it is
**not** that expression with the total density; it is

.. math::

   B_{AB} = 2 \sum_{\mu \in A} \sum_{\nu \in B}
            \left[ (D_\alpha S)_{\mu\nu} (D_\alpha S)_{\nu\mu}
                  + (D_\beta S)_{\mu\nu} (D_\beta S)_{\nu\mu} \right]

which reduces to the first when :math:`D_\alpha = D_\beta = D/2`. The
difference is silent on every closed shell and large on an open one -- triplet
O\ :sub:`2` in STO-3G is 2.00 by the right formula and 1.50 by the wrong one --
so the unrestricted case has its own test rather than an assumption.

From a deck
~~~~~~~~~~~

.. code-block:: json

   {
     "model": {"method": "hf", "basis": "6-31g"},
     "driver": "Energy",
     "molecules": [{"xyz": "ethane.xyz",
                    "molecular_charge": 0, "molecular_multiplicity": 1}],
     "properties": {"bond_orders": {"scheme": "mayer"}}
   }

The **object** is the request and ``scheme`` only says which definition, so
``"bond_orders": {}`` is a valid ask and takes Mayer -- the same shape
``properties.charges`` has. Whatever reference converged is what gets
partitioned: Hartree--Fock or Kohn--Sham, restricted or unrestricted, at no cost
beyond the SCF that was going to run anyway.

Unlike the charges, the table **is** printed, and the JSON carries the whole
matrix plus the per-atom valence :math:`V_A = \sum_{B \neq A} B_{AB}`:

.. code-block:: json

   "bond_orders": {
     "scheme": "mayer",
     "matrix": [[0.0, 0.9237, ...], ...],
     "atoms": [{"atom": 1, "valence": 3.7619}, ...]
   }

The matrix is written whole rather than as a list of bonded pairs, because where
the line between a weak bond and none falls is the reader's question and a
threshold applied here would answer it silently. The diagonal is zero: an atom
is not bonded to itself, and what the block sum would put there is a different
quantity.

**Numbers to expect**, restricted Hartree--Fock:

===============  =======  ========  ========  ========
Molecule         Basis    Bond      Order     Valence
===============  =======  ========  ========  ========
Ethane           STO-3G   C--C      1.011     C 3.968
Ethane           STO-3G   C--H      0.984     H 0.997
Ethane           6-31G    C--C      0.924     C 3.762
Ethane           6-31G    C--H      0.961     H 0.931
Benzene          6-31G    C--C      1.444     C 3.857
Benzene          6-31G    C--H      0.943     H 0.943
Water            STO-3G   O--H      0.954     O 1.908
Water            6-31G    O--H      0.803     O 1.607
Triplet O2       STO-3G   O--O      2.000     O 2.000
===============  =======  ========  ========  ========

Benzene's 1.44 is the delocalised bond-and-a-half, and its para pair -- carbons
across the ring, bonded to nothing -- comes out at 0.10, which is the honest
answer for a conjugated system rather than a zero. None of that comes from a
distance criterion; there is no distance anywhere in the definition.

Small entries can be **negative**, and that is a property of the definition and
not a bug: ethane's carbon to a hydrogen on the other carbon is +0.002 in
STO-3G and -0.015 in 6-31G. :math:`B_{AB}` is a sum of products over a
non-orthogonal basis with no positivity to it, so anything near zero is near
zero from either side. Read the sign of a hundredth as noise; a bond order is a
count of shared pairs only where there is something to share.

The basis dependence is real but mild next to a Mulliken charge's: ethane's
C--C moves 9% between STO-3G and 6-31G, where the Mulliken charge on water's
oxygen moves by a factor of two over a comparable change.

A **multiconfigurational** wave function is refused rather than skipped, for the
same reason the charges are: the 1-RDM there is in the MO basis over fractional
occupations and the AO density this needs is not formed. On the **GPU backend**
the request is refused too -- the arithmetic is backend-independent but the spin
convention is not, and an unverified factor of two in the open-shell case is
worse than a message saying to run it on the CPU.

Wiberg--Mayer bond orders from xTB
----------------------------------

Over the whole system, from one semi-empirical single point.

.. code-block:: python

   import mqc

   s = mqc.System.from_xyz("cluster.xyz")
   s.compute_bond_orders()          # one xTB single point

   s.bond_orders()                  # the full matrix, as a list of rows
   s.bond_order(0, 1)               # one pair, 0-based
   s.has_bond_orders                # whether compute has run

``compute_bond_orders(variant="gfn2", accuracy=0.0)`` takes ``"gfn2"`` or
``"gfn1"``; an accuracy of zero or less uses tblite's default.
``variant="mayer"`` runs a real RHF in ``basis`` instead and takes Mayer's
orders off its density -- closed shell only, and it costs what an SCF costs:

.. code-block:: python

   s.compute_bond_orders(variant="mayer", basis="6-31g")
   s.bond_order_scheme            # "mayer" -- ask before comparing two runs
   s.bond_order_valences()        # sum_B B_AB, per atom; the Mayer variant only

Both variants land on the same handle and are read through the same
``bond_orders()``, which is what makes the comparison a two-line script. They
are still different quantities: ``bond_order_scheme`` says which one is
currently there, and a valence is offered only for the Mayer variant because
summing the rows of the other would be a number nobody computed.

The whole system, not the monomers: the point of these is to decide where the
monomers should be, so a partition cannot be an input. A caller wanting
fragment-local orders builds a handle per fragment.

**What they are good for.** On cases with known answers, GFN2 gives 1.03, 2.03
and 3.00 for the C--C bonds of ethane, ethene and ethyne, and 0.019 for a water
dimer's hydrogen bond. A real single bond and a hydrogen bond are separated by a
factor of fifty, which is a distinction a covalent-radius rule cannot make at
all -- it sees only the distance and calls both bonds or neither.

**What they are not good for.** They do not rank cuts within a molecule. Decane's
nine C--C bonds span 1.3%, so bond order says nothing about which to break.
Treat them as a veto on unsafe cuts rather than a ranking of good ones.

Do the cheap orders rank the same?
----------------------------------

The reason both are exposed through one call. Formic acid, xTB (GFN2) against
Mayer over an RHF/6-31G density, same geometry, both read off the same handle:

=========  ==========  ==========
Pair       xTB         Mayer
=========  ==========  ==========
C--O (=O)  1.777       1.825
C--O (-O)  1.205       0.924
C--H       0.941       0.887
O--H       0.871       0.769
O...O      0.181       0.012
O...H      0.032       0.010
=========  ==========  ==========

**The four real bonds rank identically**, and in the order chemistry expects:
the carbonyl above the hydroxyl C--O, then C--H, then O--H. The magnitudes
disagree -- xTB puts the single C--O at 1.21 where Mayer puts it at 0.92 -- so
the two are not interchangeable as numbers. Over all ten pairs Spearman is
0.82; the disagreement is entirely among the non-bonded pairs, which both codes
place near zero in an order neither of them means anything by.

Ethane is the same story: the seven bonds are the top seven for both, but xTB
ranks C--C above C--H (1.031 against 0.989) where Mayer in 6-31G ranks it below
(0.924 against 0.961). Spearman over its 28 pairs is 0.93.

The conclusion to draw for fragmentation: the cheap orders are sound for
**separating bonds from non-bonds**, which is what a cut decision needs, and
they should not be trusted to order two bonds that are close together -- a
caution the xTB section already gives for a different reason.

Atomic charges
--------------

Mulliken and CHELPG partial charges. Through the Python interface, from an RHF
density; through a deck, from whatever reference the deck converged.

From the Python interface
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   s = mqc.System.from_xyz("water.xyz")
   s.compute_charges(scheme="chelpg", basis="6-31g")

   s.charges()                      # one per atom, input order
   s.charge_on(0)                   # one atom, 0-based
   s.charge_scheme                  # which scheme produced them
   s.has_charges

Unlike bond orders these cost a real SCF in the basis you name -- xTB is cheap
enough to point at anything, an RHF is not. Measured: 32 atoms in 6-31G is about
twenty seconds, nearly all of it the SCF. So the basis is the knob that matters
and the choice of scheme is not; CHELPG adds about ten percent to a calculation
that has to run anyway.

Closed shell only, and Hartree-Fock only. An odd electron count raises rather
than being quietly paired. This entry point builds its own molecule and runs its
own SCF from a geometry and a basis name, so there is nowhere for a functional
or a spin multiplicity to come from; a deck has both, which is why the deck
surface below does not share the restriction.

From a deck
~~~~~~~~~~~

``properties.charges`` partitions the density the calculation already converged.

.. code-block:: json

   {
     "model": {"method": "dft", "basis": "6-31g", "functional": "b3lyp"},
     "driver": "Energy",
     "molecules": [{"xyz": "water.xyz",
                    "molecular_charge": 0, "molecular_multiplicity": 1}],
     "properties": {"charges": {"scheme": "mulliken"}}
   }

The **object** is the request and ``scheme`` only says how, so
``"charges": {}`` is a valid ask and takes Mulliken. Both partition routines
take a density matrix and neither knows what produced it, so Hartree-Fock and
Kohn-Sham, restricted and unrestricted, all work and all cost nothing beyond the
SCF that was going to run anyway.

Note the default differs from ``properties.fukui``, which takes CHELPG. A
condensed Fukui index is a difference of two charges, where Mulliken's
basis-set sensitivity does not cancel; asked for on its own a charge is usually
wanted as the cheap population number, and Mulliken is one trace against an
overlap that already exists.

Where the numbers come out
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Not on the terminal.** Neither scheme prints a table of charges however
verbose the logger is; they are written to the JSON output beside the deck, as
``output_<deck>.json``, under ``atomic_charges``. CHELPG prints a single line
saying how many grid points it fitted, which is a progress note and not the
result -- it is easy to read that line, see no table, and conclude the request
was dropped.

.. code-block:: console

   $ mqc water.json
   ...
     chelpg: fitted 3 charges to 4358 grid points
   $ python3 -c "import json; d=json.load(open('output_water.json'));
   >   k=list(d)[0]; print(json.dumps(d[k]['atomic_charges'], indent=2))"

.. code-block:: json

   {
     "scheme": "chelpg",
     "sum": -5.55e-17,
     "atoms": [
       {"atom": 1, "charge": -0.957464},
       {"atom": 2, "charge": 0.478851},
       {"atom": 3, "charge": 0.478613}
     ]
   }

``sum`` is the total charge the partition accounts for and should equal the
molecular charge; a value that does not is the fastest way to see that a fit
went wrong, which is why it is written rather than left to be added up.

The same water in the same basis gives -0.957 on the oxygen from CHELPG and
-0.800 from Mulliken. That spread is not an error in either -- they answer
different questions, one fitting the electrostatic potential and the other
dividing the overlap -- and it is the reason a scheme is named in the output
rather than assumed by whoever reads it.

For an unrestricted reference the charges come from ``P_alpha + P_beta``, and
Mulliken additionally reports **spin populations** from ``P_alpha - P_beta``:

.. code-block:: json

   "atomic_charges": {
     "scheme": "mulliken",
     "sum": 1.0,
     "atoms": [
       {"atom": 1, "charge": -0.123936, "spin_population": 1.104060},
       {"atom": 2, "charge": 0.561968, "spin_population": -0.052030},
       {"atom": 3, "charge": 0.561968, "spin_population": -0.052030}
     ]
   }

``sum`` is written out so a consumer can check rather than trust: charges sum to
the molecular charge, and spin populations to ``n_alpha - n_beta``. There is no
CHELPG spin analogue, and its absence is meaningful rather than a gap -- that
scheme fits the electrostatic potential, which the total density alone
determines.

Two things it will not do. A **multiconfigurational** wave function is refused
rather than skipped, because the 1-RDM there is in the MO basis over orbitals
with fractional occupation and the AO density both schemes want is not formed.
And on a **fragmented** run the charges are the fragment's own, hydrogen caps
included -- which is what makes them checkable, since the column then sums to
the charge of the molecule the SCF actually saw. Dropping the caps would leave
a column summing to nothing in particular.

Which scheme
~~~~~~~~~~~~

**Mulliken** splits the density by which basis function carries it, halving every
overlap between the two atoms it spans. Cheap, and notoriously basis-set
dependent.

**CHELPG** solves for the point charges that best reproduce the molecule's own
electrostatic potential on a shell of points outside its van der Waals surface,
constrained to sum to the molecular charge. That is a physically meaningful
question -- the potential is an observable of the density -- and the answer moves
far less when the basis changes.

The same water, two bases:

===========  ===============  ===============
Basis        Mulliken, O      CHELPG, O
===========  ===============  ===============
6-31G        -0.79            -0.94
aug-cc-pVDZ  -0.30            -0.74
===========  ===============  ===============

The molecule did not change. Mulliken moves by half an electron because a diffuse
function centred on hydrogen reaches well over the oxygen while still counting as
hydrogen's; CHELPG moves less than half as far, being fitted to a field rather
than to a basis. For anything downstream that treats a charge as a physical
quantity -- embedding especially -- CHELPG is the defensible input.

One caveat worth knowing: a nearly nonpolar molecule fits badly in *relative*
terms however good the code is. Methane has no dipole and no quadrupole, so its
potential outside the van der Waals surface is ten times smaller than water's and
is mostly charge penetration, which no arrangement of point charges can
reproduce. The charges themselves are still right (C -0.45, H +0.11); it is the
fit residual as a fraction of a tiny potential that looks poor.

Requirements
------------

Charges and Mayer bond orders need an integrals backend, since they need a
density. In a build with none the Python methods raise with a message naming the
CMake option rather than failing at import. The xTB bond orders need tblite
instead, and nothing else -- the two refusals are worded differently on purpose,
because they are asking for different things to be installed.

Both backends compute Mulliken charges, from the same code: the partition is a
trace of ``D S`` against an AO-to-atom map, and only the map is arrived at
differently, so the two cannot drift apart. **CHELPG is CPU only.** It fits the
electrostatic potential at points away from the atoms, and the GPU backend
builds no integral for that; asking for it there is refused by name rather than
quietly answered with Mulliken, because the two schemes disagree by design.

.. note::

   The GPU charge path is compiled and unit-tested but has not been run against
   cuEST, which needs an sm_80 card. The arithmetic it runs is the CPU path's,
   so what is unverified is the AO-to-atom map and the overlap it is handed --
   check a molecule with inequivalent atoms against the CPU backend before
   trusting a GPU charge.
