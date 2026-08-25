Atomic charges and bond orders
==============================

Two properties exposed through the **Python interface**, and charges also
through a deck. The Python surface is deliberate for both: neither is the output
of a calculation somebody asked for, they are inputs to deciding what
calculation to run. Working out where to cut a molecule means asking about the
same system many times over many trial partitions, which is a loop, not a JSON
key.

Both follow the same shape there -- compute once onto a system handle, read many
times.

Charges have a second surface because they answer a second question. Once the
calculation you wanted has run, its charges are a property *of that
calculation*, and partitioning the density it already converged costs no second
SCF. That is ``properties.charges``, below, and it is the one to reach for
unless you are in the trial-partition loop.

Bond orders
-----------

Wiberg--Mayer bond orders from xTB, over the whole system.

.. code-block:: python

   import mqc

   s = mqc.System.from_xyz("cluster.xyz")
   s.compute_bond_orders()          # one xTB single point

   s.bond_orders()                  # the full matrix, as a list of rows
   s.bond_order(0, 1)               # one pair, 0-based
   s.has_bond_orders                # whether compute has run

``compute_bond_orders(variant="gfn2", accuracy=0.0)`` takes ``"gfn2"`` or
``"gfn1"``; an accuracy of zero or less uses tblite's default.

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

Charges need the libcint integrals backend, since they need a density. In a build
without it the Python methods raise with a message naming the CMake option rather
than failing at import. Bond orders need tblite.
