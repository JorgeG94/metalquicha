Atomic charges and bond orders
==============================

Two properties exposed through the **Python interface** rather than through a
deck. That is deliberate: neither is the output of a calculation somebody asked
for, they are inputs to deciding what calculation to run. Working out where to
cut a molecule means asking about the same system many times over many trial
partitions, which is a loop, not a JSON key.

Both follow the same shape -- compute once onto a system handle, read many times.

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

Mulliken and CHELPG partial charges, from an RHF density.

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

Closed shell only. An odd electron count raises rather than being quietly paired.

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
