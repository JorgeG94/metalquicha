Fragment Molecular Orbital method (FMO and EE-MBE)
==================================================

Two fragmentation methods that solve each fragment **in the electrostatic field of
the others**, rather than in vacuum. Both are selected in the deck by name and
both work on non-covalently bonded systems -- clusters, solvated molecules,
anything where a fragment is a whole molecule.

.. code-block:: json

   "keywords": {
     "fragmentation": {
       "method": "mbe",
       "level": 2,
       "expansion": "fmo"
     }
   }

``"expansion"`` takes ``"fmo"``, ``"ee-mbe"``, or is left out for the ordinary
many-body expansion. ``"level"`` is how many fragments at a time, so ``2`` is
FMO2, ``3`` is FMO3, and so on.

Fragments come from the ``fragments`` list the deck already declares; nothing
extra needs saying to describe them.

Why bother
----------

An ordinary many-body expansion computes each fragment in vacuum and recovers
everything else through its correction terms. FMO computes each fragment in the
field the others make, iterates that field to self-consistency, and only then
applies corrections. For hydrogen-bonded systems the difference is large -- the
monomers of a water cluster are substantially polarized, and an expansion
starting from unpolarized monomers has a great deal of ground to make up.

Measured on a water trimer against a supermolecular RHF in the same basis:

============================  ==================
Method                        Error (Hartree)
============================  ==================
FMO2                          4.4e-05
EE-MBE                        4.0e-05
Plain MBE (no embedding)      7.6e-04
============================  ==================

The embedding is worth a factor of about twenty, and that is the whole argument
for these methods over a plain expansion.

The two methods
---------------

They share all their machinery and differ in two choices.

**FMO** (``"expansion": "fmo"``) builds the field from the neighbours' actual
electron densities: nuclear attraction integrals plus the Coulomb operator of
their density matrices. It then sums *internal* energies -- each fragment's own
energy with its polarized density, not counting its interaction with the field --
and adds a term for how each n-mer's density responds to that field.

**EE-MBE** (``"expansion": "ee-mbe"``), electrostatically embedded MBE,
represents each neighbour by atomic point charges instead, and sums the *total*
embedded energies in an ordinary many-body expansion with no response term.

These are genuinely different quantities, not two spellings of one. On water
trimers they land within about 15% of each other's error, so confusing them will
not produce anything obviously wrong -- which is the reason to be clear about
which one is running.

Level, and what it costs
------------------------

``level`` truncates the expansion. Level ``n`` on ``N`` fragments computes every
combination of ``n`` fragments, so the count is ``C(N,n)``:

=========  =====================  =====================
Level      20 fragments           50 fragments
=========  =====================  =====================
2          190 n-mers             1225 n-mers
3          1140 n-mers            19600 n-mers
4          4845 n-mers            230300 n-mers
=========  =====================  =====================

Nothing refuses a high level -- the expansion is generic -- but the binomial is
the whole story and it is not gentle. Level 2 is the usual choice; level 3 is
occasionally worth it; beyond that is a research question rather than a
calculation.

Climbing the level converges on the exact answer, and when the level reaches the
number of fragments it *is* the exact answer -- the corrections telescope to the
supermolecular energy. Measured on stacked waters in STO-3G:

=================  ================  ================
Level              3 waters          4 waters
=================  ================  ================
2                  1.0e-05           2.3e-05
3                  **1.4e-13**       1.3e-06
4                  --                **8.0e-13**
=================  ================  ================

The bold entries are where the level equals the fragment count. Those are not
approximations and the agreement is SCF convergence, not chemistry.

Tuning the field
----------------

All optional, all under ``keywords.fragmentation``.

``resppc`` (default ``2.0``)
   How far away a fragment has to be before it is represented by point charges
   instead of the exact Coulomb operator. This is what makes the method scale:
   the expensive term is then needed only within a neighbourhood, so the cost per
   fragment stops growing once the system is bigger than one.

   The separation is measured as FMO measures it -- the smallest interatomic
   distance between two fragments, divided by the sum of those two atoms' van der
   Waals radii -- so it is a contact distance rather than a centre-to-centre one.
   The default matches GAMESS's ``RESPPC``.

   A negative value turns the approximation off and makes every neighbour exact.

``far_field`` (default ``"mulliken"``)
   What a distant fragment contributes: ``"mulliken"``, ``"chelpg"``, or
   ``"ignore"`` for nothing at all.

   Mulliken is the default and is what production FMO codes use. The
   approximation being made is a population one, and Mulliken populations are
   what the term being approximated reduces to.

   ``"ignore"`` is not really an approximation to the field so much as a decision
   not to have one past the cutoff. It is the honest way to ask what the
   long-range field is worth: set the cutoff where you mean to and compare.

``max_outer`` (default ``50``) and ``outer_tolerance`` (default ``1e-7``)
   Cap and convergence for the monomer self-consistency loop, the latter on the
   sum of monomer energies in Hartree.

What the cutoff costs
~~~~~~~~~~~~~~~~~~~~~

Against the same calculation with the approximation switched off, five stacked
waters:

===========  ===============  ===============
O--O sep     ``resppc`` 2.0   No cutoff
===========  ===============  ===============
2.70 A       -4.90e-04        -4.93e-04
3.20 A       1.23e-05         -1.03e-06
9.00 A       8.07e-08         1.71e-13
===========  ===============  ===============

As ratios that looks alarming. As energies it is the approximation working: where
the cutoff costs most in relative terms the absolute error it leaves is under
1e-07 Hartree, and where the error is large enough to matter the cutoff has not
engaged at all. Precision is given up where precision is not the binding
constraint.

Running it
----------

.. code-block:: bash

   ./mqc fmo_water3.json
   mpirun -np 4 ./mqc fmo_water3.json

Both phases distribute. A monomer pass is a set of independent tasks followed by
a barrier and a density exchange; the n-mers are independent of each other and of
everything else, and being the larger count are the part worth spreading. Every
rank holds the whole geometry and assembles only the fragments it was given, so
no fragment geometry is ever sent -- what crosses is the densities, energies and
charges the next pass needs, all of them small.

The answer does not depend on the rank count. That is asserted rather than
assumed: ``validation/check_fmo_mpi`` runs each method with a communicator and
again without one and compares, and they agree to about 1e-12 on one, two, three
and four ranks.

Example decks are in ``validation/inputs/cpu/mqc/fmo/``.

Limits
------

**Non-covalent fragments only, and this is enforced.** There is no capping here --
no adjusted fragment orbitals, no hybrid orbital projection, no hydrogen caps --
so a partition that severs a covalent bond is refused rather than answered, with
a message naming the two atoms and the two fragments they were put in.

The check is not a formality. Cutting a single bond leaves both fragments with an
odd electron count, which the closed-shell check catches on its own; but cutting
an even number per fragment -- a ring, a double bond -- leaves every count even.
Cyclopropane split into three CH2 groups used to come back 0.28 Hartree low,
which is 176 kcal/mol in the shape of an answer.

A hydrogen bond is not a covalent one, so clusters are unaffected: two waters
2.5 A apart, closer than anything in the validation set, pass through.

**Closed shell only.** A fragment with an odd electron count is refused rather
than quietly paired up.

**Energies only.** No gradients yet, so geometry optimization and frequencies are
not available through these.
