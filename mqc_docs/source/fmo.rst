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
supermolecular energy. That holds for both methods, and for the same reason: the
top n-mer holds every fragment, so there is nothing outside it to embed in.
Which method is running changes what a correction is made of, not that it
cancels.

Four stacked waters in STO-3G, error against the supermolecule in Hartree:

=========  ================  ================
Level      FMO_n             EE-MBE_n
=========  ================  ================
2          2.3e-05           2.0e-04
3          1.3e-06           4.1e-06
4          **6.3e-13**       **4.5e-13**
=========  ================  ================

The bold row is where the level equals the fragment count. Those are not
approximations and the agreement is SCF convergence, not chemistry.

Note the truncated rows: FMO's exact embedding is worth roughly an order of
magnitude over point charges at level 2, and the gap narrows as the level climbs
and the expansion itself does more of the work.

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

Cutting a covalent bond
-----------------------

``bond_breaking`` (default ``"none"``)
   ``"none"`` refuses a partition that severs a covalent bond, naming the two
   atoms and the two fragments they were put in. ``"afo"`` detaches the bond with
   an **adjusted frozen orbital** instead.

The refusal is not a formality, and it is still the default. Cutting a single
bond leaves both fragments with an odd electron count, which the closed-shell
check catches on its own; but cutting an even number per fragment -- a ring, a
double bond -- leaves every count even. Cyclopropane split into three CH2 groups
used to come back 0.28 Hartree low, which is 176 kcal/mol in the shape of an
answer.

A hydrogen bond is not a covalent one, so clusters are unaffected: two waters
2.5 A apart, closer than anything in the validation set, pass through and never
reach any of this.

How a bond is detached
~~~~~~~~~~~~~~~~~~~~~~

Each cut bond gets a small **model system** -- both its atoms, everything within
a radius of either, every hydrogen hanging off what that took, and a hydrogen cap
for each bond leaving the set at the standard length for the atom it hangs off.
That is solved, localized, and the orbital sitting on the bond is kept, reduced
to the coefficients on the bond-detached atom. Expressing it there is what makes
it transferable: those functions exist unchanged in any fragment containing that
atom, so putting it to work is an index map.

The two fragments then split the bond, following the assignment FMO uses. Of the
pair, one atom is the *detached* end and one the *attached* end:

============================  =========================  =========================
                              fragment of the detached   fragment of the attached
============================  =========================  =========================
nucleus                       owns it                    not owned
its basis functions           owns them                  carries them, ghosted
the bond's electron pair      none of it                 all of it
the hybrid on that atom       frozen empty               frozen occupied
electron count                ``sum(Z) - 1``             ``sum(Z) + 1``
============================  =========================  =========================

Frozen means the Fock matrix is forced block diagonal in a basis holding those
orbitals: the couplings between them and the variational space are zeroed and the
frozen virtuals held at a level shift. Zeroed rather than penalised, because a
penalty raises an orbital's energy without decoupling it and the two spaces go on
mixing however large the penalty is.

**The boundary set belongs to the n-mer, not to the fragment.** A bond cut
between two monomers is whole again inside the dimer holding both ends, so that
dimer carries no ghost, no frozen orbital and no electron shift there, while
still being cut against everything outside itself. This is worked out from each
group's own members every time.

What it costs
~~~~~~~~~~~~~

Propane in STO-3G, split at both C-C bonds into three fragments, against
ordinary RHF on the whole molecule:

=================================  ==================
Expansion                          Error, Hartree
=================================  ==================
Two fragments, one bond, MBE(2)    exact
Three fragments, MBE(3)            1.3e-13
Three fragments, MBE(2)            0.180
=================================  ==================

The middle row is the statement worth reading. An expansion carried to the
fragment count is exact by inclusion and exclusion whatever the partition did, so
landing on the whole molecule to 1e-13 says the bookkeeping is right across every
group -- three monomers with boundaries, three dimers, and one of those dimers
the pair of end fragments, which are not bonded to each other and whose group
carries a ghost of a carbon belonging to neither.

The last row is the three-body term, and across covalent bonds it is large.
Expect that: the same quantity is a rounding error for a water cluster and
110 kcal/mol here. Truncating at pairs is not advisable across detached bonds.

Restrictions
~~~~~~~~~~~~

``bond_breaking = "afo"`` requires ``embedding = "none"``. A frozen orbital and
an embedding field both describe the bond region -- the field already supplies
the neighbour's nucleus and density where the frozen orbital supplies the bond --
so the detached atom's share has to be removed from the field before the two can
be used together. That is clean for point charges and not defined for an exact
density, and neither is built yet. Asking for both is refused with that reason.

Refused by name, rather than answered badly:

* a cut through a ring, where two fragments meet in more than one place
* a bond detached at a hydrogen, which has nothing left to hybridise
* a bond carrying more than one localized orbital, which is not a single bond

Limits
------

**Closed shell only.** A fragment with an odd electron count is refused rather
than quietly paired up. A detached bond moves an electron between the two
fragments it joins, and the count checked here is the one after that move: ethane
split into two methyls is 9 and 9 before it and 8 and 10 after.

**Energies only.** No gradients yet, so geometry optimization and frequencies are
not available through these.
