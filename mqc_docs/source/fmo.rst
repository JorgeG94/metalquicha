Fragment Molecular Orbital method (FMO and EE-MBE)
==================================================

Two fragmentation methods that solve each fragment **in the electrostatic field of
the others**, rather than in vacuum. Both are selected in the deck by name and
both work on non-covalently bonded systems -- clusters, solvated molecules,
anything where a fragment is a whole molecule.

.. code-block:: json

   "keywords": {
     "fragmentation": {
       "method": "fmo",
       "level": 2
     }
   }

``"method"`` takes ``"fmo"``, ``"ee-mbe"``, or ``"mbe"`` for the ordinary
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

**FMO** (``"method": "fmo"``) builds the field from the neighbours' actual
electron densities: nuclear attraction integrals plus the Coulomb operator of
their density matrices. It then sums *internal* energies -- each fragment's own
energy with its polarized density, not counting its interaction with the field --
and adds a term for how each n-mer's density responds to that field.

**EE-MBE** (``"method": "ee-mbe"``), electrostatically embedded MBE,
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
and four ranks. It runs the three embeddings on a water tetramer and, since a
detached bond makes the exchanged density wider than the fragment's own basis,
propane split at both of its carbon-carbon bonds with point charges on top --
that case is bit-identical on one, two and four ranks.

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
a radius of either, every *singly bonded* atom hanging off what that took, and a
hydrogen cap for each bond leaving the set at the standard length for the atom it
hangs off.

A singly bonded neighbour comes in wholesale rather than by distance so that one
just outside the radius is not swapped for a cap hydrogen a few hundredths of an
Angstrom away, which is a discontinuity in anything that moves the geometry. For
a heavy one -- a carbonyl oxygen, above all -- there is a second reason: a cap
hydrogen closes exactly one electron pair, a carbonyl oxygen is held by two, and
capping it leaves the model a radical. Bond perception here is distance-based and
cannot report a bond order, so the order is not guessed and the atom comes in
whole instead. Without that rule every backbone cut of a peptide failed, amide
and C-alpha--C alike, because the sphere reaches a neighbouring carbonyl carbon
without reaching its oxygen.
That is solved, localized, and the orbital sitting on the bond is kept, reduced
to the coefficients on the bond-detached atom. Expressing it there is what makes
it transferable: those functions exist unchanged in any fragment containing that
atom, so putting it to work is an index map.

The two fragments then split the bond, following the assignment FMO uses. Of the
pair, one atom is the *detached* end and one the *attached* end:

============================  =========================  =========================
                              fragment of the detached   fragment of the attached
============================  =========================  =========================
nucleus                       ``Z - 1`` of it            ``+1`` of it, on a ghost
its basis functions           owns them                  carries them, ghosted
the bond's electron pair      none of it                 all of it
the hybrid on that atom       frozen empty               frozen occupied
electron count                ``sum(Z) - 1``             ``sum(Z) + 1``
net charge                    ``0``                      ``0``
============================  =========================  =========================

The nucleus is split because the electron pair is. One unit of charge crosses
the bond with the pair, so both fragments come out neutral closed shells and
the two halves add back to ``Z`` inside any n-mer holding the whole bond. That
is GAMESS's convention, and it is what a fragment *potential* needs: a fragment
carrying a unit charge puts a monopole term of order ``1/R`` on every
adjacent-residue pair of a protein.

**It happens only where there is a field**, because a split nucleus is only
defined when something supplies the other half. With ``embedding = "ptc"`` the
fragment across the bond supplies it -- it holds the ``+1`` and the bond pair,
and the group on this side feels both -- so the split is free, and the table
above is what runs. With ``embedding = "none"`` nothing supplies it, so the
nucleus stays whole with its owner and the two sides come out at about ``+1``
and ``-1``: the electron still moves, only the proton does not. Solving a
methyl group around a nucleus short by one proton is a worse model than the
cation, and the table under `What it costs`_ says by how much.

GAMESS splits unconditionally and that does not settle the question, because
GAMESS never runs this without a field at all -- its field-free reference state
is built with methyl caps, which is a third construction again.

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

A detached bond under a field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A detached atom is described by **two** fragments: its owner holds ``Z-1`` of
the nucleus with the hybrid there frozen empty, and the fragment across the bond
holds the same functions as a ghost carrying ``+1`` with the bond pair in that
hybrid. Three things follow, and all three are what ``embedding = "ptc"`` needed
before it could be allowed.

**Its population arrives twice, and the two shares are added.** Summing is the
only apportionment that leaves the atomic charges adding to the molecular charge,
and those charges are what the field is built from -- keep one share and every
fragment is embedded in a system carrying one elementary charge per cut that the
molecule does not have. Giving the whole atom to one side also conserves charge
but has to choose a side, and a Mulliken population is already an apportionment
of that kind.

**A group does not feel its own share.** The field is by definition the rest of
the system, so a group's members' contributions to a shared atom are taken back
out. What is left is not zero: the other fragment's share is outside and the
borrowed bond pair genuinely feels it. Dropping the whole atom from the field
instead -- the cheaper fix -- would take that term out too and leave every group
short by roughly an elementary charge per cut.

**Its density block travels with it.** A monomer density is kept at the size its
SCF produced, ghost block and all, and the members of an n-mer are laid out atom
by atom rather than as contiguous corners, because inside a group holding both
ends of a cut the borrowed block belongs in the other member's run.

**And the split nucleus makes no difference to any of this.** The charge a
group feels at an atom is ``q_all - own_q``, its own share taken back out, so a
unit of charge moved out of a group's own nucleus and into the field it sees
arrives at the same point with the same magnitude: the total one-electron
potential there is ``-(Z - pop_outside)/r`` either way. ``q_all`` does not move
either, because both Mulliken shares of a detached atom are summed into it. And
nuclear repulsion with per-atom charges that add up over the fragments is
pairwise-additive *in the fragments*, so an expansion reproduces it exactly at
level two whatever the assignment. **The embedded energy is therefore invariant
to the charge convention, exactly**, and measured that way -- propane's FMO(2)
error is 0.48908757203956554 with whole nuclei and 0.48908757203680864 with
split ones, which is convergence noise. The two runs are in the unit suite side
by side, as ``the_charge_convention_does_not_move_an_embedded_total``.

So the reason to split the nucleus is not the FMO energy at all. It is that a
fragment *potential* -- an EFMO far pair, an effective fragment -- is an
expansion about the fragment's own charge distribution, with no own-share
subtraction anywhere to undo a net charge. A unit charge per fragment is then a
monopole term of order ``1/R`` on every adjacent-residue pair of a protein,
about ninety kcal/mol where the interaction of interest is single digits, and
it lands in the electrostatics column.

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
Three fragments, FMO(3), ``ptc``   1.3e-13
Three fragments, FMO(2), ``ptc``   0.489
=================================  ==================

The middle row is the statement worth reading. An expansion carried to the
fragment count is exact by inclusion and exclusion whatever the partition did, so
landing on the whole molecule to 1e-13 says the bookkeeping is right across every
group -- three monomers with boundaries, three dimers, and one of those dimers
the pair of end fragments, which are not bonded to each other and whose group
carries a ghost of a carbon belonging to neither.

The three-body rows are the three-body term, and across covalent bonds it is
large. Expect that: the same quantity is a rounding error for a water cluster and
110 kcal/mol here. Truncating at pairs is not advisable across detached bonds.

**And the point charges make the truncated expansion worse, not better** --
0.489 Hartree against 0.180 with no embedding at all, where on a water cluster
the embedding is worth a factor of twenty. That is not a fault in the
bookkeeping, and it is not the fragments' charges either, which was the
standing explanation here until the two conventions were run side by side.
``embedding = "ptc"`` makes *every* fragment distant, including the one on the
other end of the cut bond, and a point-charge field is at its worst at bonding
contact, which is why FMO keeps an exact term inside ``resppc`` in the first
place. That is the whole of it. None of these rows moves with the charge
convention except the unembedded one, and the full-order rows land on the
supermolecule to 1e-13 regardless.

The unembedded row is why the split is not unconditional. Splitting the nucleus
there takes it from 0.180 to 0.304 Hartree, and on the same molecule numbered so
that one carbon is the detached end of *both* bonds -- and so presents ``Z-2``
-- from 0.125 to 1.549. With no field there is nothing holding the other half,
so each monomer is solved around a nucleus short by a proton. The embedded rows
do not move at all. So the convention is chosen per embedding, and nothing is
given up by doing so.

Restrictions
~~~~~~~~~~~~

``bond_breaking = "afo"`` runs with ``embedding = "none"`` and with
``embedding = "ptc"``.  ``embedding`` is read straight through as the field,
whichever ``method`` was named, so ``"fmo"`` with ``"ptc"`` is FMO's own
expansion over a point-charge field and is the pairing a detached bond runs in.
A spelling that is none of ``"exact"``, ``"ptc"`` or ``"none"`` is refused; it
used to pass validation and change nothing. It is refused with ``embedding = "exact"``. A frozen
orbital and an embedding field both describe the bond region, so the detached
atom's share has to come out of the field before the two can be used together.
With point charges that share is one number per atom -- the population that put
it there -- and is removed exactly. With an exact density the neighbour term is
a Coulomb contraction over a whole density matrix and has no per-atom part to
remove; inventing one would be the point-charge approximation smuggled into the
path defined by not making it.

Refused by name, rather than answered badly:

* a cut through a ring, where two fragments meet in more than one place
* a bond detached at a hydrogen, which has nothing left to hybridise
* a bond carrying more than one localized orbital, which is not a single bond
* a model system that comes out with an odd electron count, which means a cap
  hydrogen was asked to close a valence worth more than one pair

Where the refused bond is a backbone amide, the message names the
C-alpha--C(=O) bond to cut instead, by atom index.

.. _protein-backbone:

Cutting a protein backbone
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Cut the C-alpha--C(=O) bond, not the peptide bond.** This is FMO's convention
for proteins and it is not a matter of taste here: the amide C--N is conjugated
with the carbonyl, Boys localization puts two orbitals on it, and it is refused
by name as "not a single bond". The C-alpha--C(=O) bond one place along is a
nonpolar single bond outside that conjugation, and it is accepted.

So a fragment is one residue's carbonyl together with the *next* residue's
amine and C-alpha, and the peptide bond stays whole inside it. For the glycine
tripeptide in ``sample_inputs``, with atoms numbered from zero as a deck numbers
them:

.. code-block:: json

   "fragments": [
     [0, 1, 4, 5, 6, 7],
     [2, 3, 8, 9, 12, 13, 14],
     [10, 11, 15, 16, 17, 18, 19, 20, 21, 22, 23]
   ]

No ``connectivity`` is needed -- cut bonds are perceived from the geometry, not
declared. Error messages number atoms from **one**, so the bond reported as
"atoms 10 and 11" is the one between the deck's atoms 9 and 10.

**Charged residues are declared.** ``fragment_charges`` gives each fragment's net
charge as it would be with its cut bonds closed -- +1 for a lysine or arginine,
-1 for an aspartate or glutamate, the termini wherever they fall -- and the
charges have to add up to ``molecular_charge``, or the deck is refused. The
electron a detached bond moves is counted on top. A charged system that gives no
``fragment_charges`` is refused too: every fragment would be solved neutral.

The model system around a cut is closed with neutral caps, so any charge it
holds is a group it took in whole. An ammonium (a nitrogen with four
neighbours), a guanidinium (a carbon between three such nitrogens) and a
carboxylate (a carbon with two terminal oxygens) are recognised and counted --
the first C-alpha--C(=O) cut of a protein takes the N-terminal ammonium in, and
without that the model was a radical and the cut was refused. Any other charged
group in a model sphere still is, with a message that says so.

The field has to be point charges: ``"embedding": "ptc"`` alongside
``"bond_breaking": "afo"``. An exact field is refused with a detached bond for
the reason above, and no field at all leaves each side of a cut carrying about
plus or minus one elementary charge, which on a protein puts a spurious ``1/R``
monopole on every adjacent-residue pair.

The per-pair interaction energies -- the numbers a protein-ligand analysis is
for -- are in the output file and the log; see :ref:`fmo-pairs`. The two pairs a
detached bond joins are flagged there and are not interaction energies.

.. _fmo-pairs:

Pair interaction energies
-------------------------

Every two-fragment term of the expansion is kept and reported, rather than only
summed into the total. At info level the log carries a table, strongest first::

   fmo: pair interaction energies, strongest first
   fmo:     pair  R/Angstrom          dE/Hartree   dE/kcal/mol    Tr(dD u)/Hartree
   fmo:    3-4       1.893   ...
   fmo: joined by a detached bond -- each term carries the bond and is not an interaction energy
   fmo:    1-2       ...

and the output JSON an ``fmo`` object beside ``total_energy``:

.. code-block:: json

   "fmo": {
     "expansion": "fmo",
     "embedding": "ptc",
     "monomer_sum": -761.2,
     "pair_sum": -1.07,
     "response_sum": 0.0021,
     "level_sums": [-761.2, -1.07],
     "connected_pair_note": "...",
     "pairs": [
       {"fragments": [3, 4], "distance": 1.893, "connected": false,
        "delta_energy": -0.0055, "interaction_energy": -0.0055,
        "response": 0.0003},
       {"fragments": [1, 2], "distance": 1.52, "connected": true,
        "delta_energy": -17.6, "response": 0.0011}
     ]
   }

(values illustrative). The fields:

``fragments``
   The two fragments, **numbered from one** in the order of the deck's
   ``fragments`` list, lower first. Atom indices elsewhere are 0-based; these
   are fragments.
``distance``
   The closest approach of any atom of one to any atom of the other, in
   **Angstrom**, real atoms only. EFMO's ``distance`` is a unitless vdW-scaled
   ``R_IJ``; this one is not.
``delta_energy``
   The pair's term of the expansion in Hartree, after both monomers are taken
   off: ``E'_IJ - E'_I - E'_J + Tr(dD_IJ u_IJ)`` under the FMO expansion. On
   every row, and what the second entry of ``level_sums`` adds up.
``interaction_energy``
   The same number, written **only** where it is one: under the FMO expansion
   (or with no field at all) and on a pair no detached bond joins. This is the
   pair interaction energy FMO calls an IFIE.
``response``
   ``Tr(dD_IJ u_IJ)``, the pair density's response to the field of everything
   outside it. Already inside ``delta_energy``; zero with no field.
``connected``
   A detached bond joins the two fragments. Sorted after every other pair.

``pair_sum`` is every term of two or more fragments, so above level two it is
more than the pairs add up to. ``level_sums`` holds one sum per term size, the
first entry being ``monomer_sum``; at level three the pairs fall short of
``pair_sum`` by exactly the third entry, the three-body sum, and no row stands
in for the three-body terms. The total is ``monomer_sum + pair_sum`` and is not
changed by any of this.

**What these are.** The interaction of two fragments that have both already been
polarized by the whole rest of the system, *excluding* their mutual induction:
that relaxation happened in the monomer self-consistency, so it sits in the
monomer terms and not in the pair. On a glycine tripeptide with a water
hydrogen-bonded to a carbonyl (HF/6-31G, point charges, frozen orbitals) every
ligand IFIE comes out less attractive than the vacuum pair energy of a plain
expansion on the same partition, by more the closer the contact -- 3.7 kcal/mol
on the hydrogen bond. The ranking of residues agrees.

**What they are not.**

* **Not a binding energy when summed.** On that same system the three ligand
  IFIEs add to +0.60 kcal/mol, against a supermolecular binding energy of
  -6.02: the induction they exclude is most of the binding. A binding energy
  from FMO means two runs, with and without the ligand, and a difference of
  totals.
* **Not interaction energies where ``connected`` is true.** The monomers hold a
  split nucleus and a frozen orbital the pair restores, so the term carries the
  bond itself, about -17 Hartree on a peptide. Such a pair has no
  ``interaction_energy`` key, and a ``connected_pair_note`` says why. This is
  the same convention as the MBE fragment table's ``connected`` column.
* **Not interaction energies under EE-MBE.** There a monomer's energy is its
  total embedded energy and already holds its electrostatics with every other
  fragment, so the pair term takes that back out and is a correction: a
  distant ligand pair reads tens of kcal/mol from nothing. The rows are written
  with ``delta_energy`` only, and a ``pair_note`` says so. Use
  ``method: "fmo"`` for pair analysis.
* **Not decomposed.** There is no electrostatics, exchange, charge-transfer or
  dispersion split of an FMO pair; that is PIEDA, and it is not implemented.
  EFMO's far pairs do carry four terms.
* **Not counterpoise-corrected.**

The pairs are built on every rank from the reduced terms, so an MPI run reports
the same pairs as a serial one and no pair crosses a wire.

Limits
------

**Closed shell only.** A fragment with an odd electron count is refused rather
than quietly paired up. A detached bond moves an electron between the two
fragments it joins, and the count checked here is the one after that move: ethane
split into two methyls is 9 and 9 before it and 8 and 10 after.

**Hartree-Fock only, for now.** Every fragment and n-mer is solved with
restricted Hartree-Fock; other methods are not yet wired into these fragment
calculations, and any other ``model.method`` is refused by name. It used to be
ignored: a B3LYP deck ran as Hartree-Fock and reported that total.

**Energies only.** No gradients yet, so geometry optimization and frequencies are
not available through these.
