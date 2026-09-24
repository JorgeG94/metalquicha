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

``resdim`` (default ``2.0`` for FMO2, ``0`` otherwise)
   How far apart two fragments have to be before their pair is not solved at all
   but taken as the electrostatic interaction of the two monomers, measured as
   ``resppc`` is. GAMESS calls such a pair a separated or ES dimer and marks it
   ``D=S``; here its row in the pair table carries ``ES``, and ``separated`` in
   the JSON output. Its energy is

   .. math::

      \Delta E_{IJ} = \mathrm{Tr}(D^I u^J) + \mathrm{Tr}(D^J u^I)
         + \sum_{\mu\nu \in I} \sum_{\lambda\sigma \in J}
           D^I_{\mu\nu} D^J_{\lambda\sigma} (\mu\nu|\lambda\sigma)
         + E^{\rm nuc}_{IJ}

   with each monomer's converged density and its nuclei as its own SCF presented
   them, ghosts and detached-bond charges included, and no response term. The
   repulsion is exact, as GAMESS computes it by default.

   ``0`` solves every pair. The default is GAMESS's for FMO2; at level three
   GAMESS ties it to trimer trimming, which this code does not do, so a deck that
   says nothing solves every pair there. EE-MBE has no pair interaction energy for
   the approximation to stand in for, and refuses a non-zero value.

   Twenty waters in STO-3G separate 64 of their 190 pairs, and match GAMESS
   ``RESDIM=2.0`` to 5.5e-08 in the total and to 5e-09 in each separated pair; the
   separation moves the total by 6.4e-06 against solving every pair.

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

``detached_atoms`` (default: the sp3 end, else the lower-numbered one)
   0-based atoms that are the detached ends of cut bonds; see below.

``afo_localization`` (default ``"er"``)
   How the model system around each cut bond is localized: ``"er"`` for
   Edmiston-Ruedenberg, which is GAMESS's default and the paper's, or
   ``"boys"`` for Foster-Boys. Read only with ``bond_breaking = "afo"``.

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

The construction is GAMESS's adjusted frozen orbitals, as Fedorov, Jensen, Deka
and Kitaura describe it (*J. Phys. Chem. A* **112**, 11808 (2008)) and as
``fmolib.src`` implements it, and it reproduces GAMESS's fragment energies to
1e-8 Hartree on a single cut with Boys and 3e-8 with ER, as far as GAMESS's own
model-system convergence allows; see `Against GAMESS`_.

**Which end is detached.** Of the two atoms of a cut bond one is the *detached*
end (GAMESS's BDA) and one the *attached* end. ``keywords.fragmentation.
detached_atoms`` names the detached ends, 0-based, as the sign in GAMESS's
``$FMOBND`` does; naming both ends of one bond is refused. A bond neither end of
which is named is detached at its sp3 end -- the one with four neighbours when
the other has fewer, which on a protein backbone cut at C-alpha--C(=O) is the
C-alpha -- and otherwise at its lower-numbered atom.

**The model system.** Both atoms of the bond and every atom bonded to either,
by GAMESS's own bond test (Emsley radii times 1.2, hydrogen's unscaled); then,
for each atom bonded to one of those and not taken, a hydrogen comes in where it
is and anything else is replaced by a cap hydrogen at GAMESS's X-H length for the
atom it hangs off (1.09 A from carbon, 1.01 from nitrogen, 0.96 from oxygen). One
addition GAMESS does not make: a *terminal* heavy atom bonded to a taken one --
a carbonyl oxygen -- comes in whole, since a cap closes one electron pair and a
double bond is two. It never applies where GAMESS's model is closed-shell
already. Charged groups taken in are counted into the model's charge.

**The orbitals.** The model is solved and every occupied orbital
localized, Edmiston-Ruedenberg unless ``afo_localization`` says Boys. The detached atom's own orbitals are the ones with the largest
population on it, ``sum_{mu,nu on A} C_mu S_mu,nu C_nu`` -- five for a carbon,
its 1s and four sp3 -- and of those the one with the largest population on the
attached atom is the bond's. Each is kept on every real atom bonded to either
end and dropped elsewhere; a group takes whichever of those atoms it holds and
Gram-Schmidt orthonormalises the set, occupied first.

**Which localizer.** Edmiston-Ruedenberg by default, as in the paper and in
GAMESS (``$CONTRL LOCAL``, read into the model at ``fmolib.src:5783``); Boys was
this code's only choice until ER was added. The two split the detached carbon's
five orbitals differently, so the *monomers* move by a tenth of a Hartree --
0.086 and -0.105 on butane with a water, in GAMESS and here alike -- while the
pair holding the whole bond does not see it. What reaches a total is smaller but
not negligible: in GAMESS's FMO2 the choice moves butane with a water, cut once,
by 1.4e-8 Hartree, and the glycine tripeptide with a water, cut twice, by 2.5e-4
-- with ER 3.9e-4 above the molecule's energy and with Boys 6.4e-4. Both are
compared against GAMESS below.

The two fragments then split the bond:

============================  =========================  =========================
                              fragment of the detached   fragment of the attached
============================  =========================  =========================
nucleus                       ``Z - 1`` of it            ``+1`` of it, on a ghost
its basis functions           owns them                  carries them, ghosted
the bond's electron pair      none of it                 all of it
the bond orbital              frozen empty               frozen occupied
the atom's other orbitals     free                       frozen empty
electron count                ``sum(Z) - 1``             ``sum(Z) + 1``
net charge                    ``0``                      ``0``
============================  =========================  =========================

**The atom's other orbitals are projected out of the ghost** -- for a carbon
its 1s and its three other sp3, GAMESS's "1 occupied and 4 virtual frozen
LMOs". Without it the borrowed functions are free to hold anything, and the
fragment's electrons spread into the directions of the detached atom's other
bonds.

**What the old construction got wrong.** Until this was matched to GAMESS the
bond orbital was reduced to the detached atom's own functions and frozen alone.
On butane with a water that left the two cut monomers 0.127 and 0.298 Hartree
above GAMESS's. Put back one piece at a time: the bond orbital over its
neighbours as well carries all of the 0.127 and most of the 0.298; the four
extra empties are the last 0.026, and in the other direction -- without them
the ghost-holding monomer comes out *below* GAMESS's, as removing a constraint
must. With the point-charge field the old construction also over-polarised the
monomers at every cut: the monomer loop's energy change rose from 0.34 to 1.05
Hartree over the first three passes on a two-fragment glycine tripeptide and
took twenty passes. It now falls from the first pass, 5.0e-3, 1.1e-3, 1.5e-4,
and stops at seven; GAMESS takes nine on the same system.

The nucleus is split because the electron pair is. One unit of charge crosses
the bond with the pair, so both fragments come out neutral closed shells and
the two halves add back to ``Z`` inside any n-mer holding the whole bond. That
is GAMESS's convention, it is what a fragment *potential* needs -- a fragment
carrying a unit charge puts a monopole term of order ``1/R`` on every
adjacent-residue pair of a protein -- and it is the default with or without a
field. ``embedding = "none"`` used to keep the nucleus whole instead, which was
measured better when only the bond orbital was frozen; with the full frozen set
it is not (`What it costs`_), and on the glycine tripeptide with a water the
whole convention's C-terminal fragment is an anion whose SCF does not converge
with a Boys model, and with an ER one converges to an MBE(2) error of 2.5e-2
against the split convention's -3.7e-3.

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
the nucleus with the bond orbital frozen empty, and the fragment across the bond
holds the same functions as a ghost carrying ``+1`` with the bond pair in that
orbital. Three things follow, and all three are what ``embedding = "ptc"`` needed
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
total is -116.59715762085884 with whole nuclei and -116.59715762085621 with
split ones, with the monomer loop converged to 1e-10. The two runs are in the unit suite side
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
ordinary RHF on the whole molecule, error in Hartree with the model system
localized each way:

=================================  ===========  ===========
Expansion                          ER           Boys
=================================  ===========  ===========
Two fragments, one bond, MBE(2)    exact        exact
Three fragments, MBE(3)            1.4e-13      1.8e-13
Three fragments, MBE(2)            0.167        0.219
Three fragments, FMO(3), ``ptc``   4.3e-14      1.6e-13
Three fragments, FMO(2), ``ptc``   0.145        0.200
=================================  ===========  ===========

An expansion carried to the fragment count is exact by inclusion and exclusion
whatever the partition did, so landing on the whole molecule to 1e-13 says the
bookkeeping is right across every group. The two-body rows are the three-body
term, and across detached bonds on fragments this small it is large; with the
frozen bond orbital alone they were 0.180 and 0.489.

Propane's two ends are a bond apart from the middle carbon each detaches, so the
dimer of the two ends holds a ghost of that carbon beside the bond to it that it
owns. GAMESS freezes the same bond twice there, once as the owner's empty bond
orbital and once among the ghost's other orbitals, and its SCF for that dimer
does not converge. Here the ghost's copy is left out. On a protein backbone
cut at C-alpha--C two cuts are never that close.

The nucleus convention in vacuo, MBE(2) error in Hartree, ER (Boys):

=====================================  ================  ================
Propane, three fragments               split             whole
=====================================  ================  ================
carbons in chain order                 0.1666 (0.2186)   0.2418 (0.2200)
middle carbon detached twice (Z-2)     0.809 (0.868)     0.032 (0.046)
=====================================  ================  ================

``"whole"`` is still better where one atom is detached from two neighbours,
which GAMESS refuses outright; it can be asked for through
``fmo_options_t%cut_nucleus``. With a field the two conventions give the same
total, to 2.6e-12 on propane with the monomer loop converged to 1e-10.

The glycine tripeptide with a water hydrogen-bonded to its middle carbonyl
(``gly3_water_pair.xyz``), cut at both C-alpha--C bonds into four fragments,
RHF/STO-3G, against the molecule's -762.311670856, error in Hartree:

=================================  ===========  ===========
Expansion                          ER           Boys
=================================  ===========  ===========
MBE(4), in vacuo                   2.0e-12      2.3e-12
MBE(2), in vacuo                   -3.7e-3      -3.1e-3
FMO(2), ``ptc``                    -4.5e-4      -6.2e-5
FMO(2), ``exact``                  +3.9e-4      +6.4e-4
GAMESS FMO2, same field            +3.9e-4      +6.4e-4
=================================  ===========  ===========

With the exact field the two codes agree to 4.2e-8 in the total with either
localizer (``RESPPC=2.0 RESDIM=0 RESPAP=0`` in GAMESS, ``resppc`` 2.0 here):
-762.311281612 against -762.311281570 with ER, -762.311027620 against
-762.311027578 with Boys. GAMESS gives the same ER total at its default
``RESDIM`` and at ``$LOCAL CVGLOC=1D-10``. Pair energies against GAMESS's PIEDA
totals, kcal/mol:

============  =================  =================  =================  =================
pair          ER, here           ER, GAMESS         Boys, here         Boys, GAMESS
============  =================  =================  =================  =================
3-4           -5.4748            -5.475             -5.5039            -5.504
1-3           1.8986             1.899              1.8558             1.856
2-4           1.2650             1.265              1.2370             1.237
1-4           0.0747             0.075              0.0740             0.074
1-2 (cut)     -9109.3723         -9109.379          -9120.1131         -9120.121
2-3 (cut)     -9105.1257         -9105.133          -9116.4139         -9116.422
============  =================  =================  =================  =================

Every pair across no cut agrees to the digit GAMESS prints. **The cut pairs
differ by 7-8e-3 kcal/mol, and it is GAMESS's printout, not the energy:**
GAMESS's printed monomers plus its printed pair totals miss its own FMO2 total
by 1.3e-5 Hartree per cut -- 2.3e-5 here, 1.3e-5 on butane, with ER and with
Boys alike -- so the printed cut-pair row is not the term its total is built
from. Ours is, and summed with our monomers it gives our total to 1e-10.
Pinned, with the totals, in
``glycine_tripeptide_and_water_match_gamess_fmo2_exact_field_er``, and without
the ``_er`` for Boys.

The response part of the water pair, ``Tr(dD u)``, is 0.718 kcal/mol in both
with Boys and 0.763 with ER. Here the monomer loop stops after 6 passes at the
default ``outer_tolerance``, GAMESS's after 9; the two stop on
different measures (the change in the monomer energy sum here, density and
energy there).

Against GAMESS
~~~~~~~~~~~~~~

Butane cut at C2-C3 with a water 4.5 A beyond C4, RHF/STO-3G, both codes solving
every fragment in vacuo (GAMESS as the in-vacuo half of an EFMO run,
``RAFO(1)=1,1,1``, ``LOCAL=BOYS``):

========================================  ===================  =========
                                          GAMESS               ours - it
========================================  ===================  =========
ethyl owning C2 at ``Z-1``                -63.7464859666       4.7e-9
ethyl holding the ghost of C2             -77.4085157350       6.1e-9
water                                     -74.9620085207       -9.0e-9
MBE(2) of the three                       -230.4152004258      -1.7e-8
========================================  ===================  =========

Pinned in ``butane_and_water_match_gamess_afo``, with ``afo_localization =
"boys"``. The same system with the model ER-localized, against GAMESS with
``LOCAL=RUEDNBRG`` and ``$LOCAL CVGLOC=1D-10`` (its FMO default, ``1D-7``, leaves
the ethyls 1e-8 short of converged), every pair an SCF in vacuo:

========================================  ===================  =========
                                          GAMESS               ours - it
========================================  ===================  =========
ethyl owning C2 at ``Z-1``                -63.6608350059       -8.5e-9
ethyl holding the ghost of C2             -77.5130688706       3.3e-8
water                                     -74.9620085207       -9.0e-9
MBE(2) of the three                       -230.4152073207      -1.8e-8
========================================  ===================  =========

Pinned in ``butane_and_water_match_gamess_afo_er``. The ghost-holding ethyl is
the loosest, and it is GAMESS's model system rather than the localizer: GAMESS
stops that SCF at a density change of 1e-6, and an ER frozen set follows the
model's orbitals about four times as closely as a Boys one -- loosening our
model from 1e-10 to 1e-7 in the energy moves these monomers 1.3e-7 with ER and
3e-8 with Boys.

The glycine tripeptide with a water, cut twice, in vacuo the same way
(``glycine_tripeptide_and_water_match_gamess_afo_er``): the four monomers agree
with GAMESS to 1.8e-7 and MBE(2), -762.3153911672 there, to 4.1e-8. That is not
an ER discrepancy -- the same run with Boys against GAMESS's Boys is 2.0e-7 and
4.1e-8 -- but the same loose model SCF, on bigger models.

**The model system converges on its own terms.** Its SCF used to take the
fragments' tolerances and derive its commutator bound from them, 1e-5 at the
FMO defaults, although what leaves the model is its orbitals, which are first
order in that bound. It now runs to 1e-11 in the energy and 1e-9 in the
commutator whatever the fragments ask for, which a dozen atoms makes free.
At the default fragment tolerances that moves the tripeptide's exact-field
monomers by up to 1.3e-7 with ER and 1.6e-7 with Boys, and its totals by
4e-10; every comparison with GAMESS above holds, because GAMESS's own model
stops at a density change of 1e-6 and that remains the residual.

Pinned in ``butane_and_water_match_gamess_afo``. FMO2 in the exact field on the
same system, ``RESPPC=2.0 RESDIM=0 RESPAP=0`` in GAMESS and ``resppc`` 2.0 here:
-230.415266010 against GAMESS's -230.415265994 with Boys, 1.6e-8 apart, and
-230.415266021 against -230.415266004 with ER, 1.7e-8 apart; pinned in
``butane_and_water_match_gamess_fmo2_exact_field`` and ``..._er``. The two
pairs across no cut match GAMESS's printed 0.302 and -0.158 kcal/mol; the cut
pair is -9122.239 against a printed -9122.246, the printout offset above.
At GAMESS's default ``RESDIM=2.0`` the water and the ethyl it does not
hydrogen-bond are a separated pair, both of them fragments with a detached bond:
-230.415265614 against -230.415265534 with Boys and -230.415265628 against
-230.415265548 with ER, 8.0e-8 apart, the separated pair within 4.4e-9 of
GAMESS's printed -0.00024981 and -0.00025100 Hartree; pinned in
``butane_and_water_match_gamess_with_a_separated_pair``. With point charges it is
-230.415277443; the molecule is -230.415265709.

Restrictions
~~~~~~~~~~~~

``bond_breaking = "afo"`` runs with every field: ``"none"``, ``"ptc"`` and
``"exact"``. A spelling that is none of those is refused; it used to pass
validation and change nothing.

**With the exact field each neighbour acts as itself**, as GAMESS builds it
(``FMOESP`` in ``fmoint.src``): a near fragment -- within ``resppc``, which a
bonded neighbour always is -- through its nuclei as its own monomer presents
them, ``Z-1`` on a detached atom it owns and ``+1`` on a ghost it holds, and
through the Coulomb operator of its whole density over its whole basis, ghost
functions included; a distant fragment through its atomic charges, ghost
included. A detached atom's two shares therefore arrive from the two fragments
that hold them, and a group's own share is never put in, so nothing has to be
taken back out. This used to be refused on the grounds that the exact
neighbour term has no per-atom part to remove, which was true and beside the
point: removing it was never necessary. At full order the expansion lands on
the molecule to 1e-13 (``three_fragments_are_exact_under_an_exact_field``).

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

Use the exact field, ``"embedding": "exact"`` alongside
``"bond_breaking": "afo"``: it is FMO's own field and what GAMESS runs, and on
the glycine tripeptide with a water it matches GAMESS's total to 4e-8 and its
pair energies to the digit it prints. Point charges at bonding distance are the
field's worst case.

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
