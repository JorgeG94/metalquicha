Interaction energy of one fragment (``driver: "InteractionEnergy"``)
====================================================================

Name one fragment of a fragmented system -- a ligand in a protein cut into
residues, a solute in a cluster of waters -- and this driver returns every
many-body interaction term that contains it: the ligand--residue pairs, the
ligand--residue--residue triples, and so on up to the fragmentation level, each
one exactly the number the full many-body expansion would give for it. It
computes only the fragments those terms need, which is far fewer than the full
expansion, and it does **not** report a total energy.

.. code-block:: json

   {
     "schema": {"name": "mqc-frag", "version": "1.0"},
     "molecules": [{
       "xyz": "gly3_water.xyz",
       "molecular_charge": 0,
       "molecular_multiplicity": 1,
       "fragments": [[0, 1, 2, 3, 4, 5, 6, 7],
                     [8, 9, 10, 11, 12, 13, 14],
                     [15, 16, 17, 18, 19, 20, 21, 22, 23],
                     [24, 25, 26]],
       "connectivity": [[2, 8, 1], [10, 15, 1]]
     }],
     "model": {"method": "hf", "basis": "6-31g"},
     "keywords": {
       "fragmentation": {
         "method": "MBE",
         "level": 3,
         "reference_fragment": 3
       }
     },
     "driver": "InteractionEnergy"
   }

Two keys do it, and each is refused without the other:

- ``"driver": "InteractionEnergy"`` says *what* is computed. ``interaction_energy``
  is accepted as the same thing; case is not significant.
- ``keywords.fragmentation.reference_fragment`` says *for which fragment*. It is
  **0-based**, like every other index a deck writes: ``0`` is the first entry of
  the molecule's ``fragments`` list. The deck above names the fourth, the water.

``fragmentation.method`` stays ``"MBE"``: that key says how the system is
partitioned, and the driver says what is done with the partition.

The geometry is ``validation/inputs/sample_inputs/gly3_water.xyz``: the glycine
tripeptide ``gly3.xyz`` -- residues at atoms 0--7, 8--14 and 15--23 -- with a water
donating a 1.94 Angstrom hydrogen bond to the middle residue's carbonyl oxygen.

.. warning::

   **Declare the cut bonds.** ``connectivity`` lists the two peptide bonds the
   partition severs, C2--N8 and C10--N15, and that is what makes each residue
   a hydrogen-capped closed shell. Leave it out and the residues are computed
   uncapped, as radicals with dangling valences, and nothing refuses the deck.
   The capping is decided per term: a bond cut between two residues is whole
   again in the pair that holds both.

What it computes
----------------

The many-body expansion writes the energy of a system of fragments as a sum of
corrections, one per term :math:`T` -- a set of up to ``level`` fragments:

.. math::

   \Delta E_T = E_T - \sum_{S \subset T} \Delta E_S ,
   \qquad
   E \approx \sum_{|T| \le L} \Delta E_T .

Every correction is a property of its term alone, so the expansion splits
cleanly into the terms that contain the reference fragment :math:`R` and the
terms that do not. This driver reports the first part:

.. math::

   E_R \quad\text{and}\quad
   \Delta E_\text{int}(R) = \sum_{T \ni R,\; 2 \le |T| \le L} \Delta E_T ,

split by level. :math:`E_R` is the reference fragment's own energy, its
one-body term.

**What the number is.** Because the corrections that do not contain :math:`R`
are the same whether :math:`R` is there or not,
:math:`E_R + \Delta E_\text{int}(R)` is the level-:math:`L` expansion of the
whole system minus the level-:math:`L` expansion of the system without
:math:`R`. :math:`\Delta E_\text{int}(R)` is therefore the energy of bringing
:math:`R` into the rest of the system, with every fragment frozen at the
geometry it has in the complex -- a binding energy without deformation, at the
truncation the level sets.

At a level equal to the number of fragments the expansion is exact, and so is
this: on the deck above at ``"level": 4`` it gives -0.003100698063 hartree, and
three unfragmented runs give
:math:`E(\text{complex}) - E(\text{gly3}) - E(\text{water}) =` -0.003100698052,
the same number to the SCF convergence.

**What it is not.** It is not a total energy, and none is reported. The sum of
the corrections this run computed also includes the terms without :math:`R`
that it needed along the way (see below), and that sum is neither the system's
energy nor anything else with a name. So the output carries no
``total_energy`` key -- a consumer that reads one gets nothing, rather than a
number that looks like an energy. Nor is it a sum of pair interactions: at
level 3 the triples that contain :math:`R` are part of it, and they are not
small.

Why it computes more than the terms that contain R
---------------------------------------------------

Computing only the terms that contain the reference gives the **wrong
answer**. A correction needs the energy of every subset of its term, and most
of those subsets do not contain :math:`R`:

- :math:`\Delta E_{RA} = E_{RA} - E_R - E_A` needs :math:`E_A`.
- :math:`\Delta E_{RAB}` needs :math:`E_{RAB}`, :math:`E_{RA}`, :math:`E_{RB}`,
  :math:`E_R`, :math:`E_A`, :math:`E_B` -- **and** :math:`E_{AB}`, a pair with
  no :math:`R` in it at all.

What has to be computed is the terms that contain :math:`R`, closed under
taking subsets. The rule that produces exactly that set: **a term is computed
when it contains** :math:`R`, **or when adding** :math:`R` **to it gives a term
that is computed.** A subset :math:`S` of an :math:`R`-term either contains
:math:`R` itself, or :math:`S \cup \{R\}` is a subset of the same :math:`R`-term;
either way the rule keeps it, and it keeps nothing an :math:`R`-term does not
need.

With no distance screening that is every term of up to :math:`L` fragments
containing :math:`R`, and every term of up to :math:`L - 1` fragments not
containing it. **The only terms skipped are the** :math:`L`-**mers without**
:math:`R`. At level 2 that is every pair not involving the reference -- all of
them but :math:`n - 1` -- and at level 3 the triples not involving it. For
:math:`n` fragments:

.. math::

   N_\text{computed} = 1 + 2 \sum_{k=1}^{L-1} \binom{n-1}{k}
   \qquad\text{against}\qquad
   N_\text{full} = \sum_{k=1}^{L} \binom{n}{k} .

The list is closed under subsets by construction, so a correction can never
look for a subset that was not computed. Every correction this driver reports
is therefore the one the full expansion computes for the same term, down to the
last bit: that is checked, not assumed (see below).

What it saves
-------------

A cluster of twenty waters, ``w20_isomer1.xyz``, Hartree--Fock/STO-3G, one
water as the reference, on one core:

=========================================  ==========  ==========  ==========
Run                                        Computed    Full        Wall
=========================================  ==========  ==========  ==========
Level 2                                    39          210         0.5 s
Level 3                                    381         1350        3.7 s
Level 3, cutoffs dimer 6.0, trimer 4.5 A   275         631         2.9 s
The ordinary level-3 energy, for scale     --          1350        16.4 s
=========================================  ==========  ==========  ==========

Level 2 computes :math:`2n - 1` terms rather than :math:`n(n+1)/2` -- linear in
the system size rather than quadratic. Level 3 skips the
:math:`\binom{n-1}{3}` triples without the reference, 969 of 1350 here; what
remains is dominated by the pairs every reference triple needs, so the saving
at level 3 is a constant factor that grows with :math:`n` rather than an order.

Distance screening composes with it. A screened term is gone before the
reduction runs, so the reduction starts from the screened list and also drops
the subsets that only a screened-out :math:`R`-term needed: a fragment too far
from the reference to form a pair with it is not computed at all. The level-3
interaction energy above moves from -0.0468899392 to -0.0467545935 hartree
under those cutoffs, from 171 reference triples to 118.

What it reports
---------------

The log prints the reduction as soon as the term list is built, and the
breakdown at the end::

   Interaction energy of reference fragment 3 (0-based): computing 13 of 14 terms, skipping 1
   ...
   Interaction energy of reference fragment 3 (0-based, as in the deck; monomer 4 in the fragment table):
     Reference fragment energy:        -75.9839402918
     2-body terms:                      0.0006784376   (3 terms)
     3-body terms:                     -0.0038855273   (3 terms)
     Interaction energy:                -0.0032070897

and ``output_<name>.json`` carries it in place of ``total_energy``:

.. code-block:: json

   "interaction_energy": {
     "reference_fragment": 3,
     "reference_monomer": 4,
     "reference_energy": -75.98394029179254,
     "total": -0.0032070897117932873,
     "total_kcal_mol": -2.012479178320585,
     "by_level": [
       {"frag_level": 2, "name": "dimers",  "count": 3, "energy":  0.0006784376225397182},
       {"frag_level": 3, "name": "trimers", "count": 3, "energy": -0.0038855273343330055}
     ],
     "terms_computed": 13,
     "terms_in_full_expansion": 14
   }

``reference_fragment`` is the deck's 0-based index; ``reference_monomer`` is
the same fragment as the fragment table and the ``levels`` indices number it,
from 1. ``by_level[].count`` is how many terms at that level contain the
reference, and ``terms_in_full_expansion`` is what the ordinary run would have
computed over the same fragments, level, screening and counterpoise. From
Python, ``MBE(..., driver="InteractionEnergy", keywords={"fragmentation":
{"reference_fragment": 3}})`` runs it and ``Result.interaction_energy`` returns
this block; the result's ``energy`` means nothing for this driver.

The fragment table, ``output_<name>_fragments.csv``, has one row per computed
term with its own correction, so the individual ligand--residue numbers are
there. **It also has the terms without the reference that were computed as
subsets**, with their own corrections; select the rows whose ``m`` columns
contain ``reference_monomer`` before summing anything. The ``levels`` block of
the JSON counts those terms too, and carries no per-level energy, because a
level sum over a reduced list is no level's energy.

Reading the numbers: what level to use
--------------------------------------

The water on gly3 above, Hartree--Fock/6-31G, at every level the four
fragments allow:

==============  ======================  ======================  ============
Level           2-body (Ha)             3-body (Ha)             kcal/mol
==============  ======================  ======================  ============
2                0.0006784376           --                       +0.4257
3                0.0006784376           -0.0038855273            -2.0125
4 (exact)        0.0006784376           -0.0038855273            -1.9457
==============  ======================  ======================  ============

(The level-4 run adds a single four-body term, +0.0001063916.) The pairs alone
get the **sign** wrong. That is not the expansion failing to converge in the
usual sense: each residue in a water--residue pair is a capped fragment, and a
cap sits where the neighbouring residue's atom was, close to the water. The
three-body terms are where the capped pairs are corrected back to the real
peptide. With covalently cut fragments, read a ligand's interaction energy at
level 3 or above -- or cut the protein into fragments large enough that the
caps are far from the ligand.

What composes with it
---------------------

- **Distance screening**, as above.
- **MPI.** The reduced list is distributed like any other; on the prism water
  hexamer at levels 2 and 3 and the three references tested, one, two and three
  ranks give bitwise-identical interaction energies.
- **Counterpoise.** Each term builds its own ghosted rows, so the reduction
  carries over unchanged. For one water of the prism, VMFC(2) computes 21 of
  the ordinary 51 rows and VMFC(3) 121 of 191, and each matches the full run's
  terms holding that water: to the last digit at level 2, and to 1e-12 hartree
  at level 3 against the Valiron--Mayer sum evaluated from the full run's
  fragment table.
- **Checkpoints.** A term's energy is the same whichever driver computed it, so
  a checkpoint written by an ``Energy`` run can seed an ``InteractionEnergy``
  run over the same system, and the other way round.

What is refused
---------------

Each of these is refused before any fragment is computed, with a message naming
the keys involved, rather than approximated or ignored:

- ``driver: "InteractionEnergy"`` without ``reference_fragment``, and
  ``reference_fragment`` with any other driver -- an ``Energy`` deck carrying
  one would otherwise run the full expansion and look as though it had
  honoured it.
- A ``reference_fragment`` that is negative, not an integer, or not less than
  the number of fragments.
- ``"level": 1``, and a system with fewer than two fragments. A one-body
  expansion has no interaction terms.
- **Gradients and Hessians**, by construction: they are other drivers. The
  derivative of the interaction energy would need the derivative of every
  skipped term's absence, and nothing here computes it.
- ``"method": "gmbe"``. GMBE's terms are intersections of overlapping
  primaries, not sets of fragments, so which of them contain the reference is
  not defined.
- ``"method"`` ``fmo``, ``ee-mbe`` and ``efmo``. Those build their own term
  lists -- and under an embedding every fragment feels the field of all the
  others, so skipping one would change the ones computed.
- EFP, SAPT and NEO, which never reach the many-body expansion.
- A multi-molecule deck. A fragment index names a different fragment -- or
  none -- in each molecule.
- A supplied term list, through the Python or C interface: the reduction is
  where the saving is, and a list closed for the full expansion is not the
  reduced one.

Bonding analysis on the reference's terms
-----------------------------------------

With ``properties.bonding_analysis`` (:doc:`bonding_analysis`) on the deck, the
quasi-atomic analysis runs in every term's SCF, and the terms that hold the
reference and at least one other fragment are reported together at the end. Only
what crosses between the reference and the other fragments of the term is
kept: that is the part that describes their contact.

.. code-block:: json

   "properties": {
     "bonding_analysis": {"type": "gms_quao", "energy_threshold": 1.0}
   }

On the gly3 and water deck above at STO-3G and level 2, after the interaction
energy::

    Bonding of reference fragment 3 with its environment (quasi-atomic, |kinetic bond order| >= 1.00 kcal/mol)
      atoms and monomers are 0-based; reference atoms are marked *

      term 2: monomers 1 3
         atom pair                  index    kcal/mol
         H 25*      - O 11         0.0147       -1.93
         orbital pair, donor first                           direction       order    kcal/mol
         O 11-C 10 pi             -- H 25*-O 24 sigma        env -> ref     0.0924       -1.32

The water's hydrogen 25 is bonded to the middle residue's carbonyl oxygen 11,
and the orbital row says how: the C=O pi orbital on oxygen 11 donates into the
water's O--H sigma orbital. The other two water--residue terms have nothing at
1 kcal/mol and print ``(none above threshold)``.

There are two tables per term, and they are not equally robust:

- **Atom pairs** sum over every orbital on the reference atom against every
  orbital on the other: the ``index`` is the sum of squared bond orders, and
  the kcal/mol is the sum of kinetic bond orders. A rotation within one atom
  changes neither, so they do not depend on how well the orbitals were
  oriented. Reference atom first; strongest bonding first.
- **Orbital pairs** are the rows of the bonds and delocalization tables with
  one end on each side, labelled as the full report labels them. A
  delocalization row puts the donor first and says which way it goes. These do
  depend on the orientation; a term whose orientation stopped at its sweep
  limit (see :doc:`bonding_analysis`) says so under its table.

Both use the deck's ``energy_threshold``. A row or pair with an end on a
hydrogen cap or a ghost atom is dropped, since neither is an atom of the
system; the count of dropped orbital pairs is printed and written. A bond the
partition cut is whole inside a term holding both of its ends, so it appears
here as an ordinary ``bond`` row.

The same content goes into the JSON output, inside ``interaction_energy``:

.. code-block:: json

   "bonding": {
     "analysis": "gms_quao",
     "threshold_kcal_mol": 1.0,
     "terms": [
       {"id": 2, "monomers": [2, 4], "fragments": [1, 3],
        "orientation_stalled": true, "omitted_orbital_pairs": 0,
        "atom_pairs": [
          {"reference_atom": 25, "reference_element": "H",
           "partner_atom": 11, "partner_element": "O",
           "bond_index": 0.014652935461, "kinetic_bond_order": -1.934518468024}],
        "orbital_pairs": [
          {"kind": "delocalization", "direction": "environment_to_reference",
           "bond_order": 0.092424801839, "kinetic_bond_order": -1.321518079113,
           "ends": [
             {"atom": 11, "element": "O", "fragment": 1, "on_reference": false,
              "orbital": 19, "type": "pi", "occupation": 1.191999144051,
              "bonded_to": 10},
             {"atom": 25, "element": "H", "fragment": 3, "on_reference": true,
              "orbital": 29, "type": "sigma", "occupation": 0.755751628455,
              "bonded_to": 24}]}]}
     ]
   }

Atoms and ``fragments`` are 0-based, as the deck numbers them. ``id`` is the
term's row in the fragment table, and ``monomers`` is the same membership
1-based, as ``levels`` writes it. ``orbital`` is numbered as that term's
molecular orbitals are. ``bonded_to`` is the atom a bonding orbital's bond goes
to, and is absent for a lone pair or for a bond to a cap. Every term holding the
reference and another fragment is listed, including one with nothing above
threshold.

The per-term numbers are not combined across terms. A bond order belongs to the
calculation it came from, and there is no many-body correction of one to take:
the dimer and the trimer each describe the contact in their own surroundings.
Serial and MPI runs give the same tables.

What it does not do yet:

- The analysis also runs on the terms that do not hold the reference (the
  subsets the interaction energy needs), and its full per-fragment tables are
  still printed to the log for every term, numbered within that fragment. That
  costs about half a second per term at STO-3G without the energy
  decomposition. With ``energy_decomposition`` on, the cost is the dense
  two-electron transformation, and it is paid on every term.
- ``energy_decomposition`` is not written for any fragmented run.
- A term restored from a checkpoint carries no bonding tables, and is left out
  of the report.

Where it is checked
-------------------

- ``test/test_mqc_term_list.f90``: the reduced list has the length the formula
  above gives, for every system of 2 to 7 fragments, every level up to 4 and
  every choice of reference; it is closed under subsets, contains every term
  holding the reference, and contains nothing else no such term needs, with and
  without distance screening.
- ``test/test_mqc_mbe.f90``: the interaction energy assembled from the reduced
  list equals the full expansion's, worked out independently by
  inclusion--exclusion, for every reference; no total energy and no dipole come
  out of it.
- ``test/test_mqc_config_roundtrip.f90`` and ``test/test_mqc_json_reader.f90``:
  every refusal above, by its message.
- ``validation/inputs/cpu/mqc/interaction_energy/``: the prism at level 2, 3,
  VMFC(2) and VMFC(3), and the gly3 and water deck, each pinned to the sum of
  the same terms' corrections from the ordinary run of the same deck.

The identity itself -- every correction containing the reference equal to the
full expansion's -- was checked term by term on the prism water hexamer
(HF/STO-3G, references first, middle and last, levels 2 and 3, one to three
ranks), the gly3 and water deck (HF/6-31G, levels 2 and 3) and the twenty-water
cluster at level 3: every correction agrees bit for bit.
