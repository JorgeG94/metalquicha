Effective Fragment Molecular Orbitals (EFMO)
============================================

EFMO computes a cluster's energy as in-vacuo fragments and near groups, effective
fragment potentials for the far pairs, and one many-body induction over every
fragment at once. Each fragment's potential is built on the fly by MAKEFP, from
the very SCF that supplies the fragment's own energy, so a run needs no ``.efp``
files and no second calculation.

The near half runs to **any many-body order** -- ``keywords.fragmentation.level``,
default 2, which is EFMO as published and as GAMESS runs it. See
*Beyond pairs* below.

The method is Steinmann, Fedorov and Jensen, *J. Phys. Chem. A* **114**, 8705
(2010), in the form of Sattasathuchana *et al.*, *J. Chem. Theory Comput.* **20**,
2445 (2024), whose eq 6 is the energy expression below.

The energy
----------

.. math::

   E = \sum_{S\ {\rm near},\ |S| \le n}
         \left( \Delta E_S^0 - \Delta E_S^{\rm pol} \right)
     + \sum_{I<J,\ R_{IJ} > R_{\rm cut}}
         \left( E_{IJ}^{\rm Coul} + E_{IJ}^{\rm disp}
              + E_{IJ}^{\rm ExRep} + E_{IJ}^{\rm CT} \right)
     + E_{\rm pol}^{\rm total}

At the default level :math:`n = 2` that is eq 6 of the paper written out,

.. math::

   E = \sum_I E_I^0
     + \sum_{I<J,\ R_{IJ} \le R_{\rm cut}}
         \left( E_{IJ}^0 - E_I^0 - E_J^0 - E_{IJ}^{\rm pol} \right)
     + \sum_{I<J,\ R_{IJ} > R_{\rm cut}}
         \left( E_{IJ}^{\rm Coul} + E_{IJ}^{\rm disp}
              + E_{IJ}^{\rm ExRep} + E_{IJ}^{\rm CT} \right)
     + E_{\rm pol}^{\rm total}

because :math:`\Delta E_I^0 = E_I^0`, :math:`\Delta E_{IJ}^0 = E_{IJ}^0 - E_I^0 -
E_J^0`, :math:`\Delta E_I^{\rm pol} = 0` and :math:`\Delta E_{IJ}^{\rm pol} =
E_{IJ}^{\rm pol}`.

:math:`E_S^0` is **in vacuo**: no embedding field, no monomer
self-consistency. That is what separates EFMO from FMO, and it is what lets
diffuse basis sets work -- there is no neighbouring point charge for a diffuse
function to collapse onto.

.. note::

   **Mutually polarized QM/EFP embedding is a separate method, not planned as
   part of this one.** Solving a fragment's SCF in the field of the other
   fragments' potentials, with the induction inside the SCF, makes a group's
   energy depend on its environment -- so the many-body differences below stop
   telescoping and the level = N identity stops holding. Nothing in this page
   embeds anything; every SCF an EFMO run performs is on an isolated group of
   fragments.

:math:`E_S^{\rm pol}` is the induction energy of the isolated group :math:`S`
-- at :math:`|S| = 2` the isolated pair. Every quantum group already contains its
own fragments' mutual induction, and :math:`E_{\rm pol}^{\rm total}` -- the
induction solved over every fragment together -- contains it again, so one copy is
removed. What is left,
:math:`E_{\rm pol}^{\rm total} - \sum_S \Delta E_S^{\rm pol}`, is the induction
the truncated expansion does not reach, and at level two it is **not small**: for
three waters at four angstrom it is 44 per cent of the total. The energy is
quadratic in the field, and the square of a sum keeps cross terms no pair has. At
level = N it is zero, exactly, and that is a test rather than a coincidence.

The cutoff
----------

.. math::

   R_{IJ} = \min_{i \in I,\ j \in J}
            \frac{|\mathbf{r}_i - \mathbf{r}_j|}{r_i^{\rm vdW} + r_j^{\rm vdW}}

**Unitless.** Each interatomic distance is divided by the two van der Waals radii
of GAMESS's ``$FMO VDWRAD`` table -- hydrogen 1.20 and oxygen 1.40 Angstrom, not
Bondi's 1.10 and 1.52. The table matters: the FMO literature quotes ``RESPPC``,
``RESDIM`` and :math:`R_{\rm cut}` against *this* one, and the two disagree by six
per cent on a water pair, which is enough to move a dimer across
:math:`R_{\rm cut} = 2.0`. On seven water clusters at three cutoffs the quantum
and effective dimer counts now match GAMESS's exactly, case for case.

Dividing by the radii is what makes it unitless, so :math:`R_{IJ} = 1` is contact
and the default :math:`R_{\rm cut} = 2.0` is twice that -- a threshold that means the same thing for a water pair and for two aromatic
rings, which an angstrom threshold does not. It is *not* comparable to
``keywords.fragmentation.cutoffs``, which MBE uses and which is in angstrom. It is
FMO's ``resppc`` measured the same way, deciding a different question.

A value at or below zero is refused: it would leave no pair quantum mechanical at
all, which is EFP with in-vacuo monomers rather than the method the deck asked for.

Beyond pairs
------------

``keywords.fragmentation.level`` -- the same key MBE and FMO read; there is no
EFMO-specific level -- truncates the many-body expansion of the near groups.
Default 2. Level 1 is the fragment sum alone; the fragment count is exact.

**Both series are differenced, and that is the whole of the generalization.**
:math:`\Delta E_S^0` is the many-body difference of the in-vacuo energies of
:math:`S` and its subsets, exactly what FMO's expansion already does.
:math:`\Delta E_S^{\rm pol}` is the *same* difference applied to
:math:`E_S^{\rm pol}`, the induction energy of that group's own effective
fragment potentials -- the existing induction solver, same screening and same
convergence, run on the group's fragments instead of on a pair. So at level
three

.. math::

   \Delta E_{IJK}^{\rm pol} = E_{IJK}^{\rm pol} - E_{IJ}^{\rm pol}
                            - E_{IK}^{\rm pol} - E_{JK}^{\rm pol},

the non-additive remainder. Nothing in the expansion was ever two-body
specific: no group's Hamiltonian depends on its environment, so the differences
telescope exactly. GAMESS stops at two because ``efmo.src`` enumerates the
levels by hand, not because the method does.

**A group is near only if every pair inside it is near.** Not a convention: the
far half of the energy is pairwise by construction -- electrostatics, exchange
repulsion, dispersion and charge transfer are all two-body, and the induction is
already all orders in :math:`E_{\rm pol}^{\rm total}` -- so there is **no far
n-body term**, and a group holding a far pair would count that pair twice, once
inside its own SCF and once in the effective-fragment sum. The criterion is also
inherited by subsets, which is what makes the filtered group list one the
difference operator can be applied to at all: every subset of an all-near group
is all-near, so every group's share is actually subtracted.

**Level = N with a huge cutoff is the unfragmented calculation.** Two series
telescope at once: the in-vacuo one to the supersystem's own SCF energy, and the
induction one to :math:`E_{\rm pol}^{\rm total}`, which therefore cancels the
last term of the expression **entirely**. Nothing is left of the polarization
correction and the answer is one RHF energy. That is the sharpest available check
on the subset-induction bookkeeping -- a wrong sign, a missing group or a group
induction solved over the wrong fragments all survive the level-two limits and
break it -- and it is asserted in ``test/test_mqc_czt_efmo.f90`` and pinned as a
validation case.

**The cost is the binomial and nothing else**, so a raised level is not a small
change: there are :math:`\binom{N}{n}` groups of size :math:`n` before the
cutoff thins them, which on twenty fragments is up to 1140 SCFs at level three
against 190 at level two. No level is refused, and a run at level three or above
logs both the enumerated count and the count the cutoff kept before it computes
any of them.

Measured on the water prism, RHF/6-31G, :math:`R_{\rm cut} = 2.0` -- where the
cutoff keeps every pair, so each level is the *complete* expansion at that order
-- against the unfragmented RHF energy of the same cluster,
-456.007613510775 Hartree:

.. list-table::
   :header-rows: 1
   :widths: 10 14 26 26 24

   * - level
     - near groups
     - total / Hartree
     - error / Hartree
     - error per fragment
   * - 2
     - 15
     - -456.005523966192
     - 2.09e-3
     - 0.219 kcal/mol
   * - 3
     - 35
     - -456.007038360748
     - 5.75e-4
     - 0.060 kcal/mol
   * - 4
     - 50
     - -456.007692948772
     - -7.94e-5
     - -0.008 kcal/mol
   * - 5
     - 56
     - -456.007611814554
     - 1.70e-6
     - 0.0002 kcal/mol
   * - 6
     - 57
     - -456.007613510751
     - 2.4e-11
     - 0
   * - unfragmented
     - --
     - -456.007613510775
     - --
     - --

Level six is the fragment count, so the expansion is exact and the 2.4e-11
residual is accumulated SCF convergence over 57 groups and nothing else.

**A raised level is not unconditionally better, and one case here shows it.**
The same prism at :math:`R_{\rm cut} = 1.0` gives 0.359 kcal/mol per fragment at
level 2 and 0.583 at level 3 -- *worse*. The reason is the near criterion, not
the expansion: at that cutoff six of the fifteen pairs are effective, so only two
of the twenty triples have all three of their pairs near and the three-body
correction that gets added is two terms out of twenty. An expansion truncated
*and* filtered is not variational in the level. Raise the cutoff along with the
level, or read the level as a refinement of the quantum region rather than of the
whole cluster.

Running one
-----------

.. code-block:: json

   {
     "schema": {"name": "efmo_prism", "version": "1.0"},
     "molecules": [{
       "xyz": "prism.xyz",
       "fragments": [[0,1,2],[3,4,5],[6,7,8],[9,10,11],[12,13,14],[15,16,17]],
       "fragment_charges": [0, 0, 0, 0, 0, 0],
       "fragment_multiplicities": [1, 1, 1, 1, 1, 1],
       "molecular_charge": 0,
       "molecular_multiplicity": 1
     }],
     "model": {"method": "hf", "basis": "6-31g"},
     "keywords": {
       "fragmentation": {"method": "efmo", "level": 2, "rcut": 2.0},
       "efmo": {"charge_transfer": true, "induction_damping": 0.6}
     },
     "driver": "Energy"
   }

.. code-block:: bash

   ./mqc efmo_prism.json

Keywords
--------

.. list-table::
   :header-rows: 1
   :widths: 30 12 58

   * - Key
     - Default
     - Meaning
   * - ``keywords.fragmentation.method``
     - --
     - ``"efmo"`` selects this method. Required.
   * - ``keywords.fragmentation.level``
     - ``2``
     - How far the many-body expansion of the near groups runs; see *Beyond
       pairs*. The same key MBE and FMO read, and defaulted to 2 here rather
       than to the shared default of 1, which for EFMO would drop every near
       pair. A level below 1 is refused; a high one is warned about and run.
   * - ``keywords.fragmentation.rcut``
     - ``2.0``
     - :math:`R_{\rm cut}`, unitless. Sits here rather than under ``efmo``
       because it decides which pairs are solved quantum mechanically, which is
       a property of the partition -- the same place FMO's ``resppc`` lives.
   * - ``keywords.fragmentation.bond_breaking``
     - ``"none"``
     - ``"none"`` refuses a partition that cuts a covalent bond, naming the two
       atoms. ``"afo"`` detaches each cut bond with an adjusted frozen orbital;
       see *Cutting covalent bonds*. The same key FMO reads.
   * - ``keywords.efmo.charge_transfer``
     - ``true``
     - Include :math:`E_{IJ}^{\rm CT}` in the far pairs. GAMESS's library default
       (``$FMO MODEFM(4)=0``) leaves it out, but every EFMO deck shipped in the
       GAMESS tree turns it on, and so does the 2024 paper; the original 2012
       method used electrostatics alone. On the water prism at
       :math:`R_{\rm cut} = 0.3` it is -0.0113 Hartree, so leaving it out is a
       choice worth making deliberately -- hence ``true`` here.
   * - ``keywords.efmo.induction_damping``
     - ``0.0``
     - :math:`a` of the Tang-Toennies-like factor
       :math:`1 - e^{-aR^2}(1 + aR^2)` that damps every induction field between
       two fragments -- the static field of the permanent multipoles and the
       field of the other induced dipoles alike, and :math:`E_{IJ}^{\rm pol}`
       and :math:`E_{\rm pol}^{\rm total}` alike, since eq 6 subtracts one
       from the other. Zero is undamped. GAMESS's EFMO runs 0.6 for a cluster
       of whole molecules and 0.1 where a fragment was cut across a bond
       (``POLAB``, set from ``$FMO SCREEN``); above 2.0 the factor is one
       again, which is GAMESS's own guard. **The default is off rather than
       0.6** -- see the induction paragraph under *Against GAMESS*, which
       measures what the key does and what it does not. **With detached bonds
       zero means 0.1**, GAMESS's value for that case, because undamped the
       induction across a cut does not converge.
   * - ``model.method``
     - --
     - ``"hf"``, ``"mp2"`` or ``"ri-mp2"``. The correlation runs on the
       orbitals each monomer and each quantum dimer already converged to, so
       :math:`E_I^0` and :math:`E_{IJ}^0` become correlated energies and
       nothing else in eq 6 moves: the fragment potentials, the far pairs and
       the induction are Hartree-Fock constructions, MAKEFP being one. A
       Kohn-Sham or coupled-cluster method is refused by name.
   * - ``model.aux_basis``
     - --
     - The correlation-fitting (RIFIT) set ``"ri-mp2"`` needs. Required with
       that method and unread otherwise.
   * - ``keywords.correlation.freeze_core``
     - ``true``
     - Whether each fragment's core orbitals sit out the MP2. The count is
       derived per fragment from its elements, so a dimer's core is the sum of
       its two monomers' and :math:`E_{IJ} - E_I - E_J` differences the same
       set of correlated orbitals on both sides.
   * - ``keywords.efp.*``
     - --
     - The MAKEFP settings -- the response solve and the screening grid -- passed
       to every fragment's potential. The same keys a ``MakeFP`` run uses; see
       :doc:`makefp`.
   * - ``keywords.scf.*``
     - --
     - How every SCF here is *driven*: accelerator, DIIS subspace, level shift,
       linear-dependence threshold, incremental Fock. Its tolerances are **not**
       read -- see below.

Every SCF in an EFMO run, monomer and dimer alike, is converged to
:math:`10^{-10}` in energy and :math:`10^{-8}` in density and orbital gradient,
which are ``make_efp_potential``'s own defaults. That is deliberately tighter than
a whole-system run and is not settable: the near-dimer correction is
:math:`E_{IJ}^0 - E_I^0 - E_J^0`, four orders smaller than any of the three, so a
looser convergence leaves it with no significant figures. A looser EFMO would not
be a cheaper one either -- the cost is MAKEFP.

Correlated fragments
--------------------

``model.method: "ri-mp2"`` (or ``"mp2"``) runs the correlation on every
monomer's and every quantum dimer's converged orbitals, which is the
EFMO/RI-MP2 of the 2024 paper. The monomer's MP2 uses the SCF
``make_efp_potential`` already ran -- there is no second determinant -- and the
dimer's uses its own. What comes back is eq 6 with correlated :math:`E_I^0` and
:math:`E_{IJ}^0` and every other term unchanged:

.. code-block:: json

   "model": {"method": "ri-mp2", "basis": "6-31g", "aux_basis": "cc-pvdz-rifit"}

Two identities pin it, both exact rather than approximate and both in
``test/test_mqc_czt_efmo.f90``: on two fragments EFMO/RI-MP2 is the dimer's own
in-vacuo RI-MP2 energy, and with every pair quantum it is the RI-MP2 many-body
pair sum plus the induction no pair holds. Switching the correlation off
reproduces the Hartree-Fock total exactly.

The run reports how much of the fragment sum and of the dimer correction is
correlation. Those are reported *inside* the two sums and not beside them: a
correlated :math:`E_I^0` is the monomer energy of eq 6, not a term added to it.

Running it on several ranks
---------------------------

.. code-block:: bash

   mpirun -np 4 ./mqc efmo_prism.json

The monomers and the quantum groups are handed out round robin -- one flat task
list whatever the level, so the trimers and above distribute exactly as the pairs
do. **The balance
is struck on the monomers**, because a monomer is a MAKEFP -- an SCF, a
localization and twelve frequency-dependent response solves -- against one SCF
for a dimer. Every rank then needs every fragment's potential, since the far
pairs and the one induction over all fragments are not decomposable by owner,
so each potential is flattened into a pair of buffers and summed across ranks.
That transfer is exact rather than nearly: what crosses is the bits, not the
eight decimals a written ``.efp`` would carry.

The far pairs and the induction stay replicated. They are milliseconds beside a
potential, and replicating them means every rank reaches the same total without
a second reduction -- so any rank could write the output file, and the leader
does.

**One rank and four are bit-identical**, measured on the water prism at
:math:`R_{\rm cut}` 1.0 and 2.0 and on the EFMO/RI-MP2 trimer, at one thread.
Across thread counts the usual OpenMP reduction-order scatter of about
2 :math:`\times` 10\ :sup:`-12` applies, and it is a thread effect and not a
rank one.

Output
------

The log carries a table of every :math:`E_I^0`, every pair with its
:math:`R_{IJ}` and its class, one row per many-body level, and the eight sums. The JSON output repeats the sums
under ``efmo``, with the pair counts:

.. code-block:: json

   "efmo": {
     "monomer_sum": -455.897884059762,
     "qm_nmer_correction": -0.094501684825,
     "induction_correction": -0.013350304484,
     "efp_electrostatics": 0.0,
     "efp_dispersion": 0.0,
     "efp_exchange_repulsion": 0.0,
     "efp_charge_transfer": 0.0,
     "polarization_total": -0.026488526089,
     "qm_dimers": 15,
     "efp_dimers": 0,
     "qm_groups": 15,
     "fragment_charges": [0, 0, 0, 0, 0, 0]
   }

``qm_nmer_correction`` is :math:`\sum_S \Delta E_S^0` over the near groups of two
or more fragments -- at level two the dimer corrections. ``induction_correction``
is :math:`\sum_S \Delta E_S^{\rm pol}`, and is reported with the sign it has as a
sum while being **subtracted** from the total, so that it can be compared against
another code's pair induction directly. ``qm_groups`` counts the near groups and
equals ``qm_dimers`` at level two; the log carries the same sums split per level,
which is what says whether a raised level was worth its binomial.

The per-pair interaction map
----------------------------

The sums above say what the halves came to. They cannot say **which pairs** the
interaction is made of, and two systems with identical sums can have entirely
different pairs carrying them. ``efmo.pairs`` is that map, one object per pair.
The two below are the strongest quantum pair and the strongest classical one from
the shipped ``efmo_prism_rcut1`` deck:

.. code-block:: json

   "pairs": [
     {
       "fragments": [5, 6],
       "distance": 0.770,
       "treatment": "quantum",
       "interaction_energy": -0.010522652
     },
     {
       "fragments": [1, 4],
       "distance": 1.313,
       "treatment": "classical",
       "interaction_energy": -0.003076750,
       "electrostatics": -0.002929635,
       "dispersion": -0.000093566,
       "exchange_repulsion": -0.000052604,
       "charge_transfer": -0.000000944
     }
   ]

``fragments`` are numbered **from one**, as the log's table and the MBE fragment
lists number them. (Atom indices elsewhere in the interfaces are 0-based; these
are not atoms.) ``distance`` is the vdW-scaled :math:`R_{IJ}` that decided the
split, so a reader can see which side of :math:`R_{\rm cut}` a pair fell on.

``interaction_energy`` is what that pair contributed to the total: for a quantum
pair :math:`E_{IJ}^0 - E_I^0 - E_J^0 - E_{IJ}^{\rm pol}`, and for a far pair the
sum of its four terms. It is the same number the log's table prints, from the
same routine, so the two cannot drift apart.

**Rows are sorted strongest first**, by magnitude. Nothing else in the JSON
output sorts; this does because the enumeration order puts the largest terms
nowhere in particular, and an interaction map at protein scale is read from the
top.

``fragment_charges`` sits beside the pairs, one per fragment in fragment order.
A pair row cannot explain its own large electrostatics term, and two adjacent
charged fragments interact through a monopole term that dwarfs everything else
in the table. An EFMO fragment's charge comes from the deck and from nowhere
else -- EFMO runs no embedding of any kind -- so nothing about a run's settings
can move these.

Three things the key names do not say
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**The four named terms appear only on a classical pair.** A quantum pair's
dimer SCF produces no such split, and four zeros would read as a decomposition
that had been computed.

**``charge_transfer`` is not a residual.** It is the explicit perturbative term
-- GAMESS's ``ECHTR``, a sum over occupied-on-A to virtual-on-B amplitudes --
computed independently of the other three. A decomposition analysis usually
means by that name what is *left over* after the others, which additionally
absorbs basis-set superposition error; this does not, and the two should not be
compared as though they were the same quantity.

**No per-pair polarization is reported.** The induction is solved
self-consistently over the whole system and has no per-pair value. A quantum
pair's ``interaction_energy`` does carry its own :math:`E_{IJ}^{\rm pol}`
subtracted, because that is how the energy expression above is written -- but
that is a term of this method, not a polarization component of a decomposition,
and relabelling it as one would be inventing a number.

What the map sums to
~~~~~~~~~~~~~~~~~~~~

**At level two** the pairs are the whole of the near and far interaction, so

.. math::

   \sum_{\rm pairs} E_{\rm int} =
     \texttt{qm\_nmer\_correction} - \texttt{induction\_correction}
     + \text{the four efp sums}

exactly, and the reported total is that plus ``monomer_sum`` and
``polarization_total``. Measured on the shipped decks: 1.4e-17 Hartree on
``efmo_prism_rcut1`` and exactly zero on ``efmo_cage_rcut2``.

**Above level two it does not**, and cannot. The map holds pairs; the near
expansion at level three also holds groups of three, which are not pairs and are
not in the list. On ``efmo_w3_level3`` the sum falls short by 1.26e-4 Hartree,
which is exactly the three-body vacuum term less its induction. The pair map is
still every pair -- it is simply no longer the whole interaction.

Against GAMESS
--------------

GAMESS runs this method as ``$FMO IEFMO=1``, and its ``RESDIM`` is the same
:math:`R_{\rm cut}`. Measured on (H\ :sub:`2`\ O)\ :sub:`3`, (H\ :sub:`2`\ O)\
:sub:`4` and five water hexamers, at :math:`R_{\rm cut}` 2.0, 1.0 and 0.3, RHF/6-31G:

* **The split is identical** in all twenty-one cases -- the same pairs quantum,
  the same pairs effective.
* **The fragment sum and the quantum dimer corrections agree to 2e-8 Hartree**,
  which is GAMESS's printed precision. So does **exchange repulsion**, to 5e-8.
* **Dispersion agrees to 0.7 per cent** once GAMESS is asked for the same model
  (``MODEFM(3)=32``, damped :math:`E_6+E_7+E_8`); its default ``MODEFM(3)=1`` is
  :math:`E_6 + \tfrac{1}{3}E_6`, a different quantity.

Two terms differ on purpose, and neither is a disagreement about the same number:

* **Electrostatics.** GAMESS's EFMO sums point multipoles with no
  charge-penetration screening (``MODEFM(1)`` has no screening bit set by
  default, and its ``SCREEN(1)=-1`` alternative is marked experimental and is not
  symmetric between equivalent pairs). This code applies the screening the
  potential itself carries, which is what an ordinary EFP2 run here does. Turning
  it off for EFMO alone would break the :math:`R_{\rm cut} \to 0` limit, where the
  energy must equal this program's own EFP -- and it is the accuracy the
  screening exists for. On the prism at :math:`R_{\rm cut} = 0.3` the two Coulomb
  sums are -0.1407 against -0.1323 Hartree.
* **Induction.** GAMESS damps the induction field with a Tang-Toennies-like
  factor at :math:`a = 0.6` for a molecular cluster (``PENSAB`` in ``FRGFLD``
  for the static field, ``P1`` in ``DIPIT`` for the induced-dipole field, both
  in ``efintb.src``), and by default this code does not, so its induction runs
  two to four per cent deeper. The *pair* induction :math:`E_{IJ}^{\rm pol}`
  and the total :math:`E_{\rm pol}^{\rm total}` are the same routine in both
  codes -- same screening, same self-consistent solve, same convergence --
  differing only in how many fragments are loaded, so the subtraction is clean
  on both sides.

  ``keywords.efmo.induction_damping`` applies exactly that factor, and running
  it settles what the difference is made of. On the prism at
  :math:`R_{\rm cut} = 1.0`, in Hartree:

  .. list-table::
     :header-rows: 1

     * -
       - undamped
       - ``induction_damping: 0.6``
       - GAMESS
     * - :math:`\sum E_{IJ}^{\rm pol}`
       - 0.013033887
       - 0.012374267
       - 0.012557400
     * - :math:`E_{\rm pol}^{\rm total}`
       - -0.026488526
       - -0.025407237
       - -0.025869694
     * - total
       - -456.004182669
       - -456.003761000
       - -456.004043215

  So the damping is real and it *overshoots*: undamped we sit 2.4 per cent
  deeper than GAMESS and damped 1.8 per cent shallower. About a third of the
  induction gap is the damping and the rest is a difference in the undamped
  induction itself, which the Tang-Toennies factor cannot be blamed for.
  Damping the induced-dipole field alone moves the total by 5e-6 -- the static
  field carries all of it. The default is therefore left undamped, which is
  also what every reference in this repository was pinned with; set the key
  when the point is to reproduce a GAMESS induction rather than to be right.

Accuracy
--------

The eighteen water hexamer isomers in ``validation/inputs/sample_inputs``, against
the unfragmented RHF/6-31G of each, per fragment:

* :math:`R_{\rm cut} = 2.0`: 0.04 to 0.68 kcal/mol
* :math:`R_{\rm cut} = 1.0`: 0.14 to 0.68 kcal/mol
* :math:`R_{\rm cut} = 0.3`: -1.7 to -2.7 kcal/mol

So the paper's target of under 1 kcal/mol per fragment is met at 2.0 and 1.0 and
missed at 0.3, where no dimer is quantum at all.

**Relative energies are the harder test and it is not passed.** The error is not
constant across isomers -- it tracks how compact the cluster is, being smallest
for the prism and largest for the ring-like chair and boat -- so a relative energy
carries up to 2.5 kcal/mol at :math:`R_{\rm cut} = 2.0`. The most stable isomer
(chair) is picked correctly at every cutoff, but the ordering of the next few is
not preserved: the prism rises from fifth to second. Read EFMO energies of
different isomers against each other with that in mind.

Cutting covalent bonds
----------------------

``keywords.fragmentation.bond_breaking: "afo"`` lets a partition cut covalent
bonds, a protein by residue being the case it exists for:

.. code-block:: json

   "keywords": {
     "fragmentation": {"method": "efmo", "level": 2, "rcut": 2.0,
                       "bond_breaking": "afo"}
   }

Each cut bond is detached exactly as FMO detaches it, by FMO's own code, which
is GAMESS's adjusted frozen orbitals: the bond orbital of a small model system
is frozen empty in the fragment owning the detached atom and frozen occupied in
the fragment across the bond, which carries the atom's functions as a ghost with
the atom's core and other bonds frozen empty. **The nucleus is always
split**, ``Z - 1`` with the owner and ``+1`` on the ghost, so every fragment is a
neutral closed shell; GAMESS splits too. The cut atoms are perceived from the
geometry; the detached end of each is its sp3 end, else its lower-numbered
atom, unless ``keywords.fragmentation.detached_atoms`` names it.

Every SCF the method runs sees the cuts: each monomer's, which is also the one
its potential is made from, and each near group's. **A cut belongs to a group,
not a fragment** -- a group holding both ends of a bond has it whole, with no
ghost, no frozen orbital and both halves of the nucleus back together.

Three rules the cuts add:

* **A covalently joined pair is always quantum**, whatever ``rcut`` says.
  Joined means a cut runs between the two, or both hold centres at one cut atom
  or a bond apart -- a ghost sits a bond length from every fragment its atom is
  bonded to. No multipole expansion describes a pair that close: propane in
  three pieces, its two ends joined only through a ghost of the middle carbon,
  gave -0.336 Hartree effective against +0.274 quantum before this rule. Along a
  peptide backbone cut at C\ :sub:`alpha`--C it joins exactly the neighbouring
  residues.
* **The split is decided over the ghosts too**, since a ghost carries a charge
  and a bond pair where it sits. The reported ``distance`` stays atom to atom.
* **The induction is damped at a = 0.1** when ``keywords.efmo.induction_damping``
  is left at zero -- GAMESS's own value for a run with detached bonds. Undamped,
  the induced dipoles of two fragments sharing a bond region do not converge at
  all. A value set in the deck is used as given.

Refused by name: MP2 and RI-MP2 across a cut (the frozen virtual would be
correlated like any other); everything FMO's detachment refuses -- a cut through
a ring, a bond detached at a hydrogen, a bond carrying more than one localized
orbital; and gradients, for EFMO as a whole.

What it gives
~~~~~~~~~~~~~

**At level equal to the fragment count, with rcut huge, the total is the
molecule's own RHF energy**, whatever was cut: 1.8e-13 Hartree on propane cut
once, 3.1e-13 cut twice, 9.9e-14 numbered so the middle carbon is detached twice,
1.5e-12 on butane in four pieces and 1.1e-13 on methylcyclopropane cut at its
exocyclic bond, RHF/6-31G, in ``test/test_mqc_czt_efmo_covalent.f90``. Each
fragment's potential sums to its charge to about 1e-14.

**Far pairs of a cut fragment are as good as those of a whole molecule** from a
few angstrom out. Propane with a water beyond its first methyl, the effective
pair against the same pair solved as a quantum dimer:

=========  ================  ================
O to C     uncut, error      cut, error
=========  ================  ================
3.5 A      -7.1e-4           -2.3e-3
4.5 A      +3.2e-5           -2.3e-5
6.0 A      -9.3e-6           -6.0e-6
=========  ================  ================

At 3.5 A the pair is quantum at the default ``rcut`` anyway. **Pairs two cuts
apart are poor at short range**: butane's two end carbons, at
:math:`R_{IJ} = 0.85`, give +0.0081 effective against +0.0519 quantum. Keep
``rcut`` at 1 or above with cuts.

**Against GAMESS.** Butane cut at its middle C-C bond with a water 4.5 A beyond
one end, RHF/STO-3G, level 2, against the unfragmented -230.415265709:

======================  =================  =================
                        this code          GAMESS
======================  =================  =================
every pair quantum      -230.415258725     -230.415264851
rcut 1.0                -230.415278030     -230.415231748
error, every pair QM    +7.0e-6            +8.6e-7
======================  =================  =================

The fragments agree: the three monomer energies to 5e-9, 6e-9 and 5e-8 and the
quantum pair corrections, summed, to 1e-8. What is left is induction. The two
codes' pair and total induction differ by about a third on this system --
-0.00152 here against -0.00218 in GAMESS for the total -- and almost all of it
cancels between the two, leaving the 6e-6. At ``rcut`` 1.0 the effective pairs
add GAMESS's unscreened electrostatics against this code's screened ones, the
difference noted for whole molecules under *Against GAMESS* above.

Before the frozen set was matched to GAMESS's (see :doc:`fmo`), the same row
read -230.416556977, an error of -1.29e-3: the two cut monomers came out 0.127
and 0.298 Hartree high.

How the potential of a cut fragment is made
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A fragment cut across a covalent bond is detached the way FMO detaches it (see
:doc:`fmo`): the detached atom keeps ``Z - 1`` of its nucleus and the bond
orbital frozen empty, and the fragment across the bond carries a centre at the
same point with the same basis, ``+1`` of charge, the bond orbital frozen
occupied and the atom's other orbitals frozen empty.
Both fragments are neutral closed shells, which is what a fragment *potential*
needs -- nothing in a far pair takes a net charge back out.

That makes the charge a centre presents and the element it is two different
things, and the effective-fragment path used to read one off the other. What
reads which, now:

.. list-table::
   :header-rows: 1
   :widths: 34 33 33

   * - Where
     - Element (``atomic_numbers``, ``element``)
     - Charge presented (``nuclear_charge``, ``charge``)
   * - ``make_efp_potential``: the molecule
     - the basis looked up
     - the nuclear attraction and repulsion; the electron count is their sum
   * - core orbitals left out of the localization
     - the frozen core of the element; **none on a ghost**, whose core stays
       with the atom's owner
     - --
   * - ``distributed_multipoles``
     - bonds, so the midpoints, and the labels
     - the nuclear monopole
   * - screening grid, masses for the centre of mass
     - radii and masses
     - --
   * - projection basis header, exchange-repulsion potential
     - --
     - the *valence* charge, charge less two per core orbital: ``+3`` on a
       detached carbon, ``+1`` on the centre across the cut
   * - a fragment's basis in the pair terms
     - the element symbol
     - the charge the molecule is built with (overlap and kinetic only)

A potential made in memory carries ``element`` and ``valence`` per atom through
to the fragment and across ranks. A ``.efp`` file carries neither, and a
fragment read from one falls back on the nearest integer to its charge, which
for a whole molecule is the same number -- so nothing read from a file moves.

**The centre across the cut is an ordinary atom of the potential.** It has
basis functions, so the distributed multipoles need an expansion point there;
it is bonded to the attached atom, so it gets a midpoint; and it carries ``+1``.
Summed over the two fragments of a cut, the monopoles put the whole nucleus back
exactly once.

The frozen occupied bond orbital is localized, polarized and written with the rest of
the valence orbitals, on the constrained determinant as it stands. That is what
GAMESS's EFMO does, including its own acknowledgement that the Fock matrix is
not diagonal in the unfrozen space.

Every quantum mechanical step an EFMO run takes, with or without cuts:

* each fragment's SCF, inside ``make_efp_potential``, which is also its
  :math:`E_I^0`;
* each near group's SCF;
* with a correlated method, MP2 on both;
* with cuts, one small SCF per cut bond, on the model system its frozen orbitals come
  from.

All of them go through the same group assembly: a cut belongs to a group, not a
fragment, so a group holding both ends of a bond has it whole.

What is not here yet
--------------------

* **Kohn-Sham fragments.** MAKEFP is a Hartree-Fock construction, so a
  density-functional :math:`E_I^0` would need a potential built from a
  different density; refused by name. So is coupled cluster, and so is
  spin-component-scaled MP2, which would be a different method from the one
  the paper runs.
* **Mutually polarized QM/EFP embedding.** A separate method rather than a
  missing feature: it puts the induction inside each group's SCF, which makes a
  group's energy depend on its environment and stops the many-body differences
  telescoping. Nothing in EFMO as implemented here embeds anything.
* **Gradients**, with or without cuts; refused by name. Across a cut they
  need the frozen orbitals' own response, since the model systems' caps move
  with the atoms.
* **Correlated fragments across a cut**; the frozen empties would be
  correlated like any other virtual.
* **The rest of the induction difference.** With
  ``keywords.efmo.induction_damping`` set to GAMESS's 0.6 the two codes'
  induction still differ by about two per cent, in the other direction; what is
  left is a difference in the *undamped* induction and it is not yet accounted
  for.
