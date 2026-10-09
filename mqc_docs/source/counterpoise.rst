Counterpoise correction (VMFC and SSFC)
=======================================

A fragment solved in its own basis is described worse than the same fragment
solved inside a larger one, because in the larger one it can borrow its
neighbour's basis functions. That difference is **basis-set superposition
error**, and in a many-body expansion it does not cancel: the pair term
subtracts monomers described in a small basis from a pair described in a large
one, so the pair looks more bound than it is. The error survives truncation, and
each higher order adds its own rather than cancelling the last.

Counterpoise removes it by putting both sides of every difference in the same
basis, so the borrowing appears on each. Two schemes differ in which basis that
is:

* ``vmfc``, the Valiron--Mayer function counterpoise, solves every subfragment
  in its **parent's** basis: a monomer inside a dimer in the dimer's basis, a
  dimer inside a trimer in the trimer's.
* ``ssfc``, the site-site function counterpoise (Wells and Wilson), solves
  every term of order two and above in the basis of the **whole system**: all
  :math:`N` monomers' functions, the monomers that are not in the term present
  as ghosts. It is the Boys--Bernardi recipe applied to every term.

Both are available for energies on the MBE, and neither changes a deck except
by the one keyword.

Running one
-----------

.. code-block:: json

   "keywords": {
     "fragmentation": {
       "method": "MBE",
       "level": 2,
       "counterpoise": "vmfc"
     }
   }

``"counterpoise"`` takes ``"vmfc"``, ``"ssfc"`` or ``"none"``, and defaults to
``"none"`` -- which is the uncorrected expansion every deck written before this
existed already asked for. Nothing else in the deck changes. For the other
scheme, ``"counterpoise": "ssfc"``:

.. code-block:: json

   "keywords": {
     "fragmentation": {
       "method": "MBE",
       "level": 3,
       "counterpoise": "ssfc"
     }
   }

The run log names the scheme once the term list has been checked against it
(``Counterpoise: ssfc. Every term of order 2 and above is solved in the basis of
all 3 monomers; ...``), and the JSON output carries ``"counterpoise": "ssfc"``
beside the fingerprint. A run without counterpoise writes neither.

What VMFC is worth
------------------

With ``vmfc``, the water dimer ``validation/inputs/sample_inputs/w2_dimer.xyz`` at
Hartree-Fock/6-31G, a basis small enough that the borrowing is large:

============================  =================  ==================
Quantity                      Hartree            kcal/mol
============================  =================  ==================
Binding energy, uncorrected   -0.0102112369      -6.408
Binding energy, VMFC          -0.0080445970      -5.048
BSSE removed                   0.0021666399       1.360
============================  =================  ==================

Twenty-one per cent of the binding was basis set. The trimer ``w3.xyz`` --
three waters stacked 2.9 Angstrom apart -- at ``"level": 3``:

===========  =================  =================  ==================
Term         Uncorrected (Ha)   VMFC (Ha)          Difference
===========  =================  =================  ==================
1-body       -227.9519233972    -227.9519233972    none, by definition
2-body         -0.0178176735      -0.0145142064    2.07 kcal/mol
3-body         -0.0007565683      -0.0006102466    0.09 kcal/mol
Total        -227.9704976390    -227.9670478502    2.16 kcal/mol
===========  =================  =================  ==================

Two things in that table are worth reading. The uncorrected total,
``-227.9704976390``, is exactly the supermolecular energy -- an MBE at level
``N`` over ``N`` fragments is exact, and this one is exact to every digit
printed, which is how you know the expansion itself is right and the difference
is entirely superposition error. And the error is not spread evenly: nearly all
of it, 2.07 of 2.16 kcal/mol, sits in the two-body term. The three-body term
moves by 0.09 kcal/mol, a fifth of its own size, and keeps its sign. That is
what Valiron--Mayer is built to do: every subset a trimer subtracts is computed
in the trimer's own basis, so the three-body correction subtracts like from
like.

The one-body term is untouched on purpose. Each monomer keeps its own basis
there, because that is the reference the interaction energy is measured *from*.
Only the corrections are ghosted.

The full-cluster basis (SSFC)
-----------------------------

Write :math:`E_T(\mathrm{FB})` for the energy of the monomers in :math:`T` in the
**full basis** -- the functions of all :math:`N` monomers, those outside
:math:`T` present as ghosts -- and :math:`E_i(i)` for monomer :math:`i` in its own
basis. ``ssfc`` keeps the one-body term in each monomer's own basis and puts
every correction of order two and above in the full basis:

.. math::

   E_{\mathrm{SSFC}} = \sum_i E_i(i) + \sum_{|T| \ge 2} \Delta E_T(\mathrm{FB}),
   \qquad
   \Delta E_T(\mathrm{FB}) = E_T(\mathrm{FB})
      - \sum_{S \subsetneq T} \Delta E_S(\mathrm{FB}),
   \qquad
   \Delta E_i(\mathrm{FB}) = E_i(\mathrm{FB}).

The recursion is the ordinary one; only the basis each energy in it is computed
in has changed.

**The identity.** Carried to level :math:`N` the sum telescopes:

.. math::

   E_{\mathrm{SSFC}}(N) = E_{\mathrm{cluster}}
      + \sum_i \left[ E_i(i) - E_i(\mathrm{FB}) \right],

the energy of the whole cluster plus the superposition error of each monomer,
which is the Boys--Bernardi counterpoise-corrected total. For a variational
method, Hartree-Fock or Kohn-Sham, each bracket is positive, because a monomer
is lower in the cluster's basis than in its own, so the correction raises the
total. Below level :math:`N` the same sum is truncated
as any many-body expansion is, with every order in the same basis.

The water trimer of the last section, ``"level": 3``, in all three forms:

===========  =================  =================  =================  ==================
Term         Uncorrected (Ha)   VMFC (Ha)          SSFC (Ha)          SSFC - uncorrected
===========  =================  =================  =================  ==================
1-body       -227.9519233972    -227.9519233972    -227.9519233972    none, by definition
2-body         -0.0178176735      -0.0145142064      -0.0146751670    1.97 kcal/mol
3-body         -0.0007565683      -0.0006102466      -0.0006102466    0.09 kcal/mol
Total        -227.9704976390    -227.9670478502    -227.9672088108    2.06 kcal/mol
===========  =================  =================  =================  ==================

The SSFC total is the uncorrected supermolecular energy plus the three
monomers' superposition errors, and ``test/test_mqc_counterpoise.f90`` checks
that identity directly. The two schemes disagree in the pair term, whose
monomers VMFC takes in a pair's basis and SSFC in the trimer's, and agree in the
three-body term, which is the same trimer-basis calculation in both. At
``"level": 2`` the totals are ``-227.9664376036`` (VMFC) and ``-227.9665985642``
(SSFC).

Choosing between them
^^^^^^^^^^^^^^^^^^^^^

* To reproduce a full-cluster-basis counterpoise, **use** ``ssfc``. That is how
  a cluster is usually corrected by hand: one basis, the cluster's, for every
  term. Nguyen and Xantheas (J. Chem. Theory Comput., 2026) correct their
  many-body terms this way, and their terms for :math:`n \ge 2` are the
  ``levels[].total_energy`` of this scheme, so the tables compare directly.
* For a large cluster at a level well below :math:`N`, **use** ``vmfc``. It is
  consistent order by order -- each n-body term corrected in its own n-mer's
  basis -- and every calculation is in a basis no larger than the n-mer's. It
  pays in rows (below), but the rows are small. Near :math:`n = N` that reverses:
  the :math:`2^n - 1` rows of each n-mer sit in nearly the full basis. On the
  hexamer assumptions below, the CCSD cost against one unfragmented run is 7.3
  (SSFC) and 0.89 (VMFC) at :math:`L = 3`, then 14.2 against 7.7 at
  :math:`L = 4`, 18.5 against 27.8 at :math:`L = 5` and 19.5 against 47.3 at
  :math:`L = 6`.
* **Use either for two monomers.** At :math:`N = 2` the schemes coincide: the
  term lists are the same five rows and the energies are identical. A test
  asserts it.
* **The identity above is SSFC's alone.** It does not hold for VMFC at
  :math:`N \ge 3`, because a VMFC pair term stays in a pair's basis whatever the
  size of the cluster. The VMFC total at level :math:`N` therefore is not the
  Boys--Bernardi total, and the SSFC total is.

One difference of convention from Nguyen and Xantheas: their one-body reference
is also in the full basis, so their total at :math:`L = N` is the uncorrected
cluster energy. Here the one-body term is in each monomer's own basis -- the
reference an interaction energy is measured from, as it is in ``vmfc`` -- and
the total at :math:`L = N` is the counterpoise-corrected one. Their terms for
:math:`n \ge 2` are equal to the ones here; only the one-body line, and so the
total, differs by :math:`\sum_i [E_i(i) - E_i(\mathrm{FB})]`.

The many-body counterpoise of Richard, Lao and Herbert (MBCP(:math:`n`)) agrees
with both schemes at :math:`N = 2` and differs at a truncated order; it is not
implemented.

What it costs
-------------

Every n-mer contributes its :math:`2^n - 2` proper subsets to ``vmfc`` as extra
terms, so the row count goes from :math:`\sum_{n=1}^{L} \binom{N}{n}` to
:math:`N + \sum_{n=2}^{L} \binom{N}{n}(2^n - 1)`. ``ssfc`` adds one row per monomer and
nothing else: each kept term has a single row, in the full basis, and the
own-basis rows of the pairs and above are dropped because nothing reads them.
The count is :math:`N + \sum_{n=1}^{L} \binom{N}{n}`, the uncorrected count plus
:math:`N`, and the same at :math:`L = N`, where the whole system is one all-real
row instead of a ghosted one:

=========  =========  ============  =============  ============
Fragments  Level      Uncorrected   VMFC           SSFC
=========  =========  ============  =============  ============
20         2                   210            590           230
20         3                  1350           8570          1370
20         4                  6195          81245          6215
50         2                  1275           3725          1325
50         3                 20875         140925         20925
50         4                251175        3595425        251225
=========  =========  ============  =============  ============

VMFC multiplies the rows by about 2.9, 6.5 and 14 at levels 2, 3 and 4. Rows are
only half of what a calculation costs, and where the two schemes differ is the
other half: the basis. A VMFC row is solved in a basis of at most the n-mer's,
and a ghosted monomer of a trimer carries the trimer's full basis. An SSFC row
is solved in the basis of **all** :math:`N` monomers, however few of them it
holds. SSFC has the fewer rows and the larger calculations, and for a large
cluster the second dominates.

A worked case: the water hexamer, CCSD(T), to three-body, in aug-cc-pVQZ. The
estimate counts operations in the two steps that set the cost, CCSD
(:math:`o^2 v^4`) and (T) (:math:`o^3 v^4`), against one unfragmented hexamer
calculation, and assumes

* 172 spherical basis functions per water (80 on oxygen, 46 on each hydrogen),
  so 1032 for the hexamer;
* five doubly occupied orbitals per water, so a term with :math:`k` real
  waters has :math:`o = 5k` and :math:`v = n_{\mathrm{basis}} - 5k`; freezing
  the oxygen cores scales every row alike and leaves the ratios where they are;
* ``ssfc`` as 47 rows -- 6 own-basis monomers, then 6 + 15 + 20 = 41 rows,
  one per monomer, pair and trimer, each in all 1032 functions;
* ``vmfc`` as 191 rows -- 6 own-basis monomers, and for each pair and trimer
  its :math:`2^n - 2` subsets and itself, each in the basis of the n-mer it
  belongs to.

.. list-table::
   :header-rows: 1
   :widths: 40 12 24 24

   * - Scheme
     - Rows
     - CCSD, :math:`o^2 v^4`
     - (T), :math:`o^3 v^4`
   * - Unfragmented hexamer
     - 1
     - 1
     - 1
   * - SSFC(3)
     - 47
     - 7.3
     - 3.3
   * - VMFC(3)
     - 191
     - 0.89
     - 0.33

The twenty trimers, each in 1032 functions, are 5.3 of SSFC's 7.3 in CCSD:
a trimer holds half the hexamer's occupied orbitals but slightly more of its
virtuals, so each costs about a quarter of the whole calculation, and there are
twenty of them. VMFC(3) comes in under one calculation because every one of its
191 rows lives in at most 516 functions. These are operation-count ratios, not
timings: they leave out the SCF, the integral transformation, density fitting
and the memory that a thousand-function virtual space needs even when the
occupied space is five orbitals. Read them as the reason to choose ``ssfc`` for
a hexamer on purpose, not by accident, and ``vmfc`` for a much larger cluster at
a level well below its size.

The added rows are ordinary work, so they distribute like everything else --
they are generated before the size sort and handed to the same task server.

Reading the output
------------------

The fragment table names the ghosted rows by signed monomer index. From
``output_<name>_fragments.csv`` for the VMFC dimer above:

.. code-block:: text

   frag_index,level,m1,m2,energy,...
   1,2, 1, 2,-1.5197495012E+02,...    the pair
   2,1, 1,-2,-7.5983256080E+01,...    monomer 1 in the pair's basis
   3,1, 2,-1,-7.5983649444E+01,...    monomer 2 in the pair's basis
   4,1, 1, 0,-7.5981632908E+01,...    monomer 1 in its own basis
   5,1, 2, 0,-7.5983105976E+01,...    monomer 2 in its own basis

A negative index is a ghosted monomer: present as basis functions, absent as
nuclei. Rows 2 and 3 are auxiliary -- they exist to be subtracted by the pair
and are never summed into the total. Subtracting row 4 from row 2 gives monomer
1's own superposition error, 1.02 kcal/mol here.

Under ``ssfc`` the rows are as wide as the system, and a row that holds two or
more monomers is itself a term. For the trimer at ``"level": 2`` (the
columns after ``delta_energy`` are cut):

.. code-block:: text

   frag_index,level,m1,m2,m3,energy,delta_energy
   4,1,1,-2,-3, -7.5984866063516492E+01, -7.5984866063516492E+01   monomer 1, full basis
   5,1,2,-1,-3, -7.5985607768766613E+01, -7.5985607768766613E+01   monomer 2, full basis
   6,1,3,-1,-2, -7.5984738393132929E+01, -7.5984738393132929E+01   monomer 3, full basis
   7,1,1,0,0, -7.5983974465724899E+01, -7.5983974465724899E+01   monomer 1, own basis
   8,1,2,0,0, -7.5983974465725055E+01, -7.5983974465725055E+01   monomer 2, own basis
   9,1,3,0,0, -7.5983974465724941E+01, -7.5983974465724941E+01   monomer 3, own basis
   3,2,2,3,-1, -1.5197705606853540E+02, -6.7099066358622395E-03   pair 2,3 with 1 ghosted
   1,2,1,2,-3, -1.5197690054615407E+02, -6.4267138709652727E-03   pair 1,2 with 3 ghosted
   2,2,1,3,-2, -1.5197114300313342E+02, -1.5385464839994256E-03   pair 1,3 with 2 ghosted

The own-basis monomers (rows 7 to 9) are the one-body term. The full-basis
monomers (rows 4 to 6) are auxiliary: each is subtracted by the pairs holding
it, and never summed. The three pairs are the two-body term, each ghosting the
monomer it does not hold. At ``"level": 3`` the whole trimer is one more row,
with nothing left to ghost, and it is summed too.

The JSON output names the scheme beside the fingerprint, ``"counterpoise":
"ssfc"``, and is otherwise as for any MBE; a ghosted row in the per-fragment
breakdown (``"system": {"fragment_breakdown": "json"}``) lists its ghosts under ``ghosts``
and its real monomers under ``indices``. See :doc:`json_output`.

Energies only
-------------

Gradients, Hessians and dipole-derivative runs are **refused** under either
scheme rather than approximated, and refused before any fragment has run. The
routine that collapses the many-body recursion into one weight per fragment
predates counterpoise and assembles those weights from unghosted subsets; run
under counterpoise it would return a derivative that is wrong without looking
wrong. Run the energy, or drop the correction. A geometry optimization is a
derivative run and is refused with them. The dipole moment, which rides the same
recursion, is corrected the same way.

What it does not combine with
-----------------------------

Counterpoise is carried by the ghosted rows of the MBE term list that the driver
generates. Four things do not read that list, and a deck combining them with
``"vmfc"`` or ``"ssfc"`` is refused before any work starts:

* **GMBE** (``"method": "gmbe"``) builds its terms by
  inclusion--exclusion over overlapping primaries instead.
* **FMO**, **EE-MBE** and **EFMO** (``"method"``) build their own term lists.
* **A supplied term list** -- the Python and C interface, which hands the
  program a list of its own -- is used as given, and the ghosted rows are built
  only when the driver generates the list. Counterpoise there used to be
  dropped silently, returning an uncorrected energy; it is refused.
* **GFN1**, **GFN2** and **EFP** have no ghost centre to construct -- see below.

The first three would have returned a valid *uncorrected* energy, which is the
number the deck asked this program not to produce, and nothing about it would
have said so. That is why they are errors and not warnings.


Why not xTB
^^^^^^^^^^^

Not because superposition error is absent. GFN1 and GFN2 carry a minimal
valence STO basis, so a monomer inside a dimer does borrow a little from its
neighbour, and that borrowing is real. The correction is refused because it is
**unconstructible** there, which is a stronger reason than being small.

A ghost centre means basis functions with no nuclear charge, and that separation
only exists if the two are separable in the first place. In an ab initio method
they are: libcint is handed a shell list and a charge list, and one can be
zeroed without touching the other. In GFN the basis is welded to the atom. The
shell exponents, the diagonal Hamiltonian elements, the electronegativity
equilibration that sets the charges, the repulsion, the dispersion coefficients
-- all of them are indexed by element. Remove the nucleus and you have removed
the parameters that define the functions, and every term built on them. There is
no partial atom left to place.

Which is why the tblite interface only ever receives ``element_numbers``: there
is nothing else it could be given. A ghosted fragment handed to it becomes a
real atom weighed against an electron count computed as though that atom were
absent -- not an approximation but a different molecule, with the wrong charge.
A water dimer at GFN2 came back ``+0.0073`` Hartree that way before this was
refused.

The same holds for EFP, which describes a fragment by a classical potential
rather than by a basis at all.

Where it is checked
-------------------

``test/test_mqc_counterpoise.f90`` holds the tests on the properties the schemes
rest on: that a ghost carries basis functions and no nucleus, that ghosting
leaves the AO space unchanged, that a monomer in the pair's basis is genuinely
lower than in its own, and that a two-fragment VMFC expansion reproduces the
counterpoise-corrected supermolecular interaction energy -- checked against
``0.009315543671`` Hartree, a number SAPT reaches through the dimer-centred
basis and none of this machinery.

That last one matters because a counterpoise correction can be wired up, run,
and produce exactly the uncorrected number. The SSFC tests in the same file are
built against that failure. A water trimer through the driver at level 3 is
compared, to 1e-10 Hartree, with :math:`E_{\mathrm{cluster}} + \sum_i [E_i(i) -
E_i(\mathrm{FB})]` assembled by hand from single calculations with explicit
ghost masks, so the expansion and the identity share no code. A dimer gives the
same energy under both schemes, and a trimer at level 2 gives a different one,
which is what stops SSFC from silently running as VMFC. An interaction energy
under SSFC is checked against the same reference's terms in an ordinary run.

``test/test_mqc_counterpoise_schemes.f90`` holds the bookkeeping: which rows each
scheme sums, that each scheme accepts its own term lists and refuses the other's,
and SSFC totals against inclusion--exclusion done directly over full-basis
energies, at the top level, truncated, with no pairs, and under a reference
fragment with screening.

``test/test_mqc_term_list.f90`` covers the term list itself. At level 3 the VMFC
rule first has depth: a trimer owes six ghosted subsets, three of them pairs, and
each must be ghosted against the trimer rather than against itself. The SSFC
rows must ghost exactly the system's complement, leave no all-real n-mer below
level :math:`N`, and keep ghosting the whole system when a reference fragment
reduces the list. Both follow screening, so a pair dropped by a cutoff does not
leave its ghosted monomers behind to be paid for and subtracted by nothing.

``test/test_mqc_wide_rows.f90`` covers rows wider than the expansion level, which
SSFC produces at every level, through the recursion, the checkpoint and the
written tables, at 12 and at 100 monomers.

``test/test_mqc_config_roundtrip.f90`` covers the refusals above.

The cases against PySCF are in ``validation/inputs/cpu/mqc/counterpoise/``:
the water trimer at Hartree-Fock/6-31G for SSFC(2) and SSFC(3), at
MP2/cc-pVDZ with the frozen core for SSFC(2), and the water prism at
Hartree-Fock/STO-3G for SSFC(2) and SSFC(3). The reference is PySCF run with
ghost atoms and this repository's own basis files. All five agree to better than
1e-10 Hartree. The MP2 case is there for the frozen core: counting the ghost
oxygens as frozen-core atoms, one more occupied orbital frozen per ghost, moves
that total by 0.66 Hartree. The
interaction-energy counterparts, ``reference_fragment: 5`` on the prism, are in
``validation/inputs/cpu/mqc/interaction_energy/`` as ``ie_prism_sto-3g_l2_ssfc.json``
and ``ie_prism_sto-3g_l3_ssfc.json``, beside the VMFC ones. Run them with
``python3 run_validation.py --manifest validation_tests_cpu.json -t ssfc`` from
inside ``validation/``.
