Counterpoise correction (VMFC)
==============================

A fragment solved in its own basis is described worse than the same fragment
solved inside a larger one, because in the larger one it can borrow its
neighbour's basis functions. That difference is **basis-set superposition
error**, and in a many-body expansion it does not cancel: the pair term
subtracts monomers described in a small basis from a pair described in a large
one, so the pair looks more bound than it is. The error survives truncation, and
each higher order adds its own rather than cancelling the last.

Counterpoise removes it by solving every subfragment in its **parent's** basis,
so the borrowing appears on both sides of each difference. This program
implements the Valiron--Mayer function counterpoise scheme, ``vmfc``.

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

``"counterpoise"`` takes ``"vmfc"`` or ``"none"``, and defaults to ``"none"`` --
which is the uncorrected expansion every deck written before this existed
already asked for. Nothing else in the deck changes.

What it is worth
----------------

The water dimer ``validation/inputs/sample_inputs/w2_dimer.xyz`` at
Hartree-Fock/6-31G, a basis small enough that the borrowing is large:

============================  =================  ==================
                              Hartree            kcal/mol
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
3-body         -0.0007565683      +0.0026639427    2.15 kcal/mol
Total        -227.9704976390    -227.9637736609    4.22 kcal/mol
===========  =================  =================  ==================

Two things in that table are the argument for the correction. The uncorrected
total, ``-227.9704976390``, is exactly the supermolecular energy -- an MBE at
level ``N`` over ``N`` fragments is exact, and this one is exact to every digit
printed, which is how you know the expansion itself is right and the difference
is entirely superposition error. And the **3-body term changes sign**. The
uncorrected three-body contribution is small -- half a kcal/mol -- while the
superposition error at that order is four times larger, so for this system a
many-body analysis reading the uncorrected numbers gets the wrong sign for the
non-additive part of the binding. The correction is not only a shift in the
total; it can change what the decomposition says.

The one-body term is untouched on purpose. Each monomer keeps its own basis
there, because that is the reference the interaction energy is measured *from*.
Only the corrections are ghosted.

What it costs
-------------

Every n-mer contributes its :math:`2^n - 2` proper subsets as extra terms, so
the row count goes from :math:`N + \sum_n \binom{N}{n}` to
:math:`N + \sum_n \binom{N}{n}(2^n - 1)`:

=========  ==================  ==================  ========
Level      20 fragments        50 fragments        Factor
=========  ==================  ==================  ========
2          210 -> 590          1275 -> 3725        ~2.9
3          1350 -> 8570        20875 -> 140925     ~6.5
4          6195 -> 81245       251175 -> 3595425   ~14
=========  ==================  ==================  ========

Wall-clock is worse than the row count suggests, because every added row is
solved in a *larger* basis than the term it corrects: a ghosted monomer of a
trimer carries the trimer's full basis. Counterpoise is affordable at level 2
and expensive above it.

The added rows are ordinary work, so they distribute like everything else --
they are generated before the size sort and handed to the same task server.

Reading the output
------------------

The fragment table names the ghosted rows by signed monomer index. From
``output_<name>_fragments.csv`` for the dimer above:

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

Energies only
-------------

Gradients, Hessians and dipole derivatives are **refused** under counterpoise
rather than approximated. The routine that collapses the many-body recursion
into one weight per fragment predates counterpoise and assembles those weights
from unghosted subsets; run under VMFC it would return a derivative that is
wrong without looking wrong. Run the energy, or drop the correction.

What it does not combine with
-----------------------------

Counterpoise is carried by the ghosted rows of the MBE term list. Three things
do not read that list, and a deck combining them with ``"vmfc"`` is refused
before any work starts:

* **GMBE** (``"allow_overlapping_fragments": true``) builds its terms by
  inclusion--exclusion over overlapping primaries instead.
* **FMO** and **EE-MBE** (``"expansion"``) build their own term lists.
* **GFN1**, **GFN2** and **EFP** have no ghost centre to construct -- see below.

The first two would have returned a valid *uncorrected* energy, which is the
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

Other schemes
-------------

``vmfc`` is the only one implemented. It is the scheme that is consistent order
by order -- each n-body term is corrected in its own n-mer basis -- which is why
it is the one worth having first, and also why it is the one that multiplies the
term count. The site-function alternatives correct every fragment in the *full
system* basis instead, which is fewer calculations but each in a basis the whole
point of fragmenting was to avoid.

Where it is checked
-------------------

``test/test_mqc_counterpoise.f90`` holds nine tests on the properties the scheme
rests on: that a ghost carries basis functions and no nucleus, that ghosting
leaves the AO space unchanged, that a monomer in the pair's basis is genuinely
lower than in its own, and that a two-fragment VMFC expansion reproduces the
counterpoise-corrected supermolecular interaction energy -- checked against
``0.009315543671`` Hartree, a number SAPT reaches through the dimer-centred
basis and none of this machinery.

That last one matters because a counterpoise correction can be wired up, run,
and produce exactly the uncorrected number.

``test/test_mqc_term_list.f90`` covers the term list itself at level 3, where
the rule first has depth: a trimer owes six ghosted subsets, three of them
pairs, and each must be ghosted against the trimer rather than against itself.
It also checks that the rows follow screening, so a pair dropped by a cutoff
does not leave its ghosted monomers behind to be paid for and subtracted by
nothing.

``test/test_mqc_config_roundtrip.f90`` covers the refusals above.
