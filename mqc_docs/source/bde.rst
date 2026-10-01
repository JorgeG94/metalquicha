==========================
Bond dissociation energies
==========================

``mqc.bde`` computes the energy to break a bond in a small molecule, for one bond
or for every bond to hydrogen at once. Like :doc:`pka` it is a Python workflow
on top of the interface described in :doc:`python_interface`: every calculation
it makes is an ordinary ``mqc.MBE`` run, and nothing in the Fortran changed for
it. It shares its thermochemistry and its conformer, optimization, Hessian and
single-point stages with ``mqc.pka``.

It is built for molecules of thirty to forty atoms, at GFN2-xTB **in the gas
phase**, with the single point as the stage expected to move to DFT; see
`Swapping in a DFT single point`_. It is a screening tool: read
`What to distrust`_ before a number is used for anything.

.. code-block:: python

   import mqc
   from mqc import bde

   with mqc.session():
       result = bde.scan("methanol.xyz", bonds="X-H")
       print(result.table())              # weakest bond first

       row = result.bonds[0]
       print(row.bde_kcal, row.bde_kj, row.d0_kcal, row.dg_kcal)
       result.to_json("methanol_bde.json")

``python/examples/bde.py`` runs this end to end on methanol and phenol.

Definitions
===========

For a parent :math:`AB` that splits into :math:`A + B`, with :math:`E` the
electronic energy at the single-point level:

.. math::

   D_e &= E(A) + E(B) - E(AB) \\
   D_0 &= D_e + \mathrm{ZPE}(A) + \mathrm{ZPE}(B) - \mathrm{ZPE}(AB) \\
   \mathrm{BDE}(298) &= H(A) + H(B) - H(AB) \\
   \Delta G &= G(A) + G(B) - G(AB)

:math:`D_e` is the bottom-of-the-well energy, with no nuclear motion. :math:`D_0`
adds the change in zero-point energy and is the 0 K enthalpy of dissociation.
The conventional "bond dissociation energy" in a table is the reaction enthalpy
at 298 K, :math:`\Delta H`, which is the third line:
:math:`H = E + \mathrm{ZPE} + E_\mathrm{vib} + E_\mathrm{rot} + E_\mathrm{trans} + RT`,
where the :math:`RT` is the :math:`pV` of an ideal gas and is what makes
:math:`\Delta H` exceed :math:`\Delta E` by about :math:`RT` more than the thermal
terms alone. The temperature
is the one in the protocol (298.15 K unless the ``frequencies`` stage's
``hessian`` block says otherwise) and is in the result as ``temperature_K``;
``BondRow.bde_kcal`` is the enthalpy at that temperature whatever its name.

:math:`\Delta G` is reported beside them because the shared thermochemistry gives
it for free. It is **not** the quantity usually called a BDE: dissociation gains
translational and rotational entropy, so it is smaller than :math:`\Delta H` by
about 8 to 10 kcal/mol at 298 K for one bond cut into two pieces. By default it is
at a 1 atm standard state for every species, which is what a gas-phase free energy
of dissociation is.

Everything is in kcal/mol, and ``BondRow`` carries ``d_e_kj``, ``d0_kj``,
``bde_kj`` and ``dg_kj`` as well.

What it does, and what it does not
==================================

**Fragments are generated, not supplied.** A bond is ``(i, j)``, 0-based like every
atom index in the API. The parent's connectivity is perceived exactly as
``System.perceive_bonds`` does it: two atoms are bonded when
:math:`d < 1.2\,(r_i + r_j)` with Cordero covalent radii (copied from
``src/core/mqc_atomic_radii.f90``, and so changing together with it by hand).
Removing the bond ``i-j`` must leave two connected pieces; the one holding ``i``
is ``side_i`` and the one holding ``j`` is ``side_j``. ``bde.split`` does only
this and returns the plan without computing anything.

**A bond in a ring is refused.** Cutting it leaves one molecule, a ring-opened
diradical, not two fragments. Naming one in a list is an error. Selecting by
keyword (``bonds="all"``) skips the ring bonds and lists them in ``result.skipped``.
A ring-opening energy is a different quantity; compute the open-chain diradical
yourself and take the difference. A parent that is not one connected molecule at
the tolerance (a complex, a salt, a hydrogen-bonded dimer) is refused as well.

**It does not choose the bonds for you.** ``bonds`` is ``"X-H"`` (every bond to a
hydrogen), ``"all"`` or a list of ``(i, j)``. Two atoms that are not bonded at the
tolerance are an error, with the distance in the message; pass ``tolerance=`` for
a stretched bond.

Everything that can be wrong about the request -- a missing bond, a ring, a
charge with nowhere to go, a spin that does not match the electron count -- is
raised before the first calculation starts.

Charges and spins
=================

**Homolytic, the default.** Each electron of the bond goes to one fragment, so a
closed-shell neutral parent gives two doublet radicals. The multiplicities are
worked out from the electron count of each fragment and checked against it.

**A charged parent is never guessed.** ``scan(..., charge=+1,
charge_on="j")`` says that the fragment holding the second atom of the bond keeps
the charge; ``"i"`` is the first atom, as written in the list (the lower index for
``"X-H"`` and ``"all"``). The question has a chemical answer the program cannot
see: in methylammonium :math:`\mathrm{CH_3NH_3^+}` cut at C-N, the
:math:`\mathrm{NH_3^{\bullet+}}` side keeps it for a homolysis, and putting it on
the methyl would leave two closed shells, which is a heterolytic split under
another name -- and is refused as such. A neutral parent needs nothing.

**Heterolytic**, ``heterolytic=True`` with ``cation="i"`` or ``"j"``: the bond
pair goes to one side, so one fragment is the cation (:math:`+1` relative to the
homolytic split) and the other the anion (:math:`-1`), both normally closed
shells. A gas-phase heterolytic dissociation is hundreds of kcal/mol, because
nothing screens the charge separation, and is not comparable with a homolytic
one; use a solvent in the protocol and read `Solvent`_.

**A parent with spin** (a radical, a triplet) needs ``multiplicities=(m_i, m_j)``:
how the spin divides between the pieces is a choice. Whatever is given is checked
against the electron count of each side.

The protocol
============

``bde.Protocol`` is ``pka.Protocol`` with different defaults; the stages, the
conformer window, quasi-RRHO and the way two of them run the ``mqc`` executable are
as described in :doc:`pka`. What differs:

* **The default is GFN2-xTB with no solvent at every stage**, where ``pka``'s is
  ALPB water. The two are not interchangeable and the difference is deliberate.
* ``standard_state`` is False: the free energy is at 1 atm for every species, not
  corrected to 1 M.
* **The conformer search is on for the parent and off for the fragments**
  (``parent_conformers=True``, ``fragment_conformers=False``). A fragment is a
  piece of the parent, started from its relaxed coordinates, small, and usually a
  radical whose conformers CREST has no more reason to rank than the parent's
  geometry gave. Turn it on for a large radical whose shape changes a lot on
  losing the bond. ``conformers=None`` switches the search off everywhere.

The parent is evaluated first. Its lowest free energy conformer, optimized, is the
structure the fragments are cut from, along the *input* connectivity; the result
warns if the relaxed parent has a different bond graph.

**Vertical and relaxed fragment energies.** When the optimization stage is on, each
fragment's energy at the geometry it had inside the relaxed parent is computed too,
one extra single point. ``SpeciesData.e_vertical_hartree`` is that and
``relaxation_kcal`` is how far it drops on relaxing: a large number says the dissociation is
far from vertical, a negative one says the optimization or the conformer choice went
somewhere worse and is flagged in the warnings. An atom has nothing to relax.

Radicals in xTB
===============

Fragments are open shells, and that is the least reliable thing the workflow asks
of GFN2. Read from ``mqc_method_xtb.f90`` rather than assumed:

* The multiplicity reaches tblite as the number of unpaired electrons,
  ``uhf = multiplicity - 1``, when the molecule is built.
* The wavefunction is created with a single spin channel
  (``new_wavefunction(..., nspin=1, ...)``), and mqc does not add the
  spin-polarization interaction to the calculator. So a doublet is an
  **unpolarized** calculation with the unpaired electron placed in a singly
  occupied orbital: it runs, converges as a closed shell does, and has none of the
  spin-polarization energy that a spin-resolved xTB calculation would add. It is
  not an unrestricted calculation and not spin contaminated in the UHF sense.
* No :math:`\langle S^2\rangle` is reported for xTB, in the output or the JSON.
  (The ab initio unrestricted code computes one and does not write it either.)
  The workflow therefore cannot tell you a radical went wrong in that way.

The consequence is the expected one for a semiempirical method that was
parametrized largely on closed shells and has no spin-polarization term here:
bond dissociation energies that are a screening quantity, trustworthy for
*trends* between similar bonds (is this C-H weaker than that one) and good to
several kcal/mol, sometimes more, in absolute value. The result puts a note in
``warnings`` whenever an open-shell species was run through xTB.

Atoms
=====

A single-atom fragment (the H of an X-H bond, a halogen) has no vibrations and
no rotations and **is never given a Hessian**. It is one single point, and the
thermochemistry is analytic:

.. math::

   H = E + \tfrac52 RT, \qquad
   S = S_\mathrm{trans} + R\ln(2S+1)

with the Sackur-Tetrode translational entropy

.. math::

   S_\mathrm{trans} = R\left[\ln\!\left(\frac{(2\pi m k_B T/h^2)^{3/2}\,k_B T}{P}\right)
   + \tfrac52\right],

26.015 cal/(mol K) for hydrogen at 298.15 K and 1 atm, and
:math:`R\ln 2 = 1.377` cal/(mol K) of spin degeneracy, 27.39 in all. :math:`5/2\,RT`
is 1.481 kcal/mol at 298.15 K. This is why the Hessian is skipped and not merely believed: the Fortran
thermochemistry decides a molecule is linear when its smallest moment of inertia
vanishes, which all three do for one atom, and would then add a linear rotor's
:math:`RT` of rotational energy that an atom does not have.
**Spin-orbit coupling is ignored**: a halogen atom's ground state is
:math:`^2P_{3/2}`, a four-fold level with a second one a few kcal/mol above it
(about 2.5 kcal/mol for Cl, 10 for Br, 22 for I), and the 2S+1 degeneracy of 2 is the
same approximation applied to those. The result says so in ``warnings``. The same
path serves ``mqc.pka`` for a one-atom microstate such as a halide.

The electronic entropy
======================

The Fortran thermochemistry has a spin term, :math:`S_\mathrm{elec} = R\ln(2S+1)`,
and **no caller gives it the multiplicity**: every call to
``compute_thermochemistry`` leaves it at its default of 1. So for a radical the
JSON block reports ``spin_multiplicity: 1`` and an electronic entropy of zero.
The workflow computes the term itself, from the species' own multiplicity, and does
not use the block's value, for every species -- molecules, radicals and atoms
alike -- so :math:`\Delta G` and each free energy carry :math:`R\ln 2` per doublet.
:math:`D_e`, :math:`D_0` and :math:`\mathrm{BDE}(298)` are enthalpies and do not
see it. For the two doublets of a homolysis it is :math:`-2RT\ln 2 = -0.82`
kcal/mol in :math:`\Delta G` at 298 K. (For a closed shell the term is zero and
nothing changed in ``mqc.pka``.) Fixing it in ``mqc_thermochemistry.f90`` and its
callers is a one-line change per call, and when it lands the Python value agrees
with it.

Scanning many bonds
===================

``bde.scan(parent, bonds, protocol, charge=0, multiplicity=1)`` evaluates the
parent once and each *distinct* fragment once, then assembles one row per bond. The
three equivalent C-H bonds of a methyl group give one pair of fragments, and every
X-H bond shares the hydrogen atom: methanol's four bonds to hydrogen are
three fragment calculations (:math:`\mathrm{CH_2OH}`, :math:`\mathrm{CH_3O}` and H)
and the parent. The rows are sorted by :math:`\mathrm{BDE}(298)`, so the weakest
bond, the likely hydrogen-abstraction site, is first.

Two fragments are the same species when their *key* agrees: the sorted element
counts, the charge, the multiplicity and a hash of the fragment's bond graph. The
hash is a colour-refinement invariant: each atom starts as its element and is
relabelled, round after round, by its own label and the sorted labels of its
neighbours, and the sorted final labels are hashed. It does not depend on the order
the atoms are listed in and separates constitutional isomers (butyl from isobutyl
radicals, :math:`\mathrm{CH_2OH}` from :math:`\mathrm{CH_3O}`). **Its limits:** it is
not a canonical form, so some regular graphs that colour refinement cannot tell
apart share a hash; it carries no geometry, no bond orders and no stereochemistry,
so cis and trans, and enantiomeric, fragments are one species; and it is a
property of the connectivity perceived at the tolerance, which is only as good as
the radii. The cost of a false merge is that one fragment's energy stands in for
another's, so for stereochemically interesting radicals compute the bonds
separately.

Solvent
=======

The default is the gas phase because that is what a tabulated BDE is. A solvent
changes the meaning. Putting ``"xtb": {"solvent": "toluene", "solvation_model":
"alpb"}`` into the stages gives a *solution-phase* dissociation enthalpy in which the
solvation free energy of every species is inside its energy, with no solvation
enthalpy separated out. For a homolysis of a neutral molecule into radicals the
effect is small and not obviously of the right sign; for a heterolysis it is the
whole point, and the reason to use a solvent at all. The result notes it whenever a
solvent is set. Set ``standard_state=True`` as well if the free energy should be
for 1 M.

Swapping in a DFT single point
==============================

The single point is the stage whose result needs to be better than xTB, and the
geometries, zero-point energies and thermal corrections stay at GFN2 -- a frequency
error that is small against the energy error it fixes. As in :doc:`pka`, it is one
entry of the protocol:

.. code-block:: python

   protocol = bde.Protocol(
       single_point={
           "method": "dft", "functional": "pbe0", "basis": "def2-svp",
           "scf": {"unrestricted": True},
       },
   )

The energy is then :math:`E_\mathrm{DFT} + (H\ \text{or}\ G)_\mathrm{corr}^\mathrm{GFN2}`
and the workflow runs the single point itself, at each conformer's geometry and at
each fragment's vertical geometry. Two things come with it:

* **Open shells need an unrestricted reference**, which is what the ``scf``
  block asks for, and an unrestricted determinant is spin contaminated: a doublet
  with :math:`\langle S^2\rangle` well above 0.75 is mixed with a quartet, its
  energy is too low, and the BDE is too small by an amount that is not small for a
  delocalized radical. The workflow does not read :math:`\langle S^2\rangle` back
  (the unrestricted code computes it and does not write it to the JSON), so check
  the output log of the radicals, and trust a BDE less where it is large.
* A hybrid functional lowers the contamination and a pure one raises it; neither
  removes it. A restricted open-shell or a spin-projected treatment is not
  available here.

What to distrust
================

**The number is as good as the radical.** See `Radicals in xTB`_. The ordering of
similar bonds is more reliable than the value, and the weakest bonds in a molecule
are weakest for chemical reasons a semiempirical method captures (conjugation, a
captodative centre), so the *first row* of a scan is the part to believe.

**The geometry may be the wrong conformer.** The parent is searched; the fragments
by default are not, and a fragment is optimized to the nearest minimum of the
parent's own shape.

**Imaginary frequencies** are counted per species and surfaced in ``warnings``. A
soft imaginary mode on a methyl rotor is routine; a large one means the structure is
not a minimum and its enthalpy is not a minimum's.

**Quasi-RRHO moves the entropy and not the enthalpy**, so :math:`\mathrm{BDE}(298)`
does not see it and :math:`\Delta G` does. The harmonic enthalpy of a soft mode is
the usual error of this level, small against the rest.

**The symmetry number** is whatever the thermochemistry block used (1); the
workflow does not detect it. For the free energy it moves a fragment by
:math:`RT\ln\sigma`, which is the one place a methyl group's three-fold symmetry
would show.

**Reading the result.** ``BDEResult`` carries ``bonds`` (``BondRow``, sorted),
``parent`` and ``species`` (``SpeciesData``: energies, ZPE, :math:`H`, :math:`G`,
imaginary-mode and conformer counts, vertical and relaxed energies, and the
per-conformer record in ``detail``), ``skipped``, ``warnings``, ``table()``,
``to_dict()`` and ``to_json()``. The ``output_*.json`` files of every run are left
in the working directory as the record of what was computed, labelled
``<prefix>_<species>_c<k>_freq``, ``_sp`` and ``_vert``.

Run it inside ``mqc.session()`` and only once, as for ``mqc.pka``.
