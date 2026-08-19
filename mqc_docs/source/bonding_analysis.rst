Quasi-atomic bonding analysis
=============================

A converged wave function is a set of delocalised molecular orbitals and a
number. It does not, on its own, say where the bonds are. This analysis
recovers that picture: which atoms are bonded, by what kind of orbital, how
strongly, and where the lone pairs are -- all from the wave function itself,
with no reference to a Lewis structure anyone drew.

.. code-block:: json

   "driver": "Energy",
   "properties": {
     "bonding_analysis": {
       "type": "gms_quao",
       "energy_threshold": 2.0
     }
   }

- ``type``: ``"gms_quao"`` (or ``"quao"``) for the Ruedenberg quasi-atomic
  analysis, ``"none"`` to switch it off. Required inside the block.
- ``energy_threshold``: Orbital pairs whose kinetic bond order is weaker than
  this, in kcal/mol, are counted rather than printed (default: 1.0).

``properties`` sits beside ``keywords``, not inside it, and the distinction is
worth keeping. ``keywords`` say how to compute the wave function and change the
number that comes out; ``properties`` ask for something further to be done with
one already determined, and change no energy. That is why the driver stays
``"energy"`` -- this is a report on the calculation, not a different one.

What it produces
----------------

Formyl chloride in 6-31G::

    atoms
       atom      population      charge
      C 1          5.720976    0.279024
      O 2          8.343536   -0.343536
      Cl 3        17.100513   -0.100513
      H 4          0.834975    0.165025
      total       32.000000

    bonds
       bond              type        order    kcal/mol   orbitals
      C 1 - O 2         sigma          0.9658     -67.81     11 / 15
      C 1 - Cl 3        sigma          0.9216     -43.29      9 / 19
      C 1 - H 4         sigma          0.9495     -35.77     10 / 20
      C 1 - O 2         pi             0.9281     -28.12      8 / 14

    delocalization
       donor             into                order    kcal/mol   orbitals
      O 2 p-lone        C 1-Cl 3 sigma        0.3780      -9.21     13 / 9
      Cl 3 p-lone       C 1-O 2 pi            0.3065      -5.38     17 / 8

The C=O double bond appears as two rows, a sigma and a pi. The lone pairs are
found rather than assumed, and the delocalization block is where a lone pair
donates into a bond -- chemically a different phenomenon from a covalent bond,
so it is a different table.

Three tables follow: **bonds**, where both orbitals are bonding and point at
each other; **delocalization**, everything else above the threshold; and
**orbitals**, one row per quasi-atomic orbital with its occupation and its
s/p/d composition.

Two numbers per interaction, and they are not interchangeable. The **bond
order** is a population and says how much density two orbitals share. The
**kinetic bond order** is an energy and says what that sharing is worth; it is
negative for a bonding interaction, and it is what the tables are sorted on.
They disagree about ordering more often than one would expect -- above, the
C--Cl bond order beats the C--H one while its kinetic bond order beats it by
more than twice as much, because the C--H kinetic integral is smaller.

The kcal/mol column is Paper II's kinetic bond order: the raw kinetic
interference energy scaled by an empirical tenth, which brings it onto the
scale of tabulated bond energies. That factor is admitted to be empirical in
the paper and should be read as such.

Where the energy goes
---------------------

The tables above say which atoms are bonded and rank the interactions. They do
not say what the molecule's energy is made of, and the two questions have
different answers: the kinetic bond order carries an empirical factor of a tenth
to reach a familiar scale, so it is a gauge rather than an energy.

Four further tables resolve the actual energy. Every term in

.. math::

   E = \sum_{pq} \gamma_{pq} T_{pq}
     + \sum_{pqA} \gamma_{pq} V^{A}_{pq}
     + \tfrac{1}{2}\sum_{pqrs} \Gamma_{pqrs} (pq|rs)
     + \sum_{A<B} \frac{Z_A Z_B}{R_{AB}}

has orbitals belonging to particular atoms, so each term belongs either to one
atom or to one pair of atoms. Nothing is modelled and nothing is discarded; the
tables are a regrouping of a sum, and they add back up.

The last of the four is the one to read::

    energy decomposition
       intra-atomic                    mhartree    kcal/mol
      O 1                             -73655.015  -46219.220
      H 2                                 -6.997      -4.391
      H 3                                 -6.997      -4.391

       interatomic                     mhartree    kcal/mol
      O 1      -- H 2              -1186.141    -744.314
      O 1      -- H 3              -1186.141    -744.314
      H 2      -- H 3                 56.511      35.461

            intra-atomic total      -73.669009 hartree
             interatomic total       -2.315771 hartree
                         total      -75.984780 hartree

**These are not bond energies.** ``intra-atomic`` is an atom as it exists in the
molecule -- deformed, and in water's case with the hydrogens stripped to about
0.64 electrons each -- not a free atom. The difference between the two is the
adaptation energy, which is a separate calculation and is not computed here. So
an O--H figure of -744 kcal/mol is the whole interaction between an oxygen and a
hydrogen sitting where they sit, nuclear repulsion included, and not the energy
of breaking the bond.

The signs are worth reading. H--H comes out positive: the two hydrogens carry
like charges and repel, which the bond-order table cannot show because there is
no bond there to rank.

The three tables the total is built from separate the physics. **Kinetic** is
where covalent binding comes from in this analysis -- the interatomic entries are
the interference energy of eq (1) unscaled, so the same quantity the kinetic bond
order is a tenth of. **Nuclear attraction** carries three atomic labels rather
than two, since the nucleus doing the attracting need not sit on either orbital;
an atom's own density in a foreign nuclear field is charged to the pair, not to
the atom. **Two-electron** carries four.

Checking it
~~~~~~~~~~~

Three lines follow the tables::

                          one-electron     -122.949561 hartree
                          two-electron       37.795213 hartree
                     nuclear repulsion        9.169569 hartree
                                 total      -75.984780 hartree
            from the orbitals directly      -75.984780 hartree

The last line is computed from the converged orbitals in the atomic-orbital
basis, with nothing quasi-atomic in it. For a single determinant the two agree
to rounding, and that agreement is the check worth having: every internal
consistency test the decomposition applies to itself would still pass if the
quasi-atomic basis failed to hold the whole density.

When it does not, a further line reports the difference, under ``outside the
quasi-atomic span``. That is the same shortfall the population sum reports for a
correlated density, and it arises for the same reason.


Correlated wave functions
~~~~~~~~~~~~~~~~~~~~~~~~~

The energy tables take an MCSCF wave function as well as a determinant, and
nothing extra is asked for. What changes is one term. A determinant's
two-particle density is fixed by its one-particle one,

.. math::

   \Gamma_{pqrs} = \gamma_{pq}\gamma_{rs} - \tfrac{1}{2}\gamma_{ps}\gamma_{rq}

and correlation is precisely what makes that false. The difference between the
true two-particle density and this expression is the cumulant, it is zero unless
all four of its indices are active, and it is added to the two-electron term.

That locality is what keeps the cost down: the determinant expression is
evaluated over the whole quasi-atomic basis whatever the wave function, and the
correction lives in the active space alone.

Two lines appear that do not for a determinant::

    active orbitals outside the quasi-atomic span   5.67E-02
    ...
       reported by the calculation     -109.090026 hartree
     outside the quasi-atomic span       -0.014857 hartree

The first says the active orbitals are not entirely inside the space the
quasi-atomic orbitals span, and they are not: the valence-virtual orbitals are
chosen to look like free-atom orbitals, while the active orbitals of a converged
MCSCF are chosen to lower an energy. Those are different choices and they do not
give the same space.

The second is what that costs. For N\ :sub:`2` in cc-pVDZ with a CAS(6,6) the
decomposition reaches -109.075169 hartree against a CASSCF that reported
-109.090026, so 14.9 of the 135.9 millihartree of correlation energy are outside
what this analysis can see -- about 89% described. The comparison is against the
energy the calculation reported, since the reference built from the orbitals
assumes a determinant and would be short by the whole correlation energy.



What is classical and what is not
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each pair interaction is split into the part an electrostatic model could
produce and the part it could not::

    interatomic interactions
       pair                  total     classical  interference
      O 1      -- H 2          -1186.141      -329.244      -856.897
      O 1      -- H 3          -1186.141      -329.244      -856.897
      H 2      -- H 3             56.511        57.617        -1.106
       millihartree

**Classical** is one atom's density sitting in the other's nuclear field, the
repulsion between the two atomic charge clouds, and the repulsion of the two
nuclei. Each is large and they very nearly cancel, because a neutral atom is a
nearly neutral thing to be near.

**Interference** is what is left: density shared between the two atoms, which no
electrostatic model has any account of. The kinetic contribution is entirely
interference by construction -- two orbitals on different atoms is what the word
means.

The claim the papers make is that covalent binding comes from interference
rather than from electrostatics, and water says so plainly. Interference carries
857 of the 1186 millihartree of the O--H interaction. And the H--H pair, where
there is no bond at all, comes out at +57.6 millihartree classical against -1.1
of interference -- two like charges repelling, with essentially nothing shared.
That the analysis finds no interference where chemistry says there is no bond is
worth more than the O--H number, since nothing in the construction was told
where the bonds are.

Energy of formation
~~~~~~~~~~~~~~~~~~~

The intra-atomic terms above are large because they contain most of each atom's
own energy. Subtracting the atom as it would be on its own turns them into
something readable::

    energy of formation
       atom          in molecule      free atom     adaptation
      O 1            -73.655015     -74.780310       1.125295
      H 2             -0.006997      -0.498233       0.491236
      H 3             -0.006997      -0.498233       0.491236

              adaptation total        2.107766 hartree
             interatomic total       -2.315771 hartree
           energy of formation       -0.208004 hartree    -130.525 kcal/mol

**Adaptation is positive and that is the point.** An atom in a molecule is
promoted into the hybridisation the bonding needs, and in a polar bond it is
stripped of charge as well; both cost energy against the free atom. Water's
oxygen pays 1.13 hartree and each hydrogen 0.49. Binding happens because the
interatomic terms are more negative than the adaptation is positive, and the
difference is what the molecule is worth.

Each free atom is solved unrestricted at its Hund's-rule ground state in
**exactly the basis functions it contributes to the molecule**, which is the same
calculation the ``sad`` and ``sac`` guesses need and shares their cache. Because
neither side has seen the other atoms' functions there is no basis-set
superposition error to correct.

One caveat for a correlated wave function: the free atoms are still solved at the
unrestricted Hartree--Fock level, so a CASSCF molecule is being measured against
uncorrelated atoms and the formation energy inherits that imbalance. The
molecular correlation is included and the atomic correlation is not.

Reading the diagnostics
-----------------------

Two numbers at the end say whether to believe the tables::

    atomic character     0.9029   (6 refinement sweeps)
    valence gap          0.9977   against rejected   0.1807

**Atomic character** is how much of each quasi-atomic orbital stayed on its own
atom, one meaning the orbitals are exactly free-atom ones. The deficit is the
deformation the molecule imposes.

**Valence gap** is the separation the valence-virtual selection cut through.
Paper I reports 0.99999 against 0.105--0.272; anything much narrower means the
valence space was not cleanly separable and the analysis rests on a choice
rather than on a fact.

For correlated wave functions
-----------------------------

The analysis works on an MCSCF density as well as a reference determinant.
Nothing extra is asked for -- request it on a ``casscf`` deck and the correlated
density is what gets analysed.

The bond orders drop when it does, which is the point. N\ :sub:`2` in cc-pVDZ
gives a sigma and two pi bonds of order exactly 1.0000 at Hartree--Fock, and
0.9828, 0.9425, 0.9425 at CAS(6,6), because correlation puts occupation into the
antibonding orbitals.

One line appears that does not for a reference determinant::

    outside the valence space   0.003729 of 14 electrons

This is not an error. The quasi-atomic orbitals span the occupied valence space
plus the valence-virtual one and nothing else. A single determinant puts every
valence electron inside that span, so its populations sum exactly. A correlated
density does not, because the valence-virtual orbitals are chosen to look like
free-atom orbitals rather than to be natural orbitals. The number is the honest
measure of how much of the wave function the analysis is describing, and it is
printed rather than hidden.

Limits
------

- Closed-shell restricted references only.
- A correlated energy decomposition describes only the part of the wave
  function lying inside the quasi-atomic span, and reports the rest rather than
  hiding it. N\ :sub:`2` in cc-pVDZ with a CAS(6,6) accounts for about 89% of
  the correlation energy.
- The two-electron transformation goes through the dense ``n_ao**4`` integral
  array, which is the practical ceiling on molecule size for the energy tables.
- Hydrogen through xenon. Past that the free-atom minimal basis this projects
  onto would need a relativistic treatment that does not exist here, and the
  analysis refuses rather than using a basis that does not describe the atom.
- The output is printed, not written to the JSON output file.

References
----------

West, Schmidt, Gordon and Ruedenberg, *J. Chem. Phys.* **139**, 234107 (2013) --
the quasi-atomic orbitals and the valence-virtual space.

West, Schmidt, Gordon and Ruedenberg, *J. Phys. Chem. A* **119**, 10368 (2015)
-- the kinetic bond order and the orientation of the orbitals.

Del Angel Cruz, Gordon and Ruedenberg, *J. Am. Chem. Soc.* **147**, 42262
(2025) -- the intrinsic energy decomposition the energy tables implement.

The orbital labels and the thresholds that assign them follow GAMESS's
implementation rather than the papers, which define neither; the quantities
being labelled are the papers'. Validated against GAMESS to 1e-6 on formyl
chloride in 6-31G.
