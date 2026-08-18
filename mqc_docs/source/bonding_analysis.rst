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

The orbital labels and the thresholds that assign them follow GAMESS's
implementation rather than the papers, which define neither; the quantities
being labelled are the papers'. Validated against GAMESS to 1e-6 on formyl
chloride in 6-31G.
