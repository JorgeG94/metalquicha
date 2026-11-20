Intrinsic atomic and bond orbitals
==================================

Knizia's intrinsic atomic orbitals (IAOs) are a minimal set of atom-centred
orbitals -- one per function of a free-atom minimal basis -- that span the
occupied space of the molecule exactly. They are polarised by the molecule,
unlike the free-atom functions they start from, and they are orthonormal. Two
things follow that a Mulliken analysis cannot offer: there is no overlap
population to divide between atoms, and the numbers barely move with the
basis set of the calculation.

The intrinsic bond orbitals (IBOs) are the occupied orbitals rotated among
themselves so that each sits on as few atoms as possible, measured in IAO
charges. They come out as cores, lone pairs and bonds without anyone drawing a
Lewis structure first.

G. Knizia, *J. Chem. Theory Comput.* **9**, 4834 (2013).

.. code-block:: json

   "driver": "Energy",
   "properties": {
     "intrinsic_orbitals": {"exponent": 4}
   }

- ``exponent``: the power :math:`p` in the localisation functional
  :math:`\sum_i \sum_A (Q_A^i)^p`, where :math:`Q_A^i` is the IAO population of
  orbital :math:`i` on atom :math:`A`. ``4`` (the default) is Knizia's choice
  and gives the more local orbitals; ``2`` is Pipek-Mezey's functional on IAO
  charges. Nothing else is accepted.

Like everything under ``properties`` it changes no energy. It needs a
Hartree-Fock or Kohn-Sham reference, restricted or unrestricted.

What it reports
---------------

HCN at Hartree-Fock/cc-pVDZ::

    ==== intrinsic atomic and bond orbitals ====

       atom      IAO charge
      H 1         0.213857
      C 2        -0.013318
      N 3        -0.200539

      bond orbitals, exponent 4, 6 sweeps; electrons per atom, where above 0.10:
         1   1-centre:  N3 2.000
         2   1-centre:  C2 2.000
         3   2-centre:  C2 0.889  N3 1.099
         4   2-centre:  H1 0.768  C2 1.231
         5   1-centre:  N3 1.990
         6   2-centre:  C2 0.945  N3 1.055
         7   2-centre:  C2 0.945  N3 1.055

Two cores, the C--H bond, the nitrogen lone pair, and three orbitals shared
between carbon and nitrogen: the triple bond. The shares say how polar each
one is -- the C--H pair leans toward carbon, the C--N ones a little toward
nitrogen.

The output file carries the charges and the whole table under
``intrinsic_orbitals``: ``charges``, one per atom, and
``bond_orbital_populations``, one row per orbital with one entry per atom,
each row summing to one. The ``1-centre`` and ``2-centre`` labels in the log
count the atoms holding more than 0.05 of an orbital; the file has the numbers
and leaves the threshold to the reader. ``orbital_occupation`` says how many
electrons a row stands for: 2 here, 1 for an unrestricted reference.

An unrestricted reference
-------------------------

The IAOs are built from an occupied space, and an unrestricted wave function
has two. So with ``model.unrestricted`` the analysis is done once per spin:
alpha's occupied orbitals give alpha's IAOs and are localised among
themselves, and beta's likewise. The charges are the two spins' populations
added and taken from the nuclear charge; their difference is a spin
population per atom, which a restricted reference does not have.

The hydroxyl radical at Hartree-Fock/cc-pVDZ::

       atom      IAO charge        spin
      O 1        -0.362761    1.043722
      H 2         0.362761   -0.043722

      bond orbitals, exponent 4, 3 sweeps; electrons per atom, where above 0.05:
      alpha
         1   1-centre:  O1 1.000
         2   2-centre:  O1 0.703  H2 0.297
         3   1-centre:  O1 1.000
         4   1-centre:  O1 1.000
         5   1-centre:  O1 1.000
      beta
         1   1-centre:  O1 1.000
         2   2-centre:  O1 0.660  H2 0.340
         3   1-centre:  O1 1.000
         4   1-centre:  O1 1.000

Each orbital now holds one electron. Alpha has one more orbital on oxygen than
beta, which is the unpaired electron, and the O--H bond is not the same
orbital in the two spins: beta's leans further toward hydrogen, which is where
the small negative spin population on hydrogen comes from.

The output file gains ``spin_populations`` and
``bond_orbital_populations_beta``; ``bond_orbital_populations`` holds alpha's.
Alpha and beta orbitals are localised independently and listed in no common
order, so row *n* of one is not the partner of row *n* of the other.

A closed shell run unrestricted gives the restricted charges, zero spin
populations and the same orbitals twice.

The charges do not move with the basis
--------------------------------------

Water at Hartree-Fock, charge on oxygen, same geometry throughout:

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 20

   * - Basis
     - Mulliken
     - Lowdin
     - IAO
   * - 6-31G
     - -0.800
     - -0.587
     - -0.703
   * - cc-pVDZ
     - -0.314
     - -0.483
     - -0.704

Mulliken's charge more than halves between the two bases and Lowdin's moves by
a tenth. The intrinsic charge moves by 0.001.

The minimal basis
-----------------

The IAOs are defined against a free-atom minimal basis, and the numbers
depend on which. This implementation uses **MINAO**, the leading contracted
functions of cc-pVTZ, which is the set Knizia's paper uses and the default of
``pyscf.lo.iao``. It is ``basis_sets/minao/iao.json``, written by
``tools/minao/gen_iao_minimal_basis.py``.

It is *not* the basis the ``minao`` initial guess uses
(``basis_sets/minao/minao.json``, the leading ANO-RCC contractions). Two
different minimal bases share that name upstream, and they give different
intrinsic charges: on the water of the unit test, in 6-31G, oxygen is -0.699
against MINAO and -0.668 against the ANO-RCC set.

Hydrogen through krypton, except potassium, which the upstream table does not
carry.

Validation
----------

``test_mqc_czt_iao`` pins the IAO charges, the localisation functional at both
exponents and the oxygen share of the O--H bond orbitals of water in 6-31G to
``pyscf.lo.iao`` and ``pyscf.lo.ibo``. Through the program, six cases agree
with PySCF: water, ammonia, methane and HCN at Hartree-Fock, water at B3LYP,
and water at exponent 2, with charges within :math:`10^{-8}` and the functional
within :math:`5\times10^{-9}`. The orbitals themselves are determined only up
to order and sign, and two lone pairs on one atom only up to a rotation
between them, so it is the functional and the populations that are compared.

Unrestricted, the suite pins the water cation's charges, spin populations and
each spin's functional to the same two PySCF routines given one spin's
occupied orbitals at a time, and checks that a closed shell run unrestricted
reproduces the restricted analysis. Through the program, the hydroxyl and
methyl radicals, triplet oxygen and the water cation at Hartree-Fock, methyl at
B3LYP and hydroxyl at exponent 2 agree with PySCF to
:math:`1.2\times10^{-8}` in the charges and spin populations and
:math:`4\times10^{-9}` in the functional.

The same suite holds what needs no reference: the IAOs are orthonormal under
the overlap, the occupied orbitals lie inside their span, localising leaves
the density unchanged and does not lower the functional, and every bond
orbital's atomic shares sum to one.

Limits
------

- **Hartree-Fock or Kohn-Sham only.** A correlated wave function has no set
  of occupied orbitals to build the IAOs from, and is refused.
- **No effective core potentials.** The minimal basis would still describe the
  core electrons the potential replaced.
- **Hydrogen to krypton, without potassium.**
- **Populations only.** The orbitals are computed and not written anywhere: no
  cube files yet.
- **CPU backend.**
