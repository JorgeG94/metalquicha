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
closed-shell Hartree-Fock or Kohn-Sham reference.

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
and leaves the threshold to the reader.

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

The same suite holds what needs no reference: the IAOs are orthonormal under
the overlap, the occupied orbitals lie inside their span, localising leaves
the density unchanged and does not lower the functional, and every bond
orbital's atomic shares sum to one.

Limits
------

- **Closed-shell Hartree-Fock or Kohn-Sham only.** An unrestricted reference
  would need one set of IAOs per spin, and a correlated one has no single set
  of occupied orbitals; both are refused.
- **No effective core potentials.** The minimal basis would still describe the
  core electrons the potential replaced.
- **Hydrogen to krypton, without potassium.**
- **Populations only.** The orbitals are computed and not written anywhere: no
  cube files yet.
- **CPU backend.**
