Natural orbitals
================

The natural orbitals of a wave function are the orbitals that diagonalise its
one-particle density, and their eigenvalues are the occupations. A single
closed-shell determinant has occupations of exactly two and zero and its
natural orbitals say nothing its canonical ones did not. Every other wave
function has a density that is not idempotent, and for those the occupations
are the compact answer to two questions: how far is this from one determinant,
and in which orbitals.

.. code-block:: json

   "driver": "Energy",
   "properties": {
     "natural_orbitals": {}
   }

The object is the request and takes no settings. Like every entry under
``properties`` it changes no energy.

Which density
-------------

The most complete one-particle density the method has, and the output names it:

.. list-table::
   :header-rows: 1
   :widths: 30 25 45

   * - Method
     - ``density``
     - What is diagonalised
   * - Unrestricted Hartree-Fock or Kohn-Sham
     - ``unrestricted scf``
     - The sum of the two spin densities. These are the unrestricted natural
       orbitals (UNOs).
   * - ``ri-mp2``, closed shell
     - ``ri-mp2 relaxed``
     - The relaxed MP2 density: the reference's plus the correlation
       contribution, with the orbital response from the Z-vector equations.
       The density whose trace with an operator is the derivative of the MP2
       energy.
   * - Restricted Hartree-Fock or Kohn-Sham
     - ``scf``
     - The reference density. Twos and zeros; offered so a deck does not have
       to know in advance which kind of reference it will get.

Occupations of different densities are not comparable with one another, which
is why the name travels with the numbers.

What is reported
----------------

Every occupation is written to the output file, descending, under
``natural_orbitals.occupations``. The log lists the ones that are not within
0.001 of two or of zero, and three counts:

``unpaired_electrons_linear``
   :math:`\sum_i \min(n_i, 2 - n_i)`, Head-Gordon's number of effectively
   unpaired electrons. Zero for a closed-shell determinant and exactly the
   number of singly occupied orbitals for a high-spin one.

``unpaired_electrons_quadratic``
   :math:`\sum_i n_i^2 (2 - n_i)^2`. The same idea, weighted so that
   occupations near two and zero count at second order rather than first. The
   one to read when dynamic correlation is present, because the linear count
   adds up many small depletions into a number that looks like diradical
   character and is not.

``n_fractional``
   How many occupations lie between 0.02 and 1.98, the conventional window for
   choosing an active space. All occupations are written, so another window is
   a line of arithmetic on the output.

Two examples
------------

**A reference that has broken symmetry.** N\ :sub:`2` at 1.6 A in 6-31G, run
unrestricted with ``keywords.scf.second_order`` so that it leaves the
restricted solution (see :doc:`scf_convergence`)::

    natural orbitals of the unrestricted scf density
      electrons                            14.00000000
      effectively unpaired (linear)         2.90078945
      effectively unpaired (quadratic)      3.24003459
      6 orbitals with occupation between 0.02 and 1.98
      3 within 0.001 of two, 7 within 0.001 of zero; the rest:
        orbital    occupation
              4    1.99825998
              5    1.86737274
              6    1.34234745
              7    1.34234745
              8    0.65765255
              9    0.65765255
             10    0.13262726
             11    0.00174002

Six orbitals carry the fractional occupation -- a sigma pair at 1.87 and 0.13
and two degenerate pi pairs at 1.34 and 0.66 -- and each pair sums to two.
That is the triple bond coming apart, and it is the CAS(6,6) active space read
directly off a single-determinant calculation.

**Dynamic correlation only.** Water at RI-MP2/cc-pVDZ with a frozen core::

    natural orbitals of the ri-mp2 relaxed density
      effectively unpaired (linear)         0.19882036
      effectively unpaired (quadratic)      0.01620571
      5 orbitals with occupation between 0.02 and 1.98
      ...
              2    1.98701396
              3    1.97414226
              4    1.97069530
              5    1.96873704
              6    0.02285315
              7    0.02087790

Nothing is far from two or zero. The linear count of 0.20 is many small
depletions added up; the quadratic count of 0.016 is the one that says this
molecule is a closed shell.

A relaxed density is not the density of any wave function, so its occupations
are not bound to :math:`[0, 2]`. When one falls outside by more than
:math:`10^{-6}` the log says so and gives the extremes; here the highest is
2.0000025.

Validation
----------

The occupations agree with PySCF to :math:`10^{-10}`: against the generalised
eigenvalues of ``S D S`` for PySCF's unrestricted density on the water cation,
triplet O\ :sub:`2` and CH\ :sub:`3`, and for ``pyscf.mp.dfmp2_native``'s
relaxed density on water and ammonia with and without a frozen core. The
relaxed density itself is held to a finite-field derivative of the energy; see
:doc:`continuum_solvation`.

``test_mqc_czt_natural_orbitals`` pins what needs no other program: the
orbitals are orthonormal under the overlap and rebuild the density they came
from, a closed shell gives twos and zeros, and a high-spin determinant gives
one electron per open orbital.

Limits
------

- **MP2 means closed-shell** ``ri-mp2``. Conventional, spin-scaled and
  unrestricted MP2 have no relaxed density on this path and are refused, as is
  coupled cluster, rather than being handed the reference's occupations under
  a correlated heading.
- **Not under** ``keywords.pcm.correlation: "pte"``. The relaxed density's
  response equations do not carry the continuum. ``"ptd"`` is fine: its density
  is the gas-phase one.
- **Occupations only.** The orbitals are computed and not yet written
  anywhere; cube files and use as a starting guess are not wired.
- **CASSCF has its own.** An MCSCF deck already reports the natural
  occupations of its active space and does not go through this property.
- Near-linear dependence in the basis is removed before the density is
  diagonalised, so there can be fewer occupations than basis functions.
