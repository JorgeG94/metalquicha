===================
Spin-spin couplings
===================

The splitting of an NMR line into a multiplet is the indirect coupling of two
nuclear moments through the electrons: one nucleus polarises the electrons
around it and the other feels the polarisation. Ask for the coupling constants
with

.. code-block:: json

   "model": {"method": "dft", "functional": "b3lyp", "basis": "cc-pVTZ"},
   "driver": "Energy",
   "properties": {"spin_spin_coupling": {}}

``properties.spin_spin_coupling`` is a request by its presence and takes no
settings. It is separate from ``properties.nmr``, the shieldings, and either
can be asked for without the other.

What is computed
----------------

Of Ramsey's four contributions to a coupling, **the Fermi-contact term only**:
the interaction of a nuclear moment with the electron spin density at the
nucleus. It is most of a coupling between hydrogens and between a hydrogen and
a carbon. The spin-dipole term and the two spin-orbit terms are not built; they
matter most for couplings that involve fluorine or a multiple bond between
heavy atoms.

The contact operator of a nucleus acts on the two electron spins with opposite
signs, so the response to it is a triplet one. The coupled-perturbed equations
are solved against the triplet operator, once per nucleus, and

.. math::

   K_{AB} = 4\alpha^4 \sum_{ai} h^A_{ai}\, U^B_{ai},
   \qquad h^A_{ai} = \frac{4\pi}{3}\,\phi_a(\mathbf{R}_A)\,\phi_i(\mathbf{R}_A)

It runs for closed-shell Hartree-Fock and Kohn-Sham with exact integrals. An
open shell, a correlated method, density fitting, a continuum solvent, an
effective core potential and capped or ghosted fragments are refused by name.

Three things decide whether the numbers are worth having:

- **Use a functional, not Hartree-Fock.** Hartree-Fock overestimates couplings
  badly and is often unstable towards a spin-polarised solution, in which case
  the triplet response has no solution and the run says so rather than
  returning one. A hybrid functional is the usual choice.
- **The basis has to describe the density at the nucleus.** The contact term
  reads the orbitals at a single point, where an ordinary contracted basis is
  at its worst. Basis sets made for this -- the pcJ family, or a standard set
  with its s functions uncontracted and tight ones added -- change a coupling
  by tens of per cent against cc-pVDZ. **None of them is shipped in**
  ``basis_sets/`` **yet**, so until one is, treat a coupling as good for the
  pattern of a multiplet and not for its width.
- **A coupling is a property of one geometry.** Three-bond couplings depend
  strongly on the dihedral angle, and what a spectrum shows is their average
  over the conformations the molecule visits.

What is written
---------------

The log lists every pair with its coupling in hertz. The output file carries
``spin_spin_coupling``:

.. code-block:: json

   "spin_spin_coupling": {
     "terms": ["fermi_contact"],
     "atoms": [
       {"atom": 1, "element": "O", "isotope": "17O", "g_factor": -0.757516},
       {"atom": 2, "element": "H", "isotope": "1H", "g_factor": 5.5856946893}
     ],
     "j_hz": [[0.0, -43.63], [-43.63, 0.0]],
     "reduced_1e19_T2_per_J": [[0.0, 26.79], [26.79, 0.0]]
   }

**A coupling in hertz belongs to a pair of isotopes.** ``j_hz`` is quoted for
the isotope an experiment observes -- :sup:`1`\ H, :sup:`13`\ C, :sup:`15`\ N,
:sup:`17`\ O, :sup:`19`\ F, :sup:`29`\ Si, :sup:`31`\ P and so on through
chlorine -- and ``atoms`` says which, with its g factor. Carbon is
:sup:`13`\ C, not the abundant and spinless :sup:`12`\ C. For another isotope,
scale by the ratio of g factors. ``reduced_1e19_T2_per_J`` is the reduced
coupling :math:`K`, which no isotope enters. An element heavier than chlorine
has no isotope in the table: its row of ``j_hz`` is zero and its ``isotope``
blank, and its reduced couplings are still there.

Both matrices are symmetric with a zero diagonal, a row per atom in the order
of ``atoms``.

:doc:`spectra` reads these with the shieldings and draws the multiplets.

Validation
----------

The reference is ``make_fc`` of ``pyscf.prop.ssc``, fed this repository's
basis sets. ``test_mqc_czt_spin_coupling`` pins the three reduced couplings of
water in 6-31G to it at Hartree-Fock and at B3LYP, to a relative
:math:`10^{-5}`; they agree to :math:`4\times10^{-8}` on the one-bond couplings
and :math:`1.5\times10^{-6}` on the small two-bond one. PBE was checked by hand
to the same.
