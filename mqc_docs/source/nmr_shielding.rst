Nuclear magnetic shielding
==========================

The shielding of a nucleus is how much the electrons around it weaken an
applied magnetic field at that nucleus. It is what a chemical shift is the
difference of: :math:`\delta = \sigma_\mathrm{ref} - \sigma`, with
:math:`\sigma_\mathrm{ref}` the same nucleus in a reference compound computed
the same way.

.. code-block:: json

   "model": {"method": "hf", "basis": "cc-pVDZ"},
   "driver": "Energy",
   "properties": {"nmr": {}}

``properties.nmr`` is a request by its presence and takes no settings. Every
nucleus is shielded. Like everything under ``properties`` it changes no energy.
It runs for closed-shell Hartree-Fock and for Kohn-Sham with an LDA, a GGA or
a global hybrid of either:

.. code-block:: json

   "model": {"method": "dft", "functional": "b3lyp", "basis": "cc-pVTZ"},
   "driver": "Energy",
   "properties": {"nmr": {}}

What it computes
----------------

The shielding tensor is a second derivative of the energy, with respect to the
applied field :math:`\mathbf{B}` and the magnetic moment
:math:`\boldsymbol{\mu}_N` of the nucleus:

.. math::

   \sigma^N_{xy} = \frac{\partial^2 E}{\partial B_x\,\partial \mu_{N,y}}

A magnetic field enters through a vector potential, which has an origin, and
in a finite basis the result depends on where that origin is put. The basis
functions here are **gauge-including atomic orbitals** (GIAOs, or London
orbitals): each carries a phase that moves its own origin to its own centre.
The shielding is then independent of the origin, and of where the molecule sits
in space, at any basis size.

The tensor has two parts. The *diamagnetic* part is an expectation value over
the unperturbed density. The *paramagnetic* part is the response of the
orbitals to the field, contracted with the nucleus's own operator, and costs one
coupled-perturbed solve for each of the three field directions. The field is an
imaginary perturbation, so the density it induces is antisymmetric: it has no
Coulomb response, and in Hartree-Fock only exchange couples the equations.

Only the sum of the two parts is a physical quantity. How it divides between
them depends on the gauge, and the split here is the GIAO one.

Under a functional
------------------

Two things change, and one of them makes a pure functional *cheaper* than
Hartree-Fock.

**The exchange-correlation potential does not respond to the field.** It is a
functional of the density, and a magnetic field induces a current, not a
change of density, at first order. So the potential contributes to the
first-order Fock matrix only through the gauge phase of the basis functions:
for a function centred at :math:`\mathbf{R}` the field derivative at a point
:math:`\mathbf{r}` is :math:`\tfrac{i}{2}(\mathbf{R}\times\mathbf{r})\chi`,
and the matrix is the potential between one function and another's
derivative, less its transpose, integrated on the functional's grid. A GGA
adds the gradient of that derivative.

**Exact exchange is scaled by the fraction the functional keeps**, in the
first-order Fock matrix and in the response alike. For a pure functional that
fraction is zero, nothing couples the three equations, and they are solved by
the orbital-energy gaps in one step with no response iterations at all.

That a functional of the density alone ignores the current is an
approximation of the functional, not of this code: it is what every
conventional functional does, and it is why shieldings from pure functionals
are too deshielded for nuclei with low-lying excited states.

What it reports
---------------

For each nucleus, the isotropic shielding (a third of the trace) and the
anisotropy, in parts per million:

.. code-block:: text

     ==== nuclear magnetic shielding (GIAO), ppm ====

        atom       isotropic    anisotropy
       O 1           346.7137       30.6279
       H 2            31.1325       19.6033
       H 3            31.1325       19.6033

That is water at Hartree-Fock/cc-pVDZ. The anisotropy is
:math:`\sigma_{33} - (\sigma_{11} + \sigma_{22})/2`, from the principal values
of the symmetrised tensor with :math:`\sigma_{33}` the one furthest from the
isotropic value. At the ``verbose`` logger level the full tensors are printed
too.

The output file carries everything under ``nmr_shielding``: for each atom its
number ``atom`` (from one) and ``element``, ``isotropic``, ``anisotropy``, the ``tensor`` as three rows, and its
``diamagnetic`` and ``paramagnetic`` parts.

**The tensor is not symmetric.** Its first index is the field's direction and
its second the nuclear moment's; a row of ``tensor`` is one field direction.
For a hydrogen of water the two in-plane off-diagonal elements differ by about
2 ppm. Only the symmetric part reaches a spectrum, which is why the summary
symmetrises, and the rows are written as computed.

Validation
----------

The reference is ``pyscf.prop.nmr.RHF``, fed this repository's basis sets.

``test_mqc_czt_nmr`` pins all nine elements of all three tensors of water in
6-31G to it, to the six decimals PySCF prints, and checks that moving the
molecule by several bohr changes no element by more than :math:`10^{-6}` ppm
(measured: :math:`3\times10^{-8}`). Through the program, the cases in the table
below agree with PySCF element by element.

.. list-table::
   :header-rows: 1
   :widths: 30 25 45

   * - Molecule
     - Basis
     - Largest difference in any tensor element (ppm)
   * - Water
     - cc-pVDZ
     - 1.2e-6
   * - Ammonia
     - cc-pVDZ
     - 8.2e-7
   * - Methane
     - cc-pVDZ
     - 9.8e-7
   * - HCN
     - cc-pVDZ
     - 2.7e-7
   * - Water
     - 6-31G* (Cartesian d)
     - 7.2e-7
   * - HF
     - cc-pVTZ
     - 3.8e-7
   * - Ethane
     - cc-pVDZ
     - 1.5e-7

Under a functional the reference is ``pyscf.prop.nmr.rks``. Its quadrature
loop fails under the current PySCF, so the numbers come from the same formulas
run in one pass, and that copy reproduces the values in the module's own
self-test -- hydrogen fluoride at LDA/6-31G and at the B3LYP of its day in
cc-pVDZ -- to :math:`10^{-4}` ppm. ``test_mqc_czt_nmr`` pins all three tensors
of water in 6-31G at PBE and at B3LYP to it, agreeing to
:math:`1.2\times10^{-6}` ppm, and checks that moving the molecule changes no
PBE element by more than :math:`10^{-5}` ppm (measured: :math:`8\times10^{-12}`).
LDA was checked through the program to the same.

Limits
------

- **Closed shells with exact integrals only.** An open shell, a correlated
  method, a density-fitted reference and a continuum solvent are each refused
  by name.
- **LDA, GGA and their global hybrids.** A meta-GGA is refused, since its
  kinetic energy density depends on the gauge and that dependence is not
  built; so are range-separated hybrids and functionals with non-local
  correlation.
- **No effective core potentials.**
- **The two-electron GIAO term is not screened.** Every shell quartet is
  computed, so this step grows as the fourth power of the basis size and is the
  slowest part of a shielding. Fine for tens of atoms in a double-zeta basis;
  it needs integral screening before it is used on anything large.
- **Hartree-Fock shieldings are a starting point, not a prediction.** Electron
  correlation changes them by tens of ppm for nuclei with lone pairs or multiple
  bonds, which is what a functional is for, and a basis built for energies
  converges them slowly. Chemical shifts, being differences, are more
  forgiving than the shieldings themselves.
- **Shieldings only.** No chemical shifts against a named reference, no
  spin-spin couplings, no magnetisability.
- **CPU backend.** The GPU and terco backends warn and skip it.
