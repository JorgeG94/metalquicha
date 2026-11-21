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

Limits
------

- **Closed-shell Hartree-Fock with exact integrals only.** An open shell, a
  Kohn-Sham reference, a correlated method, a density-fitted reference and a
  continuum solvent are each refused by name. A Kohn-Sham shielding needs the
  gauge phase carried through the exchange-correlation quadrature, which is the
  next step and is not built.
- **No effective core potentials.**
- **The two-electron GIAO term is not screened.** Every shell quartet is
  computed, so this step grows as the fourth power of the basis size and is the
  slowest part of a shielding. Fine for tens of atoms in a double-zeta basis;
  it needs integral screening before it is used on anything large.
- **Hartree-Fock shieldings are a starting point, not a prediction.** Electron
  correlation changes them by tens of ppm for nuclei with lone pairs or multiple
  bonds, and a basis built for energies converges them slowly. Chemical shifts,
  being differences, are more forgiving than the shieldings themselves.
- **Shieldings only.** No chemical shifts against a named reference, no
  spin-spin couplings, no magnetisability.
- **CPU backend.** The GPU and terco backends warn and skip it.
