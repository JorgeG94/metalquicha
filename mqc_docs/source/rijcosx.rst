RIJCOSX: chain-of-spheres exchange
==================================

RIJCOSX builds the two halves of a Hartree-Fock or hybrid-DFT Fock matrix by
two different approximations. The Coulomb matrix comes from density fitting
(RI-J), and the exchange matrix from a seminumerical quadrature, the
chain-of-spheres algorithm (COSX). It is the default way to run hybrid DFT in
ORCA, where it originates [Neese2009]_.

The keyword reference lives with the other keywords in :doc:`input_files`
(``keywords.scf.exchange``, ``cosx_grid`` and ``cosx_final_grid``). This page
is what the method does, how accurate and how fast it is here, and what it does
not yet do.

.. code-block:: json

   {
     "model": {"method": "hf", "basis": "def2-svp",
               "aux_basis": "def2-universal-jkfit"},
     "keywords": {"scf": {"density_fitting": true, "exchange": "cosx"}},
     "driver": "Energy"
   }

What it computes
----------------

One electron of every two-electron integral is integrated on a grid, and the
other analytically:

.. math::

   (\mu\kappa|\nu\lambda) \approx \sum_g w_g\, \chi_\mu(\mathbf r_g)\,
   \chi_\kappa(\mathbf r_g)\, A^{\nu\lambda}_g,
   \qquad
   A^{\nu\lambda}_g = \int \frac{\chi_\nu(\mathbf r)\chi_\lambda(\mathbf r)}
   {|\mathbf r - \mathbf r_g|}\, d\mathbf r

so the exchange matrix is three steps for each batch of grid points:

.. math::

   F_{g\lambda} = w_g \sum_\kappa \chi_\kappa(\mathbf r_g) P_{\kappa\lambda},
   \qquad
   G_{g\nu} = \sum_\lambda A^{\nu\lambda}_g F_{g\lambda},
   \qquad
   K_{\mu\nu} = \sum_g \chi_\mu(\mathbf r_g)\, G_{g\nu}

followed by a symmetrisation, because the quadrature treats the two electrons
differently. The first and last steps are matrix products. The middle one is
where the time goes: an electrostatic-potential integral for every shell pair
at every grid point.

**Overlap fitting.** A grid small enough to be worth using integrates products
of basis functions poorly. Each grid therefore carries one matrix,
:math:`\mathbf S\,\mathbf S_\text{grid}^{-1}`, that makes the same quadrature
exact for the overlap [Izsak2011]_. Applied to :math:`\mathbf K` it reduces the
error of the exchange energy on the smallest grid several-fold.

**Two grids.** The iterations run on a small grid and the final energy is one
exchange build on a larger one, from the converged density. The energy's error
in the density is second order, so the small grid only has to converge the
density, not integrate the energy. ``cosx_final_grid: false`` iterates on the
larger grid throughout.

Accuracy
--------

An RIJCOSX energy differs from the unapproximated one by two things: the fit
of the Coulomb matrix and the quadrature of the exchange.

**The quadrature.** Exchange energy against the analytic exchange, on a
converged density, on each level's final grid. Worst case over water, the
water dimer, CH\ :sub:`4`, NH\ :sub:`3`, H\ :sub:`2`\ S, HCl, ethane and HCN in
def2-SVP and def2-TZVP:

.. list-table::
   :header-rows: 1

   * - ``cosx_grid``
     - Worst error (kcal/mol)
     - Final grid, points per atom
   * - 1
     - 0.27
     - about 1,900
   * - 2 (default)
     - 0.02
     - about 3,500
   * - 3
     - 0.015
     - about 6,200

**The total.** Total energy against Hartree-Fock with no approximation,
default grid, ``def2-universal-jkfit``:

.. list-table::
   :header-rows: 1

   * - System
     - Basis
     - RI-JK (hartree)
     - RIJCOSX (hartree)
   * - (H\ :sub:`2`\ O)\ :sub:`6`
     - def2-SVP
     - +3.8e-4
     - -1.4e-4
   * - (H\ :sub:`2`\ O)\ :sub:`6`
     - def2-TZVP
     - +9.8e-5
     - -1.7e-4
   * - (H\ :sub:`2`\ O)\ :sub:`6`
     - def2-QZVP
     - +1.1e-4
     - -1.9e-4
   * - Triglycine
     - def2-TZVP
     - +4.0e-4
     - -2.3e-4
   * - (H\ :sub:`2`\ O)\ :sub:`20`
     - def2-SVP
     - +1.2e-3
     - -4.5e-4

Most of the RIJCOSX error is the Coulomb fit, not the grid.

.. warning::

   **Use an auxiliary basis that fits the Coulomb matrix**: a ``-jfit`` or
   ``-jkfit`` set. With an ``-rifit`` set RIJCOSX is off by 2.4 millihartree on
   water in cc-pVDZ. Such a set fits the Coulomb matrix poorly, and RI-JK only
   looks better with it because its Coulomb and exchange errors come from the
   same fit and largely cancel. RIJCOSX has no exchange fit to cancel against.

The energy is not variational: it can land above or below the exact one. The
overlap fit also means the converged orbitals are stationary only to within
the quadrature error, which is inherent to fitted COSX.

Speed
-----

Wall time for the whole SCF on four threads, energies converged to 1e-8:

.. list-table::
   :header-rows: 1

   * - System
     - Basis functions
     - Direct HF
     - RI-JK
     - RIJCOSX
   * - (H\ :sub:`2`\ O)\ :sub:`6`, def2-SVP
     - 144
     - 5.7 s
     - 0.6 s
     - 11 s
   * - (H\ :sub:`2`\ O)\ :sub:`6`, def2-TZVP
     - 258
     - 29 s
     - 2.1 s
     - 24 s
   * - (H\ :sub:`2`\ O)\ :sub:`6`, def2-QZVP
     - 702
     - 559 s
     - 16 s
     - 127 s
   * - Triglycine, def2-TZVP
     - 469
     - 252 s
     - 15 s
     - 126 s
   * - (H\ :sub:`2`\ O)\ :sub:`20`, def2-SVP
     - 480
     - 153 s
     - 28 s
     - 166 s

**RIJCOSX here is not yet faster than RI-JK, and you should not pick it for
speed.** It overtakes the direct build as the basis grows -- four times faster
in def2-QZVP -- but the fully fitted build is several times faster again
wherever its three-index tensor fits in memory. About nine tenths of the
RIJCOSX time is the potential integral, evaluated through the general-purpose
routine of the integral library at roughly 55 ns per shell pair per grid
point. A kernel written for this one integral is what closes the gap, and it
is the step a GPU backend replaces.

Two things limit it further today:

* The Coulomb half still stores the full three-index fitting tensor, so
  RIJCOSX needs as much memory as RI-JK.
* Shell pairs are screened by their Schwarz bound times the density on the
  grid, at 1e-10. That is worth a factor of 2.5 on twenty waters. The threshold
  cannot be loosened much: what it drops has one sign and adds up.

What it does not do
-------------------

Chain-of-spheres exchange is implemented for a **restricted Hartree-Fock or
global-hybrid DFT energy on the CPU backend**. Each of the following is refused
by name, not run with a different exchange:

* gradients, Hessians and geometry optimisation;
* unrestricted references;
* range-separated functionals and double hybrids;
* MP2, coupled cluster and other correlated methods;
* excited states, the stability analysis and the second-order SCF;
* Fukui functions and the bonding energy decomposition;
* FMO and EFMO fragment SCFs;
* the cuEST and terco backends.

A pure functional has no exact exchange, so ``exchange: "cosx"`` changes
nothing for one.

How it was validated
--------------------

* **The exchange matrix** against the direct four-index build on the same
  density: it converges as the grid grows, to 3e-10 on water in 6-31G*. A wrong
  factor or a transposed index leaves an error no grid removes.
* **The exact parts to rounding**: the contraction against the dense potential
  integral, independence of the batch size, the fitted overlap equal to the
  analytic one.
* **A density of no symmetry** against the permutation-explicit direct build,
  which fixes the index order.
* **The SCF** against PySCF with the same fitted Coulomb matrix and exact
  exchange, so that the fit is common to both numbers and only the quadrature
  separates them. These are the ``cosx`` cases of the CPU validation suite.
* **The hybrid fraction**: PBE0's grid error is a quarter of Hartree-Fock's,
  and PBE is unchanged to 1e-13.

The grids are tuned against the analytic exchange on the molecules above. They
follow the structure of ORCA's [HelmichParis2021]_ -- Lebedev orders pruned by
region, radial counts growing with the period -- but not its optimised
parameters, which are not published.

References
----------

.. [Neese2009] F. Neese, F. Wennmohs, A. Hansen and U. Becker, *Chem. Phys.*
   **356**, 98 (2009).
.. [Izsak2011] R. Izsák and F. Neese, *J. Chem. Phys.* **135**, 144105 (2011).
.. [HelmichParis2021] B. Helmich-Paris, B. de Souza, F. Neese and R. Izsák,
   *J. Chem. Phys.* **155**, 104109 (2021).
