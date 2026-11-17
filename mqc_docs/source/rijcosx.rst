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

**Three grids.** The first iterations run on a coarse grid, the rest on a
small one, and the final energy is one exchange build on a larger one, from
the converged density. The energy's error in the density is second order, so
the first two only have to converge the density, not integrate the energy.
``cosx_final_grid: false`` iterates on the larger grid throughout.

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

Wall time for the whole SCF on four threads, energies converged to 1e-8. The
machine was shared while these were taken, so read them to about ten per cent:

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
     - 2.8 s
   * - (H\ :sub:`2`\ O)\ :sub:`6`, def2-TZVP
     - 258
     - 29 s
     - 2.1 s
     - 6.1 s
   * - (H\ :sub:`2`\ O)\ :sub:`6`, def2-QZVP
     - 702
     - 559 s
     - 16 s
     - 43 s
   * - Triglycine, def2-TZVP
     - 469
     - 252 s
     - 15 s
     - 34 s
   * - (H\ :sub:`2`\ O)\ :sub:`20`, def2-SVP
     - 480
     - 153 s
     - 28 s
     - 30 s
   * - Decaglycine, def2-SVP
     - 734
     - 397 s
     - does not fit
     - 68 s

RIJCOSX is several times faster than the direct build everywhere, and more so
as the basis grows: thirteen times in def2-QZVP. Against the fully fitted
build it is slower on small systems, level by twenty waters, and the only one
of the two that runs once the fitted tensor no longer fits: decaglycine in
def2-SVP wants 15.8 GB for RI-JK and 2.5 GB here.

Where the time goes, and what was done about it:

* **The potential integral** is about two thirds of a run. It is evaluated by
  a kernel written for this one integral (``mqc_czt_esp_batch``):
  McMurchie-Davidson, every inner loop over the points of a batch, contracted
  in place. It agrees with the integral library to 1e-11 for every shell type
  through g and is about three and a half times faster. A basis with shells
  above g falls back to the library.
* **Three grids, not two.** The iterations before the commutator falls below
  1e-2 run on a coarser grid still, with a looser screening threshold.
* **The final energy is a trace.** The last build wants an energy and not a
  matrix, so it screens each shell pair on both shells where a matrix build
  can screen on one.
* **Only the functions that reach a batch** are evaluated and multiplied.
* **The Coulomb matrix** comes from a Coulomb-only fit that keeps the raw
  three-centre integrals over the screened pairs and never forms the fitted
  tensor.

Shell pairs are screened by their Schwarz bound times the density on the
grid: 1e-6 on the coarse grid, 1e-8 on the iteration grid and 1e-10 where the
energy is taken. Looser than that on the final grid is not worth it: what the
screening drops has one sign and adds up, and 1e-9 moves the energy of twenty
waters by 9e-6 hartree for a tenth off the time.

The integral is still the step a GPU backend would replace, behind
``cosx_esp_contract`` and ``esp_batch_contract``.

What it does not do
-------------------

Chain-of-spheres exchange is implemented for a **restricted Hartree-Fock or
global-hybrid DFT energy on the CPU backend**. Each of the following is refused
by name, not run with a different exchange:

* gradients, Hessians and geometry optimisation;
* unrestricted references;
* range-separated functionals and double hybrids;
* MP2, coupled cluster and other correlated methods;
* excited states, the stability analysis, the second-order SCF and quantum
  nuclei (``keywords.neo``);
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
  separates them. These are the ``cosx`` cases of the CPU validation suite:
  nine cases, within 5e-7 to 2.3e-5 hartree of PySCF at the
  default grid and 1.6e-7 at level 3.
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
