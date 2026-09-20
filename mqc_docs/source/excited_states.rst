===============
Excited States
===============

Linear-response excitation energies out of a converged self-consistent field:
the Tamm-Dancoff approximation and the full Casida problem, singlets and
triplets, restricted and unrestricted, Hartree-Fock and Kohn-Sham. The
spectrum is asked for with ``keywords.excited_states`` on an ordinary
``"driver": "Energy"`` deck -- an excitation is not a derivative order, so it
does not change the driver.

What is computed
================

For each root: the excitation energy in Hartree and in eV, its spin, the
excited state's own total energy, the transition dipole in both the length and
the velocity gauge, the oscillator strength in both, and the leading natural
transition orbital weight. All of it lands in the JSON output under
``excited_states.states``, one object per root, and in an info-level table
while the run is going.

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Reference
     - Supported?
     - Notes
   * - Restricted Hartree-Fock
     - yes
     - TDHF; ``tda`` is CIS
   * - Restricted Kohn-Sham
     - yes
     - LDA, GGA, hybrid GGA, range-separated hybrid
   * - Unrestricted Hartree-Fock
     - yes
     - roots are not spin eigenstates; see below
   * - Unrestricted Kohn-Sham
     - yes
     - same functional coverage as the restricted path
   * - meta-GGA
     - **no**
     - refused; the kernel exists, no validated reference does
   * - VV10 non-local correlation
     - **no**
     - refused; both reference codes exclude it from the kernel
   * - Density-fitted reference
     - **no**
     - refused; the fitted build assumes an idempotent density
   * - Continuum solvation
     - **no**
     - refused; the solvent's response to a transition density is a
       term that does not exist here
   * - Hydrogen-capped fragment
     - **no**
     - refused; a cap's response has nowhere to be redistributed to
   * - MP2, coupled cluster, MCSCF reference
     - **no**
     - refused; that would be an equation-of-motion treatment
   * - cuEST (GPU) backend
     - **no**
     - refused; the solver is on the CPU backend

Every "no" above is a refusal by name, not an approximation. The reason is the
same in every case and is worth stating once: a response operator missing a
term still converges, and the spectrum it produces still looks like a
spectrum. Nothing in an excitation energy says which terms went into it.

Two methods
===========

``method: "tda"`` solves the Hermitian problem over the **A** block alone.
``method: "rpa"`` solves the full paired problem including the de-excitation
block **B**, through the Stratmann-Scuseria-Frisch reduction: the subspace
matrices ``b^T (A+B) b`` and ``b^T (A-B) b`` are formed, the second is
diagonalised and its square root taken, and the eigenvalues of
``S (A+B) S`` are :math:`\omega^2`. All the dense linear algebra is on the
subspace dimension, never on the full occupied-virtual space.

They are different approximations, not two routes to one answer -- RPA
energies sit a few millihartree below their Tamm-Dancoff partners -- so a
misspelled ``method`` is refused rather than resolved to whichever is nearer.

A negative eigenvalue of the projected ``(A-B)`` means the reference is not a
minimum along that rotation, and the solver stops and says so instead of
handing back a square root of a negative number. Triplet instabilities of a
closed-shell reference are the usual cause.

Singlets, triplets, and neither
===============================

Out of a **closed shell** the two spins are separate eigenproblems over the
same orbitals. ``spin: "singlet"`` and ``spin: "triplet"`` each solve one;
``spin: "both"`` solves both over one operator and one filled kernel cache and
**interleaves the results by energy**, so a ``both`` run with ``n_states: 3``
reports up to six roots and the ``spin`` field of each says which manifold it
came from. ``n_states`` counts roots per manifold, not in total.

The triplet two-electron part carries no Coulomb term -- it is exchange-only
-- and its exchange-correlation kernel is the spin-polarised functional
evaluated at :math:`\rho_\alpha = \rho_\beta = \rho/2`, giving
:math:`(f_{\alpha\alpha} - f_{\alpha\beta})/2`. A triplet transition moment is
written as an exact zero rather than computed: the spatial integral is
multiplied by an overlap of orthogonal spin functions.

Out of an **unrestricted** reference there is no such choice. The reference is
not a spin eigenfunction, so neither are its excitations; every root is
labelled ``"unrestricted"``, and a deck that names ``triplet`` or ``both``
over an unrestricted reference is refused rather than given the one spectrum
it has under a label it does not deserve. ``keywords.excited_states.spin`` is
simply not read there.

A doublet's response operator contains the rotation of its own half-filled
shell. In the paired problem that sits at :math:`\omega^2 = 0` and is dropped;
in the Tamm-Dancoff problem **A** alone is not singular along it, so it
survives as a small but nonzero root -- 6.7 millihartree on the OH radical.
Neither is a fault in the solver, and roots below 1e-3 Hartree are filtered as
rotations rather than reported as excitations.

Amplitude normalisation
=======================

Two conventions, one per route, and every consumer of an amplitude has to know
which it was handed.

* **Restricted:** :math:`\sum X^2 - \sum Y^2 = 1/2`, on all three of
  ``tda``, ``rpa`` and the Casida cross-check. A closed-shell excitation is two
  spin-orbital excitations of equal weight and the spatial amplitude carries
  both, so a pure single excitation prints a dominant amplitude of 0.707
  rather than 1.0.
* **Unrestricted:** :math:`\sum_\sigma (\sum X^2 - \sum Y^2) = 1`, on both
  routes. There is no closed-shell factor to place anywhere, so there is
  nothing for the two routes to disagree about.

The factor in a transition moment moves with the normalisation and not
separately. Restricted:

.. math::

   \mu = 2 \sum_{ia} \langle i | \mathbf{r} - \mathbf{R}_0 | a \rangle (X+Y)_{ia}

Unrestricted, where the spin sum the factor two stood for is written out and
**there is no factor two**:

.. math::

   \mu = \sum_{\sigma} \sum_{ia} \langle i_\sigma | \mathbf{r} - \mathbf{R}_0
         | a_\sigma \rangle (X+Y)^{\sigma}_{ia}

Halving one convention without the other leaves every excitation energy
exactly right and every oscillator strength wrong by a factor of four.

The transition-dipole origin
============================

:math:`\mathbf{R}_0` is the nuclear charge centroid,
:math:`\sum_A Z_A \mathbf{R}_A / \sum_A Z_A`, with each :math:`Z_A` the charge
the atom presents -- reduced where an effective core potential is in use.

The origin does not actually matter. A transition density carries no charge,
because :math:`\langle i | a \rangle` vanishes by orthogonality, so
:math:`\mu` is origin-independent to round-off. The centroid is a convention
rather than a choice, and it is PySCF's, which is what makes a cross-code
comparison a comparison of the same number.

Both gauges
===========

The **length** gauge is the expression above, with
:math:`f = \tfrac{2}{3}\,\omega\,|\mu|^2`.

The **velocity** gauge contracts :math:`\nabla` against :math:`X-Y` instead,
with :math:`f = \tfrac{2}{3}\,|v|^2 / \omega`. The sign is the one thing here
that is not obvious: ``int1e_ipovlp`` is
:math:`(\nabla \mu | \nu)`, the gradient on the *bra*, which is
:math:`-\langle \mu | \nabla | \nu \rangle`; the operator is anti-Hermitian
over real functions and integration by parts moves it across at the cost of a
sign. What is reported is the imaginary part of :math:`\langle 0 | p | n
\rangle`, component by component comparable with PySCF's.

The two agree only in a complete basis. On the lowest root of water in
cc-pVDZ they differ by a factor of three, which is a statement about the basis
and not about either implementation. Both are written out because the gap
between them is the diagnostic.

Natural transition orbitals
===========================

The singular value decomposition of the excitation amplitude ``X``,
renormalised to a unit vector, with ``Y`` dropped -- Martin's definition, and
PySCF's. The weights are the squared singular values and sum to one, so the
leading one says how nearly the root is a single orbital pair. Column phases
are fixed by making each column's largest-magnitude component positive.

An unrestricted root is decomposed per spin block, since the two blocks are
rectangles of different shapes over different orbitals. Each is divided by the
norm of the whole two-spin amplitude, so the two weight lists sum to one
between them rather than to one each; they are merged into one descending
column, and the reported leading weight is the largest of them, whichever
block it came from.

Which B3LYP
===========

``b3lyp`` resolves to libxc functional 402, the **VWN-RPA** flavour, which is
what PySCF and Psi4 both resolve the same name to. The VWN5 variant is a
different functional and gives different excitation energies. If a spectrum
from another program disagrees in the third decimal of an eV and nothing else
is different, this is the first thing to check.

The electronvolt printed beside each Hartree uses 27.211386245988. PySCF uses
27.21138602, which differs in the eighth decimal; every cross-code comparison
in the test suite and the validation manifest is done in Hartree for that
reason.

Convergence
===========

``tolerance`` is the residual at which a root is accepted, and it is **refused
below 1e-8** rather than clamped. On any grid a production run uses, the
exchange-correlation quadrature carries more error than that, so a tighter
request buys iterations and not accuracy; refusing rather than clamping means
a deck is never told it converged to a threshold it did not ask for.

The Davidson subspace is seeded with unit vectors on the lowest orbital-energy
gaps, degeneracy-aware: every pair within 1e-3 Hartree of the last one asked
for is included, so a degenerate manifold is not half-converged. New trial
vectors are handed to the operator as a block, which is one pass over the
integrals per iteration rather than one per vector.

Density screening is not simply disabled during the solve. A trial vector the
solver drives towards zero would make a density-keyed screen change the linear
map between one application and the next, which is not a linear operator at
all; the vector is rescaled instead, so the screen sees a density of ordinary
magnitude throughout.

Example
=======

.. code-block:: json

   {
     "schema": {"name": "water", "version": "1.0"},
     "molecules": [{"xyz": "water.xyz", "molecular_charge": 0,
                    "molecular_multiplicity": 1}],
     "model": {"method": "dft", "basis": "cc-pvdz", "functional": "b3lyp"},
     "driver": "Energy",
     "keywords": {
       "scf": {"tolerance": 1e-13, "gradient_tolerance": 1e-9},
       "dft": {"grid_level": 5},
       "excited_states": {
         "n_states": 5,
         "method": "rpa",
         "spin": "both",
         "tolerance": 1e-8
       }
     }
   }

See :doc:`input_files` for the full keyword list and
:doc:`developer_json_output` for the shape of what comes back.

What is not here
================

**Excited-state gradients.** There is no ``"driver": "Gradient"`` on top of a
spectrum, analytic or numerical, and asking for one gives a ground-state
gradient rather than an error. The Z-vector machinery it would need mostly
exists -- the coupled-perturbed solver, the second kernel derivative, the
non-symmetric Fock builds -- but the assembly does not.

**Spin-orbit coupling, magnetic dipoles, two-photon quantities, and state
averaging.** None of these are implemented and no keyword accepts them.

References
==========

Stratmann, Scuseria and Frisch, *J. Chem. Phys.* **109**, 8218 (1998), for the
subspace reduction the paired solver uses and for the stability test on
``(A-B)``.

Martin, *J. Chem. Phys.* **118**, 4775 (2003), for the natural transition
orbitals.

Casida, in *Recent Advances in Density Functional Methods*, Part I (1995), for
the response eigenproblem itself.
