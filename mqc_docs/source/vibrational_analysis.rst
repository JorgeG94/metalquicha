========================================
Vibrational Analysis and Thermochemistry
========================================

This document describes the theoretical background and implementation details for computing
vibrational frequencies, IR intensities, and thermochemistry from the Hessian matrix in Metalquicha.

From Hessian to Vibrational Frequencies
=======================================

The Cartesian Hessian
---------------------

The Hessian matrix contains second derivatives of the energy with respect to atomic displacements:

.. math::

   H_{ij} = \frac{\partial^2 E}{\partial x_i \partial x_j}

For a system with N atoms, the Hessian is a 3N x 3N matrix in Cartesian coordinates,
with units of Hartree/Bohr².

Mass-Weighting
--------------

The Hessian must be mass-weighted before diagonalization. Each element is divided by
the square root of the product of the corresponding atomic masses:

.. math::

   \tilde{H}_{ij} = \frac{H_{ij}}{\sqrt{m_i \times m_j}}

where :math:`m_i` is the mass of the atom associated with coordinate *i* (in atomic mass units
converted to atomic units via ``AMU_TO_AU = 1822.888``).

Projection of Translational and Rotational Modes
------------------------------------------------

For an isolated molecule, 6 degrees of freedom (5 for linear molecules) correspond to
overall translation and rotation. These must be projected out to obtain pure vibrational frequencies.

**Translation Vectors (3 modes)**

.. math::

   T_x = (1,0,0, 1,0,0, \ldots) / \sqrt{N}

for all atoms, and similarly for :math:`T_y` and :math:`T_z`.

Mass-weighted: :math:`\tilde{T} = T \times \sqrt{m}`

**Rotation Vectors (3 modes for nonlinear, 2 for linear)**

For rotation about axis :math:`\alpha`, the displacement of atom *i* at position :math:`r_i`
from center of mass:

.. math::

   R_{\alpha,i} = (r_i - r_{COM}) \times \hat{e}_\alpha

where :math:`\hat{e}_\alpha` is the unit vector along axis :math:`\alpha`.

**Projection Procedure**

1. Orthonormalize the translation/rotation vectors using Gram-Schmidt
2. Construct projection matrix: :math:`P = I - \sum_k |v_k\rangle\langle v_k|`
3. Apply to mass-weighted Hessian: :math:`\tilde{H}_{proj} = P \times \tilde{H} \times P`

Diagonalization
---------------

Diagonalize the projected mass-weighted Hessian to obtain eigenvalues :math:`\lambda_i`
and eigenvectors (normal modes):

.. math::

   \tilde{H}_{proj} \times L = L \times \Lambda

where :math:`\Lambda` is diagonal with eigenvalues :math:`\lambda_i`.

Conversion to Frequencies
-------------------------

Eigenvalues have units of Hartree/(Bohr² x m_e). Convert to cm⁻¹:

.. math::

   \nu_i \text{ (cm}^{-1}\text{)} = \text{sign}(\lambda_i) \times \sqrt{|\lambda_i|} \times \text{AU\_TO\_CM1}

where ``AU_TO_CM1 = 2.642461 x 10^7`` is the conversion factor.

Negative eigenvalues indicate imaginary frequencies (saddle points), reported as negative cm⁻¹ values.

Additional Vibrational Properties
---------------------------------

**Reduced Masses**

.. math::

   \mu_i = \frac{1}{\sum_j (L_{ij})^2}

where the sum is over mass-weighted normal mode components. Units: amu.

**Force Constants**

.. math::

   k_i = 4\pi^2 c^2 \nu_i^2 \mu_i

Can be converted to mDyne/Å using ``AU_TO_MDYNE_ANG = 15.569141``.

**Cartesian Displacements**

Un-mass-weight the eigenvectors:

.. math::

   d_{ij} = \frac{L_{ij}}{\sqrt{m_j}}


IR Intensities from Dipole Derivatives
======================================

Theory
------

IR intensity is proportional to the square of the transition dipole moment, which depends
on how the dipole moment changes along the normal mode coordinate:

.. math::

   I_i \propto \left|\frac{\partial \mu}{\partial Q_i}\right|^2

Dipole Derivative Matrix
------------------------

The dipole derivative matrix (3 x 3N) contains derivatives of the dipole moment with
respect to Cartesian displacements:

.. math::

   \left(\frac{\partial \mu}{\partial x}\right)_{\alpha i} = \frac{\partial \mu_\alpha}{\partial x_i}

where :math:`\alpha \in \{x, y, z\}` and *i* runs over all 3N Cartesian coordinates.

There are two routes to it and the backend picks, the same way it picks between
the two Hessians.

**Analytically**, for a restricted Hartree-Fock or Kohn-Sham reference and for
the double hybrids -- not, despite the Hessian being analytic there too, for a
plain MP2 reference. What a dipole derivative needs beyond one-electron
integrals is the response of the density, and the coupled-perturbed solve has
just produced that for the Hessian -- so the extra cost is two one-electron
integrals rather than a second solve:

.. math::

   \frac{\partial \mu_\alpha}{\partial R_{X\beta}}
     = Z_X \delta_{\alpha\beta}
     - \mathrm{Tr}\!\left(D \frac{\partial r_\alpha}{\partial R_{X\beta}}\right)
     - \mathrm{Tr}\!\left(\frac{\partial D}{\partial R_{X\beta}} r_\alpha\right)

the three terms being the nucleus moving, the basis functions riding it, and the
electrons relaxing. The last is the one that is easy to omit and hard to notice:
without it the derivative is smooth, correctly signed, obeys the translational
sum rule, and is wrong by tens of percent.

A **double hybrid** adds a fourth: its dipole is the *relaxed* one, so the
perturbative term contributes through the MP2 relaxed density, and the
intensity needs that derivative on top of the reference's -- the same split as
its Hessian. On water/6-31G that term is 1.26e-02 against a dipole of order one.

**Plain MP2 has no infrared intensities**, by either route, and the reason is
the same relaxed density said from the other side. Analytically, the
perturbative term's dipole derivative is built only for a double hybrid --
`mp2_correlation_hessian` refuses the request without a functional rather than
returning the reference's. By finite difference, the backend leaves
`has_dipole` false for a correlated method, because what it could publish is
the *reference's* dipole and differencing that would report intensities from a
density the energy and gradient do not use. Both are refusals rather than
approximations: a missing intensity is noticed and a plausible wrong one is
not.

**By finite differences of the dipole** otherwise, collected at the same
displaced geometries the semi-numerical Hessian already visits.

The two agree: on water/STO-3G with B3LYP the analytic route gives 7.98, 56.46
and 29.96 km/mol and the finite-difference one gives 7.98, 56.47 and 29.98,
which is two entirely different assemblies landing within 0.02 km/mol of each
other.

Validated against PySCF, which has no infrared module but does have the pieces
one is built from -- ``int1e_irp`` and a coupled-perturbed solve. The dipole
derivative agrees to 5.8e-08 for Hartree-Fock, 1.9e-07 on a basis with d
functions, and 3.8e-08 for B3LYP on a converged grid. See
:doc:`analytic_hessians` for which calculations reach the analytic path.

Transformation to Normal Coordinates
------------------------------------

Transform dipole derivatives from Cartesian to normal mode coordinates:

.. math::

   \frac{\partial \mu}{\partial Q_i} = \sum_j \frac{\partial \mu}{\partial x_j} \times \frac{\partial x_j}{\partial Q_i}

The transformation matrix :math:`\partial x / \partial Q` comes from the mass-weighted normal modes:

.. math::

   \frac{\partial x_j}{\partial Q_i} = \frac{L_{ji}}{\sqrt{m_j}}

IR Intensity Calculation
------------------------

For each normal mode *i*:

.. math::

   I_i = \frac{N_A \pi}{3 \times 4\pi\varepsilon_0 \times c^2} \times \left|\frac{\partial \mu}{\partial Q_i}\right|^2

In atomic units, the conversion factor to km/mol is:

.. code-block:: fortran

   AU_TO_KMMOL = 1.7770969e6

The intensity is:

.. math::

   I_i \text{ (km/mol)} = \text{AU\_TO\_KMMOL} \times \sum_\alpha \left(\frac{\partial \mu_\alpha}{\partial Q_i}\right)^2


Raman Activities from Polarizability Derivatives
================================================

A mode is Raman active to the extent that the molecule's polarizability changes
along it, as it is infrared active to the extent its dipole does. Ask for the
activities beside the frequencies with

.. code-block:: json

   "model": {"method": "hf", "basis": "cc-pVDZ"},
   "driver": "Hessian",
   "properties": {"raman": {}}

``properties.raman`` is a request by its presence and is refused on any driver
but ``Hessian``. Its one setting, ``method``, chooses how the polarizability
derivatives are obtained: ``"field"`` (the default) or ``"displacements"``.

How it is computed
------------------

What is needed is the Cartesian derivative of the static polarizability,
:math:`\partial\alpha/\partial x_j`: a third derivative of the energy, twice in
the electric field and once in a nuclear coordinate. It is taken by finite
differences of something analytic, and there are two ways to arrange that.

**By field** (``"method": "field"``, the default). In a uniform field
:math:`\mathbf{F}` the energy is
:math:`E_0 - \boldsymbol\mu\cdot\mathbf{F} - \tfrac12\mathbf{F}\alpha\mathbf{F}`,
so

.. math::

   \frac{\partial \alpha_{kl}}{\partial x_j}
     = -\frac{\partial^2 g_j}{\partial F_k\,\partial F_l}

with :math:`g` the analytic nuclear gradient. The reference is reconverged in
twelve fields of 0.001 atomic units -- plus and minus each axis, and plus and
minus each pair of axes together -- and the gradient taken in each. **Thirteen
self-consistent fields and gradients, whatever the size of the molecule**, all
at one geometry, so the molecule and its quadrature grid are built once.

**By displacement** (``"method": "displacements"``). Each coordinate is moved
by 0.005 Bohr both ways, the reference reconverged, and the coupled-perturbed
polarizability evaluated there: :math:`6N` self-consistent fields and
:math:`6N` response solves. It costs more than the field route from three
atoms up and grows with the molecule, and is kept because it shares nothing
with the other past the reference: it is what the field route is checked
against.

The two are different finite differences of the same quantity and agree to the
error of each, about :math:`3\times10^{-5}` Bohr\ :sup:`2` on water's
derivatives of order ten.

Projected onto mode :math:`i`,

.. math::

   \alpha'_i = \sum_j \frac{\partial \alpha}{\partial x_j}\,\frac{L_{ji}}{\sqrt{m_j}}

with :math:`L` the mass-weighted eigenvectors, the same projection the infrared
intensities use. From its mean :math:`a = \mathrm{tr}(\alpha')/3` and its
anisotropy

.. math::

   \gamma^2 = \tfrac12\left[(\alpha'_{xx}-\alpha'_{yy})^2 + (\alpha'_{yy}-\alpha'_{zz})^2
              + (\alpha'_{zz}-\alpha'_{xx})^2\right]
              + 3\left(\alpha'^2_{xy} + \alpha'^2_{xz} + \alpha'^2_{yz}\right)

the **Raman scattering activity** is :math:`S = 45a^2 + 7\gamma^2`, in
:math:`\text{\AA}^4/\mathrm{amu}`, and the **depolarization ratio** for
linearly polarized light is :math:`\rho = 3\gamma^2/(45a^2 + 4\gamma^2)`. The
ratio runs from 0, a totally symmetric mode with no anisotropy, to 3/4, any
mode that is not totally symmetric.

What it reports
---------------

Two more rows in the vibrational table,

.. code-block:: text

    Frequencies --   1807.9869   3874.1341   3993.3734
    IR Intens  --    105.6565      1.9052     34.0227
    Raman Act  --     12.2241     93.9477     46.8233
    Depolar    --      0.4023      0.2072      0.7500

(water at Hartree-Fock/6-31G: the antisymmetric stretch is the one at 3/4), and
in the output file ``raman_activities_A4_amu``, ``raman_depolarization_ratios``
and the static polarizability itself, ``polarizability_bohr3``, under
``vibrational_analysis``.

**An activity is not an intensity.** It is a property of the molecule. What a
spectrometer records also depends on the laser frequency :math:`\nu_0` and the
temperature,

.. math::

   I_i \propto \frac{(\nu_0 - \nu_i)^4}{\nu_i\left[1 - e^{-h c \nu_i / k_B T}\right]}\, S_i

and applying that, and a line shape, is left to whoever plots the spectrum.

Validation
----------

The reference is PySCF used twice, sharing nothing with this program: its
analytic Hessian for the modes, and its own coupled-perturbed polarizability
differenced over the same displacements. Water at Hartree-Fock, B3LYP and PBE
in 6-31G and ammonia at Hartree-Fock in cc-pVDZ agree to
:math:`2\times10^{-5}\ \text{\AA}^4/\mathrm{amu}` on activities between 6 and
124, and to the five decimals printed on every depolarization ratio.
``test_mqc_czt_raman`` pins three of water's Cartesian derivatives to that
reference, and ``test_mqc_vibrational_analysis`` holds the activity formula to
three tensors whose answer is known by hand.

The error of the central differences themselves is measured by the
translational sum rule: the derivatives along one direction must cancel over
the atoms, and do to :math:`3\times10^{-5}` Bohr\ :sup:`2` on derivatives of
order ten.

Limits
------

- **Closed-shell Hartree-Fock and Kohn-Sham with exact integrals.** An open
  shell, a correlated method, a density-fitted reference, a continuum solvent,
  an effective core potential and capped or ghosted fragments are refused by
  name.
- **It accompanies the analytic Hessian.** Where the Hessian is taken by
  differences of gradients instead, the activities are skipped with a warning.
- **Static polarizability.** The derivative is of the zero-frequency
  polarizability, which is the usual approximation far from an electronic
  absorption and wrong near one: there is no resonance Raman here.
- **Numerical, not analytic.** Either route differences an analytic quantity.
  The fully analytic derivative needs no extra self-consistent fields at all
  and is not built; for Kohn-Sham it would need the third functional
  derivative.
- **The extra calculations run one after another**, not distributed over MPI
  ranks as the displacements of a semi-numerical Hessian are.
- **Harmonic, like everything else on this page.**

Thermochemistry (RRHO Approximation)
====================================

The Rigid Rotor Harmonic Oscillator (RRHO) approximation computes thermodynamic properties by treating:

- Translation as an ideal gas
- Rotation as a rigid rotor (classical limit)
- Vibration as quantum harmonic oscillators

Input Requirements
------------------

- Molecular geometry (coordinates, atomic numbers)
- Vibrational frequencies (cm⁻¹)
- Temperature T (default: 298.15 K)
- Pressure P (default: 1 atm)
- Symmetry number σ (default: 1)
- Spin multiplicity (default: 1)

Moments of Inertia
------------------

Compute the inertia tensor in the center-of-mass frame:

.. math::

   I_{\alpha\beta} = \sum_i m_i \times (r_i^2 \delta_{\alpha\beta} - r_{i\alpha} \times r_{i\beta})

Diagonalize to get principal moments :math:`I_A, I_B, I_C` (in amu·Å²).

Linear molecule detection: :math:`I_A \approx 0` (smallest moment below threshold).

Rotational Temperatures
-----------------------

.. math::

   \Theta_{rot} = \frac{h^2}{8\pi^2 I k_B} = \frac{24.2637}{I} \text{ [K, for I in amu·Å²]}

Zero-Point Energy (ZPE)
-----------------------

.. math::

   \text{ZPE} = \frac{1}{2} h c \sum_i \nu_i = \frac{1}{2} \sum_i \nu_i \times \text{CM1\_TO\_KELVIN} \times k_B

where ``CM1_TO_KELVIN = 1.4387773538277`` K/cm⁻¹.

Only positive (real) frequencies are included. Imaginary frequencies are skipped
with a warning.

The threshold for "imaginary" here is exactly zero, so a projected rotational
mode that comes back at ``-1e-5`` cm-1 rather than at zero is counted and
skipped. On a genuine minimum the reported count is therefore often 3 (or 5 for
a linear molecule) rather than 0, and those modes contribute nothing to the sums
anyway. The frequency table printed above the thermochemistry is what to read
for whether a mode is really imaginary: it shows an ``i`` only past a ``10``
cm-1 window, which is also what
:ref:`hess_end <hess-end>` uses when it decides whether an optimized geometry is
a minimum.

Translational Contributions
---------------------------

For an ideal gas:

.. math::

   E_{trans} &= \frac{3}{2} RT \\
   H_{trans} &= \frac{5}{2} RT \quad \text{(includes PV = RT)} \\
   S_{trans} &= R \left[\frac{5}{2} + \ln(q_{trans})\right] \quad \text{(Sackur-Tetrode)} \\
   C_{p,trans} &= \frac{5}{2} R

Partition function:

.. math::

   q_{trans} = \frac{V}{\lambda^3} = \left(\frac{2\pi m k T}{h^2}\right)^{3/2} \times \frac{kT}{P}

Rotational Contributions
------------------------

**Nonlinear molecules** (3 rotational DOF):

.. math::

   E_{rot} &= \frac{3}{2} RT \\
   S_{rot} &= R \left[\frac{3}{2} + \ln(q_{rot})\right] \\
   C_{v,rot} &= \frac{3}{2} R \\
   q_{rot} &= \frac{\sqrt{\pi}}{\sigma} \times \frac{T^{3/2}}{\sqrt{\Theta_A \times \Theta_B \times \Theta_C}}

**Linear molecules** (2 rotational DOF):

.. math::

   E_{rot} &= RT \\
   S_{rot} &= R \left[1 + \ln\left(\frac{T}{\sigma \Theta_{rot}}\right)\right] \\
   C_{v,rot} &= R \\
   q_{rot} &= \frac{T}{\sigma \Theta_{rot}}

Vibrational Contributions
-------------------------

For each mode with vibrational temperature :math:`\theta_v = 1.4388 \times \nu` (K):

Define :math:`u = \theta_v / T`:

.. math::

   E_{vib} &= R \sum_i \frac{\theta_{v,i}}{e^{u_i} - 1} \quad \text{[thermal only, excludes ZPE]} \\
   S_{vib} &= R \sum_i \left[\frac{u_i}{e^{u_i}-1} - \ln(1-e^{-u_i})\right] \\
   C_{v,vib} &= R \sum_i \frac{u_i^2 e^{u_i}}{(e^{u_i}-1)^2} \\
   q_{vib} &= \prod_i \frac{1}{1 - e^{-u_i}}

Electronic Contributions
------------------------

Ground state only:

.. math::

   E_{elec} &= 0 \\
   S_{elec} &= R \ln(2S+1) \quad \text{[spin multiplicity]}

Total Thermodynamic Functions
-----------------------------

.. math::

   \text{Thermal correction to Energy:} \quad E_{corr} &= \text{ZPE} + E_{trans} + E_{rot} + E_{vib} \\
   \text{Thermal correction to Enthalpy:} \quad H_{corr} &= E_{corr} + RT \\
   \text{Thermal correction to Gibbs:} \quad G_{corr} &= H_{corr} - T \times S_{total}

.. math::

   \text{Total Enthalpy:} \quad H &= E_{elec} + H_{corr} \\
   \text{Total Free Energy:} \quad G &= E_{elec} + G_{corr}


Implementation Notes
====================

Matching Reference Programs
---------------------------

During development, the thermochemistry output was validated against xtb (extended tight-binding).
Several adjustments were needed to match reference values:

**Rotational Temperature Constant**

- Issue: Rotational partition function and entropy were off by a factor of ~1.7
- Root cause: Incorrect rotational temperature constant
- Fix: Changed from ``16.8576 / I`` to ``24.2637 / I`` (K, for I in amu·Å²)

**Heat Capacity Convention**

- Issue: TR (translational) heat capacity showed 2.981 cal/K/mol instead of 4.968
- Root cause: Reported Cv instead of Cp
- Fix: For the TR row, report :math:`C_p = C_v + R = \frac{5}{2}R` instead of :math:`C_v = \frac{3}{2}R`

**Enthalpy Table (TOT row)**

- Issue: TOT enthalpy showed ~15000 cal/mol instead of ~2372
- Root cause: ZPE was included in the TOT row of the contribution table
- Fix: TOT row shows only thermal contributions (VIB + ROT + TR), ZPE is reported separately

**Thermal Correction Display**

- Issue: H(0)-H(T)+PV value was ~6x too large
- Root cause: ZPE was included in this quantity
- Fix: H(0)-H(T)+PV represents only the thermal correction without ZPE:
  :math:`E_{trans} + E_{rot} + E_{vib} + RT`

Key Constants
-------------

.. list-table:: Physical Constants Used
   :header-rows: 1
   :widths: 25 25 20 30

   * - Constant
     - Value
     - Units
     - Description
   * - ``AMU_TO_AU``
     - 1822.888
     - m_e/amu
     - Mass conversion
   * - ``AU_TO_CM1``
     - 2.642461 x 10^7
     - cm⁻¹
     - Frequency conversion
   * - ``CM1_TO_KELVIN``
     - 1.4387773538277
     - K/cm⁻¹
     - Vibrational temperature
   * - Θ_rot constant
     - 24.2637
     - K·amu·Å²
     - Rotational temperature
   * - ``R_CALMOLK``
     - 1.98720425864
     - cal/(mol·K)
     - Gas constant
   * - ``AU_TO_KMMOL``
     - 1.7770969 x 10^6
     - km/mol
     - IR intensity conversion


References
==========
