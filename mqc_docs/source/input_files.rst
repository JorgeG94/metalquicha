.. _input_files:

==================
Input File Formats
==================

(this file was partially generated with an LLM but carefully checked by me, Jorge)

Metalquicha reads JSON input files.

.. contents::
   :local:
   :depth: 2

Overview
========

The workflow is:

1. Create a JSON file with your calculation setup
2. Run metalquicha with it

.. note::

   **Changed in 0.2.0.** Earlier versions read a separate section-based
   ``.mqc`` format, generated from your JSON by a ``mqc_prep.py`` helper.
   That intermediate step is gone: ``mqc`` reads the JSON directly, and both
   the ``.mqc`` format and ``mqc_prep.py`` have been removed. Your existing
   JSON inputs work unchanged -- drop the conversion step and pass the
   ``.json`` file where you used to pass the ``.mqc`` one. If you hand-wrote
   ``.mqc`` files, see :ref:`migrating_from_mqc` below.

JSON Input Format
=================

Complete JSON Schema
--------------------

Here is a complete example with all available options:

.. code-block:: json

   {
     "schema": {
       "name": "mqc-frag",
       "version": "1.0"
     },
     "molecules": [{
       "xyz": "path/to/geometry.xyz",
       "molecular_charge": 0,
       "molecular_multiplicity": 1,
       "fragments": [[0,1,2], [3,4,5]],
       "fragment_charges": [0, 0],
       "fragment_multiplicities": [1, 1]
     }],
     "model": {
       "method": "XTB-GFN1",
       "basis": "cc-pVDZ",
       "aux_basis": "cc-pVDZ-RIFIT"
     },
     "keywords": {
       "scf": {
         "maxiter": 300,
         "tolerance": 1e-6
       },
       "fragmentation": {
         "method": "MBE",
         "allow_overlapping_fragments": true,
         "level": 2,
         "embedding": "none",
         "cutoff_method": "distance",
         "distance_metric": "min",
         "cutoffs": {
           "dimer": 10.0,
           "trimer": 8.0
         }
       }
     },
     "system": {
       "logger": {
         "level": "Verbose"
       }
     },
     "driver": "Energy"
   }

Backend Selection
-----------------

Which integral backend runs, as a root-level key beside ``driver``:

.. code-block:: json

   "backend": "cuest"

- ``auto`` (default): cuEST when the build has it, the CPU path otherwise. This
  is the historical behaviour and what every existing deck gets.
- ``cuest``, or ``gpu``: require the GPU backend. A build without cuEST is
  **refused**, with the CMake flag that would fix it, rather than falling back --
  a deck that asked for cuEST and silently ran on the CPU would report a
  provenance and a timing that were not true. Asking for it alongside MP2 or
  coupled cluster is refused too, since those have no GPU implementation here.
- ``libcint``, or ``cpu``: force the CPU path even on a build that has cuEST,
  which is how the two are compared against each other.

An unrecognised name is refused rather than treated as ``auto``, so a typo
cannot quietly select a different implementation.

Schema Section
--------------

Required section that identifies the input format:

.. code-block:: json

   "schema": {
     "name": "mqc-frag",
     "version": "1.0"
   }

- ``name``: Must be ``"mqc-frag"``
- ``version``: Currently ``"1.0"``

Molecules Section
-----------------

Defines the molecular system(s) to calculate. Can contain multiple molecules for conformer/isomer studies.

Geometry Specification
^^^^^^^^^^^^^^^^^^^^^^

**Option 1: External XYZ file (recommended)**

.. code-block:: json

   "molecules": [{
     "xyz": "path/to/geometry.xyz",
     "molecular_charge": 0,
     "molecular_multiplicity": 1
   }]

**Option 2: Inline geometry**

.. code-block:: json

   "molecules": [{
     "geometry": {
       "symbols": ["H", "O", "H"],
       "coordinates": [
         [0.0, 0.0, 0.0],
         [0.0, 0.0, 0.96],
         [0.76, 0.0, -0.24]
       ]
     },
     "molecular_charge": 0,
     "molecular_multiplicity": 1
   }]

Fragment Definition
^^^^^^^^^^^^^^^^^^^

For fragmented calculations, specify which atoms belong to each fragment:

.. code-block:: json

   "fragments": [
     [0, 1, 2],      // Fragment 1: atoms 0, 1, 2
     [3, 4, 5],      // Fragment 2: atoms 3, 4, 5
     [6, 7, 8]       // Fragment 3: atoms 6, 7, 8
   ],
   "fragment_charges": [0, 0, 0],
   "fragment_multiplicities": [1, 1, 1]

**Notes**:

- Atom indices are **0-based** (first atom is 0)
- Fragment charges must sum to ``molecular_charge``
- Fragment multiplicities must be consistent with ``molecular_multiplicity``
- Fragments can overlap if ``allow_overlapping_fragments: true``

Connectivity (Optional)
^^^^^^^^^^^^^^^^^^^^^^^

For hydrogen capping at broken bonds:

.. code-block:: json

   "bonds": [
     {"atom_i": 2, "atom_j": 3, "bond_order": 1, "is_broken": true},
     {"atom_i": 5, "atom_j": 6, "bond_order": 1, "is_broken": true}
   ]

When a bond is marked ``is_broken: true``, metalquicha adds hydrogen caps at the break points.

Model Section
-------------

Specifies the quantum chemistry method:

.. code-block:: json

   "model": {
     "method": "XTB-GFN1",
     "basis": "cc-pVDZ",
     "aux_basis": "cc-pVDZ-RIFIT"
   }

**Supported methods**. Spelling is case-insensitive, and an ``XTB-`` prefix is
stripped before matching, so ``XTB-GFN2`` and ``gfn2`` are the same request.

Semi-empirical, through tblite:

- ``GFN1`` (also ``GFN1-xTB``, ``XTB-GFN1``)
- ``GFN2`` (also ``GFN2-xTB``, ``XTB-GFN2``)

Ab initio, on the CPU through libcint:

- ``HF`` (also ``RHF``, ``UHF``, ``Hartree-Fock``) -- an odd electron count or a
  multiplicity above one selects the unrestricted path whatever the spelling says
- ``MP2``, and ``SCS-MP2`` / ``SOS-MP2`` for the spin-component-scaled variants
- ``CCSD`` and ``CCSD(T)``
- ``RI-`` or ``DF-`` on any correlated method asks for the density-fitted route:
  ``RI-MP2``, ``RI-CCSD``, ``RI-CCSD(T)``
- ``CASSCF`` (also ``MCSCF``) and ``CASCI`` -- one method with the orbitals free
  or frozen; both need ``keywords.mcscf`` to say what the active space is

**Basis sets.** ``basis`` is required by every ab initio method and ignored by the
semi-empirical ones, which carry their own parameters.

``aux_basis`` is the **only** place an auxiliary basis is named, and it serves both
a density-fitted reference and a density-fitted correlation treatment. A basis set
belongs beside the orbital basis it fits; having two places to name one meant a
deck could set both and silently prefer one.

One set therefore covers both fits when a run does both. That is fine in the
direction it usually runs -- a RIFIT set fitting J and K is ordinary practice,
worth about 1.7 mHartree on a total energy against exact J and K and largely
cancelling in anything relative. The reverse is worth a warning and gets one: a
JKFIT set fitting a ``(ia|jb)`` block gives a correlation energy whose error is not
the RI error it is meant to be. Naming an auxiliary basis does **not** by itself
density-fit the reference -- that is ``keywords.scf.density_fitting``, asked for
rather than inferred.

``cartesian`` reads the basis in Cartesian form whatever its file declares --
6d rather than 5d, 10f rather than 7f. It defaults to false, meaning the file
decides, and it is a no-op for a basis with no shell above p, where the two
forms are the same functions.

It is not a preference. For a shell above p the two forms span different
spaces, so this changes the answer rather than its representation: water at
cc-pVDZ/B3LYP goes from 24 basis functions and -76.420342 Hartree to 25 and
-76.421536. Both are correct, and they are correct for different models.

It exists because the convention is not something a basis set *name* settles.
BSE is inconsistent about it -- 6-31G* is marked Cartesian and cc-pVDZ
spherical -- and GAMESS assumes Cartesian for the Pople sets throughout, so
reproducing a GAMESS number, or comparing against a program whose default runs
the other way, needs a way to say which was meant. Without it the run still
succeeds and still converges; it just answers a different question than the one
being compared against.

An auxiliary basis follows the orbital basis rather than its own file: a
three-centre integral is built in one angular form, so a spherical fitting
basis over a Cartesian orbital basis is not a mixture the integrals can
express.

It is a CPU-path keyword. cuEST builds its AO shells spherical whatever the
basis says, so asking for both is refused rather than quietly answered in the
other form -- ``backend: cuest`` with ``cartesian`` on is an error, and so is
leaving the backend at ``auto`` on a build that resolves it to the GPU. Ask for
``backend: libcint``, which honours the file and this keyword both. A basis
whose *file* is Cartesian, such as 6-31G*, was already refused on the GPU path
for the same reason.

Kohn-Sham DFT, on the CPU through libcint and libxc:

- ``DFT`` (also ``KS``, ``Kohn-Sham``) selects the method; **which** functional is
  ``model.functional``, not part of the method name. The two are separate fields
  because a functional names the theory the way a basis names the space it is
  solved in.

**Not yet reachable**: ``MCSCF`` and the F12 variants parse but have no
implementation. Unrestricted MP2 and unrestricted coupled cluster are refused
rather than quietly run restricted: both transforms take one set of orbitals and
an occupied count, so an open-shell reference needs separate alpha and beta
transforms. That also rules out the double hybrids on an open shell, since their
perturbative term is an MP2.

Unrestricted **Kohn-Sham** does work, over a spin-polarised functional
evaluation, and needs nothing said in the deck: a multiplicity above one or an
odd electron count selects it.

Functionals
^^^^^^^^^^^

Named to `libxc <https://libxc.gitlab.io/>`_, so a functional is available by its
libxc name -- ``gga_x_pbe``, ``hyb_gga_xc_b3lyp``, ``mgga_c_tpss`` and most of the
several thousand others. Nothing here shadows that list.

**Range-separated hybrids are supported** -- CAM-B3LYP and the ωB97 family split
their exchange into short and long range over an erf-attenuated kernel, and libcint
computes those integrals from the same entry points with a range parameter set. So
the long-range part is a second exchange pass rather than new integral code, and
the two are combined on the coefficients libxc reports for the functional itself.
One consequence is worth knowing: **a range-separated functional requires the
direct Fock build**, which is the default. Asking for one together with in-core or
density-fitted integrals is refused rather than run, because those tensors are
built for the full Coulomb kernel and the long-range exchange would simply be
absent.

**Functionals carrying non-local correlation** -- VV10, so ωB97X-V, ωB97M-V and
B97M-V -- are evaluated rather than refused. That term is a double integral over
the density rather than a functional of it at a point, so libxc supplies only the
semilocal half and reports the two parameters ``b`` and ``c`` through
``xc_nlc_coef``; the integral itself is this program's work. See
`Non-local correlation (VV10)`_ below for the grid it runs on.

**One family is still refused rather than approximated**, checked on libxc's own
report of the functional rather than a list of names kept here:

- **meta-GGAs needing the density Laplacian**, which is a second derivative of every
  basis function. Detected by the ``XC_FLAGS_NEEDS_LAPLACIAN`` flag.

The refusal names what is missing. Why it is a refusal rather than an
approximation is worth stating: a functional whose exchange coefficient is taken
at face value when it does not mean what a global hybrid's means returns a
converged energy several Hartree out -- 3.4 for CAM-B3LYP, 6.4 for ωB97X, before
range separation was handled -- and nothing about either run looks wrong. Dropping
VV10 is the same kind of error at a smaller size: on water/STO-3G it moves ωB97X-V
by 43 mHa, 27 kcal/mol, on three atoms, and the run converges and prints a
plausible number either way.

What metalquicha adds on top is friendly names for combinations libxc has the
parts for but not a name for the pair, plus double hybrids, which libxc does not
carry at all:

.. list-table::
   :header-rows: 1

   * - Name
     - Rung
     - Composition
   * - ``svwn``, ``lda``, ``lsda``
     - LDA
     - ``lda_x`` + ``lda_c_vwn``
   * - ``pbe``
     - GGA
     - ``gga_x_pbe`` + ``gga_c_pbe``
   * - ``blyp``
     - GGA
     - ``gga_x_b88`` + ``gga_c_lyp``
   * - ``b3lyp``
     - hybrid
     - ``hyb_gga_xc_b3lyp`` (libxc reports its own exchange fraction)
   * - ``pbe0``, ``pbeh``
     - hybrid
     - ``hyb_gga_xc_pbeh``
   * - ``tpss``
     - meta-GGA
     - ``mgga_x_tpss`` + ``mgga_c_tpss``
   * - ``m06-l``, ``m06l``
     - meta-GGA
     - ``mgga_x_m06_l`` + ``mgga_c_m06_l``
   * - ``wb97x``
     - range-separated hybrid
     - ``hyb_gga_xc_wb97x`` (ω = 0.3)
   * - ``cam-b3lyp``, ``camb3lyp``
     - range-separated hybrid
     - ``hyb_gga_xc_cam_b3lyp`` (ω = 0.33)
   * - ``wb97x-v``, ``wb97xv``
     - range-separated hybrid + VV10
     - ``hyb_gga_xc_wb97x_v``
   * - ``wb97m-v``, ``wb97mv``
     - range-separated meta-GGA hybrid + VV10
     - ``hyb_mgga_xc_wb97m_v``
   * - ``b97m-v``, ``b97mv``
     - meta-GGA + VV10
     - ``mgga_xc_b97m_v``
   * - ``b2plyp``
     - double hybrid
     - 0.53 exact exchange, 0.47 B88; 0.27 MP2, 0.73 LYP
   * - ``b2gp-plyp``
     - double hybrid
     - 0.65 / 0.36
   * - ``mpw2plyp``
     - double hybrid
     - 0.55 / 0.25 over mPW91

Every one of these is validated against a reference: the semilocal and hybrid
ones against ``pyscf.dft.RKS`` at the same grid level, agreeing to 7.7e-11 or
better; the three double hybrids against ``pyscf.dh.DFDH``, whose reported
coefficients match the table above exactly. ``blyp`` is the exception in kind
rather than in rigour -- it has no standalone energy case in the generated CPU
suite, and is instead checked against PySCF through the exchange-correlation
kernel and potential-derivative programs, where it is the reference functional
precisely because it is B2PLYP's semilocal part at full weight.

Anything beyond meta-GGA is refused, as is a functional needing the density
Laplacian -- on libxc's own say-so rather than a guess about which ones are safe.

Continuum Solvation
^^^^^^^^^^^^^^^^^^^

A polarizable continuum, on both ab initio backends. See
:doc:`continuum_solvation` for what the model is, the worked example, and the
cavity conventions that decide whether a comparison against another code can
succeed. This section is the keyword reference.

.. code-block:: json

   "keywords": {
     "pcm": {
       "method": "cpcm",
       "dielectric": 78.3553,
       "angular_points": 302,
       "radii_scale": 1.2
     }
   }

Naming the block is what turns solvation on -- there is no separate flag, because
a deck that states a dielectric wants solvent and two switches could disagree.

- ``method``: which continuum model the surface charges solve. ``"cpcm"``
  (conductor-like, the default) or ``"iefpcm"`` (the integral equation
  formalism) on the CPU backend. They are different models with different
  energies -- on the hydroxide anion in water they differ by 5.1e-5 Hartree --
  so the choice is named rather than implied. The cuEST backend's solver is
  fixed and accepts only the default; asking it for ``"iefpcm"`` is refused
  rather than substituted.
- ``dielectric``: the solvent's dielectric constant. **Required.** There is no
  solvent-name table on this path: tblite has one for its own CPCM, and a second
  table here that disagreed with it would make the same word mean two things.
- ``angular_points``: Lebedev points per atom on the cavity surface. 302 by
  default, and not a knob to economise on: 110 was the default until it was
  measured, and it gave a dielectric energy 21% short of the converged one. The
  surface integrand has cusps where the atomic spheres intersect. On the CPU
  backend the order must also carry a fitted SWIG exponent -- 6, 14, 26, 38,
  50, 86, 110, 146, 170, 194, 302, 350, 434, 590 or denser; the odd orders
  74, 230 and 266 exist in the grid tables but are refused, because the
  exponent table [J. Chem. Phys. 122, 194110 (2005)] does not cover them.
- ``radii_scale``: multiplies the van der Waals radii in ``mqc_pcm_radii``, which
  are Bondi's filled in from Mantina's. 1.2 is the usual convention. **The radii
  are the cavity and the cavity is most of the answer** -- see
  :doc:`continuum_solvation` before comparing against another code's defaults.
- ``zeta``: **cuEST only** -- the Gaussian switching prefactor for its smooth
  cavity surface, used as :math:`\zeta_i = \zeta \sqrt{n_{ang}} / R_i`, an
  empirical constant of that interface. The CPU backend refuses a deck that
  sets it: its exponents are the fitted per-point SWIG values and have no free
  prefactor, so the keyword would be silently ignored -- which is exactly what
  this program does not do.
- ``tolerance``, ``max_iter``: bounds on cuEST's iterative (conjugate gradient)
  surface-charge solve; a final solve that did not converge is refused rather
  than reported. The CPU backend factorizes the charge equations once and
  solves directly, which meets any tolerance these could ask for -- the keys
  are accepted and moot there.

This is **not** ``keywords.xtb.solvation_model: cpcm``. That configures tblite's
continuum, which builds its own cavity with its own defaults, and the two are
separate models rather than one keyword with two backends.

What each backend supports, refuses, and how each was validated is on
:doc:`continuum_solvation`.


DFT Options
^^^^^^^^^^^

The integration grid, and how the quadrature walks it:

.. code-block:: json

   "dft": {
     "grid_level": 3,
     "nlc_grid_level": 1,
     "screening_tolerance": 1e-12,
     "block_size": 512
   }

- ``grid_level``: 0 to 9, from the standard per-element tables -- the same tables
  PySCF uses, which is what makes a level-for-level comparison meaningful. Default
  3, and where a production calculation should start.
- ``nlc_grid_level``: the level of the *separate* grid VV10's double integral runs
  on, for a ``-V`` functional. Default 1. Negative keeps the default. Ignored
  entirely by a functional with no non-local term. See
  `Non-local correlation (VV10)`_.
- ``radial_points`` / ``angular_points``: override the level for every atom, which
  is what a convergence study wants. Supplying these takes the level out of
  charge; supplying one without the other is refused rather than half-applied.
- ``screening_tolerance``: the value below which a basis function is treated as
  absent from a block of grid points, so that the block's work skips it. Default
  ``1e-12``.

  **Zero or negative evaluates the whole basis.** That is a supported setting
  rather than an error: it is how a run measures what the screening is worth, and
  how a disagreement gets bisected without rebuilding. The two paths should agree
  to a few times the threshold -- on glycine tripeptide in 6-31G* they give
  ``-700.413700086596`` and ``-700.413700086592``.
- ``block_size``: grid points per block. Default 512; a non-positive value keeps
  it.

  This is two decisions wearing one name, which is worth knowing before turning
  it. A smaller block covers less space, so fewer shells reach it and the screen
  keeps less of the basis -- on a 638-function peptide, 4096-point blocks keep
  141 functions and 256-point blocks keep 116. But the loop over blocks *is* the
  OpenMP loop, so the block count is also the thread granularity: at 4096 points
  a level-1 grid is only 72 blocks, which on a 104-thread node leaves a third of
  the threads with nothing to do at all. Too small and the per-block bookkeeping
  starts to show instead. The useful value depends on the grid and the machine
  together, which is why it is a keyword rather than a constant.

.. note::

   Both apply to every grid loop the exchange-correlation code has: the
   Kohn-Sham potential built during the SCF, restricted and unrestricted alike;
   the kernel behind CPHF, the response properties and the Z-vector a double
   hybrid or an RI-MP2 gradient solves; and the exchange-correlation gradient,
   again for either reference and including the reference-operator term a double
   hybrid needs. A ``Gradient`` or ``Optimize`` run therefore sees the benefit in
   the derivative as well as in the SCF iterations before it. The SCF-side paths
   are threaded over blocks of grid points as well.

   One part of the gradient is deliberately outside this. The derivatives of the
   Becke partition weights depend on the nuclear positions and the grid, not on
   the basis, so there is nothing for a basis-function screen to drop and these
   keywords do not reach them.

Non-local correlation (VV10)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``-V`` functional -- ``wb97x-v``, ``wb97m-v``, ``b97m-v`` -- carries a
correlation term that is not a functional of the density at a point:

.. math::

   E_c^{nl} = \int \rho(\mathbf{r})
              \left[ \tfrac{1}{2} \int \rho(\mathbf{r}')\,
              \Phi(\mathbf{r},\mathbf{r}')\, d\mathbf{r}' + \beta \right]
              d\mathbf{r}

libxc evaluates the semilocal half and reports the functional's own ``b`` and
``c``; the double integral is evaluated here, folded into the potential, and the
SCF converges on the resulting density. It is not a correction applied to a
converged result afterwards.

Both spins are handled. VV10 depends on the total density only, with no spin
dependence anywhere in the kernel, so an unrestricted calculation evaluates it
once on :math:`\rho_\alpha + \rho_\beta` and adds the identical contribution
to each spin's matrix.

The **nuclear gradient** is available for a closed shell, so a ``-V`` functional
can be optimized. It needs no new integrals: the pair sum is spent producing
:math:`\delta E/\delta\rho` and :math:`\delta E/\delta\sigma`, the same two
quantities the Fock build consumes, after which the contraction is the one a GGA
uses. Two terms have no semilocal counterpart and are included -- the kernel
depends on *where* the grid points are, through :math:`|\mathbf{r}-\mathbf{r}'|`,
and a point's weight enters the energy twice rather than once. Leaving either
out costs 4.5e-4 on water and does not shrink when the grid is refined, since
neither is a quadrature error. The open shell is refused rather than
approximated. The second derivative is implemented as well, for a restricted
reference; see :doc:`analytic_hessians`.

Why it needs its own grid
"""""""""""""""""""""""""

The double integral costs the *product* of two point counts, per SCF iteration,
where everything else on the grid costs one. At the default ``grid_level`` of 3
that is 33,704 points on water and 84,712 on benzene -- 3.6 billion pairs an
iteration for benzene, which is not a calculation anyone would wait for. So the
non-local term runs on a grid of its own, and ``nlc_grid_level`` sets it.

Accuracy of the VV10 energy against PySCF's, water/STO-3G, by level:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - ``nlc_grid_level``
     - Points (water)
     - Error in the VV10 energy
   * - 0
     - 2,328
     - 2.95e-05 Ha
   * - 1 (default)
     - 10,128
     - 4.44e-07 Ha
   * - 2
     - 21,952
     - 2.09e-09 Ha
   * - 3
     - 33,704
     - 8.00e-12 Ha

Level 1 is the default because 4.4e-07 Ha is below what the rest of the
calculation is uncertain by, at an eleventh of level 3's pair count. Raise it if
you are comparing against a published number computed on a fine non-local grid;
the term converges quickly, so there is rarely a reason to go past 2.

.. note::

   The nuclear gradient (closed shell) and the analytic Hessian (restricted)
   both carry the non-local term -- a ``-V`` functional can be optimized and
   its frequencies computed. See :doc:`analytic_hessians` for the second
   derivative and its NLC-grid caveat.

.. note::

   The GPU path calls cuEST's own ``cuestNonlocalXCPotential{RKS,UKS}Compute``
   for this term, rather than reimplementing it. Whether the functional carries
   VV10 is asked of the XC plan (``CUEST_XCINTPLAN_IS_VV10``) rather than
   matched against a list of names, so ``nlc_grid_level`` does not apply there:
   cuEST chooses its own quadrature for the non-local term.

   .. warning::

      The GPU wiring compiles and is type-checked against cuEST's headers, but
      **has not been run against the library**, which needs an sm_80 card. Until
      it has, check a ``-V`` energy on the GPU against the CPU path before
      trusting it. The GPU gradient is refused, as on CPU.

Driver Section
--------------

Specifies the calculation type:

.. code-block:: json

   "driver": "Energy"

**Supported drivers**:

- ``Energy``: Single-point energy calculation
- ``Gradient``: Energy + analytical gradient (if method supports it)
- ``Hessian``: Second derivatives by finite difference, plus vibrational analysis
- ``Optimize``: Minimize the geometry, calling the method for a gradient at each
  step. ``Optimization`` and ``Opt`` are accepted for the same thing. Needs a build
  with the optimizer, and a method that has gradients. See
  :doc:`geometry_optimization`.
- ``MakeFP``: Build an effective fragment potential and write it as a ``.efp``
  file, computing no energy. ``MakeEFP`` is accepted for the same thing. See
  :doc:`makefp`.

Keywords Section
----------------

SCF Options
^^^^^^^^^^^

.. code-block:: json

   "keywords": {
     "scf": {
       "maxiter": 300,
       "tolerance": 1e-6,
       "guess": "auto",
       "unrestricted": false,
       "density_fitting": false,
       "level_shift": 0.0,
       "allow_crap_scf": false,
       "linear_dependence_threshold": 1e-7
     }
   }

- ``maxiter``: Maximum SCF iterations (default: 300)
- ``tolerance``: Convergence tolerance (default: 1e-6)
- ``guess``: Initial orbital guess (default: ``auto``). One of:

  - ``core`` -- the core Hamiltonian, ``F = H``. Cheapest and worst; on a
    46-atom peptide cation it does not converge at all in 100 cycles.
  - ``gwh`` -- generalized Wolfsberg-Helmholz, needing no atomic calculation.
  - ``sad`` -- superposition of spherically averaged atomic densities. Each
    element is solved once as a free atom and cached.
  - ``sac`` -- superposition of the free atoms' own spin densities. Arrives with
    its spatial symmetry already broken, which is what a radical wants and a
    closed shell does not.
  - ``auto`` -- let the backend choose, since the best starting point is a
    property of the backend: the CPU path resolves it to ``sad`` and the GPU path
    to ``gwh``, each having measured its own.

- ``unrestricted``: Force the unrestricted path (default: false). An odd electron
  count or a multiplicity above one forces it regardless.
- ``density_fitting``: Fit J and K in the reference (default: false). Asked for
  explicitly rather than inferred from ``aux_basis`` being present.
- ``level_shift``: Hartree added to the virtual orbitals before each
  diagonalisation, to widen the gap the next density is built through (default:
  0.0, meaning off). Use it on an SCF that oscillates rather than one that
  crawls; 0.2 to 1.0 is the useful range. A negative value is refused, since it
  narrows the gap and makes convergence worse.

  The shift is tapered off before convergence, so the orbital energies reported
  at exit belong to the unshifted operator and stay usable as correlation
  denominators. See :ref:`level-shifting` for when it will not help, which is a
  shorter list than it sounds.
- ``allow_crap_scf``: Accept a non-converged SCF instead of stopping (default:
  false). Off because the energy of an SCF that ran out of iterations has the
  right magnitude and nothing downstream can tell.
- ``linear_dependence_threshold``: Overlap eigenvalues at or below this are
  dropped when the orthogonaliser is built (default: 1e-7).

  Every SCF here runs in the canonical orthogonal basis :math:`X = U s^{-1/2}`,
  and a basis with diffuse functions on a compact or large system produces
  eigenvalues of the overlap near zero. Those modes are combinations the basis
  cannot really distinguish; dividing by their square root would put noise into
  every iteration, so they are discarded and the SCF runs in fewer orbitals than
  the basis set defines.

  **A run says which of three things happened, and the middle one is why this is
  worth setting.** Functions dropped: a warning naming the count, the eigenvalue
  and the orbital count everything downstream will use -- the basis is no longer
  the one named in the input, so energies do not compare with a run that kept
  them all. Nothing dropped but the smallest eigenvalue below 1e-5: also a
  warning, because :math:`X` then carries a large :math:`1/\sqrt{s}` into every
  iteration and the SCF can converge to different solutions from different
  guesses without ever failing. Otherwise a single line, at verbose level only.

  Raise it to shed more of a diffuse basis when the second case appears; lower it
  to keep modes the default discards. Both warnings are printed regardless of
  logger level, since a run that silently solves a smaller problem than the one
  asked for is the thing this exists to prevent.

  In a fragmented run this is often the only way to finish at all -- a few
  fragments out of millions will not converge, and stopping on the first one
  wastes the other million. What makes that safe rather than merely tolerable
  is that the fragments which failed are named in the output, with the monomers
  each was built from, so the run can be followed up rather than trusted. See
  :ref:`unconverged-fragments`.

Correlation Options
^^^^^^^^^^^^^^^^^^^

Shared by every post-Hartree-Fock method, and deliberately not under ``scf``: a
density-fitted reference followed by a conventional correlation treatment is a
combination someone will ask for, and one shared flag could not express it.

.. code-block:: json

   "correlation": {
     "freeze_core": true,
     "n_frozen_core": -1,
     "density_fitting": false,
     "scs": false
   }

- ``freeze_core``: Leave core orbitals uncorrelated (default: true)
- ``n_frozen_core``: How many to freeze (default: -1, counted from the elements)
- ``density_fitting``: Fit the correlation integrals (default: false, or true if
  the method name carries an ``RI-``/``DF-`` prefix -- an explicit keyword wins)
- ``scs``, ``scs_ss``, ``scs_os``: Spin-component scaling and its factors
  (defaults 1/3 and 1.2, applied only when ``scs`` is on or the method name asks)

Coupled Cluster Options
^^^^^^^^^^^^^^^^^^^^^^^

What only an iterative correlation method needs.

.. code-block:: json

   "cc": {
     "maxiter": 100,
     "tolerance": 1e-8,
     "diis": true,
     "diis_size": 8
   }

- ``maxiter``: Maximum amplitude iterations (default: 100)
- ``tolerance``: Correlation energy convergence (default: 1e-8)
- ``diis`` / ``diis_size``: Extrapolate the amplitudes (default: true, 8)
- ``triples``: Override whether (T) runs. Ordinarily the method name settles it,
  since ``ccsd`` and ``ccsd(t)`` are separate methods rather than one method with
  a flag; set this only to contradict the name.
- ``spin_adapted``: Which closed-shell formulation runs (default: true). Both are
  exact for a closed shell and agree to machine precision, so this chooses how a
  number is computed and not which number; spatial orbitals are roughly sixteen
  times smaller and several times faster.

**Open-shell systems** need no keyword. A multiplicity other than 1, an odd
electron count, or ``keywords.scf.unrestricted`` puts the reference in an
unrestricted SCF, and coupled cluster follows it there -- ``spin_adapted`` is
ignored, since those equations are derived for a closed shell and have no beta
orbitals to be given. Density fitting is the one combination refused: the fitted
three-index block has no spin blocks, so ``ri-ccsd`` over an open shell errors
rather than returning the alpha-only answer.

Active Space Options
^^^^^^^^^^^^^^^^^^^^

For ``casscf``, ``casci`` and ``mcscf``. A complete active space calculation
starts from a closed-shell SCF, partitions its orbitals into doubly occupied
*inactive*, an *active* set the CI distributes electrons over in every possible
way, and empty *virtual* ones, and then -- for CASSCF -- optimises the orbitals
alongside the CI coefficients.

.. code-block:: json

   "mcscf": {
     "n_active_electrons": 6,
     "n_active_orbitals": 6,
     "n_inactive_orbitals": -1,
     "optimize_orbitals": true,
     "max_macro_iter": 100,
     "orbital_convergence": 1e-6
   }

- ``n_active_electrons`` / ``n_active_orbitals``: The active space, written
  CAS(e,o). **Required unless** ``avas`` chooses it -- there is no default,
  because the right active space is a property of the chemistry rather than of
  the molecule, and a guess would produce a converged energy for a calculation
  nobody asked for.
- ``n_inactive_orbitals``: Doubly occupied orbitals below the active space
  (default: -1, derived as ``(nelec - n_active_electrons) / 2``). Set it only to
  ask for a partition other than the obvious one; the electrons still have to
  add up.
- ``optimize_orbitals``: Move the orbitals as well as the CI coefficients
  (default: true, or false if the method is spelled ``casci``). All three
  spellings are one method type, so this is the only thing that distinguishes
  CASSCF from CASCI; set it only to contradict the name.
- ``max_macro_iter``: Orbital optimisation iterations (default: 100)
- ``orbital_convergence``: Largest orbital gradient element accepted as
  converged (default: 1e-6)

The spin split is not a keyword. ``molecular_multiplicity`` settles it: every
inactive orbital is doubly occupied and contributes nothing to Ms, so the whole
of the excess alpha population sits in the active space and
``n_alpha - n_beta = multiplicity - 1`` exactly. An open-shell *state* on an
even-electron molecule is reachable this way; an odd electron count is not,
since the reference SCF is restricted.

Choosing the Active Space Automatically
"""""""""""""""""""""""""""""""""""""""

Deciding which orbitals belong in the active space is the hardest part of using
CASSCF, and getting it wrong gives a converged, plausible, wrong answer. AVAS
does it from a description of the chemistry instead: name the *atomic* orbitals
the interesting physics lives in, and the projection works out which molecular
orbitals carry that character.

.. code-block:: json

   "mcscf": {
     "avas": {"orbitals": ["N 2s", "N 2p"], "threshold": 0.2},
     "max_macro_iter": 300
   }

- ``orbitals``: Atomic orbital labels, each an element symbol, a space, a
  principal quantum number and a subshell letter -- ``"N 2p"``, ``"Cr 3d"``.
  Required inside the block. Only shells a *free atom* occupies exist to be
  asked for, so there is no ``"N 3d"``.
- ``threshold``: How much of the requested atomic character a molecular orbital
  needs to join the active space, between 0 and 1 (default: 0.2).

Every molecular orbital is scored by how much of it is the atomic orbitals
named, and everything above the threshold becomes active. For N\ :sub:`2`
asking for ``"N 2p"``, the occupied orbitals score 0.000, 0.000, 0.000, 0.779,
0.959, 0.991, 0.991 -- three with no nitrogen 2p character at all and four made
almost entirely of it, with nothing in between. Four occupied plus the matching
antibonding orbitals gives CAS(8,7): the triple bond and its antibonds, which is
what a careful person would have chosen by hand.

That gap is why ``threshold`` rarely needs touching. Anywhere in it gives the
same answer. A request that produces a *continuum* of scores instead is telling
you the question was badly posed, not that the threshold needs tuning.

Naming an active space twice -- an ``avas`` block *and* ``n_active_electrons``
-- is refused rather than resolved by precedence, since the counts would be
silently discarded in favour of whatever the projection decided.

Reference: Sayfutyarova, Sun, Chan and Knizia, *J. Chem. Theory Comput.* **13**,
4063 (2017). The projection here uses this code's own free-atom minimal basis
rather than the MINAO set of the paper; the two select the same spaces, which is
what the width of that gap predicts.

Taking the Whole Valence Shell
""""""""""""""""""""""""""""""

The third way of naming an active space, and the one that asks nothing of the
reader:

.. code-block:: json

   "mcscf": {"full_valence": true}

Every occupied orbital that is not core, plus the valence-virtual orbitals that
complete the free-atom minimal basis. There is no threshold and no list of
labels -- the size is ``n_mbs - n_occupied`` from counting the elements, so the
same molecule gives the same space in any basis set. N\ :sub:`2` comes out
CAS(10,8): two 1s cores inactive, five valence occupied, three valence virtual.

The valence-virtual orbitals are the same ones :doc:`bonding_analysis` extracts;
there they are a basis to analyse in, here they are the empty half of an active
space. The run prints the diagnostic that says whether the split is clean --
for N\ :sub:`2`, a smallest-kept singular value of 0.999 against a
largest-rejected 0.201. A narrow gap there means the minimal basis is not
finding a valence space and only the counting is holding the answer together.

It grows with the molecule rather than with what is interesting in it, so a
complete expansion over the valence shell is out of reach past small systems.
That is what ``ormas`` above is for, and the two are meant to be used together.

Naming the space twice -- ``full_valence`` with counts, or with an ``avas``
block -- is refused rather than resolved by precedence, for the reason given
there.

Restricting the Occupations
""""""""""""""""""""""""""""

A complete active space distributes its electrons over its orbitals in every
way there is, and the determinant count grows factorially: CAS(14,14) is 11.8
million and CAS(16,16) is 166 million. Most of that is excitations nobody
believes in. ORMAS -- occupation restricted multiple active space -- cuts the
active orbitals into consecutive subspaces and puts a window on how many
electrons each may hold, which keeps the part of the expansion that matters and
discards the rest.

.. code-block:: json

   "mcscf": {
     "n_active_electrons": 6,
     "n_active_orbitals": 6,
     "optimize_orbitals": false,
     "ormas": {
       "subspaces": [1, 4],
       "min_electrons": [4, 0],
       "max_electrons": [6, 2]
     }
   }

- ``subspaces``: The **active** orbital each subspace starts at, ascending. The
  first entry is always 1, and the subspaces run to the end of the active space,
  so ``[1, 4]`` on six active orbitals means 1--3 and 4--6. These are positions
  within the active space, not molecular orbital numbers: an inactive orbital is
  not part of the partition.
- ``min_electrons`` / ``max_electrons``: The fewest and most electrons a
  subspace may hold, counting **both spins together**. This is what lets the
  restriction be independent of how the spins are arranged.

All three lists have the same length, and a deck where they do not is refused
with the key names rather than an error from further down.

The example is singles and doubles: at least four of the six electrons stay in
the lower three orbitals, so at most two are promoted. What that costs is worth
seeing -- for N\ :sub:`2` in cc-pVDZ it is 118 determinants where the complete
space is 400, and the energy rises accordingly, because a restricted space is
variationally above the space that contains it.

Some shapes worth knowing:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Wanted
     - How to write it
   * - Truncated CI (singles and doubles)
     - Two subspaces, ``max_electrons`` of 2 on the upper one
   * - RAS1/RAS2/RAS3
     - Three subspaces; ``min_electrons[0] = 2*n1 - holes`` and
       ``max_electrons[2] = particles``
   * - Two non-communicating active spaces
     - ``min_electrons`` equal to ``max_electrons`` everywhere, so no charge
       moves between them
   * - Fragments with limited charge transfer
     - One subspace per fragment, windows one either side of neutral

The windows are tightened before use. A subspace cannot hold more than the rest
of the partition can spare, so a maximum that promises more is reduced, and the
numbers reported back may not be the numbers written. This removes promises that
could not have been kept and never changes which determinants exist.

**Orbital optimisation is not implemented for a restricted space** and is
refused rather than approximated -- set ``optimize_orbitals`` to false and the
CI runs on the reference orbitals. The reason is not effort: rotating one active
orbital into another stops being redundant the moment the space is not complete,
so the orbital gradient acquires a block a CASSCF has no term for, and running
one anyway converges quietly to something that is not the answer.

Reference: Ivanic, *J. Chem. Phys.* **119**, 9364 (2003). The implementation
here is checked against GAMESS, which has had ORMAS since that paper: for water
in 3-21G with one frozen core and singles and doubles out of the valence, both
give -75.7103507602.

State averaging and CASPT2/NEVPT2 are not implemented, and no keyword accepts
them -- a deck asking for either is refused rather than quietly given a
ground-state energy. Derivatives are refused for the same reason: there is no
CASSCF gradient here, analytic or numerical.

Fragmentation Options
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

   "fragmentation": {
     "method": "MBE",
     "allow_overlapping_fragments": false,
     "level": 2,
     "max_intersection_level": 3,
     "expansion": "mbe",
     "counterpoise": "none",
     "embedding": "none",
     "cutoff_method": "distance",
     "distance_metric": "min",
     "cutoffs": {
       "dimer": 10.0,
       "trimer": 8.0
     }
   }

**Parameters**:

- ``method``: ``"MBE"`` (Many-Body Expansion) or ``"GMBE"`` (Generalized MBE for overlapping fragments)
- ``allow_overlapping_fragments``: ``true`` for GMBE, ``false`` for standard MBE (default: ``false``)
- ``level``: Maximum fragment size (1=monomers only, 2=up to dimers, 3=up to trimers, etc.)
- ``max_intersection_level``: For GMBE only - maximum k-way intersection depth (default: level + 1)
- ``expansion``: ``"mbe"`` (default), ``"fmo"`` or ``"ee-mbe"`` - see :doc:`fmo`
- ``counterpoise``: ``"none"`` (default) or ``"vmfc"`` for the basis-set
  superposition correction - see :doc:`counterpoise`
- ``embedding``: Fragment embedding scheme (currently only ``"none"`` supported)
- ``cutoff_method``: How to include fragments (``"distance"``, ``"all"``)
- ``distance_metric``: For distance cutoffs: ``"min"``, ``"max"``, ``"com"`` (center of mass)
- ``cutoffs``: Distance thresholds (in Angstroms) for including dimers, trimers, etc.

Properties Section
------------------

Analyses to run once the wave function exists, beside ``keywords`` rather than
inside it. The distinction is worth keeping: ``keywords`` say how to compute the
wave function and change the number that comes out, while ``properties`` ask for
something further to be done with one already determined and change no energy.
The driver stays ``"energy"``.

.. code-block:: json

   "properties": {
     "bonding_analysis": {
       "type": "gms_quao",
       "energy_threshold": 1.0
     }
   }

- ``bonding_analysis``: The quasi-atomic bonding picture -- which atoms are
  bonded, by sigma or pi, how strongly, and where the lone pairs are. See
  :doc:`bonding_analysis`.

- ``fukui``: Where the molecule reacts.

  .. code-block:: json

     "properties": {
       "fukui": {}
     }

  Condensed Fukui indices, by difference in the electron count: the molecule is
  run again with one electron added and once with one removed, at the same
  geometry and in the same basis, and the atomic charges are differenced.
  ``f+`` says where the molecule accepts charge, ``f-`` where it gives charge
  up, and the dual descriptor ``f+ - f-`` separates electrophilic sites from
  nucleophilic ones in a single column. The three energies also give the
  ionisation potential, electron affinity, chemical potential, hardness and
  electrophilicity index, which cost nothing extra once the ions have been run.

  Computed by difference rather than from the frontier orbitals. Approximating
  ``f+`` by the shape of the LUMO costs one SCF instead of three and discards
  the relaxation of every other orbital in response to the added charge -- which
  on a polar molecule is most of the answer.

  Works with ``hf`` and with ``dft``, **double hybrids included**. Under DFT
  the two ions are run as unrestricted Kohn-Sham with the same functional as
  the neutral, so the three energies are comparable and the ionisation
  potential and electron affinity mean what they say.

  For a double hybrid the perturbative term is evaluated for all three states,
  the neutral's recomputed here rather than taken from the energy above it.
  That is deliberate: IP and EA are differences of total energies, so the term
  has to be present in all three or in none, and sourcing one of them
  elsewhere would make the result depend on two code paths agreeing about
  frozen cores and auxiliary bases. It costs one restricted MP2 beside the two
  unrestricted ones the ions already need.

  Naming an ``aux_basis`` fits that correlation, exactly as it does for the
  energy's own perturbative term -- the two follow the same choice, so the IP
  printed in the report and the total energy printed above it are never
  computed different ways.

  ``population`` is optional and defaults to ``chelpg``. It chooses how the
  density difference is condensed onto atoms:

  - ``chelpg`` -- charges fitted to the molecule's own electrostatic potential.
    The default choice, and the better behaved one.
  - ``mulliken`` -- **can be poor, and is offered because it is what much of
    the literature used.** It divides the overlap between two atoms straight
    down the middle regardless of what the atoms are, so it is sensitive to the
    basis set in a way that has nothing to do with chemistry, and it produces
    negative indices more readily. On formaldehyde it spreads ``f+`` as
    0.27/0.28/0.23/0.23 over the four atoms, ranking nothing, where CHELPG puts
    0.59 on the carbonyl carbon.

  .. warning::

     **Negative indices are spurious. Take them with a grain of salt.**

     The exact Fukui function is a derivative of the density with respect to
     electron count, so ``f+`` and ``f-`` cannot be negative: no site gives up
     charge because the molecule gained an electron. A negative condensed value
     is an artefact of dividing a continuous density among atoms -- the two
     states are fitted independently, and they can disagree by more than the
     one electron that moved between them.

     They are printed rather than hidden, with a warning, because silently
     clamping them would hide the fact that the partitioning struggled. Read
     the *ranking* and not the number: a small negative index means the site is
     unreactive in that channel, not that it repels charge. Do not quote the
     value, and do not build anything further on it -- ``f0`` and the dual
     descriptor inherit the artefact from whichever index carried it.

     A larger basis usually shrinks it. If it survives that, the atom simply
     carries very little of that channel. If you saw it under ``mulliken``,
     rerun with ``chelpg`` before concluding anything. The CHELPG fitting grid
     is not exposed as a deck setting, so it is not a knob available here.

  **Use a finer integration grid than you would for the energy.** ``grid_level``
  defaults to 3, which is right for a total energy but marginal here: the Fukui
  descriptors are *differences* between three separately converged states, so
  the quadrature error does not cancel the way it does within one SCF. On
  formaldehyde in 6-31G, ``m06-l`` at level 3 gives an ionisation potential
  2.7e-5 hartree away from the converged value; level 4 brings that to 4e-6 and
  level 5 reaches it, with level 6 no different. **Level 4 or 5 is the
  recommendation**, and it costs three grids rather than one because all three
  states are integrated.

  How much this matters depends on the functional rather than on its cost. The
  LDA, GGA and hybrid functionals are already within about 1e-6 at level 3;
  what moves are the meta-GGAs, which sample the kinetic-energy density and are
  steeper on the grid -- ``m06-l`` is the most sensitive of the set -- and
  ``wb97x``, whose range separation adds its own grid dependence.

  Two limits. Only a closed-shell neutral is accepted, so that both ions are
  doublets; anything else is refused rather than having its multiplicities
  guessed. And a warning appears when the anion comes out *above* the neutral,
  which means nothing bound the added electron -- either the basis has no
  diffuse functions or the molecule has no bound anion, and the two are
  indistinguishable from here. Water is the second kind: adding diffuse
  functions moves its affinity from -0.19 to -0.04 hartree and never makes it
  positive.

  Unfragmented calculations only. Which fragment of a many-body expansion
  receives the electron is not a well-posed question.

System Section
--------------

Logger Configuration
^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

   "system": {
     "logger": {
       "level": "Verbose"
     }
   }

**Supported log levels** (in order of verbosity):

- ``debug``: Most verbose, includes debug information
- ``verbose``: Detailed output
- ``info``: Standard output (default)
- ``performance``: Performance timing only
- ``warning``: Warnings only
- ``error``: Errors only
- ``knowledge``: Special knowledge-level output

Running Calculations
====================

Basic Usage
-----------

.. code-block:: bash

   # Run calculation (serial)
   ./mqc input.json

   # Run calculation (parallel with MPI)
   mpirun -np 4 ./mqc input.json

   # Run on multiple nodes
   mpirun -np 64 --map-by ppr:32:node ./mqc input.json

Output Files
------------

Metalquicha generates:

- **Console output**: Human-readable calculation summary
- **output_<basename>.json**: Machine-readable results with fragment energies, PIE coefficients, and total energy

Examples
========

Unfragmented Calculation
-------------------------

JSON input (``h3o.json``):

.. code-block:: json

   {
     "schema": {"name": "mqc-frag", "version": "1.0"},
     "molecules": [{
       "xyz": "h3o.xyz",
       "molecular_charge": 1,
       "molecular_multiplicity": 1
     }],
     "model": {"method": "XTB-GFN1"},
     "driver": "Energy"
   }

Run:

.. code-block:: bash

   ./mqc h3o.json

Fragmented MBE Calculation
---------------------------

JSON input (``prism.json``):

.. code-block:: json

   {
     "schema": {"name": "mqc-frag", "version": "1.0"},
     "molecules": [{
       "xyz": "prism.xyz",
       "molecular_charge": 0,
       "molecular_multiplicity": 1,
       "fragments": [
         [0,1,2], [3,4,5], [6,7,8],
         [9,10,11], [12,13,14], [15,16,17]
       ],
       "fragment_charges": [0, 0, 0, 0, 0, 0],
       "fragment_multiplicities": [1, 1, 1, 1, 1, 1]
     }],
     "model": {"method": "XTB-GFN1"},
     "keywords": {
       "fragmentation": {
         "method": "MBE",
         "level": 2
       }
     },
     "driver": "Energy"
   }

Run:

.. code-block:: bash

   ./mqc prism.json

GMBE with Overlapping Fragments
--------------------------------

JSON input (``overlapping_gly3.json``):

.. code-block:: json

   {
     "schema": {"name": "mqc-frag", "version": "1.0"},
     "molecules": [{
       "xyz": "gly3.xyz",
       "molecular_charge": 0,
       "molecular_multiplicity": 1,
       "fragments": [
         [0,1,2,3,4,5,6,7],
         [5,6,7,8,9,10,11,12,13],
         [10,11,12,13,14,15,16,17,18]
       ],
       "fragment_charges": [0, 0, 0],
       "fragment_multiplicities": [1, 1, 1]
     }],
     "model": {"method": "XTB-GFN1"},
     "keywords": {
       "fragmentation": {
         "method": "GMBE",
         "allow_overlapping_fragments": true,
         "level": 1,
         "max_intersection_level": 2
       }
     },
     "driver": "Energy"
   }

Note: Fragments 1-2 share atoms 5,6,7 and fragments 2-3 share atoms 10,11,12,13.

Run:

.. code-block:: bash

   ./mqc overlapping_gly3.json

Multi-Molecule Calculation
---------------------------

For calculating multiple conformers or isomers in one input:

.. code-block:: json

   {
     "schema": {"name": "mqc-frag", "version": "1.0"},
     "molecules": [
       {
         "xyz": "conformer1.xyz",
         "molecular_charge": 0,
         "molecular_multiplicity": 1
       },
       {
         "xyz": "conformer2.xyz",
         "molecular_charge": 0,
         "molecular_multiplicity": 1
       }
     ],
     "model": {"method": "XTB-GFN1"},
     "driver": "Energy"
   }

Each molecule is calculated independently, and results are organized by molecule index.

Best Practices
==============

1. **Keep your JSON decks in version control**: they are the input, not a
   generated artifact
2. **Start with small systems**: Test fragmentation schemes on small molecules first
3. **Check JSON output**: Verify fragment energies are reasonable
4. **Use appropriate log levels**: ``verbose`` for debugging, ``info`` for production
5. **Validate results**: Use the validation test suite (see :ref:`validation`)

Troubleshooting
===============

**"Invalid input file extension"**
   Ensure the file ends with ``.json``. Versions before 0.2.0 took ``.mqc``;
   see :ref:`migrating_from_mqc`.

**"Could not parse JSON input file"**
   The document is not valid JSON. Common causes are a trailing comma after
   the last entry of an object or array, and unquoted keys.

**"Missing required key"**
   ``schema.name``, ``schema.version``, ``model.method``, ``driver`` and
   ``molecules`` are all required, as are ``molecular_charge`` and
   ``molecular_multiplicity`` on each molecule.

**"Unknown key ..."**
   The deck is checked against the schema before it is read, so a key that is
   not part of the schema is an error rather than a setting that is silently
   ignored. The message names the key and lists what is allowed in its place;
   the usual cause is a misspelling or a key written at the wrong level.

**"fragment charges sum to N but the molecular charge is M"**
   ``fragment_charges`` must add up to ``molecular_charge``. This is checked
   because the alternative is a calculation that runs to completion on the
   wrong number of electrons.

**Fragment charge/multiplicity mismatch**
   Ensure sum of fragment charges equals molecular charge

**"No fragments generated"**
   Check that ``fragments`` array is not empty and indices are valid

**Hydrogen capping not working**
   Check that the bond is listed in ``connectivity`` and that its two atoms
   fall in different fragments -- that is what makes a bond broken, and it is
   derived rather than declared

.. _migrating_from_mqc:

Migrating from .mqc
===================

**Changed in 0.2.0.** ``.mqc`` and ``mqc_prep.py`` were removed; ``mqc`` reads
JSON directly.

If you drove metalquicha from JSON, as the documented workflow did, there is
nothing to change but the command: pass the ``.json`` file where you used to
pass the ``.mqc`` one, and drop the ``mqc_prep.py`` call. The schema is
unchanged.

If you hand-wrote ``.mqc`` files, they must be rewritten as JSON. The mapping
is direct -- each ``%section`` becomes an object of the same name:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - ``.mqc``
     - JSON
   * - ``%schema`` / ``name``, ``version``
     - ``"schema": {"name": ..., "version": ...}``
   * - ``%model`` / ``method``, ``basis``, ``aux_basis``, ``functional``,
       ``cartesian``
     - ``"model": {...}``, same keys
   * - ``%driver`` / ``type = Energy``
     - ``"driver": "Energy"``
   * - ``%structure`` + ``%geometry``
     - an entry in ``"molecules"``, with either ``"xyz"`` naming a file or
       ``"symbols"`` plus a flat ``"geometry"`` list
   * - ``%fragments`` / ``%fragment`` blocks
     - ``"fragments"``, ``"fragment_charges"``, ``"fragment_multiplicities"``
   * - ``%connectivity`` rows ``i j order``
     - ``"connectivity": [[i, j, order], ...]``
   * - ``%scf``, ``%hessian``, ``%aimd``, ``%fragmentation``, ``%xtb``
     - ``"keywords": {"scf": {...}, "hessian": {...}, ...}``
   * - ``%system`` / ``log_level``
     - ``"system": {"logger": {"level": ...}}``

Two differences worth knowing:

* **Broken bonds are no longer written down.** ``.mqc`` marked each bond
  ``broken`` or ``preserved``, because ``mqc_prep.py`` worked that out before
  emitting the file. ``mqc`` now derives it: a bond is broken when its two
  atoms do not belong to the same set of fragments. Just list the bonds.
* **JSON has no comments.** ``.mqc`` accepted ``#`` and ``!``. If you
  annotated your decks, the ``.xyz`` files they reference still take a
  comment on their second line.
