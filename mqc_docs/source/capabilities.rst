============
Capabilities
============

This document provides a comprehensive overview of Metalquicha's current features and capabilities.

Fragmentation Methods
=====================

Many-Body Expansion (MBE)
--------------------------

Standard many-body expansion for non-overlapping molecular fragments:

- **Supported levels**: 1-body (monomers) through N-body expansions
- **Energy expression**:

  .. math::

     E_{MBE(n)} = \sum_{I} E_I - \sum_{I<J} \Delta E_{IJ} + \sum_{I<J<K} \Delta E_{IJK} - \ldots

- **Fragment generation**: Automatic enumeration of all n-mer combinations
- **Counterpoise**: ``vmfc`` solves every subfragment in its parent's basis, so
  the basis-set superposition error cancels within each term instead of
  accumulating in it. Energies only; see :doc:`counterpoise`
- **Use cases**: Molecular clusters, water clusters, periodic systems

Generalized Many-Body Expansion (GMBE)
---------------------------------------

Advanced method for overlapping fragments using Principle of Inclusion-Exclusion (PIE):

- **Overlapping fragments**: Handles systems where fragments share atoms
- **PIE-based correction**: Automatically computes intersection terms with proper coefficients
- **GMBE(N) variants**:

  - GMBE(1): Primaries are monomers
  - GMBE(N): Primaries are N-mers (e.g., dimers for GMBE(2))

- **Intersection depth control**: Maximum k-way intersection level configurable
- **Use cases**: Strongly interacting systems, covalently bonded clusters

Fragment Molecular Orbital method (FMO_n)
------------------------------------------

Fragments solved in the electrostatic field of the others, iterated to
self-consistency, rather than in vacuum:

- **Two nested SCFs**: an inner fragment SCF against a fixed field, and an outer
  monomer loop that converges the field itself
- **Exact embedding**: nuclear attraction integrals plus the Coulomb operator of
  the neighbours' density matrices, not fitted point charges
- **Supported levels**: FMO2, FMO3 and beyond; level equal to the fragment count
  is the full expansion and reproduces the supermolecular energy exactly
- **Point-charge cutoff**: distant fragments approximated by atomic populations
  past a configurable separation (GAMESS's ``RESPPC``), which is what makes the
  method scale
- **Use cases**: hydrogen-bonded clusters, solvated molecules -- anywhere the
  monomers are substantially polarized
- **Restrictions**: non-covalent fragments only (enforced), closed shell,
  energies only

Electrostatically Embedded MBE (EE-MBE)
----------------------------------------

Many-body expansion with fragments embedded in point charges:

- **Point-charge embedding**: each neighbour represented by Mulliken or CHELPG
  atomic charges, self-consistently updated
- **Ordinary MBE assembly**: total embedded energies, no response term, which is
  what distinguishes it from FMO with the same charges
- **Shares all machinery with FMO**, including the level, the cutoff and the
  distributed execution
- **Use cases**: the same systems as FMO, at lower cost per fragment

See :doc:`fmo` for both.

Distance-Based Screening
-------------------------

Intelligent fragment filtering based on inter-monomer distances:

- **Cutoff specification**: Per n-mer level cutoffs (dimers, trimers, etc.)
- **Distance metric**: Minimal atom-to-atom distance between constituent monomers, for now limited up to octamers
- **Screening scope**:

  - **MBE**: Screens all generated fragments before calculation
  - **GMBE**: Screens primary fragments before PIE enumeration

- **Units**: Angstroms
- **Example**: With ``2 = 5.0``, all dimers with inter-monomer distance > 5.0 Å are excluded
- **Performance**: Dramatically reduces fragment count for large systems

Quantum Chemistry Methods
==========================

Extended Tight-Binding (XTB)
-----------------------------

Semi-empirical quantum chemistry via the `tblite <https://github.com/tblite/tblite>`_ library:

- **GFN2-xTB**: Latest parametrization, general-purpose, highest accuracy
- **GFN1-xTB**: Earlier version, faster, good for large systems

Ab Initio (CPU)
---------------

Gaussian-basis methods over `libcint <https://github.com/sunqm/libcint>`_, on the
CPU. This path exists to be checked against as much as to be run: it is the second
implementation the GPU backend is compared with, and everything below is validated
against PySCF on the same geometries and the same basis data.

- **Hartree-Fock**: restricted and unrestricted, with a direct, in-core or
  density-fitted Fock build. Initial guesses are core, GWH, and superposition of
  atomic densities or coefficients.
- **MP2**: conventional and density-fitted (RI-MP2), with spin-component scaling
  reported from separately-kept same- and opposite-spin components.
- **Coupled cluster**: CCSD and CCSD(T) over a restricted or an unrestricted
  reference. Two formulations for the closed-shell case -- spin-adapted over
  spatial orbitals by default, and spin orbitals for checking it against
  (``keywords.cc.spin_adapted: false``) -- which are exact for each other and
  agree to machine precision, so the choice is how a number is computed and not
  which number; spatial is the default because it is roughly sixteen times
  smaller. An open-shell reference takes the spin-orbital path necessarily,
  since the spin-adapted equations are derived for a closed shell. Density
  fitting is available for the restricted case only: the fitted three-index
  block has no spin blocks, so an unrestricted deck asking for it is refused
  rather than quietly given the restricted answer.
- **Kohn-Sham DFT**: the whole ladder -- LDA, GGA, hybrid, meta-GGA,
  range-separated hybrid and double hybrid -- restricted and unrestricted, over
  `libxc <https://libxc.gitlab.io/>`_, so most of what libxc carries is available
  by name, with friendly names and double-hybrid compositions layered on top.
  The quadrature drops basis functions that do not reach a block of grid points,
  which is what keeps its cost in proportion to the Fock build rather than
  dominating it; the threshold and the block size are both ``keywords.dft``.
  Range separation uses libcint's erf-attenuated integrals, which puts ωB97X and
  CAM-B3LYP in reach. Non-local correlation (VV10) is evaluated as the double
  integral it is, on a coarser grid of its own, which is what ωB97X-V, ωB97M-V and
  B97M-V need; it enters the Fock build, so the energy is self-consistent rather
  than a correction applied afterwards, and its nuclear gradient is refused.
  Laplacian-dependent meta-GGAs are still refused rather than approximated. Grids
  are Treutler-Ahlrichs
  radial times Lebedev angular with a Becke partition, from the same level tables
  PySCF uses.
- **Multiconfigurational**: CASSCF and CASCI over a complete active space, with
  the space named directly or chosen from atomic orbital character by AVAS. The
  active space can also be cut into subspaces with occupation windows (ORMAS),
  which expresses a truncated CI, a RAS, several non-communicating active spaces
  or a fragment model with limited charge transfer -- and keeps the determinant
  count down where a complete space would be hopeless. Orbital optimisation is
  complete-space only; a restricted space runs its CI on the reference orbitals.
  See :doc:`input_files`.
- **SAPT**: the interaction energy of two monomers as a sum of named physical
  terms -- electrostatics, exchange, induction, dispersion and their exchange
  counterparts -- rather than as a difference of two totals. ``sapt0`` is the
  first rung, monomers at Hartree-Fock and the interaction through second order;
  ``sapt2`` is the next one that exists, there being no SAPT1. Everything is in
  the dimer-centred basis, which is what makes the terms counterpoise-corrected
  by construction. Two fragments exactly -- see :doc:`sapt`.

Basis sets come from the Basis Set Exchange data shipped in ``basis_sets/``, and
whether a set is Cartesian or spherical is taken from the file rather than assumed.

Implicit Solvation
------------------

Implicit solvation models account for solvent effects without explicit solvent molecules:

**Supported models:**

- **PCM**: a polarizable continuum on **both** backends, for Hartree-Fock and
  Kohn-Sham, restricted and unrestricted, configured through ``keywords.pcm``.
  The CPU implementation is the smooth switching/Gaussian (SWIG) discretization
  of Lange and Herbert, solved as either C-PCM or IEF-PCM
  (``keywords.pcm.method``: ``cpcm`` or ``iefpcm``), and follows
  ``pyscf.solvent.pcm`` term for term -- the same Lebedev points per sphere, the
  same switching function, the same fitted per-point Gaussian exponents -- so the
  two can be compared directly. The cuEST (GPU) path is validated against
  PySCF's C-PCM.
- **ALPB**: Analytical Linearized Poisson-Boltzmann (recommended for GFN2-xTB)
- **GBSA**: Generalized Born with Solvent-Accessible Surface Area
- **CPCM**: Conductor-like Polarizable Continuum Model

**Supported solvents:**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Category
     - Solvents
   * - Water
     - ``water``, ``h2o``
   * - Alcohols
     - ``methanol``, ``ch3oh``, ``ethanol``, ``c2h5oh``, ``1-propanol``, ``propanol``,
       ``2-propanol``, ``isopropanol``, ``1-butanol``, ``butanol``, ``2-butanol``,
       ``1-octanol``, ``octanol``, ``decanol``
   * - Polar aprotic
     - ``acetone``, ``acetonitrile``, ``ch3cn``, ``dmso``, ``dimethylsulfoxide``,
       ``dmf``, ``dimethylformamide``, ``thf``, ``tetrahydrofuran``
   * - Aromatics
     - ``benzene``, ``toluene``, ``pyridine``, ``aniline``, ``nitrobenzene``,
       ``chlorobenzene``, ``phenol``
   * - Halogenated
     - ``chloroform``, ``chcl3``, ``dichloromethane``, ``ch2cl2``, ``dcm``,
       ``carbon tetrachloride``, ``ccl4``
   * - Ethers
     - ``diethylether``, ``ether``, ``dioxane``, ``furan``
   * - Alkanes
     - ``pentane``, ``hexane``, ``n-hexane``, ``cyclohexane``, ``heptane``,
       ``n-heptane``, ``octane``, ``n-octane``, ``decane``, ``hexadecane``
   * - Other organics
     - ``nitromethane``, ``formamide``, ``cs2``, ``carbondisulfide``
   * - Esters/acids
     - ``ethyl acetate``, ``ethylacetate``, ``acetic acid``, ``aceticacid``,
       ``formic acid``, ``formicacid``
   * - Special
     - ``woctanol`` (wet octanol), ``inf`` (infinite dielectric/conductor)

**Usage example:**

.. code-block:: json

   {
     "model": {
       "method": "XTB-GFN2",
       "solvent": "water",
       "solvation_model": "alpb"
     }
   }

Calculation Types
=================

Energy Calculations
-------------------

- **Single-point energies**: Total electronic energy
- **Component breakdown**: SCF, MP2 same/opposite spin, coupled-cluster
  singles/doubles/triples, kept separately so a total can be taken apart
- **Output**: Hartrees, per-fragment energies, MBE contributions

Gradient Calculations
---------------------

- **Analytical gradients**: Energy derivatives with respect to nuclear positions
- **Hydrogen cap redistribution**: Automatic gradient mapping for capped bonds
- **Units**: Hartree/Bohr
- **Applications**: Geometry optimization, molecular dynamics

What has a gradient on the CPU backend, at a glance:

.. list-table::
   :header-rows: 1
   :widths: 34 12 54

   * - Method
     - Gradient
     - Notes
   * - Hartree-Fock, restricted and unrestricted
     - yes
     - Exact or density fitted
   * - Kohn-Sham: LDA, GGA, hybrid, meta-GGA
     - yes
     - Restricted and unrestricted; meta-GGA restricted only
   * - Kohn-Sham: non-local correlation (``-V``)
     - yes
     - ωB97X-V, ωB97M-V, B97M-V. Self-consistent, restricted and unrestricted,
       on a separate grid set by ``keywords.dft.nlc_grid_level``. The gradient
       is restricted only, like the meta-GGAs, and CPU backend only: the GPU
       path calls cuEST's own non-local entry points and is compiled but not
       yet run against the library
   * - Effective core potentials
     - yes
     - ``model.ecp``, a separate file from the basis. Energies only: every
       nuclear derivative is refused, as are MCSCF, xTB and the GPU backend,
       and an automatically counted frozen core. Needs libfint, which is the
       default -- libcint has no ECP code, and a libcint build refuses
       ``model.ecp`` rather than ignoring it
   * - Kohn-Sham: range-separated hybrid
     - yes
     - **Needs an auxiliary basis.** The exact-ERI path builds no second
       exchange derivative at the screened omega, so the fitted one is the
       only one
   * - MP2, restricted
     - yes
     - All four combinations of fitted reference and fitted correlation,
       with or without a frozen core -- ``freeze_core`` is on by default, and
       the default deck gets the gradient of the energy it computed. The
       amplitudes span the active occupied space; the relaxed density gains an
       occupied-frozen block resolved directly from the Lagrangian. (No
       virtual-frozen block exists to build: ``n_frozen_core`` freezes leading
       core orbitals and no virtuals)
   * - Double hybrids (``b2plyp``, ``b2gp-plyp``, ``mpw2plyp``)
     - yes
     - Restricted, all-electron, GGA-based; exact or fitted reference
   * - CASSCF and ORMAS-SCF
     - yes
     - No Z-vector: a converged MCSCF is stationary with respect to both the
       orbital rotations and the CI coefficients, so the response terms vanish
   * - Coupled cluster
     - no
     - Needs the Lambda amplitudes
   * - Anything unrestricted beyond HF and semilocal/hybrid KS
     - no
     - See the refusal table below

**The four MP2 combinations are four different energies**, not one energy
reached four ways, and each gets its own gradient:

.. list-table::
   :header-rows: 1
   :widths: 22 22 56

   * - Reference
     - Correlation
     - What differentiates it
   * - exact
     - exact
     - ``mp2``. Four-centre derivatives throughout
   * - exact
     - fitted
     - ``ri-mp2``. The two-particle density stays three-index and contracts
       against three- and two-centre derivatives
   * - fitted
     - fitted
     - ``ri-mp2`` plus ``keywords.scf.density_fitting``. The response operator,
       both potentials built from the relaxed density and the reference's
       two-electron derivative all move onto the auxiliary basis
   * - fitted
     - exact
     - ``mp2`` plus ``keywords.scf.density_fitting``, and what a double hybrid
       needs. **Consistent, not cheap**: the two-particle density is still
       four-index and the four-centre integrals are still built. What fitting
       buys is that the reference being differentiated is the one the SCF
       converged

One auxiliary basis serves both halves when both are fitted, which is the only
combination a deck can express -- ``model.aux_basis`` is the only place a
fitting set is named.

Double hybrids
~~~~~~~~~~~~~~

``b2plyp``, ``b2gp-plyp`` and ``mpw2plyp`` have analytic gradients, restricted
and all-electron, over an exact or a density-fitted reference. This is not the
hybrid gradient
with a correlation gradient added on, and the difference is the whole of the
work:

- The self-consistent field is over the *semilocal* part only. The perturbative
  term is evaluated once on the converged Kohn-Sham orbitals and never fed back
  into the density, so those orbitals are not stationary with respect to the
  energy being differentiated. An orbital response is therefore mandatory rather
  than an accuracy refinement -- the same Z-vector the MP2 gradient uses, but
  contracted against a *Kohn-Sham* operator.
- That operator contains :math:`V_{xc}`, so two pieces the Hartree-Fock path
  never needed had to exist first: the exchange-correlation kernel
  :math:`f_{xc}` at the GGA rung, which turns the coupled-perturbed
  Hartree-Fock equations into coupled-perturbed Kohn-Sham ones, and the
  skeleton derivative :math:`\partial_R \mathrm{Tr}(P V_{xc}[D])` for a
  :math:`P` that is symmetric and indefinite and is not a density.

Both are validated on their own before the gradient uses them -- the kernel
through a static polarizability against PySCF, the potential derivative against
a finite difference that moves the grid with the nuclei -- because a
disagreement seen only at the end of a Z-vector solve could be any of five
things.

**The perturbative term is fitted when the deck names an auxiliary basis**, and
differentiated as fitted -- an ``Energy`` run and a ``Gradient`` run of one deck
report the same energy. That is worth stating because it was not always true:
until the fitted correlation gradient learned about functionals, a gradient run
had to fall back to exact integrals to stay consistent with its own assembly,
which left the two drivers disagreeing by the fitting error and left the dense
:math:`n_{ao}^4` two-particle density in the one place fitting exists to remove
it. Naming an auxiliary basis is now what makes a double hybrid gradient
affordable, not merely comparable to the reference implementations.

So there are three combinations a deck can ask for, and all three are
differentiated as asked:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Deck
     - Correlation
     - Reference
   * - no auxiliary basis
     - exact
     - exact
   * - ``model.aux_basis``
     - fitted
     - exact
   * - ``+ keywords.scf.density_fitting``
     - fitted
     - fitted

What is refused rather than approximated
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every one of these returns an error naming what is missing, rather than a
number computed from a formula that does not apply:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Case
     - Why it is refused
   * - Meta-GGA over an unrestricted reference
     - The kinetic energy density carries a term per spin that is not built;
       the restricted case is implemented
   * - Range-separated hybrid gradients without an auxiliary basis
     - The second exchange derivative at the screened omega comes from fitting
       a second tensor against ``erf(omega r)/r``; the exact-ERI path has no
       equivalent second pass. Give an auxiliary basis
   * - Unrestricted density fitting, Hartree-Fock or Kohn-Sham
     - The fitted J and K are written for one density, so this is refused for
       the energy and therefore for the gradient. The fitted *gradient* does
       carry both spin channels already -- it is the SCF that is missing
   * - Spin-scaled MP2 (SCS, SOS)
     - The amplitudes enter the response equations, where the two spin cases
       are no longer separable, so the scaled gradient is not the unscaled one
       rescaled
   * - Double hybrid over an open shell
     - Needs an unrestricted MP2 relaxed density and a spin-resolved response,
       neither of which is built. The *energy* is closed-shell only for the
       same reason
   * - Meta-GGA or range-separated double hybrid
     - The kernel exists at the GGA rung. A meta-GGA one needs
       :math:`f_{xc}` in :math:`\tau`; a range-separated one needs the response
       operator and the potential derivative at the screened omega. The three
       double hybrids carried here are all GGA-based, so this is reachable only
       by composing your own
   * - Frozen-core double hybrid
     - Not because the blocks are missing -- plain MP2 builds them and takes
       its frozen-core gradient -- but because the perturbative term's
       assembly has not been taught, or validated with, a frozen core over a
       Kohn-Sham operator with its kernel in the response. The *energy* does
       honour ``freeze_core``, so this refuses rather than returning the
       all-electron gradient of a frozen-core energy
   * - Coupled cluster
     - Needs the Lambda amplitudes, which are not implemented

Geometry Optimization
---------------------

- **Optimizer**: DL-FIND, through libdlfind, behind ``MQC_ENABLE_DLFIND`` (off by
  default -- it is LGPL and this program is MIT)
- **Algorithms**: L-BFGS, conjugate gradient, steepest descent
- **Coordinates**: Cartesian, HDLC (internals within each monomer, Cartesians
  between them), DLC
- **Fragmented systems**: MBE and GMBE gradients drive it unchanged. The term list
  is frozen at the starting geometry so distance screening cannot change the energy
  expression mid-run
- **Backends**: any method with a gradient -- xTB, cuEST, and the CPU ab initio
  backend, whose coverage and refusals are the table above. There is no
  optimizer-level restriction on the backend: whatever refuses a gradient
  refuses it on the first step, with its own message
- **Output**: optimized structure, a trajectory, and a machine-readable record of
  the run
- See :doc:`geometry_optimization`

Hessian Calculations
--------------------

There are two, and the backend picks -- a request that cannot be honoured
analytically is not refused, it falls back.

- **Analytic**, for a restricted Hartree-Fock or Kohn-Sham reference: second
  derivatives without displacing anything, built in the pieces the standard
  decomposition uses (nuclear repulsion, the per-atom perturbation, the
  coupled-perturbed solve, the explicit second-derivative assembly) so it can be
  compared with PySCF's ``hessian.rhf`` and ``hessian.rks`` stage by stage
  rather than only at the end. Taken whenever the calculation is restricted, not
  density fitted, and has no MP2 or coupled-cluster correlation, no continuum
  solvent and no hydrogen caps. See :doc:`analytic_hessians` for which
  functionals qualify.
- **Semi-numerical** otherwise: central differences of *analytic gradients*,
  ``H[i,j] = (g_j(x_i + h) - g_j(x_i - h)) / 2h``. Only one derivative is taken
  numerically, which keeps most of the digits that differencing energies twice
  would lose.

Why the analytic one is worth having: the finite-difference Hessian costs
``6N+1`` gradient evaluations and inherits each one's convergence noise amplified
by ``1/h``. That amplification lands hardest on the low-frequency modes, which is
where the rigid-rotor harmonic-oscillator partition function is most sensitive,
so the noise ends up in the thermochemistry numbers people quote -- and it makes
a transition-state search unreliable, since that needs one negative eigenvalue
whose sign a noisy near-zero mode can flip.

- **Configurable displacement**: step size for the semi-numerical path
  (default: 0.001 Bohr)
- **Hydrogen cap redistribution**: Proper Hessian transformation for capped atoms
- **Units**: Hartree/Bohr²
- **Applications**: Vibrational frequencies, reaction path analysis

Properties
==========

Asked for under ``properties``, beside ``keywords`` rather than inside it: these
change nothing about the energy, so ``driver`` stays ``"energy"``. There are
three.

- **Bonding analysis** -- bond orders and valences over **QUAO**, the
  quasi-atomic orbitals, a basis in which every orbital belongs to one atom. Its
  ``energy_decomposition`` option adds **IEDA**, the intrinsic energy
  decomposition analysis of Del Angel Cruz, Gordon and Ruedenberg: summed over
  that same basis, it reports what each atom contributes on its own and what
  exists only because the atoms are bonded.
  ``properties.bonding_analysis``.
- **Fukui functions**, the derivative of the density with respect to the electron
  count at fixed nuclei. Three values rather than one, because that derivative
  has a kink at every integer: adding an electron and removing one are different
  questions. Computed over Hartree-Fock or any Kohn-Sham functional, with the
  two ions unrestricted, and condensed onto atoms by a population scheme
  (``population``, default ``chelpg``). ``properties.fukui``.
- **Atomic charges**, Mulliken or CHELPG, partitioning the density the
  calculation already converged -- Hartree-Fock or any Kohn-Sham functional,
  restricted or unrestricted, at no cost beyond the SCF that was going to run
  anyway. The object is the request and ``scheme`` only says how, so
  ``"charges": {}`` is a valid ask. An unrestricted reference also gets Mulliken
  spin populations from ``P_alpha - P_beta``; there is no CHELPG counterpart,
  because that scheme fits the electrostatic potential and the total density
  alone determines it. ``properties.charges``.

Charges have a second surface as well: the Python interface reaches them
directly, from a closed-shell Hartree-Fock density only, for the trial-partition
loop that :doc:`charges_and_bond_orders` describes -- along with the
Wiberg-Mayer bond orders beside them, which have no deck surface at all. The
restriction there is the entry point rather than the schemes, which take a
density matrix and do not care what produced it.

Two more exist as machinery without a user-facing surface: distributed multipole
analysis, which builds the EFP fragment potentials, and orbital localization
(Boys/Foster and Edmiston-Ruedenberg), used by AFO and the EFP potentials.

Hydrogen Capping
================

Automatic treatment of broken covalent bonds:

**How it works:**

1. System connectivity analyzed from bond list
2. Broken bonds identified when fragments are extracted
3. Hydrogen caps placed at bond-breaking positions
4. Cap positions computed from original bond geometry

**Gradient/Hessian redistribution:**

- Forces on caps redistributed to original heavy atoms
- Maintains proper derivative continuity
- Transparent to end user - handled automatically

**Supported bond types:**

- C-C, C-N, C-O single bonds
- Configurable cap distance (typically 1.09 Å for C-H), ``cap_scale``

**Or don't cap: adjusted frozen orbitals (AFO).** ``bond_breaking: "afo"`` takes
the other route FMO offers. Rather than terminating the fragment with a hydrogen,
a small model system around the cut bond is built, closed off with hydrogens,
solved and localized, and the orbital sitting on the bond is lifted out of it and
frozen -- so the bond is represented to both sides by the same orbital. The
default remains ``"caps"``.

Input/Output
============

Input Formats
-------------

A single JSON deck describes the whole calculation:

.. code-block:: json

   {
     "molecules": [{
       "xyz": "system.xyz",
       "molecular_charge": 0,
       "molecular_multiplicity": 1
     }],
     "model": {
       "method": "gfn2"
     },
     "keywords": {
       "fragmentation": {
         "method": "MBE",
         "level": 3,
         "cutoffs": {
           "dimer": 5.0,
           "trimer": 4.0
         }
       }
     },
     "driver": "Gradient"
   }

``driver`` is a string, not an object, and everything that configures a method
lives under ``keywords``. Cutoffs are named by n-mer -- ``dimer``, ``trimer`` --
rather than numbered. See :doc:`input_files` for the complete reference.

The ``.mqc`` keyword format and its ``mqc_prep.py`` generator were removed in
0.2.0; ``mqc_docs/source/input_files.rst`` carries the migration table.

Output Formats
--------------

**JSON results:**

- Per-fragment energies and properties
- MBE/GMBE breakdown by n-mer level
- Fragment distances (for screening analysis)
- Total energy, gradient, Hessian

**Log output:**

- Fragment generation statistics
- Screening statistics (fragments kept/discarded)
- Timing information
- MPI rank information

Parallelization
===============

MPI Parallelization
-------------------

**Hierarchical design:**

- **World communicator**: All MPI ranks
- **Node communicator**: Ranks within physical compute nodes
- **Distribution strategy**:

  - Fragments assigned to node leaders
  - Node leaders distribute to local ranks
  - Load balancing across nodes

**Supported modes:**

- Multi-node calculations
- Single-node multi-core
- Serial execution (for debugging)

Fragment Distribution
---------------------

**Serial mode:**

- Single rank computes all fragments
- Useful for small systems, debugging

**MPI coordinator-worker:**

- Coordinator (rank 0) generates fragments
- Workers receive fragment assignments
- Results aggregated on coordinator

**Unfragmented calculations:**

- Direct calculation on rank 0
- No MPI overhead for single-molecule systems

Configuration Options
=====================

Everything below is JSON, and everything that configures a method lives under
``keywords``. The full reference is :doc:`input_files`; this is the shape.

.. code-block:: json

   {
     "keywords": {
       "scf":      {"maxiter": 100, "tolerance": 1e-8, "guess": "sad"},
       "dft":      {"grid_level": 3, "screening_tolerance": 1e-12, "block_size": 512},
       "cc":       {"maxiter": 50, "spin_adapted": true},
       "hessian":  {"finite_difference_displacement": 0.001},
       "pcm":      {"method": "cpcm", "dielectric": 78.4},
       "fragmentation": {
         "method": "MBE",
         "level": 3,
         "cutoffs": {"dimer": 5.0, "trimer": 4.0}
       }
     }
   }

Fragmentation Cutoffs
---------------------

Per-level distance cutoffs in Angstrom, under
``keywords.fragmentation.cutoffs``. Keys are either the n-mer name or the level
as a decimal string -- ``"dimer"`` and ``"2"`` mean the same level, and one
object may mix the two spellings:

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Level
     - Name
     - Note
   * - 2
     - ``dimer``
     -
   * - 3
     - ``trimer``
     -
   * - 4
     - ``tetramer``
     -
   * - 5
     - ``pentamer``
     -
   * - 6
     - ``hexamer``
     -
   * - 7
     - ``heptamer``
     -
   * - 8
     - ``octamer``
     - The highest level distance screening covers

**Behavior:**

- Negative or zero cutoff = no screening for that level
- Missing cutoff = no screening for that level
- Monomers (1-body) always included regardless of cutoffs

.. note::

   The ``%block ... end`` keyword format these options were once written in went
   away with the ``.mqc`` reader in 0.2.0. A deck written that way is refused by
   the schema validator before the reader sees it;
   ``mqc_docs/source/input_files.rst`` carries the migration table.

System Requirements
===================

**Supported platforms:**

- Linux (primary target)
- macOS (tested)
- Windows (via WSL)

**Compiler requirements:**

- Modern Fortran compiler (gfortran 11+, Intel ifort/ifx)
- C compiler for dependencies
- CMake 3.20+

**Dependencies:**

- MPI library (OpenMPI, MPICH, Intel MPI)
- BLAS/LAPACK
- tblite (for XTB methods)
- Optional: DFTD4, mctc-lib, multicharge

**Memory considerations:**

- Fragment count scales combinatorially with system size
- Distance screening reduces memory footprint
- Large systems (>20 monomers) may require HPC resources

Limitations and Future Work
============================

Current Limitations
-------------------

1. **Excited states**: none. No TD-DFT, no CIS, no linear-response excitation
   energies. The coupled-perturbed solver in ``mqc_libcint_response`` already
   handles an electric-field perturbation, which is the piece such a method
   would build on, but nothing consumes it yet
2. **Relativistic Hamiltonians**: none -- no ZORA, DKH or X2C, and no spin-orbit
   coupling operator. An effective core potential is the nearest thing
   available and is a real one: a set fitted to a relativistic reference
   carries much of the effect for the valence, which is why ECPs are used on
   heavy elements as much for that as for the cost. It is not a relativistic
   Hamiltonian and does not become one
3. **Dispersion corrections**: no D3 or D4 for the ab initio path. The xTB
   methods carry their own through tblite; a Kohn-Sham number from a functional
   that does not include dispersion itself does not get it from anywhere here.
   VV10 is the exception: the ``-V`` functionals carry their own non-local
   correlation, evaluated for the energy, the nuclear gradient and the
   analytic Hessian, so a ``-V`` functional can be optimized and its
   frequencies computed -- see :doc:`analytic_hessians`
4. **Multireference dynamic correlation**: CASSCF and ORMAS-SCF give the
   reference, and there is no NEVPT2, CASPT2 or MRCI on top of it
5. **Local correlation**: no DLPNO or equivalent, so coupled cluster is
   canonical and scales like it
6. **Periodic boundaries**: not implemented
7. **AIMD**: keywords are defined and reach the driver config, but there is no
   propagator behind them
8. **Analytic second derivatives**: restricted references only, and within
   those, everything except density fitting -- so Hartree-Fock, LDA, GGA,
   meta-GGA, hybrids, range-separated hybrids and the VV10 ``-V`` functionals
   are covered, and any open shell is not. What is not covered falls back to
   the semi-numerical path rather than failing. The grid response is omitted,
   as it is in the reference this is checked against. Agreement with PySCF is
   1e-8 at STO-3G and loosens to 2e-5 on cc-pVDZ for a functional, which is
   quadrature rather than the derivatives -- worth 0.16 cm-1 at worst. See
   :doc:`analytic_hessians`
9. **SCF convergence aids**: DIIS and level shifting, and no damping, Fermi
   smearing or second-order fallback. The shift is tapered off before
   convergence, so what it costs is iterations and not the orbital energies --
   see :doc:`scf_convergence`. It is a CPU-path feature; the GPU backend accepts
   the keyword without applying it

Planned Features
----------------

1. **Additional QC methods**:

   - Density-fitted unrestricted coupled cluster. The conventional path takes an
     unrestricted reference; the fitted one does not, because its three-index
     block carries no spin blocks, and a deck asking for both is refused rather
     than given the restricted answer built from the alpha orbitals
   - Unrestricted double hybrids, whose perturbative term keeps them closed-shell
   - F12 variants: these parse but have no implementation

2. **Second derivatives beyond restricted Hartree-Fock**: unrestricted,
   density-fitted and MP2 Hessians. The Kohn-Sham one is written and validated
   against PySCF to 1e-8 for every rung through meta-GGA, hybrids and
   range-separated hybrids included -- what is left there is wiring it to
   ``driver: Hessian``, not the derivatives. See :doc:`analytic_hessians`

3. **Advanced dynamics**:

   - AIMD implementation
   - Thermostats/barostats
   - Trajectory analysis

4. **Property calculations**:

   - Polarizabilities -- the coupled-perturbed solve exists and is used to
     validate the kernel, but no deck can ask for one
   - NMR chemical shifts

Performance Notes
=================

**Fragment screening impact:**

For a 20-water cluster with cutoffs ``2=5.0, 3=4.0``:

- Without screening: ~8,000+ fragments
- With screening: ~500-1,000 fragments (depends on geometry)
- Speedup: 5-10x reduction in compute time

**Scaling characteristics:**

- MBE(2): O(N²) fragments
- MBE(3): O(N³) fragments
- GMBE: Intersection enumeration can be expensive for large max_intersection_level
- Distance screening: Linear cost, high benefit

**Best practices:**

1. Use distance screening for systems but check for convergence of correction energies
2. Start with lower fragmentation levels (MBE(2) or MBE(3))
3. For overlapping systems, use GMBE(2) with controlled intersection depth
4. Profile with small test systems before production runs

Acknowledgements
================

Metalquicha builds on several external projects, without which much of the above
would not exist:

- **cuEST** -- the GPU quantum chemistry engine behind the ``cuest`` backend,
  `NVIDIA cuEST <https://developer.nvidia.com/cuda/cuda-x-libraries/cuest>`_.
  All GPU Hartree-Fock and Kohn-Sham energies, gradients and the PCM continuum
  are assembled from its APIs.
- **tblite** -- the semi-empirical engine behind the XTB methods,
  `tblite <https://github.com/tblite/tblite>`_, providing GFN1-xTB and GFN2-xTB.
- **DL-FIND** -- the geometry optimizer behind the ``Optimize`` driver,
  `DL-FIND <https://www.chemshell.org/dl-find>`_, reached through
  `libdlfind <https://github.com/digital-chemistry-laboratory/libdlfind>`_.

We are grateful to their authors and maintainers.
