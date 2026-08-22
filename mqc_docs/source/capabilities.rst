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
  spatial orbitals by default, and spin orbitals for checking it against -- which
  are exact for each other and agree to machine precision, so the choice is how a
  number is computed and not which number. An open-shell reference takes the
  spin-orbital path necessarily, since the spin-adapted equations are derived for
  a closed shell. Density fitting is available for the restricted case only: the
  fitted three-index block has no spin blocks, so an unrestricted deck asking for
  it is refused rather than quietly given the restricted answer.
- **Kohn-Sham DFT**: the whole ladder -- LDA, GGA, hybrid, meta-GGA,
  range-separated hybrid and double hybrid -- restricted and unrestricted, over
  `libxc <https://libxc.gitlab.io/>`_, so most of what libxc carries is available
  by name, with friendly names and double-hybrid compositions layered on top.
  Range separation uses libcint's erf-attenuated integrals, which puts ωB97X and
  CAM-B3LYP in reach. Non-local correlation (VV10) and Laplacian-dependent
  meta-GGAs are refused rather than approximated. Grids are Treutler-Ahlrichs
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
- **SAPT0**: the interaction energy of two monomers, decomposed into
  electrostatics, exchange, induction, dispersion and their exchange
  counterparts, in the dimer-centred basis. Two fragments exactly -- see
  :doc:`sapt`.

Basis sets come from the Basis Set Exchange data shipped in ``basis_sets/``, and
whether a set is Cartesian or spherical is taken from the file rather than assumed.

Implicit Solvation
------------------

Implicit solvation models account for solvent effects without explicit solvent molecules:

**Supported models:**

- **PCM**: a polarizable continuum on the cuEST (GPU) backend, for Hartree-Fock
  and Kohn-Sham, restricted and unrestricted. Validated against PySCF's C-PCM;
  configured through ``keywords.pcm``.
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
   * - Kohn-Sham: range-separated hybrid
     - yes
     - **Needs an auxiliary basis.** The exact-ERI path builds no second
       exchange derivative at the screened omega, so the fitted one is the
       only one
   * - MP2, restricted
     - yes
     - All four combinations of fitted reference and fitted correlation
   * - Double hybrids (``b2plyp``, ``b2gp-plyp``, ``mpw2plyp``)
     - yes
     - Restricted, all-electron, GGA-based; exact or fitted reference
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
   * - Frozen-core MP2 and RI-MP2
     - The relaxed density gains occupied-frozen and virtual-frozen blocks
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
     - The relaxed density gains occupied-frozen blocks that are not built.
       The *energy* does honour ``freeze_core``, so this refuses rather than
       returning the all-electron gradient of a frozen-core energy
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

- **Numerical Hessians**: Second derivatives via finite differences
- **Configurable displacement**: User-specified step size (default: 0.001 Bohr)
- **Hydrogen cap redistribution**: Proper Hessian transformation for capped atoms
- **Units**: Hartree/Bohr²
- **Applications**: Vibrational frequencies, reaction path analysis

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
- Configurable cap distance (typically 1.09 Å for C-H)

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

Hessian Keywords
----------------

.. code-block:: text

   %hessian
   finite_difference_displacement = 0.001  ! Bohr
   end

AIMD Keywords (Future)
----------------------

.. code-block:: text

   %aimd
   dt = 1.0                      ! Femtoseconds
   nsteps = 1000                 ! MD steps
   initial_temperature = 300.0   ! Kelvin
   output_frequency = 10         ! Write every N steps
   end

SCF Keywords
------------

.. code-block:: text

   %scf
   max_iterations = 100
   convergence_threshold = 1.0e-6
   use_diis = true
   end

Fragmentation Cutoffs
---------------------

.. code-block:: text

   %cutoffs
   2 = 5.0   ! Dimer cutoff (Angstrom)
   3 = 4.0   ! Trimer cutoff
   4 = 3.5   ! Tetramer cutoff
   5 = 3.0   ! Pentamer cutoff
   end

**Supported n-mers:**

- 2: Dimers
- 3: Trimers
- 4: Tetramers
- 5: Pentamers
- 6: Hexamers
- 7: Heptamers
- 8: Octamers

**Behavior:**

- Negative or zero cutoff = no screening for that level
- Missing cutoff = no screening for that level
- Monomers (1-body) always included regardless of cutoffs

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

1. **QC methods**: Only XTB currently integrated (HF, DFT planned)
2. **Periodic boundaries**: Not yet implemented
3. **AIMD**: Keywords defined but implementation pending

Planned Features
----------------

1. **Additional QC methods**:

   - Density-fitted unrestricted coupled cluster. The conventional path takes an
     unrestricted reference; the fitted one does not, because its three-index
     block carries no spin blocks, and a deck asking for both is refused rather
     than given the restricted answer built from the alpha orbitals
   - Unrestricted double hybrids, whose perturbative term keeps them closed-shell
   - F12 variants, and MCSCF: these parse but have no implementation

2. **Advanced dynamics**:

   - AIMD implementation
   - Thermostats/barostats
   - Trajectory analysis

3. **Property calculations**:

   - Dipole moments
   - Polarizabilities
   - NMR chemical shifts

4. **Analysis tools**:

   - Fragment interaction energies

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
