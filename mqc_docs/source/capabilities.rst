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

Calculation Types
=================

Energy Calculations
-------------------

- **Single-point energies**: Total electronic energy
- **Component breakdown**: SCF, MP2 (future), coupled-cluster (future)
- **Output**: Hartrees, per-fragment energies, MBE contributions

Gradient Calculations
---------------------

- **Analytical gradients**: Energy derivatives with respect to nuclear positions
- **Hydrogen cap redistribution**: Automatic gradient mapping for capped bonds
- **Units**: Hartree/Bohr
- **Applications**: Geometry optimization, molecular dynamics

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

**JSON configuration:**

.. code-block:: json

   {
     "molecules": {[
     "geometry": {
       "file": "system.xyz"
     }]},
     "model": {
       "method": "XTB-GFN2"
     },
     "driver": {
       "type": "gradient"
     },
     "fragmentation": {
       "method": "MBE",
       "level": 3,
       "cutoffs": {
         "2": 5.0,
         "3": 4.0
       }
     }
   }

**Internal .mqc format:**

- Human-readable keyword-value format
- Generated automatically from JSON
- Direct parsing for portability

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

   - Hartree-Fock (via libint)
   - DFT functionals (via libxc, libint)
   - Post-HF methods (MP2, CCSD)

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
