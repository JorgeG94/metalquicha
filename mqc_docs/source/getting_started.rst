.. _getting_started:

===============
Getting Started
===============

This guide is for developers who want to understand the metalquicha codebase, even without prior quantum chemistry or Fortran experience.

(this file was partially generated with an LLM but carefully checked by me, Jorge)

.. contents::
   :local:
   :depth: 2

What is Metalquicha?
====================

**Metalquicha** (Huastec word for "sunflower") is a quantum chemistry code that specializes in **fragmented calculations**. Think of it as a calculator for molecular energies, but smart enough to break large molecules into smaller pieces to make calculations faster.

What Does It Do?
----------------

- **Input**: Molecular structure (atoms and their positions) + calculation settings
- **Process**: Breaks molecules into fragments, calculates energies for each piece
- **Output**: Total energy (or other properties) of the system

Why Fragment?
-------------

Calculating the energy of a 1000-atom molecule is computationally expensive. But if you break it into 10-atom fragments, calculate each separately, and combine the results cleverly (using **Many-Body Expansion**), you can get accurate results much faster.

Quick Start
===========

Build the Code
--------------

.. code-block:: bash

   # Clone and build
   git clone git@github.com:JorgeG94/metalquicha.git
   cd metalquicha
   mkdir build && cd build
   cmake ..
   make -j

   # The executable will be at: build/mqc

Run Your First Calculation
---------------------------

First, you need to use the python script ``mqc_prep.py`` to
generate the ``mqc`` format, from the json files. The reason
behind not using the json directly is because the json
libraries in Fortran can be delicate with different
Fortran compilers, therefore we use a very simple
format so that we can use plain old Fortran.

To generate an input file, simply run:

.. code-block:: bash

   python mqc_prep.py validation/inputs/h3o.json
   python mqc_prep.py validation/inputs/prism.json

Then run the calculations:

.. code-block:: bash

   # Unfragmented calculation (simple hydronium ion)
   ./build/mqc validation/inputs/h3o.mqc

   # Fragmented calculation (water prism)
   ./build/mqc validation/inputs/prism.mqc

The output shows energies and computation details. To
run a parallel calculation, simply use ``mpirun``:

.. code-block:: bash

   mpirun -np nproc ./build/mqc validation/inputs/prism.mqc

Directory Structure
===================

.. code-block:: text

   metalquicha/
   ├── app/                          # Entry point
   │   └── main.f90                  # Program starts here
   ├── src/                          # Source code
   │   ├── core/                     # Core data structures
   │   │   ├── mqc_geometry.f90      # Molecular structures
   │   │   ├── mqc_elements.f90      # Periodic table data
   │   │   └── mqc_result_types.f90  # Energy/gradient results
   │   ├── io/                       # Input/output
   │   │   ├── mqc_config_parser.f90 # Parse .mqc files
   │   │   ├── mqc_xyz_reader.f90    # Read .xyz geometries
   │   │   └── mqc_json.f90          # JSON output
   │   ├── fragmentation/            # Fragment generation
   │   │   ├── mqc_combinatorics.f90      # Math for combinations
   │   │   ├── mqc_fragment_lookup.f90    # Fast fragment indexing
   │   │   ├── mqc_gmbe_utils.f90         # GMBE overlapping fragments
   │   │   ├── mqc_physical_fragment.f90  # Fragment geometry builder
   │   │   ├── mqc_mbe.f90                # Many-Body Expansion logic
   │   │   └── mqc_mbe_io.f90             # Fragment I/O
   │   ├── methods/                  # Quantum chemistry engines
   │   │   ├── mqc_method_xtb.f90    # XTB method (GFN1-xTB)
   │   │   └── mqc_method_base.f90   # Abstract method interface
   │   ├── parallel/                 # MPI parallelization
   │   │   └── mqc_mpi_tags.f90      # MPI message tags
   │   ├── cli/                      # Command line parsing
   │   └── mqc_driver.f90            # Main calculation orchestrator
   ├── test/                         # Unit tests
   ├── validation/                   # Physics validation tests
   ├── mqc_prep.py                   # Convert JSON → .mqc input
   └── run_validation.py             # Run all validation tests

What Lives Where?
-----------------

- **Core types** (``src/core/``) - Define what a molecule, atom, or energy is
- **I/O** (``src/io/``) - Reading input files, writing results
- **Fragmentation** (``src/fragmentation/``) - Breaking molecules into pieces
- **Methods** (``src/methods/``) - The actual quantum chemistry calculations
- **Driver** (``src/mqc_driver.f90``) - Ties everything together

Key Concepts
============

Fragments
---------

A **fragment** is a subset of atoms from the full molecule. For example, a protein with 1000 atoms might be split into:

- Fragment 1: atoms 1-10
- Fragment 2: atoms 11-20
- ... and so on

Many-Body Expansion (MBE)
--------------------------

Instead of calculating the whole molecule, we use:

.. code-block:: text

   E_total ≈ E(monomer1) + E(monomer2) + ...
             + E(dimer12) - E(monomer1) - E(monomer2)
             + ...

**Key idea**: Higher-order interactions (trimers, tetramers) contribute less, so we can stop at dimers or trimers and get good accuracy.

Generalized MBE (GMBE)
----------------------

Regular MBE assumes fragments don't overlap. **GMBE** handles **overlapping fragments** using the **Principle of Inclusion-Exclusion (PIE)**:

.. code-block:: text

   E_total = Σ E(fragments) - Σ E(intersections) + Σ E(triple_intersections) - ...

This is implemented in ``mqc_gmbe_utils.f90`` using a depth-first search (DFS) algorithm.

Hydrogen Capping
----------------

When you break a chemical bond (e.g., C-C), the fragment has a "dangling bond." We add a hydrogen atom at the break point to satisfy valency. This is called **hydrogen capping**.

See ``mqc_physical_fragment.f90:add_hydrogen_caps()``

System Geometry
---------------

The ``system_geometry_t`` type stores:

- All atom coordinates
- Fragment definitions (which atoms belong to which fragment)
- Fragment sizes
- Bonds (and which are broken)

See ``src/fragmentation/mqc_physical_fragment.f90``

Code Flow: How a Calculation Works
===================================

Here's what happens when you run ``./mqc input.mqc``:

Step 1: Entry Point (app/main.f90)
-----------------------------------

.. code-block:: fortran

   program main
     ! 1. Initialize MPI for parallel computation
     call pic_mpi_init()

     ! 2. Parse command line argument (input file)
     call get_command_argument(1, input_file)

     ! 3. Read .mqc input file
     call read_mqc_file(input_file, mqc_config, error)

     ! 4. Convert config to internal representation
     call config_to_system_geometry(mqc_config, sys_geom, error)

     ! 5. Run the calculation
     call run_calculation(world_comm, node_comm, config, sys_geom, bonds)
   end program

Step 2: Driver (src/mqc_driver.f90)
------------------------------------

The driver decides: fragmented or unfragmented?

.. code-block:: fortran

   subroutine run_calculation(...)
     if (sys_geom%fragmented) then
       ! Run MBE or GMBE
       call run_fragmented_calculation(...)
     else
       ! Run single calculation on whole molecule
       call run_unfragmented_calculation(...)
     end if
   end subroutine

Step 3: Fragment Generation
----------------------------

For fragmented calculations:

1. **Generate combinations** (``mqc_combinatorics.f90``)

   - Use binomial coefficients to determine number of fragments
   - Generate all k-mers (monomers, dimers, trimers, etc.)

2. **Build fragment geometries** (``mqc_physical_fragment.f90``)

   - Extract atom coordinates for each fragment
   - Add hydrogen caps if bonds are broken
   - Store in ``system_geometry_t``

3. **For GMBE: Find intersections** (``mqc_gmbe_utils.f90``)

   - Use DFS to enumerate all overlapping cliques
   - Compute PIE coefficients (+1 for odd cliques, -1 for even)
   - Deduplicate intersection atom sets

Step 4: Energy Evaluation
--------------------------

For each fragment:

1. **Select method** (``mqc_method_xtb.f90`` currently)

   - Initialize XTB calculator
   - Convert geometry to XTB format

2. **Run calculation**

   .. code-block:: fortran

      call method%calculate_energy(fragment_geometry, result)

3. **Store result** in ``calculation_result_t`` type

Step 5: Combine Results
------------------------

**MBE** (``mqc_mbe.f90:compute_mbe_energy``):

.. code-block:: fortran

   total_energy = sum(monomers)
                + sum(dimers) - 2*sum(monomers)
                + sum(trimers) - 3*sum(dimers) + 3*sum(monomers)
                ...

**GMBE** uses PIE coefficients:

.. code-block:: fortran

   total_energy = sum(coefficient(i) * energy(i))

Step 6: Output
--------------

Results are written to:

- **Console**: Human-readable summary
- **JSON file** (``output_*.json``): Structured data for validation

Understanding the Modules
==========================

Core Modules (Start Here!)
---------------------------

mqc_geometry.f90
^^^^^^^^^^^^^^^^

Defines ``molecular_geometry_t`` - stores atoms, positions, charges.

**Key functions**:

- ``geometry_from_arrays()`` - Create geometry from atom data
- ``get_atom_symbol()`` - Get element symbol for an atom

mqc_physical_fragment.f90
^^^^^^^^^^^^^^^^^^^^^^^^^^

Defines ``system_geometry_t`` - the master structure for everything.

**Key data**:

- ``total_atoms`` - Total number of atoms
- ``xyz(:,:)`` - Atom coordinates (3, total_atoms)
- ``atomic_numbers(:)`` - Element for each atom
- ``fragment_atoms(:,:)`` - Which atoms are in which fragment
- ``fragment_sizes(:)`` - Size of each fragment
- ``fragmented`` - Boolean: is this a fragmented calculation?

**Key functions**:

- ``build_fragment_from_indices()`` - Extract fragment geometry
- ``add_hydrogen_caps()`` - Add H atoms at broken bonds

Fragmentation Modules (The Core Algorithm!)
--------------------------------------------

mqc_combinatorics.f90
^^^^^^^^^^^^^^^^^^^^^

Pure math - no chemistry here!

**Key functions**:

- ``binomial(n, r)`` - C(n,r) = n choose r
- ``get_nfrags(n_monomers, max_level)`` - How many fragments total?
- ``generate_fragment_list()`` - Generate all k-mer combinations
- ``next_combination()`` - Iterator for combinations

**Example**:

.. code-block:: fortran

   ! For 4 monomers, level 2 (dimers):
   ! Generates: [1,2], [1,3], [1,4], [2,3], [2,4], [3,4]

mqc_fragment_lookup.f90
^^^^^^^^^^^^^^^^^^^^^^^

Hash table for O(1) fragment lookup.

**Why needed?** With 100 monomers and level 3, you have C(100,3) = 161,700 trimers. Finding "which index is trimer [5,23,87]?" needs to be fast!

**Key type**: ``fragment_lookup_t``

- Uses FNV-1a hash function
- Chaining for collision resolution
- Prime-sized table for good distribution

mqc_gmbe_utils.f90
^^^^^^^^^^^^^^^^^^

The most complex module - implements overlapping fragment logic.

**Key functions**:

- ``find_fragment_intersection()`` - Find shared atoms between 2 fragments
- ``generate_intersections()`` - Generate all k-way intersections
- ``gmbe_enumerate_pie_terms()`` - DFS to find all cliques and compute PIE coefficients

**Algorithm (simplified)**:

.. code-block:: text

   For each primary fragment i:
     Start DFS with clique = {i}
     For each candidate j > i that overlaps with current intersection:
       Add j to clique
       Compute new intersection
       If intersection non-empty:
         Coefficient = (+1 if |clique| odd, -1 if even)
         Recurse with expanded clique

mqc_mbe.f90
^^^^^^^^^^^

Energy combination logic for standard MBE.

**Key function**: ``compute_mbe_energy()``

- Loops through levels (monomers, dimers, trimers...)
- Applies inclusion-exclusion formula
- Returns total energy

I/O Modules
-----------

mqc_config_parser.f90
^^^^^^^^^^^^^^^^^^^^^

Parses ``.mqc`` input files (section-based format).

**Structure**:

.. code-block:: fortran

   type :: mqc_config_t
     character(len=:), allocatable :: method
     integer :: nlevel                    ! Max MBE level
     logical :: allow_overlapping_fragments
     type(geometry_t) :: geometry
     integer, allocatable :: fragment_indices(:,:)
     ! ... many more fields
   end type

mqc_xyz_reader.f90
^^^^^^^^^^^^^^^^^^

Reads ``.xyz`` geometry files (standard format):

.. code-block:: text

   3
   Water molecule
   O   0.0   0.0   0.0
   H   0.0   0.0   0.96
   H   0.93  0.0  -0.24

Method Modules
--------------

mqc_method_base.f90
^^^^^^^^^^^^^^^^^^^

Abstract interface - defines what a "quantum chemistry method" must provide.

.. code-block:: fortran

   type, abstract :: method_base_t
   contains
     procedure(calculate_energy_interface), deferred :: calculate_energy
     procedure(calculate_gradient_interface), deferred :: calculate_gradient
   end type

mqc_method_xtb.f90
^^^^^^^^^^^^^^^^^^

Concrete implementation using the tblite library (GFN1-xTB method).

**Key function**: ``calculate_energy()`` converts geometry → calls tblite → returns energy

Where to Start Reading
======================

For Understanding the Big Picture
----------------------------------

1. **app/main.f90** (100 lines) - See the entry point
2. **src/mqc_driver.f90** (500 lines) - See how calculations are orchestrated
3. **src/fragmentation/mqc_physical_fragment.f90** (400 lines) - Understand fragments

For Understanding Fragmentation
--------------------------------

1. **mqc_combinatorics.f90** - Pure math, easiest to understand
2. **mqc_fragment_lookup.f90** - Data structure for indexing
3. **mqc_mbe.f90** - See how energies are combined
4. **mqc_gmbe_utils.f90** - Most complex, save for last

For Understanding I/O
---------------------

1. **mqc_config_parser.f90** - How we parse ``.mqc`` files
2. **mqc_xyz_reader.f90** - How we read geometries

For Understanding Quantum Chemistry Methods
--------------------------------------------

1. **mqc_method_base.f90** - Interface definition
2. **mqc_method_xtb.f90** - Concrete implementation

Running Examples
================

Example 1: Unfragmented Calculation
------------------------------------

.. code-block:: bash

   # Simple hydronium ion
   ./build/mqc validation/inputs/h3o.mqc

**What happens**:

- Reads 4 atoms (H₃O⁺)
- No fragmentation
- Calls XTB on the whole molecule
- Prints energy: -5.773131 Ha

Example 2: Fragmented Water Prism
----------------------------------

.. code-block:: bash

   ./build/mqc validation/inputs/prism.mqc

**What happens**:

- Reads 18 atoms (6 water molecules)
- Fragments into 6 monomers
- Generates 6 monomers + 15 dimers = 21 fragments
- Calculates energy for each fragment in parallel
- Combines using MBE(2) formula
- Prints total energy: -34.673668 Ha

Example 3: GMBE with Overlapping Fragments
-------------------------------------------

.. code-block:: bash

   ./build/mqc validation/inputs/overlapping_gly3.mqc

**What happens**:

- Reads glycine tripeptide with overlapping fragments
- Uses GMBE(1) with PIE enumeration
- Finds intersections between overlapping fragments
- Applies inclusion-exclusion principle
- Prints total energy

Prepare Your Own Input
-----------------------

.. code-block:: bash

   # Create a JSON input
   cat > my_input.json << EOF
   {
     "schema": {"name": "mqc-frag", "version": "1.0"},
     "molecules": [{
       "xyz": "water.xyz",
       "molecular_charge": 0,
       "molecular_multiplicity": 1,
       "fragments": [[0,1,2]]
     }],
     "model": {"method": "XTB-GFN1"},
     "driver": "Energy"
   }
   EOF

   # Convert to .mqc format
   python3 mqc_prep.py my_input.json

   # Run calculation
   ./build/mqc my_input.mqc

Adding New Features
===================

Example: Add a New Quantum Chemistry Method
--------------------------------------------

1. **Create new method file**: ``src/methods/mqc_method_mynew.f90``

.. code-block:: fortran

   module mqc_method_mynew
     use mqc_method_base, only: method_base_t
     implicit none

     type, extends(method_base_t) :: mynew_method_t
     contains
       procedure :: calculate_energy => mynew_calculate_energy
       procedure :: calculate_gradient => mynew_calculate_gradient
     end type

   contains

     subroutine mynew_calculate_energy(this, geom, result)
       class(mynew_method_t), intent(inout) :: this
       type(molecular_geometry_t), intent(in) :: geom
       type(calculation_result_t), intent(out) :: result

       ! Your implementation here
       result%energy = 0.0_dp  ! Replace with actual calculation
     end subroutine

   end module

2. **Register in driver**: Edit ``src/mqc_driver.f90``

.. code-block:: fortran

   use mqc_method_mynew, only: mynew_method_t

   ! In select_method():
   case ('MYNEW')
     allocate(mynew_method_t :: method)

3. **Add to CMakeLists.txt**: ``src/methods/CMakeLists.txt``

.. code-block:: cmake

   target_sources(
     ${main_lib}
     PRIVATE mqc_method_mynew.f90
             ...
   )

4. **Test**: Create validation test in ``validation/inputs/test_mynew.mqc``

Example: Add New Fragment Selection Algorithm
----------------------------------------------

Want to select fragments based on distance cutoffs instead of combinatorial generation?

1. **Add to** ``mqc_physical_fragment.f90``:

.. code-block:: fortran

   subroutine generate_distance_based_fragments(sys_geom, cutoff, fragments)
     type(system_geometry_t), intent(in) :: sys_geom
     real(dp), intent(in) :: cutoff
     integer, allocatable, intent(out) :: fragments(:,:)

     ! Your algorithm:
     ! - Loop through all monomer pairs
     ! - Compute distance between centers of mass
     ! - If distance < cutoff, include dimer
   end subroutine

2. **Call from driver** in the fragmentation setup section

3. **Add config option** in ``mqc_config_parser.f90`` to enable it

Testing
=======

Unit Tests
----------

Located in ``test/`` directory. Each module has a corresponding test.

.. code-block:: bash

   # Run all tests
   cd build
   ctest

   # Run specific test
   ./test/test_mqc_frag_utils

Validation Tests
----------------

Located in ``validation/``. These test physics correctness.

.. code-block:: bash

   # Run all validation tests
   python3 run_validation.py

   # Run specific test
   ./build/mqc validation/inputs/prism.mqc

**What validation tests check**:

- Energy matches expected value (within tolerance)
- JSON output is well-formed
- All fragments were evaluated
- PIE coefficients are correct (for GMBE)

Adding a New Test
-----------------

1. **Create input file**: ``validation/inputs/my_test.mqc``

2. **Add to** ``validation_tests.json``:

.. code-block:: json

   {
     "name": "My test case",
     "input": "validation/inputs/my_test.mqc",
     "expected_energy": -123.456789,
     "type": "fragmented"
   }

3. **Run**: ``python3 run_validation.py``

The script will:

- Run the calculation
- Parse JSON output
- Compare energy to expected value
- Report PASS/FAIL

Common Patterns in the Code
============================

Error Handling
--------------

.. code-block:: fortran

   type(error_t) :: error

   call some_function(..., error)
   if (error%has_error()) then
     call logger%error(error%get_message())
     return
   end if

Logging
-------

.. code-block:: fortran

   use pic_logger, only: logger => global_logger

   call logger%info("Starting calculation...")
   call logger%debug("Fragment "//to_char(i)//" has "//to_char(n_atoms)//" atoms")
   call logger%error("Invalid input: "//error_msg)

Memory Management
-----------------

Fortran 2003+ has automatic deallocation, but we explicitly clean up:

.. code-block:: fortran

   type(system_geometry_t) :: sys_geom

   ! ... use sys_geom ...

   call sys_geom%destroy()  ! Explicitly deallocate

MPI Parallelization
-------------------

.. code-block:: fortran

   ! Each MPI rank gets a subset of fragments
   do i = world_comm%rank() + 1, n_fragments, world_comm%size()
     call calculate_fragment_energy(fragments(i), results(i))
   end do

   ! Gather results from all ranks
   call world_comm%allgatherv(local_results, all_results)

Glossary
========

- **MBE**: Many-Body Expansion - method to approximate total energy using fragments
- **GMBE**: Generalized MBE - handles overlapping fragments
- **PIE**: Principle of Inclusion-Exclusion - combinatorial method for counting overlaps
- **Monomer**: 1-body fragment (single unit)
- **Dimer**: 2-body fragment (pair of units)
- **Trimer**: 3-body fragment (triple of units)
- **k-mer**: k-body fragment (k units)
- **XTB**: Extended Tight Binding - fast semi-empirical quantum chemistry method
- **GFN1**: Geometry, Frequency, Non-covalency version 1 (XTB parameter set)
- **Hartree (Ha)**: Unit of energy in atomic units (1 Ha ≈ 27.2 eV)

Additional Resources
====================

- **pic library**: https://github.com/JorgeG94/pic - Fortran utilities used throughout
- **tblite**: https://github.com/tblite/tblite - XTB quantum chemistry engine
- **README.md**: Build instructions and overview

Getting Help
============

- **Check existing validation tests** in ``validation/inputs/`` for examples
- **Read module documentation** - each module has detailed comments
- **Run with verbose logging**: Set ``log_level = Verbose`` in your ``.mqc`` file
- **Check JSON output**: Contains detailed fragment energies and coefficients

Next Steps
==========

1. **Build and run**: Get the code compiling and run validation tests
2. **Read** ``app/main.f90``: Understand the entry point
3. **Trace a simple calculation**: Follow the code path for ``h3o.mqc`` (unfragmented)
4. **Trace a fragmented calculation**: Follow the code path for ``prism.mqc`` (MBE)
5. **Understand GMBE**: Read ``mqc_gmbe_utils.f90`` (most complex part)
6. **Add a feature**: Pick something small to modify and experiment

Happy coding! 🌻
