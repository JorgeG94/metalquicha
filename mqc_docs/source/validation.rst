.. _validation:

===========================
Physics Validation Testing
===========================

(this file was partially generated with an LLM but carefully checked by me, Jorge)

This document describes the physics validation test suite for metalquicha, which ensures that code changes don't introduce numerical regressions.

.. contents::
   :local:
   :depth: 2

Overview
========

The validation suite runs quantum chemistry calculations and compares the computed total energies against reference values. All tests must pass before merging code changes.

**Key features**:

- Automated testing of fragmented and unfragmented calculations
- Energy comparison with configurable tolerance (default: 1.0e-6 Hartree)
- Support for multi-molecule calculations (conformers, isomers)
- CI/CD integration ready
- Detailed logging for debugging failures

Directory Structure
===================

.. code-block:: text

   validation/
   └── inputs/
       ├── *.json          # Input JSON files (geometry, fragments, parameters)

   validation_tests.json   # Test manifest with expected energies
   run_validation.py       # Validation test runner script
   validation_logs/        # Generated: stdout/stderr from each calculation

Output Files
------------

The validation script generates:

- ``output_*.json`` - Calculation results for each test
- ``validation_logs/*.log`` - Stdout/stderr from each calculation (for debugging)

Workflow
========

Adding a New Test
-----------------

1. Create Input File
^^^^^^^^^^^^^^^^^^^^

Input decks live under ``validation/inputs/`` in a tree by hardware, engine and
method, so that listing a directory answers "what is covered" rather than
scrolling past two hundred files::

   validation/inputs/
   |-- sample_inputs/          geometries, shared by every backend
   |-- cpu/mqc/rhf/            the libcint CPU backend, one directory per method
   |   |-- df-hf/  uhf/  mp2/  ri-mp2/  ccsd/  ccsd-t/  ri-ccsd/  ...
   |-- cpu/tblite/gfn1/        the semi-empirical backend
   `-- gpu/cuest/dft/          the GPU backend

The geometries deliberately stay at the top, shared: a deck three levels down
reaches one with ``../../../sample_inputs/``. Paths inside a deck are resolved
relative to that deck, so nothing here is absolute and the whole tree can be
moved or copied.

Decks for the CPU suite are **generated**, so a new one goes in
``tools/cpu_validation/gen_cpu_validation.py`` rather than here -- a deck added
by hand under ``cpu/mqc/`` is deleted by the next regeneration. Add a deck by hand only for a backend the generator does not
cover:

.. code-block:: bash

   cat > validation/inputs/cpu/tblite/gfn1/my_test.json << EOF
   {
     "schema": {"name": "mqc-frag", "version": "1.0"},
     "molecules": [{
       "xyz": "../../../sample_inputs/my_molecule.xyz",
       "molecular_charge": 0,
       "molecular_multiplicity": 1
     }],
     "model": {"method": "XTB-GFN1"},
     "driver": "Energy"
   }
   EOF

2. Generate Reference Energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the calculation to get the reference energy:

.. code-block:: bash

   # Run calculation
   ./build/mqc validation/inputs/cpu/tblite/gfn1/my_test.json

   # Extract energy from output JSON
   cat output_my_test.json | grep '"total_energy"' | head -1

3. Update Test Manifest
^^^^^^^^^^^^^^^^^^^^^^^^

Edit ``validation_tests.json`` to add your test with expected energy:

.. code-block:: json

   {
     "description": "Physics validation tests for metalquicha",
     "tolerance": 1.0e-6,
     "tests": [
       {
         "name": "My test case",
         "input": "inputs/cpu/tblite/gfn1/my_test.json",
         "expected_energy": -123.456789,
         "type": "unfragmented"
       }
     ]
   }

**Test types**:

- ``unfragmented``: Single molecule, no fragmentation
- ``fragmented``: Many-body expansion calculation
- ``multi_molecule``: Multiple independent molecules (conformers, isomers)

4. Run Validation
^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Run all tests
   python3 run_validation.py

   # Verbose output
   python3 run_validation.py -v


   # Use custom executable
   python3 run_validation.py --exe ./my_build/mqc

Running the Validation Suite
=============================

Basic Usage
-----------

.. code-block:: bash

   # Run all validation tests
   python3 run_validation.py

Output:

.. code-block:: text

   ============================
   Metalquicha Validation Suite
   ============================

   ✓ Converted validation/inputs/cpu/tblite/gfn1/h3o.json
   ✓ Converted validation/inputs/cpu/tblite/gfn1/prism.json
   ...

   Running 10 validation tests...

   [ 1/10] Unfragmented glycine tripeptide ... ✓ PASSED
   [ 2/10] Charged cluster MBE ............... ✓ PASSED
   [ 3/10] Glycine 10-mer MBE ................ ✓ PASSED
   ...

   ============================
   Summary: 10/10 tests passed
   ============================

Advanced Options
----------------

.. code-block:: bash

   # Verbose mode (shows all output)
   python3 run_validation.py -v



   # Use custom executable
   python3 run_validation.py --exe /path/to/mqc

   # Specify custom test manifest
   python3 run_validation.py --tests my_tests.json

Test Types
==========

Unfragmented Tests
------------------

Single molecule calculations without fragmentation.

**Example**: ``validation/inputs/cpu/tblite/gfn1/h3o.json``

.. code-block:: json

   {
     "name": "Hydronium ion unfragmented",
     "input": "inputs/cpu/tblite/gfn1/h3o.json",
     "expected_energy": -5.773131213617977,
     "type": "unfragmented"
   }

**JSON output structure**:

.. code-block:: json

   {
     "h3o": {
       "total_energy": -5.773131213617977
     }
   }

Fragmented Tests (MBE/GMBE)
---------------------------

Many-body expansion calculations with fragment-based energies.

**Example**: ``validation/inputs/cpu/tblite/gfn1/prism.json`` (MBE(2) on water prism)

.. code-block:: json

   {
     "name": "Water prism MBE",
     "input": "inputs/cpu/tblite/gfn1/prism.json",
     "expected_energy": -34.6736678571,
     "type": "fragmented"
   }

**JSON output structure** (abbreviated):

.. code-block:: text

   {
     "prism": {
       "total_energy": -34.6736678571,
       "levels": [
         {
           "level": 1,
           "n_fragments": 6,
           "energies": [...]
         },
         {
           "level": 2,
           "n_fragments": 15,
           "energies": [...]
         }
       ]
     }
   }

Multi-Molecule Tests
--------------------

Multiple independent molecules (conformers, isomers) in one input.

**Example**: ``validation/inputs/cpu/tblite/gfn1/multi_frag.json``

.. code-block:: json

   {
     "name": "Multi-fragment calculation",
     "input": "inputs/cpu/tblite/gfn1/multi_frag.json",
     "expected_energies": {
       "molecule_1": -34.6736678571,
       "molecule_2": -34.6736678571
     },
     "type": "multi_molecule"
   }

**JSON output structure**:

.. code-block:: json

   {
     "multi_frag": {
       "molecule_1": {
         "total_energy": -34.6736678571
       },
       "molecule_2": {
         "total_energy": -34.6736678571
       }
     }
   }

Energy Extraction
=================

The validation script extracts the top-level ``total_energy`` field for each calculation type.

Unfragmented Calculations
--------------------------

.. code-block:: json

   {
     "basename": {
       "total_energy": -123.456789
     }
   }

Extracted value: ``-123.456789``

Fragmented Calculations (MBE/GMBE)
-----------------------------------

.. code-block:: text

   {
     "basename": {
       "total_energy": -123.456789,
       "levels": [...]
     }
   }

Extracted value: ``-123.456789`` (intermediate level energies are not validated)

Multi-Molecule Calculations
----------------------------

.. code-block:: json

   {
     "basename": {
       "molecule_1": {"total_energy": -123.45},
       "molecule_2": {"total_energy": -124.56}
     }
   }

Extracted values: ``molecule_1: -123.45``, ``molecule_2: -124.56``

Tolerance Configuration
=======================

Default energy tolerance is ``1.0e-6`` Hartree, configurable in ``validation_tests.json``:

.. code-block:: text

   {
     "description": "Physics validation tests for metalquicha",
     "tolerance": 1.0e-6,
     "tests": [...]
   }

The test passes if:

.. code-block:: python

   abs(computed_energy - expected_energy) < tolerance

Validation Test Examples
=========================

The validation suite includes several representative test cases:

GMBE(1) Test (Overlapping Fragments)
-------------------------------------

**Test**: ``overlapping_gly3.json``

- **System**: Glycine tripeptide with overlapping fragments
- **Method**: GMBE(1) with inclusion-exclusion principle
- **Expected energy**: -47.0192718920 Ha
- **Purpose**: Validates PIE coefficient enumeration and overlapping fragment logic

GMBE(3) Test (Higher-Order)
----------------------------

**Test**: ``nlevel_3_ov_decane.json``

- **System**: Decane molecule with overlapping fragments
- **Method**: GMBE(3) (trimers)
- **Expected energy**: -33.0506139740 Ha
- **Purpose**: Validates higher-order GMBE calculations

Standard MBE(2) Test
--------------------

**Test**: ``prism.json``

- **System**: Water prism (6 water molecules)
- **Method**: MBE(2) (dimers)
- **Expected energy**: -34.6736678571 Ha
- **Purpose**: Validates standard many-body expansion

Large System Test
-----------------

**Test**: ``w20_isomer.json``

- **System**: Water 20-mer isomer
- **Method**: MBE(2)
- **Expected energy**: -115.6850246841 Ha
- **Purpose**: Validates performance on larger systems

Charged Cluster Test
--------------------

**Test**: ``charged_cluster.json``

- **System**: Charged molecular cluster
- **Method**: MBE(2)
- **Expected energy**: -45.7161383790 Ha
- **Purpose**: Validates fragment charge handling

CI Integration
==============

For GitHub Actions or similar:

.. code-block:: yaml

   name: Validation Tests

   on: [push, pull_request]

   jobs:
     validate:
       runs-on: ubuntu-latest
       steps:
         - uses: actions/checkout@v2

         - name: Build metalquicha
           run: |
             mkdir build && cd build
             cmake ..
             make -j

         - name: Run validation tests
           run: |
             python3 run_validation.py --exe ./build/mqc

The script exits with code 1 on any test failure, causing CI to fail.

Troubleshooting
===============

Common Issues
-------------

**"No JSON files found"**
   Create ``validation/inputs/`` directory and add JSON files


**Energy mismatch**
   Check for:

   - Different compiler optimization levels (``-O2`` vs ``-O3``)
   - Different dependency versions (tblite, LAPACK)
   - Numerical instabilities in the method
   - Incorrect reference energy (re-run and update)

**Calculation timeout**
   Default is 300s per test. Edit ``run_validation.py`` if needed:

   .. code-block:: python

      timeout = 600  # 10 minutes

**Test fails intermittently**
   May indicate:

   - A bug in xtb caused by a compiler version, math library, etc. (this has been observed in MacOS)
   - Numerical instabilities in SCF convergence
   - Compiler-specific floating-point behavior

Debugging Failed Tests
----------------------

1. **Check the log file**:

   .. code-block:: bash

      cat validation_logs/my_test.log

2. **Run test manually with verbose output**:

   .. code-block:: bash

      ./build/mqc validation/inputs/cpu/tblite/gfn1/my_test.json

3. **Check JSON output**:

   .. code-block:: bash

      cat output_my_test.json | python3 -m json.tool

4. **Compare energies**:

   .. code-block:: bash

      python3 -c "import json; print(json.load(open('output_my_test.json')))"

5. **Run with different compilers** to check portability:

   .. code-block:: bash

      # gfortran
      cmake -DCMAKE_Fortran_COMPILER=gfortran ..
      make
      python3 run_validation.py

      # ifort
      cmake -DCMAKE_Fortran_COMPILER=ifort ..
      make
      python3 run_validation.py

Best Practices
==============

1. **Add tests for new features**: Every new feature should have at least one validation test
2. **Use representative systems**: Include small (fast) and large (realistic) test cases
3. **Document expected behavior**: Add comments in test manifest explaining what each test validates
4. **Update reference energies carefully**: Only update when intentional physics changes occur
5. **Check all test types**: Include unfragmented, fragmented, and multi-molecule examples
6. **Run before committing**: Always run validation suite before pushing changes

Generating Reference Energies
==============================

Best Practices
--------------

1. **Use release build** for consistency:

   .. code-block:: bash

      cmake -DCMAKE_BUILD_TYPE=Release ..
      make

2. **Run multiple times** to ensure reproducibility:

   .. code-block:: bash

      for i in {1..5}; do
        ./build/mqc validation/inputs/cpu/tblite/gfn1/my_test.json
        grep total_energy output_my_test.json
      done

3. **Document the compiler and flags** used:

   .. code-block:: bash

      gfortran --version
      cmake .. -LA | grep CMAKE_Fortran_FLAGS

4. **Cross-validate** with external codes if possible (e.g., run XTB standalone)

Validation Test Manifest Schema
================================

The ``validation_tests.json`` file has the following structure:

.. code-block:: text

   {
     "description": "Physics validation tests for metalquicha",
     "tolerance": 1.0e-6,
     "tests": [
       {
         "name": "Test name (required)",
         "input": "path/to/input.json (required)",
         "expected_energy": -123.456,  // for single molecule
         "expected_energies": {...},   // for multi-molecule
         "type": "unfragmented|fragmented|multi_molecule (required)"
       }
     ]
   }

**Fields**:

- ``description``: Human-readable description of test suite
- ``tolerance``: Energy comparison tolerance (Hartree)
- ``tests``: Array of test objects

**Test object**:

- ``name``: Descriptive test name
- ``input``: Path to the JSON input deck
- ``expected_energy``: Expected total energy (for single molecule tests)
- ``expected_energies``: Dictionary of expected energies (for multi-molecule tests)
- ``type``: Test type (``unfragmented``, ``fragmented``, or ``multi_molecule``)

Contributing New Tests
======================

When contributing new test cases:

1. **Use descriptive names**: ``"GMBE(2) water dimer with H-caps"`` not ``"Test 1"``
2. **Add comments in JSON**: Explain what the test validates
3. **Include .xyz files**: Commit geometry files to ``validation/inputs/``
4. **Document the system**: Add a comment describing the molecular system
5. **Test edge cases**: Include unusual charges, multiplicities, or geometries
6. **Verify portability**: Run on different compilers before submitting

Example contribution:

.. code-block:: json

   {
     "name": "Charged peptide with broken bonds",
     "input": "inputs/cpu/tblite/gfn1/charged_peptide.json",
     "expected_energy": -78.1234567890,
     "type": "fragmented",
     "comment": "Tests H-capping on charged system with multiple broken C-C bonds"
   }
