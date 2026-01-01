.. _input_files:

==================
Input File Formats
==================

(this file was partially generated with an LLM but carefully checked by me, Jorge)

Metalquicha supports two input file formats: JSON (for user convenience) and ``.mqc`` (the native format that the Fortran code reads).

.. contents::
   :local:
   :depth: 2

Overview
========

The recommended workflow is:

1. Create a JSON file with your calculation setup (easier to write/read)
2. Use ``mqc_prep.py`` to convert JSON → ``.mqc`` format
3. Run metalquicha with the ``.mqc`` file

The ``.mqc`` format is a simple section-based format that avoids JSON parsing complexity in Fortran while remaining human-readable.

Generating .mqc Files
=====================

Use the ``mqc_prep.py`` Python script to convert JSON inputs to ``.mqc`` format:

.. code-block:: bash

   python3 mqc_prep.py input.json

This will generate ``input.mqc`` (or ``{title}.mqc`` if a title is specified in the JSON).

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

**Supported methods**:

- ``XTB-GFN1``: GFN1-xTB semi-empirical method (default)
- ``XTB-GFN2``: GFN2-xTB semi-empirical method

**Note**: Basis sets (``basis``, ``aux_basis``) are currently only used for future ab initio methods. XTB methods ignore these fields.

Driver Section
--------------

Specifies the calculation type:

.. code-block:: json

   "driver": "Energy"

**Supported drivers**:

- ``Energy``: Single-point energy calculation
- ``Gradient``: Energy + analytical gradient (if method supports it)

Keywords Section
----------------

SCF Options
^^^^^^^^^^^

.. code-block:: json

   "keywords": {
     "scf": {
       "maxiter": 300,
       "tolerance": 1e-6
     }
   }

- ``maxiter``: Maximum SCF iterations (default: 300)
- ``tolerance``: Convergence tolerance (default: 1e-6)

Fragmentation Options
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: json

   "fragmentation": {
     "method": "MBE",
     "allow_overlapping_fragments": false,
     "level": 2,
     "max_intersection_level": 3,
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
- ``embedding``: Fragment embedding scheme (currently only ``"none"`` supported)
- ``cutoff_method``: How to include fragments (``"distance"``, ``"all"``)
- ``distance_metric``: For distance cutoffs: ``"min"``, ``"max"``, ``"com"`` (center of mass)
- ``cutoffs``: Distance thresholds (in Angstroms) for including dimers, trimers, etc.

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

.mqc File Format
================

The ``.mqc`` format is section-based with ``%section ... end`` delimiters.

Complete Example
----------------

.. code-block:: text

   %schema
   name = mqc-frag
   version = 1.0
   index_base = 0
   units = angstrom
   end

   %model
   method = XTB-GFN1
   basis = cc-pVDZ
   aux_basis = cc-pVDZ-RIFIT
   end

   %driver
   type = Energy
   end

   %system
   log_level = Verbose
   end

   %structure
   charge = 0
   multiplicity = 1
   end

   %geometry
   6

   O   0.000   0.000   0.000
   H   0.000   0.000   0.960
   H   0.757   0.000  -0.240
   end

   %fragments
   nfrag = 2

   %fragment
   charge = 0
   multiplicity = 1
   %indices
   0 1 2
   end
   end

   %fragment
   charge = 0
   multiplicity = 1
   %indices
   3 4 5
   end
   end

   end

   %connectivity
   nbonds = 1

   0 1 1 broken

   nbroken = 1
   end

   %scf
   maxiter = 300
   tolerance = 1e-06
   end

   %fragmentation
   method = MBE
   allow_overlapping_fragments = false
   level = 2
   embedding = none
   cutoff_method = distance
   distance_metric = min

   %cutoffs
   dimer = 10.0
   trimer = 8.0
   end
   end

Section Descriptions
--------------------

``%schema``
^^^^^^^^^^^

Identifies format version and conventions:

- ``name``: Format identifier (``mqc-frag``)
- ``version``: Format version (``1.0``)
- ``index_base``: 0-based (0) or 1-based (1) atom indexing
- ``units``: Coordinate units (``angstrom`` or ``bohr``)

``%model``
^^^^^^^^^^

Quantum chemistry method specification (same as JSON).

``%driver``
^^^^^^^^^^^

Calculation type (``Energy``, ``Gradient``).

``%system``
^^^^^^^^^^^

System configuration, primarily logging level.

``%structure``
^^^^^^^^^^^^^^

Molecular charge and multiplicity.

``%geometry``
^^^^^^^^^^^^^

Atomic coordinates in XYZ format:

.. code-block:: text

   %geometry
   <number_of_atoms>
   <blank line>
   <symbol> <x> <y> <z>
   <symbol> <x> <y> <z>
   ...
   end

``%fragments``
^^^^^^^^^^^^^^

Fragment definitions with nested ``%fragment`` sections:

.. code-block:: text

   %fragments
   nfrag = <number_of_fragments>

   %fragment
   charge = <fragment_charge>
   multiplicity = <fragment_multiplicity>
   %indices
   <atom1> <atom2> <atom3> ...
   end
   end

   end

``%connectivity``
^^^^^^^^^^^^^^^^^

Bond definitions, especially for hydrogen capping:

.. code-block:: text

   %connectivity
   nbonds = <number_of_bonds>

   <atom_i> <atom_j> <bond_order> [broken]

   nbroken = <number_of_broken_bonds>
   end

The ``broken`` keyword marks bonds where hydrogen caps should be added.

``%scf``
^^^^^^^^

SCF convergence parameters (same as JSON).

``%fragmentation``
^^^^^^^^^^^^^^^^^^

Fragmentation parameters with nested ``%cutoffs`` section:

.. code-block:: text

   %fragmentation
   method = MBE
   level = 2
   ...

   %cutoffs
   dimer = 10.0
   trimer = 8.0
   end
   end

Running Calculations
====================

Basic Usage
-----------

.. code-block:: bash

   # Generate .mqc from JSON
   python3 mqc_prep.py input.json

   # Run calculation (serial)
   ./mqc input.mqc

   # Run calculation (parallel with MPI)
   mpirun -np 4 ./mqc input.mqc

   # Run on multiple nodes
   mpirun -np 64 --map-by ppr:32:node ./mqc input.mqc

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

   python3 mqc_prep.py h3o.json
   ./mqc h3o.mqc

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

   python3 mqc_prep.py prism.json
   ./mqc prism.mqc

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

   python3 mqc_prep.py overlapping_gly3.json
   ./mqc overlapping_gly3.mqc

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

1. **Use JSON for input creation**: Easier to write and validate
2. **Keep .mqc files**: Version control friendly, human-readable
3. **Start with small systems**: Test fragmentation schemes on small molecules first
4. **Check JSON output**: Verify fragment energies are reasonable
5. **Use appropriate log levels**: ``verbose`` for debugging, ``info`` for production
6. **Validate results**: Use the validation test suite (see :ref:`validation`)

Troubleshooting
===============

**"Invalid input file extension"**
   Ensure file ends with ``.mqc`` (not ``.json``)

**"Error reading .mqc file"**
   Check for:
   - Missing ``end`` delimiters
   - Mismatched section names
   - Invalid values (e.g., negative atom count)

**Fragment charge/multiplicity mismatch**
   Ensure sum of fragment charges equals molecular charge

**"No fragments generated"**
   Check that ``fragments`` array is not empty and indices are valid

**Hydrogen capping not working**
   Verify bonds are marked as ``is_broken: true`` in JSON or have ``broken`` keyword in ``.mqc``
