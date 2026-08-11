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

**Basis sets.** ``basis`` is required by every ab initio method and ignored by the
semi-empirical ones, which carry their own parameters. ``aux_basis`` is the
auxiliary set for a density-fitted *reference*; the auxiliary set for a fitted
*correlation* treatment is ``keywords.correlation.aux_basis``, which is separate
because the two can legitimately differ.

**Not yet reachable**: ``DFT``, ``MCSCF`` and the F12 variants parse but have no
CPU implementation. ``RI-CCSD`` needs a restricted reference; unrestricted
coupled cluster is refused rather than quietly run restricted.

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
       "tolerance": 1e-6,
       "guess": "auto",
       "unrestricted": false,
       "density_fitting": false,
       "allow_crap_scf": false
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
- ``allow_crap_scf``: Accept a non-converged SCF instead of stopping (default:
  false). Off because the energy of an SCF that ran out of iterations has the
  right magnitude and nothing downstream can tell.

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
     "aux_basis": "cc-pVDZ-RIFIT",
     "scs": false
   }

- ``freeze_core``: Leave core orbitals uncorrelated (default: true)
- ``n_frozen_core``: How many to freeze (default: -1, counted from the elements)
- ``density_fitting``: Fit the correlation integrals (default: false, or true if
  the method name carries an ``RI-``/``DF-`` prefix -- an explicit keyword wins)
- ``aux_basis``: Auxiliary set for the fit; falls back to ``model.aux_basis``
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
   * - ``%model`` / ``method``, ``basis``, ``aux_basis``, ``functional``
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
