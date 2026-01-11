.. _json_output:

==================
JSON Output Format
==================

Metalquicha produces JSON output files containing calculation results. The output format is standardized using the json-fortran library with scientific notation (ES format) for all real numbers.

.. contents::
   :local:
   :depth: 2

Overview
========

All JSON output files follow a common structure:

.. code-block:: json

   {
     "<basename>": {
       ...calculation results...
     }
   }

Where ``<basename>`` is derived from the input filename (e.g., ``water.mqc`` produces output with key ``water``).

Output files are named ``output_<basename>.json`` and placed in the working directory.

Common Fields
=============

Several fields appear across multiple output types:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``total_energy``
     - float
     - Total electronic energy in Hartree
   * - ``gradient_norm``
     - float
     - Euclidean norm of the nuclear gradient (Hartree/Bohr)
   * - ``hessian_frobenius_norm``
     - float
     - Frobenius norm of the Hessian matrix
   * - ``dipole``
     - object
     - Dipole moment with ``x``, ``y``, ``z`` components (a.u.) and ``magnitude_debye``

Unfragmented Calculations
=========================

For single-point calculations without fragmentation (``driver: Energy``, ``Gradient``, or ``Hessian``).

Schema
------

.. code-block:: json

   {
     "<basename>": {
       "total_energy": -76.123456789012,
       "dipole": {
         "x": 0.0,
         "y": 0.0,
         "z": 0.756123456789,
         "magnitude_debye": 1.923456789012
       },
       "gradient_norm": 0.012345678901,
       "hessian_frobenius_norm": 1.234567890123
     }
   }

Fields
------

All fields are optional depending on the calculation type:

- ``total_energy``: Present for Energy, Gradient, and Hessian calculations
- ``dipole``: Present when dipole moments are computed
- ``gradient_norm``: Present for Gradient and Hessian calculations
- ``hessian_frobenius_norm``: Present for Hessian calculations

MBE Detailed Breakdown
======================

For Many-Body Expansion (MBE) calculations, provides detailed fragment energies organized by n-mer level.

Schema
------

.. code-block:: json

   {
     "<basename>": {
       "total_energy": -152.234567890123,
       "levels": [
         {
           "frag_level": 1,
           "name": "monomers",
           "count": 10,
           "total_energy": -150.123456789012,
           "fragments": [
             {
               "indices": [1],
               "energy": -15.012345678901,
               "distance": 0.0
             },
             {
               "indices": [2],
               "energy": -15.012345678901,
               "distance": 0.0
             }
           ]
         },
         {
           "frag_level": 2,
           "name": "dimers",
           "count": 45,
           "total_energy": -2.111111111111,
           "fragments": [
             {
               "indices": [1, 2],
               "energy": -30.034567890123,
               "distance": 3.456789012345,
               "delta_energy": -0.009876543210
             }
           ]
         }
       ],
       "dipole": { ... },
       "gradient_norm": 0.012345678901,
       "hessian_frobenius_norm": 1.234567890123
     }
   }

Fields
------

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``levels``
     - array
     - Array of n-mer levels (monomers, dimers, trimers, etc.)
   * - ``levels[].frag_level``
     - integer
     - The n-mer order (1 = monomer, 2 = dimer, etc.)
   * - ``levels[].name``
     - string
     - Human-readable name (monomers, dimers, ..., decamers, or "N-mers")
   * - ``levels[].count``
     - integer
     - Number of fragments at this level
   * - ``levels[].total_energy``
     - float
     - Sum of all fragment energies at this level
   * - ``levels[].fragments``
     - array
     - Individual fragment data
   * - ``fragments[].indices``
     - array[int]
     - Monomer indices that compose this fragment
   * - ``fragments[].energy``
     - float
     - Raw fragment energy (Hartree)
   * - ``fragments[].distance``
     - float
     - Characteristic distance for n-mers (n > 1)
   * - ``fragments[].delta_energy``
     - float
     - MBE correction energy for n-mers (n > 1)

GMBE Output
===========

For Generalized Many-Body Expansion with overlapping fragments.

Schema
------

.. code-block:: json

   {
     "<basename>": {
       "total_energy": -152.234567890123,
       "monomers": {
         "count": 3,
         "fragments": [
           {
             "index": 1,
             "energy": -50.123456789012
           },
           {
             "index": 2,
             "energy": -50.123456789012
           }
         ]
       },
       "intersections": {
         "total_count": 3,
         "levels": [
           {
             "level": 2,
             "count": 3,
             "fragments": [
               {
                 "indices": [1, 2],
                 "energy": -25.012345678901
               }
             ]
           }
         ]
       }
     }
   }

Fields
------

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``monomers.count``
     - integer
     - Number of primary fragments
   * - ``monomers.fragments``
     - array
     - Array of monomer data with ``index`` and ``energy``
   * - ``intersections``
     - object
     - Present only when fragments overlap
   * - ``intersections.total_count``
     - integer
     - Total number of intersection terms
   * - ``intersections.levels``
     - array
     - Intersections grouped by overlap level

GMBE PIE Output
===============

For Pairwise Interaction Energy decomposition within GMBE.

Schema
------

.. code-block:: json

   {
     "<basename>": {
       "total_energy": -152.234567890123,
       "gradient_norm": 0.012345678901,
       "hessian_frobenius_norm": 1.234567890123,
       "pie_terms": {
         "count": 100,
         "terms": [
           {
             "atom_indices": [0, 1, 2],
             "coefficient": 1,
             "energy": -50.123456789012,
             "weighted_energy": -50.123456789012
           },
           {
             "atom_indices": [0, 1],
             "coefficient": -1,
             "energy": -30.012345678901,
             "weighted_energy": 30.012345678901
           }
         ]
       }
     }
   }

Fields
------

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field
     - Type
     - Description
   * - ``pie_terms.count``
     - integer
     - Number of non-zero PIE terms
   * - ``pie_terms.terms``
     - array
     - Array of PIE term data
   * - ``terms[].atom_indices``
     - array[int]
     - Atom indices in this term (0-indexed)
   * - ``terms[].coefficient``
     - integer
     - Inclusion-exclusion coefficient (+1 or -1)
   * - ``terms[].energy``
     - float
     - Raw subsystem energy (Hartree)
   * - ``terms[].weighted_energy``
     - float
     - ``coefficient * energy`` contribution to total

Vibrational Analysis
====================

For Hessian calculations with thermochemistry analysis.

Schema
------

.. code-block:: json

   {
     "<basename>": {
       "total_energy": -76.123456789012,
       "dipole": { ... },
       "gradient_norm": 0.000012345678,
       "hessian_frobenius_norm": 1.234567890123,
       "vibrational_analysis": {
         "n_modes": 9,
         "frequencies_cm1": [1650.123, 3800.456, ...],
         "reduced_masses_amu": [1.023, 1.045, ...],
         "force_constants_mdyne_ang": [5.123, 8.456, ...],
         "ir_intensities_km_mol": [45.67, 89.01, ...]
       },
       "thermochemistry": {
         "temperature_K": 298.15,
         "pressure_atm": 1.0,
         "molecular_mass_amu": 18.015,
         "symmetry_number": 2,
         "spin_multiplicity": 1,
         "is_linear": false,
         "n_real_frequencies": 3,
         "n_imaginary_frequencies": 0,
         "moments_of_inertia_amu_ang2": {
           "Ia": 0.612,
           "Ib": 1.156,
           "Ic": 1.768
         },
         "rotational_constants_GHz": {
           "A": 825.234,
           "B": 437.123,
           "C": 285.789
         },
         "partition_functions": {
           "translational": 3.456e+06,
           "rotational": 34.567,
           "vibrational": 1.000
         },
         "contributions": {
           "translational": {
             "energy_hartree": 0.001416,
             "entropy_cal_mol_K": 34.608,
             "Cv_cal_mol_K": 2.981
           },
           "rotational": { ... },
           "vibrational": { ... },
           "electronic": { ... }
         },
         "contribution_table": {
           "VIB": {
             "H_cal_mol": 0.0,
             "Cp_cal_mol_K": 0.0,
             "S_cal_mol_K": 0.0,
             "S_J_mol_K": 0.0
           },
           "ROT": { ... },
           "INT": { ... },
           "TR": { ... },
           "TOT": { ... }
         },
         "zero_point_energy_hartree": 0.021234,
         "zero_point_energy_kcal_mol": 13.325,
         "thermal_corrections_hartree": {
           "to_energy": 0.024567,
           "to_enthalpy": 0.025511,
           "to_gibbs": 0.002345
         },
         "total_energies_hartree": {
           "electronic": -76.123456789012,
           "electronic_plus_zpe": -76.102222789012,
           "electronic_plus_thermal_E": -76.098889789012,
           "electronic_plus_thermal_H": -76.097945789012,
           "electronic_plus_thermal_G": -76.121111789012
         }
       }
     }
   }

Key Thermochemistry Fields
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Field
     - Description
   * - ``zero_point_energy_hartree``
     - Zero-point vibrational energy (ZPE)
   * - ``thermal_corrections_hartree.to_gibbs``
     - Thermal correction to Gibbs free energy
   * - ``total_energies_hartree.electronic_plus_thermal_G``
     - Total Gibbs free energy (E + thermal G correction)

The ``contribution_table`` provides a summary matching common quantum chemistry output formats:

- **VIB**: Vibrational contributions
- **ROT**: Rotational contributions
- **INT**: Internal (VIB + ROT)
- **TR**: Translational contributions
- **TOT**: Total thermodynamic properties

Number Format
=============

All real numbers are output in scientific notation (ES format) using json-fortran's default machine precision. The format is controlled by the ``JSON_REAL_FORMAT`` parameter in ``mqc_program_limits.f90``.
