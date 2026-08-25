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

Where ``<basename>`` is derived from the input filename (e.g., ``water.json`` produces output with key ``water``).

Output files are named ``output_<basename>.json`` and placed in the working directory.

A ``driver: Optimize`` run writes this document for a final single point at the
optimized geometry, and writes the record of the optimization itself -- convergence,
step counts, the final structure and the trajectory -- to a separate
``output_<basename>_optimization.json``. See :doc:`geometry_optimization`.

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
   * - ``gradient``
     - array
     - Nuclear gradient, one ``[x, y, z]`` triple per atom in input order (Hartree/Bohr)
   * - ``gradient_units``
     - string
     - Units the components are written in; always ``hartree/bohr``
   * - ``gradient_norm``
     - float
     - Euclidean norm of the nuclear gradient (Hartree/Bohr)
   * - ``hessian_frobenius_norm``
     - float
     - Frobenius norm of the Hessian matrix
   * - ``dipole``
     - object
     - Dipole moment with ``x``, ``y``, ``z`` components (a.u.) and ``magnitude_debye``

Fukui indices
=============

Written when ``properties.fukui`` asked for them::

    "fukui": {
      "population_scheme": "chelpg",
      "ionisation_potential": 0.402050,
      "electron_affinity": -0.194424,
      "chemical_potential": -0.103813,
      "hardness": 0.298237,
      "electrophilicity": 0.018068,
      "anion_bound": false,
      "atoms": [
        {"atom": 1, "f_plus": -0.1617, "f_minus": 0.7804,
         "f_zero": 0.3093, "dual": -0.9420},
        {"atom": 2, "f_plus": 0.5797, "f_minus": 0.1096,
         "f_zero": 0.3447, "dual": 0.4702}
      ]
    }

Per atom rather than reduced to a most-reactive site, because ranking sites is
what the caller is doing and which index to rank on depends on the reaction:
``f_plus`` for an incoming nucleophile, ``f_minus`` for an electrophile, and
``dual`` to tell the two kinds of site apart by sign.

``atom`` is the index in the input geometry, and there is no element symbol
because the output payload does not carry them -- a consumer that read the
geometry already knows which atom is which.

The ions are run with the deck's own method: Hartree-Fock for ``hf``, and
unrestricted Kohn-Sham with the neutral's functional for ``dft``. The energies
here are therefore differences of that method's total energies and are not
comparable with a run that used another one -- nothing in this payload records
which, because the deck that produced it does.

.. warning::

   **A negative ``f_plus`` or ``f_minus`` is spurious** -- the example above
   has one. The exact Fukui function cannot be negative, so a negative
   condensed value is an artefact of partitioning a continuous density onto
   atoms, not a site that repels charge.

   A consumer should rank on these, not quote them, and should treat a small
   negative as "unreactive in that channel". ``f_zero`` and ``dual`` are
   derived, so they inherit the artefact from whichever index carried it. A
   script that sorts ascending on ``f_minus`` to find the least reactive site
   will find the most artefactual one instead.

**Check ``anion_bound`` before using ``f_plus``.** When it is false the anion
came out above the neutral, nothing bound the added electron, and that column
describes whatever orbital the basis had left over. Nothing about the values
themselves says so, and a script sorting on ``f_plus`` would rank sites
confidently off a fiction. The same flag makes ``electron_affinity`` and
anything derived from it -- the chemical potential, the hardness, the
electrophilicity -- unreliable in the same breath.

Both indices sum to one over the molecule by construction, so a consumer can
assert it as a cheap check that it read the right section.

.. _unconverged-fragments:

Fragments that did not converge
===============================

A fragmented run with ``allow_crap_scf`` finishes rather than stopping at the
first failure, and the total it reports is built partly from fragments that
never converged. Those fragments are listed::

    "unconverged": {
      "count": 2,
      "energy_at_stake": -0.0031,
      "largest_contribution": 0.0024,
      "fragments": [
        {"id": 12, "level": 2, "monomers": [3, 7], "delta_energy": -0.0024},
        {"id": 45, "level": 1, "monomers": [9], "delta_energy": -0.0007}
      ],
      "monomers_involved": [
        {"monomer": 9, "fragments": 2},
        {"monomer": 3, "fragments": 1},
        {"monomer": 7, "fragments": 1}
      ]
    }

``monomers`` is what makes this actionable. An identifier on its own cannot be
re-run: a dimer is only reconstructible if you know which two monomers it was.
With them, a follow-up job that revisits exactly the misbehaving fragments --
tighter thresholds, a different guess, a smaller trust radius -- writes itself
from the output, which is what the Python interface is for.

``energy_at_stake`` answers the question the count does not: whether to care.
It is the net contribution of the failed fragments to the total, so a run that
lost a hundred fragments which screening had already made negligible is a run to
accept rather than repeat. It is signed, because the failures are a subset of a
sum with cancellation in it; ``largest_contribution`` sits beside it because a
net near zero can still be two large terms that happen to oppose.

``monomers_involved`` answers the other one: why. A monomer with a wrecked
geometry, a mis-assigned charge or an accidental radical drags down every
fragment it belongs to, so failures cluster on their cause. Four hundred failed
dimers sharing one monomer is one problem, not four hundred, and the list is
ordered by how many failures each monomer appears in so that the culprit is the
first row rather than one of four hundred that look alike.

Three things about how it is written.

**The section appears whenever the method reports convergence at all**, with a
``count`` of zero when nothing failed. That is deliberate: an absent section
means "this method never said", which is a different claim from "everything
converged", and a consumer should be able to tell them apart. Methods that do
not report convergence omit the section entirely rather than claiming success
for every fragment.

**Fragments whose status is unknown are not listed.** They are not failures, and
including them would fill the list with the entire calculation in exactly the
runs where it matters.

**The list is not truncated.** Ten thousand unconverged fragments is a large
list and also the situation this exists for; keeping the first hundred would
produce a follow-up job that looked complete and was not. The log still names
only the first ten, because that is the right length for something a person
reads.

Currently written for the many-body expansion. The overlapping-fragment path
does not record per-fragment convergence, so the section is absent there rather
than empty.


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
       "gradient_units": "hartree/bohr",
       "gradient": [
         [0.0, 0.0, -0.067617391054],
         [0.0, -0.016791268486, 0.033808695649],
         [0.0, 0.016791268169, 0.033808695405]
       ],
       "hessian_frobenius_norm": 1.234567890123
     }
   }

Fields
------

All fields are optional depending on the calculation type:

- ``total_energy``: Present for Energy, Gradient, and Hessian calculations
- ``dipole``: Present when dipole moments are computed
- ``gradient``, ``gradient_units``, ``gradient_norm``: Present for Gradient and
  Hessian calculations
- ``hessian_frobenius_norm``: Present for Hessian calculations

The components are what a validation case compares. A norm cannot tell a
correct gradient from one with a sign error, a swapped pair of atoms, or a
wrong component that happens to preserve the magnitude, so
``validation/run_validation.py`` reads ``gradient`` and compares element by
element against ``expected_gradient`` in the manifest -- and reports the
translational residual, :math:`\sum_A \partial E/\partial R_A`, which should
vanish, while it is there.

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
