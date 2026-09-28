# EE-MBE PySCF reference

`eembe_pyscf.py` is an independent PySCF reimplementation of metalquicha's
electrostatically-embedded many-body expansion (`keywords.fragmentation.method:
"ee-mbe"`), built by reading the Fortran rather than by fitting numbers. The
full derivation of the semantics -- which routines say what, and why -- is in
the script's own module docstring; read that before trusting or extending
this file.

It exists for two things:

1. An independent check of the EE-MBE Hartree-Fock energy against
   `build/mqc`, for the water-trimer deck named below.
2. An independent MP2 / RI-MP2 / SCS-MP2 reference for the *same* system:
   correlation added to each embedded monomer and dimer on top of its
   embedded Hartree-Fock, the charges staying Hartree-Fock, as
   `mqc_czt_fragment_solver.f90` does it. `test/test_mqc_fmo_mp2.f90` holds
   these numbers.

## Running it

```
<python-with-pyscf-2.14> tools/fmo_validation/eembe_pyscf.py --method hf
<python-with-pyscf-2.14> tools/fmo_validation/eembe_pyscf.py --method mp2 [--all-electron]
<python-with-pyscf-2.14> tools/fmo_validation/eembe_pyscf.py --method ri-mp2 [--all-electron] [--aux-basis NAME]
<python-with-pyscf-2.14> tools/fmo_validation/eembe_pyscf.py --method scs-mp2 [--all-electron]
```

It reads `validation/inputs/sample_inputs/w3.xyz` and reproduces exactly
`validation/inputs/cpu/mqc/fmo/eembe_water3.json`: three waters (file order,
one fragment per water), 6-31G, EE-MBE level 2, restricted closed-shell. It
takes no other arguments for the system itself -- this script checks one
system, not a family of them. Basis and auxiliary-basis data are read out of
this repository's own `basis_sets/*.json`, through the same loader
(`bse_to_pyscf`/`molecule_form`) `tools/cpu_validation/gen_cpu_validation.py`
uses, imported from there rather than reimplemented, so the two never
silently disagree on how a Pople set's coefficients round.

## Hartree-Fock: agreement with the binary

Run identically:

```
cp build/mqc <scratch>/mqc_ref_hf
cd <scratch>; <copy of validation/inputs/cpu/mqc/fmo/eembe_water3.json, with
  ../../../sample_inputs/w3.xyz resolving to a copy of the real file>
OMP_NUM_THREADS=4 <scratch>/mqc_ref_hf eembe_water3.json
```

| quantity | mqc (binary) | this script | difference |
|---|---:|---:|---:|
| total_energy | -227.970457333730 | -227.970457333698 | 3.2e-11 |
| monomer_sum | -227.923794994429 | -227.923794995408 | 9.8e-10 |
| pair_sum | -0.046662339301 | -0.046662338290 | 1.0e-9 |
| pair 1-2 dE | -0.028243500848 | -0.028243508700 | 7.9e-9 |
| pair 2-3 dE | -0.015613045517 | -0.015613036390 | 9.1e-9 |
| pair 1-3 dE | -0.002805792937 | -0.002805793199 | 2.6e-10 |

All within the ~1e-8 Eh target; the total agrees to 3e-11. The per-term
residuals track the *order of magnitude* mqc's own inner-SCF tolerances
(`fmo_scf_energy_tol = 1e-9`, `fmo_scf_density_tol = 1e-7`) can leave behind
against this script's much tighter `conv_tol = 1e-12` -- consistent with
both codes converging to the *same* fixed point, one more tightly than the
other, rather than computing two different things. The outer SCC loop
converges in 5 passes in both (traced to 12 decimals in the script's `--method
hf` output), which is itself part of the agreement: a different embedding
field, sign convention, or nucleus/point-charge energy term would not
reproduce mqc's own iteration trajectory this closely, only its endpoint at
best.

## MP2 / RI-MP2 / SCS-MP2 (Hartree, 12 decimals)

The embedding charges and the SCC loop are Hartree-Fock throughout (the
converged Mulliken charges and the dimers' embedding fields are identical to
the Hartree-Fock run above); MP2 correlation is added on top of each
already-embedded monomer's and dimer's own converged HF orbitals, and the
EE-MBE total is assembled from the correlated `E'_I = E_HF,embedded + E_corr`
values by the same `sum_I E'_I + sum_pairs(E'_IJ - E'_I - E'_J)` formula as
Hartree-Fock. Frozen core is `core_orbital_count` per group (1 for a monomer,
2 for a dimer -- each water's own O 1s); RI-MP2's auxiliary basis is
`cuest_scf_settings_t`'s own default, `def2-universal-jkfit`.

| method | core | total |
|---|---|---:|
| MP2 | frozen | -228.354417062659 |
| MP2 | all-electron | -228.357563922013 |
| RI-MP2 (def2-universal-jkfit) | frozen | -228.354367521908 |
| RI-MP2 (def2-universal-jkfit) | all-electron | -228.357514271563 |
| SCS-MP2 | frozen | -228.353055036906 |
| SCS-MP2 | all-electron | -228.356012217249 |

RI-MP2 disagrees with canonical MP2 by ~5e-5 Eh here -- larger than a proper
`-rifit` correlation basis would give. That is expected and is not a bug in
this script: `correlation_aux_basis` in `mqc_czt_bridge.f90` itself warns
that a JKFIT set is "not a correlation-fitting (RIFIT) set" whenever nothing
better is named, and uses it anyway, which is exactly the default this script
reproduces.

Per-monomer and per-pair breakdowns (correlation energy alongside each) are
printed by the script itself; they are not duplicated here since they are
regenerated, not hand-copied, by every run.

## What was not, and could not be, checked

MP2/RI-MP2/SCS-MP2 could not be checked against `build/mqc`:
`mqc_fragment_capabilities.f90`'s `fragment_capabilities` only lets
`METHOD_TYPE_HF`/`METHOD_TYPE_DFT` run under `FRAGMENT_SCHEME_FMO`/
`FRAGMENT_SCHEME_EE_MBE` as of this branch, so a deck asking for MP2 under
`ee-mbe` is refused before it reaches the backend. The numbers above instead
follow `mqc_czt_fragment_solver.f90::solve_fragment_method`/
`fragment_correlation` -- the method-agnostic solver EFMO's own MP2 already
uses, and the one a FMO/EE-MBE MP2 branch is expected to call unchanged --
read directly rather than approximated, so once that branch lands a real
`ee-mbe` MP2 deck on this same system should reproduce them to about the same
~1e-8 to 1e-7 Eh mqc's SCF tolerances leave for Hartree-Fock above.
