# The validation suite

Where the decks live, how they are generated, how to run them, and what a
manifest entry can pin. Read from `AGENTS.md`, which says when.

## Validation Tests

Decks live under `validation/inputs/<hardware>/<engine>/<method>/`, with the
geometries shared at `validation/inputs/sample_inputs/` and reached from a deck
with `../../../`. The CPU suite is **generated** - edit
`tools/cpu_validation/gen_cpu_validation.py` and rerun it; a deck added by hand
under `cpu/mqc/` is deleted by the next regeneration.

Validation categories under `cpu/mqc/` are `bond_orders ccsd ccsd-t charges cosx counterpoise
df-dft df-hf df-uhf dft dh ecp efmo efp fmo fukui gdm gradient hessian interaction_energy
makefp mcscf mp2 neo nmr optimize pcm quao raman rhf ri-ccsd ri-ccsd-t ri-mp2 sapt scs-mp2
rohf soscf sos-mp2 tddft udft uhf` -- worth checking against before assuming
something is untested. That list is `ls validation/inputs/cpu/mqc/`; it went stale twice
because it was maintained by hand, so check the directory rather than this
line if something is missing from it.

Examples, all under `validation/inputs/cpu/tblite/gfn1/`:
- `h3o.json` - Unfragmented hydronium
- `prism.json` - Water prism MBE(2)
- `overlapping_gly3.json` - Glycine tripeptide GMBE(1)
- `w20_isomer.json` - Water 20-mer MBE(2)

The CPU ab initio suite is the bulk of 355 decks under `cpu/mqc/`: RHF across H-Ar
and eight basis sets, UHF, density fitting, MP2, RI-MP2, CCSD(T) and RI-CCSD(T),
and one fragmented case. References are PySCF fed this repository's own basis
JSON, not PySCF's internal tables - those differ in the eighth decimal on Pople
sets, which looks exactly like a bug in whichever code you are checking.

### Running it

`validation/run_validation.py`, and it has to be run **from inside
`validation/`** -- `--exe` defaults to `../build/mqc` and `--validation-dir` to
`.`, so it finds nothing from the repository root. There are two manifests and
the default is the smaller one:

```bash
cd validation
python3 run_validation.py                                       # validation_tests.json, 23 cases
python3 run_validation.py --manifest validation_tests_cpu.json  # 347 cases, the ab initio suite
python3 run_validation.py -t water_sto-3g                       # substring filter on the test name
python3 run_validation.py --mpi --np 4                          # fragmented on 4 ranks
```

**A case that pins only the energy says nothing about a property.** What
`properties` asks for changes no energy, so a deck that requests one and a
manifest entry holding `expected_energy` alone pass whether the property is
right, wrong or absent: Hartree-Fock Raman activities were 33.8 for 45.8 for
a day with the suite green. `expected_output` maps a path into the output to
its value -- `atomic_charges.atoms[*].charge`, `stability.lowest_curvature`
-- and a missing path is a failure; `expected_nmr_isotropic` and
`expected_raman_activities` are the two with their own keys. Every property
deck that writes its result pins it. Those pins are this program's own
numbers, there to catch drift and absence; the comparison with another
program is in the unit test each `reference_note` names. The runner checks
every entry's keys before it runs anything, so a misspelt key is an error
and not a check that never happens. **The quasi-atomic bonding analysis
writes nothing to the output file**, so its seven decks still pin the energy
alone.

**A build without tblite skips most of the default manifest, and skips are not
failures.** `-t h3o` on such a build runs zero tests and prints a green summary.
Check the "Total tests" line, not the colour. `validation_tests_cpu.json` is the
one that exercises a build with the CPU integrals and no tblite.

The runner deletes every `output_*.json` in `validation/` before it starts.
Those are working outputs rather than references, but do not leave anything you
want under that name.

## Running under MPI

- With `run_validation.py --mpi`, unfragmented tests still use `np=1` while fragmented tests use several ranks. To test a parallel Hessian on an unfragmented system, mark the test `"type": "fragmented"`: its displacements are distributed across ranks.
- GitHub Actions has 2 cores, so CI's MPI validation uses `--np 2`. Run it by hand with a higher `--np` locally.
