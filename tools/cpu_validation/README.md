# CPU validation table

`validation/validation_tests_cpu.json` and the `validation/inputs/cpu_*.json`
decks it points at are generated, not written. One script produces all of them:

```bash
.venv/bin/python tools/cpu_validation/gen_cpu_validation.py
```

It needs PySCF, which is a development dependency only -- the build does not
use it. Run it from anywhere; paths are resolved relative to the repository.
`--dry-run` prints the energies and writes nothing.

Then check the code against the table:

```bash
cd validation
python3 run_validation.py --manifest validation_tests_cpu.json --exe ../build/mqc
```

To change the coverage, edit `MOLECULES` or `SWEEPS` at the top of the script
and rerun it. Do not edit the manifest or the `cpu_*` decks: the script
overwrites the manifest and deletes any `cpu_*.json` deck it did not just
write.

## Where the numbers come from

Every `expected_energy` is a PySCF RHF energy, never a number metalquicha
printed -- a manifest of our own output would detect change, not error.

PySCF is fed the **same basis JSON that metalquicha reads**, out of
`basis_sets/`, converted by `bse_to_pyscf`. That matters: PySCF's own basis
library is rounded differently from the Basis Set Exchange tables, which is
worth about 2e-8 Hartree on STO-3G water -- indistinguishable from a real bug
at the 1e-6 tolerance this suite uses, and enough to waste an afternoon. Using
one table for both sides means a disagreement is a disagreement in the
integrals or the SCF, which is what is being tested.

With the tables shared this way, all 113 cases agree to 8e-12 Hartree or
better, which is where the 1e-6 manifest tolerance gets its headroom: it is
loose enough to survive a change of compiler or BLAS and still six orders of
magnitude tighter than any real integral or contraction bug.

The shell-construction rule in `bse_to_pyscf` mirrors
`src/basis/mqc_json_basis_reader.f90`: one shell per coefficient column, with
each column taking its own angular momentum only when there are as many
angular momenta as columns (the SP case). Reading the same file the same way
is what makes a dropped coefficient column show up as an energy difference.

Geometries live in `MOLECULES` as Python literals and are written out to
`validation/inputs/sample_inputs/*.xyz`, so PySCF and the decks cannot drift
apart. Water is the exception: it reuses the pre-existing `sample_inputs/w1.xyz`
that `hf_water_*.json` already points at.

## What the sweeps cover

| sweep | what it is there for |
| --- | --- |
| STO-3G, all of H-Ar | minimal, SP-contracted shells, one closed-shell species per element |
| 3-21G, 6-31G, all of H-Ar | segmented split valence, SP-contracted |
| cc-pVDZ, all of H-Ar | general contraction (several coefficient columns on one l), d shells |
| cc-pVTZ, subset | general contraction carried to f |
| aug-cc-pVDZ, subset | diffuse functions |
| def2-SVP, subset | Karlsruhe contraction pattern, d shells |
| def2-TZVP, all of H-Ar | Karlsruhe, f shells on B-Ne and Al-Ar |
| 6-31G\*, 6-31G\*\*, all of H-Ar | the Cartesian convention |
| four density-fitted cases | the three- and two-centre path, both conventions |

Closed-shell singlets and energies only: the libcint CPU path is RHF and has no
gradients, so anything else fails by design.

## The angular form

Every other sweep here is spherical. The polarised Pople sets are not: BSE marks
their polarisation shells `function_type: "gto_cartesian"` on every element from
lithium up, which is six d functions rather than five.

The generator does not assume either way. `molecule_form` reads the same
`function_type` entries `mqc_json_basis_reader.f90` reads, by the same rules --
a shell at or below p is undecided whatever it is labelled, and above p anything
not spelled exactly `gto_cartesian` is spherical -- and sets `mol.cart` from the
answer. Hardcoding it is the thing these sweeps exist to catch. It was hardcoded
spherical here, and that is right for most of the table and silently wrong for
these two sets, so the reference agreed with the bug: water came out
-76.008427 in 18 functions where 6-31G\* means -76.009809 in 19, both converging
tidily, neither reporting anything.

H2 is worth keeping in both sweeps for a reason that is not obvious. 6-31G\*
leaves hydrogen unpolarised, so H2 is the one case in either sweep with no shell
above p anywhere in it -- the undecided form, which has to be a third state
rather than a boolean, or water would be refused for disagreeing with its own
hydrogens.

## Density fitting

Four cases rather than a sweep: DF is a different code path, not a wider range
of the same one, and one case per convention plus a second of each is enough to
tell a broken three-centre integral from a broken fitting solve.

The two Cartesian entries fit a basis with itself, which is a poor fitting set
and a deliberate choice. All three centres of a fitting integral are built in
one angular form, every JKFIT set shipped here is spherical, and so a Cartesian
orbital basis has no shipped auxiliary it can legally pair with. metalquicha
refuses that pairing by name; `pyscf_rhf` refuses it too, rather than quietly
generating a number for a calculation the code will not run.

**Transition metals.** They need ECPs, which are unimplemented.

**cc-pVTZ on the second row.** Correct but slow: 30 s for SiH4, 20 s for PH3,
18 s for AlH3, against a 4 s budget per case. Ar keeps one second-row general
contraction with f in the suite, and def2-TZVP covers Al-Ar's f shells cheaply.
