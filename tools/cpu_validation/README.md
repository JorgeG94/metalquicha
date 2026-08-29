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
worth about 2e-8 Hartree on STO-3G water -- an order of magnitude outside the
1e-9 tolerance this suite uses, and enough to waste an afternoon. Using
one table for both sides means a disagreement is a disagreement in the
integrals or the SCF, which is what is being tested.

With the tables shared this way, the Hartree-Fock and DFT cases agree to 8e-12
Hartree or better. The correlated methods do not, and not because they are worse:
a correlation energy is not stationary in the density, so where a variational
energy hides a loose SCF quadratically, MP2 and coupled cluster expose it
linearly. At an SCF tolerance of 1e-8 that put MP2 4e-9 from PySCF and CCSD 6e-9
-- entirely the SCF stopping point. The decks converge to 1e-12 instead, which
brings those to 9e-11 and 5e-10.

That is where the 1e-9 manifest tolerance gets its headroom: two orders above
the worst case, still loose enough to survive a change of compiler or BLAS, and
tight enough that a real integral or contraction bug cannot hide under it. The
suite ran at 1e-6 until the correlated cases were measured; at that level an
error four orders larger than the numerical noise passed silently.

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
| nineteen ECP cases | effective core potentials; see below |

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

**Transition metals.** They need ECPs, which the all-electron sweeps above do
not use. They have their own section below.

## Effective core potentials

`ECP_CASES` and the three tables after it. Every deck names `model.ecp:
def2-ecp`, the one potential file shipped, and lives under
`inputs/cpu/mqc/ecp/`.

The reference is PySCF fed **our** `basis_sets/def2-ecp.json`, converted by
`bse_ecp_to_pyscf`, for the same reason the orbital bases are: a potential our
own reader had mangled would otherwise reach PySCF intact, and the disagreement
would be chased through the integrals rather than through the file. The
converter reproduces PySCF's own def2 tables element for element where PySCF
ships one -- it does not for the lanthanides, which is why Yb has to be fed
ours.

The heavy elements are otherwise interchangeable to look at, so what each is
there for:

| case | what it covers |
| --- | --- |
| Sr | one ECP centre and nothing else -- the atomic term alone |
| I2 | two potentials at bonding distance: the two-centre term, which vanishes for an atom and which HI cannot check because only one end carries a potential |
| HI in def2-SVP and def2-TZVP | 31 functions against 56, `s p d` against `s p d f` -- the TZVP entry is where an f orbital shell meets the projectors |
| RbH | the lightest element def2-ECP covers |
| CsH | the 46-electron core |
| AuH, PbH4, HgCl2 | the 60-electron core |
| HgCl2, SnH4, TeH2 | a potential among all-electron neighbours |
| iodine atom | unrestricted, doublet by electron count once the core is gone |
| water | an ECP named where **no** element carries one |

Three core sizes appear on purpose: def2-ECP removes 28 electrons from Rb-Lu, 46
from Cs-La and 60 from Hf-Rn, and those are separate tabulations. A count taken
from the wrong one still produces a converging SCF, with the wrong number of
electrons in it and nothing in the output to say so.

The water case is the one worth understanding. def2-ECP carries nothing below
krypton, so naming it on water must reproduce the all-electron def2-SVP number
in the RHF sweep *exactly* -- and it does, to all twelve digits the manifest
carries. It is the case that fails if naming a potential perturbs a molecule
that has none.

**Frozen cores are written down, never counted.** The automatic count works from
the atomic number, and an ECP has already removed a core, so counting would
freeze *valence* orbitals instead. metalquicha refuses that outright, which is
why the two correlated entries name `n_frozen_core` explicitly. Four is iodine's
4s and 4p: the potential took everything through 3d, leaving 4s 4p 4d 5s 5p, and
4d is shallow enough to want correlating.

**Cost, and the two cases that are not here.** An ECP integral is steep in the
*orbital* angular momentum, not only in the potential's. Seventeen of these run
in under 1.5 s; HI/def2-TZVP takes 12 s, which is the price of the only f-shell
case and worth paying once.

Two candidates were too expensive to keep. Xenon is the bare-atom case in
principle, but def2-SVP hands xenon the same f-containing table def2-TZVP does
-- 50 functions, 11 s -- so strontium takes that slot at 23 functions and 0.4 s,
over the same 28-electron core. Ytterbium is the only closed-shell lanthanide
and so the only practical **l = 5 local channel**, at 81 functions with a g
shell and half a minute per energy. That one moved to
`test/test_mqc_ecp_matrix.f90`, which checks the same angular path against the
same PySCF reference in 0.9 s. Neither the manifest nor CI is the right place
for a 36-second case when a one-electron integral answers the question.

**What is not here.** Nuclear derivatives, because they are refused -- libfint
carries the energy integrals only, not libcint's `nr_ecp_deriv.c` -- along with
MCSCF, xTB and the GPU backend. Those refusals are unit tests
(`test/test_mqc_json_schema.f90` and `test/test_mqc_ecp_refusals.f90`) rather
than manifest entries, because `run_validation.py` compares energies and has no
way to express a deck that is supposed to fail.

**cc-pVTZ on the second row.** Correct but slow: 30 s for SiH4, 20 s for PH3,
18 s for AlH3, against a 4 s budget per case. Ar keeps one second-row general
contraction with f in the suite, and def2-TZVP covers Al-Ar's f shells cheaply.
