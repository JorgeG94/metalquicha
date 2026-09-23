# Per-pair interaction energies for protein–ligand analysis in metalquicha

Measured on `origin/main` @ `31e4e0bd73`, build `cmake --preset default -DMQC_ENABLE_TBLITE=OFF`,
gfortran, `OMP_NUM_THREADS=1`. Every number below was produced by the runs described here.

**What works today:** plain MBE at level 2. You declare each residue as a fragment and the
ligand as one more, and you read the ligand's pair rows out of a CSV that is written by
default. Each ligand–residue number is *exactly* the supermolecular interaction energy of
that capped residue with the ligand — verified to 1e-13 Hartree.

**What does not work:** `fmo` and `ee-mbe` refuse a partition that cuts the backbone.
`efmo` does not refuse it, does not cap it, and returns `NaN` with a success banner.
Do not use them for a protein.

---

## 1. The deck

```json
{
  "schema": { "name": "pair-analysis", "version": "1.0" },
  "molecules": [{
    "xyz": "complex.xyz",
    "molecular_charge": 0,
    "molecular_multiplicity": 1,

    "fragments": [
      [0,1,2,3,4,5,6,7],            // residue 1
      [8,9,10,11,12,13,14],         // residue 2
      [15,16,17,18,19,20,21,22,23], // residue 3
      [24,25,26]                    // THE LIGAND -- its own fragment, listed last
    ],
    "fragment_charges":         [0, 0, 0, 0],
    "fragment_multiplicities":  [1, 1, 1, 1],

    "connectivity": [
      [2, 8, 1],                    // the C(=O)-N peptide bond cut between residues 1 and 2
      [10, 15, 1]                   // and between 2 and 3
    ]
  }],

  "model": { "method": "hf", "basis": "6-31g" },

  "keywords": {
    "scf": { "maxiter": 200, "tolerance": 1e-10 },
    "fragmentation": { "method": "mbe", "level": 2 }
  },

  "driver": "Energy",
  "system": { "logger": { "level": "Info" } }
}
```

Run it: `./mqc complex.json`. Atom indices are **0-based**.

### The three things that must be right

1. **Every atom belongs to exactly one fragment.** An atom in no fragment is refused
   (`"atom N belongs to no monomer"`); the ligand's atoms must be in the ligand's list and
   nowhere else.
2. **Every bond that crosses a fragment boundary must appear in `connectivity`.** That is
   what triggers hydrogen capping. You do *not* need to list bonds inside a fragment — for a
   backbone cut one residue per residue, listing just the N−1 peptide bonds is enough and is
   what the shipped `gly10` deck does. A bond is marked broken automatically from the
   fragment lists; the third number in each triple is a bond order and is not used here.
3. **A ligand is not bonded to anything**, so it contributes no `connectivity` entries and
   gets no caps. That is the whole reason its rows are the trustworthy ones.

### If you forget `connectivity`

You get **no warning**. The audit that would catch it is gated on a field the JSON path never
fills, so it never fires from a deck. What happens instead is luck:

- if a fragment ends with an odd electron count, the run dies with
  `"UHF: an electron count and multiplicity that cannot be paired -- their parities disagree"`;
- if every fragment happens to come out even (two cuts per fragment, a ring, a double bond),
  **it runs to completion and prints a wrong number.** Measured: a 3-fragment cut of the same
  complex gave −772.07972 uncapped against −772.08161 capped, with zero diagnostics.

Declare the bonds.

---

## 2. Where the output lands, and how to read it

Two files are written into the **current working directory**, named from the deck's basename:

| file | content |
|---|---|
| `output_<basename>.json` | totals, per-level sums, thermochemistry |
| `output_<basename>_fragments.csv` | **one row per term — this is the pair table** |

The CSV is written by default. It is controlled by `"system": {"fragment_breakdown": "csv"}`,
whose values are `csv` (default), `json` (same data folded into the JSON instead) or `none`.
At `"system": {"logger": {"level": "verbose"}}` the same breakdown is also printed to stdout.

### Columns

```
frag_index,level,m1,m2,energy,delta_energy,distance,scf,homo,lumo,gap_ev,charge,mult
```

| column | meaning | units |
|---|---|---|
| `frag_index` | row counter | — |
| `level` | how many fragments in this term: 1 = monomer, 2 = pair | — |
| `m1 … m<level>` | which fragments, **1-based**, in the order of the deck's `fragments` list, zero-padded | — |
| `energy` | total energy of that term as computed (capped fragment, or capped pair) | Hartree |
| **`delta_energy`** | **for `level` 2 this is exactly `E_IJ − E_I − E_J`: the pair interaction energy** | Hartree |
| `distance` | **minimum interatomic distance between the two fragments' atoms**; `0` on monomer rows (a placeholder, not a measurement); caps are not included | **Ångström** |
| `scf` | `yes`/`no` — whether that term's SCF converged. Check it. | — |
| `homo`, `lumo`, `gap_ev` | frontier orbitals of that term | Hartree, Hartree, eV |
| `charge`, `mult` | of that term | — |

On a **monomer** row `delta_energy` is just a copy of `energy`, not a correction. Ignore it.

`frag_index` is not sorted by level — the dimers came first in every run here.

### Finding the ligand's row

The ligand is fragment number `len(fragments)` if you listed it last. Its pair rows are every
row with `level == 2` whose `m1`/`m2` includes that number.

```python
import csv, re
LIG = 4                                   # 1-based index of the ligand fragment
H = 627.5094740631                        # Hartree -> kcal/mol
for r in csv.DictReader(open("output_complex_fragments.csv")):
    m = sorted(int(r[k]) for k in r if re.fullmatch(r"m\d+", k or "") and r[k] != "0")
    if len(m) == 2 and LIG in m:
        other = [x for x in m if x != LIG][0]
        print(f"residue {other:3d}  d={float(r['distance']):6.2f} A  "
              f"{float(r['delta_energy'])*H:9.4f} kcal/mol  scf={r['scf']}")
```

Beware: naive `startswith("m")` also catches the `mult` column. Match `m\d+`.

### A worked ligand row

10 glycine residues + one water, HF/6-31G, MBE(2), 66 terms, 163.6 s on one thread:

| residue | distance / Å | ΔE / Hartree | ΔE / kcal mol⁻¹ |
|---:|---:|---:|---:|
| 1 | 12.349 | −0.000004050 | −0.0025 |
| 2 | 9.323 | −0.000026744 | −0.0168 |
| 3 | 5.355 | +0.000308454 | +0.1936 |
| 4 | 3.032 | +0.000888723 | +0.5577 |
| **5** | **1.790** | **−0.010207339** | **−6.4052** |
| 6 | 5.200 | +0.000165015 | +0.1035 |
| 7 | 7.558 | +0.000141243 | +0.0886 |
| 8 | 11.605 | −0.000074485 | −0.0467 |
| 9 | 14.750 | +0.000031421 | +0.0197 |
| 10 | 18.677 | −0.000004867 | −0.0031 |
| **sum** | | −0.008782629 | **−5.5112** |

Residue 5 carries the hydrogen bond; the rest decay with distance as they should. Negative is
attractive.

---

## 3. What is trustworthy, with numbers

### Each individual ligand pair energy is exact

For each ligand–residue pair I rebuilt the capped residue by hand and ran three ordinary
supermolecular calculations (`E(residue+ligand) − E(residue) − E(ligand)`) in the same basis:

| pair | distance / Å | MBE `delta_energy` / Ha | supermolecular / Ha | difference / Ha |
|---|---:|---:|---:|---:|
| res1–ligand | 4.929 | 0.0003328674 | 0.0003328674 | 1.1e-13 |
| res2–ligand | 1.790 | −0.0074571485 | −0.0074571485 | −2.8e-13 |
| res3–ligand | 3.576 | 0.0011896577 | 0.0011896577 | 6.8e-13 |

This is not an approximation that happens to be good — it is an identity, and the residual is
SCF convergence. **A row of this table is the supermolecular interaction energy of that capped
residue with the ligand, in vacuum, with no counterpoise correction.** If that is the quantity
you want, it is exact.

Internal consistency also holds: the total energy the run reports equals the sum of the whole
`delta_energy` column to 1.1e-13 Hartree.

### The pair energies do NOT sum to the binding energy

This is the important caveat and it is large.

| | kcal mol⁻¹ |
|---|---:|
| sum of the three ligand pair energies | **−3.724** |
| true supermolecular ligand binding energy | **−6.054** |
| **error** | **+2.330 (38 % of the binding energy, underbinding)** |

Adding the three-body terms (MBE level 3) recovers almost all of it — sum becomes −6.117,
error −0.063 kcal/mol — and level 4, which equals the fragment count, is exact.

**Most of that deficit is the capping, not ordinary many-body physics.** The control: the same
water bound to the middle of three *separate* glycine molecules, no covalent cuts anywhere,
same basis, same level:

| system | 2-body ligand sum | true | error |
|---|---:|---:|---:|
| 3 residues, backbone cut and capped | −3.724 | −6.054 | **+2.330** |
| 3 whole glycines, nothing cut | −5.047 | −4.737 | **−0.310** |

Capping inflates the two-body deficit by roughly **7.5×**. The missing energy is concentrated
in one three-body term, `(res2, res3, ligand)` at −2.698 kcal/mol — residues 2 and 3 are
covalently joined to each other and residue 2 is the one the ligand binds. The cap sitting
where residue 3's nitrogen belongs distorts exactly the carbonyl the ligand is hydrogen-bonded
to, and the three-body term is what repairs it. In a real protein every residue has two such
caps, so this applies everywhere, not just at the chain ends.

### Cap placement moves the numbers

`keywords.fragmentation.cap_scale` places the cap hydrogen at `R_H = R_kept + s·(R_gone − R_kept)`.
**The default is `s = 1.0`, which puts the cap hydrogen exactly on the position of the heavy
atom it replaces** — a ~1.37 Å N–H. That is deliberate in the code, not a bug, but it is not a
chemical geometry. Setting `s = 0.71` (roughly a real X–H length) moves the ligand pair
energies by:

| pair | s = 1.00 | s = 0.71 | shift |
|---|---:|---:|---:|
| res1–ligand | +0.2089 | +0.2907 | +0.082 |
| res2–ligand | −4.6794 | −5.0453 | **−0.366** |
| res3–ligand | +0.7465 | +0.6975 | −0.049 |

kcal/mol. About 8 % of the hydrogen bond. It is a real sensitivity, it is not negligible at
chemical accuracy, and there is no "right" answer the code will pick for you — so if you quote
pair energies, quote the `cap_scale` you used.

### So what can you defend?

- **Ranking residues by contribution, and the shape of the row.** Solid. The near residue
  dominates by two orders of magnitude and the tail decays monotonically.
- **An individual pair energy, quoted as "the vacuum interaction energy of this capped residue
  with the ligand."** Exact, to SCF convergence.
- **The absolute binding energy from summing pairs.** Do not. It was 38 % short here.
- **Chemical accuracy (1 kcal/mol) on any single pair adjacent to a cut.** Not established;
  cap placement alone moves it 0.37 kcal/mol.

Two further caveats that no keyword fixes: there is **no counterpoise correction** on these
numbers (`keywords.fragmentation.counterpoise: "vmfc"` exists but corrects the total, not the
pair rows you are reading), and every pair is computed **in vacuum** — no polarization from
the rest of the protein, because the embedded methods that would supply it refuse this
partition (§5).

---

## 4. Cost, and whether a protein is feasible

HF/6-31G, one thread, on a shared workstation:

| system | fragments | terms | wall |
|---|---:|---:|---:|
| 3 residues + ligand, MBE(2) | 4 | 10 | 16.2 s |
| 3 residues + ligand, MBE(3) | 4 | 14 | 81.6 s |
| 10 residues + ligand, MBE(2) | 11 | 66 | 163.6 s |
| 10 residues + ligand, MBE(2), `cutoffs.dimer = 5.0` | 11 | 30 | 67.3 s |
| 3 residues + ligand, one supermolecular SCF | — | 1 | 53.4 s |

**The scaling problem is that MBE(2) computes all N(N−1)/2 pairs, and you only want the N that
contain the ligand.** For 100 residues that is 5050 terms where 100 would do — a 50× waste.
There is a code path that accepts an explicit term list, but it is not reachable from a JSON
deck (C API and the geometry optimizer only), so today the only lever is distance screening:

```json
"fragmentation": {
  "method": "mbe", "level": 2,
  "cutoff_method": "distance",
  "cutoffs": { "dimer": 5.0 }
}
```

On the 10-residue case that cut 66 terms to 30 and the wall time from 163.6 s to 67.3 s, and
changed the ligand sum by 0.34 kcal/mol. It screens on the same minimum-interatomic-distance
the `distance` column reports, so it prunes distant residue–residue pairs and distant
ligand pairs alike — it cannot be told "keep everything touching the ligand". Choose the
cutoff against the `distance` column of a small run.

MPI works and the terms are independent, so ranks help roughly linearly (but see the
`cap_scale` bug in §6). A few hundred residues at MBE(2) in a small basis is plausible on a
cluster; it is the N² term count, not any single SCF, that will decide.

---

## 5. The other fragmentation methods, on a peptide

| `fragmentation.method` | on a backbone-cut partition |
|---|---|
| `mbe` | **works** — hydrogen-caps unconditionally. This recipe. |
| `gmbe` | runs, but **writes no `_fragments.csv` at all** and has no `delta_energy`; its JSON carries PIE terms with inclusion–exclusion coefficients, which are not pair interaction energies. Not usable for this. |
| `fmo` | **refuses**, clearly: `"fmo: the partition cuts a covalent molecule -- atoms 3 and 9 are covalently connected but were put in fragments 1 and 2 ... fragment on whole molecules"` |
| `ee-mbe` | **refuses**, same message. |
| `fmo` + `bond_breaking: "afo"` + `embedding: "none"` | **fails on a peptide bond**: `"fmo: the model system for the bond between atoms 3 and 9 could not be built"`. The frozen-orbital route works on propane (there is a shipped validation deck) but did not get through an amide C–N here. |
| `efmo` | **has no covalent-cut check at all, despite the documentation saying it refuses one.** See below. |

Accepted spellings of `method`, case-insensitive, `_` treated as `-`: `mbe`, `gmbe`, `fmo`,
`ee-mbe` / `eembe`, `efmo`.

### EFMO is the dangerous one

`mqc_docs/source/efmo.rst` states that EFMO refuses a partition that cuts a covalent bond.
**No such check exists in the code.** What happens instead:

- If a fragment ends with an odd electron count, it dies inside the fragment SCF with
  `"EFMO: RHF needs an even electron count; this system has an odd one and wants an
  unrestricted method"` — a message that says "system" when it means one fragment.
- If every fragment happens to come out even, **it runs.** I built exactly that case
  (three fragments, four severed peptide bonds, all electron counts even). EFMO ran for
  176 s, printed its full pair table and its success banner, wrote
  `output_*.json`, and **exited 0 with a total energy of `NaN`** and no error or warning
  anywhere in the log.

EFMO's per-pair table, for the water clusters it is actually valid on, is printed to **stdout
only** (not the JSON, not a CSV) with the header
`pair      R_IJ  class      contribution / Hartree`. The four named terms — Coulomb,
dispersion, exchange repulsion, charge transfer — are summed into that single `contribution`
column for far pairs and only appear separately as **system-wide totals** in the summary block
and in the JSON's `efmo` object. There is no per-pair four-term breakdown to read. `R_IJ` is
unitless (a van-der-Waals-contact-scaled separation), not an Ångström distance.

---

## 6. Traps on current main

1. **`embedding` is not an override.** Only `"none"` does anything. `"ptc"`, `"exact"` and
   anything else pass schema validation and are **silently ignored**; the embedding scheme is
   chosen by `method` (`fmo` = exact density, `ee-mbe` = point charges). A deck naming an
   embedding gets whichever one its method implies.
2. **`distance_metric` is dead.** It is parsed, allowed by the schema, and covered by a reader
   unit test — and then never read by anything. Every shipped deck says `"min"`, which happens
   to be the hard-coded behaviour. `"centroid"` would be accepted and ignored.
3. **`cap_scale` is not broadcast over MPI** (already flagged by a `TODO` in the source). On a
   multi-rank run, worker ranks cap at the default 1.0 whatever the deck asked for, so the
   fragments are inconsistent across ranks. **Either leave `cap_scale` at its default or run
   on one rank.**
4. **A refusal exits 0.** `fmo` and `ee-mbe` refusing the partition, and EFMO returning `NaN`,
   all returned exit status 0. Do not detect failure by exit code — grep the log, and check
   the `scf` column of the CSV.
5. **`level` 1 with EFMO** prints garbage QM rows (`0-0 0.000 ...`) in its pair table while the
   energy stays correct.
6. The CSV's `homo`/`lumo`/`gap_ev` columns have a known latent bug (a presence mask that is
   overwritten); harmless on this path today, but script against `delta_energy` and `distance`,
   not the orbital columns.

---

## 7. Checklist

- [ ] Ligand is its own fragment; note its 1-based index.
- [ ] Every atom in exactly one fragment.
- [ ] Every fragment-crossing bond in `connectivity`.
- [ ] `"fragmentation": {"method": "mbe", "level": 2}`.
- [ ] `cap_scale` left at default, or set *and* running on one rank.
- [ ] After the run: every `scf` column says `yes`.
- [ ] Read `delta_energy` (Hartree) from the `level == 2` rows containing the ligand.
- [ ] Quote them as vacuum, non-counterpoise-corrected pair interaction energies of *capped*
      residues — not as a decomposition of the binding energy.
