# Per-pair interaction energies for protein–ligand analysis in metalquicha

Measured on `origin/main` @ `31e4e0bd73`, build `cmake --preset default -DMQC_ENABLE_TBLITE=OFF`,
gfortran, `OMP_NUM_THREADS=1`. Every number below was produced by the runs described here.

**What works today:** plain MBE at level 2. You declare each residue as a fragment and the
ligand as one more, and you read the ligand's pair rows out of a CSV that is written by
default. Each ligand–residue number is *exactly* the supermolecular interaction energy of
that capped residue with the ligand — verified to 1e-13 Hartree.

**What does not work:** `fmo` and `ee-mbe` refuse a partition that cuts the backbone
unless it is detached with a frozen orbital, which is Part II. `efmo` refuses it.

### Since Part I was measured

Part I is a record of `31e4e0bd73`. These have changed on main since, and the sections
below say so where they are affected:

- a covalent partition that declares no `connectivity` is refused, rather than run uncapped;
- `efmo` refuses a covalent cut, rather than returning `NaN`;
- `cap_scale` reaches every rank;
- the CSV carries a `connected` column, and the pair table's reading is documented in
  `mqc_docs/source/json_output.rst` -- read that rather than §2's column list;
- EFMO's per-pair table is written to the JSON and a CSV sidecar with its four terms;
- `driver: "InteractionEnergy"` with `keywords.fragmentation.reference_fragment` computes
  only the terms holding the ligand (`mqc_docs/source/interaction_energy.rst`), which is the
  answer to §4's N² problem;
- `tools/sapt_scan` decomposes the ligand's interaction with each residue by SAPT.

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

*At `31e4e0bd73`* you got **no warning**: the audit that would catch it was gated on a
field the JSON path never filled. Main now runs it and refuses the partition. What used to
happen instead was luck:

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
*At `31e4e0bd73`* the only lever was distance screening. Main now has
`driver: "InteractionEnergy"`, which keeps only the terms holding one reference fragment and
is the right tool here; screening still works, and still prunes blindly:

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

MPI works and the terms are independent, so ranks help roughly linearly. A few hundred residues at MBE(2) in a small basis is plausible on a
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
| `efmo` | *at `31e4e0bd73`* had no covalent-cut check. **Main now refuses it.** |

Accepted spellings of `method`, case-insensitive, `_` treated as `-`: `mbe`, `gmbe`, `fmo`,
`ee-mbe` / `eembe`, `efmo`.

### EFMO was the dangerous one

Fixed on main; kept as the record of what an unguarded path did. `mqc_docs/source/efmo.rst` stated that EFMO refuses a partition that cuts a covalent bond.
**No such check exists in the code.** What happens instead:

- If a fragment ends with an odd electron count, it dies inside the fragment SCF with
  `"EFMO: RHF needs an even electron count; this system has an odd one and wants an
  unrestricted method"` — a message that says "system" when it means one fragment.
- If every fragment happens to come out even, **it runs.** I built exactly that case
  (three fragments, four severed peptide bonds, all electron counts even). EFMO ran for
  176 s, printed its full pair table and its success banner, wrote
  `output_*.json`, and **exited 0 with a total energy of `NaN`** and no error or warning
  anywhere in the log.

EFMO's per-pair table, for the water clusters it is actually valid on, was then printed to
**stdout only** (main now writes it to the JSON and a CSV sidecar, four terms per far pair) with the header
`pair      R_IJ  class      contribution / Hartree`. The four named terms — Coulomb,
dispersion, exchange repulsion, charge transfer — are summed into that single `contribution`
column for far pairs and only appear separately as **system-wide totals** in the summary block
and in the JSON's `efmo` object. There is no per-pair four-term breakdown to read. `R_IJ` is
unitless (a van-der-Waals-contact-scaled separation), not an Ångström distance.

---

## 6. Traps on current main

1. **`embedding` was not an override** at `31e4e0bd73` -- fixed by this branch, see §13.
2. **`distance_metric` is dead.** It is parsed, allowed by the schema, and covered by a reader
   unit test — and then never read by anything. Every shipped deck says `"min"`, which happens
   to be the hard-coded behaviour. `"centroid"` would be accepted and ignored.
3. **`cap_scale` was not broadcast over MPI** at `31e4e0bd73`. Fixed on main.
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
- [ ] The `cap_scale` you used, noted.
- [ ] After the run: every `scf` column says `yes`.
- [ ] Read `delta_energy` (Hartree) from the `level == 2` rows containing the ligand.
- [ ] Quote them as vacuum, non-counterpoise-corrected pair interaction energies of *capped*
      residues — not as a decomposition of the binding energy.

---

# Part II — FMO on a peptide, with frozen orbitals

Measured on `feat/fmo-peptide` (off `origin/feat/fmo-neutral-cut`, since merged), same
build and same single thread. Everything in Part I still holds; this part says what changed and what the
embedded methods give you that the plain expansion does not.

**What changed.** §5 said `fmo` + `bond_breaking: "afo"` "did not get through an amide C–N".
That was true and the diagnosis was wrong twice over. The refusal that fired was not about
the amide at all — it was `build_afo_model` reporting an **odd electron count in the model
system**, and it fired on the Cα–C bond just as readily. The model system is a sphere of
2.5 Å around the cut, closed off with one cap hydrogen per bond leaving it; on a peptide
that sphere reaches a neighbouring carbonyl **carbon** without reaching its **oxygen**, and
a cap hydrogen closes one electron pair where a C=O needs two. Every backbone model came
back a radical. A singly bonded heavy neighbour is now taken in whole instead of capped —
the same wholesale rule hydrogens already had, for the same continuity reason plus that one.

With that fixed the amide *is* refused, by the check that should have been reporting all
along:

```
fmo: 2 localized orbitals sit on the bond between atoms 3 and 9 (numbered from one), so it
is not a single bond. One frozen orbital stands in for one electron pair; cut at a single
bond. Atoms 3 and 9 are a peptide bond -- an amide C(=O)-N, conjugated with the carbonyl
and not a plain single bond. The FMO convention for a protein is to leave it whole and cut
the C-alpha--C(=O) bond one place along the backbone instead, which here is the bond
between atoms 2 and 3.
```

**Error messages number atoms from one. Decks number them from zero.** "Atoms 2 and 3"
above is the deck's atoms 1 and 2.

The default refusal — `bond_breaking` left at `"none"`, which is where a chemist starts —
used to name atoms 1 and 9, two nitrogens six bonds apart that are not a bond at all, and
close with "fragment on whole molecules", which is not advice you can take on a protein. It
now names the severed bond itself and both ways out:

```
fmo: the partition cuts a covalent molecule -- atoms 3 and 9 (numbered from one) are
covalently connected but were put in fragments 1 and 2. Nothing represents a cut bond here
unless it is asked to, so either fragment on whole molecules or set
keywords.fragmentation.bond_breaking to 'afo', which detaches the bond with a frozen
orbital. Atoms 3 and 9 are a peptide bond -- ...
```

---

## 8. The deck

Cut the **Cα–C(=O)** bond, not the peptide bond. A fragment is then one residue's carbonyl
together with the *next* residue's amine and Cα, and the peptide bond stays whole inside it.

```json
{
  "schema": { "name": "peptide-fmo", "version": "1.0" },
  "molecules": [{
    "xyz": "gly3_water_pair.xyz",
    "molecular_charge": 0,
    "molecular_multiplicity": 1,

    "fragments": [
      [0, 1, 4, 5, 6, 7],                            // N-terminal NH2-CH2
      [2, 3, 8, 9, 12, 13, 14],                      // C(=O) + NH + CH2
      [10, 11, 15, 16, 17, 18, 19, 20, 21, 22, 23],  // C(=O) + NH + CH2-COOH
      [24, 25, 26]                                   // THE LIGAND
    ],
    "fragment_charges":        [0, 0, 0, 0],
    "fragment_multiplicities": [1, 1, 1, 1]
  }],

  "model": { "method": "hf", "basis": "6-31g" },

  "keywords": {
    "scf": { "maxiter": 200, "tolerance": 1e-10 },
    "fragmentation": {
      "method": "fmo",
      "level": 2,
      "embedding": "ptc",
      "bond_breaking": "afo"
    }
  },

  "driver": "Energy",
  "system": { "logger": { "level": "verbose" } }
}
```

The geometry is `validation/inputs/sample_inputs/gly3_water_pair.xyz` — the shipped glycine
tripeptide with one water hydrogen-bonded to the residue-2 carbonyl at 1.893 Å. It is not
`gly3_water.xyz`, which main uses for the InteractionEnergy example with the water in a
different place; the numbers below are for this one.

### Four things that are different from Part I

1. **No `connectivity`.** Cut bonds are perceived from the geometry on this path. Declaring
   them is harmless and does nothing.
2. **`"embedding": "ptc"` is required, and now does something.** It used to pass validation
   and be ignored — the field followed `method` alone. It is now read straight through, and
   an unrecognised spelling is refused instead of dropped. `bond_breaking: "afo"` is refused
   with `"exact"`: a frozen orbital and an exact density both describe the bond region, and
   only a per-atom field can have the detached atom's share taken back out of it. Leaving
   the field off entirely (`"none"`) is worse than useless here — each side of a cut then
   carries about ±1 elementary charge.
3. **`"level": "verbose"`.** There is no `_fragments.csv` on this path. The pair energies
   are printed.
4. **The partition must be the Cα–C one.** An amide cut is refused, by name, with the bond
   to cut instead.

## 9. Reading the pair energies

*Superseded on `feat/fmo-pair-export`:* the pairs are now an info-level table sorted by
strength, with distances and connected pairs set apart, and an `fmo.pairs` array in the
output JSON -- see `mqc_docs/source/fmo.rst`, "Pair interaction energies". What this branch
printed at verbose, after the n-mer loop, was:

```
  fmo: n-mer interaction energies, Hartree
  fmo:   fragments                     dE
  fmo:                  1-2        -17.618384926866
  fmo:                  1-3         -0.012532415302
  fmo:                  1-4          0.000233279954
  fmo:                  2-3        -17.045007186136
  fmo:                  2-4          0.006221628160
  fmo:                  3-4         -0.005498823292
```

`1-4`, `2-4`, `3-4` are the ligand rows: fragment 4 is the water. For the FMO expansion each
of these is `E'_IJ − E'_I − E'_J` plus that pair's response to the field — the pair
interaction energy FMO calls an **IFIE**.

Two rows to ignore, and the reasons are different:

- **`1-2` and `2-3` are not interaction energies.** Those two pairs share a detached bond.
  Fragment 1's monomer holds a nucleus at `Z−1` and an empty frozen hybrid; the dimer holds
  the whole carbon and the whole bond pair. The difference carries the bond itself, which is
  why it is −17 Hartree and not −17 kcal/mol. Any pair whose two fragments are covalently
  joined reads like this.
- **On `method: "ee-mbe"` *every* row is a correction rather than an interaction energy**,
  including the ligand's. EE-MBE sums total embedded energies, and a monomer's already
  carries its electrostatics with every other fragment, so the pair term takes that back out
  with the opposite sign. Read as interaction energies they are nonsense: the same ligand
  rows come out `+0.140`, `+0.433`, `+0.374` Hartree — 88, 272 and 235 kcal/mol from one
  water. The header line says which of the two you are looking at. **Use `method: "fmo"` for
  pair analysis.** The two give the same *total* (they agree to 3e-8 Hartree here); they
  disagree about how it is split up, and only one of the splits is an interaction energy.

## 10. The numbers, against the plain expansion

Same geometry, same partition, same basis, one thread. HF/6-31G, `cap_scale` default.

| ligand pair | min. dist / Å | MBE(2) vacuum / kcal mol⁻¹ | FMO(2) IFIE / kcal mol⁻¹ | FMO − MBE |
|---|---:|---:|---:|---:|
| residue 1 – ligand | 4.728 | −0.0135 | +0.1464 | +0.160 |
| residue 2 – ligand | 3.159 | +1.4006 | +3.9041 | +2.504 |
| **residue 3 – ligand** | **1.893** | **−7.1906** | **−3.4506** | **+3.740** |
| sum | | −5.8035 | +0.5999 | +6.404 |

**They differ, and the sign is the same everywhere: every FMO IFIE is less attractive than
the vacuum pair energy, by more the closer the contact.** That is what the two quantities
are. In the plain expansion each pair is computed from two unpolarized fragments, so the
pair energy contains the induction — each fragment's density relaxing in the other's field.
In FMO that relaxation has already happened in the monomer self-consistency, against the
whole rest of the system, so the pair term is the interaction of two *already polarized*
fragments and the induction is not in it. The gap is therefore an induction energy: ~3.7
kcal/mol on the hydrogen bond, ~0.16 kcal/mol on a residue 4.7 Å away. Both the magnitude
and the distance dependence are what an induction term looks like, and the ordering of the
residues is unchanged.

**The IFIE sum is not a binding energy and is not meant to be.** The supermolecular ligand
binding energy on this geometry is **−6.0174 kcal/mol** (HF/6-31G, no counterpoise). The
MBE(2) vacuum pair sum lands at −5.8035, 0.21 short — better than Part I's 38 % because
the ligand is nowhere near a cut here. The FMO IFIE sum is +0.5999, 6.6 kcal/mol out, and
that is not an error in it: the induction it left in the monomer terms is most of the
binding. A binding energy out of FMO means running the peptide without the ligand and
subtracting totals, not summing IFIEs.

### The total, and what the truncation costs

| run | total / Hartree | error vs supermolecular |
|---|---:|---:|
| supermolecular RHF | −772.085165653 | — |
| MBE(2), hydrogen caps | −772.084073679 | +0.00109 Ha, **+0.69 kcal/mol** |
| FMO(2), `ptc` + `afo` | −772.040115669 | +0.04505 Ha, **+28.3 kcal/mol** |
| EE-MBE(2), `ptc` + `afo` | −772.040115640 | +0.04505 Ha, +28.3 kcal/mol |
| FMO(2), `afo`, no field | −772.029258666 | +0.05591 Ha, +35.1 kcal/mol |

**On the total, the embedded methods lose badly to capped MBE at level 2, and that is
expected rather than a fault.** It is the same result the propane table in `fmo.rst`
reports: across a detached bond the three-body term is large, truncating at pairs is not
advisable, and a point-charge field is at its worst at bonding contact — which is exactly
where `bond_breaking: "afo"` puts it, since an exact field is refused. Hydrogen caps happen
to be a good local model for a saturated backbone and the plain expansion benefits. The
embedding is worth a factor of twenty on a water cluster and costs you a factor of forty
here; neither figure travels.

### At full order it is exact, which is the point

STO-3G, same geometry, same partitions, level 4 = the fragment count:

| run | total / Hartree | error |
|---|---:|---:|
| supermolecular RHF | −762.311670855662 | — |
| MBE(4), hydrogen caps | −762.311670855659 | 2.5e-12 |
| FMO(4), `ptc` + `afo` | −762.311670855577 | 8.5e-11 |

Both telescope onto the supermolecule. That is the identity the whole boundary construction
is gated on — split nuclei, ghosts, frozen orbitals, group-by-group boundary sets — and it
holds on a protein backbone, not just on propane. The shipped regression for it is
`validation/inputs/cpu/mqc/fmo/afo3_gly3.json`, the tripeptide alone at level 3.

## 11. Cost

HF/6-31G, one thread, 27 atoms, 4 fragments:

| run | wall |
|---|---:|
| MBE(2), hydrogen caps | 19.6 s |
| FMO(2), `ptc` + `afo` | 622 s |
| EE-MBE(2), `ptc` + `afo` | 626 s |
| one supermolecular SCF | 53.3 s |

**FMO(2) costs about twelve supermolecular SCFs on a system this size**, nearly all of it in
the monomer self-consistency: every outer pass re-solves every fragment, and a fragment
carrying a ghost is solved in a basis wider than its own atoms. The N² pair count of Part I
§4 still applies on top. This is not a method to reach for on a four-fragment system; it
starts to pay where the monomer loop's cost is amortized over many more pairs than
fragments.

## 12. What to defend, on a peptide

- **Ranking residues by their IFIE.** Solid, and it agrees with the vacuum ranking.
- **An IFIE quoted as "the interaction of this residue with the ligand, both polarized by
  the rest of the peptide, excluding their mutual induction."** That is what it is. It is
  *not* the vacuum interaction energy and is not comparable to Part I's numbers.
- **The total at full expansion order.** Exact.
- **The total at level 2.** Do not. 28 kcal/mol out here, and worse than the plain
  expansion.
- **Summing IFIEs to a binding energy.** No. Off by the induction, which was 6.6 of 6.0
  kcal/mol here.

Unchanged from Part I: no counterpoise, and `cap_scale` does not apply on this path at all —
the frozen-orbital route has no caps outside its own model systems.

## 13. Traps, revised

§6's list still holds except for these:

1. **`embedding` is now an override**, and `"exact"`/`"ptc"`/`"none"` all do what they say.
   An unknown spelling is refused instead of ignored.
2. **A refusal still exits 0.** Grep the log.
3. **`efmo` refuses a covalent cut** on main. Do not point EFMO at a protein.
4. **`fmo` writes no fragment CSV**, so `system.fragment_breakdown` does nothing on this
   path. The pairs are in the output JSON under `fmo.pairs` (from `feat/fmo-pair-export`),
   and in an info-level table in the log.
