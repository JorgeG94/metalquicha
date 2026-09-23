# Per-residue interaction energies for a protein–ligand complex

**Read this first.** This is an experimental branch, `exp/protein-ligand`. It is not
`main`, it has not been reviewed, and it will not be merged in this form. It exists so
that one analysis — *which residues does my ligand interact with, and by how much* — can
be run today. Treat its numbers as you would an unpublished method: good enough to rank
residues and decide where to look next, and worth spot-checking against an ordinary
calculation before anything goes in a paper. §9 tells you which numbers are exact and
which are not, with measurements.

Everything you need is in this file. Beside it in the repository root is
`PAIR_ANALYSIS_RECIPE.md`, the measurement record this document is built from — read that
if you want the underlying tables rather than the conclusions.

---

## 1. What is on this branch

Development branches merged on top of `main` (`31e4e0bd73`). Three of them matter:

| branch | what it does for you |
|---|---|
| `fix/fragment-silent-failures` | Closes three ways of getting a wrong number with no warning. |
| `feat/fmo-peptide` | Makes the FMO method work on a protein backbone at all, and makes the `embedding` keyword do what it says. |
| `feat/efmo-pair-export` | Writes the EFMO method's pair table into the JSON output. Not used here; EFMO refuses a protein (§13). |

(`git log` shows five merges rather than three. Two of them,
`feat/fmo-afo-embedding` and `feat/fmo-neutral-cut`, are ancestors of `feat/fmo-peptide`
and were merged before it arrived; they are the frozen-orbital groundwork the peptide work
sits on, and merging the third on top of them brought nothing new. Nothing was dropped.)

The first two are both load-bearing for you, in different ways.

`fix/fragment-silent-failures` is the safety net. Before it, a deck that cut the peptide
backbone but forgot to declare the backbone bonds ran happily and printed a wrong number;
now it stops and tells you. (It caught a real mistake while this document was being
written — see trap 1.)

`feat/fmo-peptide` is what gives you a **second route**. Until it, the embedded FMO method
failed on every backbone cut with an unhelpful message. It now runs, which means you can
compare two different and complementary answers to the same question. §10 onward.

---

## 2. Building it

```bash
cmake --preset default -B build -DMQC_ENABLE_TBLITE=OFF
cmake --build build -j3
```

The first configure downloads dependencies, so it needs network. The build takes roughly
half an hour on a workstation; most of that is one integral library compiling several
thousand generated routines, and it is not a sign anything is wrong. The result is
`build/mqc`, a single executable that takes one argument: your input file.

`-DMQC_ENABLE_TBLITE=OFF` leaves out a semi-empirical engine this workflow does not use.

Run single-threaded:

```bash
export OMP_NUM_THREADS=1
```

Not because threading is broken, but because the pair energies are small differences of
large numbers, and a fixed thread count makes a run reproducible to the last digit. Raise
it if you want the speed — the energies move around the eleventh decimal, far below
anything that matters here.

---

## 3. What this computes, and what it cannot tell you

You declare each piece of the protein as a *fragment* and the ligand as one more. The code
runs an ordinary Hartree–Fock calculation on every fragment alone and on every **pair** of
fragments, and reports for each pair

```
ΔE(I,J) = E(I together with J) − E(I alone) − E(J alone)
```

For a pair made of **the ligand and one residue**, that is the interaction energy of that
residue with the ligand. Negative is attractive.

A residue is not a molecule — it is a piece cut out of the backbone — so something has to
stand in for the bond that was cut. That is the one real choice in this workflow, and it
is what separates the two routes:

| | **the plain route** (§4–§9) | **the FMO route** (§10–§12) |
|---|---|---|
| `method` | `mbe` | `fmo` |
| cut bond represented by | a hydrogen cap | a frozen orbital |
| each fragment sees | vacuum | a point-charge field of all the others |
| a pair energy is | the vacuum interaction of two *unpolarized* fragments | the interaction of two fragments *already polarized* by the whole system (an "IFIE") |
| output | a CSV, one row per term | a table printed to the log |
| cost, 3 residues + ligand | **16 s** | **~10 min** |

**Start with the plain route.** It is forty times cheaper, its numbers are exact in a
sense §9 makes precise, and it answers the ranking question. Reach for FMO when you want
the polarized picture, or want a second opinion on the ranking.

**What both give you:** a ranking. Which residues the ligand interacts with, in order,
with numbers attached. This is the solid part and it is what the method is for. The two
routes agree on the ranking (§11) while disagreeing on the values, which is a good sign.

**What neither can give you.** They will not tell you *why* a residue matters — whether a
given residue's −5 kcal/mol is electrostatics, or induction, or dispersion, or charge
transfer. That decomposition does exist in this code: it is the EFMO method's four-term
pair table, and one of the merged branches exports it. But it is only defined for a system
with **no covalent cuts** — a cluster of whole molecules — and a protein cut into residues
is precisely the case it refuses (§13). So for a peptide you get *how much*, not *of
what*. If you need the decomposition, you need a different partition, such as the whole
ligand against the whole protein as two intact molecules, which is a different calculation
from this one.

The one exception, and it is worth knowing: running **both** routes on the same partition
separates out a single component. The difference between a plain pair energy and the
corresponding FMO IFIE is that residue's **induction** — its polarization by the ligand and
vice versa — because FMO has already spent that in the monomer step and the plain route has
not. §11 shows it: 3.7 kcal/mol on the hydrogen bond, 0.16 kcal/mol on a residue 4.7 Å
away. That is one term out of four, obtained by subtraction rather than by a real
decomposition, and it is all you get.

---

## 4. The plain route: the deck

Input is a single JSON file. Here is a complete working one, for a glycine tripeptide with
one water molecule hydrogen-bonded to it as a ligand stand-in. The geometry is shipped:
`validation/inputs/sample_inputs/gly3_water.xyz`. Every number in this document was
produced from it.

```json
{
  "schema": { "name": "pair-analysis", "version": "1.0" },

  "molecules": [{
    "xyz": "complex.xyz",
    "molecular_charge": 0,
    "molecular_multiplicity": 1,

    "fragments": [
      [0,1,2,3,4,5,6,7],             // residue 1
      [8,9,10,11,12,13,14],          // residue 2
      [15,16,17,18,19,20,21,22,23],  // residue 3
      [24,25,26]                     // THE LIGAND, listed last
    ],
    "fragment_charges":        [0, 0, 0, 0],
    "fragment_multiplicities": [1, 1, 1, 1],

    "connectivity": [
      [2, 8, 1],
      [10, 15, 1]
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

JSON does not officially allow `//` comments — strip them before running. Field by field:

| field | what it is |
|---|---|
| `schema.name` | A label. It does **not** name the output file; the deck's filename does. |
| `molecules[0].xyz` | Path to an ordinary XYZ file, in Ångström, resolved relative to the deck. One file holding the whole complex, protein and ligand together. |
| `molecular_charge`, `molecular_multiplicity` | Of the whole complex. `0` and `1` for a neutral closed-shell system. |
| `fragments` | The partition. One list per fragment, of **0-based** atom indices into the XYZ file — the first atom is `0`. Every atom must appear in exactly one list. List the ligand last so its number is easy to remember. |
| `fragment_charges` | One integer per fragment, same order. A charged residue goes here (`-1` for a deprotonated aspartate). They must add up to `molecular_charge`. |
| `fragment_multiplicities` | One per fragment. `1` unless a fragment is genuinely open-shell. |
| `connectivity` | **The bonds you are cutting.** One `[atom_i, atom_j, order]` triple per bond crossing a fragment boundary, 0-based. Here the two peptide C–N bonds: atom 2 (residue 1's carbonyl carbon) to atom 8 (residue 2's nitrogen), and atom 10 to atom 15. **This is what tells the code to cap.** You do not need to list bonds *inside* a fragment, and the ligand — bonded to nothing — contributes no entries. The third number is a bond order and is not used. |
| `model.method` | `"hf"`. `"mp2"` also works and is better for dispersion-dominated contacts, at more cost. |
| `model.basis` | `"6-31g"` is the small, fast choice used throughout. |
| `keywords.scf.tolerance` | Convergence of each fragment calculation. Keep it tight: you are subtracting numbers near −200 Hartree to get answers near −0.007, so loose convergence lands directly in your answer. |
| `keywords.fragmentation.method` | `"mbe"` — the plain many-body expansion. |
| `keywords.fragmentation.level` | `2` — monomers and pairs. This is what produces the pair table. |
| `driver` | `"Energy"`. |
| `system.logger.level` | `"Info"`. `"verbose"` also prints the table to screen. |

### The three things that must be right

1. **Every atom in exactly one fragment.** An atom left out is refused.
2. **Every bond crossing a fragment boundary listed in `connectivity`.** Forget it and the
   run now stops (trap 1) — but a *wrong* bond is still not caught, so get it right.
3. **The ligand is its own fragment and bonded to nothing.** That is exactly why its rows
   are the trustworthy ones: no caps anywhere near them.

### Where to cut the backbone — this is worth 2 kcal/mol

The deck above puts **one residue in each fragment**, cutting the peptide C–N bond. That
is the obvious partition and it works. It is also not the best one, and the difference was
measured on this system:

| partition | cut bond | ligand pair sum | short of the true binding energy by | total energy error |
|---|---|---:|---:|---:|
| one residue per fragment | amide C–N | −3.540 kcal/mol | **2.48 kcal/mol (41 %)** | 3.94 kcal/mol |
| **Cα–C(=O)** | Cα–C | −5.804 kcal/mol | **0.21 kcal/mol (3.6 %)** | 0.69 kcal/mol |

The second partition groups one residue's carbonyl with the *next* residue's amine and Cα,
leaving each peptide bond whole inside a fragment:

```json
"fragments": [
  [0, 1, 4, 5, 6, 7],                            // N-terminal NH2-CH2
  [2, 3, 8, 9, 12, 13, 14],                      // C(=O) + NH + CH2
  [10, 11, 15, 16, 17, 18, 19, 20, 21, 22, 23],  // C(=O) + NH + CH2-COOH
  [24, 25, 26]                                   // the ligand
],
"connectivity": [ [1, 2, 1], [9, 10, 1] ]
```

It is better here for a concrete reason: this ligand hydrogen-bonds to a **carbonyl
oxygen**, and the one-residue partition puts a cap hydrogen on the nitrogen right next to
that carbonyl, distorting the very group the ligand is binding. The Cα–C cut puts the caps
a bond further away. It is also the convention the protein FMO literature uses, and — not
a coincidence — **the only partition the FMO route will accept** (§10).

So: if your ligand binds a backbone carbonyl or amide, prefer the Cα–C partition. If it
binds a side chain, either works and the one-residue partition is easier to read. Note
that "residue 3" then means the third *fragment*, not the third residue; keep a note of
which is which.

---

## 5. Running it, and where the output goes

```bash
cd /wherever/your/deck/is
export OMP_NUM_THREADS=1
/path/to/build/mqc complex.json
```

Two files appear **in the directory you ran from**, named after the deck's filename:

| file | what is in it |
|---|---|
| `output_complex.json` | Total energy and summary information. |
| `output_complex_fragments.csv` | **One row per term. This is your pair table.** |

The CSV is written by default. `"system": {"fragment_breakdown": "csv"}` controls it;
`"json"` folds the same data into the JSON instead, `"none"` suppresses it.

The tripeptide-plus-water example took **16 s** on one thread on a shared workstation. A
ten-residue peptide plus the same ligand — 11 fragments, 66 terms — took **144 s**.

**Check the exit status and read the log.** A genuine input error exits 1. A refusal by one
of the *other* fragmentation methods exits 0 anyway (trap 4).

---

## 6. Reading the CSV

```
frag_index,level,m1,m2,energy,delta_energy,distance,scf,homo,lumo,gap_ev,charge,mult
```

| column | meaning | units |
|---|---|---|
| `frag_index` | Row counter. Rows are not sorted by level; pairs came first in every run here. | — |
| `level` | `1` = one fragment alone, `2` = a pair. | — |
| `m1`, `m2` | Which fragments, **1-based**, in the order of your `fragments` list. The ligand, listed fourth, is `4`. On a `level` 1 row `m2` is `0` — filler, not a fragment. | — |
| `energy` | Total energy of that term. | Hartree |
| **`delta_energy`** | On a `level` 2 row, `E_IJ − E_I − E_J`: **the pair interaction energy. This is the column you want.** | Hartree |
| `distance` | Closest approach between the two fragments' atoms. `0` on `level` 1 rows — a placeholder, not a measurement. Cap hydrogens excluded. | **Ångström** |
| `scf` | `yes`/`no`, whether that term converged. **Check every row says `yes`.** | — |
| `homo`, `lumo`, `gap_ev` | Frontier orbitals of that term. Not needed here, and these three have a known latent bug — script against `delta_energy` and `distance` only. | Hartree, Hartree, eV |
| `charge`, `mult` | Of that term. | — |

On a `level` 1 row `delta_energy` is a copy of `energy`. Ignore it there.

**1 Hartree = 627.5095 kcal/mol.**

### Pulling out the ligand's rows

```python
import csv, re

CSV = "output_complex_fragments.csv"
LIG = 4                        # 1-based index of the ligand fragment
H   = 627.5094740631           # Hartree -> kcal/mol

for r in csv.DictReader(open(CSV)):
    members = sorted(int(r[k]) for k in r
                     if k and re.fullmatch(r"m\d+", k) and r[k].strip() != "0")
    if len(members) != 2 or LIG not in members:
        continue
    other = [x for x in members if x != LIG][0]
    d, de = float(r["distance"]), float(r["delta_energy"])
    print(f"residue {other:3d}   d = {d:6.3f} A   "
          f"dE = {de: .9f} Ha = {de*H: 9.4f} kcal/mol   scf={r['scf'].strip()}")
```

One gotcha there is deliberate: match the member columns with `m\d+`, not with
`startswith("m")`, which would also pick up the `mult` column.

---

## 7. A worked example, with the answer

Glycine tripeptide, one water hydrogen-bonded to residue 2's carbonyl oxygen at 1.893 Å.
HF/6-31G, `method: "mbe"`, `level: 2`, one residue per fragment, ligand as fragment 4. The
script above prints:

```
residue   1   d =  4.728 A   dE =  0.001111190 Ha =    0.6973 kcal/mol   scf=yes
residue   2   d =  1.893 A   dE = -0.007282852 Ha =   -4.5701 kcal/mol   scf=yes
residue   3   d =  3.645 A   dE =  0.000530023 Ha =    0.3326 kcal/mol   scf=yes
```

**Residue 2 is the answer.** It is the only attractive term, it is six times larger in
magnitude than anything else, and it is the residue at contact distance. That is the
result, and it is robust.

### Read the tail with care

Notice that residue 1, at 4.73 Å, has a *larger* number than residue 3 at 3.65 Å. That is
not a bug and it is worth understanding before you over-read a table like this. Terms this
small are electrostatic interactions with oriented backbone dipoles, and a dipole's field
depends on which way it points as much as on how far away it is. **Below about 1 kcal/mol
the ordering is not physically meaningful.** Rank the residues that matter; do not rank the
ones that do not.

The envelope does behave, once there are enough residues to see it. The same water against
a ten-residue glycine peptide, bound at residue 5 (11 fragments, 66 terms, 144 s):

```
residue   1   d = 14.661 A   dE =    0.0208 kcal/mol
residue   2   d = 11.461 A   dE =   -0.0338 kcal/mol
residue   3   d =  7.471 A   dE =    0.1125 kcal/mol
residue   4   d =  5.093 A   dE =    0.1702 kcal/mol
residue   5   d =  1.790 A   dE =   -4.3545 kcal/mol      <-- the hydrogen bond
residue   6   d =  3.372 A   dE =    0.6326 kcal/mol
residue   7   d =  5.469 A   dE =    0.1926 kcal/mol
residue   8   d =  9.306 A   dE =   -0.0058 kcal/mol
residue   9   d = 12.536 A   dE =    0.0284 kcal/mol
residue  10   d = 16.285 A   dE =   -0.0118 kcal/mol
```

One residue carries the interaction, its two neighbours carry tenths of a kcal/mol, and
past about 7 Å everything is hundredths with no meaningful order. That is the shape to
expect, and a run that does not look like this is a run to be suspicious of.

---

## 8. Which rows mean something, and which are artefacts

The example CSV has ten rows. Only some are interpretable.

**Meaningful — the three ligand pairs.** Every `level` 2 row containing the ligand's
fragment number. These are the answer.

**Artefacts — pairs of fragments that are covalently joined to each other.** In the example
these are rows `1,2` and `2,3`:

```
1,2,1,2, ... ,  1.0109126178E+00 ,  1.3669269846E+00 , yes, ...
4,2,2,3, ... ,  1.0135188652E+00 ,  1.3671965991E+00 , yes, ...
```

**+1.011 and +1.014 Hartree** — about 635 kcal/mol, three orders of magnitude larger than
anything else in the table. That is not an interaction energy. When two fragments that
share a backbone bond are put together the bond is restored and the two cap hydrogens that
stood in for it are removed, so `delta_energy` is dominated by that chemical bookkeeping.
**Ignore these rows entirely.**

They are easy to spot by two signs together: their `distance` column reads a *bond* length
rather than a contact — about 1.37 Å for the amide partition, 1.50 Å for the Cα–C one —
and their energy is in whole Hartree while everything else is in thousandths.

> If you were told these rows read "around minus fourteen Hartree": they do not, on this
> branch or this system. They read **+1.01 Hartree**, positive, for the amide cut and
> **+1.02** for the Cα–C cut, and the recipe these notes come from records no −14 figure
> either. The instruction to ignore them is right; the value and sign quoted for them were
> not. Identify them by the bond-length `distance` and the three-orders-of-magnitude jump,
> never by matching a specific number. (On the FMO route the same rows really are large and
> negative, around −17 Hartree — see §11. That is a different quantity on a different path,
> and may be where the figure came from.)

**Artefacts — pairs of *un*joined fragments.** Row `1,3` here. An honest residue–residue
interaction, but both partners carry caps, so it is the least trustworthy number in the
table. You do not need it for a ligand analysis.

**Not an interaction at all — the `level` 1 rows.** Fragment total energies. Ignore.

---

## 9. How accurate is it

### Each ligand pair energy, on its own: exact

Not "accurate" — *exact*, and this was checked rather than assumed. For each ligand pair
the capped residue was rebuilt by hand as its own XYZ file and three ordinary calculations
were run, `E(residue+ligand) − E(residue) − E(ligand)`, same basis, no fragmentation
anywhere:

| pair | `delta_energy` from the CSV | plain supermolecular | difference |
|---|---:|---:|---:|
| residue 1 – ligand | 0.001111190029 | 0.001111190029 | 0 |
| residue 2 – ligand | −0.007282851744 | −0.007282851745 | 1.1e-12 |
| residue 3 – ligand | 0.000530023141 | 0.000530023141 | 0 |

The residual is SCF convergence. So a row of your table **is** the vacuum interaction
energy of that capped residue with the ligand, with no counterpoise correction and no
polarization from the rest of the protein. If that is the quantity you want, you have it
to machine precision.

### Do not add the column up and call it a binding energy

This is the one thing to take away.

| | kcal/mol |
|---|---:|
| sum of the three ligand pair energies (one residue per fragment) | **−3.540** |
| the actual binding energy, from a plain calculation of the whole complex | **−6.017** |
| error | **+2.477 — 41 % short** |

The recipe this document is built on measured the same quantity independently on a
slightly differently oriented water and got 38 %. So: **roughly 40 % short, underbinding.**

Two things cause it and only one is ordinary physics. Two-body truncation is the ordinary
part; adding three-body terms (`"level": 3`) recovers nearly all of it, at about five times
the cost. The capping is the other part and it is the larger one. `PAIR_ANALYSIS_RECIPE.md`
§3 has the control: the same water against three *separate*, uncut glycine molecules is
only 0.3 kcal/mol out at the same level. Cutting the backbone inflates the two-body deficit
about sevenfold, and it does so right where it hurts — the cap that stands in for the
neighbouring residue sits on the group the ligand is bound to.

Moving the cut to Cα–C (§4) reduces the shortfall from 2.48 kcal/mol to **0.21**, which is
the practical lesson: the deficit is mostly about *where you cut relative to the ligand*,
not about the method.

So:

- **Rank residues by contribution, and read the shape of the decay** — defensible.
- **Quote one pair energy as "the vacuum interaction energy of this capped residue with the
  ligand"** — exact.
- **Add the column up and call it the binding energy** — no. It was 41 % short here. If you
  need a binding energy, compute it the ordinary way: the whole complex minus the protein
  minus the ligand, three calculations.
- **Claim 1 kcal/mol accuracy on a pair next to a cut** — not established; the `cap_scale`
  trap below moves the hydrogen bond by 0.37 kcal/mol on its own.

---

## 10. The FMO route: the deck

This is the route `feat/fmo-peptide` made possible, and it did not work before this branch.
It is worth running when you want each residue's interaction with the ligand *as it is in
the protein* — polarized by everything around it — rather than in vacuum.

```json
{
  "schema": { "name": "peptide-fmo", "version": "1.0" },
  "molecules": [{
    "xyz": "complex.xyz",
    "molecular_charge": 0,
    "molecular_multiplicity": 1,

    "fragments": [
      [0, 1, 4, 5, 6, 7],
      [2, 3, 8, 9, 12, 13, 14],
      [10, 11, 15, 16, 17, 18, 19, 20, 21, 22, 23],
      [24, 25, 26]
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

Four differences from the plain deck, all of them mandatory:

1. **The partition must be the Cα–C one.** An amide cut is refused (§13), because a frozen
   orbital stands in for one electron pair and a conjugated amide C–N is not a plain single
   bond. The refusal names the bond to cut instead.
2. **No `connectivity`.** Cut bonds are perceived from the geometry on this path. Declaring
   them is harmless and does nothing.
3. **`"embedding": "ptc"` is required.** This keyword used to be silently ignored; on this
   branch it is honoured, and an unrecognised spelling is refused rather than dropped.
   `"exact"` is refused together with `bond_breaking: "afo"` — a frozen orbital and an
   exact density both describe the bond region, and only a per-atom field can have the
   detached atom's share taken back out. `"none"` runs but is worse than useless here, as
   each side of a cut then carries about ±1 elementary charge.
4. **`"level": "verbose"`.** There is no CSV on this path; the pair energies are printed to
   the log and nowhere else.

`bond_breaking: "afo"` is the frozen-orbital machinery — "adjusted frozen orbital". It
replaces the hydrogen cap with an orbital frozen in the shape the real bond has, which is
why there is no `cap_scale` to worry about on this route.

**This route costs about forty times the plain one** on a system this size — about ten
minutes against sixteen seconds for three residues plus a ligand — nearly all of it in the
monomer self-consistency, which re-solves every fragment on every outer pass.

---

## 11. Reading the FMO pair table

At `"level": "verbose"` the run prints this after the n-mer loop:

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

**Read the header line every time.** `n-mer interaction energies` is what you want. If it
says something else, you are not looking at interaction energies — see §12.

`1-4`, `2-4` and `3-4` are the ligand rows; fragment 4 is the water. Each is an **IFIE** —
inter-fragment interaction energy — the standard FMO quantity.

Two rows to ignore. **`1-2` and `2-3` are not interaction energies**: those pairs share a
detached bond, and the difference between the pair and the two fragments carries the bond
itself. That is why they read **−17 Hartree** rather than −17 kcal/mol. Same rule as the
plain route (§8) with a different sign and a much larger magnitude, because a frozen
orbital and a hydrogen cap fail to cancel in different ways. Any pair whose two fragments
are covalently joined reads like this.

### The FMO numbers next to the plain ones

Same geometry, same partition, same basis, one thread:

| ligand pair | distance / Å | plain MBE(2) / kcal mol⁻¹ | FMO(2) IFIE / kcal mol⁻¹ | FMO − MBE |
|---|---:|---:|---:|---:|
| fragment 1 – ligand | 4.728 | −0.0135 | +0.1464 | +0.160 |
| fragment 2 – ligand | 3.159 | +1.4006 | +3.9041 | +2.504 |
| **fragment 3 – ligand** | **1.893** | **−7.1906** | **−3.4506** | **+3.740** |
| sum | | −5.8035 | +0.5999 | +6.403 |

**The two routes disagree, they disagree in one direction, and that is not an error in
either.** Every FMO value is less attractive than the vacuum pair energy, by more the
closer the contact. That is what the two quantities *are*. In the plain expansion a pair is
computed from two unpolarized fragments, so the pair energy contains the induction — each
fragment's density relaxing in the other's field. In FMO that relaxation has already
happened during the monomer self-consistency, against the whole rest of the system, so the
pair term is the interaction of two *already polarized* fragments and the induction is not
in it.

The gap between the columns is therefore an induction energy: ~3.7 kcal/mol on the hydrogen
bond, ~0.16 kcal/mol on a fragment 4.7 Å away. Both the size and the distance dependence
are what induction looks like.

**The ranking is identical.** Fragment 3 dominates on both, fragment 1 is negligible on
both. That is the result you use, and the two routes agreeing on it is the reassurance you
get from running both.

### Do not sum the IFIEs either — and for a different reason

| | kcal/mol |
|---|---:|
| plain MBE(2) ligand pair sum | −5.804 |
| FMO(2) IFIE sum | **+0.600** |
| true supermolecular binding energy | −6.017 |

The plain sum lands 0.21 short, which is good (§4 explains why this partition does so much
better than one-residue-per-fragment). The FMO sum is **+0.60 — the wrong sign entirely,
6.6 kcal/mol out.**

**That is not a bug and it is not an accuracy problem. It is a category error.** The
induction that FMO moved out of the pair terms and into the monomer terms *is* most of the
binding here, so the IFIEs are missing it by construction. Summing them is asking the wrong
question of the right numbers. A binding energy out of FMO means running the peptide with
and without the ligand and subtracting totals — not summing IFIEs. §9's warning applies to
the plain route; this one applies to FMO and is sharper.

### One more thing FMO is not good at here

| run | total energy / Hartree | error vs the supermolecule |
|---|---:|---:|
| supermolecular RHF | −772.085165653 | — |
| plain MBE(2), hydrogen caps, Cα–C cut | −772.084073679 | +0.69 kcal/mol |
| FMO(2), `ptc` + `afo` | −772.040116 | **+28.3 kcal/mol** |

At level 2 the embedded method loses badly to plain capped MBE on the *total* energy. That
is expected rather than a fault: across a detached bond the three-body term is large, and a
point-charge field is at its worst at bonding contact — which is exactly where
`bond_breaking: "afo"` puts it. At full expansion order both are exact to about 1e-11
(`PAIR_ANALYSIS_RECIPE.md` Part II §10). **Use FMO for the polarized pair picture, not for
a total energy at level 2.**

---

## 12. EE-MBE: a wrong answer that looks like a right one

There is a third method, `"method": "ee-mbe"`, which takes the same deck and the same
keywords and runs to completion. **Its pair rows are not interaction energies and reading
them as such will give you a badly wrong figure.**

EE-MBE sums total *embedded* energies. A monomer's energy already contains its
electrostatic interaction with every other fragment, so the pair term exists to take that
double-counting back out again, with the opposite sign. It is a correction inside a
bookkeeping scheme, not a physical interaction.

Here is what that looks like on the same deck, with `"method": "ee-mbe"` substituted:

```
  fmo: n-mer terms of the embedded expansion, Hartree -- corrections, not interaction energies: a monomer's own energy already holds its electrostatics with every other fragment
  fmo:   fragments                     dE
  fmo:                  1-2        -22.065510665003
  fmo:                  1-3          0.020915467313
  fmo:                  1-4          0.140503657230
  fmo:                  2-3        -20.723112967118
  fmo:                  2-4          0.432996260993
  fmo:                  3-4          0.373709084274
```

The header line says what they are, which is the check to make. The ligand rows are
`+0.1405`, `+0.4330` and `+0.3737` Hartree — which is **+88, +272 and +235 kcal/mol**. One
water molecule. Those are not interaction energies by any reading; they are five to seventy
times the entire binding energy, all repulsive, from a molecule that is in fact bound.

The trap is that nothing goes wrong. The run finishes, exits 0, prints a table in the same
shape as FMO's with the same column headings, and the total energy it produces is correct —
EE-MBE and FMO agree on the total here to 2.9e-8 Hartree. **The two methods differ only in
how they split that total up, and only one of the splits is an interaction energy.**

**Use `"method": "fmo"` for pair analysis.** If you see a table whose header does not say
`interaction energies`, or whose ligand rows are in the hundreds of kcal/mol, that is this.

---

## 13. Traps

**1. A deck that forgets `connectivity`.** On `main` this ran to completion and printed a
wrong number with no diagnostic. On this branch it stops:

```
invalid system: the geometry implies 2 bond(s) crossing monomer boundaries that were
never declared, starting with atoms 2 and 8; those fragments would have uncapped
valences. No connectivity was declared at all, so no bond is marked broken and no
fragment is capped; list the bonds in the molecule's connectivity, or set
system.unchecked_input if the partition is deliberate
```

and exits 1. Atom numbers in *this* message are 0-based, matching your deck. Closing this
hole is the main reason this branch exists — and it earned its place while this document
was being written, by catching a Cα–C deck here that had been copied from the FMO one and
had lost its `connectivity` in the process. Do **not** set `system.unchecked_input` to get
past it; that switch is for deliberately unusual partitions and turns the refusal back into
a warning you will not read. Note the exemption: `bond_breaking: "afo"` decks are not
audited, because they perceive their own cuts.

**2. `embedding` used to be ignored, and now is not.** On `main`,
`keywords.fragmentation.embedding` accepted `"ptc"` and `"exact"`, passed validation, and
was then **silently discarded** — the field was chosen by `method` alone, which made the
one configuration a detached bond can run in unreachable from a deck. On this branch it is
read, and an unknown spelling is refused instead of dropped. If you find older notes
saying "embedding does nothing unless it says none", they described the old behaviour.

**3. `cap_scale` is worth about 8 %, on the plain route.**
`keywords.fragmentation.cap_scale` sets where the cap hydrogen goes:
`R_H = R_kept + s·(R_gone − R_kept)`. **The default is `s = 1.0`, which puts the cap
hydrogen exactly on the position of the heavy atom it replaces** — a ~1.37 Å N–H, not a
chemical bond length. That is deliberate in the code, not a bug. Setting `s = 0.71`, about
a real X–H distance, moves the ligand pair energies by up to 0.37 kcal/mol, roughly 8 % of
the hydrogen bond. There is no right answer the code will pick for you, so **if you quote
pair energies, quote the `cap_scale` you used.** (On `main` this keyword was not sent to
other MPI ranks, so a parallel run silently ignored it and ranks disagreed about the same
fragment; that is fixed on this branch and `cap_scale` is now safe on more than one rank.)
It does not apply to the FMO route, which has no caps.

**4. A refusal can still exit 0.** The fragmentation methods refuse a partition and return
status 0 while doing it. Read the log; check the `scf` column.

**5. `distance_metric` does nothing.** Parsed, accepted, never read. The distance is always
closest approach.

**6. Cost grows as the square of the fragment count.** Level 2 computes *all* N(N−1)/2
pairs when you only want the N containing the ligand. Ten residues plus a ligand is 66
terms; a hundred residues would be 5050 where 100 would do. There is no way to ask for just
the ligand pairs from a deck today. The lever available is distance screening:

```json
"fragmentation": {
  "method": "mbe", "level": 2,
  "cutoff_method": "distance",
  "cutoffs": { "dimer": 5.0 }
}
```

On a ten-residue case that cut 66 terms to 30 and halved the wall time, moving the ligand
sum by 0.34 kcal/mol. It screens on the same closest-approach distance the CSV reports, so
it prunes distant residue–residue pairs and distant *ligand* pairs alike — it cannot be
told "keep everything touching the ligand". Pick the cutoff by looking at the `distance`
column of a small run first.

### What the other methods do on a peptide

All of these were run on the example system on this branch.

| `method` | what happens |
|---|---|
| `mbe` | **Works.** Caps the cut backbone. §4–§9. |
| `fmo` + `afo` + `ptc` | **Works**, on the Cα–C partition. §10–§11. |
| `gmbe` | Runs, but writes no `_fragments.csv` at all. No pair table. Not usable here. |
| `ee-mbe` | Runs, and its pair rows are not interaction energies. §12. |
| `fmo` with `bond_breaking` left at its default | Refuses, and tells you both ways out: `fmo: the partition cuts a covalent molecule -- atoms 3 and 9 (numbered from one) are covalently connected but were put in fragments 1 and 2. Nothing represents a cut bond here unless it is asked to, so either fragment on whole molecules or set keywords.fragmentation.bond_breaking to 'afo', which detaches the bond with a frozen orbital. Atoms 3 and 9 are a peptide bond -- an amide C(=O)-N, conjugated with the carbonyl and not a plain single bond. The FMO convention for a protein is to leave it whole and cut the C-alpha--C(=O) bond one place along the backbone instead, which here is the bond between atoms 2 and 3.` |
| `fmo` + `afo` on an **amide** partition | Refuses, naming the bond to cut instead: `fmo: 2 localized orbitals sit on the bond between atoms 3 and 9 (numbered from one), so it is not a single bond. One frozen orbital stands in for one electron pair; cut at a single bond. Atoms 3 and 9 are a peptide bond -- ... The FMO convention for a protein is to leave it whole and cut the C-alpha--C(=O) bond one place along the backbone instead, which here is the bond between atoms 2 and 3.` **This is the message that tells you to switch to the Cα–C partition.** |
| `efmo` | Refuses: `efmo: the partition cuts a covalent molecule -- atoms 1 and 9 are covalently connected but were put in fragments 1 and 2. A hydrogen cap's multipoles would act on its partner across the cut and no frozen-orbital route is wired in here, so this method cannot answer for that partition; fragment on whole molecules` — **this refusal is new on this branch.** On `main`, EFMO ran this system for three minutes, printed a success banner, and reported a total energy of `NaN`. If you find notes saying EFMO has no covalent-cut check, they describe `main`, not this branch. |

Atom numbers in the `fmo`/`ee-mbe`/`efmo` messages count from **1**, not from 0 like your
deck. The `efmo` message additionally names an arbitrary pair of atoms in the same molecule
rather than the bond itself. Only the connectivity-audit message in trap 1 uses your deck's
own numbering.

---

## 14. Checklist before you trust a run

- [ ] The ligand is its own fragment; you know its 1-based number.
- [ ] Every atom appears in exactly one fragment.
- [ ] Every bond crossing a fragment boundary is in `connectivity` (plain route), or you are
      on `bond_breaking: "afo"` and deliberately declaring none (FMO route).
- [ ] If the ligand binds a backbone carbonyl or amide, you are cutting at Cα–C, not at the
      peptide bond.
- [ ] The run exited 0 and the log has no refusal in it.
- [ ] Every `scf` column in the CSV says `yes`.
- [ ] You are reading `delta_energy` from `level == 2` rows containing the ligand.
- [ ] You are ignoring the covalently-joined pairs — `distance` at a bond length, energy in
      whole Hartree.
- [ ] You are not adding the column up and calling it a binding energy.
- [ ] You are quoting the `cap_scale` you used (plain route), or saying you used FMO with
      `afo`/`ptc` (FMO route), and never mixing numbers from the two.
