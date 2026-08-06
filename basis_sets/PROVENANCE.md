# Basis set provenance

Everything in this directory comes from one file: the
[Basis Set Exchange](https://www.basissetexchange.org) bundle
`basis_sets-json-0.12.tar.bz2`, downloaded from
[their downloads page](https://www.basissetexchange.org/downloads).

The bundle holds 922 basis sets in BSE's JSON format — 572 MB unpacked, far
too much to keep in git. It is tracked as the single compressed file, and
CMake extracts only the ones listed in `MQC_BASIS_SETS` at configure time.
The extracted `*.json` files are gitignored.

## Adding a basis set

Add its BSE name to `MQC_BASIS_SETS` in the top-level `CMakeLists.txt`, or
override the list at configure time:

```bash
cmake -B build -DMQC_BASIS_SETS="sto-3g;def2-svp;def2-universal-jkfit;def2-tzvpp"
```

BSE names are lowercase and spell `*` as `_st_`, so `6-31G**` is
`6-31g_st__st_`. Configure fails with the offending names listed if any is not
in the bundle.

Files keep their BSE name on disk — only the `.<revision>` suffix is stripped,
and where several revisions exist the highest wins. `normalize_basis_name`
produces that same spelling from whatever a user types, so there is no mapping
between the two to keep in step:

| input | resolves to |
|---|---|
| `def2-SVP` | `def2-svp.json` |
| `def2-SV(P)` | `def2-sv(p).json` |
| `6-31G**` | `6-31g_st__st_.json` |
| `cc-pVDZ` | `cc-pvdz.json` |

Parentheses are kept deliberately. def2-SV(P) and def2-SVP are different basis
sets — SV(P) leaves hydrogen unpolarized — and an earlier scheme that dropped
parentheses collapsed them onto one filename, along with their two RIFIT
auxiliaries.

`basis_sets/*.json` is derived entirely from `MQC_BASIS_SETS`: configure
deletes any JSON there that the list does not name, so dropping a name drops
the file.

## Where the code looks

`find_basis_file` searches, in order:

1. every directory in `$MQC_BASIS_PATH` (colon-separated)
2. `./basis_sets`, relative to the working directory
3. the `basis_sets/` of the source tree the binary was configured from

So a build tree works from anywhere with no setup, and `MQC_BASIS_PATH` points
at your own collection when you have one. A failed lookup names every path it
tried.

## Only JSON

BSE JSON is the sole supported format. The Gaussian94 (`.gbs`) and GAMESS
`$DATA` (`.txt`) readers were removed: JSON keeps a combined shell's
coefficient sets separate, so `SP` shells split cleanly with no
`uncontract_spdf=1` download flag to remember, it carries ECP data in the same
file for when that gets wired up, and one format means one parser to trust.

## Auxiliary bases are mandatory for the cuEST backend

cuEST exposes no conventional four-index ERI path, so J and K are always
density-fitted and every HF or DFT run needs a JKFIT auxiliary basis alongside
the orbital basis. `def2-universal-jkfit` is the default and is designed to
work with any orbital basis.

Basis set data is distributed under the terms of its respective source; see the
Basis Set Exchange for citation requirements.
