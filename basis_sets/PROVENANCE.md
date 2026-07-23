# Basis set provenance

| File | Source | Notes |
|---|---|---|
| `6-31G.txt`, `6-31Gs.txt` | pre-existing | GAMESS `$DATA` format |
| `sto-3g.gbs`, `cc-pvdz.gbs` | [Basis Set Exchange](https://www.basissetexchange.org) | Gaussian94 format, downloaded with `uncontract_spdf=1` |
| `def2-svp.gbs`, `def2-universal-jkfit.gbs` | NVIDIA cuEST sample data | Gaussian94 format, originally from the Basis Set Exchange |

## Adding a basis set

`find_basis_file` looks for `basis_sets/<normalized-name>.gbs` first and falls
back to `<normalized-name>.txt`, so either format works.

When downloading Gaussian94 (`.gbs`) files from the Basis Set Exchange, **ask
for `uncontract_spdf=1`**. Many basis sets (STO-3G, 6-31G, cc-pVDZ, ...) ship
combined `SP` shells in that format, and neither `mqc_gbs_reader` nor cuEST's
own helper handles them; the reader rejects them with a clear message rather
than guessing. The BSE web UI calls this option "Uncontract SPDF".

    curl "https://www.basissetexchange.org/api/basis/<name>/format/gaussian94/?elements=H,C,N,O&uncontract_spdf=1" \
        -o basis_sets/<name>.gbs

## Auxiliary bases are mandatory for the cuEST backend

cuEST exposes no conventional four-index ERI path, so J and K are always
density-fitted and every HF or DFT run needs a JKFIT auxiliary basis alongside
the orbital basis. `def2-universal-jkfit` is the default and is designed to
work with any orbital basis.

Basis set data is distributed under the terms of its respective source; see the
Basis Set Exchange for citation requirements.
