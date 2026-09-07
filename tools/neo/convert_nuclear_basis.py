#!/usr/bin/env python3
"""Turn PySCF-NEO's nuclear basis files into the BSE JSON this code reads.

Source: pyscf/neo/basis/*.dat in github.com/theorychemyang/pyscf, the
even-tempered proton sets of Yu, Pavosevic and Hammes-Schiffer, J. Chem.
Phys. 152, 244123 (2020). Each line pair is one uncontracted shell on H.

    python tools/neo/convert_nuclear_basis.py <pyscf-neo>/pyscf/neo/basis basis_sets/neo
"""
import json
import re
import sys
from pathlib import Path

L_OF = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4, "H": 5, "I": 6}
NAMES = {  # file stem -> the name a deck uses
    "pb4d": "pb4-d", "pb4f1": "pb4-f1", "pb4f2": "pb4-f2", "pb4f2a": "pb4-f2a",
    "pb5d": "pb5-d", "pb5f": "pb5-f", "pb5g": "pb5-g",
    "pb6d": "pb6-d", "pb6f": "pb6-f", "pb6g": "pb6-g", "pb6h": "pb6-h",
    "dzsnb": "dzsnb", "1s1p": "1s1p", "2s1p": "2s1p",
}


def convert(dat: Path, out_dir: Path) -> Path:
    name = NAMES.get(dat.stem, dat.stem)
    lines = [ln.strip() for ln in dat.read_text().splitlines()]
    description = next((ln.lstrip("#").strip() for ln in lines if ln.startswith("#")), "")
    shells = []
    i = 0
    body = [ln for ln in lines if ln and not ln.startswith("#")]
    while i < len(body):
        m = re.match(r"^H\s+([SPDFGHI])$", body[i])
        if not m:
            raise ValueError(f"{dat}: unexpected line {body[i]!r}")
        l = L_OF[m.group(1)]
        i += 1
        exps, coefs = [], []
        while i < len(body) and not re.match(r"^H\s+[SPDFGHI]$", body[i]):
            e, c = body[i].split()[:2]
            exps.append(f"{float(e):.10E}")
            coefs.append(f"{float(c):.10E}")
            i += 1
        shells.append({"function_type": "gto", "region": "",
                       "angular_momentum": [l], "exponents": exps,
                       "coefficients": [coefs]})
    doc = {
        "molssi_bse_schema": {"schema_type": "component", "schema_version": "0.1"},
        "name": name,
        "description": f"Nuclear (proton) basis {name}: {description}",
        "role": "orbital",
        "family": "neo",
        "function_types": ["gto"],
        "elements": {"1": {"electron_shells": shells,
                           "references": [{"reference_description":
                                           "J. Chem. Phys. 152, 244123 (2020); file from pyscf/neo",
                                           "reference_keys": []}]}},
    }
    out = out_dir / f"{name}.json"
    out.write_text(json.dumps(doc, indent=1) + "\n")
    return out


if __name__ == "__main__":
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    dst.mkdir(parents=True, exist_ok=True)
    for dat in sorted(src.glob("*.dat")):
        print(convert(dat, dst))
