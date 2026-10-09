#!/usr/bin/env python3
"""PySCF's minao guess density for the molecules test_mqc_czt_minao pins.

The orbital basis is read from this repository's own BSE JSON rather than
PySCF's tables, which differ in the last digits on some sets. The files used
are already sorted by angular momentum, so PySCF's AO order is mqc's.

    python3 tools/minao/minao_reference.py            # the pinned fingerprints
    python3 tools/minao/minao_reference.py --dump X   # also write each full D to X_<name>.npy
"""
import argparse
import json
import pathlib

import numpy as np
from pyscf import gto
from pyscf.scf import hf

ROOT = pathlib.Path(__file__).resolve().parents[2]
ANG = 1.8897261254578281

MOLECULES = {
    # name: (basis, [(Z, symbol, x, y, z in Angstrom)])
    "water_ccpvdz": ("cc-pvdz", [
        (8, "O", 0.0, 0.0, 0.0),
        (1, "H", 0.0, 0.0, 0.9584),
        (1, "H", 0.9268, 0.0, -0.2400)]),
    "methanol_def2svp": ("def2-svp", [
        (6, "C", -0.0467, 0.6635, 0.0000),
        (8, "O", -0.0467, -0.7581, 0.0000),
        (1, "H", -1.0868, 0.9790, 0.0000),
        (1, "H", 0.4377, 1.0739, 0.8920),
        (1, "H", 0.4377, 1.0739, -0.8920),
        (1, "H", 0.8752, -1.0659, 0.0000)]),
}


def bse_to_pyscf(name, symbols):
    data = json.loads((ROOT / "basis_sets" / (name + ".json")).read_text())
    out = {}
    for z, sym in symbols:
        shells = []
        for sh in data["elements"][str(z)]["electron_shells"]:
            exps = [float(e) for e in sh["exponents"]]
            cols = [[float(c) for c in col] for col in sh["coefficients"]]
            moms = sh["angular_momentum"]
            if len(moms) == 1:
                shells.append([moms[0]] + [[e] + [c[i] for c in cols]
                                           for i, e in enumerate(exps)])
            else:
                for l, col in zip(moms, cols):
                    shells.append([l] + [[e, col[i]] for i, e in enumerate(exps)])
        ls = [s[0] for s in shells]
        assert ls == sorted(ls), (name, sym, "not l-sorted; AO order would differ")
        out[sym] = shells
    return out


def molecule(name):
    basis, atoms = MOLECULES[name]
    syms = sorted({(a[0], a[1]) for a in atoms})
    return gto.M(atom=[(a[1], (a[2], a[3], a[4])) for a in atoms],
                 basis=bse_to_pyscf(basis, syms), unit="Angstrom",
                 cart=False, verbose=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump")
    args = ap.parse_args()
    for name in MOLECULES:
        mol = molecule(name)
        dm = np.asarray(hf.init_guess_by_minao(mol))
        s = mol.intor("int1e_ovlp")
        print(name, "nao", mol.nao)
        print("  tr(DS)   %.12f" % np.einsum("ij,ji", dm, s))
        print("  |D|_F    %.12f" % np.linalg.norm(dm))
        for i, j in [(0, 0), (1, 1), (0, 1), (2, 5), (mol.nao - 1, mol.nao - 1), (3, mol.nao - 2)]:
            print("  D(%d,%d) %.12f" % (i + 1, j + 1, dm[i, j]))
        if args.dump:
            np.save("%s_%s.npy" % (args.dump, name), dm)


if __name__ == "__main__":
    main()
