#!/usr/bin/env python3
"""Compare our electrostatic potential integrals against PySCF's.

    ./build/check_esp && python3 validation/check_esp.py

Elementwise over every AO pair at every point, not on a contracted potential.
libcint's grid integral nests the grid index *fastest*, inside the shell pair --
the opposite of its multipole blocks, where the component is slowest -- so a
transposed unpack produces a matrix that is right at one point and wrong at the
rest. Contracting first would average that into a single number that could easily
look plausible.
"""

import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tools" / "cpu_validation"))
from gen_cpu_validation import bse_to_pyscf, molecule_form, CARTESIAN  # noqa: E402

ATOMS = [("O", (0.0, 0.0, 0.0)),
         ("H", (0.0, 0.0, 0.9584)),
         ("H", (0.9268, 0.0, -0.2400))]

#: **Do not add a Pople basis with sp shells to this comparison without permuting.**
#: PySCF orders a Cartesian 6-31G* oxygen as `1s 2s 3s 2px..2pz 3px..3pz d`, grouping
#: all s functions then all p. Our reader hands libcint the s and p of each shared-
#: exponent group adjacently, giving `1s 2s 2px..2pz 3s 3px..3pz d`, which is also
#: GAMESS's order. So an elementwise comparison of AO-indexed matrices against PySCF
#: is only valid for a basis without sp shells -- the cases here are all like that,
#: which is why they pass.
#:
#: This was found the hard way: using PySCF's overlap with our own localized orbitals
#: reported them as non-orthonormal (`c^T S c` = 0.97, 0.75, 0.74, 0.43) when
#: `check_localize` verifies orthonormality to 2e-15 against our own overlap. The
#: metric was in the wrong order, not the orbitals.

#: Both codes call the same libcint routine on the same basis, so this is the
#: floor set by our BSE table against PySCF's own, not a physics tolerance.
TOL = 1e-8


def main():
    from pyscf import gto

    path = pathlib.Path("/tmp/mqc_esp.txt")
    if not path.exists():
        print(f"  MISSING {path} -- run ./build/check_esp first")
        return 1
    values = [float(x) for x in path.read_text().split()]
    nao, n_grid = int(values[0]), int(values[1])
    at = 2
    grids = np.array(values[at:at + 3*n_grid]).reshape((3, n_grid), order="F")
    at += 3*n_grid
    ours = np.array(values[at:at + nao*nao*n_grid]).reshape((nao, nao, n_grid),
                                                            order="F")

    symbols = {a[0] for a in ATOMS}
    mol = gto.Mole()
    mol.atom = ATOMS
    mol.unit = "Angstrom"
    mol.basis = {s: bse_to_pyscf("cc-pvdz", s) for s in symbols}
    mol.cart = molecule_form("cc-pvdz", symbols) == CARTESIAN
    mol.verbose = 0
    mol.build()
    theirs = mol.intor("int1e_grids", grids=grids.T)      # (n_grid, nao, nao)

    if theirs.shape != (n_grid, nao, nao):
        print(f"  FAIL PySCF returned {theirs.shape}, expected "
              f"{(n_grid, nao, nao)}")
        return 1

    worst = max(np.abs(ours[:, :, g] - theirs[g]).max() for g in range(n_grid))
    print(f"  cc-pvdz  nao {nao}  points {n_grid}")
    print(f"        worst element over {n_grid*nao*nao} integrals: {worst:.2e}")
    for g in range(n_grid):
        gap = np.abs(ours[:, :, g] - theirs[g]).max()
        print(f"        ({grids[0, g]:6.2f},{grids[1, g]:6.2f},{grids[2, g]:6.2f})"
              f"   {gap:.2e}")

    print()
    if worst >= TOL:
        print("[esp] the ESP integrals do not match PySCF's")
        return 1
    print("[esp] our electrostatic potential integrals match PySCF's")
    return 0


if __name__ == "__main__":
    sys.exit(main())
