"""Compare the assembled molecular grid against PySCF's, and integrate on both.

    compare_grid.py <check_grid binary dir>

Two comparisons. The grids themselves are matched as sets -- both codes loop
atoms then radial then angular, but our Lebedev points come out in orbit-
expansion order rather than tabulated order, so positions within a sphere
differ even though the sets agree.

The integrals are the calibration that matters: they say what accuracy a given
level actually buys, which is what a tolerance in the Fortran check should be
set from rather than guessed.
"""
import sys

import numpy as np
from pyscf import gto
from pyscf.dft import gen_grid, radi

MOLECULES = {
    "water": ([8, 1, 1], [(0.0, 0.0, 0.0), (0.0, -1.4308, 1.1078), (0.0, 1.4308, 1.1078)]),
    "auh": ([79, 1], [(0.0, 0.0, 0.0), (0.0, 0.0, 2.8)]),
}
GAUSS_ALPHA, SLATER_ZETA = 1.3, 1.7


def pyscf_grid(label, level):
    charges, coords = MOLECULES[label]
    spec = [(f"{gto.mole._std_symbol(gto.mole._symbol(int(c)))}{k+1}", xyz)
            for k, (c, xyz) in enumerate(zip(charges, coords))]
    dummy = {sym: [[0, [1.0, 1.0]]] for sym, _ in spec}
    mol = gto.M(atom=spec, unit="Bohr", basis=dummy, spin=None, verbose=0)

    g = gen_grid.Grids(mol)
    g.level = level
    g.prune = None                                   # we do not prune yet
    g.radi_method = radi.treutler_ahlrichs
    g.becke_scheme = gen_grid.original_becke
    g.radii_adjust = radi.treutler_atomic_radii_adjust
    g.alignment = 0                                  # no padding points
    g.build()
    return mol, g.coords, g.weights


def integrals(coords, weights, atom_coords):
    gnorm = (GAUSS_ALPHA / np.pi) ** 1.5
    snorm = SLATER_ZETA ** 3 / np.pi
    fg = fs = 0.0
    for R in atom_coords:
        d = coords - R
        r2 = np.einsum("ij,ij->i", d, d)
        fg = fg + gnorm * np.exp(-GAUSS_ALPHA * r2)
        fs = fs + snorm * np.exp(-2 * SLATER_ZETA * np.sqrt(r2))
    return (weights * fg).sum(), (weights * fs).sum()


print("PySCF, same settings, unpruned -- what each level actually buys:")
print(f"{'molecule':>8} {'level':>6} {'points':>9} {'gaussians':>14} {'slaters':>14}")
for label in MOLECULES:
    for level in range(1, 6):
        mol, coords, weights = pyscf_grid(label, level)
        exact = mol.natm
        fg, fs = integrals(coords, weights, mol.atom_coords())
        print(f"{label:>8} {level:>6} {len(weights):>9} "
              f"{abs(fg-exact)/exact:>14.6e} {abs(fs-exact)/exact:>14.6e}")

print()
print("full point-set comparison at level 3 (ours vs PySCF):")
print(f"{'molecule':>8} {'points':>9} {'max |dx|':>12} {'max rel |dw|':>14}   ")

griddir = sys.argv[1] if len(sys.argv) > 1 else "."
failed = []
for label in MOLECULES:
    path = f"{griddir}/grid_{label.strip()}.txt"
    ours = np.loadtxt(path)
    mol, coords, weights = pyscf_grid(label, 3)

    if len(ours) != len(weights):
        print(f"{label:>8} {len(ours):>9}   count mismatch: pyscf has {len(weights)}")
        failed.append(label)
        continue

    # Both codes loop atoms then radial then angular, but our Lebedev points
    # come out in orbit-expansion order, so sort before comparing.
    def canon(xyz, w):
        rows = np.column_stack([xyz, w])
        key = np.round(rows, 9)
        return rows[np.lexsort((key[:, 3], key[:, 2], key[:, 1], key[:, 0]))]

    a, b = canon(ours[:, :3], ours[:, 3]), canon(coords, weights)
    dx = np.abs(a[:, :3] - b[:, :3]).max()
    # Weights are compared against the grid's own scale, not point by point.
    # A pointwise relative test is meaningless here: far from the molecule the
    # Becke cell product underflows, and PySCF's C routine flushes to exactly
    # 0.0 where our Fortran keeps a denormal-scale value around 1e-16. Both are
    # zero for any purpose a quadrature has -- dividing by theirs is what
    # produces an infinity, not any disagreement about the grid.
    dw = np.abs(a[:, 3] - b[:, 3]).max() / np.abs(b[:, 3]).max()

    n_zero_ours, n_zero_theirs = int((a[:, 3] == 0).sum()), int((b[:, 3] == 0).sum())

    ok = dx < 1e-10 and dw < 1e-12
    if not ok:
        failed.append(label)
    print(f"{label:>8} {len(ours):>9} {dx:>12.3e} {dw:>14.3e}   {'ok' if ok else 'FAIL'}"
          f"   (underflowed to exactly 0: ours {n_zero_ours}, pyscf {n_zero_theirs})")

print()
if failed:
    print("FAILED:", failed)
    sys.exit(1)
print("every grid point and weight matches PySCF")
