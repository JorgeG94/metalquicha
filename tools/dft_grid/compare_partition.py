"""Compare the Fortran Becke partition against PySCF's, point by point.

    compare_partition.py partition_weights.txt

PySCF's get_partition folds the partition into the atomic grid volume, so to
isolate it we hand it our probe points with a volume weight of 1: what comes
back is then the partition weight alone. Atoms are given distinct labels
(O1, H2, H3) so each gets its own entry in the atom grid table and we can read
off w_A at every probe for every A, rather than only for the owning atom.

Note which PySCF code path each case exercises. original_becke is a marker that
routes to the C routine VXCgen_grid; stratmann falls through to the Python
loop. So this compares against both of their implementations, not one.
"""
import sys
from collections import defaultdict

import numpy as np
from pyscf import gto
from pyscf.dft import gen_grid, radi

MOLECULES = {
    "water": ([8, 1, 1], [(0.0, 0.0, 0.0), (0.0, -1.4308, 1.1078), (0.0, 1.4308, 1.1078)]),
    "auh": ([79, 1], [(0.0, 0.0, 0.0), (0.0, 0.0, 2.8)]),
}
SCHEME = {1: gen_grid.original_becke, 2: gen_grid.stratmann}
ADJUST = {0: None, 1: radi.becke_atomic_radii_adjust, 2: radi.treutler_atomic_radii_adjust}
SCHEME_NAME = {1: "Becke", 2: "Stratmann"}
ADJUST_NAME = {0: "none", 1: "becke", 2: "treutler"}

rows = defaultdict(dict)
for line in open(sys.argv[1]):
    label, scheme, adjust, ia, x, y, z, w = line.split()
    rows[(label, int(scheme), int(adjust), int(ia))][(float(x), float(y), float(z))] = float(w)

print(f"{'molecule':>9} {'scheme':>10} {'adjust':>9} {'atom':>5} {'max |dw|':>12}   ")
worst = 0.0
bad = []

for (label, scheme, adjust, ia), pts in sorted(rows.items()):
    charges, coords = MOLECULES[label]
    atom_spec = [(f"{gto.mole._std_symbol(gto.mole._symbol(int(c)))}{k+1}", xyz)
                 for k, (c, xyz) in enumerate(zip(charges, coords))]
    # The partition depends only on nuclear positions and Z, never on the
    # basis -- but gto.M insists on one, and sto-3g has no gold. A single s
    # function per atom keeps every element constructible.
    dummy = {sym: [[0, [1.0, 1.0]]] for sym, _ in atom_spec}
    mol = gto.M(atom=atom_spec, unit="Bohr", basis=dummy, spin=None, verbose=0)

    probes = np.array(list(pts.keys()))
    ours = np.array(list(pts.values()))

    atm = mol.atom_coords()
    tab = {mol.atom_symbol(k): (probes - atm[k], np.ones(len(probes)))
           for k in range(mol.natm)}

    _, weights = gen_grid.get_partition(mol, tab, radii_adjust=ADJUST[adjust],
                                        becke_scheme=SCHEME[scheme], concat=False)
    theirs = np.asarray(weights[ia - 1])

    d = np.abs(ours - theirs).max()
    worst = max(worst, d)
    ok = d < 1e-13
    if not ok:
        bad.append((label, scheme, adjust, ia, d))
    print(f"{label:>9} {SCHEME_NAME[scheme]:>10} {ADJUST_NAME[adjust]:>9} {ia:>5} {d:>12.3e}   "
          f"{'ok' if ok else 'FAIL'}")

print()
print(f"worst partition-weight deviation : {worst:.3e}")
if bad:
    print("\nFAILURES:")
    for b in bad:
        print("  ", b)
    sys.exit(1)
print("\nevery partition weight matches PySCF")
