"""Compare the Fortran Treutler-Ahlrichs mesh against PySCF's, node by node.

    compare_radial.py radial_mesh.txt

Unlike the Lebedev comparison this one is element-wise: both implementations
produce nodes in ascending order, so there is no ordering ambiguity to sort
away, and a mismatch in position is exactly what we want to catch.
"""
import sys
from collections import defaultdict

import numpy as np
from pyscf.dft import radi

rows = defaultdict(list)
for line in open(sys.argv[1]):
    z, n, r, dr = line.split()
    rows[(int(z), int(n))].append((float(r), float(dr)))

print(f"{'Z':>4} {'n':>5} {'max |dr_node|':>15} {'max rel |dw|':>14} {'tol':>11}   ")
worst_r = worst_w = 0.0
bad = []

for (z, n), vals in sorted(rows.items()):
    ours = np.array(vals)
    r_ref, dr_ref = radi.treutler_ahlrichs(n, z)

    dr_node = np.abs(ours[:, 0] - r_ref).max()
    # Weights span many orders of magnitude across the mesh, so compare
    # relatively; an absolute difference would be dominated by the outermost
    # shell and hide a bad inner one.
    rel_w = (np.abs(ours[:, 1] - dr_ref) / np.abs(dr_ref)).max()

    # The M4 mapping is worst-conditioned at its outer end, where x -> 1 and
    # the weight carries a 1/(1-x) term. A one-ulp difference in cos between
    # gfortran's libm and glibc's -- which is all that separates the two
    # implementations -- is amplified there by eps/(1-cos(pi/(n+1))). That is
    # the floor this comparison can reach, and it grows with n, so the
    # tolerance is derived from it rather than fixed. The factor of 10 is
    # headroom for the couple of ulps the rest of the expression contributes.
    x_outer = np.cos(np.pi / (n + 1))
    floor = np.spacing(abs(x_outer)) / abs(1 - x_outer)
    tol_w = 10 * floor

    worst_r = max(worst_r, dr_node)
    worst_w = max(worst_w, rel_w)
    ok = dr_node < 1e-11 and rel_w < tol_w
    if not ok:
        bad.append((z, n, dr_node, rel_w))
    print(f"{z:>4} {n:>5} {dr_node:>15.3e} {rel_w:>14.3e} {tol_w:>11.3e}   {'ok' if ok else 'FAIL'}")

print()
print(f"worst node deviation   : {worst_r:.3e} Bohr")
print(f"worst relative weight  : {worst_w:.3e}")
if bad:
    print("\nFAILURES:", bad)
    sys.exit(1)
print("\nevery node and weight matches PySCF")
