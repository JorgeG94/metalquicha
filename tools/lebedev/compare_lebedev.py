"""Compare the Fortran Lebedev grids against PySCF's, order by order.

The two implementations enumerate points in different orders -- ours expands
orbits by permutation, theirs lists them explicitly -- so this matches the two
as sets: sort both canonically, then compare element-wise at full precision.

A set match is the stronger test anyway. Enumeration order is an implementation
detail; which points exist and what weight each carries is the physics.

    compare_lebedev.py lebedev_points.txt lebedev.json
"""
import json
import sys

import numpy as np


def canonical(points, weights):
    """Sort (x, y, z, w) rows into a reproducible order."""
    rows = np.column_stack([points, weights])
    # Round only for the sort key, so near-ties order consistently in both
    # sets; the comparison itself uses the unrounded values.
    key = np.round(rows, 10)
    idx = np.lexsort((key[:, 3], key[:, 2], key[:, 1], key[:, 0]))
    return rows[idx]


def main(fortran_path, reference_path):
    ours = {}
    for line in open(fortran_path):
        parts = line.split()
        order = int(parts[0])
        ours.setdefault(order, []).append([float(x) for x in parts[1:]])

    theirs = json.load(open(reference_path))["reference"]

    print(f"{'order':>6} {'points':>7} {'max |dx|':>12} {'max |dw|':>12}  {'':>6}")
    worst_xyz = worst_w = 0.0
    failures = []

    for order in sorted(ours):
        a = np.array(ours[order])
        b = np.array(theirs[str(order)])

        if a.shape != b.shape:
            failures.append(f"order {order}: shape {a.shape} vs {b.shape}")
            continue

        ca = canonical(a[:, :3], a[:, 3])
        cb = canonical(b[:, :3], b[:, 3])

        dxyz = np.abs(ca[:, :3] - cb[:, :3]).max()
        dw = np.abs(ca[:, 3] - cb[:, 3]).max()
        worst_xyz = max(worst_xyz, dxyz)
        worst_w = max(worst_w, dw)

        ok = dxyz < 1e-13 and dw < 1e-15
        if not ok:
            failures.append(f"order {order}: max |dx| {dxyz:.3e}, max |dw| {dw:.3e}")
        print(f"{order:>6} {len(a):>7} {dxyz:>12.3e} {dw:>12.3e}  {'ok' if ok else 'FAIL':>6}")

    print()
    print(f"orders compared : {len(ours)}")
    print(f"points compared : {sum(len(v) for v in ours.values())}")
    print(f"worst coordinate deviation : {worst_xyz:.3e}")
    print(f"worst weight deviation     : {worst_w:.3e}")

    if failures:
        print("\nFAILURES:")
        for f in failures:
            print("  " + f)
        return 1
    print("\nevery point and weight matches PySCF")
    return 0


sys.exit(main(sys.argv[1], sys.argv[2]))
