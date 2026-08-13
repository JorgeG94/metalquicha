#!/usr/bin/env python3
"""Compare our dynamic polarizabilities against a GAMESS MAKEFP potential.

    ./build/check_dynamic && python3 validation/check_dynamic.py

Milestone 6. Three separate claims, checked separately because they fail for
different reasons:

1. **The frequencies.** `casimir_polder_frequencies` builds twelve points from
   Gauss-Legendre under `nu = 0.3 (1+t)/(1-t)`. The reference stamps each of its
   dynamic blocks with the frequency it used, so this is checkable exactly rather
   than by recognising the formula.
2. **The tensors**, summed over the localized orbitals at each frequency. GAMESS's
   dynamic points exclude the core, as its static ones do, so ours are compared
   with the same exclusion.
3. **The static limit**, which needs no reference at all and is checked in the
   Fortran: `nu = 0` must reproduce the static polarizability exactly.

The tolerance is GAMESS's precision rather than ours. Its static per-orbital
tensors agree with ours to 2e-5 and its own valence sum carries 4.6e-6 of an
asymmetry it should not have, so 1e-4 here is loose enough to accommodate that and
tight enough that a wrong frequency, a wrong projection or a wrong sign could not
pass.
"""

import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tools" / "efp_validation"))
from efp_format import parse_efp  # noqa: E402

REFERENCE = REPO / "tools" / "efp_validation" / "reference" / "water_6-31gs_boys.efp"

#: See the module docstring: this is GAMESS's precision, not ours.
TENSOR_TOL = 1e-4

#: How closely the quadrature points must match the ones GAMESS used. They are
#: printed to six decimals, so this is the printing.
FREQUENCY_TOL = 1e-6

POL_ORDER = [(0, 0), (1, 1), (2, 2), (1, 0), (2, 0), (2, 1), (0, 1), (0, 2), (1, 2)]


def read_dump(path):
    tokens = path.read_text().split()
    values = [float(x) for x in tokens]
    n_lmo, n_freq = int(values[0]), int(values[1])
    at = 3
    frequencies = np.array(values[at:at + n_freq]); at += n_freq
    centroids = np.array(values[at:at + 3*n_lmo]).reshape((3, n_lmo), order="F")
    at += 3*n_lmo
    tensors = np.array(values[at:at + 9*n_lmo*n_freq]).reshape(
        (3, 3, n_lmo, n_freq), order="F")
    return dict(n_lmo=n_lmo, n_freq=n_freq, frequencies=frequencies,
                centroids=centroids, tensors=tensors)


def reference_blocks():
    """GAMESS's dynamic points, grouped by the frequency stamped on each block."""
    points = parse_efp(REFERENCE.read_text())["sections"]["DYNAMIC POLARIZABLE POINTS"]
    frequencies, blocks = [], []
    for point in points:
        if "frequency" in point:
            frequencies.append(point["frequency"])
            blocks.append([])
        blocks[-1].append(point)
    return np.array(frequencies), blocks


def unpack(values):
    m = np.zeros((3, 3))
    for c, (i, j) in enumerate(POL_ORDER):
        m[i, j] = values[c]
    return m


def main():
    path = pathlib.Path("/tmp/mqc_dyn.txt")
    if not path.exists():
        print(f"  MISSING {path} -- run ./build/check_dynamic first")
        return 1
    ours = read_dump(path)
    ref_freq, blocks = reference_blocks()
    failures = 0

    # Our dump leads with nu = 0, which GAMESS does not tabulate; drop it for the
    # comparison and use it for the static-limit check instead.
    static = ours["tensors"][:, :, :, 0]
    dynamic = ours["tensors"][:, :, :, 1:]
    our_freq = ours["frequencies"][1:]

    # Pair by centroid, not by position. Water's two lone pairs are degenerate
    # under the Boys functional -- mirror images -- so which comes first follows the
    # Jacobi sweep path and moves with the OpenMP reduction order in the SCF
    # underneath. Comparing per-orbital tensors by index passed for several runs and
    # then failed on an identical localization when the pair came out swapped.
    ref_centroids = np.array([p["xyz"] for p in blocks[0]]).T
    order = []
    for j in range(ref_centroids.shape[1]):
        distances = [np.linalg.norm(ours["centroids"][:, i] - ref_centroids[:, j])
                     for i in range(ours["n_lmo"])]
        order.append(int(np.argmin(distances)))
    if sorted(order) != list(range(ours["n_lmo"])):
        print(f"  FAIL the centroid pairing is not one-to-one: {order}")
        return 1
    static = static[:, :, order]
    dynamic = dynamic[:, :, order, :]
    if order != list(range(ours["n_lmo"])):
        print(f"        our LMOs pair to CT1..CT{ours['n_lmo']} as "
              f"{[i + 1 for i in order]} -- degenerate lone pairs, arbitrary order")

    print(f"  {ours['n_lmo']} localized orbitals, {len(ref_freq)} frequencies")
    gap = np.abs(our_freq - ref_freq).max()
    print(f"        quadrature points agree with GAMESS's to {gap:.2e}")
    if gap > FREQUENCY_TOL:
        print("        FAIL the frequencies are not the ones GAMESS tabulated at")
        failures += 1
        return 1

    print(f"        {'nu':>10} {'GAMESS iso':>13} {'our iso':>13} {'worst elem':>11}")
    worst_overall = 0.0
    for k, nu in enumerate(ref_freq):
        theirs = np.array([unpack(p["tensor"]) for p in blocks[k]])
        mine = np.array([dynamic[:, :, i, k] for i in range(ours["n_lmo"])])
        worst = np.abs(mine - theirs).max()
        worst_overall = max(worst_overall, worst)
        print(f"        {nu:10.6f} {np.trace(theirs.sum(0))/3:13.8f} "
              f"{np.trace(mine.sum(0))/3:13.8f} {worst:11.2e}")

    print(f"        worst element over all frequencies and orbitals: "
          f"{worst_overall:.2e}")
    if worst_overall >= TENSOR_TOL:
        print("        FAIL the dynamic tensors do not match GAMESS's")
        failures += 1

    # The static limit, orbital by orbital, against the same projection at nu = 0.
    # An exact identity: the frequency-dependent equations reduce to the static
    # ones, so anything but agreement means the two projections disagree.
    print(f"        static limit: sum_iso at nu=0 is "
          f"{np.trace(static.sum(2))/3:.8f}")

    print()
    if failures:
        print(f"[dynamic] {failures} FAILURE(S)")
        return 1
    print("[dynamic] our dynamic polarizabilities match GAMESS's MAKEFP")
    return 0


if __name__ == "__main__":
    sys.exit(main())
