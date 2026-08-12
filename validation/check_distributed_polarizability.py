#!/usr/bin/env python3
"""Compare our distributed LMO polarizabilities against a GAMESS MAKEFP potential.

    ./build/check_distributed_polarizability && \
        python3 validation/check_distributed_polarizability.py

Milestone 5, and the first check whose reference is GAMESS rather than PySCF.

**Two failures are possible and they are separated here.** A localization can put
its orbitals somewhere else than GAMESS's while each tensor is computed correctly,
or the centroids can agree while the tensors are wrong. So orbitals are paired by
nearest centroid first, that pairing is reported, and only then are the tensors
compared. A single number would confound the two.

**Boys localization has many local maxima, so the pairing is the honest part of
this.** Nothing guarantees our optimizer and GAMESS's climb to the same one; what
is checked is that they did, and if some future basis lands elsewhere the pairing
distance says so rather than the tensor comparison failing mysteriously.

Two invariants do *not* depend on the localization at all, and are checked
regardless: the sum over the retained orbitals, and its asymmetry. Both are
properties of the retained subspace only.
"""

import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tools" / "efp_validation"))
from efp_format import parse_efp  # noqa: E402

REFERENCE = REPO / "tools" / "efp_validation" / "reference" / "water_6-31gs_boys.efp"

#: `POLARIZABLE POINTS` component order, as (row, column) of `alpha_kl =
#: d mu_k / d F_l`.
#:
#: **The off-diagonals are the transpose of what GAMESS's own labels suggest, and
#: this was measured rather than read.** The source-derived order was
#: `xx yy zz xy xz yz yx zx zy`; taking that at face value put every per-orbital
#: tensor exactly one transpose away from ours -- the worst component gap came out
#: at 1.72e-01, which is precisely the tensor's own asymmetry, the signature of a
#: transpose and of nothing else. Reading the two off-diagonal triples the other
#: way round drops the gap to 2.0e-05, which is GAMESS's own CPHF precision. So
#: GAMESS's labels are (field, dipole) where ours are (dipole, field).
#:
#: Nothing else could have caught it. The sum over orbitals is symmetric, so the
#: M4 sum rule is blind to it; so is every isotropic average; so is the asymmetry
#: magnitude, which is transpose-invariant. It takes a per-orbital comparison
#: against a reference, which is exactly what the plan warned this milestone would
#: need.
POL_ORDER = [(0, 0), (1, 1), (2, 2), (1, 0), (2, 0), (2, 1), (0, 1), (0, 2), (1, 2)]

#: How far a paired centroid may sit from GAMESS's, in Bohr.
#:
#: Loose on purpose. Two Boys optimizers agreeing on a local maximum should agree
#: far better than this, so the value is not a precision claim -- it is the line
#: between "the same localization" and "a different one", and being generous makes
#: a genuine mismatch unmistakable rather than marginal.
CENTROID_TOL = 0.05

#: Agreement demanded on each tensor component, Bohr^3.
#:
#: The reference is printed to 10 decimals but GAMESS's own CPHF is the limit here,
#: not the printing: its water sum carries 4.6e-6 of asymmetry where ours carries
#: 3.7e-7, and a solver that leaves 5e-6 in a symmetry it should satisfy exactly
#: cannot be held to better than that in the tensors themselves.
TENSOR_TOL = 5e-5


def read_dump(path):
    tokens = path.read_text().split()
    basis = tokens[0]
    values = [float(x) for x in tokens[1:]]
    nao, n_occ, n_core = int(values[0]), int(values[1]), int(values[2])
    n_lmo = int(values[3])
    at = 4
    energy = values[at]; at += 1
    alpha = np.array(values[at:at + 9]).reshape((3, 3), order="F"); at += 9
    centroids = np.array(values[at:at + 3*n_lmo]).reshape((3, n_lmo), order="F")
    at += 3*n_lmo
    tensors = np.array(values[at:at + 9*n_lmo]).reshape((3, 3, n_lmo), order="F")
    return dict(basis=basis, nao=nao, n_occ=n_occ, n_core=n_core, n_lmo=n_lmo,
                energy=energy, alpha=alpha, centroids=centroids, tensors=tensors)


def reference_points():
    s = parse_efp(REFERENCE.read_text())["sections"]
    out = []
    for rec in s["POLARIZABLE POINTS"]:
        m = np.zeros((3, 3))
        for c, (i, j) in enumerate(POL_ORDER):
            m[i, j] = rec["tensor"][c]
        out.append(dict(label=rec["label"], xyz=np.array(rec["xyz"]), tensor=m))
    return out


def main():
    if not REFERENCE.exists():
        print(f"  MISSING {REFERENCE}")
        return 1
    ref = reference_points()
    failures = 0

    for index in (1, 2):
        path = pathlib.Path(f"/tmp/mqc_distributed_{index}.txt")
        if not path.exists():
            print(f"  MISSING {path} -- run ./build/check_distributed_polarizability")
            failures += 1
            continue
        ours = read_dump(path)

        total = ours["tensors"].sum(axis=2)
        asym = np.abs(total - total.T).max()
        print(f"  {ours['basis']}  n_core {ours['n_core']}  {ours['n_lmo']} LMOs   "
              f"sum alpha_iso {np.trace(total)/3:.6f}   sum asymmetry {asym:.2e}")

        # The whole occupied space must reproduce the total exactly; that is a
        # unit-tested identity, repeated here because it costs nothing and makes
        # the valence case's small asymmetry interpretable rather than alarming.
        if ours["n_core"] == 0:
            gap = np.abs(total - ours["alpha"]).max()
            print(f"        vs our own total alpha: {gap:.2e}   "
                  f"(exact identity, no localization dependence)")
            if gap > 1e-10:
                print("        FAIL the complete sum is not the total")
                failures += 1
            if asym > 1e-10:
                print("        FAIL the complete sum is not symmetric")
                failures += 1
            continue

        # Only the valence case is comparable to GAMESS, which excludes the core.
        if ours["n_lmo"] != len(ref):
            print(f"        FAIL we have {ours['n_lmo']} tensors, GAMESS has "
                  f"{len(ref)}")
            failures += 1
            continue

        ref_total = sum(p["tensor"] for p in ref)
        iso_gap = abs(np.trace(total)/3 - np.trace(ref_total)/3)
        print(f"        GAMESS sum alpha_iso {np.trace(ref_total)/3:.6f}   "
              f"asymmetry {np.abs(ref_total - ref_total.T).max():.2e}   "
              f"our gap {iso_gap:.2e}")

        # Pair by nearest centroid, and require the pairing to be a bijection --
        # two of our orbitals claiming the same reference point would otherwise
        # pass unnoticed and make the tensor comparison meaningless.
        taken = {}
        for i in range(ours["n_lmo"]):
            c = ours["centroids"][:, i]
            distances = [np.linalg.norm(c - p["xyz"]) for p in ref]
            j = int(np.argmin(distances))
            taken.setdefault(j, []).append((i, distances[j]))

        if len(taken) != len(ref):
            print("        FAIL the centroid pairing is not one-to-one; our "
                  "localization found a different solution")
            failures += 1
            continue

        worst_centroid = 0.0
        worst_tensor = 0.0
        for j, (pairs) in sorted(taken.items()):
            i, dist = pairs[0]
            worst_centroid = max(worst_centroid, dist)
            gap = np.abs(ours["tensors"][:, :, i] - ref[j]["tensor"]).max()
            worst_tensor = max(worst_tensor, gap)
            mine = ours["tensors"][:, :, i]
            print(f"        {ref[j]['label']} <- our LMO {i + 1}   "
                  f"centroid {dist:.2e}   tensor {gap:.2e}   "
                  f"iso ours {np.trace(mine)/3:8.5f} theirs "
                  f"{np.trace(ref[j]['tensor'])/3:8.5f}   "
                  f"asym ours {np.abs(mine - mine.T).max():.4f} theirs "
                  f"{np.abs(ref[j]['tensor'] - ref[j]['tensor'].T).max():.4f}")

        ok = True
        if worst_centroid > CENTROID_TOL:
            print(f"        FAIL a centroid is {worst_centroid:.3f} Bohr from "
                  f"GAMESS's, so the two localizations disagree")
            ok = False
        if worst_tensor > TENSOR_TOL:
            print(f"        FAIL a tensor differs by {worst_tensor:.2e}")
            ok = False
        if not ok:
            failures += 1

    print()
    if failures:
        print(f"[distributed] {failures} FAILURE(S)")
        return 1
    print("[distributed] our per-orbital polarizabilities match GAMESS's MAKEFP")
    return 0


if __name__ == "__main__":
    sys.exit(main())
