#!/usr/bin/env python3
"""Compare our LMO Fock matrix against a GAMESS MAKEFP potential.

    ./build/check_projection && python3 validation/check_projection.py

The Fock matrix in the localized basis is the one piece of a potential's
exchange-repulsion data that carries **no basis-function ordering**. The orbital
coefficients and the projection basis are indexed by AO, and libcint's ordering is
not GAMESS's, so those need a permutation worked out before they can be emitted;
this needs nothing but the LMOs themselves, which the polarizability comparison
already showed line up with GAMESS's CT1..CT4 in order.

It is also a sharper check on the localization than the centroids are. Two
localizations can put their orbitals in nearly the same places and still differ in
how much each one mixes with its neighbours, and it is the off-diagonal Fock
elements that say so -- they are what the exchange-repulsion energy is built from
at use time.
"""

import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tools" / "efp_validation"))
from efp_format import parse_efp  # noqa: E402

REFERENCE = REPO / "tools" / "efp_validation" / "reference" / "water_6-31gs_boys.efp"

#: The reference prints ten decimals and both codes converge to 1e-12, so this is
#: the file's precision rather than the physics'.
FOCK_TOL = 1e-8


def read_dump(path):
    tokens = path.read_text().split()
    basis = tokens[0]
    values = [float(x) for x in tokens[1:]]
    nao, n_occ, n_lmo = int(values[0]), int(values[1]), int(values[2])
    at = 3
    fock = np.array(values[at:at + n_lmo*n_lmo]).reshape((n_lmo, n_lmo), order="F")
    at += n_lmo*n_lmo
    centroids = np.array(values[at:at + 3*n_lmo]).reshape((3, n_lmo), order="F")
    return dict(basis=basis, nao=nao, n_occ=n_occ, n_lmo=n_lmo,
                fock=fock, centroids=centroids)


def reference_fock():
    s = parse_efp(REFERENCE.read_text())["sections"]
    raw = s["FOCK MATRIX ELEMENTS"]["raw"]
    values = [float(x) for line in raw for x in line.replace(">", " ").split()]
    # A lower triangle, row by row: (1,1) (2,1) (2,2) (3,1) ...
    n = int((np.sqrt(8*len(values) + 1) - 1)/2)
    if n*(n + 1)//2 != len(values):
        raise SystemExit(f"{len(values)} values is not a triangle")
    fock = np.zeros((n, n))
    k = 0
    for i in range(n):
        for j in range(i + 1):
            fock[i, j] = fock[j, i] = values[k]
            k += 1
    centroids = np.array([r["values"][:3] for r in s["LMO CENTROIDS"]]).T
    return fock, centroids


def main():
    path = pathlib.Path("/tmp/mqc_projection_1.txt")
    if not path.exists():
        print(f"  MISSING {path} -- run ./build/check_projection first")
        return 1
    ours = read_dump(path)
    theirs, ref_centroids = reference_fock()

    if ours["n_lmo"] != theirs.shape[0]:
        print(f"  FAIL we have {ours['n_lmo']} LMOs, GAMESS has {theirs.shape[0]}")
        return 1

    centroid_gap = np.abs(ours["centroids"] - ref_centroids).max()

    # **Up to orbital phase, which is not a convention we can or should match.**
    # A localized orbital is defined only to within a sign, and flipping orbital i
    # flips row and column i of the Fock matrix while leaving the diagonal alone.
    # Our LMOs 2 and 4 come out with the opposite phase to GAMESS's, so every
    # off-diagonal element agrees in magnitude and disagrees in sign, and an
    # elementwise comparison reports a 0.38 discrepancy on a matrix that is
    # physically identical.
    #
    # The exchange-repulsion energy is phase invariant: it is built from products
    # like `S_ij S_jk F_ki`, in which flipping orbital i changes two factors and
    # cancels. What is *not* optional is that the phases be consistent between the
    # Fock matrix and the orbital coefficients, since both are emitted -- so this
    # must never be "fixed" by flipping signs in one of them alone.
    #
    # So the phases are deduced from the first row and then required to explain
    # every other element. That is a stronger check than comparing magnitudes,
    # which would accept a matrix whose sign pattern was inconsistent.
    phases = np.ones(ours["n_lmo"])
    for j in range(1, ours["n_lmo"]):
        a, b = ours["fock"][0, j], theirs[0, j]
        if abs(b) > 1e-10:
            phases[j] = np.sign(a/b)
    aligned = ours["fock"]*np.outer(phases, phases)

    diag = np.abs(np.diag(aligned) - np.diag(theirs)).max()
    offdiag = np.abs(aligned - theirs)
    np.fill_diagonal(offdiag, 0.0)
    worst = np.abs(aligned - theirs).max()
    flipped = [i + 1 for i, p in enumerate(phases) if p < 0]

    print(f"  {ours['basis']}  {ours['n_lmo']} LMOs")
    print(f"        centroids agree to {centroid_gap:.2e} Bohr "
          f"(so the LMOs are paired in order, not just in a set)")
    print(f"        LMO Fock: diagonal {diag:.2e}   off-diagonal "
          f"{offdiag.max():.2e}   worst {worst:.2e}")
    print(f"        orbital phases opposite to GAMESS's: "
          f"{flipped if flipped else 'none'}  (physically irrelevant, see the "
          f"note in this file)")
    print("        ours, phase-aligned:")
    for row in aligned:
        print("          " + "".join(f"{v:14.8f}" for v in row))
    print("        GAMESS:")
    for row in theirs:
        print("          " + "".join(f"{v:14.8f}" for v in row))

    failures = 0
    if centroid_gap > 1e-8:
        print("        FAIL the LMOs are not in the same places, so comparing "
              "their Fock matrices elementwise is not meaningful")
        failures += 1
    if worst >= FOCK_TOL:
        print(f"        FAIL the LMO Fock matrices differ by {worst:.2e}")
        failures += 1

    print()
    if failures:
        print(f"[projection] {failures} FAILURE(S)")
        return 1
    print("[projection] our LMO Fock matrix matches GAMESS's MAKEFP")
    return 0


if __name__ == "__main__":
    sys.exit(main())
