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
from efp_format import parse_efp, _numbers  # noqa: E402

REFERENCE = REPO / "tools" / "efp_validation" / "reference" / "water_6-31gs_boys.efp"

#: The reference prints ten decimals and both codes converge to 1e-12, so this is
#: the file's precision rather than the physics'.
FOCK_TOL = 1e-8

#: libcint's basis-function index for each of GAMESS's six Cartesian d slots.
#:
#: **The shell ordering already agrees and only the d shell differs.** libcint gets
#: water/6-31G* as ten shells -- s(6), s(3), p(3), s(1), p(1), d(1) on oxygen, then
#: s(3), s(1) per hydrogen -- and GAMESS as S(6), L(3), L(1), D(1) then S(3), S(1).
#: Different shell *types*, identical basis-function positions, because our split
#: s and p shells land exactly where GAMESS's combined L shell puts its s and p.
#: So no atom or shell permutation is needed at all; verified by the s and p
#: coefficients agreeing to 1.4e-9 before any reordering.
#:
#: Within the d shell libcint runs `xx xy xz yy yz zz` and GAMESS `XX YY ZZ XY XZ
#: YZ`.
D_FROM_LIBCINT = [0, 3, 5, 1, 2, 4]

#: And the Cartesian d normalization, which differs as well as the order.
#:
#: GAMESS's diagonal components carry `sqrt(pi a)` relative to libcint's and its
#: off-diagonal ones `sqrt(pi a / 3)` -- the ratio being exactly sqrt(3), which is
#: the familiar difference between normalising every Cartesian component
#: individually and normalising the `(l,0,0)` one. Measured here as 1.585330892 and
#: 0.915290921 on the oxygen d shell at exponent 0.8, whose square is 0.8*pi to
#: eight digits.
#:
#: **The sqrt(3) ratio is established; the overall sqrt(pi a) is inferred from one
#: shell.** A basis with more than one d exponent, or with f functions, would test
#: it -- and emitting `PROJECTION BASIS SET` needs the general rule anyway, since
#: GAMESS prints contraction coefficients in its own normalisation (1.11382493 for
#: a BSE coefficient of 1.0 at that exponent), so a basis emitted with ours would
#: be read as a different basis.
D_NORMALIZATION = 1.585330892


def read_dump(path):
    tokens = path.read_text().split()
    basis = tokens[0]
    values = [float(x) for x in tokens[1:]]
    nao, n_occ, n_lmo = int(values[0]), int(values[1]), int(values[2])
    at = 3
    fock = np.array(values[at:at + n_lmo*n_lmo]).reshape((n_lmo, n_lmo), order="F")
    at += n_lmo*n_lmo
    centroids = np.array(values[at:at + 3*n_lmo]).reshape((3, n_lmo), order="F")
    at += 3*n_lmo
    orbitals = np.array(values[at:at + nao*n_lmo]).reshape((nao, n_lmo), order="F")
    return dict(basis=basis, nao=nao, n_occ=n_occ, n_lmo=n_lmo,
                fock=fock, centroids=centroids, orbitals=orbitals)


def reference_coefficients(n_lmo, nao):
    """GAMESS's localized orbital coefficients, one row per LMO."""
    s = parse_efp(REFERENCE.read_text())["sections"]
    values = []
    for line in s["PROJECTION WAVEFUNCTION"]["raw"]:
        # Past the "  i  j" index pair. The numbers run together where a minus
        # sign abuts the previous exponent, so they cannot simply be split.
        values += _numbers(line[5:])
    return np.array(values).reshape((n_lmo, nao))


def to_gamess_order(coefficients, d_offset):
    """Our coefficients in GAMESS's basis-function order and normalization."""
    out = coefficients.copy()
    scale = [D_NORMALIZATION]*3 + [D_NORMALIZATION/np.sqrt(3.0)]*3
    for slot, (src, factor) in enumerate(zip(D_FROM_LIBCINT, scale)):
        out[d_offset + slot, :] = coefficients[d_offset + src, :]*factor
    return out


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

    # The orbital coefficients, which is what needs the basis-function map. The
    # phase is taken from the s and p block alone, so the d block is compared
    # against a phase it did not help choose.
    theirs_c = reference_coefficients(ours["n_lmo"], ours["nao"])
    mine = to_gamess_order(ours["orbitals"], d_offset=9)
    sp = [i for i in range(ours["nao"]) if not 9 <= i < 15]
    worst_sp = worst_d = 0.0
    for i in range(ours["n_lmo"]):
        phase = 1.0 if np.dot(ours["orbitals"][sp, i], theirs_c[i][sp]) >= 0 else -1.0
        gap = np.abs(phase*mine[:, i] - theirs_c[i])
        worst_sp = max(worst_sp, gap[sp].max())
        worst_d = max(worst_d, gap[9:15].max())
    print(f"        coefficients: s/p {worst_sp:.2e} before any reordering, "
          f"d {worst_d:.2e} after the permutation and normalization")
    if max(worst_sp, worst_d) >= 1e-7:
        print("        FAIL the localized orbital coefficients do not map onto "
              "GAMESS's basis-function order")
        failures += 1

    print()
    if failures:
        print(f"[projection] {failures} FAILURE(S)")
        return 1
    print("[projection] our LMO Fock matrix and coefficients match GAMESS's MAKEFP")
    return 0


if __name__ == "__main__":
    sys.exit(main())
