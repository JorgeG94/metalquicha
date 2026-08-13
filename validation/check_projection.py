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
    base = 3 + n_lmo*n_lmo + 3*n_lmo + nao*n_lmo
    orbital_energies = np.array(values[base:base + n_occ])
    fock = np.array(values[at:at + n_lmo*n_lmo]).reshape((n_lmo, n_lmo), order="F")
    at += n_lmo*n_lmo
    centroids = np.array(values[at:at + 3*n_lmo]).reshape((3, n_lmo), order="F")
    at += 3*n_lmo
    orbitals = np.array(values[at:at + nao*n_lmo]).reshape((nao, n_lmo), order="F")
    return dict(basis=basis, nao=nao, n_occ=n_occ, n_lmo=n_lmo, fock=fock,
                centroids=centroids, orbitals=orbitals,
                orbital_energies=orbital_energies)


def gamess_primitive_norm(l, exponent):
    """The factor GAMESS folds into a printed contraction coefficient.

    The standard normalization of the `(l,0,0)` Cartesian primitive,

        N = (2a/pi)^(3/4) (4a)^(l/2) / sqrt((2l-1)!!)

    **Derived from the reference rather than assumed, and exact.** GAMESS prints
    1.11382493 where the Basis Set Exchange table has 1.0, and 0.83172368 where it
    has 0.0018311, so a `PROJECTION BASIS SET` written with our coefficients would
    be read by GAMESS as a different basis. Checked against all fifteen primitives
    of water/6-31G* -- six s, one sp triple, one sp single and one d, spanning
    exponents from 0.27 to 5484 -- and reproduces every printed coefficient to
    eight digits.

    This is what makes the section emittable. It is separate from, and not to be
    confused with, the d-shell factor in `D_NORMALIZATION`: that one relates the
    two codes' *basis functions*, and so appears in the orbital coefficients
    inversely.
    """
    double_factorial = 1.0
    n = 2*l - 1
    while n > 1:
        double_factorial *= n
        n -= 2
    return ((2.0*exponent/np.pi)**0.75 * (4.0*exponent)**(l/2.0)
            / np.sqrt(double_factorial))


def check_basis_normalization():
    """Every printed coefficient in the reference, from our table plus the rule."""
    s = parse_efp(REFERENCE.read_text())["sections"]
    raw = s["PROJECTION BASIS SET"]["raw"]

    # (l, exponent, printed coefficient) for every primitive GAMESS listed. An "L"
    # shell prints two coefficient columns, s then p.
    listed = []
    angular = None
    for line in raw:
        token = line.split()
        if not token:
            continue
        if token[0] in ("S", "P", "D", "L"):
            angular = token[0]
            continue
        if angular is None or not token[0].isdigit():
            continue
        numbers = _numbers(line)
        if len(numbers) < 3:
            continue
        exponent = numbers[1]
        if angular == "L":
            listed.append((0, exponent, numbers[2]))
            listed.append((1, exponent, numbers[3]))
        else:
            listed.append(({"S": 0, "P": 1, "D": 2}[angular], exponent, numbers[2]))

    import json
    table = json.loads((REPO / "basis_sets" / "6-31g_st_.json").read_text())
    ours = []
    for element in ("8", "1"):
        for shell in table["elements"][element]["electron_shells"]:
            for column, l in enumerate(shell["angular_momentum"]):
                for exponent, coefficient in zip(shell["exponents"],
                                                 shell["coefficients"][column]):
                    ours.append((l, float(exponent), float(coefficient)))

    worst = 0.0
    matched = 0
    for l, exponent, printed in listed:
        for our_l, our_exponent, our_coefficient in ours:
            if our_l == l and abs(our_exponent - exponent) < 1e-6*max(1.0, exponent):
                predicted = our_coefficient*gamess_primitive_norm(l, exponent)
                worst = max(worst, abs(predicted - printed)/max(1e-8, abs(printed)))
                matched += 1
                break
    return matched, len(listed), worst


def check_charge_transfer(orbital_energies):
    """`CTFOK` is the canonical occupied orbital energies, and nothing more.

    Identified by recognising the numbers: -20.5605, -1.3414, -0.7066, -0.5709,
    -0.4979 for water in 6-31G* are its five occupied canonical eigenvalues, not
    anything charge-transfer specific. So one of the two charge-transfer sections
    needs no new physics at all -- it is the SCF's own output, serialized.

    **`CTVEC` is reproduced too, and this docstring used to say it was not.** It has
    two forms. By default GAMESS writes the occupied orbitals plus a set of
    quasi-atomic valence virtuals, which for water in 6-31G* is 5 + 2 = 7 vectors --
    the seven that made the block look unidentifiable. Asked for the canonical form
    instead, it writes the whole canonical MO matrix, and that is exactly our
    canonical orbitals in its AO ordering, agreeing column by column to 3e-08 up to
    each column's arbitrary sign.

    The earlier conclusion here -- that DAF 15 holds neither our canonical nor our
    localized orbitals -- was reached by comparing against the default form's seven
    vectors, five of which *are* our occupied canonical orbitals to 2e-10. The
    comparison missed it because it was made in a foreign AO ordering, which is the
    same trap this file documents below for the overlap matrix. The two virtuals are
    genuinely different objects, and needing a quasi-atomic construction to make them
    is why the canonical form is the one worth emitting.
    """
    s = parse_efp(REFERENCE.read_text())["sections"]
    values = _numbers(" ".join(s["CTFOK"]["raw"]).replace(">", " "))
    n = len(values)
    gap = float(np.abs(np.array(values) - orbital_energies[:n]).max())
    return n, gap


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

    # **Pair the orbitals by centroid rather than assuming a shared order.**
    # Water's two lone pairs are mirror images, so they are degenerate under the
    # Boys functional and their order is arbitrary -- it follows the Jacobi sweep
    # path, which moves with the OpenMP reduction order in the SCF underneath. This
    # check originally compared elementwise and passed for several runs before the
    # pair came out swapped, reporting a 1.0 Bohr centroid disagreement on an
    # identical localization. `check_distributed_polarizability.py` pairs by
    # centroid for exactly this reason; this one now does too.
    order = []
    for j in range(ref_centroids.shape[1]):
        distances = [np.linalg.norm(ours["centroids"][:, i] - ref_centroids[:, j])
                     for i in range(ours["n_lmo"])]
        order.append(int(np.argmin(distances)))
    if sorted(order) != list(range(ours["n_lmo"])):
        print(f"  FAIL the centroid pairing is not one-to-one: {order}")
        return 1
    ours["centroids"] = ours["centroids"][:, order]
    ours["fock"] = ours["fock"][np.ix_(order, order)]
    ours["orbitals"] = ours["orbitals"][:, order]
    if order != list(range(ours["n_lmo"])):
        print(f"        our LMOs pair to GAMESS's CT1..CT{ours['n_lmo']} as "
              f"{[i + 1 for i in order]} -- degenerate lone pairs, so the order is "
              f"arbitrary and not a disagreement")

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

    # CTFOK, which turns out to be the SCF's own occupied eigenvalues.
    n_ct, ct_gap = check_charge_transfer(ours["orbital_energies"])
    print(f"        CTFOK: {n_ct} values, our occupied orbital energies to "
          f"{ct_gap:.2e}")
    if ct_gap > 1e-6:
        print("        FAIL CTFOK is not the canonical occupied orbital energies")
        failures += 1

    matched, total, worst = check_basis_normalization()
    print(f"        basis normalization: {matched}/{total} printed coefficients "
          f"reproduced from our table, worst relative gap {worst:.2e}")
    if matched != total or worst > 1e-7:
        print("        FAIL our basis table plus the normalization rule does not "
              "reproduce what GAMESS printed, so PROJECTION BASIS SET cannot be "
              "emitted from it")
        failures += 1

    print()
    if failures:
        print(f"[projection] {failures} FAILURE(S)")
        return 1
    print("[projection] our LMO Fock matrix and coefficients match GAMESS's MAKEFP")
    return 0


if __name__ == "__main__":
    sys.exit(main())
