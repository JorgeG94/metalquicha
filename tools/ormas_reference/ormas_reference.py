"""Exact reference energies for occupation-restricted active spaces.

PySCF has no ORMAS and no RASSCF, so there is nothing to call. What there is,
for a space small enough, is something better than a call: the determinants can
be enumerated, the Hamiltonian built in that basis by the Slater-Condon rules,
and the whole thing diagonalised. That is exact -- not iterative, not truncated
-- and it works for any partition rather than only the ones some other program
ships an input deck for.

Nothing here is shared with the Fortran under test. Determinants are spin-orbital
sets in whatever order Python produces them; there are no occupation classes, no
compatibility grid, no addressing. That is deliberate: an energy does not care
how the determinants were numbered, so a reference that numbers them differently
checks the space and the physics without being able to agree by accident about
the layout.

The reference needs its own reference. With windows that exclude nothing, an
ORMAS space *is* the full space, so `pyscf.fci` is the check -- and it is a sharp
one, because a sign error in the off-diagonal Slater-Condon rules survives every
internal consistency test and shows up immediately here.

Run this file to see both: the self-check against PySCF, then reference energies
for the partitions the Fortran tests use.
"""

from __future__ import annotations

import itertools
import math

import numpy as np

NCHOL = 3


def model_integrals(norb: int, dominant: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """A Hamiltonian with the right symmetries and no physics.

    ``h_pq = -1/(p+q)`` and ``(pq|rs) = sum_L B_pqL B_rsL`` with
    ``B_pqL = 1/(p+q+L)``, on zero-based indices. Factorising rather than
    filling entries at random gives ``(pq|rs)`` the eightfold permutational
    symmetry a real integral has.

    The same formula is built in ``test/test_mqc_ormas_ci.f90``. Every value is
    the reciprocal of a small integer, so the two languages agree bit for bit
    and a number computed here can be transcribed there without rounding.
    """
    p = np.arange(norb)
    h1e = -1.0 / (p[:, None] + p[None, :] + 2.0)
    if dominant:
        # The plain model puts every state within a whisker of every other, so
        # a determinant's diagonal says almost nothing about where its energy
        # sits. That is fine for checking a Hamiltonian and useless for
        # checking a solver that preconditions with exactly that diagonal.
        # Separating the orbitals in energy separates the determinants, which
        # is the regime any real CI is in.
        np.fill_diagonal(h1e, -2.0 * (norb - p))

    b = np.empty((norb, norb, NCHOL))
    for lcho in range(NCHOL):
        b[:, :, lcho] = 1.0 / (p[:, None] + p[None, :] + lcho + 3.0)
    eri = np.einsum("pql,rsl->pqrs", b, b)
    return h1e, eri


def subspace_widths(first_orbital: list[int], norb: int) -> list[int]:
    """Orbitals in each subspace, from the orbital each one starts at."""
    bounds = list(first_orbital) + [norb + 1]
    return [bounds[k + 1] - bounds[k] for k in range(len(first_orbital))]


def occupation(orbitals: frozenset[int], first_orbital: list[int], norb: int) -> list[int]:
    """How many of ``orbitals`` fall in each subspace."""
    bounds = list(first_orbital) + [norb + 1]
    return [
        sum(1 for p in orbitals if bounds[k] <= p < bounds[k + 1])
        for k in range(len(first_orbital))
    ]


def determinants(
    norb: int,
    nalpha: int,
    nbeta: int,
    first_orbital: list[int],
    min_e: list[int],
    max_e: list[int],
) -> list[tuple[frozenset[int], frozenset[int]]]:
    """Every determinant the windows allow.

    Brute force over all ``C(norb, nalpha) * C(norb, nbeta)`` pairs, testing the
    combined occupation of each subspace directly. Slow and obviously correct,
    which is the trade this file exists to make.

    Orbitals are 1-based, matching the Fortran.
    """
    orbitals = range(1, norb + 1)
    out = []
    for a in itertools.combinations(orbitals, nalpha):
        alpha = frozenset(a)
        occ_a = occupation(alpha, first_orbital, norb)
        for b in itertools.combinations(orbitals, nbeta):
            beta = frozenset(b)
            occ_b = occupation(beta, first_orbital, norb)
            if all(
                lo <= x + y <= hi
                for x, y, lo, hi in zip(occ_a, occ_b, min_e, max_e)
            ):
                out.append((alpha, beta))
    return out


def spin_orbitals(alpha: frozenset[int], beta: frozenset[int], norb: int) -> list[int]:
    """A determinant as one ascending list of spin orbitals.

    Alpha orbital ``p`` is ``p``; beta orbital ``p`` is ``p + norb``. Ascending
    order then means every alpha creation operator stands to the left of every
    beta one, which is the convention that makes alpha and beta phases multiply
    with no extra sign between them.
    """
    return sorted(list(alpha) + [p + norb for p in beta])


def _spatial(i: int, norb: int) -> int:
    return i - 1 if i <= norb else i - norb - 1


def _same_spin(i: int, j: int, norb: int) -> bool:
    return (i <= norb) == (j <= norb)


def _eri_so(eri: np.ndarray, i: int, j: int, k: int, l: int, norb: int) -> float:
    """``(ij|kl)`` over spin orbitals, chemist notation."""
    if not _same_spin(i, j, norb) or not _same_spin(k, l, norb):
        return 0.0
    return eri[_spatial(i, norb), _spatial(j, norb), _spatial(k, norb), _spatial(l, norb)]


def _h_so(h1e: np.ndarray, i: int, j: int, norb: int) -> float:
    if not _same_spin(i, j, norb):
        return 0.0
    return h1e[_spatial(i, norb), _spatial(j, norb)]


def _excite(occ: set[int], remove: int, add: int) -> int:
    """Apply ``a+_add a_remove`` in place, returning the sign.

    The sign is the parity of how many occupied spin orbitals the operator has
    to step over, counted strictly between the two positions -- ``remove`` is an
    endpoint and ``add`` is empty, so neither counts itself.
    """
    lo, hi = (remove, add) if remove < add else (add, remove)
    n = sum(1 for x in occ if lo < x < hi)
    occ.remove(remove)
    occ.add(add)
    return -1 if n % 2 else 1


def matrix_element(
    left: tuple[frozenset[int], frozenset[int]],
    right: tuple[frozenset[int], frozenset[int]],
    norb: int,
    h1e: np.ndarray,
    eri: np.ndarray,
) -> float:
    """``<left|H|right>`` by the Slater-Condon rules."""
    occ_l = set(spin_orbitals(*left, norb))
    occ_r = set(spin_orbitals(*right, norb))

    removed = sorted(occ_l - occ_r)  # occupied on the left, not the right
    added = sorted(occ_r - occ_l)
    if len(removed) != len(added) or len(removed) > 2:
        return 0.0

    common = sorted(occ_l & occ_r)

    if not removed:
        total = sum(_h_so(h1e, i, i, norb) for i in common)
        for i, j in itertools.combinations(common, 2):
            total += _eri_so(eri, i, i, j, j, norb) - _eri_so(eri, i, j, j, i, norb)
        return total

    if len(removed) == 1:
        i, a = removed[0], added[0]
        work = set(occ_l)
        sign = _excite(work, i, a)
        total = _h_so(h1e, i, a, norb)
        for j in common:
            total += _eri_so(eri, i, a, j, j, norb) - _eri_so(eri, i, j, j, a, norb)
        return sign * total

    i, j = removed
    a, b = added
    work = set(occ_l)
    sign = _excite(work, i, a) * _excite(work, j, b)
    return sign * (
        _eri_so(eri, i, a, j, b, norb) - _eri_so(eri, i, b, j, a, norb)
    )


def hamiltonian(
    dets: list[tuple[frozenset[int], frozenset[int]]],
    norb: int,
    h1e: np.ndarray,
    eri: np.ndarray,
) -> np.ndarray:
    """The full Hamiltonian in the determinant basis, dense."""
    n = len(dets)
    matrix = np.zeros((n, n))
    for row in range(n):
        for col in range(row, n):
            value = matrix_element(dets[row], dets[col], norb, h1e, eri)
            matrix[row, col] = value
            matrix[col, row] = value
    return matrix


def solve(
    norb: int,
    nalpha: int,
    nbeta: int,
    first_orbital: list[int],
    min_e: list[int],
    max_e: list[int],
    n_roots: int = 3,
    dominant: bool = False,
) -> tuple[int, np.ndarray, float]:
    """Determinant count, lowest eigenvalues, and the trace of H.

    The trace is worth returning on its own. It is the sum of the diagonal
    elements, so it tests the space and the diagonal without needing an
    eigensolver on either side -- and being a sum, it does not care in what
    order the determinants were enumerated. That makes it the one number a
    diagonal-only implementation can already be held to.
    """
    h1e, eri = model_integrals(norb, dominant)
    dets = determinants(norb, nalpha, nbeta, first_orbital, min_e, max_e)
    if not dets:
        return 0, np.array([]), 0.0
    matrix = hamiltonian(dets, norb, h1e, eri)
    values = np.linalg.eigvalsh(matrix)
    return len(dets), values[:n_roots], float(np.trace(matrix))


def check_against_pyscf() -> None:
    """With windows that exclude nothing, this must be PySCF's full CI.

    The check that matters most. Everything else here is self-consistent by
    construction; only an outside answer can catch a sign convention that is
    wrong in the same way everywhere.
    """
    from pyscf import fci

    print("self-check: unrestricted windows against pyscf.fci")
    for norb, nalpha, nbeta in [(4, 2, 2), (5, 3, 2), (6, 3, 3), (6, 2, 1)]:
        h1e, eri = model_integrals(norb)
        nelec = nalpha + nbeta
        reference, _ = fci.direct_spin1.kernel(h1e, eri, norb, (nalpha, nbeta))

        # One subspace with the windows wide open is the full space.
        count, values, _ = solve(norb, nalpha, nbeta, [1], [nelec], [nelec], n_roots=1)
        expected = math.comb(norb, nalpha) * math.comb(norb, nbeta)

        delta = abs(values[0] - reference)
        status = "ok" if delta < 1e-11 and count == expected else "MISMATCH"
        print(
            f"  ({norb}o, {nalpha}a, {nbeta}b): {count:6d} dets "
            f"(expected {expected:6d})  E = {values[0]:.12f}  "
            f"pyscf {reference:.12f}  d = {delta:.2e}  {status}"
        )
        if status != "ok":
            raise SystemExit("the reference disagrees with pyscf; fix it before using it")


def check_partition_spellings() -> None:
    """Two subspaces with open windows must equal one subspace."""
    print()
    print("self-check: an unrestricted partition is spelling-independent")
    one = solve(4, 2, 2, [1], [4], [4], n_roots=3)
    two = solve(4, 2, 2, [1, 3], [0, 0], [4, 4], n_roots=3)
    delta = float(np.max(np.abs(one[1] - two[1])))
    print(f"  one subspace : {one[0]:4d} dets  E0 = {one[1][0]:.12f}")
    print(f"  two subspaces: {two[0]:4d} dets  E0 = {two[1][0]:.12f}  d = {delta:.2e}")
    if one[0] != two[0] or delta > 1e-12:
        raise SystemExit("the two spellings disagree")


PARTITIONS = [
    ("two subspaces, min 2 in the first", 4, 2, 2, [1, 3], [2, 0], [4, 2]),
    ("three subspaces", 6, 3, 3, [1, 3, 5], [2, 0, 0], [4, 4, 2]),
    ("unequal widths, truncated", 7, 3, 3, [1, 4], [4, 0], [6, 2]),
    ("unequal spins", 6, 3, 1, [1, 3, 5], [1, 0, 0], [4, 4, 2]),
    ("windows that get tightened", 6, 2, 2, [1, 4], [0, 3], [6, 6]),
]


def reference_table() -> None:
    print()
    print("reference energies for the partitions the Fortran tests use")
    for label, norb, na, nb, first, lo, hi in PARTITIONS:
        count, values, trace = solve(norb, na, nb, first, lo, hi, n_roots=3)
        roots = "  ".join(f"{v:.12f}" for v in values)
        print(f"  {label}")
        print(f"    norb={norb} na={na} nb={nb} start={first} min={lo} max={hi}")
        print(f"    {count:6d} determinants   trace(H) = {trace:.12f}")
        print(f"    E = {roots}")


def dominant_table() -> None:
    print()
    print("reference energies with a separated one-electron diagonal")
    print("(the regime a preconditioned solver is meant for)")
    for label, norb, na, nb, first, lo, hi in PARTITIONS:
        count, values, _ = solve(norb, na, nb, first, lo, hi, n_roots=3, dominant=True)
        roots = "  ".join(f"{v:.12f}" for v in values)
        gaps = "  ".join(f"{values[i + 1] - values[i]:.2e}" for i in range(len(values) - 1))
        print(f"  {label}")
        print(f"    {count:6d} determinants   E = {roots}")
        print(f"    gaps {gaps}")


if __name__ == "__main__":
    check_against_pyscf()
    check_partition_spellings()
    reference_table()
    dominant_table()
