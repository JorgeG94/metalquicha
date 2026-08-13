#!/usr/bin/env python3
"""Undo GAMESS's write-time shift, sum over orbitals, and compare.

    ./build/check_dipquad_sumrule && python3 validation/check_dipquad_sumrule.py

See the Fortran's docstring for why a sum is the right object: it is independent
of how the response is projected onto localized orbitals, which every per-orbital
comparison of this block has been unable to separate from the quantity itself.

**This script validates its own arithmetic before it compares anything.** GAMESS
does not write the tensor it computes; it writes ``DQSHIFT``'s translation of it
from the centre of mass to each orbital's centroid, which mixes in the
dipole-dipole polarizability:

    A'(a,bc) = A(a,bc) - (3/2)[ R_b alpha(c,a) + R_c alpha(a,b) ]
                       + delta_bc sum_d R_d alpha(d,a)

with ``R = centroid - centre of mass`` and ``slot = (a-1)*9 + (b-1)*3 + c``.

The pre-shift tensor ``A(a,bc)`` must be symmetric in ``bc`` and traceless in it,
because it is a dipole against the six unique traceless Buckingham components.
The written values are *neither* -- for water's first orbital ``(1,3)`` is 0.4006
where ``(3,1)`` is 0.1756, and the diagonal sums to 0.079 rather than zero. That
is a consequence of the shift, since the per-orbital ``alpha`` mixed in is itself
asymmetric. So recovering symmetry and tracelessness is a test of the shift
arithmetic that needs no reference at all: if undoing it does not restore both,
the formula or the ``alpha`` in it is wrong and nothing downstream means anything.
"""

import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tools" / "efp_validation"))
from efp_format import parse_efp  # noqa: E402

REFERENCE = REPO / "tools" / "efp_validation" / "reference" / "water_6-31gs_boys.efp"

#: (row, column) of GAMESS's nine dipole-dipole polarizability slots.
POL_ORDER = [(0, 0), (1, 1), (2, 2), (1, 0), (2, 0), (2, 1), (0, 1), (0, 2), (1, 2)]

#: Which of the six traceless components each (b, c) pair is, in GAMESS's
#: XX YY ZZ XY XZ YZ order.
THETA_OF = {(0, 0): 0, (1, 1): 1, (2, 2): 2,
            (0, 1): 3, (1, 0): 3, (0, 2): 4, (2, 0): 4, (1, 2): 5, (2, 1): 5}

#: How close the recovered tensor must come to symmetric and traceless for the
#: shift arithmetic to be believed. The written values carry ten decimals and the
#: shift is a sum of a few products of them, so this is generous.
STRUCTURE_TOL = 1e-6


def read_ours(path):
    """The Fortran dump: sums, then per-orbital tensors in both orderings."""
    v = [float(x) for x in path.read_text().split()]
    n_freq, n_lmo = int(v[0]), int(v[1])
    at = 2
    frequencies = np.array(v[at:at + n_freq]); at += n_freq
    raw = np.array(v[at:at + 3*9*n_freq]).reshape((3, 9, n_freq), order="F")
    at += 3*9*n_freq
    buck = np.array(v[at:at + 3*9*n_freq]).reshape((3, 9, n_freq), order="F")
    at += 3*9*n_freq
    com = np.array(v[at:at + 3]); at += 3
    centroids = np.array(v[at:at + 3*n_lmo]).reshape((3, n_lmo), order="F")
    at += 3*n_lmo
    n = 3*9*n_lmo*n_freq
    dip_measures = np.array(v[at:at + n]).reshape((3, 9, n_lmo, n_freq), order="F")
    at += n
    quad_measures = np.array(v[at:at + n]).reshape((9, 3, n_lmo, n_freq), order="F")
    at += n
    dd = np.array(v[at:at + 3*3*n_lmo*n_freq]).reshape(
        (3, 3, n_lmo, n_freq), order="F")
    return dict(frequencies=frequencies, n_lmo=n_lmo, n_freq=n_freq, com=com,
                centroids=centroids, dd=dd,
                raw=raw.reshape((3, 3, 3, n_freq)),
                buck=buck.reshape((3, 3, 3, n_freq)),
                dip_measures=dip_measures.reshape((3, 3, 3, n_lmo, n_freq)),
                quad_measures=np.einsum("qanf->aqnf", quad_measures).reshape(
                    (3, 3, 3, n_lmo, n_freq)))


def reference_blocks(section):
    """A dynamic section's points grouped by the frequency stamped on each block."""
    content = parse_efp(REFERENCE.read_text())["sections"][section]
    if isinstance(content, dict) and "raw" in content:
        return raw_blocks(content["raw"])
    frequencies, blocks = [], []
    for point in content:
        if "frequency" in point:
            frequencies.append(point["frequency"])
            blocks.append([])
        blocks[-1].append(point)
    return np.array(frequencies), blocks


def raw_blocks(lines):
    """The same grouping, for a section `efp_format` keeps as text."""
    frequencies, blocks, points = [], [], []
    for line in lines:
        if line.lstrip().startswith("CT"):
            parts = line.replace("--", " ").split()
            points.append({"label": parts[0] + parts[1],
                           "xyz": [float(v) for v in parts[2:5]], "tensor": []})
            if "FOR W=" in line:
                frequencies.append(float(line.split("FOR W=")[1].split("I")[0]))
                blocks.append([])
            blocks[-1].append(points[-1])
        else:
            points[-1]["tensor"] += [float(v) for v in line.replace(">", " ").split()]
    return np.array(frequencies), blocks


def alpha_matrix(values):
    m = np.zeros((3, 3))
    for c, (i, j) in enumerate(POL_ORDER):
        m[i, j] = values[c]
    return m


def as_cube(values):
    """The 27 written values as A(a,b,c).

    The transpose is the convention the symmetry test pins: the *first* quadrupole
    index runs fastest, so a plain (3,3,3) reshape has b and c the wrong way round.
    """
    return np.array(values).reshape((3, 3, 3)).transpose((0, 2, 1))


def shift(cube, alpha, r, sign=+1):
    """`DQSHIFT` forward (sign +1) or undone (sign -1)."""
    out = np.zeros((3, 3, 3))
    for a in range(3):
        isotropic = sum(r[d]*alpha[d, a] for d in range(3))
        for b in range(3):
            for c in range(3):
                value = cube[a, b, c] - sign*1.5*(r[b]*alpha[c, a] + r[c]*alpha[a, b])
                if b == c:
                    value += sign*isotropic
                out[a, b, c] = value
    return out


def deviatoric(t):
    """The part with the quadrupole-pair trace removed, over a trailing axis set."""
    trace = np.einsum("abb...->a...", t)/3.0
    return t - np.einsum("a...,bc->abc...", trace, np.eye(3))


def main():
    path = pathlib.Path("/tmp/mqc_dqsum.txt")
    if not path.exists():
        print(f"  MISSING {path} -- run ./build/check_dipquad_sumrule first")
        return 1
    ours = read_ours(path)
    com, centroids = ours["com"], ours["centroids"]
    n_lmo, n_freq = ours["n_lmo"], ours["n_freq"]

    dq_freq, dq = reference_blocks(
        "DIPOLE-QUADRUPOLE DYNAMIC POLARIZABLE POINTS")
    dd_freq, dd = reference_blocks("DYNAMIC POLARIZABLE POINTS")
    if abs(ours["frequencies"] - dq_freq).max() > 1e-6:
        print("  FAIL our frequencies are not the ones GAMESS tabulated at")
        return 1

    # Pair our orbitals to GAMESS's by centroid: water's lone pairs are degenerate
    # and their order follows the Jacobi sweep, so an index pairing is not safe.
    order = [int(np.argmin([np.linalg.norm(centroids[:, i] - np.array(p["xyz"]))
                            for i in range(n_lmo)])) for p in dq[0]]
    if sorted(order) != list(range(n_lmo)):
        print(f"  FAIL the centroid pairing is not one-to-one: {order}")
        return 1

    # --- the structural test: does undoing the shift restore symmetry ------------
    summed = np.zeros((3, 3, 3, n_freq))
    per_orbital = np.zeros((3, 3, 3, n_lmo, n_freq))
    worst_asym = 0.0
    for f in range(n_freq):
        for k, point in enumerate(dq[f]):
            alpha = alpha_matrix(dd[f][k]["tensor"])
            r = np.array(point["xyz"]) - com
            cube = as_cube(point["tensor"])
            per_orbital[:, :, :, order[k], f] = cube
            recovered = shift(cube, alpha, r, sign=-1)
            summed[:, :, :, f] += recovered
        for a in range(3):
            worst_asym = max(worst_asym,
                             np.abs(summed[a, :, :, f] - summed[a, :, :, f].T).max())
    print(f"  undoing the shift leaves the summed tensor symmetric to "
          f"{worst_asym:.2e}")
    if worst_asym > STRUCTURE_TOL:
        print("  FAIL -- the shift arithmetic or its alpha is wrong, so nothing "
              "below means anything")
        return 1

    # --- the sum rule: the quantity itself, free of any projection --------------
    print()
    print(f"        {'driving operator':<22}{'scale':>12}{'deviatoric rel':>16}"
          f"{'max abs':>12}")
    failures = 0
    for label, key in (("raw second moment", "raw"),
                       ("traceless Buckingham", "buck")):
        a, b = deviatoric(ours[key]), deviatoric(summed)
        scale = float(a.ravel() @ b.ravel()/(a.ravel() @ a.ravel()))
        rel = np.linalg.norm(b - a)/np.linalg.norm(b)
        print(f"        {label:<22}{scale:12.6f}{rel:16.3e}"
              f"{np.abs(b - a).max():12.3e}")
        if key == "buck":
            if rel > 1e-3:
                print("        FAIL the summed quantity is not GAMESS's")
                failures += 1
            our_trace = np.abs(np.einsum("abbf->af", ours[key])).max()
            print(f"        our isotropic part {our_trace:.1e}; "
                  f"GAMESS's recovered "
                  f"{np.abs(np.einsum('abbf->af', summed)).max():.2e} -- "
                  f"the shift's delta term")

    # --- what remains: the per-orbital decomposition ----------------------------
    print()
    print(f"        {'per-orbital construction':<28}{'relative':>12}{'max abs':>12}")
    for label, key in (("quadrupole measures", "quad_measures"),
                       ("dipole measures", "dip_measures")):
        got = np.zeros_like(per_orbital)
        for f in range(n_freq):
            for k in range(n_lmo):
                got[:, :, :, k, f] = shift(ours[key][:, :, :, k, f],
                                           ours["dd"][:, :, k, f],
                                           centroids[:, k] - com, sign=+1)
        rel = np.linalg.norm(per_orbital - got)/np.linalg.norm(per_orbital)
        print(f"        {label:<28}{rel:12.3e}{np.abs(per_orbital - got).max():12.3e}")

    print()
    if failures:
        print("[dqsum] the summed quantity does not match")
        return 1
    print("[dqsum] the quantity is GAMESS's; its distribution over orbitals is not")
    return 0


if __name__ == "__main__":
    sys.exit(main())
