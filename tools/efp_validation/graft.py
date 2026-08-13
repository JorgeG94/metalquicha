#!/usr/bin/env python3
"""Put our computed parameters into a GAMESS potential, so GAMESS can judge them.

    ./build/check_dma && ./build/check_distributed_polarizability
    python3 tools/efp_validation/graft.py > /tmp/ours.efp

**Why graft rather than emit a whole potential.** Comparing our numbers against
GAMESS's printed ones -- which M4 and M5 already do, to 1e-9 and 2e-5 -- shows we
compute the same parameters. It does not show that GAMESS, reading a file we
wrote, gets the same *energy* from them. Those are different claims, and only the
second one is end to end: it exercises our component ordering, our sign
conventions and our units all the way through a program with no reason to be kind.

The sections we cannot yet produce -- the dynamic polarizabilities, the projection
basis and wavefunction, the screening parameters, charge transfer -- are kept from
GAMESS's own file. So an energy difference localizes to what we substituted, which
is the point. This is a measurement instrument, not a step toward shipping.

Substituted: MONOPOLES, DIPOLES, QUADRUPOLES, OCTUPOLES (M4), POLARIZABLE POINTS
and LMO CENTROIDS (M5). Everything else is GAMESS's.
"""

import argparse
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "validation"))
from efp_format import parse_efp, render_efp  # noqa: E402

REFERENCE = HERE / "reference" / "water_6-31gs_boys.efp"

#: (row, column) of each of the nine polarizability components, in GAMESS's order.
#: The off-diagonal triples are the transpose of what its labels suggest -- see
#: validation/check_distributed_polarizability.py, where that was measured.
POL_ORDER = [(0, 0), (1, 1), (2, 2), (1, 0), (2, 0), (2, 1), (0, 1), (0, 2), (1, 2)]


#: libcint's index for each of GAMESS's Cartesian d slots, and the normalization
#: between the two codes' d functions. Both established in
#: validation/check_projection.py, which checks them against GAMESS's own
#: coefficients.
D_FROM_LIBCINT = [0, 3, 5, 1, 2, 4]
D_NORMALIZATION = 1.585330892


def gamess_primitive_norm(l, exponent):
    """The factor GAMESS folds into a printed contraction coefficient."""
    double_factorial = 1.0
    n = 2*l - 1
    while n > 1:
        double_factorial *= n
        n -= 2
    return ((2.0*exponent/np.pi)**0.75 * (4.0*exponent)**(l/2.0)
            / np.sqrt(double_factorial))


def projection_basis_lines(atoms, points, basis_name):
    """PROJECTION BASIS SET, in GAMESS's own columns and normalization.

    Shell types are named the way GAMESS names them, including "L" for a shared-
    exponent sp pair, because that is what its reader expects -- even though our
    own basis reader hands libcint the s and p separately. The two land on the same
    basis functions either way, which is what
    validation/check_projection.py establishes.
    """
    import json

    # basis_sets/ sanitises "*" to "_st_", as the CMake extraction names them.
    path = REPO / "basis_sets" / (basis_name.replace("*", "_st_") + ".json")
    table = json.loads(path.read_text())
    lines = []
    primitive = 0
    for index, (symbol, z, xyz) in enumerate(atoms, start=1):
        valence = z - (0 if z <= 2 else 2 if z <= 10 else 10)
        # Columns measured off a real MAKEFP file: label in ten, three
        # coordinates in fifteen each, then the valence electron count in seven.
        lines.append(f"{points[index - 1][0]:<10}" +
                     "".join(f"{v:15.10f}" for v in xyz) +
                     f"{float(valence):7.1f}")
        for shell in table["elements"][str(z)]["electron_shells"]:
            angular = shell["angular_momentum"]
            name = "L" if angular == [0, 1] else {0: "S", 1: "P", 2: "D",
                                                  3: "F"}[angular[0]]
            nprim = len(shell["exponents"])
            lines.append(f"   {name}{nprim:11d}")
            for j, exponent in enumerate(shell["exponents"]):
                primitive += 1
                exponent = float(exponent)
                values = ""
                for column, l in enumerate(angular):
                    coefficient = float(shell["coefficients"][column][j])
                    values += f"{coefficient*gamess_primitive_norm(l, exponent):15.8f}"
                lines.append(f"{primitive:6d}{exponent:21.10f}{values}")
        lines.append("  ")
    return lines


def projection_wavefunction_lines(orbitals, shells):
    """PROJECTION WAVEFUNCTION: our LMO coefficients in GAMESS's basis order."""
    nao, n_lmo = orbitals.shape
    mapped = orbitals.copy()
    for atom, l, offset, dim in shells:
        if l == 2 and dim == 6:
            scale = [D_NORMALIZATION]*3 + [D_NORMALIZATION/np.sqrt(3.0)]*3
            for slot, (src, factor) in enumerate(zip(D_FROM_LIBCINT, scale)):
                mapped[offset + slot, :] = orbitals[offset + src, :]*factor
        elif l >= 2:
            raise SystemExit(f"no basis-function map for l={l} yet")
    lines = []
    for i in range(n_lmo):
        chunk = 0
        for start in range(0, nao, 5):
            chunk += 1
            values = "".join(f"{v:15.8E}" for v in mapped[start:start + 5, i])
            lines.append(f"{i + 1:2d}{chunk:3d}{values}")
    return lines


def fock_lines(fock):
    """FOCK MATRIX ELEMENTS: the lower triangle, four to a line."""
    values = [fock[i, j] for i in range(fock.shape[0]) for j in range(i + 1)]
    lines = []
    for start in range(0, len(values), 4):
        chunk = values[start:start + 4]
        marker = " >" if start + 4 < len(values) else ""
        lines.append("".join(f"{v:16.10f}" for v in chunk) + marker)
    return lines


def load_ours():
    import check_dma
    import check_distributed_polarizability as cdp
    import check_projection as cproj

    dma = check_dma.read_dump(pathlib.Path("/tmp/mqc_dma_1.txt"))
    pol = cdp.read_dump(pathlib.Path("/tmp/mqc_distributed_1.txt"))
    proj = read_projection(pathlib.Path("/tmp/mqc_projection_1.txt"))
    return dma, pol, proj


def read_projection(path):
    """The Fock matrix, centroids, LMO coefficients and shell layout."""
    lines = path.read_text().split("\n")
    tokens = " ".join(lines).split()
    basis = tokens[0]
    values = [float(x) for x in tokens[1:]]
    nao, n_occ, n_lmo = int(values[0]), int(values[1]), int(values[2])
    at = 3
    fock = np.array(values[at:at + n_lmo*n_lmo]).reshape((n_lmo, n_lmo), order="F")
    at += n_lmo*n_lmo
    centroids = np.array(values[at:at + 3*n_lmo]).reshape((3, n_lmo), order="F")
    at += 3*n_lmo
    orbitals = np.array(values[at:at + nao*n_lmo]).reshape((nao, n_lmo), order="F")
    at += nao*n_lmo
    nbas = int(values[at])
    at += 1
    shells = [tuple(int(v) for v in values[at + 4*k:at + 4*k + 4]) for k in range(nbas)]
    return dict(basis=basis, nao=nao, n_occ=n_occ, n_lmo=n_lmo, fock=fock,
                centroids=centroids, orbitals=orbitals, shells=shells)


def graft(doc, dma, pol, proj):
    """Replace what we can compute, leave the rest, and say what was touched."""
    s = doc["sections"]
    touched = []

    labels = [r["label"] for r in s["COORDINATES"]]
    if labels != dma["labels"]:
        raise SystemExit(f"expansion points differ: {labels} vs {dma['labels']}")

    for k in range(len(labels)):
        s["MONOPOLES"][k]["values"] = [float(dma["electronic"][k]),
                                       float(dma["nuclear"][k])]
        s["DIPOLES"][k]["values"] = [float(v) for v in dma["dipole"][:, k]]
        s["QUADRUPOLES"][k]["values"] = [float(v) for v in dma["quadrupole"][:, k]]
        s["OCTUPOLES"][k]["values"] = [float(v) for v in dma["octopole"][:, k]]
    touched += ["MONOPOLES", "DIPOLES", "QUADRUPOLES", "OCTUPOLES"]

    # Polarizable points, paired to the reference by nearest centroid the same way
    # the M5 comparison does, so a reordering cannot silently mismatch them.
    ref_xyz = [np.array(p["xyz"]) for p in s["POLARIZABLE POINTS"]]
    used = set()
    for j, target in enumerate(ref_xyz):
        distances = [np.linalg.norm(pol["centroids"][:, i] - target)
                     for i in range(pol["n_lmo"])]
        i = int(np.argmin(distances))
        if i in used:
            raise SystemExit("two of our LMOs claim the same polarizable point")
        used.add(i)
        tensor = [float(pol["tensors"][a, b, i]) for a, b in POL_ORDER]
        s["POLARIZABLE POINTS"][j]["tensor"] = tensor
        s["POLARIZABLE POINTS"][j]["xyz"] = [float(v) for v in pol["centroids"][:, i]]
        s["LMO CENTROIDS"][j]["values"] = [float(v) for v in pol["centroids"][:, i]]
    touched += ["POLARIZABLE POINTS", "LMO CENTROIDS"]

    # The projection data. Raw sections, so they are built as text in GAMESS's own
    # columns; the basis-function map and the two normalizations these need are
    # established in validation/check_projection.py against GAMESS's own values.
    atoms = [("O", 8, dma["points"][:, 0]), ("H", 1, dma["points"][:, 1]),
             ("H", 1, dma["points"][:, 2])]
    labels = [(r["label"],) for r in s["COORDINATES"]]
    s["PROJECTION BASIS SET"] = {
        "raw": projection_basis_lines(atoms, labels, "6-31g*")}
    s["PROJECTION WAVEFUNCTION"] = {
        "raw": projection_wavefunction_lines(proj["orbitals"], proj["shells"])}
    s["FOCK MATRIX ELEMENTS"] = {"raw": fock_lines(proj["fock"])}
    touched += ["PROJECTION BASIS SET", "PROJECTION WAVEFUNCTION",
                "FOCK MATRIX ELEMENTS"]

    return touched


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default=None,
                    help="rename the fragment group (GAMESS defaults it to FRAGNAME)")
    ap.add_argument("--untouched", action="store_true",
                    help="re-render GAMESS's own parameters, substituting nothing; "
                         "the control that separates our numbers from our renderer")
    args = ap.parse_args()

    doc = parse_efp(REFERENCE.read_text())
    if args.untouched:
        touched = []
    else:
        dma, pol, proj = load_ours()
        touched = graft(doc, dma, pol, proj)

    if args.name:
        doc["name"] = args.name
    print(render_efp(doc), end="")
    print(f"# substituted: {', '.join(touched) if touched else 'nothing'}",
          file=sys.stderr)
    print(f"# kept from GAMESS: "
          f"{', '.join(k for k in doc['sections'] if k not in touched)}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
