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


def load_ours():
    import check_dma
    import check_distributed_polarizability as cdp

    dma = check_dma.read_dump(pathlib.Path("/tmp/mqc_dma_1.txt"))
    pol = cdp.read_dump(pathlib.Path("/tmp/mqc_distributed_1.txt"))
    return dma, pol


def graft(doc, dma, pol):
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
        dma, pol = load_ours()
        touched = graft(doc, dma, pol)

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
