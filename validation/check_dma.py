#!/usr/bin/env python3
"""Compare our distributed multipoles against a GAMESS MAKEFP potential.

    ./build/check_dma && python3 validation/check_dma.py

Milestone 4. The sum rule is checked in the Fortran, where it needs no reference
and holds to machine precision. What is left for here is the part it cannot see:
**where the density actually went**, and in which order the components are
written. Move a dipole from one atom to another, or transpose a quadrupole, and
the sum rule is unmoved -- which is exactly how the polarizability transpose in
M5 survived every internal check until a per-point comparison found it.

So this compares point by point, and reports each moment order separately, since
a partition error and a packing error look nothing alike: a wrong partition moves
one point's numbers and leaves its neighbour's, while a wrong packing permutes
components uniformly everywhere.
"""

import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tools" / "efp_validation"))
from efp_format import parse_efp  # noqa: E402

REFERENCE = REPO / "tools" / "efp_validation" / "reference" / "water_6-31gs_boys.efp"

#: Agreement demanded, in atomic units, per moment order.
#:
#: The reference is printed to ten decimals for the monopoles through quadrupoles
#: and nine for the octopoles, and both codes solve the same SCF to 1e-12, so
#: these are set by the file and not by the physics. Loosened only for the
#: octopole, where the printed precision is a digit worse and the magnitudes are
#: larger.
TOL = {"monopole": 1e-8, "dipole": 1e-8, "quadrupole": 1e-8, "octopole": 1e-7}


def read_dump(path):
    lines = path.read_text().split("\n")
    basis, molecule = lines[0].strip(), lines[1].strip()
    n = int(lines[2])
    labels = [lines[3 + i].strip() for i in range(n)]
    values = [float(x) for x in " ".join(lines[3 + n:]).split()]
    at = 0

    def take(count, shape):
        nonlocal at
        block = np.array(values[at:at + count]).reshape(shape, order="F")
        at += count
        return block

    points = take(3*n, (3, n))
    electronic = take(n, (n,))
    nuclear = take(n, (n,))
    dipole = take(3*n, (3, n))
    quadrupole = take(6*n, (6, n))
    octopole = take(10*n, (10, n))
    return dict(basis=basis, molecule=molecule, labels=labels, points=points,
                electronic=electronic, nuclear=nuclear, dipole=dipole,
                quadrupole=quadrupole, octopole=octopole)


def reference_potential():
    s = parse_efp(REFERENCE.read_text())["sections"]
    n = len(s["COORDINATES"])
    return dict(
        labels=[r["label"] for r in s["COORDINATES"]],
        points=np.array([r["values"][:3] for r in s["COORDINATES"]]).T,
        electronic=np.array([r["values"][0] for r in s["MONOPOLES"]]),
        nuclear=np.array([r["values"][1] for r in s["MONOPOLES"]]),
        dipole=np.array([r["values"][:3] for r in s["DIPOLES"]]).T,
        quadrupole=np.array([r["values"][:6] for r in s["QUADRUPOLES"]]).T,
        octopole=np.array([r["values"][:10] for r in s["OCTUPOLES"]]).T,
        n=n)


def main():
    if not REFERENCE.exists():
        print(f"  MISSING {REFERENCE}")
        return 1
    ref = reference_potential()
    failures = 0

    for index in (1, 2):
        path = pathlib.Path(f"/tmp/mqc_dma_{index}.txt")
        if not path.exists():
            print(f"  MISSING {path} -- run ./build/check_dma first")
            failures += 1
            continue
        ours = read_dump(path)

        print(f"  {ours['molecule']}/{ours['basis']}: {len(ours['labels'])} points  "
              f"{' '.join(ours['labels'])}")

        # Only water has a GAMESS reference: the HCN MAKEFP run aborts in its own
        # CPHF before writing the polarizabilities, and although its multipoles
        # are present in the partial file it is not a potential we generated
        # end to end. See tools/efp_validation/reference/README.md.
        if ours["molecule"] != "water":
            print("        (no complete GAMESS reference for this one; sum rule "
                  "checked in the Fortran)")
            continue

        if ours["labels"] != ref["labels"]:
            print(f"        FAIL labels differ: ours {ours['labels']}  "
                  f"GAMESS {ref['labels']}")
            failures += 1
            continue

        gap = np.abs(ours["points"] - ref["points"]).max()
        print(f"        expansion points agree to {gap:.2e} Bohr")
        if gap > 1e-9:
            print("        FAIL the expansion points are not in the same places")
            failures += 1
            continue

        worst = {}
        for name, key in (("monopole", "electronic"), ("dipole", "dipole"),
                          ("quadrupole", "quadrupole"), ("octopole", "octopole")):
            diff = np.abs(ours[key] - ref[key])
            worst[name] = diff.max()
            # Which point is worst, since a partition error is localized.
            where = int(np.unravel_index(diff.argmax(), diff.shape)[-1])
            flag = "ok  " if worst[name] < TOL[name] else "FAIL"
            print(f"        {flag} {name:11s} worst {worst[name]:.2e}  "
                  f"at {ours['labels'][where]}")
            if worst[name] >= TOL[name]:
                failures += 1

        if any(worst[k] >= TOL[k] for k in worst):
            print()
            print("        ours vs GAMESS, per point:")
            for k, label in enumerate(ours["labels"]):
                print(f"          {label:6s} q {ours['electronic'][k]:14.8f} "
                      f"{ref['electronic'][k]:14.8f}")
                print(f"                 mu {np.array2string(ours['dipole'][:, k], precision=6)}")
                print(f"                    {np.array2string(ref['dipole'][:, k], precision=6)}")

    print()
    if failures:
        print(f"[dma] {failures} FAILURE(S)")
        return 1
    print("[dma] our distributed multipoles match GAMESS's MAKEFP")
    return 0


if __name__ == "__main__":
    sys.exit(main())
