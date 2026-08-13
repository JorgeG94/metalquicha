#!/usr/bin/env python3
"""Try MAKEFP on every closed-shell geometry in the validation set.

    python3 tools/efp_validation/makefp_sweep.py

Everything about the emitter has been checked on water in 6-31G*: the multipoles,
the polarizabilities, the projection data, the surface a hexamer relaxes on. One
molecule with two elements and two identical bonds is thin ground for a claim about
fragment potentials in general, and the parts most likely to be wrong are exactly
the ones water cannot exercise -- a second row atom, a lone atom with no bonds at
all, a fragment carrying a charge, a bond between two heavy atoms.

So this runs the whole chain on each geometry and reports what came out. It is a
coverage sweep, not an accuracy check: there is no GAMESS reference for most of
these, and what it establishes is that the machinery does not fall over. Where it
does fall over, the failure is the useful part.

Open-shell cases are skipped -- `o2` and `oh` are radicals and EFP2 is defined over a
closed-shell reference.
"""

import argparse
import json
import pathlib
import re
import subprocess
import sys
import tempfile

REPO = pathlib.Path(__file__).resolve().parent.parent.parent
GEOM = REPO / "validation" / "inputs" / "sample_inputs"
MQC = REPO / "build" / "mqc"

#: geometry stem -> (charge, what it is there to cover)
CASES = [
    ("h2", 0, "two electrons, one bond, no core"),
    ("lih", 0, "an alkali, and a very polar bond"),
    ("beh2", 0, "linear, two bonds"),
    ("bh3", 0, "trigonal planar"),
    ("ch4", 0, "the standard organic centre"),
    ("nh3", 0, "a lone pair that is not water's"),
    ("hf", 0, "the most polar diatomic here"),
    ("ne", 0, "a lone atom: no bonds, so no midpoints"),
    ("nah", 0, "second row alkali"),
    ("mgh2", 0, "second row alkaline earth"),
    ("alh3", 0, "second row, three bonds"),
    ("sih4", 0, "second row, four bonds"),
    ("ph3", 0, "second row lone pair"),
    ("h2s", 0, "second row, water's shape"),
    ("hcl", 0, "second row halide"),
    ("ar", 0, "a lone atom with a filled second shell"),
    ("h3op", 1, "a cation, and a charged fragment at all"),
]


def deck(xyz, charge, basis):
    return {
        "schema": {"name": "makefp-sweep", "version": "1.0"},
        "molecules": [{"xyz": str(xyz), "molecular_charge": charge,
                       "molecular_multiplicity": 1}],
        "model": {"method": "hf", "basis": basis},
        "driver": "MakeFP",
    }


def sections(text):
    """Section headers in an emitted potential, so a partial file is visible."""
    out = []
    for line in text.split("\n"):
        s = line.strip()
        if s and s == s.upper() and re.match(r"^[A-Z][A-Z0-9 \-()=.]*$", s) \
                and not s.startswith("STOP") and "$" not in s:
            out.append(s.split("(")[0].strip())
    return out


def run(stem, charge, basis, keep):
    xyz = GEOM / f"{stem}.xyz"
    if not xyz.exists():
        return dict(ok=False, note=f"no geometry {xyz.name}")
    work = pathlib.Path(tempfile.mkdtemp(prefix=f"makefp_{stem}_"))
    path = work / f"{stem}.json"
    path.write_text(json.dumps(deck(xyz, charge, basis), indent=2))
    try:
        done = subprocess.run([str(MQC), path.name], cwd=work, capture_output=True,
                              text=True, timeout=1800)
        log = done.stdout + done.stderr
    except subprocess.TimeoutExpired:
        return dict(ok=False, note="timed out")
    efp = work / f"{stem}.efp"
    if not efp.exists():
        fail = re.search(r"MAKEFP failed: (.*)", log)
        return dict(ok=False, note=fail.group(1).strip() if fail else "no potential written",
                    log=log if keep else None, work=work if keep else None)
    text = efp.read_text()
    got = sections(text)
    npts = len([l for l in text.split("COORDINATES (BOHR)")[1].split("STOP")[0]
                .strip().split("\n") if l.strip()])
    energy = re.search(r"RHF energy\s*(-?\d+\.\d+)", log)
    return dict(ok=True, sections=len(got), points=npts, bytes=len(text),
                energy=float(energy.group(1)) if energy else None,
                work=work if keep else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--basis", default="6-31g*",
                    help="Cartesian s/p/d only -- a spherical basis is refused")
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--only", default=None, help="one stem, for debugging")
    args = ap.parse_args()

    if not MQC.exists():
        print(f"  MISSING {MQC} -- build first")
        return 1
    cases = [c for c in CASES if args.only is None or c[0] == args.only]

    print(f"  MAKEFP over {len(cases)} closed-shell geometries, {args.basis}")
    print()
    print(f"  {'molecule':<10}{'covers':<38}{'result':<34}")
    good = 0
    failures = []
    for stem, charge, why in cases:
        r = run(stem, charge, args.basis, args.keep)
        if r["ok"]:
            good += 1
            note = (f"{r['sections']} sections, {r['points']} points, "
                    f"{r['bytes']//1024} kB")
        else:
            note = "FAILED: " + r["note"][:60]
            failures.append((stem, r["note"]))
        print(f"  {stem:<10}{why:<38}{note:<34}")
        sys.stdout.flush()

    print()
    print(f"  {good} of {len(cases)} wrote a potential")
    if failures:
        print()
        print("  failures in full:")
        for stem, note in failures:
            print(f"    {stem}: {note}")
    return 0 if good == len(cases) else 1


if __name__ == "__main__":
    sys.exit(main())
