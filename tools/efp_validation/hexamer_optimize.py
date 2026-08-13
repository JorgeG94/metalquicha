#!/usr/bin/env python3
"""Optimize a water hexamer in GAMESS from our potential and from its own.

    ./build/check_makefp
    python3 tools/efp_validation/hexamer_optimize.py

**Why an optimization and not another single point.** `dimer_energy.py` compares
energy terms at one fixed geometry, which tests the parameters. It does not test
whether they describe a *surface*: a potential can give the right energy at one
separation and the wrong gradient, and nothing in a single-point comparison
notices. Letting GAMESS relax six fragments -- thirty-six degrees of freedom, a few
dozen steps, every term active at every step -- tests the whole surface, and the
comparison is the structure it settles into.

Both runs use the same monomer potential geometry and the same starting hexamer, so
the only difference is which file the parameters came from. What is compared:

* the final energy;
* the optimized structure, through its fifteen oxygen-oxygen distances, which are
  rotation and translation invariant and so need no alignment;
* the step count, since a potential with a subtly wrong gradient converges
  differently even when it lands somewhere similar.

The starting geometry is the water prism from the validation inputs.

**What it found.** With the screening taken from the reference and every other
parameter ours, the two potentials relax to the same structure -- 4.7e-05 Angstrom
across the fifteen distances, 1.3e-06 Hartree in energy, the same number of steps.
With our own screening fit instead, 1.2e-02 Angstrom and 7.4e-01 kcal/mol. So the
surface is right and the charge-penetration fit is the one thing that differs, which
is what `--screen-from-gamess` separates. That is also the first physical size for
that difference: it is worth three quarters of a kcal/mol on a water hexamer, which
is not negligible for fragment work.
"""

import argparse
import pathlib
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent.parent
GAMESS = REPO.parent / "mgga" / "gamess"
PRISM = REPO / "validation" / "inputs" / "sample_inputs" / "prism.xyz"
REFERENCE = HERE / "reference" / "water_6-31gs_cmo.efp"

#: Renamed so neither run can pick up GAMESS's built-in WATER fragment instead of
#: the group in the input -- see dimer_energy.py, where that trap is documented.
FRAGMENT = "EFPTEST"


def waters(path):
    """The hexamer as a list of (label, xyz) triples, three atoms each."""
    lines = path.read_text().strip().split("\n")
    atoms = []
    for line in lines[2:]:
        parts = line.split()
        if len(parts) == 4:
            atoms.append((parts[0], np.array([float(v) for v in parts[1:]])))
    if len(atoms) % 3:
        raise SystemExit(f"{path.name}: {len(atoms)} atoms is not whole waters")
    return [atoms[i:i + 3] for i in range(0, len(atoms), 3)]


def fragment_labels(potential_text):
    """The atom labels the potential uses, in order -- A01O, A02H, A03H."""
    body = potential_text.split("COORDINATES (BOHR)")[1].split("STOP")[0]
    out = []
    for line in body.strip().split("\n"):
        parts = line.split()
        if len(parts) >= 6 and not parts[0].startswith("BO"):
            out.append(parts[0])
    return out


def deck(potential_text, groups, name):
    """An all-EFP optimization of however many fragments `groups` holds."""
    labels = fragment_labels(potential_text)
    lines = [
        " $CONTRL SCFTYP=RHF RUNTYP=OPTIMIZE COORD=FRAGONLY UNITS=ANGS $END",
        " $SYSTEM MWORDS=100 $END",
        " $STATPT NSTEP=200 OPTTOL=1.0D-05 $END",
        " $EFRAG",
        "COORD=CART",
    ]
    for group in groups:
        lines.append(f"FRAGNAME={FRAGMENT}")
        for label, (_symbol, xyz) in zip(labels, group):
            lines.append(f"{label:<8}" + "".join(f"{v:15.10f}" for v in xyz))
    lines.append(" $END")
    start = potential_text.index(f" ${name}")
    body = potential_text[start:].rstrip().replace(f" ${name}", f" ${FRAGMENT}", 1)
    lines.append(body)
    return "\n".join(lines) + "\n"


def fragment_name(text):
    match = re.search(r"^ \$(\w+)\s*$", text, re.M)
    if not match:
        raise SystemExit("cannot find the $FRAGNAME group in the potential")
    return match.group(1)


def graft_screening(ours, theirs):
    """Our potential with the reference's screening blocks and nothing else changed."""
    out = ours
    for kind in ("SCREEN2", "SCREEN"):
        pattern = r"(^" + kind + r"\s+\(FROM.*?\n)(.*?)(^STOP)"
        src = re.search(pattern, theirs, re.S | re.M)
        if src is None:
            raise SystemExit(f"no {kind} block in the reference")
        out, n = re.subn(pattern, lambda m: m.group(1) + src.group(2) + m.group(3),
                         out, count=1, flags=re.S | re.M)
        if n != 1:
            raise SystemExit(f"no {kind} block to replace in ours")
    return out


def run_gamess(text, tag, keep):
    if not (GAMESS / "misc" / "automation" / "rungms").exists():
        raise SystemExit(f"no GAMESS build at {GAMESS}")
    job = f"efpopt_{tag}"
    (GAMESS / f"{job}.inp").write_text(text)
    try:
        done = subprocess.run(["./misc/automation/rungms", job, "00", "1", "1"],
                              cwd=GAMESS, capture_output=True, text=True,
                              timeout=7200)
        log = done.stdout + done.stderr
    finally:
        (GAMESS / f"{job}.inp").unlink(missing_ok=True)
        for leftover in GAMESS.glob(f"{job}.*"):
            leftover.unlink(missing_ok=True)
        for leftover in (GAMESS / "restart").glob(f"{job}*"):
            shutil.rmtree(leftover, ignore_errors=True)
    if keep:
        out = pathlib.Path(f"/tmp/{job}.log")
        out.write_text(log)
        print(f"        kept {out}")
    return log


def final_state(log):
    """Energy, oxygen positions and step count from an optimization log."""
    energies = [float(x) for x in
                re.findall(r"TOTAL ENERGY\s*=\s*(-?\d+\.\d+)", log)]
    steps = len(re.findall(r"NSERCH:\s*\d+", log)) or len(energies)
    converged = "EQUILIBRIUM GEOMETRY LOCATED" in log
    # The last printed set of fragment atom coordinates is the optimized one. The
    # block runs to the first line that is neither a FRAGNAME marker nor an atom.
    lines = log.split("\n")
    starts = [i for i, l in enumerate(lines)
              if "COORDINATES OF FRAGMENT MULTIPOLE CENTERS" in l]
    oxygens = []
    if starts:
        for line in lines[starts[-1] + 3:]:
            parts = line.split()
            if line.strip().startswith("FRAGNAME"):
                continue
            if len(parts) != 4:
                if oxygens:
                    break
                continue
            try:
                xyz = [float(v) for v in parts[1:]]
            except ValueError:
                if oxygens:
                    break
                continue
            if parts[0].upper().startswith("A") and parts[0].upper().endswith("O"):
                oxygens.append(xyz)
    return dict(energy=energies[-1] if energies else None,
                steps=steps, converged=converged,
                oxygens=np.array(oxygens) if oxygens else None)


def oo_distances(xyz):
    n = len(xyz)
    return np.array(sorted(np.linalg.norm(xyz[i] - xyz[j])
                           for i in range(n) for j in range(i + 1, n)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", default="/tmp/mqc_water.efp")
    ap.add_argument("--theirs", default=str(REFERENCE))
    ap.add_argument("--geometry", default=str(PRISM))
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--screen-from-gamess", action="store_true",
                    help="take the charge-penetration screening from the reference "
                         "and leave every other parameter ours, so an energy "
                         "difference is attributable")
    args = ap.parse_args()

    ours = pathlib.Path(args.ours)
    if not ours.exists():
        print(f"  MISSING {ours} -- run ./build/check_makefp first")
        return 1
    groups = waters(pathlib.Path(args.geometry))
    print(f"  {len(groups)} water fragments from {pathlib.Path(args.geometry).name},"
          f" all-EFP geometry optimization")

    got = {}
    theirs_text = pathlib.Path(args.theirs).read_text()
    for tag, path in (("ours", ours), ("theirs", pathlib.Path(args.theirs))):
        text = path.read_text()
        if tag == "ours" and args.screen_from_gamess:
            text = graft_screening(text, theirs_text)
            print("        (screening from the reference; everything else is ours)")
        log = run_gamess(deck(text, groups, fragment_name(text)), tag, args.keep)
        if "EXECUTION OF GAMESS TERMINATED NORMALLY" not in log:
            print(f"  FAIL GAMESS did not finish on {tag} ({path.name})")
            for line in log.strip().split("\n")[-12:]:
                print(f"        {line}")
            return 1
        got[tag] = final_state(log)
        state = got[tag]
        print(f"        {tag:7s} {path.name:28s} E = {state['energy']:.9f}"
              f"  {state['steps']:3d} steps"
              f"  {'converged' if state['converged'] else 'NOT CONVERGED'}")

    print()
    a, b = got["ours"], got["theirs"]
    gap = abs(a["energy"] - b["energy"])
    print(f"        final energy differs by {gap:.2e} Hartree "
          f"({gap*627.509:.2e} kcal/mol)")
    failures = 0
    if not (a["converged"] and b["converged"]):
        print("        FAIL at least one optimization did not converge")
        failures += 1
    if a["oxygens"] is not None and b["oxygens"] is not None:
        da, db = oo_distances(a["oxygens"]), oo_distances(b["oxygens"])
        if len(da) == len(db):
            worst = np.abs(da - db).max()
            print(f"        {len(da)} oxygen-oxygen distances agree to "
                  f"{worst:.2e} Angstrom")
            print(f"        range {da.min():.3f} to {da.max():.3f} Angstrom")
            if worst > 1.0e-3:
                print("        FAIL the optimized structures differ")
                failures += 1
        else:
            print(f"        FAIL parsed {len(da)} oxygens against {len(db)}")
            failures += 1
    else:
        # Not a note. A pass here would be claiming an agreement that was never
        # measured, which is worse than a failure.
        print("        FAIL could not parse the optimized coordinates, so the "
              "structures were not compared")
        failures += 1

    print()
    if failures:
        print("[hexamer] the two potentials do not give the same surface")
        return 1
    print("[hexamer] our potential and GAMESS's relax to the same structure")
    return 0


if __name__ == "__main__":
    sys.exit(main())
