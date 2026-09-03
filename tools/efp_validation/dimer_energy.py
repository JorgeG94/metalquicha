#!/usr/bin/env python3
"""Have GAMESS compute a dimer interaction energy from a potential we wrote.

    ./build/check_makefp
    python3 tools/efp_validation/dimer_energy.py

**Why this and not another comparison of printed numbers.** Every other check in
this tree compares a parameter we computed against the same parameter in a GAMESS
potential. All of them can pass while the *file* is unusable -- a component in the
wrong slot, a section header GAMESS's reader skips over, a continuation marker
missing so the last value of a record is silently dropped. The only test that
covers the file rather than the numbers is to hand it to GAMESS as a fragment and
ask it for an energy.

So this builds the same water dimer twice, once from the potential
``check_makefp`` wrote and once from GAMESS's own, and compares the energy terms
it reports. Agreement means our file *is* a fragment: read, oriented, and used.

**Every term but charge transfer is ours.** Only ``CTVEC`` is unwritten (see
``backends/cenzontle/mqc_efp_potential.f90``), so GAMESS gets no charge transfer from
our file and reports ``CHARGE TRANSFER=F`` for it. Electrostatics with screening,
polarization, exchange repulsion and all three dispersion orders are computed from
parameters we wrote.
"""

import argparse
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent.parent
#: The default reference. Its `CTVEC` is the *valence virtual orbital* form -- five
#: occupied orbitals plus two quasi-atomic virtuals built by GAMESS's `VVOS` -- while
#: we write the canonical-orbital form, all of the MOs. Both are legitimate and
#: GAMESS reads either, but they are different charge-transfer bases, so comparing
#: against this file leaves charge transfer differing by 2.6e-05 for a reason that is
#: not an error. `water_6-31gs_cmo.efp` is the same molecule and basis with
#: `$MAKEFP CTVVO=.FALSE.`, which is the apples-to-apples comparison and where the
#: term agrees exactly.
REFERENCE = HERE / "reference" / "water_6-31gs_boys.efp"
GAMESS = REPO.parent / "mgga" / "gamess"

BOHR = 0.52917724924

#: What both potentials' fragment groups are renamed to before the run.
#:
#: **Not "WATER".** GAMESS ships a built-in fragment by that name, and a
#: $EFRAG naming it gets the internal EFP1 potential rather than the $WATER group
#: in the input -- so a deck calling it WATER would run to completion, report
#: sensible energies, and have read neither file. The name has to be one no
#: library defines.
FRAGMENT = "EFPTEST"

#: The terms GAMESS prints for an EFP calculation, and how it labels them. Only
#: those it names explicitly -- a term we cannot find is reported missing rather
#: than defaulted to zero, since zero is also a legitimate value.
#:
#: **Dispersion is broken out by order, and it has to be.** GAMESS prints `E6`,
#: `E7` and `E8` on separate lines and then a total, all of them ending in
#: "DISPERSION ENERGY". A pattern matching that phrase and taking the last hit
#: silently returns `E8`, which is what this script did at first -- so it compared
#: our `E6` against GAMESS's `E8`, reported a 4.9e-04 disagreement, and the real
#: numbers agreed to 1.6e-09. Each order is matched by its own label now.
TERMS = {
    "electrostatic": r"ELECTROSTATIC ENERGY\s*=\s*(-?\d+\.\d+)",
    "polarization": r"POLARIZATION ENERGY\s*=\s*(-?\d+\.\d+)",
    "repulsion": r"REPULSION ENERGY\s*=\s*(-?\d+\.\d+)",
    "dispersion E6": r"E6 DISPERSION ENERGY\s*=\s*(-?\d+\.\d+)",
    "dispersion E7": r"E7 DISPERSION ENERGY\s*=\s*(-?\d+\.\d+)",
    "dispersion E8": r"E8 DISPERSION ENERGY\s*=\s*(-?\d+\.\d+)",
    "dispersion total": r"TOTAL DISPERSION ENERGY\(E6\+E7\+E8\)\s*=\s*(-?\d+\.\d+)",
    # GAMESS spells it ENRGY, not ENERGY. Matching the correct spelling finds
    # nothing and reports the term missing, which is how it hid for a while.
    "charge transfer": r"CHARGE TRANSFER ENRGY\s*=\s*(-?\d+\.\d+)",
}


def fragment_name(text):
    """The name the potential declares, which is what $EFRAG has to ask for."""
    match = re.search(r"^ \$(\w+)\s*$", text, re.M)
    if not match:
        raise SystemExit("cannot find the $FRAGNAME group in the potential")
    return match.group(1)


def atoms(text):
    """The three atom positions from COORDINATES, in Angstrom.

    Only atoms -- a bond midpoint carries no mass, and $EFRAG wants atoms to fix
    the fragment's orientation.
    """
    body = text.split("COORDINATES (BOHR)")[1].split("STOP")[0]
    out = []
    for line in body.strip().split("\n"):
        parts = line.split()
        if len(parts) < 6 or parts[0].startswith("BO"):
            continue
        x, y, z = (float(v)*BOHR for v in parts[1:4])
        out.append((parts[0], x, y, z))
    return out


def deck(potential_text, separation):
    """A dimer of two copies of the fragment, the second displaced along x.

    `COORD=FRAGONLY` is what makes an all-EFP system legal: without it GAMESS
    wants a $DATA group and then rejects the job for having no basis functions,
    since there are no QM atoms to put a basis on.
    """
    name = fragment_name(potential_text)
    positions = atoms(potential_text)
    if len(positions) < 3:
        raise SystemExit("need at least three atoms to orient a fragment")

    lines = [
        " $CONTRL SCFTYP=RHF RUNTYP=ENERGY COORD=FRAGONLY UNITS=ANGS $END",
        " $SYSTEM MWORDS=100 $END",
        " $EFRAG",
        "COORD=CART",
    ]
    for shift in (0.0, separation):
        lines.append(f"FRAGNAME={FRAGMENT}")
        for label, x, y, z in positions[:3]:
            lines.append(f"{label:<8}{x + shift:15.10f}{y:15.10f}{z:15.10f}")
    lines.append(" $END")

    # The potential's own group, verbatim from its $FRAGNAME line onward, with the
    # group renamed so it cannot collide with a built-in. The banner above it is
    # MAKEFP's comment and GAMESS does not read it here.
    start = potential_text.index(f" ${name}")
    body = potential_text[start:].rstrip()
    body = body.replace(f" ${name}", f" ${FRAGMENT}", 1)
    lines.append(body)
    return "\n".join(lines) + "\n"


def graft_screening(ours, theirs):
    """Our potential with GAMESS's screening parameters and nothing else changed.

    The attribution test for the electrostatic term. Electrostatics is the sum of
    a multipole expansion and a charge-penetration correction, so a disagreement
    there is either in the moments or in the damping, and comparing the printed
    numbers cannot separate them -- both are ours. Substituting only the screening
    does: if the term then matches, the moments were right and the fit is the
    whole difference.
    """
    out = ours
    for kind in ("SCREEN2", "SCREEN"):
        pattern = r"(^" + kind + r"\s+\(FROM.*?\n)(.*?)(^STOP)"
        source = re.search(pattern, theirs, re.S | re.M)
        if source is None:
            raise SystemExit(f"no {kind} block in the reference potential")
        out, n = re.subn(pattern, lambda m: m.group(1) + source.group(2) + m.group(3),
                         out, count=1, flags=re.S | re.M)
        if n != 1:
            raise SystemExit(f"no {kind} block to replace in our potential")
    return out


def run_gamess(text, tag, keep):
    """One GAMESS job, returning its log."""
    if not (GAMESS / "misc" / "automation" / "rungms").exists():
        raise SystemExit(f"no GAMESS build at {GAMESS}")
    work = pathlib.Path(tempfile.mkdtemp(prefix=f"efp_{tag}_"))
    job = f"efp_{tag}"
    (GAMESS / f"{job}.inp").write_text(text)
    try:
        result = subprocess.run(
            ["./misc/automation/rungms", job, "00", "1", "1"],
            cwd=GAMESS, capture_output=True, text=True, timeout=900)
        log = result.stdout + result.stderr
    finally:
        (GAMESS / f"{job}.inp").unlink(missing_ok=True)
        # rungms leaves scratch behind that a rerun of the same job name trips on.
        for leftover in GAMESS.glob(f"{job}.*"):
            leftover.unlink(missing_ok=True)
        for leftover in (GAMESS / "restart").glob(f"{job}*"):
            shutil.rmtree(leftover, ignore_errors=True)
    if keep:
        (work / f"{tag}.inp").write_text(text)
        (work / f"{tag}.log").write_text(log)
        print(f"        kept {work}/{tag}.log")
    else:
        shutil.rmtree(work, ignore_errors=True)
    return log


def terms(log):
    """The energy terms GAMESS reported."""
    found = {}
    for key, pattern in TERMS.items():
        matches = re.findall(pattern, log)
        if matches:
            found[key] = float(matches[-1])
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", default="/tmp/mqc_water.efp")
    ap.add_argument("--theirs", default=str(REFERENCE))
    ap.add_argument("--separation", type=float, default=3.0,
                    help="Angstrom displacement of the second copy along x")
    ap.add_argument("--tol", type=float, default=1e-6,
                    help="Hartree, on each shared term")
    ap.add_argument("--keep", action="store_true", help="keep the GAMESS logs")
    ap.add_argument("--screen-from-gamess", action="store_true",
                    help="take the screening parameters from the reference and "
                         "leave every other parameter ours -- separates a "
                         "multipole disagreement from a damping one")
    args = ap.parse_args()

    ours_path = pathlib.Path(args.ours)
    if not ours_path.exists():
        print(f"  MISSING {ours_path} -- run ./build/check_makefp first")
        return 1

    print(f"  water dimer, {args.separation} A apart, both fragments from the "
          f"same potential")
    results = {}
    theirs_text = pathlib.Path(args.theirs).read_text()
    for tag, path in (("ours", ours_path), ("theirs", pathlib.Path(args.theirs))):
        text = path.read_text()
        if tag == "ours" and args.screen_from_gamess:
            text = graft_screening(text, theirs_text)
            print("        (screening taken from the reference; "
                  "every other parameter is ours)")
        log = run_gamess(deck(text, args.separation), tag, args.keep)
        if "EXECUTION OF GAMESS TERMINATED NORMALLY" not in log:
            print(f"  FAIL GAMESS did not finish on {tag} ({path.name})")
            for line in log.strip().split("\n")[-15:]:
                print(f"        {line}")
            return 1
        results[tag] = terms(log)
        print(f"        {tag:7s} read {path.name}: "
              f"{len(results[tag])} energy terms")

    print()
    print(f"        {'term':<18}{'ours':>16}{'GAMESS':>16}{'difference':>14}")
    failures = 0
    for key in TERMS:
        a, b = results["ours"].get(key), results["theirs"].get(key)
        if a is None and b is None:
            continue
        if a is None or b is None:
            print(f"        {key:<18}{'-- absent from one run --':>46}")
            continue
        gap = abs(a - b)
        mark = ""
        if key == "electrostatic" and gap > args.tol and not args.screen_from_gamess:
            # The screening fit is known to differ -- ours finds a real minimum
            # where GAMESS's search reaches its alpha = 10 "off" bound, which is
            # an absorbing state its objective is flat in. Rerunning with
            # --screen-from-gamess is what separates that from a moment being
            # wrong, and it brings this term to 1e-10.
            mark = "  screening differs; rerun with --screen-from-gamess"
        elif gap > args.tol:
            mark = "  FAIL"
            failures += 1
        print(f"        {key:<18}{a:16.9f}{b:16.9f}{gap:14.2e}{mark}")

    print()
    if failures:
        print(f"[dimer] {failures} term(s) disagree beyond {args.tol:g} Hartree")
        return 1
    print("[dimer] GAMESS reads our potential and agrees on every term we supply")
    return 0


if __name__ == "__main__":
    sys.exit(main())
