#!/usr/bin/env python3
"""Time this build on this machine, and say what changed since last time.

Run it once after configuring and it records a baseline for the machine; run it
again after a change and it compares against that baseline and flags anything
outside the measured noise. It is a performance regression suite and a "does
this build perform the way it should here" report, and those want the same
measurements.

Four things it does that a wall-clock timer does not, each of which came from a
real investigation where the plain number was misleading:

  * **Repeats, and it reports the spread.** A single timing is not a
    measurement. Density functional energies here repeat to under one per cent
    and gradients used to vary by sixteen, and a change smaller than the spread
    is not a change. Comparisons are made against the spread, never against a
    fixed percentage.

  * **Hartree-Fock alongside every functional.** A pure functional carries no
    exact exchange and should cost *less* than Hartree-Fock. When this suite was
    written it cost twenty-six times as much, and no absolute number would have
    said so -- only the ratio to a reference that shares the same integrals.

  * **The same case at every thread count.** A stage that does not parallelise
    is invisible in a run at one thread count: it just looks expensive. The
    exchange-correlation quadrature was serial for the whole life of the code
    and this is the report that would have shown it.

  * **A fragmented case, serial and under MPI.** Fragment work pins itself to
    one OpenMP thread and parallelises with MPI instead, so a change that helps
    a single molecule can hurt the fragment path. One did: threaded BLAS is five
    times faster on one molecule and thirty-one per cent slower on four ranks.

The build configuration is recorded with the numbers, because the BLAS choice
alone moves a run by a factor of five and cannot be recovered afterwards.
"""

import argparse
import json
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
GEOM = HERE / "geometries"

STAGE_ROW = re.compile(r"^\s{4}(\S.*?)\s{2,}(\d+\.\d+) s\s+(\d+)\s+(\d+\.\d+) s\s*$")
TOTAL_LINE = re.compile(r"Total processing time:\s*([0-9.]+)")
FEATURES = re.compile(r"^features:\s*(.*)$", re.M)
BUILDINFO = re.compile(r"^build:\s*(.*)$", re.M)

# One per rung, plus the reference. Hartree-Fock is not decoration: it shares
# the Coulomb and exchange build with the hybrids and has no quadrature at all,
# so it separates "the functional is expensive" from "this molecule is big".
# (name, method, functional, geometry, basis, family)
#
# Geometry and basis are per case rather than global because the methods do not
# share a cost scale. Second-order perturbation theory goes as the fifth power
# of the basis and coupled cluster as the sixth or seventh, so a system where
# MP2 is worth timing makes CCSD(T) take an hour, and one where CCSD(T) is
# quick leaves MP2 too small to measure. Each is placed where it lands in a few
# seconds and where its own work, not the reference SCF, is what moves.
ENERGY_CASES = [
    ("hf", "hf", None, "w10", "6-31g", "reference, no XC"),
    ("svwn", "dft", "svwn", "w10", "6-31g", "LDA"),
    ("pbe", "dft", "pbe", "w10", "6-31g", "GGA"),
    ("tpss", "dft", "tpss", "w10", "6-31g", "meta-GGA"),
    ("b3lyp", "dft", "b3lyp", "w10", "6-31g", "global hybrid"),
    ("camb3lyp", "dft", "camb3lyp", "w10", "6-31g", "range-separated hybrid"),
]
CORRELATED_CASES = [
    ("mp2", "mp2", None, "w10", "cc-pvdz", "conventional MP2"),
    ("ri-mp2", "ri-mp2", None, "w10", "cc-pvdz", "density-fitted MP2"),
    ("ccsd", "ccsd", None, "w5", "6-31g", "coupled cluster, singles and doubles"),
    ("ccsd(t)", "ccsd(t)", None, "w5", "6-31g", "and perturbative triples"),
]
GRADIENT_CASES = [
    ("hf", "hf", None, "w10", "6-31g", "reference, no XC"),
    ("pbe", "dft", "pbe", "w10", "6-31g", "GGA"),
    ("tpss", "dft", "tpss", "w10", "6-31g", "meta-GGA"),
]


def deck(path, geometry, method, functional, driver, basis, fragments=None, level=None):
    model = {"method": method, "basis": basis}
    if functional:
        model["functional"] = functional
    molecule = {"xyz": str(geometry), "molecular_charge": 0, "molecular_multiplicity": 1}
    keywords = {"scf": {"maxiter": 100, "tolerance": 1e-8}, "dft": {"grid_level": 3}}
    if fragments:
        molecule["fragments"] = fragments
        keywords["fragmentation"] = {"method": "mbe", "level": level}
    path.write_text(json.dumps({
        "schema": {"name": "mqc-benchmark", "version": "1.0"},
        "molecules": [molecule],
        "model": model,
        "keywords": keywords,
        "system": {"logger": {"level": "Verbose"}},
        "driver": driver,
    }, indent=2) + "\n")
    return path


def run_once(exe, deck_path, threads, ranks=1):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    cmd = [str(exe), str(deck_path)]
    if ranks > 1:
        cmd = ["mpirun", "-np", str(ranks)] + cmd
    start = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=deck_path.parent)
    wall = time.perf_counter() - start
    if proc.returncode != 0:
        return {"failed": (proc.stdout + proc.stderr).strip().splitlines()[-3:]}
    stages = {}
    for line in proc.stdout.splitlines():
        hit = STAGE_ROW.match(line)
        if hit and hit.group(1).strip() not in ("stage", "total"):
            stages[hit.group(1).strip()] = float(hit.group(2))
    reported = TOTAL_LINE.search(proc.stdout)
    return {"wall": float(reported.group(1)) if reported else wall, "stages": stages}


def repeat(exe, deck_path, threads, repeats, ranks=1):
    runs = [run_once(exe, deck_path, threads, ranks) for _ in range(repeats)]
    bad = [r for r in runs if "failed" in r]
    if bad:
        return {"failed": bad[0]["failed"]}
    walls = [r["wall"] for r in runs]
    median = statistics.median(walls)
    stages = {k: statistics.median(r["stages"].get(k, 0.0) for r in runs)
              for k in runs[0]["stages"]}
    # Spread is what a later comparison is judged against, so it is carried
    # with the number rather than recomputed from a single run.
    spread = (max(walls) - min(walls)) / median * 100 if median else 0.0
    return {"median": median, "spread_pct": spread, "n_runs": len(walls),
            "runs": walls, "stages": stages}


def describe_build(exe):
    proc = subprocess.run([str(exe), "--version"], capture_output=True, text=True)
    text = proc.stdout + proc.stderr
    features = FEATURES.search(text)
    build = BUILDINFO.search(text)
    return {
        "features": features.group(1).strip() if features else "unknown",
        "build": build.group(1).strip() if build else "unknown (binary predates build info)",
        "host": platform.node(),
        "cpus": os.cpu_count(),
    }


def thread_ladder(cpus):
    ladder, n = [], 1
    while n <= cpus:
        ladder.append(n)
        n *= 2
    if cpus not in ladder:
        ladder.append(cpus)
    return ladder


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", default="build/mqc", help="path to the mqc binary")
    ap.add_argument("--threads", type=int, default=min(16, os.cpu_count() or 1))
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--basis", default="6-31g")
    ap.add_argument("--baseline", help="baseline file (default: per-host, beside this script)")
    ap.add_argument("--record", action="store_true", help="write the results as the new baseline")
    ap.add_argument("--quick", action="store_true", help="skip the thread ladder and the MPI case")
    ap.add_argument("--keep", help="keep the generated decks in this directory")
    args = ap.parse_args()

    exe = Path(args.exe).resolve()
    if not exe.exists():
        sys.exit(f"no binary at {exe} -- pass --exe")
    baseline_path = Path(args.baseline) if args.baseline else HERE / f"baseline-{platform.node()}.json"

    info = describe_build(exe)
    print("=" * 78)
    print(f"  host {info['host']}, {info['cpus']} cpus, running at {args.threads} threads")
    print(f"  {info['features']}")
    print(f"  {info['build']}")
    print("=" * 78)

    work = Path(args.keep) if args.keep else Path(tempfile.mkdtemp(prefix="mqc-bench-"))
    work.mkdir(parents=True, exist_ok=True)
    for xyz in GEOM.glob("*.xyz"):
        shutil.copy(xyz, work / xyz.name)
    results = {"info": info, "threads": args.threads, "basis": args.basis, "cases": {}}

    def say(*parts):
        # Unbuffered on purpose: this suite runs for many minutes and a user
        # watching a silent terminal cannot tell it from a hang.
        print(*parts, flush=True)

    def record(name, outcome, note=""):
        results["cases"][name] = outcome
        if "failed" in outcome:
            say(f"  {name:28s} FAILED  {outcome['failed']}")
        else:
            spread = (f"spread {outcome['spread_pct']:4.1f}%" if outcome.get("n_runs", 0) > 1
                      else "1 run, no spread")
            say(f"  {name:28s} {outcome['median']:8.2f} s  ({spread}) {note}")

    say("\nenergies, w10, one per rung")
    hf_time = None
    for stem, method, func, geom, basis, family in ENERGY_CASES:
        d = deck(work / f"e_{stem}.json", work / f"{geom}.xyz", method, func, "Energy", basis)
        out = repeat(exe, d, args.threads, args.repeats)
        note = f"[{family}]"
        if stem == "hf" and "median" in out:
            hf_time = out["median"]
        elif hf_time and "median" in out:
            note += f"  {out['median']/hf_time:.1f}x HF"
        record(f"energy/{stem}", out, note)

    say("\ncorrelated methods -- sized per method, see the table above")
    for stem, method, func, geom, basis, family in CORRELATED_CASES:
        safe = stem.replace("(", "_").replace(")", "")
        d = deck(work / f"c_{safe}.json", work / f"{geom}.xyz", method, func, "Energy", basis)
        record(f"correlated/{stem}", repeat(exe, d, args.threads, args.repeats),
               f"[{family}, {geom}/{basis}]")

    say("\ngradients, w10")
    for stem, method, func, geom, basis, family in GRADIENT_CASES:
        d = deck(work / f"g_{stem}.json", work / f"{geom}.xyz", method, func, "Gradient", basis)
        record(f"gradient/{stem}", repeat(exe, d, args.threads, max(2, args.repeats - 1)), f"[{family}]")

    if not args.quick:
        say("\nthread ladder, PBE energy -- a stage flat in the thread count is serial")
        d = deck(work / "ladder.json", work / "w10.xyz", "dft", "pbe", "Energy", args.basis)
        one = None
        for n in thread_ladder(info["cpus"] or 1):
            out = repeat(exe, d, n, max(2, args.repeats - 1))
            if "median" in out:
                one = one or out["median"]
                record(f"ladder/{n:02d}", out, f"{one/out['median']:5.1f}x on {n:2d} threads")
            else:
                record(f"ladder/{n:02d}", out)

        say("\nfragmented, w20 MBE(2) -- fragment work pins to one thread and spreads over ranks")
        frags = [[3 * i, 3 * i + 1, 3 * i + 2] for i in range(20)]
        d = deck(work / "frag.json", work / "w20.xyz", "dft", "pbe", "Energy", args.basis,
                 fragments=frags, level=2)
        record("fragmented/serial", repeat(exe, d, args.threads, 2))
        if shutil.which("mpirun"):
            record("fragmented/mpi4", repeat(exe, d, max(1, args.threads // 4), 2, ranks=4))
        else:
            print("  fragmented/mpi4              skipped, no mpirun")

    if baseline_path.exists():
        compare(results, json.loads(baseline_path.read_text()))
    else:
        print(f"\nno baseline for this host yet; --record writes one to {baseline_path.name}")
    if args.record:
        baseline_path.write_text(json.dumps(results, indent=2) + "\n")
        print(f"\nbaseline written to {baseline_path}")
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)


def band_for(outcome, old):
    """The width inside which nothing is claimed.

    Twice the measured spread, and never under five per cent, because the
    spread understates the noise it is standing in for. Repeats inside one
    invocation run back to back on a warm cache at a settled clock; a baseline
    recorded minutes or days earlier saw none of those conditions. Measured
    here, cases that repeat to under one per cent within a run still drift
    around three between runs -- enough to report a regression that is nothing
    but the machine, which is how a suite teaches people to ignore it.

    Real regressions are not marginal. Everything this suite was built from was
    tens of per cent or a factor.
    """
    return max(5.0, 2.0 * (outcome.get("spread_pct", 0.0) + old.get("spread_pct", 0.0)))


def compare(now, before):
    """Judge against the measured spread, never a fixed percentage.

    A case that repeats to 0.3 per cent and moves by 3 has regressed; one that
    repeats to 16 and moves by 3 has not moved at all. Using one threshold for
    both is how a suite becomes noise people stop reading.
    """
    print("\n" + "=" * 78)
    if before.get("info", {}).get("build") != now["info"]["build"]:
        print("  NOTE: the baseline was recorded from a different build --")
        print(f"    was: {before.get('info', {}).get('build')}")
        print(f"    now: {now['info']['build']}")
    print("  against baseline")
    print("=" * 78)
    verdicts = []
    for name, outcome in now["cases"].items():
        old = before.get("cases", {}).get(name)
        if not old or "median" not in outcome or "median" not in old:
            continue
        change = (outcome["median"] - old["median"]) / old["median"] * 100
        # Two runs, two spreads: the noise on a difference is the sum, and
        # doubling that is the band inside which nothing is claimed.
        band = band_for(outcome, old)
        if outcome.get("n_runs", 1) < 2 or old.get("n_runs", 1) < 2:
            # No spread was measured on one side, so there is nothing to judge
            # against. Ten per cent is a guess, and saying so beats reporting a
            # regression that is one machine hiccup.
            band = max(band, 10.0)
        if abs(change) <= band:
            continue
        verdicts.append((name, old["median"], outcome["median"], change,
                         "SLOWER" if change > 0 else "faster"))

    for name, outcome in now["cases"].items():
        old = before.get("cases", {}).get(name)
        if not old:
            continue
        for stage, seconds in outcome.get("stages", {}).items():
            was = old.get("stages", {}).get(stage)
            # Stages under a second are noise dressed as a measurement.
            if not was or was < 1.0:
                continue
            change = (seconds - was) / was * 100
            if abs(change) <= max(10.0, band_for(outcome, old)):
                continue
            verdicts.append((f"{name} :: {stage}", was, seconds, change,
                             "SLOWER" if change > 0 else "faster"))
    if not verdicts:
        print("  nothing outside the measured spread")
    for name, was, now_s, change, word in sorted(verdicts, key=lambda v: -abs(v[3])):
        print(f"  {name:28s} {was:8.2f} -> {now_s:8.2f} s  {change:+6.1f}%  {word}")
    if any(v[4] == "SLOWER" for v in verdicts):
        print("\n  regressions above; re-run to confirm before believing a single result")


if __name__ == "__main__":
    main()
