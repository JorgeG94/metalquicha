#!/usr/bin/env python3
"""Pre-commit hook: an SCF call must forward the settings its caller holds.

metalquicha has one SCF driver and about two dozen callers, and every
caller passes its arguments by hand.  Forgetting one compiles cleanly,
passes every test, and produces a run that behaves as though the deck had
never mentioned the setting -- the schema accepted it, the reader read it,
the config carried it, and then it stopped at a call site.

That failure is silent by construction, and it has happened repeatedly:

  * the Fukui deltaSCF ions dropped six of thirteen documented keys;
  * the MakeFP SCF dropped six settings and hardcoded its iteration cap;
  * `climb_basis_ladder` never received the settings at all, so every rung
    of a basis-set-projection guess ran on defaults -- including
    `linear_dependence`, on the one path that exists to climb into a large
    diffuse basis where the overlap goes near-singular.

Each was found by reading code, which does not scale and did not catch the
next one.  This does: it reads every `call run_libcint_rhf` / `_uhf` in the
tree and reports which of the optional settings the call forwards.

## What counts as a violation

A call site that forwards fewer settings than the ones listed in SETTINGS,
unless it is either

  * marked with `! mqc: scf-subset -- <reason>` on the call or just above
    it, for a site where forwarding genuinely makes no sense; or
  * listed in the baseline beside this file, which is the ratchet.

## The baseline

The tree was NOT at zero when this was written -- FMO, SAPT and AFO drop
everything, and fixing them is its own work.  So known sites are recorded
in `scf_call_sites_baseline.txt` with a count, and this fails only on a
site that is new or that got worse.  The baseline is meant to shrink: when
a site is fixed, drop its line.  It is not a place to add to.

Receives candidate files as args from pre-commit.  Also runs manually:

    tools/lint/scf_call_sites.py --all
    tools/lint/scf_call_sites.py --all --report   # the inventory, no verdict
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

SCAN_DIRS = ("src", "app", "backends")

#: The optional arguments that describe how an SCF is driven.  A caller
#: holding a configuration should forward all of them; `max_iter` and the
#: two tolerances are positional and always passed, so they are not here.
SETTINGS = (
    "level_shift",
    "linear_dependence",
    "accelerator",
    "diis_vectors",
    "incremental_fock",
    "grad_tol",
    "convergence",
)

CALL = re.compile(r"\bcall\s+run_libcint_(?:rhf|uhf)\s*\(", re.IGNORECASE)
EXEMPT = re.compile(r"!\s*mqc:\s*scf-subset\b", re.IGNORECASE)

BASELINE = Path(__file__).with_name("scf_call_sites_baseline.txt")


def statement(text: str, start: int) -> tuple[str, int]:
    """The whole call, from `call` to its matching paren, and its end offset."""
    depth, i = 0, start
    while i < len(text):
        c = text[i]
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
            if depth == 0:
                return text[start : i + 1], i
        i += 1
    return text[start:], len(text)


def exempt_near(lines: list[str], line_no: int, seg: str) -> bool:
    """A marker on the call itself, or on any comment line just above it."""
    if EXEMPT.search(seg):
        return True
    i = line_no - 2  # 0-based, the line above
    while i >= 0:
        stripped = lines[i].strip()
        if not stripped.startswith("!"):
            return False
        if EXEMPT.search(stripped):
            return True
        i -= 1
    return False


def scan(paths: list[Path]) -> list[tuple[str, int, list[str]]]:
    """Every call site, as (path, line, missing settings)."""
    found = []
    for p in sorted(paths):
        try:
            text = p.read_text()
        except (OSError, UnicodeDecodeError):
            continue
        lines = text.splitlines()
        for m in CALL.finditer(text):
            seg, _ = statement(text, m.start())
            line_no = text.count("\n", 0, m.start()) + 1
            if exempt_near(lines, line_no, seg):
                continue
            # `scf=` is the whole configuration in one argument, and is what a
            # call site should pass. It satisfies the rule on its own: the
            # individual keywords predate it and remain valid, but repeating
            # six of them at every site is the shape that lost settings.
            if re.search(r"\bscf\s*=", seg):
                continue
            missing = [s for s in SETTINGS if not re.search(rf"\b{s}\s*=", seg)]
            if missing:
                found.append((p.as_posix(), line_no, missing))
    return found


def read_baseline() -> dict[str, int]:
    """path -> how many incomplete call sites it is allowed to still have."""
    allowed: dict[str, int] = {}
    if not BASELINE.exists():
        return allowed
    for raw in BASELINE.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        path, _, count = line.rpartition(" ")
        allowed[path.strip()] = int(count)
    return allowed


def main(argv: list[str]) -> int:
    report = "--report" in argv
    argv = [a for a in argv if a != "--report"]

    if "--all" in argv:
        paths = [
            p
            for d in SCAN_DIRS
            for p in Path(d).rglob("*")
            if p.suffix in (".f90", ".F90")
        ]
    else:
        # Only the directories that hold production call sites. A test or a
        # validation program IS the caller -- it constructs exactly the SCF it
        # means to run and has no configuration to forward -- so requiring it
        # to pass a full set would be noise, and marking every one of them
        # exempt would be worse.
        paths = [
            p
            for a in argv[1:]
            if a.endswith((".f90", ".F90"))
            for p in [Path(a)]
            if p.parts and p.parts[0] in SCAN_DIRS
        ]
    if not paths:
        return 0

    found = scan(paths)

    if report:
        if not found:
            print("every SCF call site forwards its caller's settings")
            return 0
        print(f"{len(found)} SCF call site(s) forwarding an incomplete configuration:\n")
        for path, line, missing in found:
            print(f"  {path}:{line}")
            print(f"      missing: {', '.join(missing)}")
        print("\nCounts per file, in baseline format:\n")
        per: dict[str, int] = {}
        for path, _, _ in found:
            per[path] = per.get(path, 0) + 1
        for path in sorted(per):
            print(f"{path} {per[path]}")
        return 0

    # Ratchet: a file may keep the count the baseline records, no more.
    allowed = read_baseline()
    per: dict[str, list[tuple[int, list[str]]]] = {}
    for path, line, missing in found:
        per.setdefault(path, []).append((line, missing))

    failed = False
    for path in sorted(per):
        hits = per[path]
        budget = allowed.get(path, 0)
        if len(hits) <= budget:
            continue
        failed = True
        print(
            f"{path}: {len(hits)} SCF call site(s) forward an incomplete "
            f"configuration, but only {budget} are recorded in the baseline."
        )
        for line, missing in hits:
            print(f"  {path}:{line}  missing: {', '.join(missing)}")

    if failed:
        print(
            "\nAn SCF call that does not forward its caller's settings makes those\n"
            "settings silently inert: the deck is read, validated, carried, and then\n"
            "dropped here. Forward them, or mark the call\n"
            "\n"
            "    ! mqc: scf-subset -- <why this SCF should not inherit them>\n"
            "\n"
            f"The baseline in {BASELINE.name} is meant to shrink. Do not add to it."
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
