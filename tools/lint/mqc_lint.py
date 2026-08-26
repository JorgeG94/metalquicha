#!/usr/bin/env python3
"""House rules for metalquicha's Fortran that no other linter covers.

fprettify enforces layout and fortitude enforces language-level hygiene.
Neither has an opinion about the two things this checks:

  MQC001  a comment marker other than ``!`` or ``!!``
  MQC002  a parameter whose name already exists in mqc_physical_constants
  MQC003  a parameter whose value already exists there under another name

MQC002 and MQC003 exist because mqc_physical_constants says in its own header
that it is the only place a constant may be written down -- the Bohr radius had
already drifted to three different values across the tree before that rule
existed, which is a class of bug no test catches, because every copy is
self-consistent.

Run over specific files (this is how pre-commit calls it)::

    tools/lint/mqc_lint.py src/scf/mqc_scf_common.f90

or over the whole tree::

    tools/lint/mqc_lint.py --all

Suppress a single finding by putting ``mqc-lint: ignore MQC003`` in a comment
on the offending line or on the line above it, with a reason.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

CONSTANTS_MODULE = Path("src/core/mqc_physical_constants.F90")

SCAN_DIRS = ("src", "app", "backends", "test", "validation", "tools")
SUFFIXES = (".f90", ".F90", ".f", ".F")

# Compiler directives, not comments. `!$` alone is a conditional-compilation
# sentinel, `!$omp`/`!$acc` are OpenMP/OpenACC, `!dir$`/`!gcc$` are vendor.
DIRECTIVE = re.compile(r"^[$]|^(dir|gcc|ibm|dec)[$]", re.I)

# FORD's alternative doc markers. `>` is its predocmark, `*` and `|` its
# alternates. The house style uses none of them.
BAD_MARKERS = {">": "predocmark", "|": "alternate predocmark"}

SUPPRESS = re.compile(r"mqc-lint:\s*ignore\s+(MQC\d{3}(?:\s*,\s*MQC\d{3})*)", re.I)

# Values too common to be evidence of a copied constant. A parameter that
# happens to equal 2.0 is not a duplicated physical constant.
BORING = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 10.0, 100.0, 1000.0}

DECL_START = re.compile(
    r"""^\s*(
        (real|integer|complex|logical|character|double\s+precision)\s*[(,:]
      | (type|class)\s*\(
      | (module|program|submodule)\s+\w
      | (pure\s+|elemental\s+|recursive\s+|impure\s+|module\s+)*(subroutine|function)\s
      | [\w()=*, ]*\bfunction\s+\w+\s*\(
      | interface\b | abstract\s+interface\b
      | type\s*(,|::|\s+\w) | end\s+type\b | contains\b | enum\b
      | use\s+\w | implicit\s+none | procedure\b | generic\b | enumerator\b
    )""",
    re.X | re.I,
)


class Finding:
    def __init__(self, path, line, rule, message):
        self.path, self.line, self.rule, self.message = path, line, rule, message

    def render(self, fmt):
        if fmt == "github":
            return (
                f"::error file={self.path},line={self.line},"
                f"title={self.rule}::{self.message}"
            )
        return f"{self.path}:{self.line}: {self.rule} {self.message}"


def comment_column(line: str) -> int:
    """Index of the ``!`` that opens a comment, or -1.

    Quotes are tracked so that a bang inside a string literal -- a format
    string, an error message -- is not mistaken for one.
    """
    quote = None
    i = 0
    while i < len(line):
        c = line[i]
        if quote is not None:
            if c == quote:
                if i + 1 < len(line) and line[i + 1] == quote:
                    i += 2
                    continue
                quote = None
        elif c in "'\"":
            quote = c
        elif c == "!":
            return i
        i += 1
    return -1


def code_before_comment(line: str) -> str:
    col = comment_column(line)
    return line if col < 0 else line[:col]


def suppressed(lines, idx, rule) -> bool:
    for probe in (idx, idx - 1):
        if 0 <= probe < len(lines):
            m = SUPPRESS.search(lines[probe])
            if m and rule.upper() in [r.strip().upper() for r in m.group(1).split(",")]:
                return True
    return False


# --------------------------------------------------------------------------
# MQC001 -- comment markers
# --------------------------------------------------------------------------


def next_code_line(lines, i):
    """The first line after ``i`` that is neither blank, comment nor directive.

    Preprocessor lines are stepped over: a docstring in a ``.F90`` file often
    sits above an ``#if`` that selects between two declarations, and it
    documents the declaration rather than the conditional.
    """
    j = i
    while j < len(lines):
        s = lines[j].strip()
        if s and not s.startswith("!") and not s.startswith("#"):
            return lines[j]
        j += 1
    return ""


def check_comment_markers(path, lines, allow_predoc):
    findings = []
    for i, line in enumerate(lines):
        col = comment_column(line)
        if col < 0:
            continue
        rest = line[col:].lstrip("!")
        bangs = len(line[col:]) - len(line[col:].lstrip("!"))
        if bangs and DIRECTIVE.match(rest):
            continue
        if bangs > 2:
            marker = "!" * bangs
            if not suppressed(lines, i, "MQC001"):
                findings.append(
                    Finding(
                        path,
                        i + 1,
                        "MQC001",
                        f"comment marker '{marker}': use '!' for an inline comment "
                        f"or '!!' to document an entity",
                    )
                )
            continue
        if not rest or rest[0] not in BAD_MARKERS:
            continue

        marker = "!" * bangs + rest[0]
        # A block of these documents whatever follows it; find that.
        block_end = i
        while block_end < len(lines):
            c = comment_column(lines[block_end])
            if c < 0 or lines[block_end][c:].lstrip("!")[:1] not in ("", rest[0]):
                break
            block_end += 1
        documents_decl = bool(DECL_START.match(next_code_line(lines, block_end)))

        if allow_predoc and documents_decl:
            continue
        if suppressed(lines, i, "MQC001"):
            continue

        if documents_decl:
            fix = "move the text below the declaration and use '!!'"
        else:
            fix = "this comment documents no entity, so use '!'"
        findings.append(
            Finding(path, i + 1, "MQC001", f"comment marker '{marker}': {fix}")
        )
        # Report the block once, not once per line.
        for skip in range(i + 1, block_end):
            lines[skip] = lines[skip]  # no-op, kept for clarity
    return dedupe_blocks(findings)


def fix_comment_markers(lines):
    """Rewrite inline ``!>`` to ``!``; return (new_lines, count).

    Only the inline case is rewritten. A marker that documents a declaration
    cannot be repaired by swapping characters -- ``!!`` above a declaration
    attaches to whatever precedes it -- so that one is reported for a human to
    move, never auto-fixed.
    """
    out, changed = list(lines), 0
    i = 0
    while i < len(out):
        col = comment_column(out[i])
        if col < 0:
            i += 1
            continue
        rest = out[i][col:].lstrip("!")
        bangs = len(out[i][col:]) - len(out[i][col:].lstrip("!"))
        if bangs != 1 or not rest or rest[0] not in BAD_MARKERS or DIRECTIVE.match(rest):
            i += 1
            continue
        block_end = i
        while block_end < len(out):
            c = comment_column(out[block_end])
            if c < 0 or out[block_end][c:].lstrip("!")[:1] not in ("", rest[0]):
                break
            block_end += 1
        if DECL_START.match(next_code_line(out, block_end)):
            i = block_end
            continue
        for k in range(i, block_end):
            c = comment_column(out[k])
            head, tail = out[k][:c], out[k][c:]
            body = tail.lstrip("!")[1:]
            out[k] = head + "!" + (" " + body.lstrip() if body.strip() else "")
            changed += 1
        i = block_end
    return out, changed


def dedupe_blocks(findings):
    """Collapse consecutive MQC001 lines into the first of each block."""
    out, prev = [], {}
    for f in findings:
        key = (f.path, f.rule)
        if key in prev and f.line == prev[key] + 1:
            prev[key] = f.line
            continue
        prev[key] = f.line
        out.append(f)
    return out


# --------------------------------------------------------------------------
# MQC002 / MQC003 -- constants
# --------------------------------------------------------------------------

NUM = re.compile(
    r"^[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?(_\w+)?$"
)

PARAM_DECL = re.compile(
    r"^\s*(?:real|integer|complex|double\s+precision)\b[^:]*\bparameter\b[^:]*::\s*(.+?)\s*$",
    re.I,
)


def to_float(text: str):
    t = text.strip()
    t = re.sub(r"_\w+$", "", t)
    t = re.sub(r"([\d.])[dD]([+-]?\d)", r"\1e\2", t)
    if not NUM.match(t + ""):
        return None
    try:
        return float(t)
    except ValueError:
        return None


def split_top_level(text: str, sep=","):
    """Split on ``sep`` outside parentheses and string literals."""
    parts, depth, quote, cur = [], 0, None, []
    for c in text:
        if quote is not None:
            if c == quote:
                quote = None
        elif c in "'\"":
            quote = c
        elif c in "([":
            depth += 1
        elif c in ")]":
            depth -= 1
        elif c == sep and depth == 0:
            parts.append("".join(cur))
            cur = []
            continue
        cur.append(c)
    parts.append("".join(cur))
    return [p.strip() for p in parts if p.strip()]


def load_constants(root: Path):
    """Name -> value for every public parameter of mqc_physical_constants."""
    path = root / CONSTANTS_MODULE
    consts = {}
    if not path.exists():
        return consts
    for line in path.read_text(errors="replace").split("\n"):
        code = code_before_comment(line)
        m = re.match(
            r"^\s*(?:real|integer)\b[^:]*\bparameter\b[^:]*\bpublic\b[^:]*::\s*(.+)$",
            code,
            re.I,
        )
        if not m:
            continue
        for item in split_top_level(m.group(1)):
            if "=" not in item:
                continue
            name, _, value = item.partition("=")
            consts.setdefault(name.strip().upper(), to_float(value))
    return consts


def significant_digits(text: str) -> int:
    digits = re.sub(r"_\w+$", "", text.strip())
    digits = re.sub(r"[eEdD][+-]?\d+$", "", digits)
    return len(re.sub(r"[^0-9]", "", digits).lstrip("0"))


def check_constants(path, lines, consts, is_constants_module):
    if is_constants_module:
        return []
    findings = []
    by_value = [(n, v) for n, v in consts.items() if v is not None and v not in BORING]

    for i, line in enumerate(lines):
        code = code_before_comment(line)

        # PI written as an expression rather than a literal.
        if re.search(r"\bacos\s*\(\s*-\s*1", code, re.I) and re.search(
            r"\bparameter\b", code, re.I
        ):
            if not suppressed(lines, i, "MQC003"):
                findings.append(
                    Finding(
                        path,
                        i + 1,
                        "MQC003",
                        "acos(-1) is PI; use mqc_physical_constants, only: PI",
                    )
                )
            continue

        m = PARAM_DECL.match(code)
        if not m:
            continue
        for item in split_top_level(m.group(1)):
            if "=" not in item:
                continue
            name, _, value = item.partition("=")
            name = name.strip().upper()
            value = value.strip()

            # A name match short-circuits the value check whether or not it
            # was reported, so that silencing MQC002 does not leave MQC003
            # firing on the same declaration.
            if name in consts:
                if not suppressed(lines, i, "MQC002"):
                    findings.append(
                        Finding(
                            path,
                            i + 1,
                            "MQC002",
                            f"'{name}' is already in mqc_physical_constants; "
                            f"use it from there instead of redeclaring it",
                        )
                    )
                continue

            num = to_float(value)
            if num is None or num in BORING or significant_digits(value) < 5:
                continue
            for cname, cval in by_value:
                if cval == 0.0:
                    continue
                if abs(num - cval) <= 1e-9 * abs(cval):
                    if not suppressed(lines, i, "MQC003"):
                        findings.append(
                            Finding(
                                path,
                                i + 1,
                                "MQC003",
                                f"'{name}' has the value of {cname} in "
                                f"mqc_physical_constants; take it from there so the "
                                f"two cannot drift apart",
                            )
                        )
                    break
    return findings


# --------------------------------------------------------------------------


def gather(root: Path, explicit):
    if explicit:
        return [Path(p) for p in explicit]
    files = []
    for d in SCAN_DIRS:
        base = root / d
        if not base.is_dir():
            continue
        for p in sorted(base.rglob("*")):
            if p.suffix in SUFFIXES and "build" not in p.parts and "_deps" not in p.parts:
                files.append(p)
    return files


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("files", nargs="*", help="files to check")
    ap.add_argument("--all", action="store_true", help="check the whole tree")
    ap.add_argument(
        "--rules",
        default="MQC001,MQC002,MQC003",
        help="comma-separated rules to apply (default: all)",
    )
    ap.add_argument(
        "--allow-predoc",
        action="store_true",
        help="under MQC001, permit '!>' where it documents a declaration and "
        "flag it only where it precedes executable code",
    )
    ap.add_argument(
        "--fix",
        action="store_true",
        help="rewrite inline '!>' to '!' in place; markers that document a "
        "declaration are reported, not touched",
    )
    ap.add_argument("--format", choices=("text", "github"), default="text")
    ap.add_argument(
        "--stats", action="store_true", help="print a per-rule count and exit 0"
    )
    args = ap.parse_args(argv)

    root = Path(os.environ.get("MQC_ROOT", ".")).resolve()
    # Walk up to the repo root so the constants module is found no matter
    # where the hook is invoked from.
    probe = root
    while not (probe / CONSTANTS_MODULE).exists() and probe != probe.parent:
        probe = probe.parent
    if (probe / CONSTANTS_MODULE).exists():
        root = probe

    rules = {r.strip().upper() for r in args.rules.split(",") if r.strip()}
    consts = load_constants(root)
    if not consts:
        print(
            f"mqc-lint: could not read {CONSTANTS_MODULE}; MQC002/MQC003 disabled",
            file=sys.stderr,
        )

    findings = []
    for path in gather(root, [] if args.all else args.files):
        try:
            lines = path.read_text(errors="replace").split("\n")
        except OSError:
            continue
        rel = os.path.relpath(path, root)
        if args.fix and "MQC001" in rules:
            fixed, n = fix_comment_markers(lines)
            if n:
                path.write_text("\n".join(fixed))
                lines = fixed
                print(f"{rel}: rewrote {n} inline comment marker(s)")
        if "MQC001" in rules:
            findings += check_comment_markers(rel, lines, args.allow_predoc)
        if rules & {"MQC002", "MQC003"}:
            same = Path(rel).as_posix() == CONSTANTS_MODULE.as_posix()
            findings += [
                f for f in check_constants(rel, lines, consts, same) if f.rule in rules
            ]

    if args.stats:
        counts = {}
        for f in findings:
            counts[f.rule] = counts.get(f.rule, 0) + 1
        for rule in sorted(set(list(counts) + sorted(rules))):
            print(f"{rule}: {counts.get(rule, 0)}")
        return 0

    findings.sort(key=lambda f: (f.path, f.line))
    for f in findings:
        print(f.render(args.format))
    if findings:
        print(f"\nmqc-lint: {len(findings)} finding(s)", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
