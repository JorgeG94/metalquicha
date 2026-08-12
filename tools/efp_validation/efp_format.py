#!/usr/bin/env python3
"""Parse and render GAMESS ``.efp`` effective fragment potential files.

Test scaffolding, not a shipped product: it exists so a reimplementation of
GAMESS's ``RUNTYP=MAKEFP`` can be compared against GAMESS parameter by
parameter. Nothing in ``python/`` imports it and nothing packages it.

Two directions, wanted at different times:

* **parse** -- a ``.efp`` into a plain dict. This is the half that does the work:
  comparing our JSON output against GAMESS means putting both into the same
  structure and diffing numbers with per-quantity tolerances, which is a very
  different thing from diffing fixed-width text where a trailing zero is a
  failure and nothing tells you which parameter moved.
* **render** -- a dict back into ``.efp`` text, so GAMESS can be handed a
  potential we generated. Only needed for the final end-to-end check, and only
  has to be good enough for GAMESS's reader.

Run it directly to round-trip every reference potential GAMESS ships::

    python3 tools/efp_validation/efp_format.py [gamess_root]

The round-trip is the test: parse, render, parse again, and require the two
structures to be identical, plus require that every non-blank input line was
consumed by some section. The first catches a renderer that drops a field; the
second catches a *parser* that skipped a section it did not recognise, which is
the failure that would otherwise look like success.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

#: Sections holding one record per expansion point, and how many numbers each
#: record carries. The counts are the packed GAMESS conventions, not the full
#: Cartesian tensors: 6 quadrupole components and 10 octopole ones, against the
#: 9 and 27 libcint returns. That packing map is a real source of silent
#: transposition and belongs written down somewhere.
RECORD_SECTIONS = {
    "COORDINATES": 5,      # x y z mass nuclear_charge
    "MONOPOLES": 2,        # electronic, nuclear
    "DIPOLES": 3,
    "QUADRUPOLES": 6,
    "OCTUPOLES": 10,
    "LMO CENTROIDS": None,  # x y z (+ trailing fields in some files)
}

#: Sections whose records are a labelled point followed by a 3x3 tensor.
TENSOR_SECTIONS = ("POLARIZABLE POINTS", "DYNAMIC POLARIZABLE POINTS")

#: Sections kept as raw text. The projection basis and the wavefunction blocks
#: are bulk data we only ever pass through, so parsing them into structure would
#: be work with no consumer -- until M7 wants to emit them, and even then a
#: faithful copy is what is required.
RAW_SECTIONS = ("PROJECTION BASIS SET", "PROJECTION WAVEFUNCTION",
                "FOCK MATRIX ELEMENTS",
                # Present in real MAKEFP output but *not* in the shipped library,
                # because the library files were curated. A default run has
                # DISP7, DISP8 and CHTR on, so it emits these -- and a parser
                # built only against auxdata/EFP silently dropped 1347 of 2247
                # lines of a real file. Kept raw: bulk data with no consumer yet.
                "DIPOLE-QUADRUPOLE DYNAMIC POLARIZABLE POINTS",
                "LMOQQPOL DYNAMIC POLARIZABLE POINTS",
                "LMODOPOL DYNAMIC POLARIZABLE POINTS",
                "MOLECULAR DIPOLE POLARIZABILITY",
                "MOLECULAR DIP-QUAD POLARIZABILITY",
                "MOLECULAR DIP-OCTO POLARIZABILITY",
                "MOLECULAR QUAD-QUAD POLARIZABILITY",
                "CTVEC", "CTFOK")

#: The imaginary frequency is written into each dynamic record as a trailing
#: comment: ``-- FOR W= 0.002792I A.U.``. That is how the twelve Casimir-Polder
#: quadrature points can be read off a reference file rather than guessed.
_FREQ = re.compile(r"--\s*FOR\s+W=\s*([0-9.]+)I", re.I)

_NUM = re.compile(r"[-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?|[-+]?\d+")


def _numbers(text: str) -> list[float]:
    """Every number in ``text``, with Fortran's D exponent accepted."""
    return [float(t.replace("D", "E").replace("d", "e")) for t in _NUM.findall(text)]


def _split_label(line: str) -> tuple[str, str]:
    """A record's label and its numeric remainder.

    Labels contain digits (``A01N1``, ``BO21``), so the numbers cannot simply be
    scraped off the whole line -- ``A01N1`` would contribute a 01 and a 1. The
    label is therefore taken as the leading token, with one exception: the
    dynamic polarizability records write theirs as ``CT  2``, with the index
    split off by whitespace, where every other section writes ``CT2``.
    """
    parts = line.split()
    if not parts:
        return "", ""
    if parts[0].upper() == "CT" and len(parts) > 1 and parts[1].isdigit():
        return f"CT{parts[1]}", line.split(parts[1], 1)[1]
    return parts[0], line[len(parts[0]):]


def _records(lines: list[str]) -> list[tuple[str, str]]:
    """Group a section's lines into (label, payload) records.

    A record begins on a line whose first column is non-blank; anything indented
    continues it. That rule is what the file actually obeys, and it is more
    reliable than the trailing ``>`` continuation marker -- the polarizability
    sections put a point's tensor on *unmarked* following lines, so a
    ``>``-driven parser silently truncates them to their coordinates.
    """
    out: list[tuple[str, str]] = []
    for line in lines:
        if not line.strip():
            continue
        payload = line.rstrip().rstrip(">").rstrip()
        if line[0] not in " \t":
            label, rest = _split_label(payload)
            out.append((label, rest))
        elif out:
            out[-1] = (out[-1][0], out[-1][1] + " " + payload)
    return out


def parse_efp(text: str) -> dict:
    """A GAMESS ``.efp`` file as a dict."""
    lines = text.split("\n")
    doc: dict = {"sections": {}}
    consumed = 0

    # Header: a $-prefixed fragment name, then a free-text comment naming the
    # bases. The comment matters for provenance -- the shipped potentials use
    # one basis for electrostatics and a larger one for exchange-repulsion, and
    # a comparison against parameters made in a different basis is meaningless.
    # A generated file opens with two banner lines ("RUNTYP=MAKEFP ... DATA
    # FOLLOWS", "<NAME> EFP GENERATED AT <time>") before the $-name; the curated
    # library files start at the $-name directly. Skip anything up to it rather
    # than assuming either shape.
    i = 0
    banner: list[str] = []
    while i < len(lines) and not lines[i].strip().startswith("$"):
        if lines[i].strip():
            banner.append(lines[i].rstrip())
            consumed += 1
        i += 1
    if banner:
        doc["banner"] = banner
    if i < len(lines) and lines[i].strip().startswith("$"):
        doc["name"] = lines[i].strip().lstrip("$").strip()
        consumed += 1
        i += 1
    if i < len(lines) and not _is_header(lines[i]):
        doc["comment"] = lines[i].strip()
        consumed += 1
        i += 1

    current: str | None = None
    body: list[str] = []

    def close() -> None:
        nonlocal current, body
        if current is not None:
            doc["sections"][current] = _finish(current, body)
        current, body = None, []

    for line in lines[i:]:
        if not line.strip():
            continue
        stripped = line.strip()
        if stripped.upper() == "STOP":
            consumed += 1
            close()
            continue
        # The document terminator. Easy to miss, since it looks like the $-name
        # on line one and sits after the last STOP where nothing else appears --
        # which is exactly how it was missed until the unconsumed-line count was
        # fixed to only credit lines actually filed somewhere.
        if stripped.upper() == "$END":
            consumed += 1
            close()
            doc["terminated"] = True
            continue
        header = _header_name(line)
        if header is not None:
            consumed += 1
            close()
            current = header
            # MULTIPLICITY and PROJECTION WAVEFUNCTION carry values on the
            # header line itself rather than below it.
            trailing = _numbers(stripped[len(header):])
            if trailing:
                doc.setdefault("header_values", {})[header] = trailing
            continue
        if current is not None:
            consumed += 1
            body.append(line)
    close()

    doc["_consumed_lines"] = consumed
    doc["_nonblank_lines"] = sum(1 for l in lines if l.strip())
    return doc


#: Header keywords, longest first so "DYNAMIC POLARIZABLE POINTS" is not
#: mistaken for "POLARIZABLE POINTS".
_HEADERS = sorted(
    list(RECORD_SECTIONS) + list(TENSOR_SECTIONS) + list(RAW_SECTIONS)
    + ["MULTIPLICITY", "SCREEN2", "SCREEN", "CHARGE TRANSFER"],
    key=len, reverse=True,
)


def _header_name(line: str) -> str | None:
    """The section this line opens, or None if it is a record.

    Matched on the keyword alone, deliberately *not* on indentation. Most
    headers are indented by one column and records start in column 1, which
    makes indentation look like a usable rule -- but ``SCREEN2`` starts in
    column 1 like a record does, and an indentation rule drops it. No expansion
    point label (``A01N1``, ``BO21``, ``CT1``) begins with a section keyword, so
    keyword matching alone is unambiguous.
    """
    upper = line.strip().upper()
    for key in _HEADERS:
        if upper.startswith(key):
            return key
    return None


def _is_header(line: str) -> bool:
    return _header_name(line) is not None


def _finish(section: str, body: list[str]) -> object:
    """Turn a section's raw lines into structure."""
    if section in RAW_SECTIONS:
        return {"raw": [l.rstrip() for l in body]}

    if section in TENSOR_SECTIONS:
        points = []
        for label, payload in _records(body):
            freq = _FREQ.search(payload)
            values = _numbers(_FREQ.sub("", payload))
            points.append({
                "label": label,
                "xyz": values[:3],
                # Nine components: the full 3x3, not a symmetrised six.
                "tensor": values[3:],
                **({"frequency": float(freq.group(1))} if freq else {}),
            })
        return points

    records = []
    for label, payload in _records(body):
        records.append({"label": label, "values": _numbers(payload)})
    return records


def render_efp(doc: dict) -> str:
    """A dict back into ``.efp`` text.

    Aimed at GAMESS's reader, not at byte-identity with its writer: free-form
    numbers with the continuation marker where the sections expect one. Byte
    fidelity would be a much larger job for no extra confidence, since the test
    that matters is whether GAMESS accepts the file and agrees with it.
    """
    out: list[str] = list(doc.get("banner", []))
    out.append(f" ${doc.get('name', 'FRAG')}")
    if "comment" in doc:
        out.append(doc["comment"])

    for section, content in doc["sections"].items():
        header = f" {section}"
        for value in doc.get("header_values", {}).get(section, []):
            header += f"   {int(value) if float(value).is_integer() else value}"
        out.append(header)

        if isinstance(content, dict) and "raw" in content:
            out.extend(content["raw"])
        elif section in TENSOR_SECTIONS:
            for point in content:
                tag = point["label"]
                suffix = ""
                if "frequency" in point:
                    tag = f"CT  {tag[2:]}"
                    suffix = f" -- FOR W= {point['frequency']:.6f}I A.U."
                coords = "".join(f"{v:16.10f}" for v in point["xyz"])
                out.append(f"{tag}{coords}{suffix}")
                out.extend(_wrap(point["tensor"]))
        else:
            for record in content:
                values = "".join(f"{v:16.10f}" for v in record["values"])
                out.extend(_wrap_record(record["label"], values, record["values"]))
        out.append(" STOP")

    if doc.get("terminated"):
        out.append(" $END")
    return "\n".join(out) + "\n"


def _wrap(values: list[float], per_line: int = 4) -> list[str]:
    """Continuation lines, ``>`` on every line but the last."""
    lines = []
    for start in range(0, len(values), per_line):
        chunk = values[start:start + per_line]
        text = "    " + "".join(f"{v:16.10f}" for v in chunk)
        lines.append(text + (" >" if start + per_line < len(values) else ""))
    return lines


def _wrap_record(label: str, _values_text: str, values: list[float],
                 per_line: int = 4) -> list[str]:
    lines = []
    for start in range(0, len(values), per_line):
        chunk = values[start:start + per_line]
        text = "".join(f"{v:16.10f}" for v in chunk)
        prefix = f"{label:<10s}" if start == 0 else " " * 10
        lines.append(prefix + text + (" >" if start + per_line < len(values) else ""))
    return lines


# --------------------------------------------------------------------------
# the round-trip test
# --------------------------------------------------------------------------

def _compare(a: object, b: object, path: str = "") -> list[str]:
    """Differences between two parsed documents, ignoring bookkeeping keys."""
    problems = []
    if isinstance(a, dict) and isinstance(b, dict):
        keys = {k for k in set(a) | set(b) if not k.startswith("_")}
        for key in sorted(keys):
            if key not in a or key not in b:
                problems.append(f"{path}/{key}: present in only one")
            else:
                problems += _compare(a[key], b[key], f"{path}/{key}")
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            problems.append(f"{path}: length {len(a)} vs {len(b)}")
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                problems += _compare(x, y, f"{path}[{i}]")
    elif isinstance(a, float) and isinstance(b, float):
        if abs(a - b) > 1e-9:
            problems.append(f"{path}: {a} vs {b}")
    elif a != b:
        problems.append(f"{path}: {a!r} vs {b!r}")
    return problems


def main(argv: list[str]) -> int:
    root = Path(argv[1]) if len(argv) > 1 else Path("../mgga/gamess")
    efp_dir = root / "auxdata" / "EFP"
    files = sorted(efp_dir.glob("*.efp"))
    if not files:
        print(f"no .efp files under {efp_dir}")
        return 1

    failures = 0
    frequencies: set[float] = set()
    for path in files:
        text = path.read_text()
        doc = parse_efp(text)
        again = parse_efp(render_efp(doc))
        problems = _compare(doc, again)

        # Did the parser file everything somewhere? Counted only for lines
        # actually assigned to a section, which is the point: an unrecognised
        # header leaves its records belonging to nothing, and the round-trip
        # alone cannot see them because they are absent from both sides. The
        # first version of this counted every line it looked at and therefore
        # reported zero while quietly dropping the whole SCREEN2 block.
        missed = doc["_nonblank_lines"] - doc["_consumed_lines"]

        for point in doc["sections"].get("DYNAMIC POLARIZABLE POINTS", []):
            if "frequency" in point:
                frequencies.add(point["frequency"])

        status = "ok  "
        if problems or missed:
            status, failures = "FAIL", failures + 1
        print(f"  {status} {path.name:<12s} "
              f"{len(doc['sections'])} sections  "
              f"{doc['_nonblank_lines']:5d} lines  "
              f"unconsumed {missed}"
              + (f"  {problems[:2]}" if problems else ""))

    print(f"\n  {len(files) - failures}/{len(files)} round-tripped")
    if frequencies:
        ordered = sorted(frequencies)
        print(f"  {len(ordered)} distinct imaginary frequencies: "
              + ", ".join(f"{f:.6f}" for f in ordered))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
