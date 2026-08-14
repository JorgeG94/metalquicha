#!/usr/bin/env python3
"""Zero chosen sections of a potential, so one energy term can be isolated.

    python3 tools/efp_validation/zero_sections.py in.efp --zero DIPOLES QUADRUPOLES \
        OCTUPOLES SCREEN SCREEN2 > monopoles_only.efp

**Why this exists.** GAMESS reports one number for electrostatics: the sum over
every pair of expansion points of every multipole interaction, damped. Comparing
our implementation against that number tests the whole multipole expansion at
once, and when it disagrees there is nothing to say which term is wrong -- the
quadrupole component order, the octupole normalization, the damping form, a sign.

Zeroing sections turns one comparison into a ladder. With only monopoles left,
GAMESS's electrostatic energy is exactly the charge-charge sum, which has no
conventions in it to get wrong. Add dipoles and the difference is the
charge-dipole and dipole-dipole terms and nothing else. Each rung pins one thing.

**Zeroed rather than deleted.** A section a potential is expected to carry and
does not is a different test -- possibly a failing one, since GAMESS may require
it -- and that is not what is being asked here. Keeping the section present with
zero values changes the energy and nothing else.
"""

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from efp_format import parse_efp, render_efp  # noqa: E402

#: Sections whose values are numbers this script knows how to zero. Anything else
#: is refused rather than silently passed through, because a typo in a section
#: name would otherwise look exactly like a term that made no difference.
ZEROABLE = (
    "MONOPOLES", "DIPOLES", "QUADRUPOLES", "OCTUPOLES",
    "SCREEN", "SCREEN2",
    "POLARIZABLE POINTS", "DYNAMIC POLARIZABLE POINTS",
)


def zero_records(doc, section):
    """Every numeric value in one section set to zero, its labels untouched."""
    body = doc["sections"].get(section)
    if body is None:
        raise SystemExit(f"the potential has no {section} section")
    count = 0
    for record in body:
        if isinstance(record, dict) and "values" in record:
            record["values"] = [0.0 for _ in record["values"]]
            count += 1
        elif isinstance(record, dict) and "tensor" in record:
            record["tensor"] = [[0.0 for _ in row] for row in record["tensor"]]
            count += 1
    return count


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("potential", type=pathlib.Path)
    ap.add_argument("--zero", nargs="*", default=[],
                    help=f"sections to zero; one or more of {', '.join(ZEROABLE)}")
    ap.add_argument("--delete", nargs="*", default=[],
                    help="sections to remove entirely. Zeroing is wrong for the "
                         "screening sections: their value is a damping exponent, "
                         "and zero means the strongest possible damping rather "
                         "than none, so switching penetration off means removing "
                         "the section.")
    ap.add_argument("-o", "--output", type=pathlib.Path, default=None,
                    help="write here rather than to stdout")
    args = ap.parse_args()

    unknown = [s for s in args.zero if s not in ZEROABLE]
    if unknown:
        raise SystemExit(f"cannot zero {unknown}; known: {list(ZEROABLE)}")
    if not args.zero and not args.delete:
        raise SystemExit("nothing to do: pass --zero and/or --delete")

    doc = parse_efp(args.potential.read_text())
    for section in args.delete:
        if section not in doc["sections"]:
            raise SystemExit(f"the potential has no {section} section")
        del doc["sections"][section]
        doc.get("stop_after", set()).discard(section) if isinstance(
            doc.get("stop_after"), set) else None
        for key in ("header_lines", "stop_lines"):
            if isinstance(doc.get(key), dict):
                doc[key].pop(section, None)
        print(f"  removed {section}", file=sys.stderr)
    for section in args.zero:
        n = zero_records(doc, section)
        print(f"  zeroed {n} records in {section}", file=sys.stderr)
    text = render_efp(doc)
    if args.output:
        args.output.write_text(text)
        print(f"  wrote {args.output}", file=sys.stderr)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
