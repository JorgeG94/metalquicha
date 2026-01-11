#!/usr/bin/env python3
"""
Strict MQC JSON -> v1 Fortran-input emitter.

This script converts JSON input files to the MQC v1 format.
Output file is automatically named based on input JSON or 'title' field.
"""

from __future__ import annotations
import argparse
from pathlib import Path

from utils import parse_json_to_input, emit_v1, die


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Strict MQC JSON -> v1 Fortran-input emitter",
        epilog="Output file is automatically named based on input JSON or 'title' field."
    )
    parser.add_argument("json_file", help="Input JSON file")
    args = parser.parse_args()

    json_path = Path(args.json_file)
    if not json_path.exists():
        die(f"Input file not found: {json_path}")

    inp = parse_json_to_input(json_path)
    text, out_path = emit_v1(inp, json_path)

    print(f"Generated: {out_path}")
    print(f"  Lines: {len(text.splitlines())}")
    print(f"  Size: {len(text)} bytes")


if __name__ == "__main__":
    main()
