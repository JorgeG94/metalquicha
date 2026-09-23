#!/usr/bin/env python3
"""Check a SAPT pair against a counterpoise-corrected supermolecular energy.

If the terms of a decomposition are going to be quoted, the total they sum to
needs a measured error rather than an assumed one. This measures it.

The comparison is against `e_int_hf_cp`, the counterpoise-corrected
Hartree-Fock interaction energy SAPT reports in its own table. That is the
term with an independent definition: it is what a supermolecular calculation
of the same pair, in the dimer basis, computes directly. The VMFC
counterpoise path over the same two fragments produces exactly that quantity
as its two-body correction, so the two numbers are the same physical thing
arrived at by different routes.

`total` is deliberately **not** the thing compared. It carries dispersion,
which Hartree-Fock does not have, so it cannot equal a supermolecular HF
interaction and a disagreement there would mean nothing.

Usage::

    python3 validate_against_counterpoise.py <pair-dir> --exe ./build/mqc

where the pair directory is one the scan wrote, holding `pair.xyz` and the
SAPT deck and output.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import subprocess
import sys

HARTREE_TO_KCAL = 627.509474


def counterpoise_interaction(pair_dir, exe, basis, tolerance, timeout):
    """The VMFC two-body correction for the pair, in Hartree.

    Under counterpoise every monomer is recomputed in the dimer basis, so the
    two-body correction is the basis-set-superposition-free interaction.
    """
    # Not "*_ligand.json": that also matches `output_..._ligand.json`, the
    # result of the run rather than its input.
    sapt_deck = [p for p in glob.glob(os.path.join(pair_dir, "*_ligand.json"))
                 if not os.path.basename(p).startswith("output_")]
    if not sapt_deck:
        raise SystemExit(f"{pair_dir}: no SAPT deck here; is this a pair directory?")
    with open(sapt_deck[0]) as handle:
        fragments = json.load(handle)["molecules"][0]["fragments"]
    n_a, n_b = len(fragments[0]), len(fragments[1])

    deck = {
        "schema": {"name": "mqc-frag", "version": "1.0"},
        "molecules": [{
            "xyz": "pair.xyz",
            "molecular_charge": 0, "molecular_multiplicity": 1,
            "fragments": [list(range(n_a)), list(range(n_a, n_a + n_b))],
            "fragment_charges": [0, 0], "fragment_multiplicities": [1, 1],
        }],
        "model": {"method": "HF", "basis": basis},
        "driver": "Energy",
        "keywords": {
            "fragmentation": {"method": "mbe", "level": 2, "embedding": "none",
                              "counterpoise": "vmfc"},
            "scf": {"tolerance": tolerance},
        },
    }
    with open(os.path.join(pair_dir, "cp_check.json"), "w") as handle:
        json.dump(deck, handle, indent=1)

    subprocess.run([os.path.abspath(exe), "cp_check.json"], cwd=pair_dir,
                   capture_output=True, text=True, timeout=timeout)

    table = os.path.join(pair_dir, "output_cp_check_fragments.csv")
    if not os.path.exists(table):
        raise SystemExit(f"{pair_dir}: the counterpoise run wrote no table. "
                         "Read cp_check's log -- a refusal exits zero.")
    with open(table) as handle:
        rows = [r for r in csv.DictReader(
            line for line in handle if not line.startswith("#"))]
    return sum(float(r["delta_energy"]) for r in rows if int(r["level"]) == 2)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("pair_dir")
    parser.add_argument("--exe", default="./build/mqc")
    parser.add_argument("--basis", default="sto-3g")
    parser.add_argument("--scf-tolerance", type=float, default=1e-10)
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument("--threshold", type=float, default=1e-8,
                        help="Hartree; the two routes should agree to SCF convergence")
    args = parser.parse_args(argv)

    out = glob.glob(os.path.join(args.pair_dir, "output_*_ligand.json"))
    if not out:
        raise SystemExit(f"{args.pair_dir}: no SAPT output here")
    with open(out[0]) as handle:
        document = json.load(handle)
    sapt = document[list(document)[0]]["sapt"]

    cp = counterpoise_interaction(args.pair_dir, args.exe, args.basis,
                                  args.scf_tolerance, args.timeout)
    ref = sapt["e_int_hf_cp"]
    difference = abs(cp - ref)

    print("counterpoise-corrected supermolecular HF : %.12f Ha  (%+8.4f kcal/mol)"
          % (cp, cp * HARTREE_TO_KCAL))
    print("SAPT e_int_hf_cp                         : %.12f Ha  (%+8.4f kcal/mol)"
          % (ref, ref * HARTREE_TO_KCAL))
    print("difference                               : %.2e Ha  (%.2e kcal/mol)"
          % (difference, difference * HARTREE_TO_KCAL))
    print("SAPT0 total (carries dispersion, not comparable to HF) : %+8.4f kcal/mol"
          % (sapt["total"] * HARTREE_TO_KCAL))

    if difference > args.threshold:
        print("FAILED: above %.1e Ha" % args.threshold, file=sys.stderr)
        return 1
    print("agrees to SCF convergence")
    return 0


if __name__ == "__main__":
    sys.exit(main())
