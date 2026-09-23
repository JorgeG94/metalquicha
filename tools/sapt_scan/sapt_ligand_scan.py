#!/usr/bin/env python3
"""Decompose a ligand's interaction with each nearby residue, by SAPT.

One SAPT calculation per residue-ligand pair, collected into one table:
residue, distance, and the named SAPT terms, sorted by total interaction.

Why this works where the fragmented methods do not. EFMO has the named terms
and cannot cut a protein at all -- it has no bond-breaking option. FMO can cut
one but its pair energy is a single number. SAPT decomposes the interaction of
*exactly two* monomers, and a ligand is not covalently bonded to anything, so
every residue-ligand pair is a legitimate SAPT problem as the program stands.

The input is an ordinary fragmented deck -- the same one an MBE run would take,
with ``fragments`` and ``connectivity`` already declared -- plus which fragment
is the ligand. Nothing about the partition is restated here.

**The capping is done by this script, and that is not a detail.** The SAPT
route in ``mqc_driver.f90`` slices its two monomers straight out of the
geometry: it never calls ``build_fragment_from_indices``, so it adds no
hydrogen caps and reads no ``connectivity``. Declaring bonds in a SAPT deck
does nothing. A residue handed over uncapped is a radical with a dangling
valence, which either dies on an electron-count parity error or, worse,
converges to a wrong number. So the caps are placed here, by the same rule the
Fortran uses, and written into the geometry as real atoms.

Run with ``--help`` for the options. See ``mqc_docs/source/sapt_ligand_scan.rst``
for what the resulting numbers do and do not include.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

# `R_H = R_kept + s (R_gone - R_kept)`, and s = 1 puts the cap exactly on the
# atom it replaces. That is `place_caps` in `mqc_physical_fragment.f90` and it
# is the default there, so it is the default here: a scan that silently used a
# different rule from the rest of the program would not be comparable with it.
DEFAULT_CAP_SCALE = 1.0

BOHR_PER_ANGSTROM = 1.8897261254578281

# Slot order of `SAPT_TERM_NAMES` / `SAPT2_TERM_NAMES` in `mqc_program_limits`.
# Named here so the table's columns are the literature's names rather than
# whatever order a JSON object happened to serialise in.
SAPT0_TERMS = [
    "elst10", "exch10_s2", "exch10", "ind20_u", "ind20_r",
    "exch_ind20_u", "exch_ind20_r", "disp20", "exch_disp20",
    "delta_hf", "e_int_hf_cp", "total",
]
SAPT2_EXTRA = ["elst12", "exch11", "exch12", "ind22", "exch_ind22", "total_sapt2"]

HARTREE_TO_KCAL = 627.509474


class ScanError(Exception):
    """Something the caller can act on, as opposed to a traceback."""


def read_xyz(path):
    """Elements and Angstrom coordinates from a plain XYZ file."""
    with open(path) as handle:
        lines = handle.read().splitlines()
    n = int(lines[0].split()[0])
    elements, coords = [], []
    for line in lines[2:2 + n]:
        parts = line.split()
        elements.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return elements, coords


def write_xyz(path, elements, coords, comment=""):
    with open(path, "w") as handle:
        handle.write(f"{len(elements)}\n{comment}\n")
        for element, xyz in zip(elements, coords):
            handle.write("%-3s %18.10f %18.10f %18.10f\n" % (element, *xyz))


def load_deck(path):
    """The geometry, fragments and bonds of an ordinary fragmented deck.

    Coordinates come back in Angstrom. A deck may carry its geometry inline as
    ``symbols`` plus a flat ``geometry`` list, or point at an ``xyz`` file
    relative to itself; both are read here because both are shipped.
    """
    with open(path) as handle:
        deck = json.load(handle)
    molecule = deck["molecules"][0]

    if "xyz" in molecule:
        xyz_path = os.path.join(os.path.dirname(os.path.abspath(path)), molecule["xyz"])
        elements, coords = read_xyz(xyz_path)
    elif "symbols" in molecule:
        elements = list(molecule["symbols"])
        flat = molecule["geometry"]
        coords = [flat[i:i + 3] for i in range(0, len(flat), 3)]
    else:
        raise ScanError(f"{path}: the molecule declares neither 'xyz' nor 'symbols'")

    fragments = molecule.get("fragments")
    if not fragments:
        raise ScanError(
            f"{path}: no 'fragments'. This scan needs the residue partition the "
            "deck would use for a fragmented run; it does not invent one."
        )
    bonds = [(b[0], b[1]) for b in molecule.get("connectivity", [])]
    charges = molecule.get("fragment_charges", [0] * len(fragments))
    return elements, coords, [list(f) for f in fragments], bonds, charges


def min_distance(coords, atoms_a, atoms_b):
    """Closest approach between two atom sets, in Angstrom.

    The minimum, matching what the fragmented path screens on -- a centroid
    distance would call a long residue far away while one end of it sits on
    the ligand.
    """
    best = float("inf")
    for i in atoms_a:
        xi, yi, zi = coords[i]
        for j in atoms_b:
            xj, yj, zj = coords[j]
            d = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2)
            if d < best:
                best = d
    return best


def capped_residue(elements, coords, residue, bonds, cap_scale):
    """A residue's atoms plus a hydrogen for every bond leaving it.

    `R_H = R_kept + s (R_gone - R_kept)`, the rule `place_caps` uses, applied
    to each bond with exactly one end inside the residue. Returns elements and
    Angstrom coordinates, the caps appended after the residue's own atoms.
    """
    members = set(residue)
    out_elements = [elements[i] for i in residue]
    out_coords = [list(coords[i]) for i in residue]
    n_caps = 0
    for atom_i, atom_j in bonds:
        in_i, in_j = atom_i in members, atom_j in members
        if in_i == in_j:
            continue
        kept, gone = (atom_i, atom_j) if in_i else (atom_j, atom_i)
        rk, rg = coords[kept], coords[gone]
        out_elements.append("H")
        out_coords.append([rk[k] + cap_scale * (rg[k] - rk[k]) for k in range(3)])
        n_caps += 1
    return out_elements, out_coords, n_caps


def sapt_deck(xyz_name, n_a, n_b, basis, method, charge_a, charge_b, scf_tolerance):
    """A two-fragment SAPT deck over a geometry of residue-then-ligand.

    No ``connectivity``: the SAPT route never reads it, and writing one would
    suggest this deck gets capped when the caps are already in the geometry.
    """
    return {
        "schema": {"name": "mqc-frag", "version": "1.0"},
        "molecules": [{
            "xyz": xyz_name,
            "molecular_charge": charge_a + charge_b,
            "molecular_multiplicity": 1,
            "fragments": [list(range(n_a)), list(range(n_a, n_a + n_b))],
        }],
        "model": {"method": method, "basis": basis},
        "keywords": {"scf": {"tolerance": scf_tolerance}},
        "driver": "Energy",
    }


def run_one(job):
    """One residue-ligand SAPT calculation, in its own directory.

    Failure is detected by reading the output, never by exit status: a refusal
    in this program exits zero, so an exit code says nothing about whether a
    number came back.
    """
    label = job["label"]
    workdir = os.path.join(job["scratch"], label)
    os.makedirs(workdir, exist_ok=True)

    write_xyz(os.path.join(workdir, "pair.xyz"), job["elements"], job["coords"],
              comment=f"{label}: capped residue then ligand")
    deck_path = os.path.join(workdir, f"{label}.json")
    with open(deck_path, "w") as handle:
        json.dump(job["deck"], handle, indent=1)

    completed = subprocess.run(
        [job["exe"], f"{label}.json"], cwd=workdir,
        capture_output=True, text=True, timeout=job["timeout"],
    )
    log = completed.stdout + completed.stderr
    with open(os.path.join(workdir, "run.log"), "w") as handle:
        handle.write(log)

    out_path = os.path.join(workdir, f"output_{label}.json")
    if not os.path.exists(out_path):
        return {"label": label, "error": _first_refusal(log) or "no output file was written"}
    with open(out_path) as handle:
        document = json.load(handle)
    body = document.get(label)
    if body is None or "sapt" not in body:
        return {"label": label, "error": _first_refusal(log) or "the output carries no sapt section"}

    terms = body["sapt"]
    bad = [k for k, v in terms.items()
           if isinstance(v, float) and not math.isfinite(v)]
    if bad:
        return {"label": label, "error": f"non-finite SAPT terms: {', '.join(sorted(bad))}"}
    return {"label": label, "terms": terms}


def _first_refusal(log):
    """The line a refusal printed, if one did.

    Refusals name themselves and then exit zero, so this is the only signal
    that distinguishes 'declined' from 'crashed'.
    """
    for line in log.splitlines():
        stripped = line.strip()
        if re.match(r"^(SAPT|Error|ERROR|.*refus)", stripped) and len(stripped) > 12:
            return stripped[:300]
    return None


def scan(args):
    elements, coords, fragments, bonds, charges = load_deck(args.deck)
    if not 1 <= args.ligand <= len(fragments):
        raise ScanError(
            f"--ligand {args.ligand} is out of range; the deck has "
            f"{len(fragments)} fragments, numbered from one"
        )
    ligand = fragments[args.ligand - 1]
    ligand_charge = charges[args.ligand - 1] if args.ligand - 1 < len(charges) else 0

    # A ligand covalently bonded to the protein is not a ligand for this
    # purpose, and its pair would be the one case SAPT cannot answer.
    ligand_set = set(ligand)
    crossing = [(i, j) for i, j in bonds
                if (i in ligand_set) != (j in ligand_set)]
    if crossing:
        raise ScanError(
            f"fragment {args.ligand} has {len(crossing)} covalent bond(s) to the rest "
            "of the system, so it is not an unbonded ligand. SAPT between two "
            "halves of a molecule is not defined; this scan refuses rather than "
            "capping the ligand and reporting a number for it."
        )

    jobs = []
    for index, residue in enumerate(fragments, start=1):
        if index == args.ligand:
            continue
        distance = min_distance(coords, residue, ligand)
        if distance > args.cutoff:
            continue
        res_elements, res_coords, n_caps = capped_residue(
            elements, coords, residue, bonds, args.cap_scale)
        pair_elements = res_elements + [elements[i] for i in ligand]
        pair_coords = res_coords + [list(coords[i]) for i in ligand]
        label = f"res{index:03d}_ligand"
        residue_charge = charges[index - 1] if index - 1 < len(charges) else 0
        jobs.append({
            "label": label, "index": index, "distance": distance, "n_caps": n_caps,
            "elements": pair_elements, "coords": pair_coords,
            "deck": sapt_deck("pair.xyz", len(res_elements), len(ligand),
                              args.basis, args.method, residue_charge,
                              ligand_charge, args.scf_tolerance),
            "exe": args.exe, "scratch": args.scratch, "timeout": args.timeout,
        })

    if not jobs:
        raise ScanError(
            f"no fragment lies within {args.cutoff} A of fragment {args.ligand}. "
            "Nothing to do -- widen --cutoff or check the ligand number."
        )

    print(f"{len(jobs)} residue-ligand pairs within {args.cutoff} A, "
          f"{args.method} / {args.basis}, cap_scale {args.cap_scale}", file=sys.stderr)

    results = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        for outcome in pool.map(run_one, jobs):
            results[outcome["label"]] = outcome
            mark = "ok" if "terms" in outcome else "FAILED"
            print(f"  {outcome['label']}  {mark}", file=sys.stderr)

    rows = []
    for job in jobs:
        outcome = results[job["label"]]
        row = {"residue": job["index"], "distance": job["distance"],
               "n_caps": job["n_caps"]}
        if "terms" in outcome:
            row.update(outcome["terms"])
        else:
            row["error"] = outcome["error"]
        rows.append(row)
    return rows, args


def write_table(rows, args, handle):
    """The scan as CSV, strongest interaction first.

    The provenance lines are not decoration. Cap placement moves a
    hydrogen-bonded pair energy by about eight per cent, and nothing in the
    numbers themselves says which placement produced them, so two tables are
    not comparable unless both say.
    """
    term_names = SAPT0_TERMS + (SAPT2_EXTRA if args.method.lower() == "sapt2" else [])
    total_key = "total_sapt2" if args.method.lower() == "sapt2" else "total"

    handle.write(f"# method,{args.method}\n")
    handle.write(f"# basis,{args.basis}\n")
    handle.write(f"# cap_scale,{args.cap_scale}\n")
    handle.write(f"# cutoff_angstrom,{args.cutoff}\n")
    handle.write("# monomers are isolated: no protein environment is felt\n")
    handle.write("# energies in kcal/mol\n")

    done = [r for r in rows if total_key in r]
    failed = [r for r in rows if total_key not in r]
    done.sort(key=lambda r: r[total_key])

    handle.write("residue,distance,n_caps," + ",".join(term_names) + "\n")
    for row in done:
        cells = ["%d" % row["residue"], "%.3f" % row["distance"], "%d" % row["n_caps"]]
        cells += ["%.6f" % (row.get(name, float("nan")) * HARTREE_TO_KCAL)
                  for name in term_names]
        handle.write(",".join(cells) + "\n")
    for row in failed:
        handle.write("%d,%.3f,%d,%s\n" % (row["residue"], row["distance"],
                                          row["n_caps"], "ERROR: " + row["error"]))
    return len(done), len(failed)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="SAPT decomposition of a ligand's interaction with each nearby residue.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("deck", help="a fragmented deck: geometry, fragments, connectivity")
    parser.add_argument("--ligand", type=int, required=True,
                        help="which fragment is the ligand, numbered from one")
    parser.add_argument("--cutoff", type=float, default=6.0,
                        help="include residues within this many Angstrom of the ligand")
    parser.add_argument("--basis", default="sto-3g")
    parser.add_argument("--method", default="sapt0", choices=["sapt0", "sapt2"])
    parser.add_argument("--cap-scale", type=float, default=DEFAULT_CAP_SCALE,
                        dest="cap_scale",
                        help="R_H = R_kept + s (R_gone - R_kept); 1.0 puts the cap on "
                             "the atom it replaces, as the rest of the program does")
    parser.add_argument("--scf-tolerance", type=float, default=1e-10, dest="scf_tolerance")
    parser.add_argument("--exe", default="./build/mqc")
    parser.add_argument("--workers", type=int, default=1,
                        help="pairs to run at once; each is one small serial job")
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument("--scratch", default=None,
                        help="where the per-pair decks and logs go; a temporary "
                             "directory that is kept, if not given")
    parser.add_argument("-o", "--output", default=None, help="CSV out; stdout if absent")
    args = parser.parse_args(argv)

    args.exe = os.path.abspath(args.exe)
    if not os.path.exists(args.exe):
        raise ScanError(f"no executable at {args.exe}; build it or pass --exe")
    if args.scratch is None:
        args.scratch = tempfile.mkdtemp(prefix="sapt_scan_")
        print(f"scratch: {args.scratch}", file=sys.stderr)
    os.makedirs(args.scratch, exist_ok=True)

    rows, args = scan(args)
    if args.output:
        with open(args.output, "w") as handle:
            done, failed = write_table(rows, args, handle)
        print(f"wrote {args.output}", file=sys.stderr)
    else:
        done, failed = write_table(rows, args, sys.stdout)

    print(f"{done} pairs decomposed, {failed} failed", file=sys.stderr)
    # A failed pair is a missing row, not a smaller table: exit non-zero so a
    # caller that does check exit status is told, even though the program this
    # drives does not.
    return 1 if failed else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ScanError as exc:
        print(f"sapt_ligand_scan: {exc}", file=sys.stderr)
        sys.exit(2)
