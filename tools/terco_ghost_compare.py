"""Row-by-row terco vs libcint for the VMFC water-dimer B3LYP decks.

    python tools/terco_ghost_compare.py output_cp_b3lyp_libcint_fragments.csv output_cp_b3lyp_terco_fragments.csv

A bare monomer row ([1] or [2]) has no ghosts, so a difference there is plain
B3LYP disagreeing between the two codes (grid or functional definition). A
ghosted row ([1] ghosts [2]) adds the ghost; if it differs by much more than
the bare rows, the ghost handling is what disagrees.
"""
import csv
import sys


def rows(path):
    out = {}
    with open(path) as handle:
        for r in csv.DictReader(handle):
            mons = [int(r[k]) for k in r if k and k[0] == "m" and k[1:].isdigit() and r[k].strip()]
            real = tuple(m for m in mons if m > 0)
            ghosts = tuple(-m for m in mons if m < 0)
            out[(real, ghosts)] = float(r["energy"])
    return out


cpu, gpu = rows(sys.argv[1]), rows(sys.argv[2])
print(f"{'row':<22}{'libcint':>20}{'terco':>20}{'terco-cpu':>12}")
for key in sorted(cpu, key=lambda k: (len(k[0]), len(k[1]), k)):
    real, ghosts = key
    label = str(list(real)) + (f" ghosts {list(ghosts)}" if ghosts else "")
    e_c, e_g = cpu[key], gpu.get(key)
    shown = f"{e_g:20.12f}" if e_g is not None else f"{'missing':>20}"
    diff = f"{e_g - e_c:+12.2e}" if e_g is not None else ""
    print(f"{label:<22}{e_c:20.12f}{shown}{diff}")
for (real, ghosts), e_c in cpu.items():
    if ghosts and (real, ()) in cpu and (real, ()) in gpu and (real, ghosts) in gpu:
        b_c = e_c - cpu[(real, ())]
        b_g = gpu[(real, ghosts)] - gpu[(real, ())]
        print(f"ghost lowering of {list(real)}: libcint {b_c:+.6e}  terco {b_g:+.6e}  "
              f"diff {b_g - b_c:+.2e}")
