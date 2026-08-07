#!/bin/bash
# Fragmented (MBE) test of the cuEST backend: water hexamer, MBE(2), one rank.
#
#   ./validation/run_fragment_test.sh
#
# 6 monomers + 15 dimers = 21 subsystems queued through a single rank. Two
# things get exercised that a single-molecule run cannot:
#
#   * the cuEST handle is created once and reused for all 21 fragments --
#     re-creating it per fragment would dominate the runtime;
#   * monomers have 24 AO functions and dimers 48, so the device scratch pool
#     has to GROW mid-run. A distributed-Hessian test could never reach that
#     branch, since every displacement is the same size.
set -u
cd "$(dirname "$0")/.."
MQC=./build/mqc

for tag in hf pbe; do
    echo "=== water hexamer MBE(2), $tag/def2-svp ==="
    # Anchored patterns, not substrings: a bare case-insensitive "MBE" also
    # matches "nuMBEr of atoms", which buries the result in per-fragment noise.
    /usr/bin/time -f "  wall: %e s   peak RSS: %M kB" \
        $MQC "validation/inputs/frag_water6_$tag.json" 2>&1 \
        | grep -E "Total fragments:|Computing Many-Body|^ *wall:|ERROR:|Total processing time"

    # The MBE totals live in the JSON, not the log.
    python3 -c "
import json
d = json.load(open('output_frag_water6_$tag.json'))['frag_water6_$tag']
one, tot = d['levels'][0]['total_energy'], d['total_energy']
print(f'  1-body           : {one:.6f} Ha')
print(f'  MBE(2) total     : {tot:.6f} Ha')
print(f'  2-body           : {tot-one:+.6f} Ha = {(tot-one)*627.5094740631:+.1f} kcal/mol')
"
    echo
done
