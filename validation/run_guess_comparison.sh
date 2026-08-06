#!/bin/bash
# Compare SCF initial guesses: core, GWH and SAC.
#
#   ./validation/run_guess_comparison.sh [path/to/mqc]
#
# Each system is run with all three guesses. What to look for:
#
#   * The ENERGY should be identical across guesses. A guess changes the path
#     to convergence, not the answer -- unless it changes which solution you
#     land on, which is exactly what happens to OH below.
#   * The ITERATION count should fall from core to GWH to SAC. That is the
#     whole point of a better guess.
#
# The OH row is the cautionary one: with a core guess it converges cleanly,
# with a respectable <S^2>, onto the A2Sigma+ excited state about 4.3 eV above
# the ground state. Nothing about the convergence gives it away.
set -u
cd "$(dirname "$0")/.."
MQC=${1:-./build/mqc}

if [ ! -x "$MQC" ]; then
    echo "No mqc binary at $MQC -- build with -DMQC_ENABLE_CUEST=ON first." >&2
    exit 1
fi

SYSTEMS="hf_water_def2-svp dft_water_pbe uhf_oh uhf_o2"
GUESSES="core gwh sac"
# Temporaries live beside the originals: the `xyz` path inside each JSON is
# relative to the JSON's own directory, so moving them elsewhere breaks it.
WORK=validation/inputs
PREFIX=_guess_tmp_
trap 'rm -f "$WORK/${PREFIX}"*' EXIT

printf '%-22s %-6s %22s %6s %10s\n' SYSTEM GUESS "TOTAL ENERGY (Ha)" ITERS "<S^2>"
echo "---------------------------------------------------------------------------"

for sys in $SYSTEMS; do
    first_energy=""
    for guess in $GUESSES; do
        # Rewrite the input with this guess, keeping everything else identical.
        stem="${PREFIX}${sys}_${guess}"
        python3 - "$sys" "$guess" "$WORK/$stem.json" <<'PY'
import json, sys
name, guess, dest = sys.argv[1], sys.argv[2], sys.argv[3]
d = json.load(open(f"validation/inputs/{name}.json"))
d.setdefault("keywords", {}).setdefault("scf", {"maxiter": 100, "tolerance": 1e-8})["guess"] = guess
json.dump(d, open(dest, "w"), indent=4)
PY
        python3 mqc_prep.py "$WORK/$stem.json" > /dev/null 2>&1 || { echo "prep failed"; continue; }

        out=$($MQC "$WORK/$stem.mqc" 2>&1)
        energy=$(printf '%s' "$out" | grep "Final energy:" | awk '{print $NF}')
        # Iteration lines are "<n> <energy> <dE> <diis>"; the last one wins.
        iters=$(printf '%s' "$out" | awk '/^ *[0-9]+ +-?[0-9]+\.[0-9]+ +-?[0-9]/ {n=$1} END {print n+0}')
        s2=$(printf '%s' "$out" | grep "<S\^2>" | awk '{print $3}')
        [ -z "$s2" ] && s2="-"

        if [ -z "$energy" ]; then
            energy=$(printf '%s' "$out" | grep -i "ERROR:" | head -1 | sed 's/.*: //' | cut -c1-22)
            iters="-"
        fi
        printf '%-22s %-6s %22s %6s %10s\n' "$sys" "$guess" "$energy" "$iters" "$s2"

        # Flag a guess that changed the answer rather than just the path.
        if [ -n "$energy" ] && [ "$iters" != "-" ]; then
            if [ -z "$first_energy" ]; then
                first_energy=$energy
            else
                python3 -c "
import sys
a,b='$first_energy','$energy'
try:
    d=abs(float(a)-float(b))
    if d > 1e-6: print(f'    ^^ differs from the first guess by {d:.4f} Ha = {d*27.2114:.2f} eV')
except ValueError: pass
"
            fi
        fi
    done
    echo
done
