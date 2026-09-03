#!/bin/bash
# Compare SCF initial guesses on the CPU path: core, GWH, SAC and SAD.
#
#   ./validation/run_cpu_guess_comparison.sh [path/to/mqc]
#
# The libcint twin of run_guess_comparison.sh, which does the same for cuEST.
# What to look for:
#
#   * The ENERGY must be identical across guesses. A guess changes the path to
#     convergence, not the answer -- unless it changes which solution is found,
#     which is what happens to OH below.
#   * The ITERATION count is the point. cpu_peptide46 is the case that motivated
#     this: 46 atoms, charge +1, 134 functions, and a core guess does not
#     converge at all inside 100 cycles.
#
# Two results worth not being surprised by:
#
#   * SAC is usually *worse* than SAD on a closed-shell system, and can be worse
#     than GWH. It keeps the free atom's own orbitals, so an open-shell atom
#     hands over a density that singled out one p direction for reasons the
#     molecule knows nothing about, and the SCF spends iterations undoing it.
#     SAC earns its place on radicals, where that broken symmetry is the point.
#   * A better guess is not always fewer iterations by much. Most of what SAD
#     buys on well-behaved systems is the first two or three cycles; what it
#     really buys is the systems where the core guess does not converge.
set -u
cd "$(dirname "$0")/.."
MQC=${1:-./build/mqc}

# Decks live under inputs/<hardware>/<engine>/<method>/ now, so a deck name is
# resolved rather than assumed. Keeping the lookup here means the layout can move
# again without touching this script.
find_deck() { find validation/inputs -name "$1.json" -print -quit; }


if [ ! -x "$MQC" ]; then
    echo "No mqc binary at $MQC -- build with -DMQC_ENABLE_CZT=ON first." >&2
    exit 1
fi

# The peptide is the headline case. The density-fitted water is here because the
# fitted exchange is the one path that cannot take a density directly -- it goes
# through the pseudo-orbital factorisation instead -- so it is the only case that
# exercises it end to end. OH is the open-shell end, and the one where a guess can
# change the answer rather than the path.
SYSTEMS="cpu_peptide46_sto-3g cpu_water_cc-pvdz_df cpu_oh_def2-svp_m2"
GUESSES="core gwh sac sad"
PREFIX=_cpu_guess_tmp_
trap 'find validation/inputs -name "${PREFIX}*" -delete' EXIT

printf '%-24s %-6s %24s %7s\n' SYSTEM GUESS "TOTAL ENERGY (Ha)" ITERS
echo "-----------------------------------------------------------------------"

for sys in $SYSTEMS; do
    if [ -z "$(find_deck "$sys")" ]; then
        printf '%-24s %s\n' "$sys" "(no such input, skipped)"
        continue
    fi
    first_energy=""
    for guess in $GUESSES; do
        src=$(find_deck "$sys")
        tmp="$(dirname "$src")/${PREFIX}${sys}_${guess}.json"
        python3 - "$src" "$guess" "$tmp" <<'PY'
import json, sys
src, guess, dest = sys.argv[1], sys.argv[2], sys.argv[3]
d = json.load(open(src))
scf = d.setdefault("keywords", {}).setdefault("scf", {"maxiter": 100, "tolerance": 1e-8})
scf["guess"] = guess
d.setdefault("system", {}).setdefault("logger", {})["level"] = "Verbose"
json.dump(d, open(dest, "w"), indent=4)
PY

        out=$($MQC "$tmp" 2>&1)
        energy=$(printf '%s' "$out" | grep -i "Final energy" | awk '{print $NF}')
        iters=$(printf '%s' "$out" | grep -cE '^ *iter +[0-9]+')

        if [ -z "$energy" ]; then
            printf '%-24s %-6s %24s %7s\n' "$sys" "$guess" "did not converge" "$iters"
            continue
        fi
        printf '%-24s %-6s %24s %7s\n' "$sys" "$guess" "$energy" "$iters"

        # Flag a guess that changed the answer rather than just the path.
        if [ -z "$first_energy" ]; then
            first_energy=$energy
        else
            python3 -c "
a, b = '$first_energy', '$energy'
try:
    d = abs(float(a) - float(b))
    if d > 1e-6:
        print(f'    ^^ differs from the first converged guess by {d:.6f} Ha = {d*27.2114:.2f} eV')
except ValueError:
    pass
"
        fi
    done
    echo
done
