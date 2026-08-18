#!/bin/bash
# Check that spreading one system over several GPUs does not change its answer.
#
#   ./validation/run_multigpu_test.sh [n_ranks] [path/to/mqc]
#
# The claim the multi-GPU mode makes is exact, not approximate: the radial
# quadrature is partitioned across ranks and the pieces are summed, and a Becke
# weight does not depend on which other grid points exist, so the reduced
# integral is the single-GPU integral term for term. This script is what turns
# that claim into a number.
#
# One deck is run twice. At -np 1 the team is dormant -- `gpu_team_active` is
# false whatever the deck asked for -- so that run IS the single-GPU reference,
# built by the same code path rather than by a second deck that might drift
# away from it. At -np N the grid is sliced N ways.
#
# What to expect: agreement to the last printed digit. This is not a tolerance
# test in spirit. The two runs differ only in the ORDER an identical set of
# quadrature contributions is summed, which in floating point is worth an ulp
# or two on the total energy -- so a disagreement in the sixth decimal is not
# "close enough", it means the partition dropped or double-counted points.
#
# Needs a GPU of compute capability 8.0 or newer per rank, and at least as many
# visible GPUs as ranks -- ranks past the device count share a device, which
# still gives the right answer but measures nothing about memory.
set -u
cd "$(dirname "$0")/.."
NP=${1:-2}
MQC=${2:-./build/mqc}

if [ ! -x "$MQC" ]; then
    echo "No mqc binary at $MQC -- build with -DMQC_ENABLE_CUEST=ON first." >&2
    exit 1
fi

find_deck() { find validation/inputs -name "$1.json" -print -quit; }

# The energy and the gradient norm both, deliberately. The energy alone would
# miss a partition that is wrong only in the XC derivative, and the gradient is
# the calculation this mode exists for.
extract() {
    printf '%s' "$1" | awk '
        /Final energy:/      { e = $NF }
        /Gradient norm/      { g = $NF }
        END { printf "%s %s", (e == "" ? "NONE" : e), (g == "" ? "NONE" : g) }'
}

run_case() {   # $1 = label, $2 = ranks, $3 = deck stem
    local deck out
    deck=$(find_deck "$3")
    if [ -z "$deck" ]; then
        printf '%-28s %s\n' "$1" "deck '$3' not found"
        return 1
    fi
    if [ "$2" -eq 1 ]; then
        out=$("$MQC" "$deck" 2>&1)
    else
        out=$(mpirun -np "$2" "$MQC" "$deck" 2>&1)
    fi
    extract "$out"
}

echo "Multi-GPU consistency: one system, 1 rank vs $NP ranks"
echo "-------------------------------------------------------------------"

status=0
for stem in multigpu_gly3_b3lyp_grad multigpu_gly3_b3lyp_grad_direct; do
    one=$(run_case "$stem (1 rank)" 1 "$stem") || { status=1; continue; }
    many=$(run_case "$stem ($NP ranks)" "$NP" "$stem") || { status=1; continue; }

    printf '%-38s\n' "$stem"
    printf '  %-10s energy=%-22s |grad|=%s\n' "1 rank"  ${one}
    printf '  %-10s energy=%-22s |grad|=%s\n' "$NP ranks" ${many}

    if [ "$one" = "$many" ]; then
        echo "  MATCH"
    else
        echo "  DIFFER -- the radial partition is not summing back to the whole grid."
        status=1
    fi
    echo
done

# The second deck adds df_integral_direct and a recompute gradient policy. It
# must reach the same numbers as the first: those keywords trade time for
# device memory and are not allowed to change the answer. If deck one matches
# across ranks and deck two does not, suspect the memory policy rather than the
# grid.
echo "Both decks should also agree with EACH OTHER: integral-direct and the"
echo "gradient memory policy trade time for memory and must not move the energy."

exit $status
