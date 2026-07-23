#!/bin/bash
# Validate the cuEST-backed SCF (HF and DFT) against reference energies for water.
#
#   ./validation/run_cuest_validation.sh [path/to/mqc]
#
# Requires a GPU of compute capability 8.0 or newer. On anything older -- the
# Gadi login node's V100 -- every row reports UNSUPPORTED_ARCHITECTURE, which
# is expected rather than a failure of this code.
#
# Geometry: r(OH) = 0.959 A, angle 107.3 deg. Every row fits J and K with
# def2-universal-JKFIT, so each carries a small (~1e-4 Ha) RI error relative to
# the un-fitted reference.
#
# On the references: def2-SVP sits ~0.066 Ha ABOVE cc-pVDZ for the same method.
# Both contract to [3s2p1d], but oxygen's primitive sets differ -- (7s4p1d) vs
# (9s4p1d) -- and total energy is dominated by the 1s core, so the two missing
# tight s functions cost real energy while valence properties stay comparable.
# The def2-SVP references below already include that offset.
set -u
cd "$(dirname "$0")/.."
MQC=${1:-./build/mqc}

if [ ! -x "$MQC" ]; then
    echo "No mqc binary at $MQC -- build with -DMQC_ENABLE_CUEST=ON first." >&2
    exit 1
fi

INPUTS="hf_water_sto-3g hf_water_cc-pvdz hf_water_def2-svp
        dft_water_svwn5 dft_water_pbe dft_water_b3lyp dft_water_pbe0"
for name in $INPUTS; do
    python3 mqc_prep.py "validation/inputs/$name.json" > /dev/null || exit 1
done

run_one() {   # $1 = label, $2 = input stem, $3 = reference
    local out energy
    out=$($MQC "validation/inputs/$2.mqc" 2>&1)
    energy=$(printf '%s' "$out" | grep "Final energy:" | awk '{print $NF}')
    if [ -z "$energy" ]; then
        energy=$(printf '%s' "$out" | grep -i "ERROR:" | head -1 | sed 's/.*: //' | cut -c1-40)
    fi
    printf '%-22s %24s %12s\n' "$1" "$energy" "$3"
}

printf '%-22s %24s %12s\n' "CALCULATION" "TOTAL ENERGY (Ha)" "REFERENCE"
echo "--- Hartree-Fock: basis-set convergence ------------------------------"
run_one "HF/STO-3G"      hf_water_sto-3g   "-74.962"
run_one "HF/cc-pVDZ"     hf_water_cc-pvdz  "-76.027"
run_one "HF/def2-SVP"    hf_water_def2-svp "-75.961"
echo "--- DFT: functional ladder, def2-SVP ---------------------------------"
run_one "SVWN5/def2-SVP"  dft_water_svwn5  "-75.795"
run_one "PBE/def2-SVP"    dft_water_pbe    "-76.272"
run_one "B3LYP/def2-SVP"  dft_water_b3lyp  "-76.321"
run_one "PBE0/def2-SVP"   dft_water_pbe0   "-76.276"
