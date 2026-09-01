#!/usr/bin/env bash
# Launch mqc on a Perlmutter GPU node (64 physical / 128 logical cores).
set -euo pipefail

BINARY=./build/mqc
NTASKS=1
NTHREADS=64
CPUS_PER_TASK=
DRY_RUN=0
INPUT=

usage() {
    cat <<'USAGE'
Usage: ./run.sh [options] <input.json>

Options:
  -nt, --nthreads N       OpenMP threads per rank        (default: 64)
  -n,  --ntasks N         MPI ranks                      (default: 1)
  -c,  --cpus-per-task N  srun -c value                  (default: logical cores / ntasks)
  -b,  --binary PATH      mqc executable                 (default: ./build/mqc)
       --dry-run          print the command, do not run it
  -h,  --help             this message

Examples:
  ./run.sh fukui.json
  ./run.sh -nt 32 validation/inputs/cpu/mqc/rhf/h2o.json
  ./run.sh -n 4 -nt 16 prism.json

mqc takes the deck and nothing else -- a second argument is a hard error.
Verbosity is "system": {"logger": {"level": "verbose"}} in the deck itself;
levels are debug, verbose, info (default), performance, warning, error.
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -nt|--nthreads)      NTHREADS="$2"; shift 2 ;;
        -n|--ntasks)         NTASKS="$2"; shift 2 ;;
        -c|--cpus-per-task)  CPUS_PER_TASK="$2"; shift 2 ;;
        -b|--binary)         BINARY="$2"; shift 2 ;;
        --dry-run)           DRY_RUN=1; shift ;;
        -h|--help)           usage; exit 0 ;;
        -*)                  echo "run.sh: unknown option '$1'" >&2; usage >&2; exit 2 ;;
        *)
            if [[ -n "$INPUT" ]]; then
                echo "run.sh: only one input file may be given (got '$INPUT' and '$1')" >&2
                exit 2
            fi
            INPUT="$1"; shift ;;
    esac
done

if [[ -z "$INPUT" ]]; then
    echo "run.sh: no input file given" >&2
    usage >&2
    exit 2
fi
[[ -f "$INPUT" ]]   || { echo "run.sh: input '$INPUT' not found" >&2; exit 1; }
[[ -x "$BINARY" ]]  || { echo "run.sh: '$BINARY' is missing or not executable -- build first?" >&2; exit 1; }

for pair in "threads:$NTHREADS" "ranks:$NTASKS"; do
    [[ "${pair#*:}" =~ ^[1-9][0-9]*$ ]] || { echo "run.sh: ${pair%%:*} must be a positive integer" >&2; exit 2; }
done

# Logical CPUs available to the job; the node has 128 when SLURM does not say.
TOTAL_CPUS=${SLURM_CPUS_ON_NODE:-128}
if [[ -z "$CPUS_PER_TASK" ]]; then
    CPUS_PER_TASK=$(( TOTAL_CPUS / NTASKS ))
    (( CPUS_PER_TASK > 0 )) || CPUS_PER_TASK=1
fi

export OMP_NUM_THREADS="$NTHREADS"
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

declare -a CMD
if command -v srun >/dev/null 2>&1; then
    CMD=(srun -n "$NTASKS" -c "$CPUS_PER_TASK" --cpu-bind=cores "$BINARY" "$INPUT")
elif (( NTASKS > 1 )); then
    CMD=(mpirun -np "$NTASKS" "$BINARY" "$INPUT")
else
    CMD=("$BINARY" "$INPUT")
fi
echo "# ${NTASKS} rank(s) x ${NTHREADS} thread(s), -c ${CPUS_PER_TASK}"
echo "# ${CMD[*]}"
(( DRY_RUN )) && exit 0

exec "${CMD[@]}"
