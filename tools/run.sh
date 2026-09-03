#!/usr/bin/env bash
# Launch mqc into whatever allocation you already hold: srun inside a SLURM job,
# mpirun outside one, or the binary directly for a single rank. Rank and thread
# placement is worked out from the node rather than left to whoever types the
# command.
set -euo pipefail

REPO_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)

BINARY=${MQC_BINARY:-}
NTASKS=1
NTHREADS=
CPUS_PER_TASK=
GPUS_PER_TASK=
DRY_RUN=0
INPUT=

die() { printf 'run.sh: %s\n' "$*" >&2; exit 2; }           # the command line is wrong
fail() { printf 'run.sh: %s\n' "$*" >&2; exit 1; }          # the machine is not ready
warn() { printf 'run.sh: warning: %s\n' "$*" >&2; }

usage() {
    cat <<'USAGE'
Usage: tools/run.sh [options] <input.json>

Options:
  -nt, --nthreads N       OpenMP threads per rank    (default: physical cores per rank)
  -n,  --ntasks N         MPI ranks                  (default: 1)
  -c,  --cpus-per-task N  srun -c value              (default: logical cores / ntasks)
  -g,  --gpus-per-task N  srun --gpus-per-task value (default: 1 in a GPU job, else none)
  -b,  --binary PATH      mqc executable             (default: $MQC_BINARY, else the
                                                      first of build/mqc, build/*/mqc)
       --dry-run          print the command and its environment, do not run it
       --version          ask mqc for its version
  -h,  --help             this message

Examples:
  tools/run.sh fukui.json
  tools/run.sh -nt 32 validation/inputs/cpu/mqc/rhf/h2o.json
  tools/run.sh -n 4 prism.json            # threads sized to the ranks, not repeated

This runs in the allocation you already hold; it does not submit one. On a
Perlmutter GPU node, take one first:

  salloc -N 1 -C gpu -q interactive -t 60 -A <account>

mqc takes the deck and nothing else -- a second argument is a hard error.
Verbosity is "system": {"logger": {"level": "verbose"}} in the deck itself;
levels are debug, verbose, info (default), performance, warning, error.
USAGE
}

require_positive_int() { # label, value
    [[ $2 =~ ^[1-9][0-9]*$ ]] || die "$1 must be a positive integer (got '$2')"
}

require_nonnegative_int() { # label, value
    [[ $2 =~ ^(0|[1-9][0-9]*)$ ]] || die "$1 must be a non-negative integer (got '$2')"
}

while [[ $# -gt 0 ]]; do
    # Every option in this list takes a value. Catching a missing one here is
    # what keeps `run.sh -nt` from dying as a bare "$2: unbound variable".
    case "$1" in
        -nt|--nthreads|-n|--ntasks|-c|--cpus-per-task|-g|--gpus-per-task|-b|--binary)
            (( $# >= 2 )) || die "option '$1' needs a value" ;;
    esac
    case "$1" in
        -nt|--nthreads)      NTHREADS="$2"; shift 2 ;;
        -n|--ntasks)         NTASKS="$2"; shift 2 ;;
        -c|--cpus-per-task)  CPUS_PER_TASK="$2"; shift 2 ;;
        -g|--gpus-per-task)  GPUS_PER_TASK="$2"; shift 2 ;;
        -b|--binary)         BINARY="$2"; shift 2 ;;
        --dry-run)           DRY_RUN=1; shift ;;
        --version)           INPUT="--version"; shift ;;
        -h|--help)           usage; exit 0 ;;
        -*)                  printf "run.sh: unknown option '%s'\n" "$1" >&2; usage >&2; exit 2 ;;
        *)
            [[ -z $INPUT ]] || die "only one input file may be given (got '$INPUT' and '$1')"
            INPUT="$1"; shift ;;
    esac
done

# Argument checks before environment checks, so a mistyped count is reported as
# a mistyped count and not as a missing executable.
require_positive_int "--ntasks" "$NTASKS"
[[ -z $NTHREADS ]] || require_positive_int "--nthreads" "$NTHREADS"
[[ -z $CPUS_PER_TASK ]] || require_positive_int "--cpus-per-task" "$CPUS_PER_TASK"
[[ -z $GPUS_PER_TASK ]] || require_nonnegative_int "--gpus-per-task" "$GPUS_PER_TASK"

if [[ -z $INPUT ]]; then
    printf 'run.sh: no input file given\n' >&2
    usage >&2
    exit 2
fi
[[ $INPUT == "--version" || -f $INPUT ]] || fail "input '$INPUT' not found"

find_binary() {
    # `cmake --preset default` lands in build/, every other preset in
    # build/<preset>/, so one hard-coded path is wrong on whichever machine did
    # not use the default preset -- perlmutter among them.
    local candidate
    for candidate in "$REPO_ROOT"/build/mqc "$REPO_ROOT"/build/*/mqc; do
        if [[ -x $candidate ]]; then
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    return 1
}

if [[ -z $BINARY ]]; then
    BINARY=$(find_binary) || BINARY=
fi
if [[ -z $BINARY || ! -x $BINARY ]]; then
    # A dry run is most useful before a build exists, so it may not insist on one.
    if (( ! DRY_RUN )); then
        fail "no mqc executable at ${BINARY:-$REPO_ROOT/build/mqc or $REPO_ROOT/build/*/mqc} -- build first?"
    fi
    warn "no executable found; showing the command for a nominal path"
    BINARY=${BINARY:-$REPO_ROOT/build/mqc}
fi

# Logical CPUs this job may use, and how many of them share a physical core.
TOTAL_CPUS=${SLURM_CPUS_ON_NODE:-$(nproc 2>/dev/null || echo 1)}
[[ $TOTAL_CPUS =~ ^[1-9][0-9]*$ ]] || TOTAL_CPUS=1
THREADS_PER_CORE=$(lscpu 2>/dev/null | awk -F: '/^Thread\(s\) per core/ {gsub(/ /, "", $2); print $2; exit}')
[[ $THREADS_PER_CORE =~ ^[1-9][0-9]*$ ]] || THREADS_PER_CORE=1
PHYSICAL_CPUS=$(( TOTAL_CPUS / THREADS_PER_CORE ))
(( PHYSICAL_CPUS > 0 )) || PHYSICAL_CPUS=1

if [[ -z $CPUS_PER_TASK ]]; then
    CPUS_PER_TASK=$(( TOTAL_CPUS / NTASKS ))
    # srun binds by core, so a -c that is not a whole number of cores splits a
    # hardware-thread pair between two ranks.
    CPUS_PER_TASK=$(( CPUS_PER_TASK - CPUS_PER_TASK % THREADS_PER_CORE ))
    (( CPUS_PER_TASK > 0 )) || CPUS_PER_TASK=$THREADS_PER_CORE
fi

if [[ -z $NTHREADS ]]; then
    # One thread per physical core in this rank's share. Sizing to the node
    # instead is how `-n 4` on a 64-core node asks for 256 threads.
    NTHREADS=$(( CPUS_PER_TASK / THREADS_PER_CORE ))
    (( NTHREADS > 0 )) || NTHREADS=1
fi

if (( NTASKS * NTHREADS > PHYSICAL_CPUS )); then
    warn "$NTASKS rank(s) x $NTHREADS thread(s) = $(( NTASKS * NTHREADS )) threads on $PHYSICAL_CPUS physical core(s)"
fi

if [[ -z $GPUS_PER_TASK ]]; then
    # cuEST wants a device per rank. Without this every rank in the job lands on
    # the same one, which looks like a slow run rather than a misconfigured one.
    if [[ ${SLURM_GPUS_ON_NODE:-0} =~ ^[1-9][0-9]*$ ]]; then
        GPUS_PER_TASK=1
    else
        GPUS_PER_TASK=0
    fi
fi

export OMP_NUM_THREADS="$NTHREADS"
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_STACKSIZE="${OMP_STACKSIZE:-64M}"
# The BLAS is sequential on purpose: fragment paths pin themselves to one OpenMP
# thread and parallelise over MPI, so a threaded BLAS competes with the ranks
# instead of helping them -- measured at 31 per cent slower on four ranks. See
# AGENTS.md, "Performance". Both are overridable for a single-molecule run.
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"

declare -a CMD
if [[ -n ${SLURM_JOB_ID:-} ]] && command -v srun >/dev/null 2>&1; then
    # SLURM_JOB_ID rather than just "srun exists": srun is on the login nodes
    # too, where it queues a fresh allocation instead of running anything here.
    CMD=(srun -n "$NTASKS" -c "$CPUS_PER_TASK" --cpu-bind=cores)
    if (( GPUS_PER_TASK > 0 )); then
        CMD+=(--gpus-per-task="$GPUS_PER_TASK")
    fi
    CMD+=("$BINARY" "$INPUT")
elif (( NTASKS > 1 )); then
    command -v mpirun >/dev/null 2>&1 ||
        fail "$NTASKS ranks asked for, but there is no SLURM job and no mpirun"
    CMD=(mpirun -np "$NTASKS" "$BINARY" "$INPUT")
else
    CMD=("$BINARY" "$INPUT")
fi

printf '# %s rank(s) x %s thread(s), -c %s' "$NTASKS" "$NTHREADS" "$CPUS_PER_TASK"
if (( GPUS_PER_TASK > 0 )); then
    printf ', %s gpu(s)/task' "$GPUS_PER_TASK"
fi
printf '\n# %s\n' "${CMD[*]}"

if (( DRY_RUN )); then
    printf '# OMP_NUM_THREADS=%s OMP_PROC_BIND=%s OMP_PLACES=%s OMP_STACKSIZE=%s\n' \
        "$OMP_NUM_THREADS" "$OMP_PROC_BIND" "$OMP_PLACES" "$OMP_STACKSIZE"
    printf '# MKL_NUM_THREADS=%s OPENBLAS_NUM_THREADS=%s\n' \
        "$MKL_NUM_THREADS" "$OPENBLAS_NUM_THREADS"
    exit 0
fi

exec "${CMD[@]}"
