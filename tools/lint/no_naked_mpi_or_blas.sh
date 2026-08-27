#!/usr/bin/env bash
# Pre-commit hook: ban direct MPI and direct BLAS/LAPACK in metalquicha.
#
# Every MPI call goes through `pic_mpi_lib` and every BLAS/LAPACK call
# through `pic_blas`.  The MPI one is the load-bearing rule: pic-mpi is
# the single layer that knows whether the build has `mpi_f08` or the
# legacy `mpi` module, and some compilers ship only one of the two.  A
# naked `use mpi_f08` compiles fine here and breaks the build on someone
# else's machine, which is why no test catches it.  pic-blas is the same
# argument one step down -- it hides which BLAS implementation is
# underneath, and gives a call that reads as what it computes.
#
# Both are FORTRAN_STYLE.md, "No Naked MPI or BLAS Calls".  The tree was
# at zero violations when this was written, so the rule is enforced with
# no grandfathered exemptions.
#
# The documented exception is pic-mpi or pic-blas genuinely lacking the
# functionality.  Take it with a reason on the line, or the line above:
#
#     ! mqc: allow-naked-blas -- pic_blas has no ?syevr, and this needs
#     ! the selected-eigenvalue form
#     call dsyevr(...)
#
# Receives candidate files as args from pre-commit.  Also runs manually:
#
#     bash tools/lint/no_naked_mpi_or_blas.sh src/scf/*.f90
#     bash tools/lint/no_naked_mpi_or_blas.sh --all
set -euo pipefail

SCAN_DIRS=(src app backends test validation tools)

if [[ "${1:-}" == "--all" ]]; then
    shift
    mapfile -t files < <(find "${SCAN_DIRS[@]}" \
        \( -name '*.f90' -o -name '*.F90' -o -name '*.f' -o -name '*.F' \) \
        -type f 2>/dev/null | sort)
    set -- "${files[@]}"
fi

[[ $# -eq 0 ]] && exit 0

# BLAS/LAPACK routine names, without the type prefix.  Kept to routines
# pic_blas actually covers plus the LAPACK drivers that turn up here; a
# name pic-blas has no wrapper for is a contribution upstream, not a
# silent direct call.
BLAS_NAMES='gemm|gemv|ger|symm|syrk|syr2k|trmm|trsm|trsv|axpy|copy|scal|swap|dot|dotu|dotc|nrm2|asum|potrf|potri|potrs|getrf|getri|getrs|gesv|posv|syev|syevd|syevr|syevx|heev|heevd|geev|gesvd|gesdd|gelss|gels|orgqr|geqrf'

# Strip full-line and trailing comments before matching, so a "bad
# example" written into a docstring is not a violation.  A `!` inside a
# string literal takes the rest of the line with it -- that can only
# hide a violation inside a string, which is not code.
strip_comments() { sed 's/!.*$//'; }

# A line carrying the escape hatch, or preceded by one.
suppressed() {
    local file=$1 line=$2 rule=$3
    sed -n "$((line > 1 ? line - 1 : 1)),${line}p" "$file" \
        | grep -qiE "mqc:[[:space:]]*allow-naked-${rule}\b"
}

status=0

report() {
    local rule=$1 pattern=$2 message=$3
    shift 3
    local hits=""
    for f in "$@"; do
        [[ -f "$f" ]] || continue
        while IFS=: read -r lineno _; do
            [[ -n "$lineno" ]] || continue
            suppressed "$f" "$lineno" "$rule" && continue
            hits+="$f:$lineno: $(sed -n "${lineno}p" "$f" | sed 's/^[[:space:]]*//')"$'\n'
        done < <(strip_comments < "$f" | grep -nEi "$pattern" || true)
    done
    if [[ -n "$hits" ]]; then
        printf '%s' "$hits"
        echo
        echo "ERROR: $message"
        echo
        status=1
    fi
}

report mpi \
    '(^[[:space:]]*use[[:space:]]+(mpi|mpi_f08)\b|^[[:space:]]*include[[:space:]]+.[[:space:]]*mpif\.h|\bcall[[:space:]]+MPI_[A-Za-z_]+)' \
    "direct MPI is banned. Use pic_mpi_lib (comm_t, send, recv, bcast,
       allgather, isend, irecv, wait, iprobe, abort_comm). pic-mpi is what
       lets the mpi_f08 and legacy mpi backends be swapped; a direct call
       builds here and fails on a compiler that ships only the other one.
       See FORTRAN_STYLE.md, 'No Naked MPI or BLAS Calls'." \
    "$@"

report blas \
    "(\bcall[[:space:]]+[sdcz]($BLAS_NAMES)[[:space:]]*\(|[^[:alnum:]_%]([sdcz](dot|dotu|dotc|nrm2|asum)|i[sdcz]amax)[[:space:]]*\()" \
    "direct BLAS/LAPACK is banned. Use pic_blas (pic_gemm, pic_dot, ...).
       If pic-blas has no wrapper for what you need, add one upstream
       rather than calling through.
       See FORTRAN_STYLE.md, 'No Naked MPI or BLAS Calls'." \
    "$@"

exit $status
