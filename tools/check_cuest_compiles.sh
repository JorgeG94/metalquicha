#!/usr/bin/env bash
#
# Compile the cuEST backend with gfortran, without cuEST and without a GPU.
#
# Nothing else in the tree does this. MQC_ENABLE_CUEST is OFF by default and
# turning it on needs CUEST_ROOT and an sm_80 card, so on any machine without
# one -- which is every CI runner and most developer boxes -- backends/cuest/
# is never handed to a compiler at all. It can therefore rot silently: a
# renamed field or a changed signature elsewhere in the tree breaks it and
# nothing says so until someone with a card tries to build.
#
# The bindings in backends/cuest/bindings/ are plain `bind(C)` interfaces, so
# the library is needed to *link* but not to *compile*. That is what makes this
# possible: gfortran type-checks every call against those interfaces, which is
# where the mistakes actually are -- a wrong integer kind, an argument in the
# wrong position, a missing `use`.
#
# What this does NOT do: link, run, or say anything about whether cuEST returns
# the right numbers. It is a compile check.
#
# Usage:  tools/check_cuest_compiles.sh [build-dir]
#
# The build directory only supplies the .mod files for the rest of the program,
# so any existing CPU build will do.
set -euo pipefail

BUILD=${1:-build}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

if [[ ! -d $BUILD/modules ]]; then
   echo "error: $BUILD/modules not found." >&2
   echo "Build the program first, e.g. cmake -B $BUILD && cmake --build $BUILD" >&2
   exit 2
fi

FC=${FC:-gfortran}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

INC=(-I"$BUILD/modules" -I"$BUILD/modules_shared" -I"$WORK")

# The dependencies that carry Fortran modules, named one by one.
#
# A CMake target never needs this: pic, pic-blas and pic-mpi each export their
# module directory as a BUILD_INTERFACE include, so linking `pic::pic` carries
# it and cannot be got wrong. This script compiles by hand, outside CMake, and
# so has to say where they are.
#
# Spelled out rather than globbed over `_deps/*/modules`. A glob would keep
# working when a dependency moved its modules, which sounds like a virtue and
# is not: it would also keep working when one went missing, and quietly put
# some other library's directory on the include path. A name that stops being
# right should stop the build.
#
# pic-mpi is deliberately absent. v0.6.2 still writes its .mod files into
# `${CMAKE_BINARY_DIR}/modules` -- this build tree's own, not its -- so
# `pic_mpi_lib.mod` arrives through the first `-I` above. Its main has since
# moved to the same convention pic and pic-blas use, but does not yet create
# the directory at configure time, so a consumer cannot generate against it.
# When a release carries that too, pic-mpi joins this loop.

for _dep in pic pic-blas; do
   _moddir="$BUILD/_deps/${_dep}-build/modules"
   if [[ ! -d $_moddir ]]; then
      echo "error: $_moddir not found." >&2
      echo "$_dep either did not build or has moved its .mod files; this script" >&2
      echo "names the directory and has to be told when it changes." >&2
      exit 2
   fi
   INC+=(-I"$_moddir")
done

if [[ -d $BUILD/_deps/jsonfortran-build/include ]]; then
   INC+=(-I"$BUILD/_deps/jsonfortran-build/include")
fi

# The C-interface layer first: everything else uses it.
BINDINGS=(cuda_runtime cublas cusolver cuest cuda_helpers cuest_helpers)

# Backend modules, roughly in dependency order. The order is not load-bearing:
# the loop below retries until it stops making progress, so a new file can be
# appended anywhere.
BACKEND=(mqc_cuest_runtime mqc_cuest_basis mqc_cuest_grid mqc_cuest_functionals
         mqc_cuest_context mqc_diis_device mqc_cuest_integrals mqc_cuest_scf
         mqc_cuest_atomic_guess mqc_cuest_gradient mqc_cuest_driver
         mqc_cuest_bridge)

fail=0
# Either extension: `cuda_runtime` is `.F90` because it carries the
# `#if MQC_CUDA_VERSION_MAJOR` guards for the entry points CUDA 13 renamed, and
# the capital F is what gets it preprocessed. This used to assume `.f90` for
# everything, so the rename turned into "Cannot open file
# backends/cuest/bindings/cuda_runtime.f90" -- a message that reads like a
# missing file rather than a script that looked in one place.
#
# The guard is left undefined here on purpose. This is a compile check with no
# CUDA toolkit in reach, so an undefined macro takes the pre-13 branch, which is
# the one that has to keep compiling for every build that is not CUDA 13.
src_for() { # dir, name -> path
   if [[ -f "$1/$2.F90" ]]; then
      printf '%s/%s.F90' "$1" "$2"
   else
      printf '%s/%s.f90' "$1" "$2"
   fi
}

compile() { # dir, name
   local src
   src=$(src_for "$1" "$2")
   if ! $FC -c "${INC[@]}" -J"$WORK" "$src" -o "$WORK/$2.o" 2>"$WORK/$2.err"; then
      return 1
   fi
   return 0
}

for m in "${BINDINGS[@]}"; do
   if ! compile backends/cuest/bindings "$m"; then
      echo "FAILED: $(src_for backends/cuest/bindings "$m")"
      sed -n '1,25p' "$WORK/$m.err"
      fail=1
   fi
done
[[ $fail -eq 0 ]] || { echo; echo "cuEST bindings do not compile."; exit 1; }

todo=("${BACKEND[@]}")
for _ in 1 2 3 4; do
   left=()
   for m in "${todo[@]}"; do
      compile backends/cuest/backend "$m" || left+=("$m")
   done
   todo=("${left[@]:-}")
   [[ ${#todo[@]} -eq 0 || -z ${todo[0]} ]] && break
done

if [[ ${#todo[@]} -gt 0 && -n ${todo[0]} ]]; then
   for m in "${todo[@]}"; do
      echo "FAILED: backends/cuest/backend/$m.f90"
      sed -n '1,25p' "$WORK/$m.err"
      echo
   done
   exit 1
fi

echo "cuEST backend compiles: ${#BINDINGS[@]} bindings, ${#BACKEND[@]} backend modules."
