#!/usr/bin/env bash
# Compile and run one throwaway driver against an existing build tree.
#
#   ./scratch/build.sh scratch/dump_link_table.f90 [build-dir]
#
# Nothing here is part of the CMake build. These are for looking at something
# by hand -- printing a table, timing a loop, checking a number against PySCF --
# in cases where a test-drive assertion is not what you want yet.
#
# Links against the shared library rather than the static one, deliberately: the
# static archive carries only the Fortran objects, and the C halves of libcint
# and libxc live in their own archives whose link order is not worth
# reconstructing here. libmqc.so already has all of it resolved.
set -euo pipefail
src="${1:?usage: build.sh <source.f90> [build-dir]}"
build="${2:-build-dh}"
out="${src%.f90}"
abs_build="$(cd "$build" && pwd)"

gfortran -O2 -fopenmp -o "$out" "$src" \
    -I"$abs_build/modules" \
    -L"$abs_build" -lmqc \
    -Wl,-rpath,"$abs_build"

echo "built $out"
"$out"
