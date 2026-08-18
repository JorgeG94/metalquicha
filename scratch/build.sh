#!/usr/bin/env bash
# Compile and run one throwaway driver against an existing build tree.
#
#   ./scratch/build.sh scratch/dump_link_table.f90 [build-dir]
#
# Nothing here is part of the CMake build. These are for looking at something
# by hand -- printing a table, timing a loop, checking a number against PySCF --
# in cases where a test-drive assertion is not what you want yet.
set -euo pipefail
src="${1:?usage: build.sh <source.f90> [build-dir]}"
build="${2:-build-dh}"
out="${src%.f90}"

gfortran -O2 -fopenmp -o "$out" "$src" \
    -I"$build/modules" \
    "$build/libmetalquicha.a" \
    $(find "$build/_deps" -name '*.a' | tr '\n' ' ') \
    -llapack -lblas

echo "built $out"
"$out"
