#!/usr/bin/env bash
# Pre-commit hook: the four places that state metalquicha's version agree.
#
#   CMakeLists.txt            project(... VERSION x.y.z ...)
#   fpm.toml                  version = "x.y.z"
#   src/core/mqc_version.f90  MQC_VERSION_STR = "x.y.z"  (what `mqc` prints)
#   mqc_docs/source/conf.py   release = 'x.y.z'          (what the docs say)
#
# They are kept in step by hand, and they have drifted before: the v0.1.0
# and v0.1.1 binaries reported themselves as 0.2.0, because CMake and the
# version module had been bumped and fpm and the docs had not.  A release
# changes all four in one commit; this refuses a commit that changes some.
#
# Run from the repository root, which is where pre-commit runs hooks.
set -euo pipefail

extract() {
    # $1 file, $2 sed expression that prints the version and nothing else
    local v
    v=$(sed -n -E "$2" "$1" | head -n 1)
    if [ -z "$v" ]; then
        echo "versions_agree: no version found in $1" >&2
        exit 1
    fi
    printf '%s' "$v"
}

cmake=$(extract CMakeLists.txt 's/^[[:space:]]*VERSION[[:space:]]+([0-9]+\.[0-9]+\.[0-9]+).*/\1/p')
fpm=$(extract fpm.toml 's/^version[[:space:]]*=[[:space:]]*"([^"]+)".*/\1/p')
fortran=$(extract src/core/mqc_version.f90 's/.*MQC_VERSION_STR[[:space:]]*=[[:space:]]*"([^"]+)".*/\1/p')
docs=$(extract mqc_docs/source/conf.py "s/^release[[:space:]]*=[[:space:]]*'([^']+)'.*/\\1/p")

if [ "$cmake" = "$fpm" ] && [ "$cmake" = "$fortran" ] && [ "$cmake" = "$docs" ]; then
    exit 0
fi

cat >&2 <<EOF
versions_agree: the version is not the same everywhere
  CMakeLists.txt            $cmake
  fpm.toml                  $fpm
  src/core/mqc_version.f90  $fortran
  mqc_docs/source/conf.py   $docs
Change all four together.
EOF
exit 1
