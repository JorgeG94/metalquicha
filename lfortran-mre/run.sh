#!/usr/bin/env bash
# Compile each MRE and report the outcome against LFortran.
#
# 01-05 are OpenMP bugs and need --openmp to show up at all; 06-07 are not
# OpenMP and are compiled the way the serial build compiles this project.
# Every one of these files is accepted by `gfortran -fopenmp`.
#
# Expected on LFortran 0.64.0:
#   01 semantic error   02 SEGFAULT   02b ok   03 semantic error
#   04 semantic error   05 semantic error   06 semantic error   07 ICE
for f in 0*.f90; do
   case $f in 0[1-5]*) flags="--openmp" ;; *) flags="" ;; esac
   printf '%-40s ' "$f"
   out=$(lfortran $flags -c "$f" -o /dev/null 2>&1); rc=$?
   # shellcheck disable=SC2001  # stripping ANSI escapes needs a regex class
   clean=$(sed 's/\x1b\[[0-9;]*m//g' <<<"$out")
   if   [ $rc -eq 0 ];   then echo "ok"
   elif [ $rc -eq 139 ]; then echo "SEGFAULT"
   elif grep -q LCompilersException <<<"$clean"; then
        echo "ICE: $(grep -m1 LCompilersException <<<"$clean")"
   else echo "exit $rc: $(grep -m1 -oP '(semantic|syntax) error: .*' <<<"$clean")"
   fi
done
