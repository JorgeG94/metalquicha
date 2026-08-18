# scratch

Throwaway Fortran drivers, for poking at the library by hand.

Nothing in here is compiled by the CMake build or run by `ctest`, and
everything except this file and `build.sh` is gitignored. Write whatever you
like; it will not be committed and it cannot break CI.

```
./scratch/build.sh scratch/whatever.f90          # uses build-dh/
./scratch/build.sh scratch/whatever.f90 build     # or another build tree
```

The build tree has to exist already -- this compiles one file against it, it
does not configure anything.

## What this is for, and what it is not

Use it when you want to *look* at something: print a link table, dump a matrix,
time a loop, check one number against PySCF before deciding what the assertion
should be.

Once you know what the answer should be, the check belongs in `test/` as a
test-drive unit test, not here. A driver in this directory proves nothing to
anyone else and is not run again after you close the terminal.
