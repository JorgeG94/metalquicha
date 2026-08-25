#!/usr/bin/env python3
"""Check that a Fortran-written amplitude file is the right tensor to C.

Run `check_amplitude_layout` first; it writes `amplitude_layout.h5` and
asserts nothing about it.

WHY THIS IS NOT A FORTRAN TEST
------------------------------
It cannot be one. HDF5 declares a dataset in C order -- first extent
slowest-varying -- while Fortran's leftmost index varies fastest. Passing
`shape(values)` straight into `H5Screate_simple` therefore stores the correct
bytes under an incorrect declared shape: `t1(2,3)` becomes a 2x3 dataset whose
rows interleave the two indices. Reading it back through the same bindings is
exact, because the same convention is wrong on both sides, so every assertion
a Fortran test could make passes. The file is only wrong to somebody else --
h5py, h5dump, MATLAB, any C consumer -- and being readable by somebody else is
the entire reason to use HDF5 rather than a stream of unformatted records.

This ran as a manual cross-check once and immediately found the bug it is
named for; it is here so that it runs every time instead.

WHAT IS CHECKED, AND WHY BOTH HALVES ARE NEEDED
-----------------------------------------------
`h5dump -H` would show the declared extents and is free -- it ships inside
the HDF5 package. It is not enough on its own: a "fix" that transposed the
array in Fortran before writing, rather than reversing the extents, prints the
right header over the wrong numbers. So the shapes are checked *and* every
element is checked against the indices it encodes.
"""
import os
import sys

import h5py
import numpy as np

PATH = "amplitude_layout.h5"
N1, N2, N3, N4 = 2, 3, 4, 5
ITERATION = 3
ENERGY = -0.5


def main():
    failures = []

    with h5py.File(PATH, "r") as f:
        # Attributes first: `complete` is what a resumed run trusts, and a
        # file that says 0 here would be silently ignored rather than used.
        got = f.attrs.get("complete")
        if got != 1:
            failures.append(f"complete is {got!r}, expected 1")
        got = f.attrs.get("iteration")
        if got != ITERATION:
            failures.append(f"iteration is {got!r}, expected {ITERATION}")
        got = f.attrs.get("energy")
        if got != ENERGY:
            failures.append(f"energy is {got!r}, expected {ENERGY}")

        t1 = f["t1"][...]
        t2 = f["t2"][...]

    # Shape: reversed, because that is what C order means. All four extents
    # of t2 differ, so a partial or wrong-axis reversal cannot pass.
    if t1.shape != (N2, N1):
        failures.append(f"t1 shape is {t1.shape}, expected {(N2, N1)} "
                        f"(Fortran t1({N1},{N2}) reversed)")
    if t2.shape != (N4, N3, N2, N1):
        failures.append(f"t2 shape is {t2.shape}, expected {(N4, N3, N2, N1)} "
                        f"(Fortran t2({N1},{N2},{N3},{N4}) reversed)")

    # Values: the half a header dump cannot see. t1_h5[a, i] must be the
    # element Fortran calls t1(i+1, a+1).
    if t1.shape == (N2, N1):
        want = np.fromfunction(
            lambda a, i: 100 * (i + 1) + (a + 1), (N2, N1), dtype=float)
        if not np.array_equal(t1, want):
            failures.append("t1 values do not match t1_h5[a,i] == 100*i + a; "
                            f"got\n{t1}\nwanted\n{want}")

    if t2.shape == (N4, N3, N2, N1):
        want = np.fromfunction(
            lambda b, j, a, i: 1000 * (i + 1) + 100 * (a + 1) + 10 * (j + 1) + (b + 1),
            (N4, N3, N2, N1), dtype=float)
        if not np.array_equal(t2, want):
            bad = np.argwhere(t2 != want)
            failures.append(
                f"t2 values do not match t2_h5[b,j,a,i] == 1000*i + 100*a + 10*j + b; "
                f"{len(bad)} of {t2.size} elements differ, first at {tuple(bad[0])}")

    if failures:
        print(f"FAIL: {PATH} is not the tensor a C-order reader expects")
        for line in failures:
            print(f"  - {line}")
        return 1

    print(f"OK: {PATH} reads correctly from h5py "
          f"(HDF5 {h5py.version.hdf5_version}, h5py {h5py.__version__})")
    print(f"  t1 {t1.shape} and t2 {t2.shape}, both reversed as C order requires")
    # Removed only on success, so a failing run leaves the evidence to h5dump.
    os.remove(PATH)
    return 0


if __name__ == "__main__":
    sys.exit(main())
