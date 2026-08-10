# Lebedev grid data

`mqc_lebedev_data.f90` is generated, not written. These three scripts produce
and check it. They need PySCF, which is a development dependency only -- the
build does not use them.

```bash
python extract_lebedev.py > lebedev.json          # orbit parameters + reference points
python gen_lebedev_fortran.py lebedev.json > ../../src/methods/mqc_lebedev_data.f90
```

Then build `check_lebedev`, run it, and compare:

```bash
./check_lebedev                                    # self-consistency, writes lebedev_points.txt
python compare_lebedev.py lebedev_points.txt lebedev.json
```

## What is taken from where

The parameters are the published Lebedev-Laikov values [Lebedev & Laikov,
Dokl. Math. 59, 477 (1999)]; PySCF is used as a convenient carrier of them, and
`extract_lebedev.py` records them by instrumenting `SphGenOh` rather than
parsing source. The orbit expansion in `mqc_lebedev.f90` is written from the
octahedral symmetry, not transliterated -- so the comparison is between two
independent expansions of the same parameters, and agreement confirms both.

PySCF's own Lebedev code arrived there as Fortran -> C++ -> Python; porting it
back to Fortran would have been a round trip, and PySCF is Apache-2.0 against
this project's MIT.

## Result

32 orders, 46976 points. Coordinates agree with PySCF to 9.0e-16, weights
exactly. Orders 74, 230 and 266 have negative weights in both -- a real
property of those grids, not an error.
