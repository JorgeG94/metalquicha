#!/usr/bin/env python3
"""Compare our basis-function values against PySCF's.

    ./build/check_ao && python3 validation/check_ao.py

libcint does not evaluate basis functions, so the values `check_ao` writes come
from our own code and no integral call can confirm them. This compares them with
`pyscf.dft.numint.eval_ao`, fed this repository's own basis JSON.

**What is compared, and why it is not the raw numbers.** An SCF energy is
invariant under a diagonal rescaling of the basis, so two codes can normalise
individual basis functions differently, agree on every energy, and still disagree
elementwise here. PySCF does exactly that: it renormalises each contracted
function to unit self-overlap, while libcint takes the contraction coefficients as
given. So the test is on the *ratio*:

  * a **constant** ratio per function means the two agree on ordering, on sign and
    on the angular transform, and differ only by that scale -- which is the whole
    of what this comparison can prove and the only part that matters, since our
    density matrix is expressed in the same convention as our integrals;
  * a ratio that **varies across points** means the functions themselves differ:
    a permutation, or a different angular convention. That is a real failure.

**Cartesian bases are excluded, with a reason.** PySCF's `mol.cart=True`
convention is not libcint's -- the ratio there varies by two orders of magnitude
across points and changes sign, which is a permutation rather than a scale. Our
Cartesian path is checked instead by the overlap identity in
`test/test_mqc_libcint_ao.f90`, which compares against *our* integrals and is the
comparison that actually has to hold. Making PySCF's Cartesian ordering line up
would be a separate piece of work with no benefit to the DFT path.
"""
import sys
import pathlib
import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tools" / "cpu_validation"))
from gen_cpu_validation import bse_to_pyscf, molecule_form, CARTESIAN  # noqa: E402

# Spherical only; see the module docstring for why Cartesian is not comparable.
BASES = ["sto-3g", "cc-pvdz", "cc-pvtz"]

# How constant the per-function ratio has to be. This is the assertion that
# catches a permutation or a sign, and it is tight because a real ordering error
# moves it to order one.
SPREAD_TOL = 1e-7
# How far the scale itself may sit from one. PySCF's contraction renormalisation
# is a part in a million on these sets; anything larger is not that.
SCALE_TOL = 1e-5


def read_dump(path):
    v = path.read_text().split()
    n_points, n_ao, cart = int(v[0]), int(v[1]), v[2] in ("T", "t")
    rest = [float(x) for x in v[3:]]
    coords = np.array(rest[: 3 * n_points]).reshape(n_points, 3)
    ao = np.array(rest[3 * n_points: 3 * n_points + n_points * n_ao])
    return coords, ao.reshape(n_points, n_ao), cart


def main():
    from pyscf import gto
    from pyscf.dft import numint

    atoms = [("O", (0.0, 0.00000000009155, 0.10077199490609)),
             ("H", (0.0, 0.77250895271063, -0.46780199741728)),
             ("H", (0.0, -0.77250895280218, -0.46780199748881))]
    symbols = {a[0] for a in atoms}

    failures = 0
    for basis in BASES:
        path = pathlib.Path(f"/tmp/mqc_ao_{basis}.txt")
        if not path.exists():
            print(f"  MISSING {path} -- run ./build/check_ao first")
            failures += 1
            continue
        coords, ours, cart = read_dump(path)

        mol = gto.Mole()
        mol.atom = atoms
        mol.unit = "Angstrom"
        mol.basis = {s: bse_to_pyscf(basis, s) for s in symbols}
        mol.cart = molecule_form(basis, symbols) == CARTESIAN
        mol.verbose = 0
        mol.build()

        if mol.cart:
            print(f"  skip {basis}: Cartesian, not comparable (see docstring)")
            continue
        if mol.nao != ours.shape[1]:
            print(f"  FAIL {basis}: {ours.shape[1]} functions, PySCF has {mol.nao}")
            failures += 1
            continue

        theirs = numint.eval_ao(mol, coords)
        spreads, scales = [], []
        for mu in range(mol.nao):
            keep = np.abs(theirs[:, mu]) > 1e-8
            if keep.sum() < 2:
                continue
            r = ours[keep, mu] / theirs[keep, mu]
            spreads.append(r.max() - r.min())
            scales.append(r.mean())
        spread = max(spreads)
        off = max(abs(s - 1.0) for s in scales)
        ok = spread <= SPREAD_TOL and off <= SCALE_TOL
        print(f"  {'ok  ' if ok else 'FAIL'} {basis:10s} nao={mol.nao:3d}  "
              f"ratio spread {spread:.2e} (<= {SPREAD_TOL:.0e})  "
              f"|scale-1| {off:.2e} (<= {SCALE_TOL:.0e})")
        if not ok:
            failures += 1
            worst = int(np.argmax(spreads))
            print(f"       function {worst} has a non-constant ratio: the functions "
                  f"themselves differ, not just their normalisation")

    print()
    if failures:
        print(f"[ao] {failures} basis set(s) disagree")
        return 1
    print("[ao] ordering, signs and the angular transform match PySCF on every "
          "spherical basis")
    return 0


if __name__ == "__main__":
    sys.exit(main())
