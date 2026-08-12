#!/usr/bin/env python3
"""Compare our Boys localization against PySCF.

    ./build/check_localize && python3 validation/check_localize.py

**Comparing localizations needs care, and the obvious comparison is wrong.** A
Boys solution is defined only up to the ordering and sign of the localized
orbitals, so coefficient matrices cannot be compared elementwise. Worse, the
functional has many local maxima and the PySCF optimizer is strongly start
dependent: on water/STO-3G its own cost function gives 8.641 for the canonical
orbitals, 6.271 from its default init_guess="atomic", and 7.699 from the canonical
orbitals with that guess disabled. So "does our number equal their number" is not
a test of correctness -- it is a test of whose starting point was luckier.

What is a test, and is the sharp one:

1. **Stationarity under the PySCF cost function.** Hand PySCF our orbitals with
   its initial guess disabled and let it optimize. If it cannot improve on them,
   they are a genuine stationary point of the objective PySCF itself defines,
   which is a stronger statement than agreeing with one of its local maxima.
2. **Quality.** Our spread should be no worse than what PySCF reaches on its own.
   On these cases it is better, which is the right way round -- but that is a
   property of the starting point rather than of correctness, so it is reported
   and not asserted.

The PySCF cost function is r2 - sum_i |<i|r|i>|^2, the sum of orbital spreads,
which it minimizes. That is the same objective this module maximizes in the
equivalent form sum_i |<i|r|i>|^2, since sum_i <i|r^2|i> is invariant under a
rotation among the occupied orbitals.
"""

import sys
import pathlib

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tools" / "cpu_validation"))
from gen_cpu_validation import bse_to_pyscf, molecule_form, CARTESIAN  # noqa: E402

ATOMS = [("O", (0.0, 0.0, 0.0)),
         ("H", (0.0, 0.0, 0.9584)),
         ("H", (0.9268, 0.0, -0.2400))]

#: How much PySCF may improve on our orbitals before we call them unconverged.
#: The Jacobi sweep stops on a 1e-10 angle, so the spread is stationary to far
#: better than this.
STATIONARY_TOL = 1e-7

def read_dump(path):
    tokens = path.read_text().split()
    basis = tokens[0]
    values = [float(x) for x in tokens[1:]]
    nao, n_occ = int(values[0]), int(values[1])
    at = 2
    l_before, l_after = values[at], values[at + 1]
    at += 2
    centroids = np.array(values[at:at + 3 * n_occ]).reshape((3, n_occ), order="F")
    at += 3 * n_occ
    orbitals = np.array(values[at:at + nao * n_occ]).reshape((nao, n_occ), order="F")
    return dict(basis=basis, nao=nao, n_occ=n_occ, l_before=l_before,
                l_after=l_after, centroids=centroids, orbitals=orbitals)


def main():
    from pyscf import gto, scf, lo

    failures = 0
    for index in (1, 2):
        path = pathlib.Path(f"/tmp/mqc_localize_{index}.txt")
        if not path.exists():
            print(f"  MISSING {path} -- run ./build/check_localize first")
            failures += 1
            continue
        ours = read_dump(path)

        symbols = {a[0] for a in ATOMS}
        mol = gto.Mole()
        mol.atom = ATOMS
        mol.unit = "Angstrom"
        mol.basis = {s: bse_to_pyscf(ours["basis"], s) for s in symbols}
        mol.cart = molecule_form(ours["basis"], symbols) == CARTESIAN
        mol.verbose = 0
        mol.build()

        mf = scf.RHF(mol)
        mf.conv_tol = 1e-11
        mf.kernel()
        n_occ = ours["n_occ"]

        def spread(coeff):
            """The PySCF Boys cost function: the sum of orbital spreads."""
            return float(lo.Boys(mol, coeff).cost_function(np.eye(coeff.shape[1])))

        ours_spread = spread(ours["orbitals"])

        # 1. Can PySCF improve on our orbitals? init_guess must be disabled --
        #    by default Boys discards whatever it was given in favour of an
        #    atomic guess, which is exactly why a naive comparison looks like a
        #    disagreement when there is none.
        refine = lo.Boys(mol, ours["orbitals"])
        refine.init_guess = None
        refine.conv_tol = 1e-12
        refine.max_cycle = 1000
        refined_spread = spread(refine.kernel())
        gain = ours_spread - refined_spread

        # 2. And what does PySCF reach by itself, for context?
        theirs = lo.Boys(mol, mf.mo_coeff[:, :n_occ])
        theirs.conv_tol = 1e-12
        theirs.max_cycle = 1000
        their_spread = spread(theirs.kernel())

        ok = gain < STATIONARY_TOL
        status = "ok  " if ok else "FAIL"
        print(f"  {status} {ours['basis']:8s} n_occ {n_occ}  spreads: "
              f"ours {ours_spread:.8f}   PySCF refining ours {refined_spread:.8f}"
              f"   PySCF alone {their_spread:.8f}")
        if not ok:
            failures += 1
            print(f"       PySCF improved our spread by {gain:.2e}, so ours is not "
                  f"a stationary point")
        elif their_spread < ours_spread - 1e-8:
            print(f"       note: PySCF localized better by "
                  f"{ours_spread - their_spread:.2e}")

        # Water should be recognisable: an oxygen core, two O-H bonds, and two
        # equivalent lone pairs off the molecular plane. Printed for a human; the
        # assertion above is what a machine checks.
        for i in range(n_occ):
            c = ours["centroids"][:, i]
            print(f"        centroid {i + 1}: "
                  f"({c[0]:8.4f},{c[1]:8.4f},{c[2]:8.4f})  |r| = "
                  f"{np.linalg.norm(c):7.4f}")

    print()
    if failures:
        print(f"[localize] {failures} FAILURE(S)")
        return 1
    print("[localize] our orbitals are stationary under the PySCF Boys objective")
    return 0


if __name__ == "__main__":
    sys.exit(main())
