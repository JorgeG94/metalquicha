#!/usr/bin/env python3
"""Compare our exchange-correlation energy against PySCF's.

    ./build/check_dft && python3 validation/check_dft.py

`check_dft` converges RHF, evaluates the density on our own DFT grid, and
integrates an LDA exchange-correlation energy through libxc. It dumps the grid and
the density so this script can separate two questions that a single number would
confound:

  1. **Is the functional called correctly?** Given *our* rho and *our* weights,
     PySCF evaluates the same two libxc functionals and integrates them the same
     way. Both codes are calling the same library on the same numbers, so this
     must agree to roundoff -- it is a test of our quadrature and our factor of
     rho, nothing else.
  2. **Is the density right?** PySCF builds its own converged density and
     evaluates it on *our* grid points with *our* weights. The grid therefore
     cancels out of the comparison entirely and what is left is the SCF and the AO
     evaluation, which is why this agrees to 4e-10 rather than to grid accuracy --
     a much stronger statement than comparing two independently built grids, and
     the reason the tolerance is 1e-8 and not 1e-5.

Splitting them matters because a factor of two in the density and a missing rho in
the quadrature both move a single end-to-end number, in opposite directions.
"""
import sys
import pathlib
import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tools" / "cpu_validation"))
from gen_cpu_validation import bse_to_pyscf, molecule_form, CARTESIAN  # noqa: E402

CASES = [("cc-pvdz", 3), ("cc-pvdz", 5), ("sto-3g", 4)]
FUNCTIONAL_TOL = 1e-10   # same library, same numbers: roundoff only
CHAIN_TOL = 1e-8         # the grid cancels; see below


def read_dump(path):
    v = path.read_text().split()
    n_points, n_ao = int(v[0]), int(v[1])
    rest = [float(x) for x in v[2:]]
    n_elec, e_x, e_c = rest[0], rest[1], rest[2]
    body = np.array(rest[3:3 + 4 * n_points]).reshape(n_points, 4)
    rho = np.array(rest[3 + 4 * n_points: 3 + 4 * n_points + n_points])
    return dict(coords=body[:, :3], weights=body[:, 3], rho=rho,
                n_elec=n_elec, e_x=e_x, e_c=e_c, n_ao=n_ao)


def main():
    from pyscf import gto, scf
    from pyscf.dft import libxc, numint

    atoms = [("O", (0.0, 0.00000000009155, 0.10077199490609)),
             ("H", (0.0, 0.77250895271063, -0.46780199741728)),
             ("H", (0.0, -0.77250895280218, -0.46780199748881))]
    symbols = {a[0] for a in atoms}

    failures = 0
    for basis, level in CASES:
        path = pathlib.Path(f"/tmp/mqc_dft_{basis}_L{level}.txt")
        if not path.exists():
            print(f"  MISSING {path} -- run ./build/check_dft first")
            failures += 1
            continue
        d = read_dump(path)

        # ---- 1. the functional, on our own density and weights --------------
        ex = libxc.eval_xc("LDA_X", d["rho"])[0]
        ec = libxc.eval_xc("LDA_C_VWN", d["rho"])[0]
        px = float(np.dot(d["weights"] * d["rho"], ex))
        pc = float(np.dot(d["weights"] * d["rho"], ec))
        dx, dc = abs(px - d["e_x"]), abs(pc - d["e_c"])
        ok1 = dx < FUNCTIONAL_TOL and dc < FUNCTIONAL_TOL

        # ---- 2. the whole chain, PySCF's density on our grid ----------------
        mol = gto.Mole()
        mol.atom = atoms
        mol.unit = "Angstrom"
        mol.basis = {s: bse_to_pyscf(basis, s) for s in symbols}
        mol.cart = molecule_form(basis, symbols) == CARTESIAN
        mol.verbose = 0
        mol.build()
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-12
        mf.kernel()
        ao = numint.eval_ao(mol, d["coords"])
        rho_p = numint.eval_rho(mol, ao, mf.make_rdm1())
        n_p = float(np.dot(d["weights"], rho_p))
        exp_ = libxc.eval_xc("LDA_X", rho_p)[0]
        ecp = libxc.eval_xc("LDA_C_VWN", rho_p)[0]
        cx = float(np.dot(d["weights"] * rho_p, exp_))
        cc = float(np.dot(d["weights"] * rho_p, ecp))
        ex_chain, ec_chain = abs(cx - d["e_x"]), abs(cc - d["e_c"])
        n_chain = abs(n_p - d["n_elec"])
        ok2 = ex_chain < CHAIN_TOL and ec_chain < CHAIN_TOL and n_chain < CHAIN_TOL

        print(f"  {basis:9s} L{level}  points {len(d['weights']):6d}  N={d['n_elec']:.7f}")
        print(f"      functional  dEx {dx:.2e}  dEc {dc:.2e}   "
              f"{'ok' if ok1 else 'FAIL'} (<= {FUNCTIONAL_TOL:.0e})")
        print(f"      whole chain dEx {ex_chain:.2e}  dEc {ec_chain:.2e}  dN {n_chain:.2e}  "
              f"{'ok' if ok2 else 'FAIL'} (<= {CHAIN_TOL:.0e})")
        if not (ok1 and ok2):
            failures += 1

    print()
    if failures:
        print(f"[dft] {failures} case(s) disagree")
        return 1
    print("[dft] exchange-correlation energies match PySCF, functional and chain")
    return 0


if __name__ == "__main__":
    sys.exit(main())
