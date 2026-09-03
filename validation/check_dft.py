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
# Self-consistent Kohn-Sham: functional, PySCF's name for it, basis, grid level,
# and whether this is the direct, Schwarz-screened Fock build rather than the
# in-core one. Both are compared, because the two builds write separate files --
# screening that changed an answer would show up only in the direct rows.
KS_CASES = [("lda_x", "LDA_X", "cc-pvdz", 3, False), ("svwn", "SVWN", "cc-pvdz", 3, False),
            ("gga_x_pbe", "GGA_X_PBE", "cc-pvdz", 3, False), ("pbe", "PBE", "cc-pvdz", 3, False),
            ("b3lyp", "B3LYP", "cc-pvdz", 3, False), ("pbe0", "PBE0", "cc-pvdz", 3, False),
            ("tpss", "TPSS", "cc-pvdz", 3, False), ("m06-l", "M06-L", "cc-pvdz", 3, False),
            ("svwn", "SVWN", "cc-pvdz", 3, True), ("b3lyp", "B3LYP", "cc-pvdz", 3, True),
            ("tpss", "TPSS", "cc-pvdz", 3, True),
            # Range-separated hybrids, and the sharpest test of the attenuated
            # integrals: omega reaches libcint only through the long-range
            # exchange pass, and the total is wrong by Hartrees if it fails to
            # arrive. Direct only -- the in-core path holds one full-range
            # exchange tensor and has nothing to attenuate.
            ("wb97x", "wb97x", "cc-pvdz", 3, True),
            ("cam-b3lyp", "camb3lyp", "cc-pvdz", 3, True)]

# Unrestricted Kohn-Sham: molecule, functional, PySCF's name, basis, grid level.
#
# The methyl radical for the open-shell cases -- a doublet whose ground state is
# unambiguous, which matters because two codes' UKS can converge to different
# stationary points and then disagree for a reason that is not a bug. <S^2> is
# compared alongside the energy so that a disagreement of that kind says so.
#
# Closed-shell water is here for a different reason: it must reproduce the
# *restricted* energy exactly. That is the sharper test of the interleaved arrays
# and the cross-spin gradient term, and it is asserted in test_mqc_czt_uks too.
UKS_CASES = [("ch3", "svwn", "SVWN", "cc-pvdz", 3),
             ("ch3", "pbe", "PBE", "cc-pvdz", 3),
             ("ch3", "b3lyp", "B3LYP", "cc-pvdz", 3),
             ("ch3", "tpss", "TPSS", "cc-pvdz", 3),
             ("ch3", "wb97x", "wb97x", "cc-pvdz", 3),
             ("water", "pbe", "PBE", "cc-pvdz", 3),
             ("water", "b3lyp", "B3LYP", "cc-pvdz", 3)]
UKS_S2_TOL = 1e-6

# Geometries the Fortran side uses, repeated because there is no shared reader.
CH3_ATOMS = [("C", (0.0, 0.0, 0.0)),
             ("H", (1.079, 0.0, 0.0)),
             ("H", (-0.5395, 0.93444, 0.0)),
             ("H", (-0.5395, -0.93444, 0.0))]
CH3_SPIN = 1   # PySCF counts unpaired electrons, not 2S+1
KS_TOL = 1e-9   # measured 2e-11; the SCF thresholds set the floor


def check_pyscf_grid_convention():
    """Refuse a PySCF whose radial grid is not the one we integrate on.

    PySCF changed its default. From roughly 2.7 it applies the per-element xi
    scale factors of the Treutler-Ahlrichs mesh (JCP 102, 346), gated on
    `pyscf.dft.radi.ATOM_SPECIFIC_TREUTLER_GRIDS`; 2.6.2 and earlier use xi = 1
    for every element. We apply them -- src/methods/mqc_dft_radial.f90 -- so
    against an older PySCF the two codes integrate on different grids.

    That failure is silent and convincing. Every Kohn-Sham case misses by 1e-8
    to 1e-6 while the two checks above still pass, because those are the ones
    where the grid cancels; the size ordering follows grid sensitivity, worst
    for M06-L and the B97 series. It reads as a regression in the
    exchange-correlation code and is not one. It cost an afternoon once.

    The usual cause is the interpreter rather than the install: `python3.10`
    picks up ~/.local packages, which are visible to any Python 3.10 whether or
    not a venv is active. Use the venv's `python3`.
    """
    import pyscf
    from pyscf.dft import radi

    on = getattr(radi, "ATOM_SPECIFIC_TREUTLER_GRIDS", False)
    if not on:
        raise SystemExit(
            f"[dft] refusing to compare against PySCF {pyscf.__version__} at "
            f"{pyscf.__file__}\n"
            "      Its Treutler-Ahlrichs radial mesh has no per-element xi scale "
            "factors\n"
            "      (pyscf.dft.radi.ATOM_SPECIFIC_TREUTLER_GRIDS is absent or "
            "false), so it\n"
            "      integrates on a different grid than we do and every "
            "Kohn-Sham case would\n"
            "      appear to fail by 1e-8 to 1e-6. Use the project venv's "
            "python3 (PySCF\n"
            "      >= 2.7), not python3.10. See mqc_docs/source/validation.rst, "
            "'Density\n"
            "      Functional Energies'."
        )

# Double hybrids, against pyscf-forge's DFDH on the same geometry and basis.
#
# The tolerance is 5e-5 and that is not slack -- it is DFDH's own error. Its SCF is
# density-fitted with an auto-generated even-tempered auxiliary basis while ours is
# exact, which is worth ~1.9e-5 here; forcing DFDH onto cc-pVDZ-RIFIT for J/K makes
# it *worse* by 1.7e-3, because RIFIT fits correlation and not J/K. So the residual
# is the reference's approximation, not ours, and tightening this bound would mean
# reproducing an approximation rather than an answer.
#
# What is checked tightly instead: the perturbative term on its own, where both
# codes fit and the only difference is which auxiliary basis. And the composition,
# in test_mqc_xc_spec, against the coefficients DFDH reports for itself.
DH_CASES = [("b2plyp", "B2PLYP", -76.3533921732, -0.0652456679),
            ("b2gp-plyp", "B2GPPLYP", -76.3369776292, None),
            ("mpw2plyp", "mPW2PLYP", -76.3528724926, None)]
DH_TOL = 5e-5
DH_PT2_TOL = 1e-5
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
    check_pyscf_grid_convention()
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

    # ---- self-consistent Kohn-Sham, against PySCF's own RKS -----------------
    #
    # A different kind of check from the two above: not a fixed density and a
    # functional, but the whole loop, with the potential in the Fock matrix
    # deciding the orbitals. Both codes build their own grid at the same level,
    # so this is the one comparison where the grids do *not* cancel -- and it
    # still agrees to 2e-11, because the grid tables are the same tables.
    from pyscf import dft
    for tag, pyscf_xc, basis, level, direct in KS_CASES:
        mode = "_direct" if direct else ""
        path = pathlib.Path(f"/tmp/mqc_ks_{tag}_{basis}_L{level}{mode}.txt")
        if not path.exists():
            print(f"  MISSING {path} -- run ./build/check_dft first")
            failures += 1
            continue
        v = path.read_text().split()
        ours, iters = float(v[0]), int(v[1])

        mol = gto.Mole()
        mol.atom = atoms
        mol.unit = "Angstrom"
        mol.basis = {s: bse_to_pyscf(basis, s) for s in symbols}
        mol.cart = molecule_form(basis, symbols) == CARTESIAN
        mol.verbose = 0
        mol.build()
        mf = dft.RKS(mol)
        mf.xc = pyscf_xc
        mf.grids.level = level
        mf.conv_tol = 1e-11
        theirs = mf.kernel()
        d = abs(theirs - ours)
        ok = d < KS_TOL
        print(f"  {'ok  ' if ok else 'FAIL'} KS {tag + mode:16s}/{basis} L{level}  "
              f"ours {ours:.10f}  PySCF {theirs:.10f}  diff {d:.2e} "
              f"(<= {KS_TOL:.0e})  {iters} iterations")
        if not ok:
            failures += 1

    # ---- unrestricted Kohn-Sham, against PySCF's own UKS ---------------------
    #
    # <S^2> as well as the energy. Two codes' UKS can land on different stationary
    # points, and then the energies differ for a reason that is not a bug in either
    # -- so a disagreement should be able to say which kind it is.
    for molecule, tag, pyscf_xc, basis, level in UKS_CASES:
        path = pathlib.Path(f"/tmp/mqc_uks_{molecule}_{tag}_{basis}_L{level}.txt")
        if not path.exists():
            print(f"  MISSING {path} -- run ./build/check_dft first")
            failures += 1
            continue
        v = path.read_text().split()
        ours, iters, s2 = float(v[0]), int(v[1]), float(v[2])

        case_atoms = CH3_ATOMS if molecule == "ch3" else atoms
        spin = CH3_SPIN if molecule == "ch3" else 0
        case_symbols = {a[0] for a in case_atoms}
        mol = gto.Mole()
        mol.atom = case_atoms
        mol.unit = "Angstrom"
        mol.basis = {t: bse_to_pyscf(basis, t) for t in case_symbols}
        mol.cart = molecule_form(basis, case_symbols) == CARTESIAN
        mol.spin = spin
        mol.verbose = 0
        mol.build()
        mf = dft.UKS(mol)
        mf.xc = pyscf_xc
        mf.grids.level = level
        mf.conv_tol = 1e-11
        theirs = mf.kernel()
        their_s2 = mf.spin_square()[0]
        d = abs(theirs - ours)
        ds2 = abs(their_s2 - s2)
        ok = d < KS_TOL and ds2 < UKS_S2_TOL
        print(f"  {'ok  ' if ok else 'FAIL'} UKS {molecule:5s} {tag:9s}/{basis} L{level}  "
              f"ours {ours:.10f}  PySCF {theirs:.10f}  diff {d:.2e} "
              f"(<= {KS_TOL:.0e})  <S^2> {s2:.6f} vs {their_s2:.6f}  {iters} iterations")
        if not ok:
            failures += 1

    # ---- double hybrids -----------------------------------------------------
    for tag, name, ref_total, ref_pt2 in DH_CASES:
        path = pathlib.Path(f"/tmp/mqc_dh_{tag}_cc-pvdz_ri.txt")
        if not path.exists():
            print(f"  MISSING {path} -- run ./build/check_dft first")
            failures += 1
            continue
        v = [float(x) for x in path.read_text().split()]
        total, ks, pt2 = v[0], v[1], v[2]
        d = abs(total - ref_total)
        ok = d < DH_TOL
        line = (f"  {'ok  ' if ok else 'FAIL'} DH {tag:10s} ours {total:.10f}  "
                f"DFDH {ref_total:.10f}  diff {d:.2e} (<= {DH_TOL:.0e})")
        if ref_pt2 is not None:
            dp2 = abs(pt2 - ref_pt2)
            ok = ok and dp2 < DH_PT2_TOL
            line += f"\n       PT2 {pt2:.10f} vs {ref_pt2:.10f}  diff {dp2:.2e}"
        print(line)
        if not ok:
            failures += 1

    # The fitting error in the perturbative term, measured against our own
    # conventional transform rather than against anyone else's auxiliary basis.
    ri = pathlib.Path("/tmp/mqc_dh_b2plyp_cc-pvdz_ri.txt")
    conv = pathlib.Path("/tmp/mqc_dh_b2plyp_cc-pvdz_conv.txt")
    if ri.exists() and conv.exists():
        a = [float(x) for x in ri.read_text().split()]
        b = [float(x) for x in conv.read_text().split()]
        # Same Kohn-Sham part by construction; only the PT2 term may move.
        assert abs(a[1] - b[1]) < 1e-12, "the KS part should not depend on the MP2 route"
        print(f"  ok   DH b2plyp RI vs conventional PT2: {abs(a[2] - b[2]):.2e} "
              f"(the fitting error, ours, measured internally)")

    print()
    if failures:
        print(f"[dft] {failures} case(s) disagree")
        return 1
    print("[dft] exchange-correlation energies and Kohn-Sham totals match PySCF")
    return 0


if __name__ == "__main__":
    sys.exit(main())
