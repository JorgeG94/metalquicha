"""Mayer bond orders from PySCF, as the reference the Fortran tests are pinned to.

Run it, and it prints the matrices that ``test/test_mqc_mayer_bond_orders.f90``
carries as literals, plus the geometries those tests build.

Two things this script exists to get right, both of which have bitten us:

**The basis comes out of this repository's own ``basis_sets/``**, converted with
``bse_to_pyscf`` from the validation generator, not from PySCF's internal
tables. The two differ in the eighth decimal of the exponents on Pople sets,
which fakes a disagreement that looks exactly like a bug in whichever code you
are checking.

**The open-shell formula is not the closed-shell one with the total density.**
It is

    B_AB = 2 sum_{mu in A, nu in B} [ (Da S)_mu,nu (Da S)_nu,mu
                                    + (Db S)_mu,nu (Db S)_nu,mu ]

and it reduces to the closed-shell expression only when Da = Db. So the
open-shell case below prints *both* numbers: the right one, and the one the
closed-shell expression gives on the same density. The Fortran test is pinned
tightly enough to tell them apart, and this is where you check that the gap is
real before trusting that.

Usage::

    source ~/dev/mqc_worktrees/mqc_env.sh
    python tools/cpu_validation/mayer_bond_orders_reference.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gen_cpu_validation import bse_to_pyscf  # noqa: E402

BOHR = 0.52917721092

#: Triplet O2 at its experimental bond length, in Angstrom. The open-shell
#: case the tests are pinned to, because it is where the two formulas are
#: furthest apart: the open-shell one gives exactly 2.00, the closed-shell one
#: applied to the same total density gives 1.50.
OXYGEN = [("O", (0.0, 0.0, 0.0)), ("O", (0.0, 0.0, 1.208))]

#: Water, the geometry `test_mqc_czt_charges.f90` already uses, in Angstrom.
WATER = [
    ("O", (0.0, 0.0, 0.0)),
    ("H", (0.0, 0.0, 0.9584)),
    ("H", (0.9268, 0.0, -0.2400)),
]


def ethane():
    """Staggered ethane: C-C 1.526 A, C-H 1.088 A, H-C-C 111.2 degrees.

    Built rather than tabulated so the Fortran test and this script cannot
    drift apart in the fourth decimal -- the numbers this prints are what the
    test carries.
    """
    r_cc, r_ch, theta = 1.526, 1.088, np.radians(111.2)
    atoms = [("C", (0.0, 0.0, 0.0)), ("C", (0.0, 0.0, r_cc))]
    s, c = np.sin(theta), np.cos(theta)
    for k in range(3):
        phi = 2.0 * np.pi * k / 3.0
        atoms.append(
            ("H", (r_ch * s * np.cos(phi), r_ch * s * np.sin(phi), r_ch * c))
        )
    for k in range(3):
        phi = 2.0 * np.pi * k / 3.0 + np.pi / 3.0
        atoms.append(
            (
                "H",
                (
                    r_ch * s * np.cos(phi),
                    r_ch * s * np.sin(phi),
                    r_cc - r_ch * c,
                ),
            )
        )
    return atoms


def build(atoms, basis, charge=0, spin=0):
    from pyscf import gto

    mol = gto.Mole()
    mol.atom = [(sym, tuple(xyz)) for sym, xyz in atoms]
    mol.unit = "Angstrom"
    mol.charge = charge
    mol.spin = spin
    mol.basis = {sym: bse_to_pyscf(basis, sym) for sym, _ in atoms}
    mol.build()
    return mol


def mayer(mol, densities):
    """B_AB from a list of density matrices: one entry closed shell, two open.

    The closed-shell entry is the total density and carries its own factor of
    two; the open-shell pair is alpha and beta, each squared and the sum
    doubled.
    """
    s = mol.intor("int1e_ovlp")
    slices = mol.aoslice_by_atom()
    natm = mol.natm
    b = np.zeros((natm, natm))
    factor = 1.0 if len(densities) == 1 else 2.0
    for d in densities:
        ds = d @ s
        for a in range(natm):
            ia, ja = slices[a][2], slices[a][3]
            for c in range(natm):
                ic, jc = slices[c][2], slices[c][3]
                b[a, c] += factor * np.sum(ds[ia:ja, ic:jc] * ds[ic:jc, ia:ja].T)
    np.fill_diagonal(b, 0.0)
    return b


def report(name, mol, b, note=""):
    print(f"\n=== {name} {note}")
    labels = [f"{mol.atom_symbol(i)}{i + 1}" for i in range(mol.natm)]
    for i in range(mol.natm):
        for j in range(i + 1, mol.natm):
            if b[i, j] < 1.0e-4:
                continue
            print(f"  {labels[i]:>4} -- {labels[j]:<4} {b[i, j]: .12f}")
    print("  valences:")
    for i in range(mol.natm):
        print(f"  {labels[i]:>4} {b[i].sum(): .12f}")


def closed_shell(atoms, basis, name, charge=0):
    from pyscf import scf

    mol = build(atoms, basis, charge=charge)
    mf = scf.RHF(mol)
    mf.conv_tol = 1.0e-12
    mf.kernel()
    b = mayer(mol, [mf.make_rdm1()])
    report(name, mol, b)
    return b


def open_shell(atoms, basis, name, charge, spin):
    from pyscf import scf

    mol = build(atoms, basis, charge=charge, spin=spin)
    mf = scf.UHF(mol)
    mf.conv_tol = 1.0e-12
    mf.kernel()
    da, db = mf.make_rdm1()
    right = mayer(mol, [da, db])
    wrong = mayer(mol, [da + db])
    report(name, mol, right, "(open-shell formula, the reference)")
    report(name, mol, wrong, "(closed-shell formula on the total density -- WRONG)")
    gap = np.abs(right - wrong).max()
    print(f"\n  largest element the wrong formula moves: {gap:.6f}")
    return right


if __name__ == "__main__":
    print("geometry, ethane, Angstrom:")
    for sym, xyz in ethane():
        print(f"  {sym} {xyz[0]: .6f} {xyz[1]: .6f} {xyz[2]: .6f}")
    closed_shell(WATER, "sto-3g", "water sto-3g")
    closed_shell(WATER, "6-31g", "water 6-31g")
    closed_shell(ethane(), "sto-3g", "ethane sto-3g")
    closed_shell(ethane(), "6-31g", "ethane 6-31g")
    open_shell(OXYGEN, "sto-3g", "triplet O2 sto-3g", charge=0, spin=2)
    open_shell(WATER, "sto-3g", "water cation doublet sto-3g", charge=1, spin=1)
