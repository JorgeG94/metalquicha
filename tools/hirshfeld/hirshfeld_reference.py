#!/usr/bin/env python3
"""Hirshfeld charges by this repository's recipe, in PySCF, for test_mqc_hirshfeld.

The recipe is mqc_czt_hirshfeld's: each element's free atom is the spherically
averaged, spin-polarised wB97X atom in def2-TZVP -- the polarized case of
tools/spherical_atom/spherical_atom_reference.py, imported from there so the two
cannot drift -- and a point of space is shared between the atoms in proportion
to their free densities there (Hirshfeld, Theor. Chim. Acta 44, 129 (1977)).
The approach follows the PySCF extension github.com/frobnitzem/hirshfeld.

Two things are done differently from mqc on purpose, so that agreement means
something:

* each free atom's density is evaluated exactly, from its own orbitals at every
  molecular grid point, where mqc interpolates a radial table;
* the molecular integral runs on PySCF's own grid, at a size well beyond mqc's
  default.

Nothing is renormalised: the charges sum to the molecular charge to the
accuracy of the grid, and that sum is printed.

    python3 tools/hirshfeld/hirshfeld_reference.py
"""
import pathlib
import sys

import numpy as np
from pyscf import dft, gto, scf

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "spherical_atom"))
from spherical_atom_reference import bse_to_pyscf, solve  # noqa: E402

FREE_ATOM_BASIS = "def2-tzvp"
FREE_ATOM_FUNCTIONAL = "wb97x"
MOLECULAR_GRID = (150, 974)

# Bohr: the water of test_mqc_fukui, C2v with the hydrogens exactly equivalent.
WATER = [("O", (0.0, 0.0, 0.0)),
         ("H", (0.0, 1.4308, 1.1078)),
         ("H", (0.0, -1.4308, 1.1078))]

_free_atoms = {}


def free_atom(symbol):
    """The converged polarized free atom, as (mol, total density matrix)."""
    if symbol not in _free_atoms:
        mol, mf, dm, _ = solve(symbol, FREE_ATOM_BASIS, FREE_ATOM_FUNCTIONAL, True)
        _free_atoms[symbol] = (mol, dm[0] + dm[1])
    return _free_atoms[symbol]


def molecule(atoms, basis, charge=0, spin=0):
    symbols = sorted({a[0] for a in atoms})
    return gto.M(atom=atoms, unit="Bohr", charge=charge, spin=spin,
                 basis={s: bse_to_pyscf(basis, gto.charge(s)) for s in symbols},
                 cart=False, verbose=0)


def converge(mol):
    mf = scf.RHF(mol) if mol.spin == 0 else scf.UHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.run()
    assert mf.converged
    dm = np.asarray(mf.make_rdm1())
    if dm.ndim == 2:
        return mf, dm, None
    return mf, dm[0] + dm[1], dm[0] - dm[1]


def populations(mol, total, spin=None):
    grids = dft.Grids(mol)
    grids.atom_grid = MOLECULAR_GRID
    grids.build()
    ni = dft.numint.NumInt()
    ao = ni.eval_ao(mol, grids.coords)
    rho = ni.eval_rho(mol, ao, total)
    spin_rho = None if spin is None else ni.eval_rho(mol, ao, spin)

    free = np.empty((mol.natm, len(grids.weights)))
    for a in range(mol.natm):
        amol, adm = free_atom(mol.atom_symbol(a))
        shifted = grids.coords - mol.atom_coord(a)
        free[a] = ni.eval_rho(amol, ni.eval_ao(amol, shifted), adm)
    total_free = free.sum(axis=0)
    share = np.divide(free, total_free, out=np.zeros_like(free), where=total_free > 0)

    electrons = share @ (grids.weights * rho)
    charges = mol.atom_charges() - electrons
    spins = None if spin_rho is None else share @ (grids.weights * spin_rho)
    return charges, spins


def show(label, values):
    print("  %-10s" % label + " ".join("%17.12f" % v for v in values))


def main():
    neutral = molecule(WATER, "cc-pvdz")
    _, d0, _ = converge(neutral)
    q0, _ = populations(neutral, d0)
    print("water_ccpvdz RHF")
    show("charges", q0)
    print("  sum        %.3e" % q0.sum())

    cation = molecule(WATER, "cc-pvdz", charge=1, spin=1)
    _, dp, sp = converge(cation)
    qp, spins = populations(cation, dp, sp)
    print("water_ccpvdz_cation UHF")
    show("charges", qp)
    show("spins", spins)
    print("  sum        %.3e" % (qp.sum() - 1.0))

    anion = molecule(WATER, "cc-pvdz", charge=-1, spin=1)
    _, dm, _ = converge(anion)
    qm, _ = populations(anion, dm)
    print("water_ccpvdz_fukui")
    show("f+", q0 - qm)
    show("f-", qp - q0)


if __name__ == "__main__":
    main()
