#!/usr/bin/env python3
"""PySCF's spherically averaged free atoms, for test_mqc_czt_spherical_atom.

`pyscf.scf.atom_hf.AtomSphericAverageRHF` and `pyscf.scf.atom_ks`'s RKS twin
solve one atom spin-restricted with every angular momentum's electrons spread
evenly over its 2l+1 components, the configuration taken from
`pyscf.data.elements.NRSRHF_CONFIGURATION` (Lehtola, Phys. Rev. A 101, 012516
(2020)). The same solver underlies the Hirshfeld analysis of the PySCF
extension github.com/frobnitzem/hirshfeld, which is where the approach this
repository follows was taken from. Those are the "restricted" cases.

PySCF has no spin-polarised counterpart, so the "polarized" cases are its
UHF and UKS with the same two pieces replaced: `eig` averages each spin's Fock
matrix over components exactly as `atom_hf` does, and `get_occ` splits each l
by Hund's rule -- closed shells to both spins, the open shell to alpha first.

The orbital basis is read from this repository's own BSE JSON rather than
PySCF's tables, which differ in the last digits on some sets. The files used
are sorted by angular momentum, so PySCF's AO order is mqc's.

    python3 tools/spherical_atom/spherical_atom_reference.py
"""
import json
import pathlib

import numpy as np
from pyscf import gto
from pyscf.data import elements
from pyscf.dft import uks
from pyscf.scf import atom_hf, atom_ks, uhf

ROOT = pathlib.Path(__file__).resolve().parents[2]

# name: (symbol, basis, functional or None for Hartree-Fock, polarized)
ATOMS = {
    "h_ccpvdz_hf_restricted": ("H", "cc-pvdz", None, False),
    "c_ccpvdz_hf_restricted": ("C", "cc-pvdz", None, False),
    "o_ccpvdz_hf_restricted": ("O", "cc-pvdz", None, False),
    "ne_ccpvdz_hf_restricted": ("Ne", "cc-pvdz", None, False),
    "fe_def2svp_hf_restricted": ("Fe", "def2-svp", None, False),
    "o_ccpvdz_pbe_restricted": ("O", "cc-pvdz", "pbe", False),
    "h_ccpvdz_hf_polarized": ("H", "cc-pvdz", None, True),
    "c_ccpvdz_hf_polarized": ("C", "cc-pvdz", None, True),
    "o_ccpvdz_hf_polarized": ("O", "cc-pvdz", None, True),
    "fe_def2svp_hf_polarized": ("Fe", "def2-svp", None, True),
    "o_ccpvdz_pbe_polarized": ("O", "cc-pvdz", "pbe", True),
    "o_ccpvdz_b3lyp_polarized": ("O", "cc-pvdz", "b3lyp", True),
    "o_ccpvdz_wb97x_polarized": ("O", "cc-pvdz", "wb97x", True),
}


def spin_electrons(symbol):
    """Electrons per l in each spin, Hund's rule within each l."""
    alpha, beta = [], []
    for l, n_l in enumerate(elements.NRSRHF_CONFIGURATION[gto.charge(symbol)]):
        per_spin = 2 * l + 1
        n_double, n_open = divmod(n_l, 2 * per_spin)
        alpha.append(per_spin * n_double + min(n_open, per_spin))
        beta.append(per_spin * n_double + max(0, n_open - per_spin))
    return alpha, beta


def channel_occupation(mol, electrons):
    """One spin's occupation in atom_hf's orbital order: l, then column, then m."""
    occ = []
    for l in range(4):
        degen = 2 * l + 1
        nbas_l = mol._bas[mol._bas[:, gto.ANG_OF] == l, gto.NCTR_OF].sum()
        occ_l = np.zeros(nbas_l)
        per_component = electrons[l] / degen
        for k in range(nbas_l):
            if per_component <= 1e-12:
                break
            occ_l[k] = min(1.0, per_component)
            per_component -= occ_l[k]
        assert per_component <= 1e-12, (mol.atom_symbol(0), l, "basis too small")
        occ.append(np.repeat(occ_l, degen))
    return np.hstack(occ)


class _SphericalPolarized:
    """eig and get_occ for an unrestricted solver, spherically averaged."""

    def eig(self, fock, s, overwrite=False, x=None):
        ea, ca = atom_hf.AtomSphericAverageRHF.eig(self, fock[0], s)
        eb, cb = atom_hf.AtomSphericAverageRHF.eig(self, fock[1], s)
        return np.array((ea, eb)), np.array((ca, cb))

    def get_occ(self, mo_energy=None, mo_coeff=None):
        alpha, beta = spin_electrons(self.mol.atom_symbol(0))
        return np.array((channel_occupation(self.mol, alpha),
                         channel_occupation(self.mol, beta)))

    def get_grad(self, mo_coeff, mo_occ, fock=None):
        return commutator(self, mo_coeff, mo_occ, fock)


def commutator(mf, mo_coeff, mo_occ, fock):
    """F D S - S D F for each spin, standing in for the orbital gradient.

    atom_hf returns zero here, which leaves convergence to the energy alone and
    the density good only to the square root of conv_tol. The commutator
    vanishes exactly when D is built from F's own eigenvectors, so it can be
    held to conv_tol_grad instead.
    """
    dm = np.asarray(mf.make_rdm1(mo_coeff, mo_occ))
    if fock is None:
        fock = mf.get_fock(dm=dm)
    s = mf.get_ovlp()
    fock, dm = np.asarray(fock), np.asarray(dm)
    if dm.ndim == 2:
        fock, dm = fock[None], dm[None] / 2
    return np.hstack([(f @ d @ s - s @ d @ f).ravel() for f, d in zip(fock, dm)])


class AtomSphericalRHF(atom_hf.AtomSphericAverageRHF):
    get_grad = _SphericalPolarized.get_grad


class AtomSphericalRKS(atom_ks.AtomSphericAverageRKS):
    get_grad = _SphericalPolarized.get_grad


class AtomSphericalUHF(_SphericalPolarized, uhf.UHF):
    pass


class AtomSphericalUKS(_SphericalPolarized, uks.UKS):
    pass


def bse_to_pyscf(name, z):
    data = json.loads((ROOT / "basis_sets" / (name + ".json")).read_text())
    shells = []
    for sh in data["elements"][str(z)]["electron_shells"]:
        exps = [float(e) for e in sh["exponents"]]
        cols = [[float(c) for c in col] for col in sh["coefficients"]]
        moms = sh["angular_momentum"]
        if len(moms) == 1:
            shells.append([moms[0]] + [[e] + [c[i] for c in cols]
                                       for i, e in enumerate(exps)])
        else:
            for l, col in zip(moms, cols):
                shells.append([l] + [[e, col[i]] for i, e in enumerate(exps)])
    ls = [s[0] for s in shells]
    assert ls == sorted(ls), (name, z, "not l-sorted; AO order would differ")
    return shells


def solve(symbol, basis, functional, polarized):
    z = gto.charge(symbol)
    alpha, beta = spin_electrons(symbol)
    mol = gto.M(atom=[(symbol, (0.0, 0.0, 0.0))],
                basis={symbol: bse_to_pyscf(basis, z)},
                spin=(sum(alpha) - sum(beta)) if polarized else z % 2,
                cart=False, verbose=0)
    if polarized:
        mf = AtomSphericalUHF(mol) if functional is None else AtomSphericalUKS(mol)
    elif functional is None:
        mf = AtomSphericalRHF(mol)
    else:
        mf = AtomSphericalRKS(mol)
    if functional is not None:
        mf.xc = functional
        # Far finer than either code's default, so what is left between the
        # two is the quadrature scheme and not its resolution.
        mf.grids.atom_grid = (200, 1202)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-10
    mf.max_cycle = 200
    mf.run()
    assert mf.converged, (symbol, basis, functional)
    dm = np.asarray(mf.make_rdm1())
    if dm.ndim == 2:
        dm = np.array((dm / 2, dm / 2))
    s = mol.intor("int1e_ovlp")
    return mol, mf, dm, s


def main():
    for name, (symbol, basis, functional, polarized) in ATOMS.items():
        mol, mf, dm, s = solve(symbol, basis, functional, polarized)
        print(name, "nao", mol.nao)
        print("  energy      %.12f" % mf.e_tot)
        print("  tr(DaS)     %.12f" % np.einsum("ij,ji", dm[0], s))
        print("  tr(DbS)     %.12f" % np.einsum("ij,ji", dm[1], s))
        print("  Da(1,1)     %.12f" % dm[0][0, 0])
        print("  Db(1,1)     %.12f" % dm[1][0, 0])
        print("  |Da+Db|_F   %.12f" % np.linalg.norm(dm[0] + dm[1]))


if __name__ == "__main__":
    main()
