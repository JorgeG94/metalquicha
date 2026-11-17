#!/usr/bin/env python3
"""A reference RI-CCSD analytic nuclear gradient, in numpy over PySCF pieces.

PySCF has a conventional CCSD gradient and a density-fitted CCSD energy, and no
gradient for the second. This is that missing reference, written before the
Fortran for the reason `ri_mp2_gradient.py` gives: a factor resolved here costs
minutes.

What is differentiated, precisely: an *exact-ERI* restricted Hartree-Fock
reference with a *fitted* CCSD correlation energy, frozen core allowed. That is
what an `ri-ccsd` deck asks for. The Fock matrix in the correlation energy is
the reference's, so its derivative is never fitted; only the fluctuation
integrals are.

The amplitudes, the Lambda amplitudes and the densities are PySCF's
(`dfccsd`, `ccsd_lambda`, `ccsd_rdm`). What is written here is everything after
them, and it is written for a general non-variational correlation energy

    E_corr = sum_pq gamma_pq f_pq + 1/2 sum_pqrs Gamma_pqrs (pq|rs)_fit

over all molecular orbitals, frozen ones included, with `gamma` and `Gamma`
symmetrised. Nothing below knows the method is coupled cluster.

Orbital rotations. Under C -> C U the energy changes by 2 sum_tp U_tp X_tp with

    X_tp = (f gamma)_tp + sum_qrs (tq|rs) Gamma_pqrs + 2 V[gamma]_tp [p occupied]

and V[gamma] = J[gamma] - K[gamma]/2 from exact integrals. The orbitals are
held canonical by multipliers z on f_ai = 0 and, with a frozen core, on
f_ci = 0 between core and active occupied orbitals; z is chosen so that X built
from gamma + z is symmetric in those blocks. What is left of U is -S^x/2, and

    dE/dx = dE_HF/dx + sum gamma_rel f^(x) + 1/2 sum Gamma (pq|rs)_fit^(x)
            - sum S^(x) X_rel.

Run it directly to check against finite differences::

    python tools/cpu_validation/ri_ccsd_gradient.py
"""

from __future__ import annotations

import numpy as np
from pyscf import df, gto, scf
from pyscf.cc import ccsd_rdm, dfccsd


# --------------------------------------------------------------------------
# energy, so the finite difference differentiates exactly what the analytic
# gradient claims to
# --------------------------------------------------------------------------

def run_ri_ccsd(mol, auxbasis, frozen=0, with_lambda=False):
    """Exact-ERI RHF, then CCSD over fitted integrals."""
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-13
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("RHF did not converge")
    mycc = dfccsd.RCCSD(mf, frozen=frozen)
    mycc.with_df = df.DF(mol, auxbasis=auxbasis)
    mycc.conv_tol = 1e-12
    mycc.conv_tol_normt = 1e-10
    mycc.kernel()
    if not mycc.converged:
        raise RuntimeError("CCSD did not converge")
    if with_lambda:
        mycc.conv_tol_normt = 1e-10
        mycc.solve_lambda()
        if not mycc.converged_lambda:
            raise RuntimeError("the Lambda equations did not converge")
    return mf, mycc


def ri_ccsd_energy(mol, auxbasis, frozen=0):
    mf, mycc = run_ri_ccsd(mol, auxbasis, frozen)
    return mycc.e_tot


# --------------------------------------------------------------------------
# densities
# --------------------------------------------------------------------------

def correlation_densities(mycc):
    """gamma (n,n) and Gamma (n,n,n,n) over all MOs, symmetrised.

    The parts of the one- and two-particle densities that multiply the Fock
    matrix and the fluctuation integrals in the correlation energy: no
    reference density, and none of its products with gamma.
    """
    t1, t2, l1, l2 = mycc.t1, mycc.t2, mycc.l1, mycc.l2
    d1 = ccsd_rdm._gamma1_intermediates(mycc, t1, t2, l1, l2)
    d2 = ccsd_rdm._gamma2_intermediates(mycc, t1, t2, l1, l2)
    gamma = ccsd_rdm._make_rdm1(mycc, d1, with_frozen=True, with_mf=False)
    big = ccsd_rdm._make_rdm2(mycc, d1, d2, with_dm1=False, with_frozen=True)
    gamma = 0.5 * (gamma + gamma.T)
    big = 0.25 * (big + big.transpose(1, 0, 2, 3) + big.transpose(0, 1, 3, 2)
                  + big.transpose(1, 0, 3, 2))
    big = 0.5 * (big + big.transpose(2, 3, 0, 1))
    return gamma, big


# --------------------------------------------------------------------------
# the gradient
# --------------------------------------------------------------------------

def generalised_fock(fock, g_mo, gamma, big, v_gamma, n_occ):
    """X_tp, whose antisymmetric part is the orbital gradient."""
    x = fock @ gamma + np.einsum("tqrs,pqrs->tp", g_mo, big)
    x[:, :n_occ] += 2.0 * v_gamma[:, :n_occ]
    return x


def ri_ccsd_gradient(mol, auxbasis, frozen=0, checks=None):
    mf, mycc = run_ri_ccsd(mol, auxbasis, frozen, with_lambda=True)
    c = mf.mo_coeff
    n = c.shape[0]
    n_occ = mol.nelectron // 2
    n_core = frozen
    eps = mf.mo_energy
    fock = np.diag(eps)
    dm_hf = mf.make_rdm1()

    auxmol = df.addons.make_auxmol(mol, auxbasis)
    n_aux = auxmol.nao
    three = df.incore.aux_e2(mol, auxmol, "int3c2e")          # (mu nu|P)
    metric = auxmol.intor("int2c2e")
    fit = np.linalg.solve(metric, three.reshape(n * n, n_aux).T)  # c^P_{mu nu}
    fit = fit.reshape(n_aux, n, n)
    three_mo = np.einsum("mnP,mp,nq->pqP", three, c, c)
    fit_mo = np.einsum("Pmn,mp,nq->Ppq", fit, c, c)
    g_mo = np.einsum("pqP,Prs->pqrs", three_mo, fit_mo)

    gamma, big = correlation_densities(mycc)

    def v_of(gam_mo):
        """J - K/2 of an MO-basis symmetric matrix, back in the MO basis."""
        p_ao = c @ gam_mo @ c.T
        vj, vk = mf.get_jk(mol, p_ao, hermi=1)
        return c.T @ (vj - 0.5 * vk) @ c

    # The correlation energy from these densities must be the CCSD one, or
    # nothing after this line differentiates the right thing.
    e_corr = np.einsum("pq,pq", gamma, fock) + 0.5 * np.einsum("pqrs,pqrs", big, g_mo)
    if checks is not None:
        checks["energy"] = abs(e_corr - mycc.e_corr)

    # ----- the multipliers that keep the orbitals canonical -----
    pairs = [(a, i) for i in range(n_occ) for a in range(n_occ, n)]
    pairs += [(i, k) for k in range(n_core) for i in range(n_core, n_occ)]

    def z_matrix(z):
        zs = np.zeros((n, n))
        for value, (p, q) in zip(z, pairs):
            zs[p, q] += 0.5 * value
            zs[q, p] += 0.5 * value
        return zs

    def residual(x):
        return np.array([x[p, q] - x[q, p] for p, q in pairs])

    x0 = generalised_fock(fock, g_mo, gamma, big, v_of(gamma), n_occ)
    if checks is not None:
        # Rotations among active occupied, and among virtual, orbitals leave
        # the energy alone at convergence.
        act = slice(n_core, n_occ)
        vir = slice(n_occ, n)
        checks["occ_occ"] = abs(x0[act, act] - x0[act, act].T).max()
        checks["vir_vir"] = abs(x0[vir, vir] - x0[vir, vir].T).max()
    rhs = residual(x0)
    n_z = len(pairs)
    op = np.zeros((n_z, n_z))
    zero_big = np.zeros_like(big)
    for k in range(n_z):
        unit = np.zeros(n_z)
        unit[k] = 1.0
        zs = z_matrix(unit)
        op[:, k] = residual(generalised_fock(fock, g_mo, zs, zero_big, v_of(zs), n_occ))
    z = np.linalg.solve(op, -rhs)
    gamma_rel = gamma + z_matrix(z)
    x_rel = generalised_fock(fock, g_mo, gamma_rel, big, v_of(gamma_rel), n_occ)
    if checks is not None:
        checks["stationary"] = abs(residual(x_rel)).max()
    x_sym = 0.5 * (x_rel + x_rel.T)

    # ----- back to the atomic orbitals -----
    p_ao = c @ gamma_rel @ c.T
    w_ao = c @ x_sym @ c.T
    big_ao = np.einsum("pqrs,mp,nq,lr,ks->mnlk", big, c, c, c, c, optimize=True)
    gamma3 = np.einsum("mnlk,Plk->Pmn", big_ao, fit)            # Gamma^P_{mu nu}
    gamma2 = np.einsum("Pmn,Qmn->PQ", gamma3, fit)
    gamma2 = 0.5 * (gamma2 + gamma2.T)
    fock_density = (np.einsum("mn,lk->mnlk", p_ao, dm_hf)
                    - 0.5 * np.einsum("ml,nk->mnlk", p_ao, dm_hf))

    # ----- derivative integrals -----
    grad_hf = mf.nuc_grad_method()
    grad_hf.verbose = 0
    gradient = grad_hf.kernel().copy()
    hcore_deriv = grad_hf.hcore_generator(mol)
    s1 = mol.intor("int1e_ipovlp")                              # (nabla mu|nu)
    e1 = mol.intor("int2e_ip1").reshape(3, n, n, n, n)          # (nabla mu nu|l k)
    ip1 = df.incore.aux_e2(mol, auxmol, "int3c2e_ip1")          # (nabla mu nu|P)
    ip2 = df.incore.aux_e2(mol, auxmol, "int3c2e_ip2")          # (mu nu|nabla P)
    v1 = auxmol.intor("int2c2e_ip1")                            # (nabla P|Q)
    ao_slices = mol.aoslice_by_atom()
    aux_slices = auxmol.aoslice_by_atom()

    for atom in range(mol.natm):
        p0, p1 = ao_slices[atom, 2:]
        q0, q1 = aux_slices[atom, 2:]
        g = np.einsum("xmn,mn->x", hcore_deriv(atom), p_ao)
        # The exact four-centre derivative, one slot at a time.
        g -= np.einsum("xijkl,ijkl->x", e1[:, p0:p1], fock_density[p0:p1])
        g -= np.einsum("xjikl,ijkl->x", e1[:, p0:p1], fock_density[:, p0:p1])
        g -= np.einsum("xklij,ijkl->x", e1[:, p0:p1], fock_density[:, :, p0:p1])
        g -= np.einsum("xlkij,ijkl->x", e1[:, p0:p1], fock_density[:, :, :, p0:p1])
        # The fitted fluctuation integrals.
        g -= 2.0 * np.einsum("xmnP,Pmn->x", ip1[:, p0:p1], gamma3[:, p0:p1])
        g -= np.einsum("xmnP,Pmn->x", ip2[:, :, :, q0:q1], gamma3[q0:q1])
        g += np.einsum("xPQ,PQ->x", v1[:, q0:q1], gamma2[q0:q1])
        # Orthonormality.
        g += 2.0 * np.einsum("xmn,mn->x", s1[:, p0:p1], w_ao[p0:p1])
        gradient[atom] += g
    return gradient


def finite_difference(mol, auxbasis, frozen, atom, axis, step=1.0e-3):
    coords = mol.atom_coords()
    values = []
    for sign in (1.0, -1.0):
        moved = coords.copy()
        moved[atom, axis] += sign * step
        shifted = mol.copy()
        shifted.set_geom_(moved, unit="Bohr")
        values.append(ri_ccsd_energy(shifted, auxbasis, frozen))
    return (values[0] - values[1]) / (2.0 * step)


def main():
    cases = [
        ("sto-3g", "cc-pvdz-ri", 0),
        ("sto-3g", "cc-pvdz-ri", 1),
        ("cc-pvdz", "cc-pvdz-ri", 1),
    ]
    for basis, auxbasis, frozen in cases:
        mol = gto.M(atom="O 0 0 0.117; H 0 0.757 -0.469; H 0.1 -0.757 -0.469",
                    basis=basis, verbose=0)
        checks = {}
        gradient = ri_ccsd_gradient(mol, auxbasis, frozen, checks)
        print(f"water {basis} / {auxbasis}, {frozen} frozen")
        for name, value in checks.items():
            print(f"  {name:12s} {value:.2e}")
        print(f"  sum over atoms {abs(gradient.sum(axis=0)).max():.2e}")
        worst = 0.0
        for atom, axis in ((0, 2), (1, 1), (2, 0)):
            fd = finite_difference(mol, auxbasis, frozen, atom, axis)
            worst = max(worst, abs(fd - gradient[atom, axis]))
            print(f"  atom {atom} axis {axis}: analytic {gradient[atom, axis]: .9f}"
                  f"  finite difference {fd: .9f}")
        print(f"  worst difference {worst:.2e}")


if __name__ == "__main__":
    main()
