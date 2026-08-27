#!/usr/bin/env python3
"""A reference RI-MP2 analytic nuclear gradient, in numpy over PySCF integrals.

PySCF has no RI-MP2 gradient -- both `pyscf.mp.dfmp2.DFMP2.nuc_grad_method` and
`dfmp2_native`'s raise NotImplementedError -- so there is nothing to check a
Fortran implementation against. This is that missing reference, written first
because a factor resolved here costs minutes and the same factor resolved in
Fortran costs an afternoon.

The formulation follows Stocks, Palethorpe and Barca, J. Chem. Theory Comput.
2024, 20, 2505 (doi:10.1021/acs.jctc.3c01424), section 3; equation numbers below
are theirs. See ``rimp2grad.md`` in the repository root.

What is differentiated, precisely: an *exact-ERI* restricted Hartree-Fock
reference with a *fitted* correlation energy. That is what our RI-MP2 decks ask
for, and it is why the four-centre coupling tensor and the Fock derivative below
are never fitted -- the paper makes the same point for the same reason.

Run it directly to check against finite differences::

    python tools/cpu_validation/ri_mp2_gradient.py
"""

from __future__ import annotations

import numpy as np
from pyscf import df, gto, lib, scf


# --------------------------------------------------------------------------
# energy, so the finite difference differentiates exactly what the analytic
# gradient claims to
# --------------------------------------------------------------------------

def ri_mp2_energy(mol, auxmol, conv_tol=1e-13, frozen=0):
    """Exact-ERI RHF plus a fitted MP2 correlation energy."""
    mf = scf.RHF(mol)
    mf.conv_tol = conv_tol
    mf.kernel()
    if not mf.converged:
        raise SystemExit("RHF did not converge")
    b, _, _ = fitted_tensor(mol, auxmol, mf.mo_coeff, mol.nelectron // 2,
                            frozen=frozen)
    x, ovov = amplitudes(b, mf.mo_energy, mol.nelectron // 2, frozen=frozen)
    e_corr = np.einsum("ijab,ijab->", x, 2 * ovov - ovov.transpose(0, 1, 3, 2))
    return mf.e_tot + e_corr, mf


def fitted_tensor(mol, auxmol, mo_coeff, nocc, frozen=0):
    """B^P_ia (eq 5), together with the metric and its inverse square root.

    The inverse square root is returned rather than rebuilt by the caller: a
    different eigenvalue threshold would fit a different energy, and the
    gradient has to differentiate the one that was computed.
    """
    nao = mol.nao
    orbo = mo_coeff[:, frozen:nocc]
    orbv = mo_coeff[:, nocc:]

    three = df.incore.aux_e2(mol, auxmol, intor="int3c2e", aosym="s1")
    three = three.reshape(nao, nao, -1)
    metric = auxmol.intor("int2c2e")

    # (ia|Q)
    ia_q = np.einsum("ui,uvQ,va->iaQ", orbo, three, orbv, optimize=True)
    jm12 = metric_inverse_sqrt(metric)
    b = np.einsum("iaQ,QP->Pia", ia_q, jm12, optimize=True)
    return b, metric, jm12


def metric_inverse_sqrt(metric, threshold=1e-10):
    """J^(-1/2) through the eigendecomposition, dropping near-null modes."""
    w, v = np.linalg.eigh(metric)
    keep = w > threshold
    return (v[:, keep] / np.sqrt(w[keep])) @ v[:, keep].T


def amplitudes(b, mo_energy, nocc, frozen=0):
    """X^ab_ij (eq 9) and the fitted integrals it came from (eq 4)."""
    e_o = mo_energy[frozen:nocc]
    e_v = mo_energy[nocc:]
    ovov = np.einsum("Pia,Pjb->iajb", b, b, optimize=True)
    ovov = ovov.transpose(0, 2, 1, 3)          # (i, j, a, b)
    denom = (e_o[:, None, None, None] + e_o[None, :, None, None]
             - e_v[None, None, :, None] - e_v[None, None, None, :])
    return ovov / denom, ovov


# --------------------------------------------------------------------------
# the gradient
# --------------------------------------------------------------------------

def ri_mp2_gradient(mol, auxmol, mf, verbose=False, frozen=0):
    """dE/dR for exact RHF plus fitted MP2, in Hartree/Bohr.

    **The assembly is not a fresh expansion of eq 7.** Everything one-particle
    here -- the relaxed density, the energy-weighted density, and the
    core-Hamiltonian, overlap and Hartree-Fock two-electron derivative terms --
    is the conventional MP2 gradient's, in the arrangement `pyscf.grad.mp2`
    uses and our own Fortran reproduces to 1e-9. Only two things are replaced:

    * the two-particle density contracted with four-centre derivatives becomes
      the three- and two-centre terms of eq 6, and
    * the Lagrangian `Imat` is built from `L'` and `L''` (eqs 13, 14) instead of
      from a four-index density.

    Re-deriving the rest was tried and produced a term-3 error of about ten
    times the answer, invisible to translational invariance. The reference
    gradient lives *inside* this block -- `D_SCF` appears in the densities --
    so the Hartree-Fock gradient is not added separately.
    """
    nocc = mol.nelectron // 2
    nao, nmo = mf.mo_coeff.shape
    nvir = nmo - nocc
    mo, e_mo = mf.mo_coeff, mf.mo_energy
    orbo, orbv = mo[:, :nocc], mo[:, nocc:]

    b, metric, jm12 = fitted_tensor(mol, auxmol, mo, nocc, frozen=frozen)
    x, _ = amplitudes(b, e_mo, nocc, frozen=frozen)

    # ---- three- and two-index densities (eqs 8, 10) ----------------------
    xbar = 2 * x - x.transpose(0, 1, 3, 2)             # 2X^ab_ij - X^ba_ij
    gamma = np.einsum("ijab,Qjb,PQ->Pia", xbar, b, jm12, optimize=True)
    gamma_pq = np.einsum("Pia,Ria,QR->PQ", gamma, b, jm12, optimize=True)

    three = df.incore.aux_e2(mol, auxmol, intor="int3c2e", aosym="s1")
    three = three.reshape(nao, nao, -1)

    # Two exact identities, and the cheapest gates in the whole routine. Both
    # follow from sum_P (ia|P) J^{-1/2}_PQ = B^Q_ia:
    #
    #   sum_{P,ia} Gamma^P_ia (ia|P) = E_corr,   sum_PQ gamma_PQ J_PQ = E_corr
    #
    # A factor slipped between the paper's `2X - X^swap` and the conventional
    # `4X - 2X^swap` shows up here and in no later quantity.
    if verbose:
        ia_p = np.einsum("ui,uvP,va->iaP", mo[:, frozen:nocc], three, orbv,
                         optimize=True)
        print("  E_corr from Gamma^P (ia|P): %.12f   from gamma_PQ J_PQ: %.12f"
              % (np.einsum("Pia,iaP->", gamma, ia_p),
                 np.einsum("PQ,PQ->", gamma_pq, metric)))

    # ---- the unrelaxed one-particle blocks (eqs 17, 18) ------------------
    #
    # In the form our conventional MP2 gradient uses, checked element-wise
    # against pyscf.mp.mp2._gamma1_intermediates; eq 17 scanned ambiguously and
    # this is the same object.
    doo = -(2 * np.einsum("kiab,kjab->ij", x, x, optimize=True)
            - np.einsum("kiab,kjba->ij", x, x, optimize=True))
    dvv = (2 * np.einsum("ijca,ijcb->ba", x, x, optimize=True)
           - np.einsum("ijca,ijbc->ba", x, x, optimize=True))

    dm1mo = np.zeros((nmo, nmo))
    dm1mo[frozen:nocc, frozen:nocc] = doo + doo.T
    dm1mo[nocc:, nocc:] = dvv + dvv.T

    # ---- the Lagrangian, from the fitted integrals (eqs 13, 14) ----------
    #
    # `(pq|P)` over all molecular orbitals, which is why the three-centre
    # integrals are transformed in full rather than only to `(ia|P)`.
    #
    # The sign is what makes this drop into the conventional assembly's slot:
    # there the occupied-virtual blocks of `Imat` enter the Z-vector
    # right-hand side as `Imat[i,a]^T - Imat[a,i]`, and eq 21 wants
    # `L''_ai - L'_ia`, so `Imat = -L`.
    pq_p = np.einsum("up,uvP,vq->pqP", mo, three, mo, optimize=True)
    lag_o = 2 * np.einsum("Pia,qaP->iq", gamma, pq_p[:, nocc:, :], optimize=True)
    lag_v = 2 * np.einsum("Pia,iqP->aq", gamma, pq_p[frozen:nocc, :, :],
                          optimize=True)

    # Transposed into the slot. `L'_iq` and `L''_aq` are indexed with the
    # Lagrangian's own index first; the conventional assembly's `Imat` carries
    # it second. Verified numerically against a saturated-auxiliary
    # `pyscf.grad.mp2` Lagrangian: the two agree to 3.8e-7, which is the
    # fitting error, and disagree by 0.012 without the transpose -- with
    # matching diagonals and matching norms either way, so neither of those is
    # evidence.
    imat = np.zeros((nmo, nmo))
    imat[:, frozen:nocc] = -lag_o.T
    imat[:, nocc:] = -lag_v.T

    # The occupied-frozen block of the relaxed density, resolved directly from
    # the Lagrangian because both orbitals are occupied -- and before the veff
    # build, so the Z-vector's right-hand side sees the whole unrelaxed
    # density. pyscf.grad.mp2's dco, in its arrangement.
    if frozen:
        gap = e_mo[:frozen, None] - e_mo[None, frozen:nocc]
        dm1mo[:frozen, frozen:nocc] = imat[:frozen, frozen:nocc] / gap
        dm1mo[frozen:nocc, :frozen] = dm1mo[:frozen, frozen:nocc].T

    # ---- the Z-vector ----------------------------------------------------
    hf_dm = mf.make_rdm1()
    dm1_ao = mo @ dm1mo @ mo.T
    vhf = veff(mf, mol, dm1_ao) * 2
    xvo = orbv.T @ vhf @ orbo + imat[:nocc, nocc:].T - imat[nocc:, :nocc]

    dvo = solve_z_vector(mol, mf, xvo, nocc, nvir)
    dm1mo[nocc:, :nocc] += dvo
    dm1mo[:nocc, nocc:] += dvo.T

    imat[nocc:, :nocc] = imat[:nocc, nocc:].T
    im1 = mo @ imat @ mo.T

    # ---- the energy-weighted density -------------------------------------
    zeta = (e_mo[:, None] + e_mo[None, :]) * 0.5
    zeta[nocc:, :nocc] = e_mo[:nocc]
    zeta[:nocc, nocc:] = e_mo[:nocc][:, None]
    zeta = mo @ (zeta * dm1mo) @ mo.T

    dm1_ao = mo @ dm1mo @ mo.T
    p_occ = orbo @ orbo.T
    vhf_s1occ = p_occ @ veff(mf, mol, dm1_ao + dm1_ao.T) @ p_occ

    dm1p = hf_dm + dm1_ao * 2
    dm1_total = hf_dm + dm1_ao
    zeta += (orbo * e_mo[:nocc]) @ orbo.T * 2.0

    # ---- assembly --------------------------------------------------------
    grad = mf.nuc_grad_method()
    hcore_deriv = grad.hcore_generator(mol)
    s1 = grad.get_ovlp(mol)
    vhf1 = hf_two_electron_deriv(mol, hf_dm)
    aoslices = mol.aoslice_by_atom()

    de = np.zeros((mol.natm, 3))
    de += gradient_three_centre(mol, auxmol, gamma, mo, nocc, frozen=frozen)
    de += gradient_two_centre(auxmol, gamma_pq, mol)
    for a in range(mol.natm):
        p0, p1 = aoslices[a, 2], aoslices[a, 3]
        de[a] += np.einsum("xij,ij->x", s1[:, p0:p1], im1[p0:p1])
        de[a] += np.einsum("xij,ij->x", s1[:, p0:p1], im1.T[p0:p1])
        de[a] += np.einsum("xij,ji->x", hcore_deriv(a), dm1_total)
        de[a] -= np.einsum("xij,ij->x", s1[:, p0:p1], zeta[p0:p1])
        de[a] -= np.einsum("xij,ij->x", s1[:, p0:p1], zeta.T[p0:p1])
        de[a] -= 2 * np.einsum("xij,ij->x", s1[:, p0:p1], vhf_s1occ[p0:p1])
        de[a] -= np.einsum("xij,ij->x", vhf1[a], dm1p)
    de += grad.grad_nuc()

    if verbose:
        print("  |gamma| %.4e  |z| %.4e  |D| %.4e  |Imat| %.4e"
              % (abs(gamma).sum(), abs(dvo).sum(), abs(dm1_ao).sum(),
                 abs(imat).sum()))
    return de


def veff(mf, mol, dm):
    """J - K/2 from a symmetric density, the SCF's own combination."""
    vj, vk = mf.get_jk(mol, dm)
    return vj - vk * 0.5


def hf_two_electron_deriv(mol, hf_dm):
    """Per-atom (3, nao, nao) matrices the reference density contracts with.

    The four accumulations are the four positions an atomic orbital on the
    differentiated centre can take in a quartet; `int2e_ip1` differentiates the
    first index only, so the other three are reached by the integral's
    permutational symmetry.
    """
    nao = mol.nao
    eri1 = mol.intor("int2e_ip1").reshape(3, nao, nao, nao, nao)
    aoslices = mol.aoslice_by_atom()
    out = np.zeros((mol.natm, 3, nao, nao))
    for a in range(mol.natm):
        p0, p1 = aoslices[a, 2], aoslices[a, 3]
        blk = eri1[:, p0:p1]
        out[a] += np.einsum("xijkl,ij->xkl", blk, hf_dm[p0:p1])
        out[a] -= np.einsum("xijkl,il->xkj", blk, hf_dm[p0:p1]) * 0.5
        out[a][:, p0:p1] += np.einsum("xijkl,kl->xij", blk, hf_dm)
        out[a][:, p0:p1] -= np.einsum("xijkl,jk->xil", blk, hf_dm) * 0.5
    return out


def solve_z_vector(mol, mf, xvo, nocc, nvir):
    """Solve (eps_a - eps_i) D_ai + sum_bj A_aibj D_bj = L_ai, densely.

    Densely on purpose. PySCF's Krylov CPHF carries error its `tol` does not
    control -- it is what put 4e-8 between our conventional MP2 gradient and
    `pyscf.grad.mp2` -- and at reference sizes the operator is small enough to
    build column by column and hand to LAPACK.
    """
    mo = mf.mo_coeff
    orbo, orbv = mo[:, :nocc], mo[:, nocc:]
    e_mo = mf.mo_energy
    gaps = (e_mo[nocc:, None] - e_mo[None, :nocc]).ravel()

    n = nvir * nocc

    # The response operator in the *conventional* assembly's normalisation:
    # `2 C_v^T (J - K/2) C_o`, not the paper's `A`, which is twice it. Mixing
    # the two doubles the coupling and leaves the answer wrong by a few percent
    # in a way translational invariance cannot see.
    def apply(vec):
        d = orbv @ vec.reshape(nvir, nocc) @ orbo.T
        vj, vk = mf.get_jk(mol, d + d.T)
        return 2 * (orbv.T @ (vj - vk * 0.5) @ orbo).ravel()

    operator = np.zeros((n, n))
    basis = np.zeros(n)
    for k in range(n):
        basis[:] = 0.0
        basis[k] = 1.0
        operator[:, k] = apply(basis)
    operator += np.diag(gaps)
    # `A z = -X`, the sign the conventional assembly's right-hand side carries.
    return np.linalg.solve(operator, -xvo.ravel()).reshape(nvir, nocc)


# --------------------------------------------------------------------------
# the four gradient terms of eq 6
# --------------------------------------------------------------------------

def gradient_three_centre(mol, auxmol, gamma, mo, nocc, frozen=0):
    """4 sum_{munuP} Gamma^P_munu (P|munu)^xi -- the first term of eq 6."""
    nao = mol.nao
    orbo, orbv = mo[:, frozen:nocc], mo[:, nocc:]
    gamma_ao = np.einsum("Pia,ui,va->Puv", gamma, orbo, orbv, optimize=True)

    de = np.zeros((mol.natm, 3))
    aoslices = mol.aoslice_by_atom()
    auxslices = auxmol.aoslice_by_atom()

    # d/dR of the orbital centres, through int3c2e_ip1. (mu nu|P) is symmetric
    # in mu and nu, so the nu derivative is the mu one transposed.
    ip1 = df.incore.aux_e2(mol, auxmol, intor="int3c2e_ip1", aosym="s1", comp=3)
    ip1 = ip1.reshape(3, nao, nao, -1)
    for a in range(mol.natm):
        p0, p1 = aoslices[a, 2], aoslices[a, 3]
        de[a] -= 4 * np.einsum("xuvP,Puv->x", ip1[:, p0:p1], gamma_ao[:, p0:p1, :])
        de[a] -= 4 * np.einsum("xuvP,Pvu->x", ip1[:, p0:p1], gamma_ao[:, :, p0:p1])

    # d/dR of the auxiliary centre, through int3c2e_ip2.
    ip2 = df.incore.aux_e2(mol, auxmol, intor="int3c2e_ip2", aosym="s1", comp=3)
    ip2 = ip2.reshape(3, nao, nao, -1)
    for a in range(mol.natm):
        q0, q1 = auxslices[a, 2], auxslices[a, 3]
        de[a] -= 4 * np.einsum("xuvP,Puv->x", ip2[:, :, :, q0:q1],
                               gamma_ao[q0:q1])
    return de


def gradient_two_centre(auxmol, gamma_pq, mol):
    """-2 sum_PQ gamma_PQ (P|Q)^xi -- the second term of eq 6."""
    de = np.zeros((mol.natm, 3))
    ip1 = auxmol.intor("int2c2e_ip1", comp=3)
    auxslices = auxmol.aoslice_by_atom()
    for a in range(mol.natm):
        q0, q1 = auxslices[a, 2], auxslices[a, 3]
        de[a] += 2 * np.einsum("xPQ,PQ->x", ip1[:, q0:q1], gamma_pq[q0:q1])
        de[a] += 2 * np.einsum("xQP,PQ->x", ip1[:, q0:q1], gamma_pq[:, q0:q1])
    return de


# --------------------------------------------------------------------------
# checking
# --------------------------------------------------------------------------

def finite_difference(atoms, basis, auxbasis, step=2.5e-3, frozen=0):
    """Central differences of the same energy the analytic gradient claims."""
    de = np.zeros((len(atoms), 3))
    for a in range(len(atoms)):
        for c in range(3):
            plus = [list(x) for x in atoms]
            minus = [list(x) for x in atoms]
            plus[a][1 + c] += step
            minus[a][1 + c] -= step
            e_plus = ri_mp2_energy(*build(plus, basis, auxbasis),
                                   frozen=frozen)[0]
            e_minus = ri_mp2_energy(*build(minus, basis, auxbasis),
                                    frozen=frozen)[0]
            de[a, c] = (e_plus - e_minus) / (2 * step)
    return de


def build(atoms, basis, auxbasis, unit="Bohr"):
    """The two molecules, both bases read out of ``basis_sets/``.

    Not PySCF's own tables. They and this repository's JSON differ around the
    eighth decimal of an exponent, which puts about 3e-8 between the two codes'
    energies -- comfortably more than the difference between two correct
    gradients, so a comparison built on PySCF's tables cannot tell a real
    discrepancy from a basis one. Reading the same file both codes read removes
    the question.
    """
    from gen_cpu_validation import CARTESIAN, bse_to_pyscf, molecule_form

    symbols = {a[0] for a in atoms}

    # `make_auxmol` inherits `mol.cart`, so PySCF will quietly build a
    # spherical auxiliary set in Cartesian form if asked. metalquicha refuses
    # that combination outright -- libcint builds all three centres of a
    # fitting integral in one form -- so a reference for it would be a
    # reference for a deck that cannot run. Refuse it here too, and for the
    # same reason.
    orbital_form = molecule_form(basis, symbols)
    aux_form = molecule_form(auxbasis, symbols)
    if CARTESIAN in (orbital_form, aux_form) and orbital_form != aux_form:
        raise SystemExit(
            f"{basis} is {orbital_form} and {auxbasis} is {aux_form}; "
            "metalquicha refuses a fitting integral whose centres disagree, so "
            "there is nothing here to be a reference for")

    mol = gto.Mole()
    mol.atom = [(s, (x, y, z)) for s, x, y, z in atoms]
    mol.unit = unit
    mol.basis = {s: bse_to_pyscf(basis, s) for s in symbols}
    mol.cart = orbital_form == CARTESIAN
    mol.verbose = 0
    mol.build()
    auxmol = df.addons.make_auxmol(
        mol, {s: bse_to_pyscf(auxbasis, s) for s in symbols})
    return mol, auxmol


def main():
    lib.num_threads(8)
    cases = [
        ("H2 / sto-3g", [["H", 0.0, 0.0, -0.7], ["H", 0.0, 0.0, 0.7]],
         "sto-3g", "cc-pvdz-rifit", 0),
        ("H2O / sto-3g", [["O", 0.0, 0.0, 0.0], ["H", 0.0, -1.4308, 1.1078],
                          ["H", 0.0, 1.4308, 1.1078]],
         "sto-3g", "cc-pvdz-rifit", 0),
        ("HCN / sto-3g", [["H", 0.0, 0.0, -2.0], ["C", 0.0, 0.0, 0.0],
                          ["N", 0.0, 0.0, 2.2]],
         "sto-3g", "cc-pvdz-rifit", 0),
        # The frozen finite difference freezes the same core: the relaxed
        # density's occupied-frozen block comes from the Lagrangian rather
        # than the amplitudes, and only a derivative of the frozen energy
        # sees it as a derivative.
        ("H2O / 6-31g, frozen 1", [["O", 0.0, 0.0, 0.0],
                                   ["H", 0.0, -1.4308, 1.1078],
                                   ["H", 0.0, 1.4308, 1.1078]],
         "6-31g", "cc-pvtz-rifit", 1),
        ("HCN / sto-3g, frozen 2", [["H", 0.0, 0.0, -2.0],
                                    ["C", 0.0, 0.0, 0.0],
                                    ["N", 0.0, 0.0, 2.2]],
         "sto-3g", "cc-pvdz-rifit", 2),
    ]

    worst_overall = 0.0
    for label, atoms, basis, auxbasis, frozen in cases:
        mol, auxmol = build(atoms, basis, auxbasis)
        energy, mf = ri_mp2_energy(mol, auxmol, frozen=frozen)
        analytic = ri_mp2_gradient(mol, auxmol, mf, verbose=True, frozen=frozen)
        numeric = finite_difference(atoms, basis, auxbasis, frozen=frozen)

        print(f"== {label}   E(RI-MP2) = {energy:.12f}")
        for a in range(mol.natm):
            print("   %2d  analytic %16.11f %16.11f %16.11f   fd %16.11f %16.11f %16.11f"
                  % (a + 1, *analytic[a], *numeric[a]))
        worst = np.abs(analytic - numeric).max()
        worst_overall = max(worst_overall, worst)
        print(f"   largest deviation from finite difference: {worst:.4e}")
        print(f"   |sum over atoms|:                         "
              f"{np.abs(analytic.sum(axis=0)).max():.4e}")

    print()
    print("worst over all cases: %.4e" % worst_overall)
    return 0 if worst_overall < 1e-5 else 1


if __name__ == "__main__":
    raise SystemExit(main())
