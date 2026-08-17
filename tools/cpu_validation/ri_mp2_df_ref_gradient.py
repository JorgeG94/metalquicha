#!/usr/bin/env python3
"""RI-MP2 over a *density-fitted* reference: the analytic gradient, in numpy.

`ri_mp2_gradient.py` next door differentiates an exact-ERI Hartree-Fock
reference with fitted correlation, which is what an ``ri-mp2`` deck means. This
one differentiates the other combination -- the reference fitted as well, which
is what ``keywords.scf.density_fitting`` asks for and what the Fortran side
refuses today.

**What actually changes.** Only the reference. The amplitudes, the three- and
two-index correlation densities, the Lagrangian and the one-particle blocks are
identical, and are imported rather than copied. What moves is every place the
*reference* two-electron interaction appears:

1. `J - K/2` for a density that is not idempotent -- the Z-vector's right-hand
   side, the `S^(1)` occupied term, and the response operator inside the solve.
   The SCF's own build routes exchange through the occupied orbitals and cannot
   be reused: a relaxed correlation density is symmetric and *indefinite* and
   has no such orbitals.
2. The skeleton derivative of that interaction, which in the exact case is a
   four-centre `int2e_ip1` contraction and here becomes three- and two-centre
   terms with their own `Gamma` and `Omega`.

The second is the whole of the new derivation, and the exchange half of it is
where the factor conventions are expected to break first.

Run it directly to check against finite differences::

    python tools/cpu_validation/ri_mp2_df_ref_gradient.py
"""

from __future__ import annotations

import numpy as np
from pyscf import df, gto, lib, scf

from ri_mp2_gradient import (amplitudes, build, fitted_tensor,
                             metric_inverse_sqrt)


# --------------------------------------------------------------------------
# the fitted reference interaction
# --------------------------------------------------------------------------

class FittedReference:
    """`(mu nu|P)`, `J^-1`, and the two things built from them.

    Held together because every consumer wants more than one of them, and
    because the gradient has to differentiate the *same* fitting the energy
    used -- a second set of integrals built with a different threshold would be
    a gradient of a different function.
    """

    def __init__(self, mol, auxmol, k_scale=1.0):
        nao = mol.nao
        self.mol = mol
        self.auxmol = auxmol
        self.k_scale = k_scale
        self.three = df.incore.aux_e2(
            mol, auxmol, intor="int3c2e", aosym="s1").reshape(nao, nao, -1)
        self.metric = auxmol.intor("int2c2e")
        self.jinv = np.linalg.inv(self.metric)
        # `Y^P = sum_Q J^-1_PQ M^Q`, the metric applied once to the three-centre
        # integrals. Both the exchange build and its derivative want it, and it
        # is n_aux n^2 to store against n_aux n^3 to rebuild.
        self.y = np.einsum("PQ,uvQ->uvP", self.jinv, self.three, optimize=True)

    def rho(self, dm):
        """`J^-1 g`, with `g_P = sum_uv D_uv (uv|P)`."""
        return self.jinv @ np.einsum("uvP,uv->P", self.three, dm)

    def veff(self, dm):
        """`J - k K/2` for an arbitrary symmetric density.

        The general form, not the SCF's. Per auxiliary function

            J += (sum_uv M^P_uv D_uv) Y^P,    K += M^P D Y^P

        at `n^3 n_aux` rather than the `n^2 n_occ n_aux` an idempotent density
        allows. Evaluated a handful of times per gradient rather than once per
        SCF iteration, so the cost is not the point -- correctness for an
        indefinite density is.
        """
        vj = np.einsum("uvP,P->uv", self.three, self.rho(dm), optimize=True)
        vk = np.einsum("ulP,ls,svP->uv", self.three, dm, self.y, optimize=True)
        return vj - 0.5 * self.k_scale * vk

    def gamma_omega(self, dm_a, dm_b):
        """`Gamma^P_uv` and `Omega_PQ` for `F = sum_uv G[A]_uv B_uv`.

        This is the object the exact path reaches through `int2e_ip1`: the
        skeleton derivative of the reference two-electron interaction, with one
        density inside the operator and another contracted against it. Both are
        held fixed as matrices; only the integrals move.

        **Coulomb.** `F_J = sum_PQ g^A_P J^-1_PQ g^B_Q`, so differentiating the
        three-centre integrals picks up one term from each density and
        `d(J^-1) = -J^-1 (dJ) J^-1` gives the two-centre one.

        **Exchange.** `F_K = -(k/2) sum_PQ J^-1_PQ Tr(B M^P A M^Q)`. Both
        `M` factors differentiate, and the two contributions are transposes of
        one another rather than equal -- which is the trap, because `Gamma` is
        symmetric either way and only the assembled gradient can tell.

        The two-centre half reuses `S^P = B Y^P A` rather than forming the
        `n_aux^2` intermediate `W_PQ = Tr(B M^P A M^Q)`: contracting `S^R`
        against `Y^S` is the same number, since `J^-1 W J^-1` applied on both
        sides is what the metric derivative wants anyway.
        """
        k = self.k_scale
        rho_a = self.rho(dm_a)
        rho_b = self.rho(dm_b)

        gamma = (np.einsum("uv,P->uvP", dm_a, rho_b, optimize=True)
                 + np.einsum("uv,P->uvP", dm_b, rho_a, optimize=True))
        omega = -np.outer(rho_a, rho_b)

        # S^P = B Y^P A, one n^3 product per auxiliary function.
        s = np.einsum("ul,lsP,sv->uvP", dm_b, self.y, dm_a, optimize=True)
        gamma = gamma - 0.5 * k * (s + s.transpose(1, 0, 2))
        omega = omega + 0.5 * k * np.einsum("uvR,uvS->RS", s, self.y,
                                            optimize=True)

        return gamma, 0.5 * (omega + omega.T)


def reference_gradient(mol, auxmol, gamma, omega):
    """`sum_P Gamma^P_uv d(uv|P) + sum_PQ Omega_PQ dJ_PQ`, per atom.

    libcint's `ip` integrals differentiate the electronic coordinate, so the
    derivative with respect to a nucleus carrying the function is minus them --
    which is the single sign in this file that cannot be read off the algebra.
    """
    nao = mol.nao
    de = np.zeros((mol.natm, 3))
    aoslices = mol.aoslice_by_atom()
    auxslices = auxmol.aoslice_by_atom()

    ip1 = df.incore.aux_e2(mol, auxmol, intor="int3c2e_ip1", aosym="s1",
                           comp=3).reshape(3, nao, nao, -1)
    ip2 = df.incore.aux_e2(mol, auxmol, intor="int3c2e_ip2", aosym="s1",
                           comp=3).reshape(3, nao, nao, -1)
    for a in range(mol.natm):
        p0, p1 = aoslices[a, 2], aoslices[a, 3]
        # `(mu nu|P)` is symmetric in mu and nu, so the ket derivative is the
        # bra one transposed rather than a second integral.
        de[a] -= np.einsum("xuvP,uvP->x", ip1[:, p0:p1], gamma[p0:p1])
        de[a] -= np.einsum("xuvP,vuP->x", ip1[:, p0:p1], gamma[:, p0:p1])
        q0, q1 = auxslices[a, 2], auxslices[a, 3]
        de[a] -= np.einsum("xuvP,uvP->x", ip2[:, :, :, q0:q1],
                           gamma[:, :, q0:q1])

    j1 = auxmol.intor("int2c2e_ip1", comp=3)
    for a in range(mol.natm):
        q0, q1 = auxslices[a, 2], auxslices[a, 3]
        de[a] -= np.einsum("xPQ,PQ->x", j1[:, q0:q1], omega[q0:q1])
        de[a] -= np.einsum("xQP,PQ->x", j1[:, q0:q1], omega[:, q0:q1])
    return de


# --------------------------------------------------------------------------
# energy and gradient
# --------------------------------------------------------------------------

def energy(mol, auxmol, conv_tol=1e-13):
    """Density-fitted RHF plus a fitted MP2 correlation energy.

    One auxiliary basis for both, which is what our decks can express:
    `model.aux_basis` is the only place a fitting set is named. PySCF would
    happily use different ones, and a comparison that let it would be comparing
    two approximations rather than two implementations.
    """
    mf = scf.RHF(mol).density_fit(auxbasis=auxmol.basis)
    mf.conv_tol = conv_tol
    mf.kernel()
    if not mf.converged:
        raise SystemExit("DF-RHF did not converge")
    nocc = mol.nelectron // 2
    b, _, _ = fitted_tensor(mol, auxmol, mf.mo_coeff, nocc)
    x, ovov = amplitudes(b, mf.mo_energy, nocc)
    e_corr = np.einsum("ijab,ijab->", x, 2 * ovov - ovov.transpose(0, 1, 3, 2))
    return mf.e_tot + e_corr, mf


def gradient(mol, auxmol, mf, verbose=False):
    """dE/dR for DF-RHF plus fitted MP2, in Hartree/Bohr."""
    nocc = mol.nelectron // 2
    nao, nmo = mf.mo_coeff.shape
    nvir = nmo - nocc
    mo, e_mo = mf.mo_coeff, mf.mo_energy
    orbo, orbv = mo[:, :nocc], mo[:, nocc:]

    ref = FittedReference(mol, auxmol)

    b, metric, jm12 = fitted_tensor(mol, auxmol, mo, nocc)
    x, _ = amplitudes(b, e_mo, nocc)

    # ---- correlation densities and the Lagrangian, unchanged -------------
    xbar = 2 * x - x.transpose(0, 1, 3, 2)
    gamma = np.einsum("ijab,Qjb,PQ->Pia", xbar, b, jm12, optimize=True)
    gamma_pq = np.einsum("Pia,Ria,QR->PQ", gamma, b, jm12, optimize=True)

    doo = -(2 * np.einsum("kiab,kjab->ij", x, x, optimize=True)
            - np.einsum("kiab,kjba->ij", x, x, optimize=True))
    dvv = (2 * np.einsum("ijca,ijcb->ba", x, x, optimize=True)
           - np.einsum("ijca,ijbc->ba", x, x, optimize=True))

    dm1mo = np.zeros((nmo, nmo))
    dm1mo[:nocc, :nocc] = doo + doo.T
    dm1mo[nocc:, nocc:] = dvv + dvv.T

    pq_p = np.einsum("up,uvP,vq->pqP", mo, ref.three, mo, optimize=True)
    lag_o = 2 * np.einsum("Pia,qaP->iq", gamma, pq_p[:, nocc:, :], optimize=True)
    lag_v = 2 * np.einsum("Pia,iqP->aq", gamma, pq_p[:nocc, :, :], optimize=True)
    imat = np.zeros((nmo, nmo))
    imat[:, :nocc] = -lag_o.T
    imat[:, nocc:] = -lag_v.T

    # ---- the Z-vector, over the fitted operator --------------------------
    hf_dm = mf.make_rdm1()
    dm1_ao = mo @ dm1mo @ mo.T
    vhf = ref.veff(dm1_ao) * 2
    xvo = orbv.T @ vhf @ orbo + imat[:nocc, nocc:].T - imat[nocc:, :nocc]

    dvo = solve_z_vector(ref, mo, e_mo, xvo, nocc, nvir)
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
    vhf_s1occ = p_occ @ ref.veff(dm1_ao + dm1_ao.T) @ p_occ

    dm1p = hf_dm + dm1_ao * 2
    dm1_total = hf_dm + dm1_ao
    zeta += (orbo * e_mo[:nocc]) @ orbo.T * 2.0

    # ---- assembly --------------------------------------------------------
    grad = scf.RHF(mol).nuc_grad_method()
    hcore_deriv = grad.hcore_generator(mol)
    s1 = grad.get_ovlp(mol)
    aoslices = mol.aoslice_by_atom()

    de = np.zeros((mol.natm, 3))
    de += gradient_three_centre(mol, auxmol, gamma, mo, nocc)
    de += gradient_two_centre(auxmol, gamma_pq, mol)
    # The reference two-electron term, which is the whole of what differs from
    # the exact-reference gradient next door.
    #
    # **The half is not a fudge.** `gamma_omega` returns the honest skeleton
    # derivative of `Tr(G[A] B)`, checked on its own against a finite difference
    # with two fixed indefinite matrices. What the assembly wants is
    # `(1/2) Tr(G^x[D] D) + Tr(G^x[D] P)`, the reference's own two-electron
    # gradient plus the cross term -- which is half of `Tr(G^x[D] (D + 2P))`.
    # The exact path reaches the same half by building only two of the four
    # differentiated positions in `hf_two_electron_deriv` and letting the
    # integral's symmetry supply the rest.
    de += 0.5 * reference_gradient(mol, auxmol, *ref.gamma_omega(hf_dm, dm1p))
    for a in range(mol.natm):
        p0, p1 = aoslices[a, 2], aoslices[a, 3]
        de[a] += np.einsum("xij,ij->x", s1[:, p0:p1], im1[p0:p1])
        de[a] += np.einsum("xij,ij->x", s1[:, p0:p1], im1.T[p0:p1])
        de[a] += np.einsum("xij,ji->x", hcore_deriv(a), dm1_total)
        de[a] -= np.einsum("xij,ij->x", s1[:, p0:p1], zeta[p0:p1])
        de[a] -= np.einsum("xij,ij->x", s1[:, p0:p1], zeta.T[p0:p1])
        de[a] -= 2 * np.einsum("xij,ij->x", s1[:, p0:p1], vhf_s1occ[p0:p1])
    de += grad.grad_nuc()

    if verbose:
        print("  |gamma| %.4e  |z| %.4e  |D| %.4e  |Imat| %.4e"
              % (abs(gamma).sum(), abs(dvo).sum(), abs(dm1_ao).sum(),
                 abs(imat).sum()))
    return de


def solve_z_vector(ref, mo, e_mo, xvo, nocc, nvir):
    """The Z-vector with the fitted response operator.

    Densely, for the reason the exact-reference prototype gives: PySCF's Krylov
    CPHF carries error its `tol` does not control, and at reference sizes the
    operator fits in memory.
    """
    orbo, orbv = mo[:, :nocc], mo[:, nocc:]
    gaps = (e_mo[nocc:, None] - e_mo[None, :nocc]).ravel()
    n = nvir * nocc

    def apply(vec):
        d = orbv @ vec.reshape(nvir, nocc) @ orbo.T
        return 2 * (orbv.T @ ref.veff(d + d.T) @ orbo).ravel()

    operator = np.zeros((n, n))
    basis = np.zeros(n)
    for k in range(n):
        basis[:] = 0.0
        basis[k] = 1.0
        operator[:, k] = apply(basis)
    operator += np.diag(gaps)
    return np.linalg.solve(operator, -xvo.ravel()).reshape(nvir, nocc)


# The correlation terms, imported by value rather than by name because the
# exact-reference module builds its own `three` and this one already has it.
from ri_mp2_gradient import (gradient_three_centre,  # noqa: E402
                             gradient_two_centre)


# --------------------------------------------------------------------------
# checking
# --------------------------------------------------------------------------

def finite_difference(atoms, basis, auxbasis, step=2.5e-3):
    """Central differences of the same energy the analytic gradient claims."""
    de = np.zeros((len(atoms), 3))
    for a in range(len(atoms)):
        for c in range(3):
            plus = [list(x) for x in atoms]
            minus = [list(x) for x in atoms]
            plus[a][1 + c] += step
            minus[a][1 + c] -= step
            e_plus = energy(*build(plus, basis, auxbasis))[0]
            e_minus = energy(*build(minus, basis, auxbasis))[0]
            de[a, c] = (e_plus - e_minus) / (2 * step)
    return de


def richardson(atoms, basis, auxbasis, step=2.5e-3):
    """Two central differences combined to kill the `h^2` truncation term.

    Worth the second pass here rather than a looser bound. A plain central
    difference at 2.5e-3 leaves 5e-6 on HCN, which is exactly the size of a
    factor-of-two error in one of the two exchange `Gamma` halves -- so the
    cheap check cannot tell the finished routine from a broken one. Measured:
    5e-3 gives 2.2e-5, 2.5e-3 gives 5.4e-6 and 1.25e-3 gives 1.3e-6, shrinking
    by 4.00 and 4.04, and `(4 D(h/2) - D(h))/3` leaves 1.9e-8.
    """
    coarse = finite_difference(atoms, basis, auxbasis, step=step)
    fine = finite_difference(atoms, basis, auxbasis, step=step / 2)
    return (4.0 * fine - coarse) / 3.0


def main():
    lib.num_threads(8)
    cases = [
        ("H2 / sto-3g", [["H", 0.0, 0.0, -0.7], ["H", 0.0, 0.0, 0.7]],
         "sto-3g", "cc-pvdz-rifit"),
        ("H2O / sto-3g", [["O", 0.0, 0.0, 0.0], ["H", 0.0, -1.4308, 1.1078],
                          ["H", 0.0, 1.4308, 1.1078]],
         "sto-3g", "cc-pvdz-rifit"),
        ("HCN / sto-3g", [["H", 0.0, 0.0, -2.0], ["C", 0.0, 0.0, 0.0],
                          ["N", 0.0, 0.0, 2.2]],
         "sto-3g", "cc-pvdz-rifit"),
    ]

    worst_overall = 0.0
    for label, atoms, basis, auxbasis in cases:
        mol, auxmol = build(atoms, basis, auxbasis)
        e_tot, mf = energy(mol, auxmol)
        analytic = gradient(mol, auxmol, mf, verbose=True)
        numeric = richardson(atoms, basis, auxbasis)

        print(f"== {label}   E(RI-MP2/DF-RHF) = {e_tot:.12f}")
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
    # Two orders below the plain central difference this used to accept, which
    # is the point of extrapolating: at 1e-5 a wrong exchange factor passes.
    return 0 if worst_overall < 1e-7 else 1


if __name__ == "__main__":
    raise SystemExit(main())
