#!/usr/bin/env python3
"""A reference RI-MP2 analytic Hessian, in numpy over PySCF integrals.

Nothing computes this, and a Fortran implementation has no one to check
against: PySCF has no RI-MP2 Hessian (it has no RI-MP2 gradient either --
``ri_mp2_gradient.py`` next door is that missing reference). This is the
second derivative of the same energy: an *exact-ERI* restricted Hartree-Fock
reference with a *fitted* correlation energy, optionally frozen-core.

Run it directly::

    python tools/cpu_validation/ri_mp2_hessian.py [--cases 0,1,2,3]

It prints a pass/fail line per check and a cost table, and exits non-zero on
any failure. The oracle for the total Hessian is five-point differences of the
analytic RI-MP2 gradient, which ``ri_mp2_gradient.py`` validates against
differences of the energy.

What is returned is the **correlation block**: the Hessian of E_corr alone, to
be added to the RHF Hessian (``mf.Hessian()``) and nothing else. Nuclear
repulsion and the reference are not in it.

THE FORMULATION
---------------

The conventional MP2 Hessian in ``mqc_czt_mp2_hessian.f90`` is pycc's
``_hessian_blocks`` transcribed. This script is that same ladder, with exactly
one change of principle:

    **the correlation functional contracts the fitted tensor
    R_pqrs = (pr|P) J^-1_PQ (Q|qs); everything that belongs to the reference
    keeps the exact one.**

Where each lives, in pycc's terms (<pq|rs> physicist, L = 2<pq|rs> - <pq|sr>):

    amplitudes, d_Y t2 numerator      R      (fitted)
    4 <pr|st> Gamma_qrst in I'        R      (fitted)
    Gamma . <pq|rs>^(X) in the sum    R      (fitted), skeleton includes J^(X)
    D_rs L_rpsq folds in I'           L      (exact)   -- the Fock operator
    orbital Hessian, d_X f, U^X       L      (exact)
    f^(X) = h^(X) + G^(X)[P_occ]      exact

The ladder, in the order the code runs it::

    E^X  = D_rel . f^(X) + Gamma . R^(X) + I' . S^(X)                    (gradient)
    H_XY = pass 1  : the same three contracted with second-derivative
                     integrals, densities fixed
         + orb(X,Y): 2 U^Y_ai X~^(X)_ai + S^(Y)_pq I''~^(X)_pq [+ P^(X)_pq f^(Y)_pq]
         + d_Y D_rel . f^(X) + d_Y Gamma . R^(X) + d_Y I' . S^(X)         (response)

``X^(X) = I'^(X)T - I'^(X)`` and ``I''^(X)`` are the *skeleton Lagrangian*'s
carriers; the bracket is the frozen core's independent core<->active rotation.
``d_Y`` is the full first derivative along Y, orbital response included, with
a **non-canonical** gauge: U_oo = U_vv = -S/2, U_ov = -S - U_vo, U_vo from the
coupled-perturbed equations. Amplitudes then respond through the full oo and
vv blocks of d_Y f, not only its diagonal. That is also what makes the
degenerate pi orbitals of HCN a non-event (no divide by an orbital-energy gap
anywhere except the frozen core's core-active block, whose gap is never small).

One batched Z-vector solve closes the response: for each Y,
``A z^Y = dX^Y - G^Y z`` with the unperturbed operator A, and the solved block
enters D_rel^Y and I'^Y linearly, which is why the four-index derivative need
not outlive the loop in the dense path (and nothing four-index exists in the
lean one).

CONVENTIONS THAT COST TIME (each was either measured or the reason a check exists)
---------------------------------------------------------------------------------

* ``D_rel`` has **minus** the Z-vector in its virtual-occupied block, where the
  Z-vector solves ``A z = I'_ia - I'_ai`` (unrelaxed density) -- pycc's ``-z^Y``.
  Reproduced to 1e-13 by the first check below.
* ``I'`` is not symmetric; only its symmetric part multiplies S^(X), and the
  vo block of ``I''`` is the transposed ov block.
* The fitted Lagrangian's two-particle term is, with Gamma^P_ia = u_ijab c^P_jb,
  c = J^-1 (jb|P), u = 2t - t^swap, and ``a`` any (n, n, naux) integral-like
  tensor:

      4 <pr|st> Gamma_qrst   ->   4 sum_P A^P_pa Gamma^P_ia      (q = i)
                                  4 sum_P A^P_ip Gamma^P_ia      (q = a)

  which is ``gp`` below. The same-shaped term in ``ri_mp2_gradient.py`` is
  ``imat = -L^T``; the transposes there and the slot here are the same thing.
* The skeleton derivative of that term is NOT ``gp(I^(X), Gamma)`` alone: the
  fitting coefficient c = J^-1 I moves too, so the derivative of
  R = I c is  I^(X) c + c I^(X) - c J^(X) c  and the extra pieces are
  ``gp(c, Theta)``, Theta^Q_ia = sum_jb u_ijab I^Q_jb - sum_R J_QR Gamma^R_ia.
  Dropping them is invisible to translational invariance.
* An orbital-basis rotation acts on an integral tensor as U^T I + I U (covariant)
  and on a density-like one as U D + D U^T (contravariant). The check on
  Gamma^P (AO) failed by 9e-2 until the check, not the Hessian, was fixed.
* ``d c`` is conditioned by the metric, and sits at 1e-8 against differences
  where everything else sits at 1e-10: expected, not a defect.
* The inverse metric is ``jm12 @ jm12`` from ``ri_mp2_gradient.metric_inverse_sqrt``
  (eigenvalue threshold 1e-10), so the energy being differentiated is the one the
  gradient differentiates. Its derivative is ``-J^-1 J^(X) J^-1`` as long as no
  mode is dropped, which none is for these auxiliary sets.

THE LEAN PATH (nothing n^4)
---------------------------

Every use of a four-index tensor reduces to one of two things.

1. The reference operator ``T1(D)_pq = sum_rs D_rs L_rpsq = (2J[D] - K[D]^T)`` in
   the MO basis, for a general (non-symmetric) MO matrix D, and its first
   derivative at fixed MO coefficients ``T1x(D)[X]``. The total derivative of
   ``T1(D)`` along Y at a fixed MO matrix D is

       T1x(D)[Y] + T1(U D + D U^T) + U^T T1(D) + T1(D) U

   (``Lean.DT1``): one skeleton derivative build per fixed density, once, and
   one ordinary Fock-type build per Y. ``d_Y f``, the core<->active rotation,
   the CPHF right-hand side, every ``D . dL`` fold and the Z-vector's ``G^Y z``
   are all this. In Fortran these are the fitted-or-exact J/K builds that
   ``make_h1_atom`` and ``cphf_solve`` already perform.
2. The fitted tensor, only ever against the cumulant, as the three-index
   objects Gamma^P_ia = u c, gamma_PQ = Gamma c and Theta above.

Pass 1 splits the same way. The reference share (D_rel h^(XY), the Fock fold
D_rel (x) P_occ, I' S^(XY)) is the conventional code's and is here *differenced
from analytic first derivatives*, not rebuilt from second-derivative
integrals. The fitted share is analytic:

    4 Gamma^P_mn (mn|P)^(XY) - 2 gamma_PQ (P|Q)^(XY)        [second derivatives]
  + 4 (d_Y Gamma) . (ia|P)^(X) - 2 (d_Y gamma) . (P|Q)^(X)   [cross terms]

    d_Y Gamma = J^-1 Theta^Y,
    d_Y gamma = d_Y Gamma c + Gamma J^-1 (I^(Y) - J^(Y) c)       (u fixed)

``contract_second_3c2c`` is the first, with libcint's ``ipip1``, ``ipvip1``,
``ip1ip2`` and ``ipip2`` (three-centre) and ``ipip1``, ``ip1ip2`` (two-centre).
Each ordered pair of centres is one array, picked out by the atom that owns
each function; the nine components of every array are ``[first derivative
direction, second derivative direction]`` in the order the integral name
reads, and a derivative with respect to an electron coordinate is minus the
nuclear one, so a product of two is plus. The (nu, mu) pair of ``ipvip1`` is
the array indexed ``[y, x]`` with the functions *not* swapped -- getting that
one wrong was the only error here, 6.6 on a random density. Checked against
differences of the first-derivative contraction to 8e-14.

WHAT IS CHECKED, AND AGAINST WHAT
---------------------------------

Per case (see ``CASES``; water at 6-31G all-electron and frozen-core with a
triple-zeta fitting set, and HCN frozen-core, whose pi orbitals are degenerate):

* the unperturbed identity: ``D_rel f + Gamma R + I' S`` equals the validated
  RI-MP2 gradient minus the RHF gradient (1e-13);
* the dense four-index path (transcription of the Fortran) against the lean
  path; with the same pass 1 they agree to rounding;
* the analytic fitted pass 1 against differences of the first-derivative
  contraction;
* each first-order quantity against its own finite difference, **at fixed U** --
  geometry displaced and orbitals taken as C0 (1 + h U), not re-converged, so
  the difference isolates one formula in one gauge: d(ia|P), d(J^-1/2)-based
  dB, d c, d(ia|jb)_RI, d f, and d t2 (the non-canonical amplitude equations
  solved at the displaced point);
* the gauge-free AO objects D_rel, I' and Gamma^P against differences across
  *re-converged* SCFs, mapped from the analytic side with the orbital rotation
  (D -> C (dD + U D + D U^T) C^T). These pin the perturbed Z-vector solution
  and the perturbed Lagrangian, including the frozen core's Sylvester
  derivative;
* the total (RHF Hessian from PySCF + this block) against five-point
  differences of ``ri_mp2_gradient.ri_mp2_gradient``;
* the fitting error: this block minus the conventional MP2 correlation block
  from the same code with R replaced by the exact tensor.

COST STRUCTURE
--------------

``cost_report`` prints counts for the last case; the shape is in
``run_case``. With n MOs, o active occupied, v virtual, P auxiliary functions,
N atoms: resident once are the 3-index MO integrals and c (n^2 P), Gamma, t2
and a handful of n^2 matrices; resident over all 3N perturbations only n^2
matrices (S^(X), f^(X), U^X, D_rel^Y, I'^Y, three T1x stacks). Per
perturbation, transient: the total-derivative 3-index integrals (only the
(n, v) and (o, n) blocks are read), d(ia|jb) and d t2 (o^2 v^2), the
three-index response dGamma~ = du c (o v P), and a few n^2 matrices. The pair
sum 4 dGamma~ . (ia|P)^(X) - 2 dgamma~ . (P|Q)^(X) does not depend on the
Z-vector, so it is formed inside the perturbation loop and nothing o v P
outlives it. The prototype keeps (ia|P)^(X) for all X; a port does not need to,
because that sum -- together with the cross terms of pass 1, which are the same
contraction with a different three-index density -- is the existing RI-MP2
*gradient's* three- and two-centre contraction handed a different density:
one gradient-style sweep per perturbation.
"""

from __future__ import annotations

import os
import sys

import numpy as np
from pyscf import df, scf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ri_mp2_gradient as grad_ref  # noqa: E402

ES = True


# ---------------------------------------------------------------------------
# first-derivative integrals, one array per atom
# ---------------------------------------------------------------------------

def owners(mol):
    sl = mol.aoslice_by_atom()
    m = np.zeros((mol.natm, mol.nao), dtype=bool)
    for a in range(mol.natm):
        m[a, sl[a, 2]:sl[a, 3]] = True
    return m


def first_derivatives(mol, auxmol, four_centre=True):
    """d/dR_A of every integral (functions on A move), per atom, (3, ...)."""
    nao, naux = mol.nao, auxmol.nao
    own, aown = owners(mol), owners(auxmol)
    out = {"S": [], "H": [], "J": [], "I": []}
    if four_centre:
        out["E"] = []
    ip_s = mol.intor("int1e_ipovlp")
    hgen = scf.RHF(mol).nuc_grad_method().hcore_generator(mol)
    ip3_1 = df.incore.aux_e2(mol, auxmol, intor="int3c2e_ip1", aosym="s1",
                             comp=3).reshape(3, nao, nao, naux)
    ip3_2 = df.incore.aux_e2(mol, auxmol, intor="int3c2e_ip2", aosym="s1",
                             comp=3).reshape(3, nao, nao, naux)
    ip2c = auxmol.intor("int2c2e_ip1", comp=3)
    if four_centre:
        ip4 = mol.intor("int2e_ip1").reshape(3, nao, nao, nao, nao)
    for a in range(mol.natm):
        m, q = own[a], aown[a]
        t = ip_s * m[None, :, None]
        out["S"].append(-t - t.transpose(0, 2, 1))
        out["H"].append(np.asarray(hgen(a)))
        t = ip3_1 * m[None, :, None, None]
        out["I"].append(-t - t.transpose(0, 2, 1, 3)
                        - ip3_2 * q[None, None, None, :])
        t = ip2c * q[None, :, None]
        out["J"].append(-t - t.transpose(0, 2, 1))
        if four_centre:
            t = ip4 * m[None, :, None, None, None]
            out["E"].append(-t - t.transpose(0, 2, 1, 3, 4)
                            - t.transpose(0, 3, 4, 1, 2)
                            - t.transpose(0, 3, 4, 2, 1))
    return out


# ---------------------------------------------------------------------------
# reference state
# ---------------------------------------------------------------------------

class Ref:
    """RHF at one geometry plus every AO and MO object the Hessian reads."""


def make_ref(atoms, basis, auxbasis, nf=0):
    r = Ref()
    r.atoms, r.basis, r.auxbasis, r.nf = atoms, basis, auxbasis, nf
    mol, auxmol = grad_ref.build(atoms, basis, auxbasis)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-13
    mf.conv_tol_grad = 1e-10
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        raise SystemExit("RHF did not converge")
    r.mol, r.auxmol, r.mf = mol, auxmol, mf
    r.C, r.eps = mf.mo_coeff, mf.mo_energy
    r.nao, r.n = mol.nao, mf.mo_coeff.shape[1]
    r.no = mol.nelectron // 2
    r.nv = r.n - r.no
    r.naux = auxmol.nao
    r.natm = mol.natm
    return r


def ao_to_mo4(a, c):
    return np.einsum("up,vq,lr,ws,uvlw->pqrs", c, c, c, c, a, optimize=ES)


def ri_objects(r, c=None):
    """AO and MO three-index integrals, metric, and fitted coefficients."""
    c = r.C if c is None else c
    three = df.incore.aux_e2(r.mol, r.auxmol, intor="int3c2e", aosym="s1")
    three = three.reshape(r.nao, r.nao, -1)
    metric = r.auxmol.intor("int2c2e")
    jm12 = grad_ref.metric_inverse_sqrt(metric)
    jinv = jm12 @ jm12
    imo = np.einsum("up,uvP,vq->pqP", c, three, c, optimize=ES)
    cmo = imo @ jinv
    return three, metric, jinv, imo, cmo


def phys(chem):
    """(pq|rs) -> <pq|rs> = (pr|qs)."""
    return chem.transpose(0, 2, 1, 3)


def setup(r, corr):
    """Unperturbed tensors. corr is 'ri' or 'exact' (what the correlation
    functional contracts); the reference operators are always exact."""
    no, nf = r.no, r.nf
    eri_ao = r.mol.intor("int2e").reshape((r.nao,) * 4)
    r.eri = phys(ao_to_mo4(eri_ao, r.C))
    r.L = 2 * r.eri - r.eri.transpose(0, 1, 3, 2)
    r.three, r.metric, r.jinv, r.imo, r.cmo = ri_objects(r)
    if corr == "ri":
        r.R = phys(np.einsum("pqP,rsP->pqrs", r.imo, r.cmo, optimize=ES))
    else:
        r.R = r.eri
    r.corr = corr
    e = r.eps
    eo, ev = e[nf:no], e[no:]
    d = (eo[:, None, None, None] + eo[None, :, None, None]
         - ev[None, None, :, None] - ev[None, None, None, :])
    r.t2 = r.R[nf:no, nf:no, no:, no:] / d
    return r


# ---------------------------------------------------------------------------
# pycc ladder, unperturbed
# ---------------------------------------------------------------------------

def cumulant(t2, nf, no, n):
    u = 2 * t2 - t2.transpose(0, 1, 3, 2)
    gam = np.zeros((n,) * 4)
    gam[nf:no, nf:no, no:, no:] = u
    gam[no:, no:, nf:no, nf:no] = u.transpose(2, 3, 0, 1)
    return gam


def unrelaxed_density(t2, nf, no, n):
    l2 = 2 * (2 * t2 - t2.transpose(0, 1, 3, 2))
    d0 = np.zeros((n, n))
    d0[nf:no, nf:no] = -np.einsum("imef,jmef->ij", t2, l2, optimize=ES)
    d0[no:, no:] = np.einsum("mnbe,mnae->ab", t2, l2, optimize=ES)
    return d0


def lagrangian(r, d, gam, rcorr=None):
    """I'_pq = -1/2 [eps_p (D+D^T) + delta_{q occ} D_rs (L_rpsq + L_rqsp)
    + 4 <pr|st> Gamma_qrst], the two-particle term over the correlation
    tensor, the rest over the reference's."""
    rc = r.R if rcorr is None else rcorr
    no = r.no
    ip = r.eps[:, None] * (d + d.T)
    ip += 4 * np.einsum("prst,qrst->pq", rc, gam, optimize=ES)
    fold = (np.einsum("rs,rpsq->pq", d, r.L, optimize=ES)
            + np.einsum("rs,rqsp->pq", d, r.L, optimize=ES))
    ip[:, :no] += fold[:, :no]
    return -0.5 * ip


def orbital_hessian(r):
    """A_{ia,jb} z_jb = (eps_a - eps_i) z_ia + (L_ajib + L_abij) z_jb."""
    no, nv = r.no, r.nv
    a = np.zeros((no, nv, no, nv))
    g = r.eps[None, no:] - r.eps[:no, None]
    for i in range(no):
        for k in range(nv):
            a[i, k, i, k] += g[i, k]
    oi, vi = slice(0, no), slice(no, r.n)
    # L[a, j, i, b] -> (i, a, j, b)
    t1 = r.L[vi, oi, oi, vi].transpose(2, 0, 1, 3)
    # L[a, b, i, j] -> (i, a, j, b)
    t2 = r.L[vi, vi, oi, oi].transpose(2, 0, 3, 1)
    a += t1 + t2
    return a.reshape(no * nv, no * nv)


def solve_unperturbed(r):
    no, nf, nv, n = r.no, r.nf, r.nv, r.n
    t2 = r.t2
    r.gam = cumulant(t2, nf, no, n)
    r.d0 = unrelaxed_density(t2, nf, no, n)
    ip0 = lagrangian(r, r.d0, r.gam)
    r.pco = None
    if nf > 0:
        r.pco = np.zeros((nf, no - nf))
        for i in range(no - nf):
            for c in range(nf):
                r.pco[c, i] = (ip0[c, nf + i] - ip0[nf + i, c]) \
                    / (r.eps[c] - r.eps[nf + i])
    r.A = orbital_hessian(r)
    r.Ainv = np.linalg.inv(r.A)
    x = (ip0[:no, no:] - ip0[no:, :no].T).copy()
    if nf > 0:
        for a in range(nv):
            for i in range(no):
                acc = 0.0
                for c in range(nf):
                    for j in range(no - nf):
                        acc += r.pco[c, j] * (r.L[no + a, nf + j, i, c]
                                              + r.L[no + a, c, i, nf + j])
                x[i, a] += acc
    r.z = (r.Ainv @ x.ravel()).reshape(no, nv)
    drel = r.d0.copy()
    if nf > 0:
        drel[nf:no, :nf] += r.pco.T
        drel[:nf, nf:no] += r.pco
    drel[no:, :no] -= r.z.T
    drel[:no, no:] -= r.z
    r.drel = drel
    r.ip = lagrangian(r, drel, r.gam)
    return r


# ---------------------------------------------------------------------------
# skeleton objects at the reference geometry
# ---------------------------------------------------------------------------

def first_order_skeletons(r, d1):
    """fx, sx and the skeleton MO integrals, stacked over X = 3*atom + comp."""
    no = r.no
    cocc = r.C[:, :no]
    pocc = cocc @ cocc.T
    fx, sx, ex, rx = [], [], [], []
    for a in range(r.natm):
        for c in range(3):
            sx.append(r.C.T @ d1["S"][a][c] @ r.C)
            e = d1["E"][a][c]
            f = r.C.T @ d1["H"][a][c] @ r.C
            f += r.C.T @ (2 * np.einsum("pqls,ls->pq", e, pocc, optimize=ES)
                          - np.einsum("plsq,ls->pq", e, pocc, optimize=ES)) @ r.C
            fx.append(f)
            ex.append(phys(ao_to_mo4(e, r.C)))
            ix = np.einsum("up,uvP,vq->pqP", r.C, d1["I"][a][c], r.C, optimize=ES)
            jx = d1["J"][a][c]
            rc = (np.einsum("pqP,rsP->pqrs", ix, r.cmo, optimize=ES)
                  + np.einsum("pqP,rsP->pqrs", r.cmo, ix, optimize=ES)
                  - np.einsum("pqP,PQ,rsQ->pqrs", r.cmo, jx, r.cmo, optimize=ES))
            rx.append(phys(rc))
    r.fx, r.sx, r.erix_e, r.erix_c0 = fx, sx, ex, rx
    if r.corr == "exact":
        r.erix_c = ex
    else:
        r.erix_c = rx
    return r


def corr_gradient_from_state(r):
    """E^X = D_rel f^(X) + Gam ERI^(X) + I' S^(X), correlation only."""
    g = np.zeros(3 * r.natm)
    for x in range(3 * r.natm):
        g[x] = (np.sum(r.drel * r.fx[x]) + np.sum(r.gam * r.erix_c[x])
                + np.sum(r.ip * r.sx[x]))
    return g.reshape(r.natm, 3)


# ---------------------------------------------------------------------------
# first-order orbitals
# ---------------------------------------------------------------------------

def solve_mo1(r):
    """U^X_ai from the coupled-perturbed equations, dense, plus -S/2 rows."""
    no, nv, n = r.no, r.nv, r.n
    mo1 = []
    for x in range(3 * r.natm):
        sx, fx = r.sx[x], r.fx[x]
        rhs = fx[no:, :no] - sx[no:, :no] * r.eps[None, :no]
        rhs -= 0.5 * (np.einsum("nm,anim->ai", sx[:no, :no], r.L[no:, :no, :no, :no],
                                optimize=ES)
                      + np.einsum("nm,amin->ai", sx[:no, :no], r.L[no:, :no, :no, :no],
                                  optimize=ES))
        u = -(r.Ainv @ rhs.T.ravel()).reshape(no, nv).T
        m = np.zeros((n, no))
        m[no:] = u
        m[:no] = -0.5 * sx[:no, :no]
        mo1.append(m)
    r.mo1 = mo1
    return r


def full_u(r, x):
    no, nf, n = r.no, r.nf, r.n
    sx, fx, mo1 = r.sx[x], r.fx[x], r.mo1[x]
    u = np.zeros((n, n))
    u[:, :no] = mo1
    u[no:, no:] = -0.5 * sx[no:, no:]
    u[:no, no:] = -sx[:no, no:] - mo1[no:, :].T
    if nf > 0:
        L, e = r.L, r.eps
        for i in range(nf, no):
            for c in range(nf):
                acc = fx[c, i] - sx[c, i] * e[i]
                acc -= 0.5 * np.einsum("nm,nm->", sx[:no, :no],
                                       L[c, :no, i, :no] + L[c, :no, i, :no].T)
                acc += np.einsum("dm,dm->", u[no:, :no],
                                 L[c, no:, i, :no] + L[c, :no, i, no:].T)
                u[c, i] = -acc / (e[c] - e[i])
                u[i, c] = -sx[c, i] - u[c, i]
    return u


def perturbed_fock(r, x, u):
    no = r.no
    sx, fx, L, e = r.sx[x], r.fx[x], r.L, r.eps
    df = fx + u * e[:, None] + u.T * e[None, :]
    # -1/2 S_nm (L_pnqm + L_pmqn)
    df -= 0.5 * (np.einsum("nm,pnqm->pq", sx[:no, :no], L[:, :no, :, :no], optimize=ES)
                 + np.einsum("nm,pmqn->pq", sx[:no, :no], L[:, :no, :, :no], optimize=ES))
    # + U_cm (L_pcqm + L_pmqc)
    df += (np.einsum("cm,pcqm->pq", u[no:, :no], L[:, no:, :, :no], optimize=ES)
           + np.einsum("cm,pmqc->pq", u[no:, :no], L[:, :no, :, no:], optimize=ES))
    return df


def perturbed_eri(erix, u, eri):
    d = np.einsum("tp,tqrs->pqrs", u, eri, optimize=ES)
    d += np.einsum("tq,ptrs->pqrs", u, eri, optimize=ES)
    d += np.einsum("tr,pqts->pqrs", u, eri, optimize=ES)
    d += np.einsum("ts,pqrt->pqrs", u, eri, optimize=ES)
    return d + erix


def perturbed_t2(r, dk, df):
    """d_X t2 from the ovov block of d_X <pq|rs>, ``dk[i, j, a, b]``."""
    nf, no = r.nf, r.no
    t2, e = r.t2, r.eps
    acc = dk.copy()
    dvv = df[no:, no:]
    doo = df[nf:no, nf:no]
    acc += np.einsum("ac,ijcb->ijab", dvv, t2, optimize=ES)
    acc += np.einsum("bc,ijac->ijab", dvv, t2, optimize=ES)
    acc -= np.einsum("ik,kjab->ijab", doo, t2, optimize=ES)
    acc -= np.einsum("jk,ikab->ijab", doo, t2, optimize=ES)
    eo, ev = e[nf:no], e[no:]
    d = (eo[:, None, None, None] + eo[None, :, None, None]
         - ev[None, None, :, None] - ev[None, None, None, :])
    return acc / d


def perturbed_onepdm(r, dt2):
    nf, no, n = r.nf, r.no, r.n
    t2 = r.t2
    l2 = 2 * (2 * t2 - t2.transpose(0, 1, 3, 2))
    la = 2 * (2 * dt2 - dt2.transpose(0, 1, 3, 2))
    dd = np.zeros((n, n))
    dd[nf:no, nf:no] = -(np.einsum("imef,jmef->ij", dt2, l2, optimize=ES)
                         + np.einsum("imef,jmef->ij", t2, la, optimize=ES))
    dd[no:, no:] = (np.einsum("mnbe,mnae->ab", dt2, l2, optimize=ES)
                    + np.einsum("mnbe,mnae->ab", t2, la, optimize=ES))
    return dd


def perturbed_lagrangian(r, df, deri_c, dl, d, dd, dgam):
    """d_X I'; the two-particle term over the correlation tensor."""
    no = r.no
    dsym = d + d.T
    dip = df @ dsym
    dip += r.eps[:, None] * (dd + dd.T)
    dip += 4 * (np.einsum("prst,qrst->pq", deri_c, r.gam, optimize=ES)
                + np.einsum("prst,qrst->pq", r.R, dgam, optimize=ES))
    fold = (np.einsum("rs,rpsq->pq", dd, r.L, optimize=ES)
            + np.einsum("rs,rqsp->pq", dd, r.L, optimize=ES)
            + np.einsum("rs,rpsq->pq", d, dl, optimize=ES)
            + np.einsum("rs,rqsp->pq", d, dl, optimize=ES))
    dip[:, :no] += fold[:, :no]
    return -0.5 * dip


def perturbed_zvector_rhs(r, dip, df, dl):
    no = r.no
    rhs = dip[:no, no:] - dip[no:, :no].T          # (i, a)
    z = r.z
    # (dL_ajib + dL_abij) z_jb
    rhs -= np.einsum("ajib,jb->ia", dl[no:, :no, :no, no:], z, optimize=ES)
    rhs -= np.einsum("abij,jb->ia", dl[no:, no:, :no, :no], z, optimize=ES)
    rhs -= np.einsum("ab,ib->ia", df[no:, no:], z, optimize=ES)
    rhs += np.einsum("ij,ja->ia", df[:no, :no], z, optimize=ES)
    return rhs




# ---------------------------------------------------------------------------
# skeleton Lagrangian and the orbital-response carriers (per X)
# ---------------------------------------------------------------------------

def skeleton_lagrangian(r, x):
    no = r.no
    fx1, sx1 = r.fx[x], r.sx[x]
    ec, ee = r.erix_c[x], r.erix_e[x]
    drel, imat = r.drel, r.ip
    ip = fx1 @ drel.T + fx1.T @ drel + sx1 @ imat.T + sx1.T @ imat
    ip += 4 * np.einsum("prst,qrst->pq", ec, r.gam, optimize=ES)
    fold = (2 * np.einsum("rs,rpsq->pq", drel, ee, optimize=ES)
            - np.einsum("rs,rpqs->pq", drel, ee, optimize=ES)
            + 2 * np.einsum("rs,rqsp->pq", drel, ee, optimize=ES)
            - np.einsum("rs,rqps->pq", drel, ee, optimize=ES))
    ip[:, :no] += fold[:, :no]
    ip *= -0.5
    xov = ip.T - ip
    i2 = ip.copy()
    i2[no:, :no] = ip[:no, no:].T
    return ip, xov, i2


def pair_rotation_augment(r, ip, xov, i2):
    nf, no, n = r.nf, r.no, r.n
    L, e = r.L, r.eps
    pof = np.zeros((no, no))
    for q in range(nf, no):
        for p in range(nf):
            pof[p, q] = (ip[p, q] - ip[q, p]) / (e[p] - e[q])
            pof[q, p] = pof[p, q]
    pf = np.zeros((n, n))
    pf[:no, :no] = pof
    xt = xov.copy()
    fold_v = (np.einsum("kl,kali->ai", pof, L[:no, no:, :no, :no], optimize=ES)
              + np.einsum("kl,kila->ai", pof, L[:no, :no, :no, no:], optimize=ES))
    xt[no:, :no] += 0.5 * fold_v
    it = i2.copy()
    fold_o = (np.einsum("kl,kplq->pq", pof, L[:no, :no, :no, :no], optimize=ES)
              + np.einsum("kl,kqlp->pq", pof, L[:no, :no, :no, :no], optimize=ES))
    it[:no, :no] -= pof * e[None, :no] + 0.5 * fold_o
    return xt, it, pf


def perturbed_response(r):
    """Transcription of mp2_perturbed_response, dense, one X at a time."""
    no, nf, nv, n, na = r.no, r.nf, r.nv, r.n, r.no - r.nf
    npert = 3 * r.natm
    e = r.eps
    dt2, ddrel, di, rhs = [], [], [], []
    for y in range(npert):
        u = full_u(r, y)
        df = perturbed_fock(r, y, u)
        deri_e = perturbed_eri(r.erix_e[y], u, r.eri)
        deri_c = deri_e if r.corr == "exact" else perturbed_eri(r.erix_c[y], u, r.R)
        dl = 2 * deri_e - deri_e.transpose(0, 1, 3, 2)
        ta = perturbed_t2(r, deri_c[nf:no, nf:no, no:, no:], df)
        ddg = perturbed_onepdm(r, ta)
        dgam = cumulant(ta, nf, no, n)
        dip = perturbed_lagrangian(r, df, deri_c, dl, r.d0, ddg, dgam)
        rhs1 = perturbed_zvector_rhs(r, dip, df, dl)         # (i, a)
        if nf > 0:
            dpco = np.zeros((nf, na))
            for i in range(na):
                for c in range(nf):
                    acc = dip[c, nf + i] - dip[nf + i, c]
                    acc -= np.dot(df[c, :nf], r.pco[:, i])
                    acc += np.dot(r.pco[c, :], df[nf:no, nf + i])
                    dpco[c, i] = acc / (e[c] - e[nf + i])
            L = r.L
            for a in range(nv):
                for i in range(no):
                    acc = 0.0
                    for c in range(nf):
                        for j in range(na):
                            acc += (dpco[c, j] * (L[no + a, nf + j, i, c]
                                                  + L[no + a, c, i, nf + j])
                                    + r.pco[c, j] * (dl[no + a, nf + j, i, c]
                                                     + dl[no + a, c, i, nf + j]))
                    rhs1[i, a] += acc
            for i in range(na):
                for c in range(nf):
                    ddg[c, nf + i] += dpco[c, i]
                    ddg[nf + i, c] += dpco[c, i]
        rhs.append(rhs1)
        dip2 = perturbed_lagrangian(r, df, deri_c, dl, r.drel, ddg, dgam)
        di.append(dip2)
        ddrel.append(ddg)
        dt2.append(ta)
    # one batched solve:   A x = -rhs   (response added to the vo block)
    for y in range(npert):
        x = -(r.Ainv @ rhs[y].ravel()).reshape(no, nv)       # (i, a)
        zm = np.zeros((n, n))
        zm[:no, no:] = x
        zm[no:, :no] = x.T
        ddrel[y] = ddrel[y] + zm
        di[y] = di[y] - 0.5 * e[:, None] * (zm + zm.T)
        fold = (np.einsum("rs,rpsq->pq", zm, r.L, optimize=ES)
                + np.einsum("rs,rqsp->pq", zm, r.L, optimize=ES))
        di[y][:, :no] -= 0.5 * fold[:, :no]
    r.dt2, r.ddrel, r.di = dt2, ddrel, di
    return r


def orbital_response(r, carriers):
    npert = 3 * r.natm
    no = r.no
    orb = np.zeros((npert, npert))
    for y in range(npert):
        for x in range(npert):
            xt, it, pf = carriers[x]
            acc = 2.0 * np.sum(r.mo1[y][no:, :no] * xt[no:, :no])
            acc += np.sum(r.sx[y] * it)
            if pf is not None:
                acc += np.sum(pf * r.fx[y])
            orb[x, y] = acc
    return orb


def carriers_all(r):
    out = []
    for x in range(3 * r.natm):
        ip, xov, i2 = skeleton_lagrangian(r, x)
        if r.nf > 0:
            xt, it, pf = pair_rotation_augment(r, ip, xov, i2)
        else:
            xt, it, pf = xov, i2, None
        out.append((xt, it, pf))
    return out


# ---------------------------------------------------------------------------
# pass 1, by differencing the analytic skeleton gradient at frozen densities
# ---------------------------------------------------------------------------

def skeleton_vector(r, atoms, parts=("ref", "corr")):
    """E^X at a displaced geometry with C, D_rel, Gamma and I' held fixed.

    ``ref`` is the reference operator's share -- D_rel against the core
    Hamiltonian and the occupied-density Fock fold, and I' against S -- and
    ``corr`` the correlation functional's. Differencing it along Y is pass 1.
    """
    mol, auxmol = grad_ref.build(atoms, r.basis, r.auxbasis)
    d1 = first_derivatives(mol, auxmol)
    no, nf = r.no, r.nf
    cocc = r.C[:, :no]
    pocc = cocc @ cocc.T
    dao = r.C @ r.drel @ r.C.T
    nx = 3 * mol.natm
    g = np.zeros(nx)
    if "corr" in parts and r.corr == "ri":
        three = df.incore.aux_e2(mol, auxmol, intor="int3c2e", aosym="s1")
        three = three.reshape(mol.nao, mol.nao, -1)
        jinv = grad_ref.metric_inverse_sqrt(auxmol.intor("int2c2e"))
        jinv = jinv @ jinv
        imo = np.einsum("up,uvP,vq->pqP", r.C, three, r.C, optimize=ES)
        cov = (imo @ jinv)[nf:no, no:]
        u = 2 * r.t2 - r.t2.transpose(0, 1, 3, 2)
        gp_ = np.einsum("ijab,jbP->iaP", u, cov, optimize=ES)
        gq = np.einsum("iaP,iaQ->PQ", gp_, cov, optimize=ES)
    for a in range(mol.natm):
        for c in range(3):
            x = 3 * a + c
            e = d1["E"][a][c]
            val = 0.0
            if "ref" in parts:
                val += np.sum(dao * d1["H"][a][c])
                val += 2 * np.einsum("pqls,pq,ls->", e, dao, pocc, optimize=ES)
                val -= np.einsum("plsq,pq,ls->", e, dao, pocc, optimize=ES)
                val += np.sum(r.ip * (r.C.T @ d1["S"][a][c] @ r.C))
            if "corr" in parts:
                if r.corr == "ri":
                    ix = np.einsum("up,uvP,vq->pqP", r.C[:, nf:no], d1["I"][a][c],
                                   r.C[:, no:], optimize=ES)
                    val += 4 * np.sum(gp_ * ix) - 2 * np.sum(gq * d1["J"][a][c])
                else:
                    val += np.sum(r.gam * phys(ao_to_mo4(e, r.C)))
            g[x] = val
    return g


def displaced(atoms, x, h):
    a, c = divmod(x, 3)
    out = [list(t) for t in atoms]
    out[a][1 + c] += h
    return out


def pass1(r, h=1e-3, parts=("ref", "corr")):
    n3 = 3 * r.natm
    p1 = np.zeros((n3, n3))
    for y in range(n3):
        gp1 = skeleton_vector(r, displaced(r.atoms, y, h), parts)
        gm1 = skeleton_vector(r, displaced(r.atoms, y, -h), parts)
        gp2 = skeleton_vector(r, displaced(r.atoms, y, 2 * h), parts)
        gm2 = skeleton_vector(r, displaced(r.atoms, y, -2 * h), parts)
        p1[:, y] = (-gp2 + 8 * gp1 - 8 * gm1 + gm2) / (12 * h)
    return p1


def correlation_hessian(r, h=1e-3):
    """H[X, Y] over 3N x 3N, correlation block only."""
    first_order_skeletons(r, first_derivatives(r.mol, r.auxmol))
    solve_mo1(r)
    perturbed_response(r)
    carriers = carriers_all(r)
    orb = orbital_response(r, carriers)
    n3 = 3 * r.natm
    resp = np.zeros((n3, n3))
    for y in range(n3):
        dg = cumulant(r.dt2[y], r.nf, r.no, r.n)
        for x in range(n3):
            resp[x, y] = (np.sum(r.ddrel[y] * r.fx[x])
                          + np.sum(dg * r.erix_c[x])
                          + np.sum(r.di[y] * r.sx[x]))
    p1 = pass1(r, h)
    r.parts = {"pass1": p1, "orb": orb, "resp": resp}
    return p1 + orb + resp


# ===========================================================================
# the lean path: no four-index MO tensor anywhere
# ===========================================================================
#
# The dense path above holds <pq|rs> (reference) and R (correlation) as n^4
# arrays because that is the shape pycc's formulas are written in. Every use
# of them reduces to one of two things, and the lean path is those two things:
#
#   * the reference's two-electron operator  T1(D)_pq = sum_rs D_rs L_rpsq,
#     L = 2<pq|rs> - <pq|sr>, which is 2J[D] - K[D]^T carried to the MO basis,
#     and its first derivative at fixed MO coefficients, T1x(D)[X]; and
#   * the correlation tensor R = (pa|P) J^-1 (P|jb), which is only ever
#     contracted against the cumulant u_ijab, so it appears as three-index
#     objects: Gamma^P_ia = u c^P and its relatives.

class Lean:
    """Reference operators and the 3-index correlation objects."""

    def __init__(self, r, d1):
        self.r = r
        self.d1 = d1
        self.eri_ao = r.mol.intor("int2e").reshape((r.nao,) * 4)

    # ---- the reference's operator ----------------------------------------
    def _t1(self, dao, e):
        j = np.einsum("pqls,ls->pq", e, dao, optimize=ES)
        k = np.einsum("plqs,ls->pq", e, dao, optimize=ES)
        return 2 * j - k.T

    def T1(self, d):
        c = self.r.C
        return c.T @ self._t1(c @ d @ c.T, self.eri_ao) @ c

    def T1x(self, d, x):
        c = self.r.C
        a, k = divmod(x, 3)
        return c.T @ self._t1(c @ d @ c.T, self.d1["E"][a][k]) @ c

    def T1x_stack(self, d):
        return [self.T1x(d, x) for x in range(3 * self.r.natm)]

    def DT1(self, d, t1d, t1x_y, u):
        """Total derivative of T1(d) for a fixed MO matrix d, along y."""
        return (t1x_y + self.T1(u @ d + d @ u.T) + u.T @ t1d + t1d @ u)


def embed_oo(r, s):
    m = np.zeros((r.n, r.n))
    m[:r.no, :r.no] = s[:r.no, :r.no]
    return m


def embed_vo(r, u):
    """The virtual-occupied block of u, alone, as an n x n matrix."""
    m = np.zeros((r.n, r.n))
    m[r.no:, :r.no] = u[r.no:, :r.no]
    return m


def embed_z(r, z):
    m = np.zeros((r.n, r.n))
    m[:r.no, r.no:] = z
    m[r.no:, :r.no] = z.T
    return m


def gp(r, a, g):
    """sum_P A^P_pa G^P_ia (-> [p,i]) and sum_P A^P_ip G^P_ia (-> [p,a]).

    ``a`` is an (n, n, naux) integral-like tensor, ``g`` a (nact, nv, naux)
    density-like one. This is the whole of every two-particle Lagrangian term.
    """
    no, nf, n = r.no, r.nf, r.n
    out = np.zeros((n, n))
    out[:, nf:no] = np.einsum("paP,iaP->pi", a[:, no:, :], g, optimize=ES)
    out[:, no:] = np.einsum("ipP,iaP->pa", a[nf:no, :, :], g, optimize=ES)
    return out


def theta(r, a_ov, j1, gam3, u):
    """Theta^Q_ia = sum_jb u_ijab A^Q_jb - sum_R J1_QR Gamma^R_ia."""
    return (np.einsum("ijab,jbQ->iaQ", u, a_ov, optimize=ES)
            - np.einsum("QR,iaR->iaQ", j1, gam3, optimize=ES))


def corr_objects(r):
    """Gamma^P_ia, gamma_PQ, c_ov, the cumulant u."""
    no, nf = r.no, r.nf
    r.u = 2 * r.t2 - r.t2.transpose(0, 1, 3, 2)
    r.c_ov = r.cmo[nf:no, no:, :]
    r.gam3 = np.einsum("ijab,jbP->iaP", r.u, r.c_ov, optimize=ES)
    r.gamma2 = np.einsum("iaP,iaQ->PQ", r.gam3, r.c_ov, optimize=ES)
    return r


def mo3(r, ao3, c=None):
    """(nao, nao, naux) -> (n, n, naux) with the reference's C."""
    c = r.C if c is None else c
    return np.einsum("up,uvP,vq->pqP", c, ao3, c, optimize=ES)


def setup_lean(r):
    """Fitted objects and amplitudes with no four-index tensor of any kind."""
    no, nf = r.no, r.nf
    r.three, r.metric, r.jinv, r.imo, r.cmo = ri_objects(r)
    r.corr = "ri"
    kov = np.einsum("iaP,jbP->ijab", r.imo[nf:no, no:], r.cmo[nf:no, no:],
                    optimize=ES)
    eo, ev = r.eps[nf:no], r.eps[no:]
    d = (eo[:, None, None, None] + eo[None, :, None, None]
         - ev[None, None, :, None] - ev[None, None, None, :])
    r.t2 = kov / d
    return r


def lean_unperturbed(r, d1):
    """The unperturbed ladder -- amplitudes to relaxed density and I' -- from
    three-index objects and the reference operator only."""
    no, nf, nv, n = r.no, r.nf, r.nv, r.n
    lean = Lean(r, d1)
    r.lean = lean
    corr_objects(r)
    r.d0 = unrelaxed_density(r.t2, nf, no, n)
    ip0 = lean_lagrangian_unperturbed(r, r.d0)
    r.pco = None
    if nf > 0:
        r.pco = np.zeros((nf, no - nf))
        for i in range(no - nf):
            for c in range(nf):
                r.pco[c, i] = (ip0[c, nf + i] - ip0[nf + i, c]) \
                    / (r.eps[c] - r.eps[nf + i])
    r.Ainv_lean = np.linalg.inv(lean_orbital_hessian(r))
    x = (ip0[:no, no:] - ip0[no:, :no].T).copy()
    if nf > 0:
        ps = np.zeros((n, n))
        ps[nf:no, :nf] = r.pco.T
        ps[:nf, nf:no] = r.pco
        r.psym = ps
        r.T1_p = lean.T1(ps)
        r.T1x_p = lean.T1x_stack(ps)
        x += r.T1_p[no:, :no].T
    r.z = (r.Ainv_lean @ x.ravel()).reshape(no, nv)
    drel = r.d0.copy()
    if nf > 0:
        drel[nf:no, :nf] += r.pco.T
        drel[:nf, nf:no] += r.pco
    drel[no:, :no] -= r.z.T
    drel[:no, no:] -= r.z
    r.drel = drel
    r.ip = lean_lagrangian_unperturbed(r, drel)
    # T1 of the unperturbed densities, and their first derivatives
    r.T1_d0 = lean.T1(r.d0)
    r.T1_drel = lean.T1(r.drel)
    r.T1x_d0 = lean.T1x_stack(r.d0)
    r.T1x_drel = lean.T1x_stack(r.drel)
    r.zsym = embed_z(r, r.z)
    r.T1_z = lean.T1(r.zsym)
    r.T1x_z = lean.T1x_stack(r.zsym)
    return lean


def lean_lagrangian_unperturbed(r, d):
    """I' from three-index objects only (checks the dense lagrangian)."""
    ip = r.eps[:, None] * (d + d.T)
    ip += 4 * gp(r, r.imo, r.gam3)
    t = r.lean.T1(d)
    ip[:, :r.no] += (t + t.T)[:, :r.no]
    return -0.5 * ip


def lean_skeletons(r, d1):
    """sx, fx, and the three-index first-derivative integrals, per X."""
    no, n = r.no, r.n
    docc = np.zeros((n, n))
    docc[:no, :no] = np.eye(no)
    r.fx, r.sx, r.Ix, r.Jx = [], [], [], []
    for a in range(r.natm):
        for c in range(3):
            x = 3 * a + c
            r.sx.append(r.C.T @ d1["S"][a][c] @ r.C)
            r.fx.append(r.C.T @ d1["H"][a][c] @ r.C + r.lean.T1x(docc, x))
            r.Ix.append(mo3(r, d1["I"][a][c]))
            r.Jx.append(d1["J"][a][c])
    r.Ix_ov = np.array([i[r.nf:no, no:, :] for i in r.Ix])
    r.Jx_all = np.array(r.Jx)
    return r


def lean_orbital_hessian(r):
    no, nv = r.no, r.nv
    a = np.zeros((no * nv, no * nv))
    g = (r.eps[None, no:] - r.eps[:no, None]).ravel()
    for k in range(no * nv):
        z = np.zeros(no * nv)
        z[k] = 1.0
        t = r.lean.T1(embed_z(r, z.reshape(no, nv)))
        a[:, k] = t[no:, :no].T.ravel()
        a[k, k] += g[k]
    return a


def lean_solve_mo1(r):
    no, nv, n = r.no, r.nv, r.n
    ainv = r.Ainv_lean
    r.mo1 = []
    for x in range(3 * r.natm):
        sx, fx = r.sx[x], r.fx[x]
        t = r.lean.T1(embed_oo(r, sx))
        rhs = fx[no:, :no] - sx[no:, :no] * r.eps[None, :no] - t[no:, :no]
        u = -(ainv @ rhs.T.ravel()).reshape(no, nv).T
        m = np.zeros((n, no))
        m[no:] = u
        m[:no] = -0.5 * sx[:no, :no]
        r.mo1.append(m)


def lean_full_u(r, x):
    no, nf, n = r.no, r.nf, r.n
    sx, fx, mo1, e = r.sx[x], r.fx[x], r.mo1[x], r.eps
    u = np.zeros((n, n))
    u[:, :no] = mo1
    u[no:, no:] = -0.5 * sx[no:, no:]
    u[:no, no:] = -sx[:no, no:] - mo1[no:, :].T
    if nf > 0:
        uvo = embed_vo(r, u)
        t = r.lean.T1(-embed_oo(r, sx) + uvo + uvo.T)
        for i in range(nf, no):
            for c in range(nf):
                acc = fx[c, i] - sx[c, i] * e[i] + t[c, i]
                u[c, i] = -acc / (e[c] - e[i])
                u[i, c] = -sx[c, i] - u[c, i]
    return u


def lean_perturbed_fock(r, x, u):
    sx, fx, e = r.sx[x], r.fx[x], r.eps
    uvo = embed_vo(r, u)
    return (fx + u * e[:, None] + u.T * e[None, :]
            + r.lean.T1(-embed_oo(r, sx) + uvo + uvo.T))


def lean_dI(r, x, u):
    """Total derivative of the MO three-index integrals along x."""
    return (r.Ix[x] + np.einsum("tp,tqP->pqP", u, r.imo, optimize=ES)
            + np.einsum("tq,ptP->pqP", u, r.imo, optimize=ES))


def lean_dK(r, di_ov, dj):
    """d_X of (ia|jb)_RI, as dk[i, j, a, b]."""
    c = r.c_ov
    return (np.einsum("iaP,jbP->ijab", di_ov, c, optimize=ES)
            + np.einsum("iaP,jbP->ijab", c, di_ov, optimize=ES)
            - np.einsum("iaP,PQ,jbQ->ijab", c, dj, c, optimize=ES))


def lean_corr_lagrangian_term(r, a_full, j1, d_gam3=None):
    """4 * [gp(A, Gamma) + gp(c, Theta(A, J1)) (+ gp(I, dGamma))].

    The correlation functional's whole contribution to a Lagrangian-shaped
    matrix: a first-derivative (or total-derivative) integral tensor ``A``
    with its metric derivative ``J1``, and optionally the response of the
    three-index density.
    """
    no, nf = r.no, r.nf
    th = theta(r, a_full[nf:no, no:, :], j1, r.gam3, r.u)
    out = gp(r, a_full, r.gam3) + gp(r, r.cmo, th)
    if d_gam3 is not None:
        out = out + gp(r, r.imo, d_gam3)
    return 4 * out


def lean_perturbed_lagrangian(r, df, di_full, dj, d_gam3, d, dd, dt1):
    dsym = d + d.T
    dip = df @ dsym
    dip += r.eps[:, None] * (dd + dd.T)
    dip += lean_corr_lagrangian_term(r, di_full, dj, d_gam3)
    t = r.lean.T1(dd) + dt1
    dip[:, :r.no] += (t + t.T)[:, :r.no]
    return -0.5 * dip


def lean_zvector_rhs(r, dip, df, dt1z):
    no = r.no
    rhs = dip[:no, no:] - dip[no:, :no].T
    rhs -= dt1z[no:, :no].T
    rhs -= r.z @ df[no:, no:].T
    rhs += df[:no, :no] @ r.z
    return rhs


def lean_response(r):
    no, nf, nv, n, na = r.no, r.nf, r.nv, r.n, r.no - r.nf
    npert = 3 * r.natm
    e, lean = r.eps, r.lean
    dt2, ddrel, di, rhs = [], [], [], []
    pair = np.zeros((npert, npert))
    for y in range(npert):
        u = lean_full_u(r, y)
        df = lean_perturbed_fock(r, y, u)
        dI = lean_dI(r, y, u)
        dJ = r.Jx[y]
        dk = lean_dK(r, dI[nf:no, no:, :], dJ)
        ta = perturbed_t2(r, dk, df)
        ddg = perturbed_onepdm(r, ta)
        du = 2 * ta - ta.transpose(0, 1, 3, 2)
        dG = np.einsum("ijab,jbP->iaP", du, r.c_ov, optimize=ES)
        # the pair sum 2 sum du_ijab R^(X)_iajb, now, while dG is in hand: it
        # does not depend on the Z-vector, so nothing o*v*naux outlives this y
        dgt = np.einsum("iaP,iaQ->PQ", dG, r.c_ov, optimize=ES)
        pair[:, y] = (4 * np.einsum("iaP,xiaP->x", dG, r.Ix_ov, optimize=ES)
                      - 2 * np.einsum("PQ,xPQ->x", dgt, r.Jx_all, optimize=ES))
        dt1_d0 = lean.DT1(r.d0, r.T1_d0, r.T1x_d0[y], u)
        dt1_z = lean.DT1(r.zsym, r.T1_z, r.T1x_z[y], u)
        dip = lean_perturbed_lagrangian(r, df, dI, dJ, dG, r.d0, ddg, dt1_d0)
        rhs1 = lean_zvector_rhs(r, dip, df, dt1_z)
        if nf > 0:
            dpco = np.zeros((nf, na))
            for i in range(na):
                for c in range(nf):
                    acc = dip[c, nf + i] - dip[nf + i, c]
                    acc -= np.dot(df[c, :nf], r.pco[:, i])
                    acc += np.dot(r.pco[c, :], df[nf:no, nf + i])
                    dpco[c, i] = acc / (e[c] - e[nf + i])
            dps = np.zeros((n, n))
            dps[nf:no, :nf] = dpco.T
            dps[:nf, nf:no] = dpco
            dt1_p = lean.DT1(r.psym, r.T1_p, r.T1x_p[y], u)
            rhs1 += (lean.T1(dps) + dt1_p)[no:, :no].T
            ddg = ddg + dps
        rhs.append(rhs1)
        dt1_dr = lean.DT1(r.drel, r.T1_drel, r.T1x_drel[y], u)
        di.append(lean_perturbed_lagrangian(r, df, dI, dJ, dG, r.drel, ddg, dt1_dr))
        ddrel.append(ddg)
        dt2.append(ta)
    for y in range(npert):
        x = -(r.Ainv_lean @ rhs[y].ravel()).reshape(no, nv)
        zm = embed_z(r, x)
        ddrel[y] = ddrel[y] + zm
        di[y] = di[y] - 0.5 * e[:, None] * (zm + zm.T)
        t = lean.T1(zm)
        di[y][:, :no] -= 0.5 * (t + t.T)[:, :no]
    r.dt2, r.ddrel, r.di, r.pair = dt2, ddrel, di, pair
    return r


def lean_carriers(r):
    no, nf, n = r.no, r.nf, r.n
    out = []
    for x in range(3 * r.natm):
        fx1, sx1 = r.fx[x], r.sx[x]
        ip = (fx1 @ r.drel.T + fx1.T @ r.drel + sx1 @ r.ip.T + sx1.T @ r.ip)
        ip += lean_corr_lagrangian_term(r, r.Ix[x], r.Jx[x])
        t = r.T1x_drel[x]
        ip[:, :no] += (t + t.T)[:, :no]
        ip *= -0.5
        xov = ip.T - ip
        i2 = ip.copy()
        i2[no:, :no] = ip[:no, no:].T
        if nf > 0:
            pof = np.zeros((no, no))
            for q in range(nf, no):
                for p in range(nf):
                    pof[p, q] = (ip[p, q] - ip[q, p]) / (r.eps[p] - r.eps[q])
                    pof[q, p] = pof[p, q]
            pf = np.zeros((n, n))
            pf[:no, :no] = pof
            tp = r.lean.T1(pf)
            xt = xov.copy()
            xt[no:, :no] += tp[no:, :no]
            it = i2.copy()
            it[:no, :no] -= pof * r.eps[None, :no] + tp[:no, :no]
            out.append((xt, it, pf))
        else:
            out.append((xov, i2, None))
    return out


def lean_correlation_hessian(r, d1, pass1_matrix):
    """pass 1 + orbital response + density response, no n^4 anywhere.

    Expects ``setup_lean``, ``lean_unperturbed`` and ``lean_skeletons`` to have
    run: pass 1 needs the first two and the fitted share of it the third.
    """
    lean_solve_mo1(r)
    lean_response(r)
    carriers = lean_carriers(r)
    orb = orbital_response(r, carriers)
    n3 = 3 * r.natm
    resp = r.pair.copy()
    for y in range(n3):
        for x in range(n3):
            resp[x, y] += (np.sum(r.ddrel[y] * r.fx[x])
                           + np.sum(r.di[y] * r.sx[x]))
    r.lean_parts = {"orb": orb, "resp": resp}
    return pass1_matrix + orb + resp


# ---------------------------------------------------------------------------
# pass 1 for the fitted part, analytically
# ---------------------------------------------------------------------------

def second_derivatives_3c2c(mol, auxmol):
    """The libfint/libcint second-derivative integrals, 9 components each."""
    nao, naux = mol.nao, auxmol.nao
    def aux3(name):
        return df.incore.aux_e2(mol, auxmol, intor=name, aosym="s1",
                                comp=9).reshape(3, 3, nao, nao, naux)
    return {
        "ipip1": aux3("int3c2e_ipip1"), "ipvip1": aux3("int3c2e_ipvip1"),
        "ip1ip2": aux3("int3c2e_ip1ip2"), "ipip2": aux3("int3c2e_ipip2"),
        "c_ipip1": auxmol.intor("int2c2e_ipip1", comp=9).reshape(3, 3, naux, naux),
        "c_ip1ip2": auxmol.intor("int2c2e_ip1ip2", comp=9).reshape(3, 3, naux, naux),
    }


def contract_second_3c2c(mol, auxmol, gam_ao, gamma2, d2=None):
    """4 Gamma^P_mn (mn|P)^{XY} - 2 gamma_PQ (P|Q)^{XY}, all (X, Y).

    ``gam_ao`` is (nao, nao, naux), ``gamma2`` (naux, naux). Both derivatives
    land on any of the three centres of a three-index integral (two of a
    two-index one); each ordered pair of centres is one libcint array,
    selected by the atoms that own the two functions. A derivative with respect
    to an electron coordinate is minus the nuclear one, so the product of two
    is plus.
    """
    if d2 is None:
        d2 = second_derivatives_3c2c(mol, auxmol)
    natm = mol.natm
    m = owners(mol)
    q = owners(auxmol)
    # symmetric image of the density: (mn|P) is symmetric in m, n
    g = gam_ao
    h = np.zeros((3 * natm, 3 * natm))
    for a in range(natm):
        for b in range(natm):
            blk = np.zeros((3, 3))
            # (mu, mu): both derivatives on the first function
            w = m[a][:, None] * m[b][:, None]
            blk += np.einsum("xyuvP,uvP->xy", d2["ipip1"], g * w[:, :, None])
            # (nu, nu)
            w = m[a][None, :] * m[b][None, :]
            blk += np.einsum("xyvuP,uvP->xy", d2["ipip1"], g * w[:, :, None])
            # (mu, nu) and (nu, mu)
            w = m[a][:, None] * m[b][None, :]
            blk += np.einsum("xyuvP,uvP->xy", d2["ipvip1"], g * w[:, :, None])
            w = m[a][None, :] * m[b][:, None]
            blk += np.einsum("yxuvP,uvP->xy", d2["ipvip1"], g * w[:, :, None])
            # (mu, P) and (P, mu)
            w = m[a][:, None, None] * q[b][None, None, :]
            blk += np.einsum("xyuvP,uvP->xy", d2["ip1ip2"], g * w)
            w = q[a][None, None, :] * m[b][:, None, None]
            blk += np.einsum("yxuvP,uvP->xy", d2["ip1ip2"], g * w)
            # (nu, P) and (P, nu)
            w = m[a][None, :, None] * q[b][None, None, :]
            blk += np.einsum("xyvuP,uvP->xy", d2["ip1ip2"], g * w)
            w = q[a][None, None, :] * m[b][None, :, None]
            blk += np.einsum("yxvuP,uvP->xy", d2["ip1ip2"], g * w)
            # (P, P)
            w = q[a][None, None, :] * q[b][None, None, :]
            blk += np.einsum("xyuvP,uvP->xy", d2["ipip2"], g * w)
            blk4 = 4.0 * blk
            # two-index: (P,P) (Q,Q) (P,Q) (Q,P)
            c = np.zeros((3, 3))
            w = q[a][:, None] * q[b][:, None]
            c += np.einsum("xyPQ,PQ->xy", d2["c_ipip1"], gamma2 * w)
            w = q[a][None, :] * q[b][None, :]
            c += np.einsum("xyQP,PQ->xy", d2["c_ipip1"], gamma2 * w)
            w = q[a][:, None] * q[b][None, :]
            c += np.einsum("xyPQ,PQ->xy", d2["c_ip1ip2"], gamma2 * w)
            w = q[a][None, :] * q[b][:, None]
            c += np.einsum("yxQP,PQ->xy", d2["c_ip1ip2"], gamma2 * w)
            h[3 * a:3 * a + 3, 3 * b:3 * b + 3] = blk4 - 2.0 * c
    return h


def fitted_pass1_analytic(r, d1):
    """The fitted correlation's share of pass 1: second derivatives of the
    3c/2c integrals against Gamma and gamma, plus the cross terms from the
    first-order change of c = J^-1 (ia|P) at fixed cumulant."""
    no, nf = r.no, r.nf
    n3 = 3 * r.natm
    gam_ao = np.einsum("ui,iaP,va->uvP", r.C[:, nf:no], r.gam3, r.C[:, no:],
                       optimize=ES)
    second = contract_second_3c2c(r.mol, r.auxmol, gam_ao, r.gamma2)
    cross = np.zeros((n3, n3))
    for y in range(n3):
        iy = r.Ix[y][nf:no, no:, :]
        jy = r.Jx[y]
        th = theta(r, iy, jy, r.gam3, r.u)
        dgam3 = np.einsum("PQ,iaQ->iaP", r.jinv, th, optimize=ES)
        dc = np.einsum("QR,iaR->iaQ", r.jinv,
                       iy - np.einsum("RS,iaS->iaR", jy, r.c_ov, optimize=ES),
                       optimize=ES)
        dgamma = (np.einsum("iaP,iaQ->PQ", dgam3, r.c_ov, optimize=ES)
                  + np.einsum("iaP,iaQ->PQ", r.gam3, dc, optimize=ES))
        for x in range(n3):
            ix = r.Ix[x][nf:no, no:, :]
            cross[x, y] = 4 * np.sum(dgam3 * ix) - 2 * np.sum(dgamma * r.Jx[x])
    return second, cross


# ---------------------------------------------------------------------------
# intermediates, each against a finite difference
# ---------------------------------------------------------------------------

def five_point(fun, h):
    """d/dh of an array-valued function at 0."""
    f = [fun(s * h) for s in (2, 1, -1, -2)]
    return (-f[0] + 8 * f[1] - 8 * f[2] + f[3]) / (12 * h)


def inverse_sqrt_pair(metric, threshold=1e-10):
    w, v = np.linalg.eigh(metric)
    keep = w > threshold
    return v[:, keep], w[keep]


def moved(r, y, h, u):
    """Geometry displaced by h along y with orbitals C0 (1 + h u).

    The orbitals are *not* re-converged and not re-orthonormalised: u is the
    first-order rotation the analytic formulas used, and a finite difference at
    fixed u differentiates exactly what those formulas claim to, one gauge,
    one quantity at a time.
    """
    no, nf = r.no, r.nf
    at = displaced(r.atoms, y, h)
    mol, auxmol = grad_ref.build(at, r.basis, r.auxbasis)
    c = r.C + h * r.C @ u
    three = df.incore.aux_e2(mol, auxmol, intor="int3c2e", aosym="s1")
    three = three.reshape(mol.nao, mol.nao, -1)
    metric = auxmol.intor("int2c2e")
    jm12 = grad_ref.metric_inverse_sqrt(metric)
    jinv = jm12 @ jm12
    imo = np.einsum("up,uvP,vq->pqP", c, three, c, optimize=ES)
    b = np.einsum("pqP,PQ->pqQ", imo, jm12, optimize=ES)
    cmo = imo @ jinv
    hcore = mol.intor("int1e_kin") + mol.intor("int1e_nuc")
    eri = mol.intor("int2e").reshape((mol.nao,) * 4)
    cocc = c[:, :no]
    dm = cocc @ cocc.T
    fock = hcore + 2 * np.einsum("pqls,ls->pq", eri, dm, optimize=ES) \
        - np.einsum("plqs,ls->pq", eri, dm, optimize=ES)
    f = c.T @ fock @ c
    kov = np.einsum("iaP,jbP->ijab", imo[nf:no, no:], cmo[nf:no, no:],
                    optimize=ES)
    # non-canonical amplitudes: fixed point of  D t = K + (off-diagonal f) t
    fo, fv = f[nf:no, nf:no], f[no:, no:]
    do, dv = np.diag(fo), np.diag(fv)
    den = do[:, None, None, None] + do[None, :, None, None] \
        - dv[None, None, :, None] - dv[None, None, None, :]
    fo_off, fv_off = fo - np.diag(do), fv - np.diag(dv)
    t = kov / den
    for _ in range(200):
        new = (kov
               + np.einsum("ac,ijcb->ijab", fv_off, t, optimize=ES)
               + np.einsum("bc,ijac->ijab", fv_off, t, optimize=ES)
               - np.einsum("ik,kjab->ijab", fo_off, t, optimize=ES)
               - np.einsum("jk,ikab->ijab", fo_off, t, optimize=ES)) / den
        if np.abs(new - t).max() < 1e-15:
            t = new
            break
        t = new
    return {"imo": imo, "b": b, "cmo": cmo, "f": f, "kov": kov, "t2": t,
            "metric": metric, "c": c}


def d_jm12(metric, djm):
    """d(J^-1/2) from d(J) by the eigenbasis formula."""
    v, w = inverse_sqrt_pair(metric)
    s = np.sqrt(w)
    m = v.T @ djm @ v
    return v @ (-m / (s[:, None] * s[None, :] * (s[:, None] + s[None, :]))) @ v.T


def intermediate_checks(r, d1, y, h=1e-3):
    """Max deviation of each analytic first-order quantity from its FD."""
    no, nf = r.no, r.nf
    u = lean_full_u(r, y) if hasattr(r, "Ix") else full_u(r, y)
    out = {}
    fd = lambda key: five_point(lambda s: moved(r, y, s, u)[key], h)  # noqa: E731
    # B = (ia|Q) J^-1/2 : the total derivative needs d(J^-1/2)
    di_ = lean_dI(r, y, u)
    jm12 = grad_ref.metric_inverse_sqrt(r.metric)
    djm = d_jm12(r.metric, r.Jx[y])
    db = (np.einsum("pqP,PQ->pqQ", di_, jm12, optimize=ES)
          + np.einsum("pqP,PQ->pqQ", r.imo, djm, optimize=ES))
    out["dB (all MO pairs)"] = abs(db - fd("b")).max()
    out["d(ia|P)"] = abs(di_ - fd("imo")).max()
    dc = np.einsum("PQ,pqQ->pqP", r.jinv,
                   di_ - np.einsum("QR,pqR->pqQ", r.Jx[y], r.cmo, optimize=ES),
                   optimize=ES)
    out["d c = d(J^-1 (ia|P))"] = abs(dc - fd("cmo")).max()
    dk = lean_dK(r, di_[nf:no, no:], r.Jx[y])
    fdk = five_point(lambda s: moved(r, y, s, u)["kov"], h)
    out["d(ia|jb)_RI"] = abs(dk - fdk).max()
    df_ = lean_perturbed_fock(r, y, u)
    out["d f (Fock, moving basis)"] = abs(df_ - fd("f")).max()
    ta = perturbed_t2(r, dk, df_)
    out["d t2 (non-canonical)"] = abs(ta - fd("t2")).max()
    return out


def ao_density_checks(r, d1, y, h=1e-3):
    """Gauge-free AO-basis quantities, differenced across re-converged SCFs.

    The analytic side is the lean path's own first-order objects mapped to the
    AO basis with the orbital rotation: D -> C (dD + U D + D U^T) C^T.
    """
    no, nf = r.no, r.nf
    u = lean_full_u(r, y)
    c = r.C
    out = {}

    def state_at(s):
        rr = make_ref(displaced(r.atoms, y, s), r.basis, r.auxbasis, r.nf)
        setup(rr, "ri")
        solve_unperturbed(rr)
        corr_objects(rr)
        gam_emb = np.zeros((rr.n, rr.n, rr.naux))
        gam_emb[nf:no, no:] = rr.gam3
        # AO-basis Gamma^P, with the phases of this geometry's own orbitals
        g_ao = np.einsum("up,pqP,vq->uvP", rr.C, gam_emb, rr.C, optimize=ES)
        return {"drel": rr.C @ rr.drel @ rr.C.T,
                "ip": rr.C @ ((rr.ip + rr.ip.T) / 2) @ rr.C.T,
                "gam_ao": g_ao,
                "gamma2": rr.gamma2}

    ref = {}
    for s in (2, 1, -1, -2):
        ref[s] = state_at(s * h)
    fd = {k: (-ref[2][k] + 8 * ref[1][k] - 8 * ref[-1][k] + ref[-2][k]) / (12 * h)
          for k in ref[1]}
    # analytic
    dd = c @ (r.ddrel[y] + u @ r.drel + r.drel @ u.T) @ c.T
    out["d D_rel (AO)"] = abs(dd - fd["drel"]).max()
    dis = (r.di[y] + r.di[y].T) / 2
    # the rotation of a non-symmetric I' is u ip + ip u^T on each side
    di_ao = c @ (dis + (u @ r.ip + r.ip @ u.T + (u @ r.ip + r.ip @ u.T).T) / 2) @ c.T
    out["d I' (AO, symmetrised)"] = abs(di_ao - fd["ip"]).max()
    # Gamma^P: u * dc  +  du * c  +  rotation
    dI = lean_dI(r, y, u)
    dc_ov = np.einsum("PQ,iaQ->iaP", r.jinv,
                      dI[nf:no, no:] - np.einsum("QR,iaR->iaQ", r.Jx[y], r.c_ov,
                                                 optimize=ES), optimize=ES)
    dt = r.dt2[y]
    dg_cum = np.einsum("ijab,jbP->iaP", 2 * dt - dt.transpose(0, 1, 3, 2),
                       r.c_ov, optimize=ES)
    dG_tot = dg_cum + np.einsum("ijab,jbP->iaP", r.u, dc_ov, optimize=ES)
    emb = np.zeros((r.n, r.n, r.naux))
    emb[nf:no, no:] = r.gam3
    demb = np.zeros_like(emb)
    demb[nf:no, no:] = dG_tot
    # a density-like object rotates contravariantly: U G + G U^T
    rot = (np.einsum("pt,tqP->pqP", u, emb, optimize=ES)
           + np.einsum("qt,ptP->pqP", u, emb, optimize=ES))
    dg_ao = np.einsum("up,pqP,vq->uvP", c, demb + rot, c, optimize=ES)
    out["d Gamma^P (AO)"] = abs(dg_ao - fd["gam_ao"]).max()
    return out


# ---------------------------------------------------------------------------
# the public entry point and the validation driver
# ---------------------------------------------------------------------------

def ri_mp2_correlation_hessian(atoms, basis, auxbasis, nf=0, h=1e-3):
    """The RI-MP2 correlation Hessian, (3N, 3N), over an exact RHF reference.

    Lean path, analytic fitted pass 1. The reference operator's share of
    pass 1 (core Hamiltonian, occupied-density Fock fold, overlap) is the
    conventional MP2 Hessian's, unchanged, and is differenced here from the
    analytic first derivatives rather than rebuilt from second-derivative
    integrals.
    """
    r = make_ref(atoms, basis, auxbasis, nf)
    setup_lean(r)
    d1 = first_derivatives(r.mol, r.auxmol)
    lean_unperturbed(r, d1)
    lean_skeletons(r, d1)
    p_ref = pass1(r, h, parts=("ref",))
    second, cross = fitted_pass1_analytic(r, d1)
    hess = lean_correlation_hessian(r, d1, p_ref + second + cross)
    return hess, r


def rhf_hessian(r):
    n3 = 3 * r.natm
    return r.mf.Hessian().kernel().transpose(0, 2, 1, 3).reshape(n3, n3)


def total_ri_gradient(atoms, basis, auxbasis, nf):
    mol, auxmol = grad_ref.build(atoms, basis, auxbasis)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-13
    mf.conv_tol_grad = 1e-10
    mf.max_cycle = 200
    mf.kernel()
    return grad_ref.ri_mp2_gradient(mol, auxmol, mf, frozen=nf).ravel()


def fd_total_hessian(atoms, basis, auxbasis, nf, step=2e-3):
    n3 = 3 * len(atoms)
    out = np.zeros((n3, n3))
    for y in range(n3):
        out[:, y] = five_point(
            lambda s: total_ri_gradient(displaced(atoms, y, s), basis, auxbasis, nf),
            step)
    return out


WATER = [["O", 0.0, 0.0, 0.0], ["H", 0.0, -1.4308, 1.1078],
         ["H", 0.0, 1.4308, 1.1078]]
HCN = [["H", 0.0, 0.0, -2.0], ["C", 0.0, 0.0, 0.0], ["N", 0.0, 0.0, 2.2]]
CASES = [
    ("H2O / sto-3g / cc-pvdz-rifit, all electron", WATER, "sto-3g",
     "cc-pvdz-rifit", 0),
    ("H2O / 6-31g / cc-pvdz-rifit, all electron", WATER, "6-31g",
     "cc-pvdz-rifit", 0),
    ("H2O / 6-31g / cc-pvtz-rifit, frozen core 1", WATER, "6-31g",
     "cc-pvtz-rifit", 1),
    ("HCN / sto-3g / cc-pvdz-rifit, frozen core 2 (degenerate pi)", HCN,
     "sto-3g", "cc-pvdz-rifit", 2),
]


def report(name, value, tol, results):
    ok = value < tol
    results.append(ok)
    print(f"   {'ok  ' if ok else 'FAIL'} {name:<58s} {value:10.3e}  (< {tol:.0e})")


def run_case(label, atoms, basis, auxbasis, nf, results):
    print(f"== {label}")
    # -- the unperturbed identity: the pycc-form gradient is the RI-MP2 one
    r = make_ref(atoms, basis, auxbasis, nf)
    setup(r, "ri")
    solve_unperturbed(r)
    first_order_skeletons(r, first_derivatives(r.mol, r.auxmol))
    g = corr_gradient_from_state(r)
    gref = grad_ref.ri_mp2_gradient(r.mol, r.auxmol, r.mf, frozen=nf) \
        - r.mf.nuc_grad_method().kernel()
    report("gradient, D_rel.f + Gamma.R + I'.S  vs  RI-MP2 gradient",
           abs(g - gref).max(), 1e-10, results)

    # -- dense (n^4) and lean (no n^4) Hessians agree
    h_dense = correlation_hessian(r)
    hess, rl = ri_mp2_correlation_hessian(atoms, basis, auxbasis, nf)
    report("lean path vs dense four-index path (correlation block)",
           abs(hess - h_dense).max(), 1e-8, results)
    report("symmetry of the correlation block, |H - H^T|",
           abs(hess - hess.T).max(), 1e-9, results)
    print(f"        correlation block: |H| = {np.linalg.norm(hess):.6f}")

    # -- the fitted pass 1, analytic against differenced
    d1 = first_derivatives(rl.mol, rl.auxmol)
    # (given the dense path's own differenced pass 1, the lean path agrees
    # with it to rounding: everything above is then the same arithmetic)
    h_same = lean_correlation_hessian(rl, d1, r.parts["pass1"])
    report("lean path vs dense, same pass 1 (response + orbital terms only)",
           abs(h_same - h_dense).max(), 1e-11, results)
    second, cross = fitted_pass1_analytic(rl, d1)
    p_fit = pass1(rl, 1e-3, parts=("corr",))
    report("fitted pass 1 (3c/2c second derivatives + cross terms) vs FD",
           abs(second + cross - p_fit).max(), 1e-8, results)

    # -- intermediates, each against its own finite difference
    for y in sorted({1, 3 * len(atoms) - 2}):
        for key, val in intermediate_checks(rl, d1, y).items():
            tol = 5e-8 if key.startswith("d c") else 1e-8
            report(f"y={y}: {key}", val, tol, results)
        # differences across re-converged SCFs carry the SCF's own error,
        # ~1e-10 in the density over a step of 1e-3
        for key, val in ao_density_checks(rl, d1, y).items():
            report(f"y={y}: {key}", val, 1e-7, results)

    # -- the total, against differences of the analytic RI-MP2 gradient
    total = hess + rhf_hessian(rl)
    fd = fd_total_hessian(atoms, basis, auxbasis, nf)
    dev = total - fd
    report("RHF Hessian + correlation block vs FD of the RI-MP2 gradient",
           abs(dev).max(), 1e-6, results)
    print(f"        |H_total| = {np.linalg.norm(total):.6f}   "
          f"|dev| = {np.linalg.norm(dev):.3e}   "
          f"|dev|/|H| = {np.linalg.norm(dev) / np.linalg.norm(total):.3e}   "
          f"|FD - FD^T| = {abs(fd - fd.T).max():.3e}")

    # -- the fitting error itself: RI minus conventional MP2 correlation block
    rx = make_ref(atoms, basis, auxbasis, nf)
    setup(rx, "exact")
    solve_unperturbed(rx)
    h_exact = correlation_hessian(rx)
    diff = hess - h_exact
    print(f"        RI - conventional correlation block: max {abs(diff).max():.3e}, "
          f"|diff|/|H_conv| = {np.linalg.norm(diff) / np.linalg.norm(h_exact):.3e}")
    return rl


def cost_report(r):
    n, o, v = r.n, r.no - r.nf, r.nv
    p = r.naux
    nn = 3 * r.natm
    print("== cost structure for the last case "
          f"(n={n}, active o={o}, v={v}, naux={p}, 3N={nn})")
    rows = [
        ("resident, once", "(pq|P) and c=J^-1(pq|P), all MO pairs", f"2 n^2 P = {2*n*n*p}"),
        ("", "Gamma^P_ia, c_ia, gamma_PQ, t2", f"o v P + P^2 + o^2 v^2 = {o*v*p + p*p + o*o*v*v}"),
        ("", "d0, D_rel, I', z, T1(.) of each", f"~8 n^2 = {8*n*n}"),
        ("resident, 3N", "S^(X), f^(X), U^X (mo1), D_rel^Y, I'^Y", f"5 n^2 x 3N = {5*n*n*nn}"),
        ("", "T1x of D0, D_rel, z (first-derivative Fock folds)",
         f"3 n^2 x 3N = {3*n*n*nn}"),
        ("", "(ia|P)^(X) and (P|Q)^(X), prototype keeps all X",
         f"3N (o v P + P^2) = {nn*(o*v*p+p*p)}"),
        ("per Y, transient", "d(pq|P), all MO pairs (blocks (n,v), (o,n) read)",
         f"n^2 P = {n*n*p}"),
        ("", "d(ia|jb)_RI, d t2", f"2 o^2 v^2 = {2*o*o*v*v}"),
        ("", "dGamma~^P_ia = du.c  and  dgamma~_PQ", f"o v P + P^2 = {o*v*p + p*p}"),
        ("", "U^Y, d f, dI', rhs, T1(.) of perturbed densities", f"~10 n^2 = {10*n*n}"),
        ("per (X, Y)", "pair sum 4 dGamma~.(ia|P)^(X) - 2 dgamma~.(P|Q)^(X)",
         f"flops ~ o v P = {o*v*p}"),
    ]
    for a, b, c in rows:
        print(f"   {a:<18s} {b:<62s} {c}")
    print("   (counts are doubles; nothing four-index is resident, and the only "
          "3N-resident\n    objects carrying P are the first-derivative integrals "
          "the prototype keeps for\n    convenience and a Fortran port recomputes "
          "per atom.)")


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cases", default="",
                    help="comma-separated indices into the case list (default all)")
    args = ap.parse_args(argv)
    chosen = [int(i) for i in args.cases.split(",") if i != ""] or range(len(CASES))
    results = []
    last = None
    for i in chosen:
        last = run_case(*CASES[i], results)
    if last is not None:
        cost_report(last)
    print()
    print(f"{sum(results)}/{len(results)} checks passed")
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
