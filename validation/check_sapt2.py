#!/usr/bin/env python3
"""SAPT2 reference implementation, conventional four-index.

SAPT2 is exactly SAPT0 plus four terms -- Elst12, Exch11, Exch12, Ind22 --
each an intramonomer-correlation correction built from monomer MP2 amplitudes
(sapt2.cc:131-176; Exch-Ind22 is not computed but scaled, ind22.cc:52).

This transcribes psi4's libsapt_solver verbatim, with one substitution: the
density-fitting factors B^P are replaced by an *exact* factorization of the
conventional two-electron integrals (an eigendecomposition of the AO ERI
matrix), so that

    sum_P B^P_pq B^P_rs = (pq|rs)      to machine precision.

Every contraction, dressing column and factor then follows psi4 line by line
-- including the three "dressed" columns appended to each three-index block,
which fold the one-electron and nuclear parts of the intermolecular operator
into the same GEMMs (utils.cc, get_*_ints) -- while the result is what a
conventional four-index implementation must produce. psi4 itself cannot make
these numbers: its closed-shell SAPT is density-fitted by construction.

Checks that need no external SAPT2 reference (the plan's section 0):
  * zeroing the monomer amplitudes makes all four terms vanish identically,
    so SAPT2 falls back to SAPT0 bit for bit;
  * the monomer MP2 correlation energy recovered from the amplitudes equals
    PySCF's MP2 on the same ghosted monomer to machine precision;
  * the dressed-contraction machinery reproduces Elst10 as 4 diagAA.diagBB
    (elst10.cc:37) against the independent SAPT0 formula;
  * the exchange-induction potential vAR assembled here equals minus the
    USAPT0-factorized potential of check_sapt0.py.

Every formula is cited to psi4 master (2026-08-12) at file:line.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_sapt0 as s0  # noqa: E402


# ---------------------------------------------------------------------------
# The context: everything psi4's SAPT2 carries, in psi4's naming
# ---------------------------------------------------------------------------

class Sapt2:
    """Mirrors the SAPT2 object's state: MO blocks of the exact three-index
    factors with psi4's dressing columns, the one-electron MO matrices, the
    omega integrals and the amplitudes."""

    def __init__(self, c, foccA=0, foccB=0):
        self.c = c
        mol = c["mol"]
        self.nao = mol.nao

        # Full MO sets, occupied first -- psi4's CA_/CB_.
        self.CA = np.hstack([c["Cocc_A"], c["Cvir_A"]])
        self.CB = np.hstack([c["Cocc_B"], c["Cvir_B"]])
        self.evalsA = np.hstack([c["eps_occ_A"], c["eps_vir_A"]])
        self.evalsB = np.hstack([c["eps_occ_B"], c["eps_vir_B"]])
        self.noccA = c["Cocc_A"].shape[1]
        self.noccB = c["Cocc_B"].shape[1]
        self.nvirA = c["Cvir_A"].shape[1]
        self.nvirB = c["Cvir_B"].shape[1]
        self.foccA, self.foccB = foccA, foccB
        self.aoccA = self.noccA - foccA
        self.aoccB = self.noccB - foccB

        # One-electron matrices in MO pairs (sapt.cc:228-271).
        S, VA, VB = c["S"], c["V_A"], c["V_B"]
        self.sAB = self.CA.T @ S @ self.CB
        self.vBAA = self.CA.T @ VB @ self.CA
        self.vABB = self.CB.T @ VA @ self.CB
        self.vAAB = self.CA.T @ VA @ self.CB
        self.vBAB = self.CA.T @ VB @ self.CB

        # NA_/NB_ are ELECTRON counts (sapt.cc:169,179), and enuc_ is the
        # intermolecular nuclear repulsion (sapt.cc:184). The dressing spreads
        # the one- and zero-electron parts of the interaction over them.
        self.NA = 2 * self.noccA
        self.NB = 2 * self.noccB
        self.enuc = c["E_nuc"]
        self.esq = np.sqrt(self.enuc / (self.NA * self.NB))

        # Exact three-index factors: eigendecompose the AO ERI matrix and keep
        # the positive spectrum. Rank ~ nao(nao+1)/2; reconstruction error is
        # O(eps * lambda_max), far below every tolerance used here.
        eri = mol.intor("int2e").reshape(self.nao**2, self.nao**2)
        w, U = np.linalg.eigh(eri)
        keep = w > w.max() * 1e-15
        B = U[:, keep] * np.sqrt(w[keep])
        self.ndf = B.shape[1]
        self.B_ao = np.ascontiguousarray(B.reshape(self.nao, self.nao, self.ndf))
        self.fact_err = np.abs(eri - B @ B.T).max()
        del eri, U, B

        # Raw MO blocks, each with the three (zeroed) dressing columns.
        def mo(X, Y):
            t = np.tensordot(self.B_ao, X, axes=([0], [0]))   # (q, P, x)
            t = np.tensordot(t, Y, axes=([0], [0]))           # (P, x, y)
            out = np.zeros((X.shape[1], Y.shape[1], self.ndf + 3))
            out[:, :, :self.ndf] = np.transpose(t, (1, 2, 0))
            return out

        oA, vA = self.CA[:, :self.noccA], self.CA[:, self.noccA:]
        oB, vB = self.CB[:, :self.noccB], self.CB[:, self.noccB:]
        self.bAA = mo(oA, oA)
        self.bAR = mo(oA, vA)
        self.bRR = mo(vA, vA)
        self.bBB = mo(oB, oB)
        self.bBS = mo(oB, vB)
        self.bSS = mo(vB, vB)
        self.bAB = mo(oA, oB)
        self.bAS = mo(oA, vB)
        self.bRB = mo(vA, oB)

        # diagAA_/diagBB_, dressed (sapt2.cc:815-833 with get_diag_*_ints(1)).
        nf = self.ndf
        self.diagAA = np.zeros(nf + 3)
        self.diagAA[:nf] = np.einsum("aaP->P", self.bAA[:, :, :nf])
        self.diagAA[nf] = self.noccA
        self.diagAA[nf + 1] = np.trace(self.vBAA[:self.noccA, :self.noccA]) / self.NB
        self.diagAA[nf + 2] = self.noccA * self.esq
        self.diagBB = np.zeros(nf + 3)
        self.diagBB[:nf] = np.einsum("bbP->P", self.bBB[:, :, :nf])
        self.diagBB[nf] = np.trace(self.vABB[:self.noccB, :self.noccB]) / self.NA
        self.diagBB[nf + 1] = self.noccB
        self.diagBB[nf + 2] = self.noccB * self.esq

        # Omega integrals, w = v + 2J in MO blocks (sapt2.cc:815-905). The
        # GEMV there runs over ndf_ only -- no dressing enters w.
        nA, nB = self.noccA, self.noccB
        self.wBAA = self.vBAA[:nA, :nA] + 2.0 * self.bAA[:, :, :nf] @ self.diagBB[:nf]
        self.wBAR = self.vBAA[:nA, nA:] + 2.0 * self.bAR[:, :, :nf] @ self.diagBB[:nf]
        self.wBRR = self.vBAA[nA:, nA:] + 2.0 * self.bRR[:, :, :nf] @ self.diagBB[:nf]
        self.wABB = self.vABB[:nB, :nB] + 2.0 * self.bBB[:, :, :nf] @ self.diagAA[:nf]
        self.wABS = self.vABB[:nB, nB:] + 2.0 * self.bBS[:, :, :nf] @ self.diagAA[:nf]
        self.wASS = self.vABB[nB:, nB:] + 2.0 * self.bSS[:, :, :nf] @ self.diagAA[:nf]

    # -- dressed blocks, transcribing utils.cc get_*_ints ------------------
    # A-monomer pairs dress as (delta, vB/NB, enuc*delta); B-monomer pairs as
    # (vA/NA, delta, enuc*delta); mixed pairs carry sAB where a same-monomer
    # pair carries delta. dress=1/2 on mixed blocks selects which side of the
    # convention the block plays.

    def AA(self, dress=1, f1=0, f2=0):
        A = self.bAA[f1:, f2:, :].copy()
        if dress:
            nf = self.ndf
            A[:, :, nf + 1] = self.vBAA[f1:self.noccA, f2:self.noccA] / self.NB
            for i in range(A.shape[0]):
                j = i + f1 - f2
                if 0 <= j < A.shape[1]:
                    A[i, j, nf] = 1.0
                    A[i, j, nf + 2] = self.esq
        return A

    def BB(self, dress=1, f1=0, f2=0):
        A = self.bBB[f1:, f2:, :].copy()
        if dress:
            nf = self.ndf
            A[:, :, nf] = self.vABB[f1:self.noccB, f2:self.noccB] / self.NA
            for i in range(A.shape[0]):
                j = i + f1 - f2
                if 0 <= j < A.shape[1]:
                    A[i, j, nf + 1] = 1.0
                    A[i, j, nf + 2] = self.esq
        return A

    def AR(self, dress=1, focc=0):
        A = self.bAR[focc:, :, :].copy()
        if dress:
            A[:, :, self.ndf + 1] = self.vBAA[focc:self.noccA, self.noccA:] / self.NB
        return A

    def BS(self, dress=1, focc=0):
        A = self.bBS[focc:, :, :].copy()
        if dress:
            A[:, :, self.ndf] = self.vABB[focc:self.noccB, self.noccB:] / self.NA
        return A

    def RR(self, dress=1):
        A = self.bRR.copy()
        if dress:
            nf = self.ndf
            A[:, :, nf + 1] = self.vBAA[self.noccA:, self.noccA:] / self.NB
            for r in range(self.nvirA):
                A[r, r, nf] = 1.0
                A[r, r, nf + 2] = self.esq
        return A

    def SS(self, dress=1):
        A = self.bSS.copy()
        if dress:
            nf = self.ndf
            A[:, :, nf] = self.vABB[self.noccB:, self.noccB:] / self.NA
            for s in range(self.nvirB):
                A[s, s, nf + 1] = 1.0
                A[s, s, nf + 2] = self.esq
        return A

    def AB(self, dress, fA=0, fB=0):
        A = self.bAB[fA:, fB:, :].copy()
        nf = self.ndf
        s = self.sAB[fA:self.noccA, fB:self.noccB]
        if dress == 1:
            A[:, :, nf] = s
            A[:, :, nf + 1] = self.vBAB[fA:self.noccA, fB:self.noccB] / self.NB
            A[:, :, nf + 2] = self.esq * s
        elif dress == 2:
            A[:, :, nf] = self.vAAB[fA:self.noccA, fB:self.noccB] / self.NA
            A[:, :, nf + 1] = s
            A[:, :, nf + 2] = self.esq * s
        return A

    def AS(self, dress, focc=0):
        A = self.bAS[focc:, :, :].copy()
        nf = self.ndf
        s = self.sAB[focc:self.noccA, self.noccB:]
        if dress == 1:
            A[:, :, nf] = s
            A[:, :, nf + 1] = self.vBAB[focc:self.noccA, self.noccB:] / self.NB
            A[:, :, nf + 2] = self.esq * s
        elif dress == 2:
            A[:, :, nf] = self.vAAB[focc:self.noccA, self.noccB:] / self.NA
            A[:, :, nf + 1] = s
            A[:, :, nf + 2] = self.esq * s
        return A

    def RB(self, dress, focc=0):
        A = self.bRB[:, focc:, :].copy()
        nf = self.ndf
        s = self.sAB[self.noccA:, focc:self.noccB]
        if dress == 1:
            A[:, :, nf] = self.vAAB[self.noccA:, focc:self.noccB] / self.NA
            A[:, :, nf + 1] = s
            A[:, :, nf + 2] = self.esq * s
        elif dress == 2:
            A[:, :, nf] = s
            A[:, :, nf + 1] = self.vBAB[self.noccA:, focc:self.noccB] / self.NB
            A[:, :, nf + 2] = self.esq * s
        return A


def ein(*args):
    return np.einsum(*args, optimize=True)


# ---------------------------------------------------------------------------
# Amplitudes -- amplitudes.cc, one monomer at a time
# ---------------------------------------------------------------------------

def antisym(t):
    """theta_{ar,a'r'} = 2 t_{ar,a'r'} - t_{a'r,ar'} (utils.cc:1326-1346)."""
    return 2.0 * t - t.transpose(2, 1, 0, 3)


def monomer_amplitudes(x, side):
    """t, theta, Theta^P, pOO, pVV, Y2, the T2 singles and t2 for one monomer.

    Transcribes SAPT2::amplitudes() (amplitudes.cc:42-105) for the labels the
    four SAPT2 terms consume: tARAR, Theta AR, pAA/pRR, Y2 AR, T2 AR, t2ARAR
    and Theta 2 AR (and the BS mirrors).
    """
    if side == "A":
        focc, nocc, nvir = x.foccA, x.noccA, x.nvirA
        evals = x.evalsA
        bOO, bOV, bVV = x.bAA, x.bAR, x.bRR
        OVd = x.AR(1, focc)
    else:
        focc, nocc, nvir = x.foccB, x.noccB, x.nvirB
        evals = x.evalsB
        bOO, bOV, bVV = x.bBB, x.bBS, x.bSS
        OVd = x.BS(1, focc)

    nf = x.ndf
    eo, ev = evals[focc:nocc], evals[nocc:]
    D = (eo[:, None, None, None] + eo[None, None, :, None]
         - ev[None, :, None, None] - ev[None, None, None, :])

    out = {}
    # tOVOV (amplitudes.cc:106-137): t = (ar|a'r')/D, denominator INCLUDED,
    # chemists' notation, NOT antisymmetrized.
    v = ein("arP,cqP->arcq", bOV[focc:, :, :nf], bOV[focc:, :, :nf])
    t = v / D
    th = antisym(t)
    out["t"] = t
    out["theta"] = th

    # pOOpVV (amplitudes.cc:138-166): the unrelaxed MP2 density blocks,
    # WITHOUT the factor -2/+2 -- those live in the consuming terms.
    out["pOO"] = ein("arcq,brcq->ab", t, th)
    out["pVV"] = ein("axcr,axcq->rq", t, th)

    # theta (amplitudes.cc:167-200): Theta^P = theta . B_OV(dressed).
    out["ThOV"] = ein("arcq,cqP->arP", th, OVd)

    # Y2 (amplitudes.cc:201-317). Y2_3 first; the T2 singles amplitude is
    # written from Y2_3 ALONE, before Y2_1/Y2_2 are added (amplitudes.cc:216).
    y = (-ein("acP,crP->ar", bOO[focc:, focc:, :], out["ThOV"])
         + ein("aqP,rqP->ar", out["ThOV"], bVV))
    out["tOV1"] = y / (eo[:, None] - ev[None, :])
    # Y2_1 (amplitudes.cc:230-262)
    X = ein("xyP,xy->P", bVV, out["pVV"])
    y = y + 2.0 * bOV[focc:, :, :] @ X
    y -= ein("qx,aqP,rxP->ar", out["pVV"], bOV[focc:, :, :], bVV)
    # Y2_2 (amplitudes.cc:263-291)
    Xp = ein("xyP,xy->P", bOO[focc:, focc:, :], out["pOO"])
    y = y - 2.0 * bOV[focc:, :, :] @ Xp
    y += ein("cx,xaP,crP->ar", out["pOO"], bOO[focc:, focc:, :],
             bOV[focc:, :, :])
    out["y2"] = y

    # t2OVOV (amplitudes.cc:318-395): the second-order doubles.
    vOOVV = ein("acP,rqP->acrq", bOO[focc:, focc:, :nf], bVV[:, :, :nf])
    N = -ein("xrcy,axqy->arcq", t, vOOVV)          # after the first swap
    N -= ein("arxy,cxqy->arcq", t, vOOVV)          # the unswapped partner
    N += ein("arP,cqP->arcq", bOV[focc:, :, :], out["ThOV"])
    vOOOO = ein("axP,cyP->axcy", bOO[focc:, focc:, :nf], bOO[focc:, focc:, :nf])
    N += 0.5 * ein("axcy,xryq->arcq", vOOOO, t)
    vVVVV = ein("rxP,qyP->rxqy", bVV[:, :, :nf], bVV[:, :, :nf])
    N += 0.5 * ein("rxqy,axcy->arcq", vVVVV, t)
    del vVVVV
    N = N + N.transpose(2, 3, 0, 1)                # symmetrize (:599)
    out["t2"] = N / D

    # Theta 2 (amplitudes.cc:100-105): same theta() on t2.
    out["Th2OV"] = ein("arcq,cqP->arP", antisym(out["t2"]), OVd)
    return out


def mp2_energy(x, amps, side):
    """The monomer MP2 correlation energy the amplitudes imply.

    E2 = sum t_{ar,a'r'} [2 (ar|a'r') - (a'r|ar')] -- must match PySCF's MP2
    on the same ghosted monomer to machine precision before anything
    downstream is trusted (plan section 0.2)."""
    focc = x.foccA if side == "A" else x.foccB
    bOV = (x.bAR if side == "A" else x.bBS)[focc:, :, :x.ndf]
    v = ein("arP,cqP->arcq", bOV, bOV)
    return float(ein("arcq,arcq->", amps["t"], 2.0 * v - v.transpose(2, 1, 0, 3)))


# ---------------------------------------------------------------------------
# The exchange-induction potentials uAR/vAR -- exch-ind20.cc:682-1145
# ---------------------------------------------------------------------------

def k2f_integrals_A(x):
    """uAR ("AR Exch12 K2f Integrals") and vAR ("AR Exch-Ind Integrals"),
    transcribing exch_ind20rA_B (exch-ind20.cc:682-911). The thirteen pieces
    are shared; only a few coefficients differ between the two outputs."""
    nA, nB = x.noccA, x.noccB
    sOcc = x.sAB[:nA, :nB]
    sVir = x.sAB[nA:, :nB]
    AB1, AB2 = x.AB(1), x.AB(2)
    RB1 = x.RB(1)
    AA1, AR1, BB1 = x.AA(1), x.AR(1), x.BB(1)

    t1 = ein("abP,rbP->ar", AB1, RB1)
    t2 = ein("ab,rbP,P->ar", sOcc, RB1, x.diagAA)
    t3 = ein("acP,cb,rbP->ar", AA1, sOcc, RB1)
    Cp = ein("ab,abP->P", sOcc, AB2)
    t4 = AR1 @ Cp
    t5 = ein("cb,abP,arP->cr", sOcc, AB2, AR1)
    t6 = ein("abP,P,rb->ar", AB1, x.diagBB, sVir)
    t7 = ein("adP,bdP,rb->ar", AB1, BB1, sVir)
    X_BB = ein("bdP,P->bd", BB1, x.diagAA)
    t8 = ein("ad,bd,rb->ar", sOcc, X_BB, sVir)
    X_AA = ein("acP,P->ac", AA1, x.diagBB)
    S_AR = sOcc @ sVir.T
    t9 = X_AA @ S_AR
    t10 = ein("cb,acP,rd,dbP->ar", sOcc, AA1, sVir, BB1)
    S_BB = sOcc.T @ sOcc
    Cp2 = ein("bdP,bd->P", BB1, S_BB)
    t11 = AR1 @ Cp2
    S_AA = sOcc @ sOcc.T
    t12 = S_AA @ ein("arP,P->ar", AR1, x.diagBB)
    t13 = ein("ab,cd,bdP,crP->ar", sOcc, x.sAB[:nA, :nB], BB1, AR1)
    # t13: C_p_BA[b,c,P] = sum_d sAB[c,d] BB1[b,d,P]; C_p_AA[a,c,P] =
    # sum_b sAB[a,b] C_p_BA[b,c,P]; then sum_c C_p_AA . AR1[c,r,P].

    uAR = (t1 + 2.0 * t2 - t3 + 4.0 * t4 - t5 + 2.0 * t6 - t7 - 4.0 * t8
           - 4.0 * t9 + 2.0 * t10 - 4.0 * t11 - 4.0 * t12 + 2.0 * t13)
    vAR = (t1 + 2.0 * t2 - t3 + 2.0 * t4 - t5 + 2.0 * t6 - t7 - 2.0 * t8
           - 2.0 * t9 + t10 - 2.0 * t11 - 2.0 * t12 + t13)
    return uAR, vAR


def k2f_integrals_B(x):
    """uBS and vBS, transcribing exch_ind20rB_A (exch-ind20.cc:913-1145)."""
    nA, nB = x.noccA, x.noccB
    sOcc = x.sAB[:nA, :nB]
    sVir = x.sAB[:nA, nB:]           # (occA, virB)
    AB1, AB2 = x.AB(1), x.AB(2)
    AS1 = x.AS(1)
    AA1, BB1, BS1 = x.AA(1), x.BB(1), x.BS(1)

    s1 = ein("abP,asP->bs", AB2, AS1)
    s2 = ein("ab,asP,P->bs", sOcc, AS1, x.diagBB)
    s3 = ein("bcP,ac,asP->bs", BB1, sOcc, AS1)
    Cp = ein("ab,abP->P", sOcc, AB1)
    s4 = ein("bsP,P->bs", BS1, Cp)
    s5 = ein("ac,abP,bsP->cs", sOcc, AB1, BS1)
    s6 = ein("abP,P,as->bs", AB2, x.diagAA, sVir)
    # s7: X_AB[a',b] accumulates AA1[a,a',P] AB2[a,b,P] over a, then dots with
    # sAB[a', noccB+s]. Note AB2 -- exch_ind20rB_A reassigns B_p_AB to
    # get_AB_ints(2) just above this term (exch-ind20.cc:995-1023), where the
    # A-side's mirror term uses dress 1.
    s7 = ein("acP,abP,cs->bs", AA1, AB2, sVir)
    X_AA = ein("acP,P->ac", AA1, x.diagBB)
    s8 = ein("ac,cb,as->bs", X_AA, sOcc, sVir)
    X_BB = ein("bcP,P->bc", BB1, x.diagAA)
    S_BS = sOcc.T @ sVir
    s9 = X_BB @ S_BS
    s10 = ein("cd,bdP,es,ecP->bs", sOcc, BB1, sVir, AA1)
    # s10: C_p_BA[b,c,P] = sum_d sAB[c,d] BB1[b,d,P]; C_p_SA[s,c,P] =
    # sum_e sAB[e,noccB+s] AA1[e,c,P]; contract over (c,P).
    S_AA = sOcc @ sOcc.T
    Cp3 = ein("acP,ac->P", AA1, S_AA)
    s11 = ein("bsP,P->bs", BS1, Cp3)
    S_BB = sOcc.T @ sOcc
    s12 = S_BB @ ein("bsP,P->bs", BS1, x.diagAA)
    # s13: C_p_AB[e,c,P] = sum_a sAB[a,c] AA1[e,a,P]; C_p_BB[b,c,P] =
    # sum_e sAB[e,b] C_p_AB[e,c,P]; then sum_b C_p_BB[b,c,P] BS1[b,s,P].
    s13 = ein("eb,ac,eaP,bsP->cs", sOcc, sOcc, AA1, BS1)

    uBS = (s1 + 2.0 * s2 - s3 + 4.0 * s4 - s5 + 2.0 * s6 - s7 - 4.0 * s8
           - 4.0 * s9 + 2.0 * s10 - 4.0 * s11 - 4.0 * s12 + 2.0 * s13)
    vBS = (s1 + 2.0 * s2 - s3 + 2.0 * s4 - s5 + 2.0 * s6 - s7 - 2.0 * s8
           - 2.0 * s9 + s10 - 2.0 * s11 - 2.0 * s12 + s13)
    return uBS, vBS


# ---------------------------------------------------------------------------
# Elst12 -- elst12.cc
# ---------------------------------------------------------------------------

def elst12(x, ampsA, ampsB, chfA, chfB):
    """Elst12,r = Elst120 + Elst102 (elst12.cc:38-58; elst120 at :61-95).

    Each half is the partner's electrostatic potential contracted against the
    monomer's MP2 density: the oo and vv blocks directly, and the ov (orbital
    relaxation) block through the CPHF coefficients already solved for Ind20r
    -- Tr(Z w) = Tr(Y x_w), so no Z-vector is ever formed."""
    fA, fB = x.foccA, x.foccB
    e120 = (-2.0 * ein("ab,ab->", ampsA["pOO"], x.wBAA[fA:, fA:])
            + 2.0 * ein("rq,rq->", ampsA["pVV"], x.wBRR)
            + 4.0 * ein("ar,ar->", ampsA["y2"], chfA[fA:, :]))
    e102 = (-2.0 * ein("ab,ab->", ampsB["pOO"], x.wABB[fB:, fB:])
            + 2.0 * ein("rq,rq->", ampsB["pVV"], x.wASS)
            + 4.0 * ein("ar,ar->", ampsB["y2"], chfB[fB:, :]))
    return e120 + e102, e120, e102


# ---------------------------------------------------------------------------
# Exch11 -- exch11.cc
# ---------------------------------------------------------------------------

def exch110(x, ThAR):
    """exch110 (exch11.cc:57-127): first-order exchange against monomer A's
    correlated (theta-level) density. Also reused by Exch12 with Theta 2
    (exch12.cc:44)."""
    fA = x.foccA
    nA, nB = x.noccA, x.noccB
    saB = x.sAB[fA:nA, :nB]
    sRB = x.sAB[nA:, :nB]

    AB2 = x.AB(2, fA, 0)
    C_p_AB = ein("rb,arP->abP", sRB, ThAR)
    e1 = -2.0 * ein("abP,abP->", C_p_AB, AB2)

    BB1 = x.BB(1)
    C_p_BB = ein("ab,acP->bcP", saB, C_p_AB)
    e2 = 4.0 * ein("bcP,bcP->", BB1, C_p_BB)

    RB1 = x.RB(1)
    C_p_AR = ein("ab,rbP->arP", saB, RB1)
    e3 = -2.0 * ein("arP,arP->", ThAR, C_p_AR)

    xAR = saB @ sRB.T
    yAR = ThAR @ x.diagBB
    e4 = -8.0 * ein("ar,ar->", xAR, yAR)
    return e1 + e2 + e3 + e4


def exch101(x, ThBS):
    """exch101 (exch11.cc:129-194), the B-side mirror."""
    fB = x.foccB
    nA, nB = x.noccA, x.noccB
    sAb = x.sAB[:nA, fB:nB]
    sAS = x.sAB[:nA, nB:]

    AB1 = x.AB(1, 0, fB)
    C_p_AB = ein("as,bsP->abP", sAS, ThBS)
    e1 = -2.0 * ein("abP,abP->", C_p_AB, AB1)

    AA1 = x.AA(1)
    C_p_AA = ein("cb,abP->acP", sAb, C_p_AB)
    e2 = 4.0 * ein("acP,acP->", AA1, C_p_AA)

    AS1 = x.AS(1)
    C_p_BS = ein("ab,asP->bsP", sAb, AS1)
    e3 = -2.0 * ein("bsP,bsP->", ThBS, C_p_BS)

    xBS = sAb.T @ sAS
    yBS = ThBS @ x.diagAA
    e4 = -8.0 * ein("bs,bs->", xBS, yBS)
    return e1 + e2 + e3 + e4


def exch11(x, ampsA, ampsB):
    """Exch11 = Exch110 + Exch101 (exch11.cc:37-55)."""
    return exch110(x, ampsA["ThOV"]) + exch101(x, ampsB["ThOV"])


# ---------------------------------------------------------------------------
# Ind22 -- ind22.cc
# ---------------------------------------------------------------------------

def ind220_side(x, amps, side):
    """One direction of Ind22 (ind220, ind22.cc:58-102; ind202 mirrors it with
    every A-label swapped for the B one, ind22.cc:104-148).

    The response here is UNCOUPLED -- iAR = wBAR / (e_a - e_r) -- and the
    seven pieces are MBPT contractions of that response with the monomer's
    amplitudes; no CPHF solve appears anywhere in psi4's Ind22.

    A correction to the plan's section 0.1: pieces 2-7 vanish identically when
    the stored amplitudes are zeroed, but piece 1 does NOT -- its perturbed
    doubles x1 are built on the fly from bare two-electron integrals times the
    uncoupled response, with the stored t entering only two of its four RHS
    terms. Zeroing t therefore turns SAPT2 into SAPT0 plus exactly that
    residual, and the t->0 test must assert it piecewise."""
    nf = x.ndf
    if side == "A":
        focc, nocc, nvir = x.foccA, x.noccA, x.nvirA
        evals = x.evalsA
        wOO, wOV, wVV = x.wBAA, x.wBAR, x.wBRR
        bOO, bOV, bVV = x.bAA, x.bAR, x.bRR
        wOV_other = x.wABS
        evals_other, focc_other, nocc_other = x.evalsB, x.foccB, x.noccB
        bOV_other = x.bBS
    else:
        focc, nocc, nvir = x.foccB, x.noccB, x.nvirB
        evals = x.evalsB
        wOO, wOV, wVV = x.wABB, x.wABS, x.wASS
        bOO, bOV, bVV = x.bBB, x.bBS, x.bSS
        wOV_other = x.wBAR
        evals_other, focc_other, nocc_other = x.evalsA, x.foccA, x.noccA
        bOV_other = x.bAR

    eo, ev = evals[focc:nocc], evals[nocc:]
    i_OV = wOV[focc:, :] / (eo[:, None] - ev[None, :])
    eoo, evo = evals_other[focc_other:nocc_other], evals_other[nocc_other:]
    i_other = wOV_other[focc_other:, :] / (eoo[:, None] - evo[None, :])

    t, th = amps["t"], amps["theta"]
    tOV1, pOO, pVV = amps["tOV1"], amps["pOO"], amps["pVV"]
    D4 = (eo[:, None, None, None] + eo[None, None, :, None]
          - ev[None, :, None, None] - ev[None, None, None, :])

    # ind220_1 (ind22.cc:150-229): the wB-perturbed first-order doubles.
    # x1[a,r,c,q] = sum_x i[a,x] (x r| c q) - sum_x i[x,r] (a x| c q)
    #             - sum_x wOO[a,x] t[x,r,c,q] + sum_x t[a,r,c,x] wVV[q,x]
    x1 = ein("ax,xrP,cqP->arcq", i_OV, bVV, bOV[focc:, :, :])
    x1 -= ein("xr,axP,cqP->arcq", i_OV, bOO[focc:, focc:, :], bOV[focc:, :, :])
    x1 -= ein("ax,xrcq->arcq", wOO[focc:, focc:], t)
    x1 += ein("arcx,qx->arcq", t, wVV)
    x1 = x1 + x1.transpose(2, 3, 0, 1)
    y1 = antisym(x1)
    e1 = float(ein("arcq,arcq->", x1 / D4, y1))

    # ind220_2 (:231-255)
    z2 = i_OV @ wVV.T - wOO[focc:, focc:] @ i_OV
    e2 = 4.0 * float(ein("ar,ar->", tOV1, z2))

    # ind220_3 (:257-296)
    xOO = i_OV @ wOV[focc:, :].T
    xVV = i_OV.T @ wOV[focc:, :]
    e3 = -2.0 * float(ein("ab,ab->", pOO, xOO)) \
         - 2.0 * float(ein("rq,rq->", pVV, xVV))

    # ind220_4 (:298-334)
    xOO4 = i_OV @ i_OV.T
    xVV4 = i_OV.T @ i_OV
    C = ein("ac,crP->arP", xOO4, bOV[focc:, :, :]) \
        + ein("rq,aqP->arP", xVV4, bOV[focc:, :, :])
    e4 = -2.0 * float(ein("arP,arP->", C, amps["ThOV"]))

    # ind220_5 (:336-366): antisym(t2) with the denominator multiplied back.
    # NOTE the amplitude is the SECOND-order t2 (ind22.cc:88 passes "t2ARAR
    # Amplitudes"), not the first-order t -- the one place Ind22 consumes t2.
    e5 = 2.0 * float(ein("arcq,arcq,ar,cq->", antisym(amps["t2"]), D4,
                         i_OV, i_OV))

    # ind220_6 (:368-421): g = 2(ar|cq) - (ac|rq).
    g = 2.0 * ein("arP,cqP->arcq", bOV[focc:, :, :nf], bOV[focc:, :, :nf]) \
        - ein("acP,rqP->arcq", bOO[focc:, focc:, :nf], bVV[:, :, :nf])
    xar = ein("arcq,cq->ar", g, i_OV)
    yar = ein("arcq,cq->ar", th, i_OV)
    e6 = -4.0 * float(ein("ar,ar->", xar, yar))

    # ind220_7 (:423-460): this monomer's density against the OTHER monomer's
    # uncoupled response, raw blocks throughout -- pure Coulomb coupling.
    W = ein("acP,ac->P", bOO[focc:, focc:, :], pOO)
    X = ein("rqP,rq->P", bVV, pVV)
    Y = ein("arP,ar->P", bOV[focc:, :, :], tOV1)
    Z = ein("bsP,bs->P", bOV_other[focc_other:, :, :], i_other)
    e7 = -8.0 * float(W @ Z) + 8.0 * float(X @ Z) + 16.0 * float(Y @ Z)

    return [e1, e2, e3, e4, e5, e6, e7]


def ind22(x, ampsA, ampsB):
    """Ind22 = Ind220 + Ind202 (ind22.cc:37-56)."""
    pA = ind220_side(x, ampsA, "A")
    pB = ind220_side(x, ampsB, "B")
    return sum(pA) + sum(pB), pA, pB


# ---------------------------------------------------------------------------
# Exch12 -- exch12.cc
# ---------------------------------------------------------------------------

def exch111(x, ThAR, ThBS):
    """exch111 (exch12.cc:99-155): both monomers' theta densities exchanging."""
    fA, fB = x.foccA, x.foccB
    nA, nB = x.noccA, x.noccB
    C_p_AB = ein("rb,arP->abP", x.sAB[nA:, fB:nB], ThAR)
    D_p_AB = ein("as,bsP->abP", x.sAB[fA:nA, nB:], ThBS)
    e1 = -4.0 * ein("abP,abP->", C_p_AB, D_p_AB)

    C_p_AS = ein("rs,arP->asP", x.sAB[nA:, nB:], ThAR)
    C_p_BS = ein("ab,asP->bsP", x.sAB[fA:nA, fB:nB], C_p_AS)
    e2 = -4.0 * ein("bsP,bsP->", ThBS, C_p_BS)
    return e1 + e2


def exch120_k2f(x, ampsA, uAR):
    """exch120_k2f (exch12.cc:157-241): the T2 singles against the K2f
    potential and its overlap-dressed relatives."""
    fA = x.foccA
    nA, nB = x.noccA, x.noccB
    tAR = ampsA["tOV1"]
    saB = x.sAB[fA:nA, :nB]
    sRB = x.sAB[nA:, :nB]

    e1 = -2.0 * ein("ar,ar->", tAR, uAR[fA:, :])

    RB2 = x.RB(2)
    AB2 = x.AB(2)
    C_p_AB = ein("ar,rbP->abP", tAR, RB2)
    e2 = -2.0 * ein("abP,abP->", AB2[fA:, :, :], C_p_AB)

    BB1 = x.BB(1)
    C_p_BB = ein("ab,acP->bcP", saB, C_p_AB)
    e3 = 2.0 * ein("bcP,bcP->", BB1, C_p_BB)

    xAB = C_p_AB @ x.diagBB
    e4 = -4.0 * ein("ab,ab->", xAB, saB)

    xAB5 = AB2[fA:, :, :] @ x.diagAA
    yAB = tAR @ sRB
    e5 = -4.0 * ein("ab,ab->", xAB5, yAB)

    AA1 = x.AA(1)
    # D_p_AB[c,b,P] = sum_{x in aocc} yAB[x,b] AA1[c, foccA+x, P]
    D_p_AB = ein("xb,cxP->cbP", yAB, AA1[:, fA:, :])
    e6 = 2.0 * ein("cbP,cbP->", AB2, D_p_AB)

    AR1 = x.AR(1)
    C_p_AA = ein("cr,arP->caP", tAR, AR1)         # c in aocc, a in nocc
    D_p_AA = ein("ab,cbP->caP", x.sAB[:nA, :nB], AB2[fA:, :, :])
    e7 = 2.0 * ein("caP,caP->", C_p_AA, D_p_AA)
    return e1 + e2 + e3 + e4 + e5 + e6 + e7


def exch102_k2f(x, ampsB, uBS):
    """exch102_k2f (exch12.cc:243-380), the B-side mirror."""
    fB = x.foccB
    nA, nB = x.noccA, x.noccB
    tBS = ampsB["tOV1"]
    sAb = x.sAB[:nA, fB:nB]
    sAS = x.sAB[:nA, nB:]

    e1 = -2.0 * ein("bs,bs->", tBS, uBS[fB:, :])

    AS2 = x.AS(2)
    AB1 = x.AB(1)
    C_p_AB = ein("bs,asP->abP", tBS, AS2)
    e2 = -2.0 * ein("abP,abP->", AB1[:, fB:, :], C_p_AB)

    AA1 = x.AA(1)
    C_p_AA = ein("cb,abP->acP", sAb, C_p_AB)
    e3 = 2.0 * ein("acP,acP->", AA1, C_p_AA)

    xAB = C_p_AB @ x.diagAA
    e4 = -4.0 * ein("ab,ab->", xAB, sAb)

    xAB5 = AB1[:, fB:, :] @ x.diagBB
    yAB = sAS @ tBS.T
    e5 = -4.0 * ein("ab,ab->", xAB5, yAB)

    BB1 = x.BB(1)
    D_p_AB = ein("ax,xbP->abP", yAB, BB1[fB:, :, :])
    e6 = 2.0 * ein("abP,abP->", AB1, D_p_AB)

    BS1 = x.BS(1)
    C_p_BB = ein("cs,bsP->cbP", tBS, BS1)         # c in aocc, b in nocc
    # D_p_BB[c,b,P] = sum_a sAB[a,b] AB1[a, foccB+c, P]
    D_p_BB = ein("ab,acP->cbP", x.sAB[:nA, :nB], AB1[:, fB:, :])
    e7 = 2.0 * ein("cbP,cbP->", C_p_BB, D_p_BB)
    return e1 + e2 + e3 + e4 + e5 + e6 + e7


def exch120_k11u_1(x, ampsA):
    """exch120_k11u_1 (exch12.cc:382-718 upper half): the pRR density."""
    nA, nB = x.noccA, x.noccB
    pRR = ampsA["pVV"]
    sRB = x.sAB[nA:, :nB]
    sOcc = x.sAB[:nA, :nB]
    RB1, RB2 = x.RB(1), x.RB(2)
    AR1, AB2, BB1, RR1 = x.AR(1), x.AB(2), x.BB(1), x.RR(1)

    E = 0.0
    yRR = ein("rbP,qbP->rq", RB1, RB2)
    E += 2.0 * ein("rq,rq->", pRR, yRR)

    D_p_RB = ein("rq,qbP->rbP", pRR, RB1)
    E_p_RB = ein("rq,qbP->rbP", pRR, RB2)
    D_p_BR = ein("ab,arP->brP", sOcc, AR1)
    E -= 2.0 * ein("brP,rbP->", D_p_BR, D_p_RB)

    E += 4.0 * ein("rb,rb->", sRB, D_p_RB @ x.diagAA)
    E += 4.0 * ein("rb,rb->", sRB, E_p_RB @ x.diagBB)

    C_p_BB = ein("rb,rcP->bcP", sRB, E_p_RB)
    E -= 2.0 * ein("bcP,bcP->", BB1, C_p_BB)

    yRB = ein("arP,abP->rb", AR1, AB2)
    zRB = pRR @ sRB
    E -= 2.0 * ein("rb,rb->", yRB, zRB)

    xBR = D_p_BR @ x.diagBB
    E -= 8.0 * ein("br,rb->", xBR, zRB)

    D_p_BB = ein("rc,brP->bcP", zRB, D_p_BR)
    E += 4.0 * ein("bcP,bcP->", BB1, D_p_BB)

    zBB = sRB.T @ zRB
    wBB = BB1 @ x.diagAA
    E -= 4.0 * ein("bc,bc->", wBB, zBB)

    X = ein("rqP,rq->P", RR1, pRR)
    xAB = AB2 @ X
    E += 4.0 * ein("ab,ab->", sOcc, xAB)

    xBB = BB1 @ X
    yBB = sOcc.T @ sOcc
    E -= 4.0 * ein("bc,bc->", xBB, yBB)
    return -E


def exch102_k11u_1(x, ampsB):
    """exch102_k11u_1 (exch12.cc:562-718), the pSS mirror."""
    nA, nB = x.noccA, x.noccB
    pSS = ampsB["pVV"]
    sAS = x.sAB[:nA, nB:]
    sOcc = x.sAB[:nA, :nB]
    AS1, AS2 = x.AS(1), x.AS(2)
    BS1, AB1, AA1, SS1 = x.BS(1), x.AB(1), x.AA(1), x.SS(1)

    E = 0.0
    ySS = ein("asP,atP->st", AS1, AS2)
    E += 2.0 * ein("st,st->", pSS, ySS)

    D_p_AS = ein("st,atP->asP", pSS, AS1)
    E_p_AS = ein("st,atP->asP", pSS, AS2)
    F_p_AS = ein("ab,bsP->asP", sOcc, BS1)
    E -= 2.0 * ein("asP,asP->", D_p_AS, F_p_AS)

    E += 4.0 * ein("as,as->", sAS, D_p_AS @ x.diagBB)
    E += 4.0 * ein("as,as->", sAS, E_p_AS @ x.diagAA)

    C_p_AA = ein("cs,asP->acP", sAS, E_p_AS)
    E -= 2.0 * ein("acP,acP->", AA1, C_p_AA)

    yAS = ein("abP,bsP->as", AB1, BS1)
    zAS = sAS @ pSS
    E -= 2.0 * ein("as,as->", yAS, zAS)

    xAS = F_p_AS @ x.diagAA
    E -= 8.0 * ein("as,as->", xAS, zAS)

    D_p_AA = ein("cs,asP->acP", zAS, F_p_AS)
    E += 4.0 * ein("acP,acP->", AA1, D_p_AA)

    zAA = zAS @ sAS.T
    wAA = AA1 @ x.diagBB
    E -= 4.0 * ein("ac,ac->", wAA, zAA)

    X = ein("stP,st->P", SS1, pSS)
    xAB = AB1 @ X
    E += 4.0 * ein("ab,ab->", sOcc, xAB)

    xAA = AA1 @ X
    yAA = sOcc @ sOcc.T
    E -= 4.0 * ein("ac,ac->", xAA, yAA)
    return -E


def exch120_k11u_2(x, ampsA):
    """exch120_k11u_2 (exch12.cc:720-876): the pAA density."""
    fA = x.foccA
    nA, nB = x.noccA, x.noccB
    paa = ampsA["pOO"]
    sOcc = x.sAB[:nA, :nB]
    saB = x.sAB[fA:nA, :nB]
    A_aB = x.AB(1, fA, 0)
    B_AB, B_aB = x.AB(2), x.AB(2, fA, 0)
    A_aA, A_Aa, A_aa = x.AA(1, fA, 0), x.AA(1, 0, fA), x.AA(1, fA, fA)
    B_BB = x.BB(1)

    E = 0.0
    x1 = ein("abP,cbP->ac", A_aB, B_aB)
    E += 2.0 * ein("ac,ac->", x1, paa)

    X = ein("abP,ab->P", B_AB, sOcc)
    E += 4.0 * ein("ac,ac->", A_aa @ X, paa)

    C_p_aB = ein("db,adP->abP", sOcc, A_aA)
    x3 = ein("abP,cbP->ac", C_p_aB, B_aB)
    E -= 2.0 * ein("ac,ac->", x3, paa)

    xaB = B_aB @ x.diagAA
    x4 = saB @ xaB.T
    E += 4.0 * ein("ac,ac->", x4, paa)

    xaB2 = ein("dcP,dbP->cb", A_Aa, B_AB)
    x5 = saB @ xaB2.T
    E -= 2.0 * ein("ac,ac->", x5, paa)

    xaB3 = A_aB @ x.diagBB
    x6 = xaB3 @ saB.T
    E += 4.0 * ein("ac,ac->", x6, paa)

    xaB4 = ein("adP,bdP->ab", A_aB, B_BB)
    x7 = xaB4 @ saB.T
    E -= 2.0 * ein("ac,ac->", x7, paa)

    xBB = sOcc.T @ sOcc
    X2 = ein("bdP,bd->P", B_BB, xBB)
    E -= 4.0 * ein("ac,ac->", A_aa @ X2, paa)

    xBB2 = B_BB @ x.diagAA
    xaB5 = saB @ xBB2
    x9 = xaB5 @ saB.T
    E -= 4.0 * ein("ac,ac->", x9, paa)

    xaA = saB @ sOcc.T
    yaA = A_aA @ x.diagBB
    x10 = xaA @ yaA.T
    E -= 8.0 * ein("ac,ac->", x10, paa)

    C_p_aB2 = ein("ad,dbP->abP", saB, B_BB)
    C_p_aA = ein("db,abP->adP", sOcc, C_p_aB2)
    x11 = ein("adP,cdP->ac", C_p_aA, A_aA)
    E += 4.0 * ein("ac,ac->", x11, paa)
    return E


def exch102_k11u_2(x, ampsB):
    """exch102_k11u_2 (exch12.cc:878-1040): the pBB mirror."""
    fB = x.foccB
    nA, nB = x.noccA, x.noccB
    pbb = ampsB["pOO"]
    sOcc = x.sAB[:nA, :nB]
    sAb = x.sAB[:nA, fB:nB]
    A_AB, A_Ab = x.AB(1), x.AB(1, 0, fB)
    B_Ab = x.AB(2, 0, fB)
    A_AA = x.AA(1)
    B_bB, B_Bb, B_bb = x.BB(1, fB, 0), x.BB(1, 0, fB), x.BB(1, fB, fB)

    E = 0.0
    x1 = ein("abP,acP->bc", B_Ab, A_Ab)
    E += 2.0 * ein("bc,bc->", x1, pbb)

    X = ein("abP,ab->P", A_AB, sOcc)
    E += 4.0 * ein("bc,bc->", B_bb @ X, pbb)

    # x3[b,c] = sum_{a,P} C_p_Ab[a,b,P] A_Ab[a,c,P], with
    # C_p_Ab[a,b,P] = sum_d sOcc[a,d] B_bB[b,d,P]
    C_p_Ab = ein("ad,bdP->abP", sOcc, B_bB)
    x3 = ein("abP,acP->bc", C_p_Ab, A_Ab)
    E -= 2.0 * ein("bc,bc->", x3, pbb)

    xAb = A_Ab @ x.diagBB
    x4 = sAb.T @ xAb
    E += 4.0 * ein("bc,bc->", x4, pbb)

    xAb2 = ein("adP,dbP->ab", A_AB, B_Bb)
    x5 = sAb.T @ xAb2
    E -= 2.0 * ein("bc,bc->", x5, pbb)

    xAb3 = B_Ab @ x.diagAA
    x6 = xAb3.T @ sAb
    E += 4.0 * ein("bc,bc->", x6, pbb)

    xAb4 = ein("daP,dbP->ab", A_AA, B_Ab)
    x7 = xAb4.T @ sAb
    E -= 2.0 * ein("bc,bc->", x7, pbb)

    xAA = sOcc @ sOcc.T
    X2 = ein("adP,ad->P", A_AA, xAA)
    E -= 4.0 * ein("bc,bc->", B_bb @ X2, pbb)

    xAA2 = A_AA @ x.diagBB
    xAb5 = xAA2 @ sAb
    x9 = xAb5.T @ sAb
    E -= 4.0 * ein("bc,bc->", x9, pbb)

    xbB = sAb.T @ sOcc
    ybB = B_bB @ x.diagAA
    x10 = xbB @ ybB.T
    E -= 8.0 * ein("bc,bc->", x10, pbb)

    C_p_Ab2 = ein("db,daP->baP", sAb, A_AA)
    C_p_bB = ein("ac,baP->bcP", sOcc, C_p_Ab2)
    x11 = ein("bcP,dcP->bd", C_p_bB, B_bB)
    E += 4.0 * ein("bd,bd->", x11, pbb)
    return E


def exch120_k11u_3(x, ampsA):
    """exch120_k11u_3 (exch12.cc:1042-1170): the theta.t "square" terms."""
    nA, nB = x.noccA, x.noccB
    t, th = ampsA["t"], ampsA["theta"]
    sRB = x.sAB[nA:, :nB]
    RB1, RR1, BB1 = x.RB(1), x.RR(1), x.BB(1)

    # thetaRBAA[r1,b] = sum_r2 S_{r2 b} theta[a1,r1,a2,r2] (:1076-1080), and
    # the analogous tRBAA. The triangular r1>=r2 loop plus the swapped-t term
    # equals the full square because the dressed RR block is symmetric.
    thS = ein("crgq,qb->crgb", th, sRB)
    tS = ein("crgq,qb->crgb", t, sRB)

    E = 0.0
    # first square (:1090-1110): (r1 r2 | r3 b)-type
    Z = ein("rxP,qbP->rxqb", RR1, RB1)
    G = ein("cxgq,crgb->xqrb", t, thS)
    E += 2.0 * ein("rxqb,xqrb->", Z, G)

    # xRR . diagBB (:1124-1134)
    xRR = ein("crgb,cqgb->rq", tS, thS)
    yRR = RR1 @ x.diagBB
    E += 4.0 * ein("rq,rq->", xRR, yRR)

    # second square (:1138-1160): (r1 r2 | b b')-type
    Z2 = ein("rxP,bdP->rxbd", RR1, BB1)
    G2 = ein("cxgb,crgd->xbrd", tS, thS)
    E -= 2.0 * ein("rxbd,xbrd->", Z2, G2)
    return -E


def exch102_k11u_3(x, ampsB):
    """exch102_k11u_3 (exch12.cc:1172-1330), the B-side mirror."""
    nA, nB = x.noccA, x.noccB
    t, th = ampsB["t"], ampsB["theta"]
    sAS = x.sAB[:nA, nB:]
    AS1, SS1, AA1 = x.AS(1), x.SS(1), x.AA(1)

    thS = ein("csgt,at->csga", th, sAS)
    tS = ein("csgt,at->csga", t, sAS)

    E = 0.0
    Z = ein("sxP,atP->sxta", SS1, AS1)     # B_p_SA[(s,a)] = AS1[a,s]
    G = ein("cxgt,csga->xtsa", t, thS)
    E += 2.0 * ein("sxta,xtsa->", Z, G)

    xSS = ein("csga,ctga->st", tS, thS)
    ySS = SS1 @ x.diagAA
    E += 4.0 * ein("st,st->", xSS, ySS)

    Z2 = ein("sxP,adP->sxad", SS1, AA1)
    G2 = ein("cxga,csgd->xasd", tS, thS)
    E -= 2.0 * ein("sxad,xasd->", Z2, G2)
    return -E


def exch120_k11u_4(x, ampsA):
    """exch120_k11u_4 (exch12.cc:1332-1424): the occupied-space square."""
    fA = x.foccA
    nA, nB = x.noccA, x.noccB
    t, th = ampsA["t"], ampsA["theta"]
    saB = x.sAB[fA:nA, :nB]

    # M[a,a1,a2,a3] = sum_{r,r1} theta[a,r,a1,r1] t[a2,r,a3,r1] (:1352-1360)
    M = ein("argq,crhq->agch", th, t)
    AA1f = x.AA(1, fA, fA)
    C_p_AA = ein("agch,ghP->acP", M, AA1f)

    AB2f = x.AB(2, fA, 0)
    D_p_AA = ein("cb,abP->acP", saB, AB2f)
    E = 2.0 * ein("acP,acP->", C_p_AA, D_p_AA)

    sAA = saB @ saB.T
    X = ein("acP,ac->P", C_p_AA, sAA)
    E += 4.0 * float(X @ x.diagBB)

    BB1 = x.BB(1)
    C_p_AB = ein("ad,dbP->abP", saB, BB1)
    E_p_AA = ein("cb,abP->acP", saB, C_p_AB)
    E -= 2.0 * ein("acP,acP->", C_p_AA, E_p_AA)
    return -E


def exch102_k11u_4(x, ampsB):
    """exch102_k11u_4 (exch12.cc:1426-1500), the B-side mirror."""
    fB = x.foccB
    nA, nB = x.noccA, x.noccB
    t, th = ampsB["t"], ampsB["theta"]
    sAb = x.sAB[:nA, fB:nB]

    M = ein("bsgt,dsht->bgdh", th, t)
    BB1f = x.BB(1, fB, fB)
    C_p_BB = ein("bgdh,ghP->bdP", M, BB1f)

    AB1f = x.AB(1, 0, fB)
    D_p_BB = ein("ab,adP->bdP", sAb, AB1f)
    E = 2.0 * ein("bdP,bdP->", C_p_BB, D_p_BB)

    sBB = sAb.T @ sAb
    X = ein("bdP,bd->P", C_p_BB, sBB)
    E += 4.0 * float(X @ x.diagAA)

    AA1 = x.AA(1)
    C_p_BA = ein("ab,acP->bcP", sAb, AA1)
    # E_p_BB[b,d,P] = sum_c sAb[c,d] C_p_BA[b,c,P]
    E_p_BB = ein("cd,bcP->bdP", sAb, C_p_BA)
    E -= 2.0 * ein("bdP,bdP->", C_p_BB, E_p_BB)
    return -E


def exch120_k11u_5(x, ampsA):
    """exch120_k11u_5 (exch12.cc:1471-1570): theta contracted with Theta^P,
    then the exch110-like overlap chains."""
    fA = x.foccA
    nA, nB = x.noccA, x.noccB
    th = ampsA["theta"]
    saB = x.sAB[fA:nA, :nB]
    sRB = x.sAB[nA:, :nB]

    T_p_AR = ein("arcq,cqP->arP", th, ampsA["ThOV"])

    RB1 = x.RB(1)
    C_p_BR = ein("ab,arP->brP", saB, T_p_AR)
    E = ein("brP,rbP->", C_p_BR, RB1)

    AB2f = x.AB(2, fA, 0)
    C_p_AB = ein("rb,arP->abP", sRB, T_p_AR)
    E += ein("abP,abP->", AB2f, C_p_AB)

    BB1 = x.BB(1)
    C_p_BB = ein("ab,acP->bcP", saB, C_p_AB)
    E -= 2.0 * ein("bcP,bcP->", BB1, C_p_BB)

    xAR = saB @ sRB.T
    yAR = T_p_AR @ x.diagBB
    E += 4.0 * ein("ar,ar->", xAR, yAR)
    return -2.0 * E


def exch102_k11u_5(x, ampsB):
    """exch102_k11u_5 (exch12.cc:1572-1660), the B-side mirror."""
    fB = x.foccB
    nA, nB = x.noccA, x.noccB
    th = ampsB["theta"]
    sAb = x.sAB[:nA, fB:nB]
    sAS = x.sAB[:nA, nB:]

    T_p_BS = ein("bsdt,dtP->bsP", th, ampsB["ThOV"])

    AS1 = x.AS(1)
    C_p_AS = ein("ab,bsP->asP", sAb, T_p_BS)
    E = ein("asP,asP->", C_p_AS, AS1)

    # C_p_BA[a,b,P] = sum_s sAS[a,s] T_p_BS[b,s,P] (psi4 stores it [b][a]).
    AB1f = x.AB(1, 0, fB)
    C_p_BA = ein("as,bsP->abP", sAS, T_p_BS)
    E += ein("abP,abP->", AB1f, C_p_BA)

    # C_p_AA[a,c,P] = sum_b sAb[a,b] C_p_BA[c,b,P] -- psi4's
    # C_p_AA[a][a'][P] = sum_b sAB[a][foccB+b] C_p_BA[b][a'][P].
    AA1 = x.AA(1)
    C_p_AA = ein("ab,cbP->acP", sAb, C_p_BA)
    E -= 2.0 * ein("acP,acP->", AA1, C_p_AA)

    xBS = sAb.T @ sAS
    yBS = T_p_BS @ x.diagAA
    E += 4.0 * ein("bs,bs->", xBS, yBS)
    return -2.0 * E


def exch120_k11u_6(x, ampsA):
    """exch120_k11u_6 (exch12.cc:1662-1786): the 3 t.t + theta.theta square."""
    fA = x.foccA
    nA, nB = x.noccA, x.noccB
    t, th = ampsA["t"], ampsA["theta"]
    saB = x.sAB[fA:nA, :nB]
    sRB = x.sAB[nA:, :nB]

    # T[a,r,g,h] = 3 sum_{c,q} t[a,r,c,q] t[g,h,c,q]
    #            +   sum_{c,q} theta[c,r,a,q] theta[c,h,g,q]  (:1668-1684)
    # (the second is theta after OVOpVp_to_OVpOpV, squared)
    T = 3.0 * ein("arcq,ghcq->argh", t, t) \
        + ein("craq,chgq->argh", th, th)

    RR1 = x.RR(1)
    T_p_AA = ein("argh,rhP->agP", T, RR1)
    AA1f = x.AA(1, fA, fA)
    T_p_RR = ein("argh,agP->rhP", T, AA1f)

    AB2f = x.AB(2, fA, 0)
    C_p_BA = ein("cb,cgP->bgP", saB, T_p_AA)
    E = -ein("bgP,gbP->", C_p_BA, AB2f)

    BB1 = x.BB(1)
    C_p_BB = ein("gd,bgP->bdP", saB, C_p_BA)
    E += ein("bdP,bdP->", C_p_BB, BB1)

    RB1 = x.RB(1)
    C_p_BR = ein("qb,qrP->brP", sRB, T_p_RR)
    E -= ein("rbP,brP->", RB1, C_p_BR)

    D_p_BB = ein("rd,brP->bdP", sRB, C_p_BR)
    E += ein("bdP,bdP->", D_p_BB, BB1)

    sAA = saB @ saB.T
    sRR = sRB @ sRB.T
    E -= 2.0 * float(ein("agP,ag->P", T_p_AA, sAA) @ x.diagBB)
    E -= 2.0 * float(ein("rhP,rh->P", T_p_RR, sRR) @ x.diagBB)
    return -E


def exch102_k11u_6(x, ampsB):
    """exch102_k11u_6 (exch12.cc:1788-1900), the B-side mirror."""
    fB = x.foccB
    nA, nB = x.noccA, x.noccB
    t, th = ampsB["t"], ampsB["theta"]
    sAb = x.sAB[:nA, fB:nB]
    sAS = x.sAB[:nA, nB:]

    T = 3.0 * ein("bsdt,ghdt->bsgh", t, t) \
        + ein("dsbt,dhgt->bsgh", th, th)

    SS1 = x.SS(1)
    T_p_BB = ein("bsgh,shP->bgP", T, SS1)
    BB1f = x.BB(1, fB, fB)
    T_p_SS = ein("bsgh,bgP->shP", T, BB1f)

    AB1f = x.AB(1, 0, fB)
    # C_p_AB[a,g,P] = sum_b sAb[a,b] T_p_BB[b,g,P]
    C_p_AB = ein("ab,bgP->agP", sAb, T_p_BB)
    E = -ein("agP,agP->", C_p_AB, AB1f)

    AA1 = x.AA(1)
    C_p_AA = ein("cb,abP->acP", sAb, C_p_AB)
    E += ein("acP,acP->", C_p_AA, AA1)

    AS1 = x.AS(1)
    C_p_AS = ein("at,tsP->asP", sAS, T_p_SS)
    E -= ein("asP,asP->", AS1, C_p_AS)

    D_p_AA = ein("cs,asP->acP", sAS, C_p_AS)
    E += ein("acP,acP->", D_p_AA, AA1)

    sBB = sAb.T @ sAb
    sSS = sAS.T @ sAS
    E -= 2.0 * float(ein("bgP,bg->P", T_p_BB, sBB) @ x.diagAA)
    E -= 2.0 * float(ein("shP,sh->P", T_p_SS, sSS) @ x.diagAA)
    return -E


def exch12(x, ampsA, ampsB, uAR, uBS):
    """Exch12 (exch12.cc:37-97): exch111 + the K2u pieces (exch110/exch101 on
    Theta 2) + K2f + the six K11u pieces per side."""
    parts = {}
    parts["exch111"] = exch111(x, ampsA["ThOV"], ampsB["ThOV"])
    parts["k2u_A"] = exch110(x, ampsA["Th2OV"])
    parts["k2u_B"] = exch101(x, ampsB["Th2OV"])
    parts["k2f_A"] = exch120_k2f(x, ampsA, uAR)
    parts["k2f_B"] = exch102_k2f(x, ampsB, uBS)
    parts["k11u_A"] = (exch120_k11u_1(x, ampsA) + exch120_k11u_2(x, ampsA)
                       + exch120_k11u_3(x, ampsA) + exch120_k11u_4(x, ampsA)
                       + exch120_k11u_5(x, ampsA) + exch120_k11u_6(x, ampsA))
    parts["k11u_B"] = (exch102_k11u_1(x, ampsB) + exch102_k11u_2(x, ampsB)
                       + exch102_k11u_3(x, ampsB) + exch102_k11u_4(x, ampsB)
                       + exch102_k11u_5(x, ampsB) + exch102_k11u_6(x, ampsB))
    return sum(parts.values()), parts


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_sapt2(basis="6-31g", geom=None, focc=(0, 0), verbose=True):
    """All four SAPT2 terms plus Exch-Ind22, on top of check_sapt0's SAPT0."""
    if geom is None:
        geom = (s0.WATER_A, s0.WATER_B)
    c, t0 = s0.run(basis, geom=geom)
    x = Sapt2(c, foccA=focc[0], foccB=focc[1])
    if verbose:
        print(f"  factorization rank {x.ndf}, error {x.fact_err:.2e}")

    # Machinery check 1: Elst10 through the dressed diagonal vectors
    # (elst10.cc:37) must equal the independent SAPT0 formula.
    e10 = 4.0 * float(x.diagAA @ x.diagBB)
    assert abs(e10 - t0["Elst10,r"]) < 1e-10, (e10, t0["Elst10,r"])

    ampsA = monomer_amplitudes(x, "A")
    ampsB = monomer_amplitudes(x, "B")

    # Machinery check 2 (plan 0.2): the amplitudes' MP2 energy against PySCF.
    from pyscf import mp
    for side, mf, mol, amps, f in (("A", c["mfA"], c["molA"], ampsA, focc[0]),
                                   ("B", c["mfB"], c["molB"], ampsB, focc[1])):
        e2 = mp2_energy(x, amps, side)
        pt = mp.MP2(mf, frozen=f if f else None)
        pt.kernel()
        if verbose:
            print(f"  E2({side}) amplitudes {e2:+.12f}   pyscf {pt.e_corr:+.12f}"
                  f"   d {e2 - pt.e_corr:+.2e}")
        assert abs(e2 - pt.e_corr) < 1e-10

    uAR, vAR = k2f_integrals_A(x)
    uBS, vBS = k2f_integrals_B(x)

    # Machinery check 3: vAR must be minus the USAPT0-factorized potential.
    exA = s0.exch_ind_potential(c, "A", "B")
    exB = s0.exch_ind_potential(c, "B", "A")
    dA = np.abs(vAR + exA).max()
    dB = np.abs(vBS + exB).max()
    if verbose:
        print(f"  vAR vs -EX_A: {dA:.2e}   vBS vs -EX_B: {dB:.2e}")
    assert dA < 1e-10 and dB < 1e-10

    # CPHF coefficients, the same ones Ind20,r used (ind20.cc:931-944).
    chfA = s0._cphf_dense(c, "A", c["Cocc_A"].T @ c["w_B"] @ c["Cvir_A"])
    chfB = s0._cphf_dense(c, "B", c["Cocc_B"].T @ c["w_A"] @ c["Cvir_B"])

    t2 = dict(t0)
    e12, e120, e102 = elst12(x, ampsA, ampsB, chfA, chfB)
    t2["Elst12,r"] = e12
    t2["Exch11"] = exch11(x, ampsA, ampsB)
    ex12, ex12_parts = exch12(x, ampsA, ampsB, uAR, uBS)
    t2["Exch12"] = ex12
    i22, i22_A, i22_B = ind22(x, ampsA, ampsB)
    t2["Ind22"] = i22
    # Exch-Ind22 is NOT computed: it is Ind22 scaled by the ratio of the
    # second-order pair (ind22.cc:52). Guard the division psi4 does not.
    if abs(t0["Ind20,r"]) > 1e-12:
        t2["Exch-Ind22"] = t2["Ind22"] * (t0["Exch-Ind20,r"] / t0["Ind20,r"])
    else:
        t2["Exch-Ind22"] = 0.0

    # sapt2.cc:262: E(SAPT2) = E(SAPT0) + Elst12 + Exch11 + Exch12
    #                          + Ind22 + Exch-Ind22.
    t2["Total SAPT2"] = (t0["Total SAPT0"] + t2["Elst12,r"] + t2["Exch11"]
                         + t2["Exch12"] + t2["Ind22"] + t2["Exch-Ind22"])

    # The t2 -> 0 identity (plan 0.1, as corrected): with the amplitudes
    # zeroed, Elst12, Exch11 and Exch12 vanish identically, as do Ind22's
    # pieces 2-7; Ind22's piece 1 keeps its integral-driven part (see
    # ind220_side's docstring), so SAPT2 - SAPT0 collapses to exactly that.
    zero = {k: np.zeros_like(v) for k, v in ampsA.items()}
    zeroB = {k: np.zeros_like(v) for k, v in ampsB.items()}
    z_elst = elst12(x, zero, zeroB, chfA, chfB)[0]
    z_x11 = exch11(x, zero, zeroB)
    # uAR does NOT vanish with t2 (it is an HF-level potential); the K2f terms
    # vanish because they contract it with amplitude-derived objects.
    z_x12 = exch12(x, zero, zeroB, uAR, uBS)[0]
    _, z_iA, z_iB = ind22(x, zero, zeroB)
    assert z_elst == 0.0 and z_x11 == 0.0 and z_x12 == 0.0
    assert all(v == 0.0 for v in z_iA[1:]) and all(v == 0.0 for v in z_iB[1:])
    if verbose:
        print("  t2 -> 0: Elst12, Exch11, Exch12 and Ind22(2-7) vanish "
              "identically; Ind22(1) keeps its integral-driven part "
              f"({z_iA[0] + z_iB[0]:+.3e})")

    extras = {"Elst120": e120, "Elst102": e102, "exch12_parts": ex12_parts,
              "ind22_A": i22_A, "ind22_B": i22_B,
              "uAR": uAR, "uBS": uBS, "chfA": chfA, "chfB": chfB}
    return c, x, ampsA, ampsB, t2, extras


#: psi4 tests/sapt6 (water dimer, aug-cc-pVDZ, freeze_core) -- DF *and*
#: MP2-natural-orbital truncated, so a conventional implementation lands
#: within ~1e-5 Eh of these, not on them. Loose sanity only, never a pin.
PSI4_SAPT6 = {
    "Elst12,r": 0.000044903392,
    "Exch11": 0.000045589076,
    "Exch12": 0.002153870009,
    "Ind22": -0.000837186437,
}

#: sapt6's geometry (psi4 tests/sapt6/input.dat), Angstrom.
SAPT6_A = """
O  -1.551007  -0.114520   0.000000
H  -1.934259   0.762503   0.000000
H  -0.599677   0.040712   0.000000
"""
SAPT6_B = """
O   1.350625   0.111469   0.000000
H   1.680398  -0.373741  -0.758561
H   1.680398  -0.373741   0.758561
"""


def emit_json(path, basis, geom, label, focc=(0, 0)):
    import json
    _, _, _, _, t2, ex = run_sapt2(basis, geom=geom, focc=focc, verbose=True)
    out = {"system": label, "basis": basis,
           "note": "conventional four-index SAPT2; psi4 cannot produce these, "
                   "its closed-shell SAPT being density-fitted by construction",
           "terms": {k: float(v) for k, v in sorted(t2.items())},
           "elst12_split": {"Elst120": ex["Elst120"], "Elst102": ex["Elst102"]},
           "exch12_split": {k: float(v) for k, v in ex["exch12_parts"].items()}}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")
    return out


if __name__ == "__main__":
    if "--json" in sys.argv:
        dest = sys.argv[sys.argv.index("--json") + 1]
        o = emit_json(dest, "6-31g", (s0.WATER_A, s0.WATER_B),
                      "water dimer, 3.0 Angstrom along x")
        print(f"wrote {dest}")
        for k, v in o["terms"].items():
            print(f"  {k:<28} {v:+.12f}")
        sys.exit(0)

    if "--sapt6" in sys.argv:
        print("== psi4 tests/sapt6 cross-check (DF + NO refs: LOOSE) ==")
        _, _, _, _, t2, _ = run_sapt2("aug-cc-pvdz", geom=(SAPT6_A, SAPT6_B),
                                      focc=(1, 1))
        for k, ref in PSI4_SAPT6.items():
            got = t2[k]
            print(f"  {k:<12} {got:+.12f}   psi4 {ref:+.12f}   d {got - ref:+.2e}")
        sys.exit(0)

    print("== SAPT2, water dimer 6-31G (the Fortran reference case) ==")
    c, x, aA, aB, t2, ex = run_sapt2()
    for k in ("Elst10,r", "Exch10", "Ind20,r", "Exch-Ind20,r", "Disp20",
              "Exch-Disp20", "dHF(2)", "Total SAPT0",
              "Elst12,r", "Exch11", "Exch12", "Ind22", "Exch-Ind22",
              "Total SAPT2"):
        print(f"  {k:<14} {t2[k]:+.12f}")
    print("\n  Exch12 pieces:")
    for k, v in ex["exch12_parts"].items():
        print(f"    {k:<10} {v:+.12f}")
