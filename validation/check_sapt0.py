#!/usr/bin/env python3
"""SAPT0 reference implementation, conventional four-index.

Exists to produce reference numbers for the Fortran, at machine precision and
with exactly the conventions the Fortran will use -- psi4 cannot do that, being
density-fitted by construction.

Every formula is cited to psi4 master (2026-08-12); see SAPT_PLAN.md section 5.
"""
import numpy as np
from pyscf import gto, scf, lib

#: Named systems. `sapt1` is psi4's tests/sapt1 and is what the psi4 reference
#: values in PSI4_SAPT1 belong to; `water6-31g` is small enough for a Fortran
#: unit test and shares the geometry every EFP number in this tree is pinned to.
SYSTEMS = {}

# Ethene-ethyne, the psi4 tests/sapt1 geometry, Angstrom.
GEOM_A = """
C     0.000000    -0.667578    -2.124659
C     0.000000     0.667578    -2.124659
H     0.923621    -1.232253    -2.126185
H    -0.923621    -1.232253    -2.126185
H    -0.923621     1.232253    -2.126185
H     0.923621     1.232253    -2.126185
"""
GEOM_B = """
C     0.000000     0.000000     2.900503
C     0.000000     0.000000     1.693240
H     0.000000     0.000000     0.627352
H     0.000000     0.000000     3.963929
"""


WATER_A = """
O     0.00000000   0.00000000   0.10077199
H     0.00000000   0.77250895  -0.46780200
H     0.00000000  -0.77250895  -0.46780200
"""
WATER_B = """
O     3.00000000   0.00000000   0.10077199
H     3.00000000   0.77250895  -0.46780200
H     3.00000000  -0.77250895  -0.46780200
"""


def _atoms(block):
    out = []
    for line in block.strip().splitlines():
        f = line.split()
        out.append((f[0], (float(f[1]), float(f[2]), float(f[3]))))
    return out


def build(basis="cc-pvdz", verbose=0, geom=None):
    """The dimer and the two monomers, all in the dimer-centred basis.

    Ghost atoms carry their basis functions and no nuclear charge. PySCF spells
    that with a `ghost-` prefix on the symbol; the charge follows from the label
    and the basis has to be named for the ghost label too.
    """
    ga, gb = geom if geom else (GEOM_A, GEOM_B)
    A, B = _atoms(ga), _atoms(gb)

    def mol_of(real):
        """The dimer's atoms in the dimer's order, with one monomer ghosted.

        **The atom order must be the dimer's, always.** Building monomer B as
        "B's atoms then A's ghosts" is the obvious thing and is wrong: the AO
        ordering then differs from the dimer's, and every matrix built from one
        molecule and contracted against another is silently permuted. It shows
        up as Tr(D_B S) != n_occ_B, which is the cheapest way to catch it.
        """
        spec, bas = [], {}
        for idx, (s, x) in enumerate(A + B):
            label = s if idx in real else "ghost-" + s
            spec.append((label, x))
            bas[label] = basis
        return gto.M(atom=spec, basis=bas, unit="Angstrom", verbose=verbose)

    dimer = gto.M(atom=A + B, basis=basis, unit="Angstrom", verbose=verbose)
    molA = mol_of(set(range(len(A))))
    molB = mol_of(set(range(len(A), len(A) + len(B))))
    return dimer, molA, molB


def monomer_scf(mol, conv=1e-11):
    mf = scf.RHF(mol)
    mf.conv_tol = conv
    mf.kernel()
    assert mf.converged
    return mf


def check_m1(basis="cc-pvdz"):
    """M1: ghosts work, and the dimer basis lowers each monomer's energy (BSSE)."""
    dimer, molA, molB = build(basis)
    A, B = _atoms(GEOM_A), _atoms(GEOM_B)

    print(f"  dimer nao   {dimer.nao}")
    print(f"  molA  nao   {molA.nao}   nelec {molA.nelectron}")
    print(f"  molB  nao   {molB.nao}   nelec {molB.nelectron}")
    assert molA.nao == dimer.nao == molB.nao, "monomers must span the dimer basis"

    # Nuclear charges: only the real atoms carry one.
    zA = molA.atom_charges()
    print(f"  molA charges {zA}")
    assert zA[len(A):].sum() == 0, "ghosts must carry no nuclear charge"
    assert zA[:len(A)].sum() == sum(gto.charge(s) for s, _ in A)

    # BSSE: the same monomer in the bigger basis must be lower.
    own = gto.M(atom=A, basis=basis, unit="Angstrom", verbose=0)  # A alone, own basis
    e_own = monomer_scf(own).e_tot
    e_dcbs = monomer_scf(molA).e_tot
    print(f"  A in own basis  {e_own:.10f}")
    print(f"  A in dimer basis {e_dcbs:.10f}")
    print(f"  BSSE            {e_dcbs - e_own:+.3e}")
    assert e_dcbs < e_own, "the dimer basis must lower the monomer energy"
    return dimer, molA, molB


def cache(basis="cc-pvdz", df_aux=None, geom=None):
    """Everything the terms are built from, in psi4's naming.

    `D` carries NO factor of two (sapt_jk_terms.py:65-67) -- every factor of 2
    and 4 downstream exists because of that.
    """
    dimer, molA, molB = build(basis, geom=geom)
    # `df_aux` reproduces psi4 rather than being the production path: DF-SCF for
    # the orbitals, then the same fitting basis for J/K. It exists so a term that
    # disagrees can be shown to disagree by exactly the fitting error.
    if df_aux:
        mfA = scf.RHF(molA).density_fit(auxbasis=df_aux[0]); mfA.conv_tol = 1e-11; mfA.kernel()
        mfB = scf.RHF(molB).density_fit(auxbasis=df_aux[0]); mfB.conv_tol = 1e-11; mfB.kernel()
    else:
        mfA, mfB = monomer_scf(molA), monomer_scf(molB)

    c = {}
    c["mol"] = dimer
    c["molA"], c["molB"] = molA, molB
    c["mfA"], c["mfB"] = mfA, mfB

    for tag, mf in (("A", mfA), ("B", mfB)):
        occ = mf.mo_occ > 0
        c["Cocc_" + tag] = mf.mo_coeff[:, occ]
        c["Cvir_" + tag] = mf.mo_coeff[:, ~occ]
        c["eps_occ_" + tag] = mf.mo_energy[occ]
        c["eps_vir_" + tag] = mf.mo_energy[~occ]
        c["D_" + tag] = c["Cocc_" + tag] @ c["Cocc_" + tag].T
        c["P_" + tag] = c["Cvir_" + tag] @ c["Cvir_" + tag].T

    # Nuclear attraction of each monomer's nuclei ALONE, in the dimer basis --
    # which is exactly what the ghosted molecule's one-electron integral is.
    c["V_A"] = molA.intor("int1e_nuc")
    c["V_B"] = molB.intor("int1e_nuc")
    c["S"] = dimer.intor("int1e_ovlp")

    # J and K from each monomer's own density.
    if df_aux:
        from pyscf import df as _df
        obj = _df.DF(dimer); obj.auxbasis = df_aux[1]; obj.build()
        c["_df"] = obj
    for tag in ("A", "B"):
        if df_aux:
            j, k = c["_df"].get_jk(c["D_" + tag], hermi=1)
        else:
            j, k = scf.hf.get_jk(dimer, c["D_" + tag])
        c["J_" + tag], c["K_" + tag] = j, k

    c["w_A"] = c["V_A"] + 2.0 * c["J_A"]
    c["w_B"] = c["V_B"] + 2.0 * c["J_B"]
    c["h_A"] = c["w_A"] - c["K_A"]
    c["h_B"] = c["w_B"] - c["K_B"]

    c["E_nuc"] = dimer.energy_nuc() - molA.energy_nuc() - molB.energy_nuc()
    return c


def elst10(c):
    """Elst10,r -- sapt_jk_terms.py:131-134."""
    e = 4.0 * np.einsum("pq,pq->", c["D_B"], c["J_A"])
    e += 2.0 * np.einsum("pq,pq->", c["D_A"], c["V_B"])
    e += 2.0 * np.einsum("pq,pq->", c["D_B"], c["V_A"])
    e += c["E_nuc"]
    return e


_DF_OBJ = [None]


def _jk(mol, dm, hermi=0):
    """J and K for a possibly non-symmetric density."""
    if _DF_OBJ[0] is not None:
        return _DF_OBJ[0].get_jk(dm, hermi=hermi)
    return scf.hf.get_jk(mol, dm, hermi=hermi)


def transition(c):
    """K_O and J_O, from the inter-monomer transition density D_A S D_B.

    Built in FISAPT's orientation directly (`fisapt.cc:3665`, K_O = K[D_A S D_B]).
    psi4's python driver feeds the JK pair reversed and transposes afterwards
    (`sapt_jk_terms.py:96-111`); doing it this way round removes the single most
    bug-prone line in that file. The density is NOT symmetric, hence hermi=0.
    """
    O = c["D_A"] @ c["S"] @ c["D_B"]
    J_O, K_O = _jk(c["mol"], O, hermi=0)
    return O, J_O, K_O


def exch10_s2(c):
    """Exch10(S^2), FISAPT's MCBS form -- fisapt.cc:2450-2465.

    Preferred over the DCBS form of sapt_jk_terms.py:215-222 because it never
    mentions P, so it does not depend on D + P = S^-1 -- an identity that fails
    silently if the monomer SCF drops linearly dependent MOs.
    """
    S, D_A, D_B = c["S"], c["D_A"], c["D_B"]
    _, _, K_O = transition(c)
    AB = D_A @ S @ D_B
    BA = D_B @ S @ D_A
    # Every contraction below is psi4's `vector_dot`: elementwise, Tr(X^T Y).
    # For the two terms whose operands are both non-symmetric that is not Tr(X Y),
    # and using the wrong one puts Exch10 55% high.
    e = -2.0 * np.einsum("pq,pq->", D_A, c["K_B"])
    e -= 2.0 * np.einsum("pq,pq->", AB, c["h_A"])
    e -= 2.0 * np.einsum("pq,pq->", BA, c["h_B"])
    e += 2.0 * np.einsum("pq,pq->", D_B @ S @ D_A @ S @ D_B, c["w_A"])
    e += 2.0 * np.einsum("pq,pq->", D_A @ S @ D_B @ S @ D_A, c["w_B"])
    e -= 2.0 * np.einsum("pq,pq->", AB, K_O)
    return e


def exch10(c):
    """Exch10 at S^inf -- sapt_jk_terms.py:169-237.

    This is the one the SAPT0 total uses (sapt0.cc:231); Exch10(S^2) only feeds
    the sSAPT0 scaling ratio.
    """
    S = c["S"]
    ca, cb = c["Cocc_A"], c["Cocc_B"]
    na, nb = ca.shape[1], cb.shape[1]

    Sab = np.zeros((na + nb, na + nb))
    SAB = ca.T @ S @ cb
    Sab[:na, na:] = SAB
    Sab[na:, :na] = SAB.T
    Sab[np.diag_indices_from(Sab)] += 1.0
    Sab = np.linalg.inv(Sab)
    Sab[np.diag_indices_from(Sab)] -= 1.0      # note the -1: sapt_jk_terms.py:177

    Tmo_AA, Tmo_BB, Tmo_AB = Sab[:na, :na], Sab[na:, na:], Sab[:na, na:]
    T_A = ca @ Tmo_AA @ ca.T
    T_B = cb @ Tmo_BB @ cb.T
    T_AB = ca @ Tmo_AB @ cb.T

    # psi4 feeds these as (C_left=Cocc_B, C_right=Cocc_A Tmo_AB), so the implied
    # density is T_AB TRANSPOSED (sapt_jk_terms.py:201-202) -- not T_AB. The `.T`
    # on KT_AB in the last term below is what converts it back. Building T_AB
    # here instead is a silent 55% error in Exch10 and nothing else changes.
    JT_A, KT_A = _jk(c["mol"], ca @ Tmo_AA @ ca.T, hermi=0)
    JT_AB, KT_AB = _jk(c["mol"], T_AB.T, hermi=0)

    e = -2.0 * np.einsum("pq,pq->", c["D_A"], c["K_B"])
    e += 2.0 * np.einsum("pq,pq->", T_A, c["h_B"])
    e += 2.0 * np.einsum("pq,pq->", T_B, c["h_A"])
    e += 2.0 * np.einsum("pq,pq->", T_AB, c["h_A"] + c["h_B"])
    e += 4.0 * np.einsum("pq,pq->", T_B, JT_AB - 0.5 * KT_AB)
    e += 4.0 * np.einsum("pq,pq->", T_A, JT_AB - 0.5 * KT_AB)
    e += 4.0 * np.einsum("pq,pq->", T_B, JT_A - 0.5 * KT_A)
    e += 4.0 * np.einsum("pq,pq->", T_AB, JT_AB - 0.5 * KT_AB.T)
    return e


def exch_ind_potential(c, this, other):
    """EX for one direction, USAPT0's factorization -- usapt0.cc:1261-1315.

    Two triplet products rather than the eighteen chains of
    sapt_jk_terms.py:288-310, and verified equal to them to 4.5e-17.

    `this`/`other` are "A"/"B" or "B"/"A". Called the second way round, the only
    other change is K_O -> K_O^T, which is what makes one routine serve both.
    """
    S = c["S"]
    D_t, D_o = c["D_" + this], c["D_" + other]
    h_t, h_o = c["h_" + this], c["h_" + other]
    w_t, w_o = c["w_" + this], c["w_" + other]
    K_o = c["K_" + other]

    _, J_O, K_O = transition(c)
    if this == "B":
        K_O = K_O.T                      # the one asymmetry between the two calls

    # J of D_other S D_this S D_other
    J_P = _jk(c["mol"], D_o @ S @ D_t @ S @ D_o, hermi=0)[0]

    W = -K_o - 2.0 * J_O + K_O + 2.0 * J_P
    # NOTE: S D_this w_other -- `this`, not `other`; the outer S D_o supplies
    # the second projector. Getting it wrong is dimensionally invisible.
    T1 = -h_t + S @ D_t @ w_o + w_t @ D_o @ S - K_O.T
    W = W + S @ D_o @ T1
    T2 = -h_o + w_o @ D_t @ S - K_O
    W = W + T2 @ D_o @ S
    return c["Cocc_" + this].T @ W @ c["Cvir_" + this]


def _cphf_dense(c, tag, rhs):
    """Solve the CPHF response densely -- pyscf.scf.cphf.solve is iterative and
    carries ~5e-7 that its tolerance does not control, which is far too loose to
    be a reference. The orbital Hessian here is at most a few hundred square."""
    from pyscf import ao2mo
    co, cv = c["Cocc_" + tag], c["Cvir_" + tag]
    eo, ev = c["eps_occ_" + tag], c["eps_vir_" + tag]
    no, nv = co.shape[1], cv.shape[1]
    mol = c["mol"]

    ovov = ao2mo.general(mol, (co, cv, co, cv), compact=False).reshape(no, nv, no, nv)
    oovv = ao2mo.general(mol, (co, co, cv, cv), compact=False).reshape(no, no, nv, nv)

    # A_{ar,bs} = (e_r - e_a) d_ab d_rs + 4(ar|bs) - (ab|rs) - (as|br)
    H = 4.0 * ovov
    H -= np.einsum("abrs->arbs", oovv)
    H -= np.einsum("asbr->arbs", ovov)
    H = H.reshape(no * nv, no * nv)
    d = (ev[None, :] - eo[:, None]).ravel()
    H[np.diag_indices_from(H)] += d
    return np.linalg.solve(H, -rhs.ravel()).reshape(no, nv)


def induction(c):
    """Ind20 and Exch-Ind20, uncoupled and coupled -- sapt_jk_terms.py:344-357, 518-535.

    Uncoupled and coupled differ ONLY in how x is obtained; the contraction after
    is identical, which is why this returns both from one code path.
    """
    out = {}
    ex = {"A": exch_ind_potential(c, "A", "B"), "B": exch_ind_potential(c, "B", "A")}
    for this, other in (("A", "B"), ("B", "A")):
        co, cv = c["Cocc_" + this], c["Cvir_" + this]
        w_mo = co.T @ c["w_" + other] @ cv
        eo, ev = c["eps_occ_" + this], c["eps_vir_" + this]

        x_u = w_mo / (eo[:, None] - ev[None, :])
        x_r = _cphf_dense(c, this, w_mo)

        out[f"Ind20,u ({this})"] = 2.0 * np.einsum("ar,ar->", x_u, w_mo)
        out[f"Ind20,r ({this})"] = 2.0 * np.einsum("ar,ar->", x_r, w_mo)
        out[f"Exch-Ind20,u ({this})"] = 2.0 * np.einsum("ar,ar->", x_u, ex[this])
        out[f"Exch-Ind20,r ({this})"] = 2.0 * np.einsum("ar,ar->", x_r, ex[this])
    for k in ("Ind20,u", "Ind20,r", "Exch-Ind20,u", "Exch-Ind20,r"):
        out[k] = out[f"{k} (A)"] + out[f"{k} (B)"]
    return out


def disp20(c):
    """Disp20 -- exch-disp20.cc:184-195 (NOT disp20.cc, which is debug/NO code).

    Disp20 = 4 sum_ar sum_bs (ar|bs)^2 / (e_a + e_b - e_r - e_s)

    Chemists' notation, no antisymmetrisation, factor 4, denominator negative.
    """
    from pyscf import ao2mo
    v = ao2mo.general(
        c["mol"], (c["Cocc_A"], c["Cvir_A"], c["Cocc_B"], c["Cvir_B"]), compact=False
    )
    na, ra = c["Cocc_A"].shape[1], c["Cvir_A"].shape[1]
    nb, rb = c["Cocc_B"].shape[1], c["Cvir_B"].shape[1]
    v = v.reshape(na, ra, nb, rb)
    d = (c["eps_occ_A"][:, None, None, None] + c["eps_occ_B"][None, None, :, None]
         - c["eps_vir_A"][None, :, None, None] - c["eps_vir_B"][None, None, None, :])
    return 4.0 * np.einsum("arbs,arbs->", v * v, 1.0 / d)


def exch_disp20(c):
    """Exch-Disp20, S^2 -- FISAPT's form, fisapt.cc:4351-4733.

    The libsapt_solver version is factorised into ~20 `h`/`q` pieces and is not
    term-comparable; FISAPT's is the readable one and agrees with it to 2.6e-8
    (the difference being libsapt's Laplace denominators).

    Four extra integral transforms per monomer with modified coefficient sets,
    plus four rank-one terms built from J/K/V. Same shape as Disp20 throughout,
    so no new bottleneck.
    """
    from pyscf import ao2mo
    S = c["S"]
    D_A, D_B, P_A, P_B = c["D_A"], c["D_B"], c["P_A"], c["P_B"]
    oA, vA, oB, vB = c["Cocc_A"], c["Cvir_A"], c["Cocc_B"], c["Cvir_B"]
    _, _, K_O = transition(c)

    I = np.eye(S.shape[0])
    Cr1 = (I - D_B @ S) @ vA
    Cs1 = (I - D_A @ S) @ vB
    Ca2 = D_B @ S @ oA
    Cb2 = D_A @ S @ oB
    Cr3 = 2.0 * (D_B @ S - D_A @ S @ D_B @ S) @ vA
    Cs3 = 2.0 * (D_A @ S - D_B @ S @ D_A @ S) @ vB
    Ca4 = -2.0 * D_A @ S @ D_B @ S @ oA
    Cb4 = -2.0 * D_B @ S @ D_A @ S @ oB

    def g(*cs):
        shape = tuple(x.shape[1] for x in cs)
        return ao2mo.general(c["mol"], cs, compact=False).reshape(shape)

    # v[a,r,b,s], assembled in FISAPT's order (fisapt.cc:4706-4716)
    v = np.transpose(g(oB, Cr1, oA, Cs1), (2, 1, 0, 3))       # B_br . B_as
    v += np.transpose(g(Cb2, vA, Ca2, vB), (2, 1, 0, 3))      # C_br . C_as
    v += g(oA, vA, oB, Cs3) + g(oA, vA, Cb4, vB)              # A_ar . F_bs
    v += g(oA, Cr3, oB, vB) + g(Ca4, vA, oB, vB)              # F_ar . A_bs

    # Orbital-space overlap blocks and the Q matrices (fisapt.cc:4374-4452)
    Sas = oA.T @ S @ vB
    Sbr = oB.T @ S @ vA
    SBar = oA.T @ S @ D_B @ S @ vA
    SAbs = oB.T @ S @ D_A @ S @ vB

    J_A, K_A, J_B, K_B = c["J_A"], c["K_A"], c["J_B"], c["K_B"]
    V_A, V_B = c["V_A"], c["V_B"]
    Qbr = (2.0 * oB.T @ J_A @ vA - oB.T @ K_A @ vA + oB.T @ K_O.T @ vA
           - 2.0 * oB.T @ S @ D_A @ J_B @ vA - 2.0 * oB.T @ J_A @ D_B @ S @ vA
           - oB.T @ S @ D_A @ V_B @ vA + oB.T @ V_A @ P_B @ S @ vA)
    Qas = (2.0 * oA.T @ J_B @ vB - oA.T @ K_B @ vB + oA.T @ K_O @ vB
           - 2.0 * oA.T @ J_B @ D_A @ S @ vB - 2.0 * oA.T @ S @ D_B @ J_A @ vB
           - oA.T @ S @ D_B @ V_A @ vB + oA.T @ V_B @ P_A @ S @ vB)
    Qar = 4.0 * oA.T @ J_B @ vA + 2.0 * oA.T @ V_B @ vA
    Qbs = 4.0 * oB.T @ J_A @ vB + 2.0 * oB.T @ V_A @ vB

    v += np.einsum("br,as->arbs", Qbr, Sas)
    v += np.einsum("br,as->arbs", Sbr, Qas)
    v += np.einsum("ar,bs->arbs", Qar, SAbs)
    v += np.einsum("ar,bs->arbs", SBar, Qbs)

    # t is the plain dispersion amplitude, reused (fisapt.cc:4697)
    amp = g(oA, vA, oB, vB)
    d = (c["eps_occ_A"][:, None, None, None] + c["eps_occ_B"][None, None, :, None]
         - c["eps_vir_A"][None, :, None, None] - c["eps_vir_B"][None, None, None, :])
    return -2.0 * np.einsum("arbs,arbs->", amp / d, v)


def delta_hf(c, terms, basis="cc-pvdz"):
    """dHF(2) = E_int^HF(CP) - (Elst10 + Exch10 + Ind20,r + Exch-Ind20,r).

    sapt0.cc:229. Defined as the residual, so the sum of those four plus dHF is
    E_int^HF(CP) *by construction* -- an identity, and the sharpest check we have
    on the four terms collectively.
    """
    dimer = c["mol"]
    mfD = monomer_scf(dimer)
    e_int = mfD.e_tot - c["mfA"].e_tot - c["mfB"].e_tot
    four = (terms["Elst10,r"] + terms["Exch10"]
            + terms["Ind20,r"] + terms["Exch-Ind20,r"])
    return e_int - four, e_int


#: psi4 tests/sapt1, ethene-ethyne cc-pVDZ. Density-fitted, so a conventional
#: implementation agrees to roughly 1e-6 absolute on these, not to 1e-9.
PSI4_SAPT1 = {
    "Elst10,r": -0.003599150578,
    "Exch10": 0.003629111562,
    "Exch10(S^2)": 0.003615540809,
    "Ind20,r": -0.001179619202,
    "Exch-Ind20,r": 0.000841886061,
    "Disp20": -0.001737715247,
    "Exch-Disp20": 0.000232291505,
    "Total HF": -0.00080141016,
    "dHF(2)": -0.00049363800,
    "Total SAPT0": -0.00230683390,
}


def report(name, got):
    ref = PSI4_SAPT1[name]
    print(f"  {name:<16} {got:+.12f}   psi4 {ref:+.12f}   d {got - ref:+.2e}")


def run(basis="cc-pvdz", geom=None):
    c = cache(basis, geom=geom)
    # The ordering guard: a monomer built in its own atom order rather than the
    # dimer's gives a silently permuted AO basis, and this is what catches it.
    for t in ("A", "B"):
        n = np.trace(c["D_" + t] @ c["S"])
        assert abs(n - round(n)) < 1e-8, f"Tr(D_{t} S) = {n}, AO ordering is wrong"

    t = {}
    t["Elst10,r"] = elst10(c)
    t["Exch10(S^2)"] = exch10_s2(c)
    t["Exch10"] = exch10(c)
    t.update(induction(c))
    t["Disp20"] = disp20(c)
    t["Exch-Disp20"] = exch_disp20(c)
    t["dHF(2)"], t["Total HF"] = delta_hf(c, t, basis)
    t["Total SAPT0"] = (
        t["Elst10,r"] + t["Exch10"] + t["Ind20,r"] + t["Exch-Ind20,r"]
        + t["dHF(2)"] + t["Disp20"] + t["Exch-Disp20"]
    )
    return c, t


def emit_json(path, basis, geom, label):
    """Write the reference values a Fortran test reads.

    The point of writing them rather than printing them is that a number pasted
    into a test drifts from the script that produced it. Regenerate with
    `python validation/check_sapt0.py --json <path>`.
    """
    import json
    _, t = run(basis, geom=geom)
    out = {"system": label, "basis": basis,
           "note": "conventional four-index SAPT0; psi4 cannot produce these, "
                   "its closed-shell SAPT being density-fitted by construction",
           "terms": {k: v for k, v in sorted(t.items())}}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")
    return out


if __name__ == "__main__":
    import sys
    if "--json" in sys.argv:
        dest = sys.argv[sys.argv.index("--json") + 1]
        o = emit_json(dest, "6-31g", (WATER_A, WATER_B), "water dimer, 3.0 Angstrom along x")
        print(f"wrote {dest}")
        for k, v in o["terms"].items():
            print(f"  {k:<28} {v:+.12f}")
        sys.exit(0)

    print("== M1: ghost atoms and monomers in the dimer basis ==")
    check_m1()
    print("  OK\n")

    c, t = run()
    print("== M2-M7: terms, conventional four-index ==")
    for k in ("Elst10,r", "Exch10(S^2)", "Exch10", "Ind20,u", "Ind20,r",
              "Exch-Ind20,u", "Exch-Ind20,r", "Disp20", "Exch-Disp20",
              "Total HF", "dHF(2)", "Total SAPT0"):
        if k in PSI4_SAPT1:
            report(k, t[k])
        else:
            print(f"  {k:<16} {t[k]:+.12f}")

    print("\n== the dHF identity (exact by construction) ==")
    lhs = (t["Elst10,r"] + t["Exch10"] + t["Ind20,r"] + t["Exch-Ind20,r"] + t["dHF(2)"])
    print(f"  four terms + dHF  {lhs:+.12f}")
    print(f"  E_int^HF(CP)      {t['Total HF']:+.12f}")
    print(f"  residual          {lhs - t['Total HF']:+.2e}")
