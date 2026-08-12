#!/usr/bin/env python3
"""Check a GAMESS MAKEFP potential against PySCF, before we can produce one.

    python3 tools/efp_validation/check_sum_rules.py

**Why this exists before M4 is written.** The verification M4 gets for free is the
sum rule: translate every distributed multipole to a common origin, add, and
recover the total molecular multipoles. But a sum rule constrains the *total* and
not the distribution -- a transposed quadrupole sums to a transposed total -- so it
is worth knowing exactly how much it does pin down, and worth pinning GAMESS's
conventions with it while the reference is the only thing under test. Every
convention confirmed here is one that cannot silently differ later.

Run against `reference/*.efp` and PySCF. Needs no metalquicha build.
"""

import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from efp_format import parse_efp  # noqa: E402

#: `QUADRUPOLES` packing, cross-checked against GAMESS's reader as well as its
#: writer. The sum rule below is what makes this testable rather than asserted:
#: water's XZ is -0.477, so a wrong off-diagonal order would break the total.
#: XY and YZ vanish by symmetry here, so those two remain unconstrained -- an
#: asymmetric molecule would be needed to pin them.
QUAD_ORDER = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]

#: `POLARIZABLE POINTS` tensor order. Nine components, not six: the per-LMO
#: tensors are genuinely asymmetric and only their sum is symmetric.
POL_ORDER = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2), (1, 0), (2, 0), (2, 1)]

CASES = [
    ("water_6-31gs_boys.efp", "6-31G*",
     [("O", (0.0, 0.0, 0.0)),
      ("H", (0.0, 0.0, 0.9584)),
      ("H", (0.9268, 0.0, -0.2400))]),
]


def unpack(values, order, symmetric):
    """Components in GAMESS's order to a 3x3 matrix.

    `symmetric` mirrors the off-diagonals, and the two callers differ on it: a
    quadrupole stores 6 unique components and *must* be mirrored, while a
    polarizability stores all 9 and must not be -- the per-LMO tensors are
    genuinely asymmetric. Forgetting the mirror on the quadrupole costs half of
    every off-diagonal contribution and shows up as a 0.17 error in the second
    moment sum rule, which is how it was caught here.
    """
    m = np.zeros((3, 3))
    for c, (i, j) in enumerate(order):
        m[i, j] = values[c]
        if symmetric:
            m[j, i] = values[c]
    return m


def reference_totals(atoms, basis):
    """Total dipole, total raw second moment and exact alpha, from PySCF."""
    from pyscf import gto, scf, ao2mo

    mol = gto.Mole()
    mol.atom = atoms
    mol.unit = "Angstrom"
    # PySCF's own basis library, not ours: this compares GAMESS against PySCF and
    # should not depend on metalquicha at all. The two tables agree only to
    # PySCF's eight significant figures, which is why the sum rule below lands at
    # 1e-7 rather than at the 1e-10 the file is printed to. Fine for confirming a
    # convention; not the tolerance to reuse for comparing *our* moments.
    mol.basis = basis
    mol.cart = True          # GAMESS leaves ISPHER at its Cartesian default
    mol.verbose = 0
    mol.build()
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-13
    mf.conv_tol_grad = 1e-11
    mf.kernel()

    dm = mf.make_rdm1()
    r = mol.intor_symmetric("int1e_r")
    rr = mol.intor_symmetric("int1e_rr").reshape(3, 3, mol.nao, mol.nao)
    z, coords = mol.atom_charges(), mol.atom_coords()
    dipole = -np.einsum("xpq,qp->x", r, dm) + np.einsum("i,ix->x", z, coords)
    second = (-np.einsum("xypq,qp->xy", rr, dm)
              + np.einsum("i,ix,iy->xy", z, coords, coords))

    # Exact alpha, all occupied orbitals, by dense solve -- see
    # validation/check_cphf.py for why this rather than pyscf.scf.cphf.
    n_occ = mol.nelectron // 2
    mo = mf.mo_coeff
    n_mo = mo.shape[1]
    n_vir = n_mo - n_occ
    occ, vir = slice(0, n_occ), slice(n_occ, n_mo)
    eri = ao2mo.kernel(mol, mo, compact=False).reshape([n_mo]*4)
    a = 4.0*eri[vir, occ, vir, occ].copy()
    a -= np.einsum("abij->aibj", eri[vir, vir, occ, occ])
    a -= np.einsum("ajib->aibj", eri[vir, occ, occ, vir])
    for x in range(n_vir):
        for i in range(n_occ):
            a[x, i, x, i] += mf.mo_energy[n_occ + x] - mf.mo_energy[i]
    a = a.reshape(n_vir*n_occ, n_vir*n_occ)
    orbv, orbo = mo[:, n_occ:], mo[:, :n_occ]
    h1 = np.einsum("xpq,pa,qi->xai", r, orbv, orbo)
    u = np.linalg.solve(a, -h1.reshape(3, n_vir*n_occ).T).T.reshape(3, n_vir, n_occ)
    alpha = -4.0*np.einsum("xai,yai->xy", h1, u)
    return dict(energy=mf.e_tot, dipole=dipole, second=second, alpha=alpha,
                n_occ=n_occ)


def check(path, basis, atoms):
    doc = parse_efp(path.read_text())
    s = doc["sections"]
    ref = reference_totals(atoms, basis)
    problems = []

    R = np.array([r["values"][:3] for r in s["COORDINATES"]])
    qe = np.array([r["values"][0] for r in s["MONOPOLES"]])
    qn = np.array([r["values"][1] for r in s["MONOPOLES"]])
    mu = np.array([r["values"][:3] for r in s["DIPOLES"]])
    quad = np.array([r["values"][:6] for r in s["QUADRUPOLES"]])
    q = qe + qn

    print(f"  {path.name}   {len(R)} expansion points, "
          f"{len(s['POLARIZABLE POINTS'])} polarizable points")

    # The electron count, which catches a partitioning that loses density.
    n_elec = -qe.sum()
    expected = float(sum(_z(a[0]) for a in atoms))
    print(f"    electrons from monopoles {n_elec:.7f}  (expect {expected:.1f})   "
          f"net charge {q.sum():+.2e}")
    if abs(n_elec - expected) > 1e-6:
        problems.append("distributed monopoles do not account for every electron")

    # Bond midpoints are plain arithmetic means of their two atoms.
    for i, rec in enumerate(s["COORDINATES"]):
        label = rec["label"]
        if not label.startswith("BO"):
            continue
        hi, lo = int(label[2]) - 1, int(label[3]) - 1
        if np.abs(R[i] - 0.5*(R[hi] + R[lo])).max() > 1e-9:
            problems.append(f"{label} is not the midpoint of its two atoms")

    # M4's sum rule: translate to a common origin and add.
    dipole = (q[:, None]*R).sum(0) + mu.sum(0)
    second = np.zeros((3, 3))
    for k in range(len(R)):
        second += (q[k]*np.outer(R[k], R[k])
                   + np.outer(mu[k], R[k]) + np.outer(R[k], mu[k])
                   + unpack(quad[k], QUAD_ORDER, True))
    d_gap = np.abs(dipole - ref["dipole"]).max()
    s_gap = np.abs(second - ref["second"]).max()
    print(f"    sum rule vs PySCF: dipole {d_gap:.2e}   second moment {s_gap:.2e}")
    if d_gap > 1e-6 or s_gap > 1e-5:
        problems.append("the distributed multipoles do not sum to the total")

    # Raw, not traceless -- the Buckingham conversion happens in the consumer.
    trace = np.trace(unpack(quad[0], QUAD_ORDER, True))
    print(f"    first quadrupole trace {trace:+.4f}  "
          f"({'raw' if abs(trace) > 1e-3 else 'TRACELESS -- convention changed'})")
    if abs(trace) < 1e-3:
        problems.append("quadrupoles look traceless; the raw-moment assumption broke")

    # M5: per-LMO tensors asymmetric, their sum symmetric, core excluded.
    total = np.zeros((3, 3))
    worst_asym = 0.0
    for rec in s["POLARIZABLE POINTS"]:
        a = unpack(rec["tensor"], POL_ORDER, False)
        worst_asym = max(worst_asym, np.abs(a - a.T).max())
        total += a
    sum_asym = np.abs(total - total.T).max()
    iso_sum = np.trace(total)/3.0
    iso_ref = np.trace(ref["alpha"])/3.0
    print(f"    per-LMO asymmetry {worst_asym:.4f}, asymmetry of the sum "
          f"{sum_asym:.2e}")
    print(f"    alpha_iso: sum of LMO tensors {iso_sum:.6f}   "
          f"all-occupied exact {iso_ref:.6f}   core {iso_ref - iso_sum:+.2e}")
    if worst_asym < 1e-3:
        problems.append("per-LMO tensors are symmetric; M5 may symmetrize after all")
    if sum_asym > 1e-4:
        problems.append("the sum of the LMO tensors is not symmetric")

    n_pol = len(s["POLARIZABLE POINTS"])
    if n_pol == ref["n_occ"]:
        print("    NOTE every occupied orbital has a polarizable point -- "
              "the core is NOT excluded here")
    else:
        print(f"    core excluded: {ref['n_occ']} occupied orbitals, "
              f"{n_pol} polarizable points")

    # Polarizable points sit on the LMO centroids.
    cen = {r["label"]: np.array(r["values"]) for r in s["LMO CENTROIDS"]}
    off = max(np.abs(np.array(r["xyz"]) - cen[r["label"]]).max()
              for r in s["POLARIZABLE POINTS"])
    print(f"    polarizable points sit on the LMO centroids to {off:.1e}")
    if off > 1e-9:
        problems.append("polarizable points are not at the LMO centroids")
    return problems


def _z(symbol):
    return {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9}[symbol]


def main():
    failures = 0
    for name, basis, atoms in CASES:
        path = HERE / "reference" / name
        if not path.exists():
            print(f"  MISSING {path}")
            failures += 1
            continue
        problems = check(path, basis, atoms)
        for p in problems:
            print(f"    FAIL {p}")
        failures += len(problems)

    print()
    if failures:
        print(f"[efp sum rules] {failures} FAILURE(S)")
        return 1
    print("[efp sum rules] GAMESS's conventions are what M4 and M5 will assume")
    return 0


if __name__ == "__main__":
    sys.exit(main())
