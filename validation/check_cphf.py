#!/usr/bin/env python3
"""Compare our static polarizability against PySCF, three ways.

    ./build/check_cphf && python3 validation/check_cphf.py

**The primary reference is a dense exact solve, not PySCF's own CPHF solver.**
The orbital Hessian is built here elementwise from PySCF's MO integrals and
PySCF's orbitals and inverted with LAPACK, so no iterative tolerance enters the
reference at all and every convention is written out where it can be read. Two
further numbers are reported for context:

  * `pyscf.scf.cphf.solve`, the same equations through PySCF's Krylov solver.
  * a Richardson-extrapolated finite field through PySCF's SCF, which shares no
    machinery with either -- no orbital rotation, no response function.

**PySCF's Krylov CPHF is the least accurate of the three, and that is worth
recording.** On aug-cc-pVDZ it sits 5.7e-7 from the dense solve while our CG sits
4e-9 from it, and its answer does not move when `tol` is taken from 1e-9 to 1e-15
or `max_cycle` from 50 to 200 -- so the discrepancy is not convergence and cannot
be tightened away from the outside. The finite field, entirely independent of
both, lands 3.5e-8 from the dense solve and thus agrees with us and not with it.
Three references disagreeing is only useful if you can say which one is right,
and here the dense solve is right by construction: it is a direct solve of a
180x180 system.

That ordering is why our tolerance is set against the dense solve. Comparing
against PySCF's solver instead would mean either accepting 1e-6 agreement or
"fixing" code that is already correct to 4e-9.
"""

import sys
import pathlib

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tools" / "cpu_validation"))
from gen_cpu_validation import bse_to_pyscf, molecule_form, CARTESIAN  # noqa: E402

ATOMS = [("O", (0.0, 0.0, 0.0)),
         ("H", (0.0, 0.0, 0.9584)),
         ("H", (0.9268, 0.0, -0.2400))]

#: Agreement demanded against the dense exact solve. Our CG drives the residual
#: to 1e-11 relative, so this is loose by three orders of magnitude and still
#: admits nothing real.
ALPHA_TOL = 1e-8

#: Largest field for the finite-difference reference; a halved step is also taken.
FIELD = 1e-3

#: How far the two context references may sit from the dense solve before this
#: script stops reporting and starts complaining.
#:
#: A finite field is limited by SCF convergence amplified by 1/F -- a
#: thousandfold here -- and PySCF defaults `conv_tol_grad` to sqrt(conv_tol),
#: around 3e-7, which is nowhere near enough to differentiate. That was the
#: dominant error at first, and the giveaway was that Richardson extrapolation
#: made the agreement *worse*: it cancelled a truncation term that was not
#: dominant while amplifying the noise that was. With the SCF gradient tightened
#: and the extrapolation in place, the field reference lands within 6e-8 of the
#: dense solve, and PySCF's Krylov CPHF within 6e-7.
CONTEXT_TOL = 2e-6


def read_dump(path):
    tokens = path.read_text().split()
    basis = tokens[0]
    values = [float(x) for x in tokens[1:]]
    nao, n_occ = int(values[0]), int(values[1])
    energy = values[2]
    alpha = np.array(values[3:12]).reshape((3, 3), order="F")
    return dict(basis=basis, nao=nao, n_occ=n_occ, energy=energy, alpha=alpha)


def build_mol(basis):
    from pyscf import gto

    symbols = {a[0] for a in ATOMS}
    mol = gto.Mole()
    mol.atom = ATOMS
    mol.unit = "Angstrom"
    mol.basis = {s: bse_to_pyscf(basis, s) for s in symbols}
    mol.cart = molecule_form(basis, symbols) == CARTESIAN
    mol.verbose = 0
    mol.build()
    return mol


def converged_scf(mol):
    from pyscf import scf

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-13
    mf.conv_tol_grad = 1e-11
    mf.kernel()
    return mf


def mo_dipole_ov(mol, mf, n_occ):
    """The occupied-virtual dipole blocks, `(3, n_vir, n_occ)`."""
    orbo = mf.mo_coeff[:, :n_occ]
    orbv = mf.mo_coeff[:, n_occ:]
    int_r = mol.intor_symmetric("int1e_r")
    return np.einsum("xpq,pa,qi->xai", int_r, orbv, orbo)


def dense_polarizability(mol, mf, n_occ):
    """alpha from a dense LAPACK solve of the orbital Hessian.

    The reference with no solver in it. `A` is written out term by term,

        A[ai,bj] = delta (eps_a - eps_i) + 4(ai|bj) - (ab|ij) - (aj|ib)

    inverted exactly, and contracted as `alpha = -4 h A^-1 h`. Cubic in the
    occupied-virtual dimension and quadratic in memory, so it is a reference for
    validation-sized cases and nothing else -- which is the entire reason the
    Fortran solves the same equations iteratively.

    Its symmetry is asserted rather than assumed: `A` symmetric is what makes
    conjugate gradients the right solver, and the claim is cheap to check here.
    """
    from pyscf import ao2mo

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
    asymmetry = np.abs(a - a.T).max()

    h1 = mo_dipole_ov(mol, mf, n_occ)
    u = np.linalg.solve(a, -h1.reshape(3, n_vir*n_occ).T).T
    u = u.reshape(3, n_vir, n_occ)
    return -4.0*np.einsum("xai,yai->xy", h1, u), asymmetry


def krylov_polarizability(mol, mf, n_occ):
    """alpha from pyscf.scf.cphf.solve, the same equations through its solver."""
    from pyscf.scf import cphf

    orbo = mf.mo_coeff[:, :n_occ]
    orbv = mf.mo_coeff[:, n_occ:]
    n_vir = orbv.shape[1]
    h1 = mo_dipole_ov(mol, mf, n_occ)
    vresp = mf.gen_response(hermi=1)

    def fvind(mo1):
        # The Krylov solver hands this however many trial vectors it currently
        # holds, which is not three -- hence the leading -1. The first-order
        # density carries the closed-shell factor of two, and `gen_response`
        # contracts a total density the way a Fock build does, so this returns
        # the two-electron term of the response operator with no further scaling.
        mo1 = mo1.reshape(-1, n_vir, n_occ)
        dm1 = np.einsum("xai,pa,qi->xpq", mo1, orbv, orbo)*2
        dm1 = dm1 + dm1.transpose(0, 2, 1)
        return np.einsum("xpq,pa,qi->xai", vresp(dm1), orbv, orbo).ravel()

    mo1 = cphf.solve(fvind, mf.mo_energy, mf.mo_occ, h1,
                     max_cycle=100, tol=1e-13)[0]
    return -4.0*np.einsum("xai,yai->xy", h1, mo1.reshape(3, n_vir, n_occ))


def field_polarizability(mol):
    """alpha by Richardson-extrapolated finite field through PySCF's SCF.

    Two central differences, at `F` and `F/2`, combined as
    `(4 D(F/2) - D(F))/3`. Each is `alpha + c F^2 + O(F^4)`, so the combination
    cancels the quadratic term -- the second hyperpolarizability, which grows
    quickly once there are diffuse functions.
    """
    from pyscf import scf

    int_r = mol.intor_symmetric("int1e_r")
    nuc = np.einsum("i,ix->x", mol.atom_charges(), mol.atom_coords())

    def dipole(field):
        mf = scf.RHF(mol)
        h0 = scf.hf.get_hcore(mol)
        mf.get_hcore = lambda *a, **k: h0 + np.einsum("x,xpq->pq", field, int_r)
        # Tight on the *gradient*, not just the energy: a finite difference
        # amplifies whatever density error the SCF leaves by 1/F.
        mf.conv_tol = 1e-14
        mf.conv_tol_grad = 1e-12
        mf.kernel()
        return -np.einsum("xpq,qp->x", int_r, mf.make_rdm1()) + nuc

    def difference(step):
        alpha = np.zeros((3, 3))
        for axis in range(3):
            shift = np.zeros(3)
            shift[axis] = step
            alpha[:, axis] = (dipole(shift) - dipole(-shift))/(2.0*step)
        return alpha

    return (4.0*difference(FIELD/2.0) - difference(FIELD))/3.0


def main():
    failures = 0
    for index in (1, 2, 3):
        path = pathlib.Path(f"/tmp/mqc_cphf_{index}.txt")
        if not path.exists():
            print(f"  MISSING {path} -- run ./build/check_cphf first")
            failures += 1
            continue
        ours = read_dump(path)

        mol = build_mol(ours["basis"])
        mf = converged_scf(mol)
        n_occ = mol.nelectron // 2

        dense, asymmetry = dense_polarizability(mol, mf, n_occ)
        krylov = krylov_polarizability(mol, mf, n_occ)
        by_field = field_polarizability(mol)

        worst = np.abs(ours["alpha"] - dense).max()
        ok = worst < ALPHA_TOL and asymmetry < 1e-10
        status = "ok  " if ok else "FAIL"
        iso = np.trace(ours["alpha"])/3.0
        print(f"  {status} {ours['basis']:12s} nao {ours['nao']:3d}  "
              f"alpha_iso {iso:9.6f}   vs dense {worst:8.2e}   "
              f"[PySCF krylov {np.abs(dense - krylov).max():8.2e}   "
              f"field {np.abs(dense - by_field).max():8.2e}]   "
              f"E {abs(ours['energy'] - mf.e_tot):8.2e}")

        if asymmetry >= 1e-10:
            print(f"       the orbital Hessian came out asymmetric by {asymmetry:.2e}, "
                  f"so conjugate gradients is the wrong solver for it")
        if not ok:
            failures += 1
            for k in range(3):
                print(f"       ours  {ours['alpha'][k]}")
                print(f"       dense {dense[k]}")

        # Not a failure -- these are context, and PySCF's Krylov solver is known
        # to sit further out than either the dense solve or the finite field.
        for label, value in (("krylov", krylov), ("field", by_field)):
            gap = np.abs(dense - value).max()
            if gap > CONTEXT_TOL:
                print(f"       note: PySCF's {label} route is {gap:.2e} from the "
                      f"dense solve, further than expected")

    print()
    if failures:
        print(f"[cphf] {failures} FAILURE(S)")
        return 1
    print("[cphf] our polarizabilities match an exact solve of the same equations")
    return 0


if __name__ == "__main__":
    sys.exit(main())
