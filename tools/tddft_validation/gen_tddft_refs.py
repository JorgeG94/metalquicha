#!/usr/bin/env python3
"""Generate the TDHF/TDDFT *unit test* gate numbers from PySCF 2.14.

What this is for, and what it is not for
----------------------------------------

``test/test_mqc_czt_tddft.f90`` pins every linear-response gate as a Fortran
parameter -- the 10x10 Tamm-Dancoff matrices of H2O/STO-3G, the excitation
energies and oscillator strengths of H2O/cc-pVDZ, the OH radical's
unrestricted roots, the water cation's Kohn-Sham ones. Those literals were
produced by five throwaway scripts, one per layer, living untracked in five
worktrees. This is the one script they became: rerun a section and the numbers
it prints paste straight into the test file.

The *validation suite's* references are a separate thing and are NOT here.
They live in ``tools/cpu_validation/gen_cpu_validation.py``, in ``pyscf_tddft``
beside ``pyscf_hessian``, because a validation reference has to be emitted
together with the deck it belongs to and the manifest entry that names it, and
every other reference driver in this repository is already in that file. The
split is by output, not by physics: manifest JSON there, Fortran literals here.
Both read the *same* basis JSON through ``bse_to_pyscf``, which is the part
that must not be duplicated and is therefore imported rather than copied.

Conventions, all of which matter at the gates these feed
---------------------------------------------------------

* **Geometries in Bohr.** PySCF 2.14 carries the CODATA 2010 Bohr radius and
  this code the 2018 one. Converting an Angstrom geometry on each side moves
  orbital energies by 5e-10, which was once the entire disagreement at a
  1e-10 gate.
* **The basis out of ``basis_sets/``**, never PySCF's own table: the two are
  rounded differently and disagree in the eighth decimal on some sets.
* **``conv_tol = 1e-15`` for Hartree-Fock**, 1e-13 for Kohn-Sham. The core
  elements of an A matrix are near 20 Hartree and move 7e-11 between 1e-13 and
  1e-15; a Kohn-Sham reference is gated at 1e-7 and asking a quadrature for
  more than that only risks hitting ``max_cycle``.
* **Grid level 5**, PySCF's, which is 90064 points on the water below against
  our 90058 -- worth 3.4e-10 in a matrix element.
* **One PySCF thread.** Reference values must not depend on the core count of
  whatever machine generated them.

Run it with the repository .venv:

    .venv/bin/python tools/tddft_validation/gen_tddft_refs.py --section all
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools" / "cpu_validation"))

from gen_cpu_validation import bse_to_pyscf  # noqa: E402

#: Bohr in one Angstrom, CODATA 2018 -- the reciprocal of what
#: ``mqc_physical_constants.F90`` calls ``BOHR_TO_ANGSTROM``. Used to carry a
#: geometry quoted in Angstrom into the Bohr the molecules below are declared
#: in, so that both codes see the identical coordinates.
BOHR_PER_ANGSTROM = 1.8897261254578281

#: The plan's water, in Bohr. C2v, no symmetry imposed, spherical.
WATER_BOHR = [
    ("O", (0.0, 0.0, 0.22259084031767759)),
    ("H", (0.0, 1.4275992706554927, -0.89036525099683572)),
    ("H", (0.0, -1.4275992706554927, -0.89036525099683572)),
]

#: The OH radical at r = 0.9697 Angstrom, in Bohr. A doublet, and a 2-Pi one:
#: see the note on ``converge_unrestricted``.
OH_BOHR = [
    ("O", (0.0, 0.0, 0.0)),
    ("H", (0.0, 0.0, 0.9697 * BOHR_PER_ANGSTROM)),
]

#: Every restricted gate runs over these four, which between them carry no
#: exact exchange, a fixed fraction of it, and a range-separated one.
FUNCTIONALS = (None, "pbe", "b3lyp", "camb3lyp")

GRID_LEVEL = 5


# --------------------------------------------------------------------------
# molecules and references
# --------------------------------------------------------------------------

def make_mol(atoms, basis, spin=0, charge=0):
    """One PySCF molecule on this repository's basis JSON, coordinates in Bohr."""
    from pyscf import gto

    mol = gto.Mole()
    mol.atom = [(s, xyz) for s, xyz in atoms]
    mol.unit = "Bohr"
    mol.basis = {s: bse_to_pyscf(basis, s) for s in sorted({s for s, _ in atoms})}
    mol.cart = False
    mol.spin = spin
    mol.charge = charge
    mol.verbose = 0
    mol.build()
    return mol


def converge(mol, functional=None, level=GRID_LEVEL):
    """A restricted reference, driven past what any gate here reads."""
    from pyscf import dft, scf

    if functional is None:
        mf = scf.RHF(mol)
        mf.conv_tol = 1e-15
        mf.conv_tol_grad = 1e-11
    else:
        mf = dft.RKS(mol)
        mf.xc = functional
        mf.grids.level = level
        mf.conv_tol = 1e-13
        mf.conv_tol_grad = 1e-9
    mf.max_cycle = 200
    mf.kernel()
    assert mf.converged, "the restricted reference did not converge"
    return mf


def converge_unrestricted(mol, functional=None, level=GRID_LEVEL, shift=True):
    """An unrestricted reference, and why the Kohn-Sham one is a two-stage solve.

    OH is a 2-Pi radical and a quadrature grid is not cylindrically symmetric,
    so the two orientations of the singly occupied pi hole are distinct
    stationary points 6e-7 Hartree apart whose soft root moves 4e-5 between
    them. Plain DIIS from ``minao`` lands in either one and stagnates near
    |g| = 1e-6 on perhaps half its runs. A first pass with a level shift picks
    the lower basin every time; the second releases the shift and converges
    from there. PySCF's own ``converged`` flag can still read False at the end
    because the gradient plateaus near 1e-9 -- what is asserted is the gradient
    reached, since that is what the orbitals inherit their error from.
    """
    from pyscf import dft, scf

    if functional is None:
        mf = scf.UHF(mol)
        mf.conv_tol = 1e-15
        mf.conv_tol_grad = 1e-11
        mf.max_cycle = 200
        mf.kernel()
        assert mf.converged, "the unrestricted Hartree-Fock reference did not converge"
        return mf

    dm0 = None
    if shift:
        first = dft.UKS(mol)
        first.xc = functional
        first.grids.level = level
        first.init_guess = "minao"
        first.level_shift = 0.5
        first.conv_tol = 1e-9
        first.conv_tol_grad = 1e-5
        first.max_cycle = 200
        first.kernel()
        dm0 = first.make_rdm1()

    mf = dft.UKS(mol)
    mf.xc = functional
    mf.grids.level = level
    mf.conv_tol = 1e-13
    mf.conv_tol_grad = 1e-9
    mf.max_cycle = 300
    mf.kernel(dm0=dm0)
    grad = np.linalg.norm(mf.get_grad(mf.mo_coeff, mf.mo_occ))
    assert grad < 5e-8, f"orbital gradient {grad:.3e} is too loose for a 1e-7 gate"
    return mf


# --------------------------------------------------------------------------
# the operators, probed rather than assembled
# --------------------------------------------------------------------------

def is_casida(mf):
    """Whether PySCF routes this reference through ``TDDFTNoHybrid``.

    A pure functional takes that path, and its ``gen_vind`` is the Casida
    matrix on a single vector rather than the paired product -- so its
    eigenvalues are omega squared outright and there is no separate B.
    """
    if not hasattr(mf, "xc"):
        return False
    return not mf._numint.libxc.is_hybrid_xc(mf.xc)


def tda_matrix(mf, singlet=True):
    """The explicit A, one column per unit vector through ``gen_vind``."""
    from pyscf import tdscf

    td = tdscf.TDA(mf)
    if singlet is not None and not isinstance(mf.mo_occ, tuple) and mf.mo_occ.ndim == 1:
        td.singlet = singlet
    vind, hdiag = td.gen_vind(mf)
    n = hdiag.size
    return vind(np.eye(n)).T


def ab_matrices(mf, singlet=True):
    """A and B out of the paired ``gen_vind``, or the Casida matrix and None.

    ``tdscf.rhf.get_ab`` ignores ``td.singlet``, so a triplet reference cannot
    come from it; the paired product does respect the flag. Feeding it
    ``(x, 0)`` returns ``(A x, -B x)``, which is where the sign on B comes from.
    """
    from pyscf import tdscf

    td = tdscf.TDDFT(mf) if hasattr(mf, "xc") else tdscf.TDHF(mf)
    restricted = not isinstance(mf.mo_occ, tuple) and mf.mo_occ.ndim == 1
    if restricted and singlet is not None:
        td.singlet = singlet
    vind, hdiag = td.gen_vind(mf)
    if is_casida(mf):
        n = hdiag.size
        return vind(np.eye(n)).T, None
    n = hdiag.size // 2
    eye = np.eye(n)
    out = vind(np.hstack([eye, np.zeros_like(eye)])).reshape(n, 2, n)
    return out[:, 0, :].T, -out[:, 1, :].T


def paired_roots(a, b):
    """The RPA excitation energies, through the same reduction the solver uses.

    Returns ``(omega, min_eig_amb)``; ``min_eig_amb`` is None for the Casida
    route, whose eigenvalues are already omega squared.
    """
    if b is None:
        return np.sqrt(np.abs(np.linalg.eigvalsh(a))), None
    wm, vm = np.linalg.eigh(a - b)
    root = vm @ np.diag(np.sqrt(np.abs(wm))) @ vm.T
    w2 = np.linalg.eigvalsh(root @ (a + b) @ root)
    return np.sqrt(np.abs(w2)), wm.min()


def iterative(mf, route, nstates=5, singlet=True):
    """An iterative solve, for the properties the dense route does not give."""
    from pyscf import tdscf

    td = tdscf.TDA(mf) if route == "tda" else tdscf.TDDFT(mf)
    restricted = not isinstance(mf.mo_occ, tuple) and mf.mo_occ.ndim == 1
    if restricted:
        td.singlet = singlet
    td.nstates = nstates
    td.conv_tol = 1e-10
    td.max_cycle = 400
    td.kernel()
    return td


# --------------------------------------------------------------------------
# emission
# --------------------------------------------------------------------------

def fortran_list(values, per_line=3, indent=26):
    """A real(dp) array literal, wrapped the way the test file writes them."""
    out = []
    for i in range(0, len(values), per_line):
        chunk = ", ".join(f"{v:.12f}_dp" for v in values[i:i + per_line])
        tail = ", &" if i + per_line < len(values) else "]"
        out.append(" " * indent + chunk + tail)
    return "\n".join(out)


def emit_matrix(name, a, indent=46):
    """One A matrix as a column-major ``reshape`` literal, plus its spectrum."""
    flat = a.reshape(-1, order="F")
    print(f"   real(dp), parameter :: {name}(N_OV, N_OV) = reshape([ &")
    lines = []
    for i in range(0, flat.size, 4):
        lines.append(" " * indent + ", ".join(f"{v:.12f}_dp" for v in flat[i:i + 4]) + ", &")
    body = "\n".join(lines)
    print(body[:-3] + "], [N_OV, N_OV])")
    w = np.linalg.eigvalsh(a)
    print(f"   real(dp), parameter :: {name}_EIG(N_OV) = [ &")
    print(fortran_list(w))


# --------------------------------------------------------------------------
# sections
# --------------------------------------------------------------------------

def section_matrices():
    """The 10x10 STO-3G Tamm-Dancoff matrices, singlet and triplet."""
    print("############ H2O / STO-3G, the explicit A matrices ############")
    mol = make_mol(WATER_BOHR, "sto-3g")
    mf = converge(mol)
    print(f"! RHF E = {mf.e_tot:.12f}")
    print("! orbital energies: " + ", ".join(f"{v:.12f}_dp" for v in mf.mo_energy))
    emit_matrix("RHF_SINGLET_A", tda_matrix(mf, singlet=True))
    emit_matrix("RHF_TRIPLET_A", tda_matrix(mf, singlet=False))
    mfp = converge(mol, "pbe")
    print(f"! PBE E = {mfp.e_tot:.12f}")
    emit_matrix("PBE_SINGLET_A", tda_matrix(mfp, singlet=True))
    emit_matrix("PBE_TRIPLET_A", tda_matrix(mfp, singlet=False))


def section_restricted():
    """The H2O/cc-pVDZ tables: four functionals, two routes, two spins."""
    print("############ H2O / cc-pVDZ, restricted ############")
    mol = make_mol(WATER_BOHR, "cc-pvdz")
    for functional in FUNCTIONALS:
        mf = converge(mol, functional)
        label = functional.upper() if functional else "RHF"
        print(f"--- {label}: E = {mf.e_tot:.12f} ---")
        for singlet in (True, False):
            spin = "singlet" if singlet else "triplet"
            a = tda_matrix(mf, singlet=singlet)
            print(f"  TDA {spin}:")
            print(fortran_list(np.linalg.eigvalsh(a)[:5]))
            omega, amb = paired_roots(*ab_matrices(mf, singlet=singlet))
            note = "" if amb is None else f"   (min eig(A-B) = {amb:.6e})"
            print(f"  RPA {spin}:{note}")
            print(fortran_list(omega[:5]))


def section_properties():
    """Transition moments, oscillator strengths, natural transition orbitals."""
    print("############ H2O / cc-pVDZ, transition properties ############")
    mol = make_mol(WATER_BOHR, "cc-pvdz")
    for functional in FUNCTIONALS:
        mf = converge(mol, functional)
        label = functional.upper() if functional else "RHF"
        for route in ("tda", "rpa"):
            td = iterative(mf, route, nstates=5)
            mu = td.transition_dipole()
            vel = td.transition_velocity_dipole()
            print(f"--- {label} {route.upper()}: E = {mf.e_tot:.12f} ---")
            print("  omega:")
            print(fortran_list(td.e))
            print("  f(length):")
            print(fortran_list(td.oscillator_strength(gauge="length")))
            print("  f(velocity):")
            print(fortran_list(td.oscillator_strength(gauge="velocity")))
            for k in range(len(td.e)):
                print(f"  S{k+1} mu  = " + ", ".join(f"{v: .12f}_dp" for v in mu[k]))
                print(f"  S{k+1} vel = " + ", ".join(f"{v: .12f}_dp" for v in vel[k]))
                print(f"  S{k+1} mu.v = {float(np.dot(mu[k], vel[k])): .12f}")
        td = iterative(mf, "tda", nstates=5)
        for state in (1, 3):
            weights, _ = td.get_nto(state=state)
            print(f"  NTO weights S{state}:")
            print(fortran_list(weights[:6]))
        break  # the NTO and moment gates are Hartree-Fock only


def section_sum_rule():
    """The whole STO-3G spectrum, and what its oscillator strengths sum to.

    Not ten. Thomas-Reiche-Kuhn holds in a complete basis and STO-3G is not
    one: the sum over every RPA root of this water is 1.946, which PySCF and
    this code agree on to 5e-11. So the gate is iterative-against-dense and
    exact, not the sum rule.
    """
    print("############ H2O / STO-3G, the whole spectrum ############")
    mol = make_mol(WATER_BOHR, "sto-3g")
    mf = converge(mol)
    nocc = int((mf.mo_occ > 0).sum())
    nvir = mf.mo_occ.size - nocc
    td = iterative(mf, "rpa", nstates=nocc * nvir)
    f = td.oscillator_strength(gauge="length")
    print(f"! RHF E = {mf.e_tot:.12f}, {len(td.e)} of {nocc*nvir} roots")
    print(f"! sum f(length) = {f.sum():.12f} against {mol.nelectron} electrons")
    print("  omega:")
    print(fortran_list(td.e))
    print("  f(length):")
    print(fortran_list(f))


def section_unrestricted():
    """OH and the water cation: roots, and the oscillator strengths of Part B.

    **The Kohn-Sham gates are not on OH, and that is measured rather than
    cautious.** Its two pi-hole orientations are separate stationary points on
    a quadrature grid, and the two codes reach different ones: PBE root 1 comes
    out 3.6e-5 from what the plan table records, against a 1e-7 gate. B3LYP is
    worse -- no guess, level shift, second-order solver or damping tried here
    brings PySCF's orbital gradient below 2.5e-6, so there is no converged
    reference to compare with at all. The water cation is the Kohn-Sham case:
    a doublet whose singly occupied orbital is non-degenerate, so there is one
    stationary point and both codes reach it from any guess. OH keeps the
    Hartree-Fock gates, where there is no grid and no basin to choose.
    """
    print("############ unrestricted ############")
    oh = make_mol(OH_BOHR, "cc-pvdz", spin=1)
    print(f"--- OH / cc-pVDZ, doublet ({oh.nao} functions) ---")
    for functional in (None,):
        mf = converge_unrestricted(oh, functional)
        label = functional.upper() if functional else "UHF"
        print(f"  {label}: E = {mf.e_tot:.12f}")
        print("  TDA roots:")
        print(fortran_list(np.linalg.eigvalsh(tda_matrix(mf))[:5]))
        omega, amb = paired_roots(*ab_matrices(mf))
        drop = int((omega <= 1e-5).sum())
        print(f"  RPA roots (dropping {drop} at the numerical zero):")
        print(fortran_list(omega[drop:drop + 5]))
        for route in ("tda", "rpa"):
            td = iterative(mf, route, nstates=6)
            print(f"  {route.upper()} iterative omega:")
            print(fortran_list(td.e))
            print(f"  {route.upper()} f(length):")
            print(fortran_list(td.oscillator_strength(gauge="length")))
            print(f"  {route.upper()} f(velocity):")
            print(fortran_list(td.oscillator_strength(gauge="velocity")))

    # The water cation, which is where the Kohn-Sham gates live: its singly
    # occupied orbital is not degenerate, so there is one stationary point and
    # both codes reach it from any guess.
    cation = make_mol(WATER_BOHR, "cc-pvdz", spin=1, charge=1)
    print(f"--- H2O+ / cc-pVDZ, doublet ({cation.nao} functions) ---")
    for functional in (None, "pbe", "b3lyp"):
        mf = converge_unrestricted(cation, functional, shift=False)
        label = functional.upper() if functional else "UHF"
        print(f"  {label}: E = {mf.e_tot:.12f}")
        print("  TDA roots:")
        print(fortran_list(np.linalg.eigvalsh(tda_matrix(mf))[:5]))
        omega, amb = paired_roots(*ab_matrices(mf))
        drop = int((omega <= 1e-5).sum())
        print(f"  RPA roots (dropping {drop} at the numerical zero):")
        print(fortran_list(omega[drop:drop + 5]))
        for route in ("tda", "rpa"):
            td = iterative(mf, route, nstates=6)
            print(f"  {route.upper()} iterative omega:")
            print(fortran_list(td.e))
            print(f"  {route.upper()} f(length):")
            print(fortran_list(td.oscillator_strength(gauge="length")))


def section_cross_check():
    """The closed shell through the unrestricted operator, A_aa +/- A_ab.

    Psi4's ``test_RU_TDA_C1``: an unrestricted operator applied to a closed
    shell has to reproduce the restricted singlet and triplet manifolds
    exactly. Stronger than any cross-code root comparison, because both sides
    come out of the same integrals.
    """
    print("############ closed shell through the unrestricted operator ############")
    mol = make_mol(WATER_BOHR, "sto-3g")
    for functional in (None, "pbe"):
        mu = converge_unrestricted(mol, functional, shift=False)
        a = tda_matrix(mu)
        half = a.shape[0] // 2
        label = functional.upper() if functional else "UHF"
        print(f"  {label}: E = {mu.e_tot:.12f}")
        print("  singlet A_aa + A_ab:")
        print(fortran_list(np.linalg.eigvalsh(a[:half, :half] + a[:half, half:])[:4]))
        print("  triplet A_aa - A_ab:")
        print(fortran_list(np.linalg.eigvalsh(a[:half, :half] - a[:half, half:])[:4]))


SECTIONS = {
    "matrices": section_matrices,
    "restricted": section_restricted,
    "properties": section_properties,
    "sum-rule": section_sum_rule,
    "unrestricted": section_unrestricted,
    "cross-check": section_cross_check,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--section", default="all",
                        choices=("all",) + tuple(SECTIONS),
                        help="which gate set to regenerate (default: all)")
    args = parser.parse_args()

    # One thread, for the reason the CPU generator pins one: a reference that
    # moves with the core count of the generating machine is not a reference.
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    from pyscf import lib as pyscf_lib

    pyscf_lib.num_threads(1)

    chosen = SECTIONS if args.section == "all" else {args.section: SECTIONS[args.section]}
    for name, run in chosen.items():
        print(f"\n======================== {name} ========================")
        run()


if __name__ == "__main__":
    main()
