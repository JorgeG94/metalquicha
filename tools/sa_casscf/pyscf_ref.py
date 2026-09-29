#!/usr/bin/env python3
"""A PySCF reference for state-averaged CASSCF energies and gradients.

Phase 0 of the SA-CASSCF gradient project (``SA_CASSCF_GRADIENT_PLAN.md`` at
the repository root). What this checks against later phases: every root's
energy and (with ``--grad``) every root's analytic nuclear gradient from
PySCF's own Lagrangian implementation (``pyscf.grad.sacasscf``), plus a
finite-difference arbiter for the cases where PySCF's own iterative solvers
are too loose to tell a real disagreement from noise.

**Basis, not PySCF's internal tables.** ``bse_to_pyscf``/``molecule_form``
are imported from ``tools/cpu_validation/gen_cpu_validation.py`` -- the same
conversion the CPU validation suite feeds PySCF -- because PySCF's own
tables differ from this repository's basis JSON around the eighth decimal of
an exponent, enough to look exactly like a bug in whichever code is being
checked. See that module's own note on the same point.

**Cartesian vs spherical d.** 6-31G* marks its d shell ``gto_cartesian`` in
the basis JSON (confirmed by reading ``basis_sets/6-31g_st_.json``), and
``molecule_form`` from ``gen_cpu_validation`` turns that into ``mol.cart =
True`` the same way the CPU validation generator does. That is the default
here (``--cart``/``--sph`` overrides it); getting this wrong for a Pople set
is a way to make two correct codes disagree in the fourth decimal for a
reason that has nothing to do with either gradient.

**Equal weights only.** PySCF's own SA-CASSCF gradient code
(``pyscf.grad.sacasscf.Gradients.__init__``) raises ``NotImplementedError``
outright once the weight spread exceeds ``1e-8`` -- read directly out of that
file, not assumed. ``--weights`` is still a free CLI argument for the energy
path (state averaging itself has no such restriction), but ``--grad`` with
unequal weights fails exactly as PySCF fails, which is useful to know before
metalquicha's own implementation has to decide the same thing (see
``mqc_docs/source/developer_sa_casscf.rst``, "Open questions").

**Finite differences.** A neighbouring geometry's MO coefficients are not
orthonormal in the displaced AO basis, so they cannot be handed to
``mc.kernel()`` as they are; ``mcscf.project_init_guess`` projects them first.
Each displaced point is then converged with the first-order solver and polished
with ``mc.newton()``: a root energy is not stationary in the orbitals, so
whatever orbital gradient is left over goes straight into it, and the default
solver's residual gave errors that grew as ``1/h``. After the Newton polish
PySCF still stalls at an orbital gradient of about 1e-7 on these systems, which
leaves about 1e-6 Hartree/Bohr of scatter in a ``h = 1e-3`` central difference,
with a sign that changes from coordinate to coordinate. Below that level the
finite difference cannot decide anything.

**Singlets.** Two electrons in two orbitals have a triplet among the lowest
states, and at planar C2H4 it is the second root. For a closed-shell molecule
the CI solver is ``direct_spin0``, which keeps the CI vector symmetric under
alpha/beta transposition and so excludes every odd-S state exactly. A
``fix_spin_`` penalty does the same job approximately, but it made the orbital
gradient stall near 1e-6 on twisted C2H4, and it is only used when
``--spin-shift`` is given.

Run it, e.g.::

    source ~/dev/mqc_worktrees/mqc_env.sh
    source /home/jorge/dev/metalquicha/.venv/bin/activate
    python tools/sa_casscf/pyscf_ref.py --system c2h4_planar --grad
    python tools/sa_casscf/pyscf_ref.py --system c2h4_twisted --grad \\
        --fd 0:2,1:2,2:1 --fd-root 0
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CPU_VALIDATION = HERE.parent / "cpu_validation"
if str(CPU_VALIDATION) not in sys.path:
    sys.path.insert(0, str(CPU_VALIDATION))

from gen_cpu_validation import CARTESIAN, bse_to_pyscf, molecule_form  # noqa: E402

ANGSTROM_PER_BOHR = 0.52917721092
BOHR_PER_ANGSTROM = 1.0 / ANGSTROM_PER_BOHR

# --------------------------------------------------------------------------
# geometries
# --------------------------------------------------------------------------

#: LiH and H2O are the fast unit-test cases: small enough that SA-2-CAS(2,2)
#: converges in a fraction of a second, so they are worth running on every
#: check rather than only once. Bond lengths are unremarkable textbook values,
#: not tuned to anything -- the point is speed, not chemistry.
INLINE_GEOMETRIES_ANGSTROM = {
    "lih": [
        ("Li", 0.0, 0.0, 0.0),
        ("H", 0.0, 0.0, 1.5949),
    ],
    "h2o": [
        ("O", 0.0, 0.0, 0.0),
        ("H", 0.0, -0.7572, 0.5865),
        ("H", 0.0, 0.7572, 0.5865),
    ],
}

#: C2H4 planar and twisted are read from the two xyz files next to this
#: script (Angstrom), written once by this tool rather than regenerated on
#: every run, so a geometry a plot or a NAC check depends on does not move
#: under it.
XYZ_SYSTEMS = {
    "c2h4_planar": HERE / "c2h4_planar.xyz",
    "c2h4_twisted": HERE / "c2h4_twisted.xyz",
}


def read_xyz(path):
    """A plain xyz file (Angstrom) into a list of (symbol, x, y, z)."""
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].strip())
    atoms = []
    for line in lines[2:2 + n]:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms


def geometry_bohr(system):
    """The named system's geometry, converted to Bohr once at load time.

    Everything downstream (molecule construction, finite differences) works
    in Bohr from here on, so a finite-difference step given in Bohr needs no
    further conversion and the analytic and numerical gradients are directly
    comparable in Hartree/Bohr.
    """
    if system in XYZ_SYSTEMS:
        atoms = read_xyz(XYZ_SYSTEMS[system])
    elif system in INLINE_GEOMETRIES_ANGSTROM:
        atoms = INLINE_GEOMETRIES_ANGSTROM[system]
    else:
        raise SystemExit(
            f"unknown --system '{system}'; choose one of "
            f"{sorted(set(XYZ_SYSTEMS) | set(INLINE_GEOMETRIES_ANGSTROM))}")
    return [[s, x * BOHR_PER_ANGSTROM, y * BOHR_PER_ANGSTROM, z * BOHR_PER_ANGSTROM]
            for s, x, y, z in atoms]


# --------------------------------------------------------------------------
# molecule / SA-CASSCF construction
# --------------------------------------------------------------------------

def build_mol(atoms_bohr, basis, cart, charge=0, spin=0, verbose=0):
    """A PySCF ``Mole`` from a Bohr geometry and this repository's basis JSON."""
    from pyscf import gto

    symbols = sorted({a[0] for a in atoms_bohr})
    mol = gto.Mole()
    mol.atom = [(s, (x, y, z)) for s, x, y, z in atoms_bohr]
    mol.unit = "Bohr"
    mol.basis = {s: bse_to_pyscf(basis, s) for s in symbols}
    mol.cart = cart
    mol.charge = charge
    mol.spin = spin
    mol.verbose = verbose
    mol.build()
    return mol


def resolve_cart(basis, atoms_bohr, cart_override):
    """Cartesian or spherical d, from the CLI or from the basis JSON itself."""
    if cart_override is not None:
        return cart_override
    symbols = {a[0] for a in atoms_bohr}
    return molecule_form(basis, symbols) == CARTESIAN


def parse_nelecas(text):
    """``"2"`` -> 2, ``"1,1"`` -> (1, 1) -- whatever ``mcscf.CASSCF`` itself takes."""
    if "," in text:
        na, nb = text.split(",")
        return (int(na), int(nb))
    return int(text)


def run_sacasscf(mol, ncas, nelecas, nroots, weights, mo_guess=None, ci_guess=None,
                  conv_tol=1e-12, conv_tol_grad=1e-6, spin_shift=None, verbose=0,
                  newton_polish=False, project_from=None):
    """RHF, then SA-CASSCF at tight convergence.

    PySCF's default CASSCF ``conv_tol`` (1e-7) and ``conv_tol_grad`` (~1e-4)
    would leave a finite-difference check testing the convergence criterion
    rather than the gradient.

    For a closed-shell molecule with ``spin_shift`` unset, the CI solver is
    ``direct_spin0`` (singlets and other even-S states only). With
    ``spin_shift`` set it is ``fix_spin_`` at that shift, doubled up to three
    times if a root still comes back spin-contaminated.

    ``mo_guess`` seeds the CASSCF orbitals and must already belong to
    ``mol``'s AO basis. ``project_from = (prev_mol, prev_mo)`` instead takes
    orbitals from another geometry through ``mcscf.project_init_guess``. RHF is always converged
    from its own default guess. With ``newton_polish`` the first-order result
    is refined by ``mc.newton()`` at the same geometry.
    """
    from pyscf import mcscf, scf
    from pyscf.fci import addons, direct_spin0
    from pyscf.fci.spin_op import spin_square

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        raise SystemExit("RHF did not converge")

    guess_mo = mo_guess if mo_guess is not None else mf.mo_coeff
    if project_from is not None:
        guess_mo = mcscf.project_init_guess(mcscf.CASSCF(mf, ncas, nelecas),
                                            project_from[1], prev_mol=project_from[0])
    shift = spin_shift

    def make_mc():
        mc = mcscf.CASSCF(mf, ncas, nelecas)
        if mol.spin == 0 and shift is None:
            mc.fcisolver = direct_spin0.FCI(mol)
        elif mol.spin == 0:
            mc.fcisolver = addons.fix_spin_(mc.fcisolver, ss=0, shift=shift)
        if nroots > 1 or len(weights) > 1:
            mc.state_average_(list(weights))
        else:
            mc.fcisolver.nroots = nroots
        mc.natorb = False  # keep the orbital ordering the active-space indices assume
        return mc

    for attempt in range(4):
        mc = make_mc()
        mc.conv_tol = conv_tol
        mc.conv_tol_grad = conv_tol_grad
        mc.max_cycle_macro = 200
        mc.kernel(guess_mo, ci0=ci_guess)
        if not mc.converged:
            raise SystemExit("SA-CASSCF did not converge")

        if newton_polish:
            # The Newton solver stalls near 1e-7 rather than reporting
            # convergence, so its flag is not checked.
            mn = make_mc().newton()
            mn.conv_tol = 1e-13
            mn.conv_tol_grad = 1e-9
            mn.kernel(mc.mo_coeff, ci0=mc.ci)
            mc = mn

        if mol.spin != 0 or shift is None:
            return mf, mc

        # A closed-shell reference asks for singlets; check every root actually
        # is one rather than trusting the penalty to have been strong enough.
        ci_list = mc.ci if isinstance(mc.ci, list) else [mc.ci]
        ss_values = [spin_square(c, ncas, nelecas)[0] for c in ci_list]
        if max(ss_values) < 0.1:
            return mf, mc
        shift = (shift or 0.1) * 2.0
        guess_mo = mc.mo_coeff

    raise SystemExit(
        f"could not converge every SA-CASSCF root to a singlet (S^2 = {ss_values}); "
        "raise --spin-shift by hand or accept a spin-mixed average")


def state_energies(mc, nroots):
    """Every root's total energy, ascending, whatever attribute PySCF used."""
    if hasattr(mc, "e_states") and mc.e_states is not None:
        return [float(e) for e in np.atleast_1d(mc.e_states)]
    return [float(mc.e_tot)] * nroots


def active_space_counts(mc):
    """``(n_core, n_active)`` CASSCF actually used, for the output record.

    Not a check on *which* orbitals ended up active -- that would need an AO
    population analysis this script does not do -- only a record of the
    counts, which for CAS(2,2) pins the active pair to exactly the HOMO and
    LUMO of the reference determinant (``ncore = (nelectron - nelecas) // 2``
    is the HOMO index already), so there is no orbital-selection ambiguity to
    check for that specific case. A larger active space would need AVAS or
    ``sort_mo`` and a real population check; this script's systems do not.
    """
    return mc.ncore, mc.ncas


# --------------------------------------------------------------------------
# gradients
# --------------------------------------------------------------------------

def sa_gradients(mc, nroots, conv_rtol=1e-11, conv_atol=1e-13, max_cycle=200):
    """The analytic nuclear gradient of every root, Hartree/Bohr.

    ``conv_rtol``/``conv_atol``/``max_cycle`` are ``pyscf.grad.lagrange``'s
    attribute names for the Z-vector (Lagrange multiplier) solve -- read out
    of ``lagrange.Gradients.__init__`` and ``solve_lagrange`` directly, not
    guessed. The defaults there are ``1e-7``/``1e-12``/``50``; tightened here
    for the same reason the CASSCF convergence is tightened above.
    """
    grads = []
    for state in range(nroots):
        g = mc.nuc_grad_method()
        g.conv_rtol = conv_rtol
        g.conv_atol = conv_atol
        g.max_cycle = max_cycle
        de = g.kernel(state=state)
        grads.append(np.asarray(de))
    return grads


def finite_difference_root(system, basis, cart, ncas, nelecas, nroots, weights,
                            state, coords, h, conv_tol, conv_tol_grad, spin_shift=None):
    """Central differences of one root's energy, re-optimised at each geometry.

    Every displaced point starts from the undisplaced solution: its orbitals
    projected onto the displaced basis with ``mcscf.project_init_guess``, and
    its CI vectors as the Davidson guess, so the same stationary point and the
    same root order are followed. Each point is polished with ``mc.newton()``;
    see the module docstring for the 1e-6 floor that is left.
    """
    base = geometry_bohr(system)
    _, mc0 = run_sacasscf(build_mol(base, basis, cart), ncas, nelecas, nroots,
                          weights, conv_tol=conv_tol, conv_tol_grad=conv_tol_grad,
                          spin_shift=spin_shift, newton_polish=True)

    results = {}
    for (iatom, comp) in coords:
        energies = []
        for sign in (1.0, -1.0):
            atoms = [list(a) for a in base]
            atoms[iatom][1 + comp] += sign * h
            _, mc = run_sacasscf(build_mol(atoms, basis, cart), ncas, nelecas,
                                 nroots, weights, ci_guess=mc0.ci,
                                 project_from=(mc0.mol, mc0.mo_coeff),
                                 conv_tol=conv_tol, conv_tol_grad=conv_tol_grad,
                                 spin_shift=spin_shift, newton_polish=True)
            energies.append(state_energies(mc, nroots)[state])
        results[(iatom, comp)] = (energies[0] - energies[1]) / (2.0 * h)
    return results


def parse_coord_list(text):
    """``"0:2,1:0"`` -> ``[(0, 2), (1, 0)]``, atom index then x/y/z index."""
    coords = []
    for token in text.split(","):
        atom_str, comp_str = token.split(":")
        coords.append((int(atom_str), int(comp_str)))
    return coords


# --------------------------------------------------------------------------
# NAC availability (for phase 7)
# --------------------------------------------------------------------------

def report_nac_availability():
    """Whether this PySCF has an SA-CASSCF NAC implementation, and where."""
    try:
        import pyscf.nac.sacasscf as nac_mod
        return {
            "available": True,
            "module": nac_mod.__file__,
            "class": "NonAdiabaticCouplings"
            if hasattr(nac_mod, "NonAdiabaticCouplings") else None,
        }
    except ImportError as exc:
        return {"available": False, "reason": str(exc)}


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--system", default="c2h4_planar",
                        choices=sorted(set(XYZ_SYSTEMS) | set(INLINE_GEOMETRIES_ANGSTROM)))
    parser.add_argument("--basis", default="6-31g*")
    parser.add_argument("--nroots", type=int, default=2)
    parser.add_argument("--weights", default=None,
                        help="comma-separated, defaults to equal weights over --nroots")
    parser.add_argument("--ncas", type=int, default=2)
    parser.add_argument("--nelecas", default="2", help='e.g. "2" or "1,1"')
    cart_group = parser.add_mutually_exclusive_group()
    cart_group.add_argument("--cart", dest="cart", action="store_true", default=None)
    cart_group.add_argument("--sph", dest="cart", action="store_false")
    parser.add_argument("--grad", action="store_true",
                        help="also compute every root's analytic gradient")
    parser.add_argument("--fd", default=None,
                        help='e.g. "0:2,1:2,2:1" -- atom:component pairs to '
                             "finite-difference (0-based atom, 0=x/1=y/2=z)")
    parser.add_argument("--fd-root", type=int, default=0)
    parser.add_argument("--fd-h", type=float, default=5.0e-4,
                        help="central-difference step, Bohr")
    parser.add_argument("--conv-tol", type=float, default=1.0e-12)
    parser.add_argument("--conv-tol-grad", type=float, default=1.0e-6,
                        help="orbital-gradient norm the first-order CASSCF "
                             "stops at, before the Newton polish")
    parser.add_argument("--spin-shift", type=float, default=None,
                        help="use a fix_spin_ penalty of this size instead of "
                             "the singlet-only direct_spin0 solver (closed "
                             "shell only)")
    parser.add_argument("--out", default=None, help="write JSON here instead of stdout")
    parser.add_argument("--nac-check", action="store_true",
                        help="report whether this PySCF has an SA-CASSCF NAC "
                             "implementation, and exit")
    args = parser.parse_args()

    if args.nac_check:
        print(json.dumps(report_nac_availability(), indent=2))
        return 0

    weights = ([1.0 / args.nroots] * args.nroots if args.weights is None
              else [float(w) for w in args.weights.split(",")])
    if len(weights) != args.nroots:
        raise SystemExit(f"--weights has {len(weights)} entries for --nroots {args.nroots}")
    nelecas = parse_nelecas(args.nelecas)

    atoms_bohr = geometry_bohr(args.system)
    cart = resolve_cart(args.basis, atoms_bohr, args.cart)
    mol = build_mol(atoms_bohr, args.basis, cart)

    mf, mc = run_sacasscf(mol, args.ncas, nelecas, args.nroots, weights,
                          conv_tol=args.conv_tol, conv_tol_grad=args.conv_tol_grad,
                          spin_shift=args.spin_shift, newton_polish=True)
    ncore, ncas_actual = active_space_counts(mc)
    energies = state_energies(mc, args.nroots)
    e_sa = float(np.dot(weights, energies))

    from pyscf.fci.spin_op import spin_square
    ci_list = mc.ci if isinstance(mc.ci, list) else [mc.ci]
    spin_squares = [float(spin_square(c, args.ncas, nelecas)[0]) for c in ci_list]

    out = {
        "system": args.system,
        "basis": args.basis,
        "cart": cart,
        "ncas": args.ncas,
        "nelecas": nelecas if isinstance(nelecas, int) else list(nelecas),
        "ncore": ncore,
        "nroots": args.nroots,
        "weights": weights,
        "e_rhf": float(mf.e_tot),
        "e_states": energies,
        "e_sa": e_sa,
        "casscf_converged": bool(mc.converged),
        "spin_square": spin_squares,
    }

    if args.grad:
        if max(weights) - min(weights) > 1.0e-8:
            raise SystemExit(
                "PySCF's SA-CASSCF gradient code "
                "(pyscf.grad.sacasscf.Gradients.__init__) refuses unequal "
                "weights outright; --grad needs equal weights")
        grads = sa_gradients(mc, args.nroots)
        out["gradients"] = {f"state_{i}": g.tolist() for i, g in enumerate(grads)}
        out["max_abs_gradient"] = {f"state_{i}": float(np.abs(g).max())
                                   for i, g in enumerate(grads)}

    if args.fd:
        coords = parse_coord_list(args.fd)
        fd_values = finite_difference_root(
            args.system, args.basis, cart, args.ncas, nelecas, args.nroots,
            weights, args.fd_root, coords, args.fd_h, args.conv_tol, args.conv_tol_grad,
            spin_shift=args.spin_shift)
        fd_entry = {"root": args.fd_root, "h_bohr": args.fd_h, "coords": []}
        analytic = None
        if args.grad:
            analytic = np.asarray(out["gradients"][f"state_{args.fd_root}"])
        for (iatom, comp), value in fd_values.items():
            row = {"atom": iatom, "component": comp, "finite_difference": value}
            if analytic is not None:
                row["analytic"] = float(analytic[iatom, comp])
                row["abs_diff"] = abs(row["analytic"] - value)
            fd_entry["coords"].append(row)
        out["finite_difference"] = fd_entry

    text = json.dumps(out, indent=2)
    if args.out:
        Path(args.out).write_text(text + "\n")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
