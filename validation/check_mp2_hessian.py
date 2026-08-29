#!/usr/bin/env python3
"""The analytic MP2 Hessian against pycc's, as frequencies and as a norm.

The MP2 sibling of ``check_hessian.py``, and it compares the same two scalars
for the same reasons: the Hessian is not in mqc's JSON output, the Frobenius
norm is one number every element contributes to, and the frequencies are what
the matrix is for. Between them a wrong element has nowhere to hide.

**The reference is pycc, not PySCF.** PySCF has no analytic MP2 Hessian at all
(``pyscf.hessian`` is rhf/rks/uhf/uks only), so the RHF script's reference is
unavailable here. pycc has one, frozen-core aware, and it is what every gate on
this implementation was built against -- see ``tools/mp2_hessian_oracle/``.

**Both codes are fed this repository's basis JSON**, never psi4's internal
tables. Psi4's Pople sets carry eight significant figures where the Basis Set
Exchange bundle under ``basis_sets/`` carries ten; measured on water/6-31G that
difference alone moves the correlation Hessian by 3e-9, which is far above the
agreement this script expects and looks exactly like a bug in whichever code you
are checking. The conversion is imported from the oracle harness rather than
rewritten, because writing it twice is how the element loop got emitted once per
*atom* the first time.

Usage::

    source ~/dev/mqc_worktrees/mqc_env.sh
    validation/check_mp2_hessian.py --binary build/mqc

Requires a configured build tree: the BSE bundle is unpacked at configure time.
"""

import argparse
import importlib.util
import json
import pathlib
import subprocess
import sys
import tempfile

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent

# Import the oracle harness's basis conversion rather than duplicating it. That
# function has a history: emitting a block once per atom instead of once per
# element does not fail, it silently inflates nao with linearly dependent
# functions the SCF then drops.
_spec = importlib.util.spec_from_file_location(
    "dump_pycc", REPO / "tools" / "mp2_hessian_oracle" / "dump_pycc.py"
)
_oracle = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_oracle)

# Angstrom, as a deck gives them.
WATER = [("O", (0.0, 0.0, 0.0)),
         ("H", (0.0, 0.757, 0.5862)),
         ("H", (0.0, -0.757, 0.5862))]
# Bent off linear on purpose: a degenerate pi pair makes the two codes'
# projected modes hard to line up, and the point here is the Hessian.
HCN = [("H", (0.033, -0.022, -1.058)),
       ("C", (0.0, 0.0, 0.0)),
       ("N", (0.0, 0.0, 1.164))]
AMMONIA = [("N", (0.0, 0.0, 0.116)),
           ("H", (0.0, 0.939, -0.271)),
           ("H", (0.813, -0.470, -0.271)),
           ("H", (-0.813, -0.470, -0.271))]

# (label, atoms, basis, freeze_core). Each row buys something the row above
# cannot: a second basis, then d functions, then a core of more than one
# orbital -- every gate on the implementation was measured at water/6-31G with
# ncore at most 1, and a validation suite that only repeats that is decoration.
CASES = [
    ("water", WATER, "sto-3g", False),
    ("water", WATER, "sto-3g", True),
    ("water", WATER, "6-31g", False),
    ("water", WATER, "6-31g", True),
    ("water", WATER, "6-31g*", True),      # Cartesian d on oxygen
    ("ammonia", AMMONIA, "6-31g", True),   # a four-atom case
    ("hcn", HCN, "6-31g", True),           # two core orbitals
]

# The frequencies are compared where both sides return an all-real set. A
# harmonic analysis projects out six modes and which six is well defined only at
# a stationary point; away from one the two codes can partition them
# differently and one hands back an imaginary frequency where the other reports
# a small real one. That is bookkeeping, not the Hessian -- so the norm is
# checked always and the frequencies when they mean something, which is the same
# arrangement `check_hessian.py` settled on.
NORM_TOL = 1.0e-7
FREQ_TOL = 1.0e-2   # cm^-1

# mqc's masses, from `src/core/mqc_elements.f90` -- standard atomic weights.
# psi4 defaults to the most abundant *isotope* instead (O at 15.9949 against
# 15.999), and frequencies go as 1/sqrt(m), so that 2.6e-4 in the mass shows up
# as 1.3e-4 in every frequency: 0.43 cm-1 on an O-H stretch, measured. The
# Hessians agree to 5e-9 while the frequencies disagree by forty times the
# tolerance, which reads as a broken derivative and is a mass table. PySCF
# defaults to the same standard weights, which is why `check_hessian.py` never
# had to do this.
MQC_MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999}


def run_mqc(binary, atoms, basis, freeze_core, workdir):
    """One MP2 `driver: Hessian` deck through mqc; frequencies and the norm."""
    xyz = workdir / "geom.xyz"
    lines = [str(len(atoms)), ""]
    lines += [f"{s}  {c[0]:.10f} {c[1]:.10f} {c[2]:.10f}" for s, c in atoms]
    xyz.write_text("\n".join(lines) + "\n")

    deck = workdir / "case.json"
    deck.write_text(json.dumps({
        "schema": {"name": "case", "version": "1.0"},
        "model": {"method": "mp2", "basis": basis},
        "keywords": {
            "scf": {"maxiter": 200, "tolerance": 1e-12},
            "correlation": {"freeze_core": freeze_core},
        },
        "driver": "Hessian",
        "molecules": [{"xyz": "geom.xyz",
                       "molecular_charge": 0,
                       "molecular_multiplicity": 1}],
    }))

    proc = subprocess.run([str(binary), "case.json"], cwd=workdir,
                          capture_output=True, text=True)
    out = workdir / "output_case.json"
    label = f"{basis}{' fc' if freeze_core else ''}"
    if proc.returncode != 0 or not out.exists():
        raise SystemExit(f"mqc failed on {label}:\n{proc.stdout[-2000:]}")

    # The analytic path announces itself. Without this the script would compare
    # the semi-numerical Hessian to pycc and pass while testing nothing it
    # exists to test -- and the frozen-core deck is exactly the one that used to
    # take that fallback.
    if "computing the analytic MP2 Hessian" not in proc.stdout:
        raise SystemExit(f"{label}: mqc did not take the analytic MP2 path, so "
                         f"there is no analytic Hessian here to check")

    data = json.load(open(out))
    block = data[next(iter(data))]
    return (np.array(block["vibrational_analysis"]["frequencies_cm1"]),
            block["hessian_frobenius_norm"])


def run_pycc(atoms, basis, freeze_core):
    """The same molecule through pycc; frequencies and the norm."""
    import psi4
    import pycc
    from pycc import vibanalysis

    psi4.core.be_quiet()
    symbols = [s for s, _ in atoms]
    basis_json = REPO / "basis_sets" / f"{_basis_stem(basis)}.json"
    if not basis_json.exists():
        raise SystemExit(f"{basis_json} not found -- configure the build first")
    gbs = _oracle.bse_gbs(basis_json, symbols)

    psi4.core.clean()
    psi4.core.clean_options()
    geom = ["units angstrom", "symmetry c1", "no_com", "no_reorient"]
    geom += [f"{s} {c[0]:.10f} {c[1]:.10f} {c[2]:.10f}" for s, c in atoms]
    molecule = psi4.geometry("\n".join(geom) + "\n")
    for index, (symbol, _) in enumerate(atoms):
        molecule.set_mass(index, MQC_MASS[symbol])
    psi4.basis_helper(gbs, name="mqcbse")
    psi4.set_options({"scf_type": "pk",
                      "freeze_core": "true" if freeze_core else "false",
                      "e_convergence": 1e-12, "d_convergence": 1e-12})
    _, wfn = psi4.energy("scf", return_wfn=True)
    mp = pycc.MPwfn(wfn, orbital_basis="spatial")
    mp.compute_energy()
    drv = pycc.MPderiv(mp)
    natom = len(atoms)
    hess = np.asarray(drv.hessian().total).reshape(3 * natom, 3 * natom)

    # `source` here is the molecule, not the driver, and projection is off by
    # pycc's default where mqc's vibrational analysis projects the six
    # rigid-body modes out. Ask for it, or this compares two different
    # conventions and blames the Hessian.
    analysis = vibanalysis.harmonic_analysis(
        molecule, hessian=hess, project_trans=True, project_rot=True,
    )
    return (np.asarray(analysis["frequencies"]), float(np.linalg.norm(hess)),
            int(analysis["n_tr"]))


def _basis_stem(basis):
    """`6-31g*` is `6-31g_st_` on disk, the way the bundle spells it."""
    return basis.lower().replace("*", "_st_").replace("(", "(").replace("+", "+")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--binary", default=str(REPO / "build" / "mqc"))
    parser.add_argument("--filter", default="", help="substring of the label")
    args = parser.parse_args()

    binary = pathlib.Path(args.binary)
    if not binary.exists():
        raise SystemExit(f"{binary} not found -- build first")

    bad = 0
    for name, atoms, basis, fc in CASES:
        label = f"{name}/{basis}{' frozen' if fc else ''}"
        if args.filter and args.filter not in label:
            continue
        with tempfile.TemporaryDirectory() as tmp:
            mqc_freq, mqc_norm = run_mqc(binary, atoms, basis, fc, pathlib.Path(tmp))
        ref_freq, ref_norm, n_tr = run_pycc(atoms, basis, fc)

        dnorm = abs(mqc_norm - ref_norm)
        line = f"  {label:24s} |H| {mqc_norm:14.9f} vs {ref_norm:14.9f}  d {dnorm:.2e}"
        ok = dnorm < NORM_TOL
        if len(mqc_freq) == len(ref_freq):
            # Compare the genuine vibrations only. Both codes project the
            # rigid-body modes out, but what is left of a projected mode is
            # residue, not a frequency: measured here the largest disagreement
            # in the whole set sat on a mode at 0.0 cm-1 while every real
            # vibration agreed far inside tolerance. Drop the `n_tr` smallest
            # by magnitude -- pycc's own count of what it projected -- rather
            # than widening the tolerance until the residue fits under it.
            a = np.sort(mqc_freq[np.argsort(np.abs(mqc_freq))][n_tr:])
            b = np.sort(ref_freq[np.argsort(np.abs(ref_freq))][n_tr:])
            worst = int(np.argmax(np.abs(a - b)))
            dfreq = float(np.abs(a - b)[worst])
            line += f"   dnu {dfreq:8.2e} cm-1 @ {b[worst]:7.1f} ({len(a)} modes)"
            ok = ok and dfreq < FREQ_TOL
        else:
            line += "   dnu (skipped, mode counts differ)"
        print(line + ("" if ok else "   <== FAIL"))
        bad += 0 if ok else 1

    print(f"\n{'all MP2 Hessian checks passed' if not bad else f'{bad} case(s) failed'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
