#!/usr/bin/env python3
"""The analytic RHF Hessian against PySCF's, as frequencies and as a norm.

Runs `mqc` on a deck with `driver: Hessian`, then builds the same molecule in
PySCF -- from **this repository's** basis JSON, not PySCF's internal tables,
which differ in the eighth decimal on some sets and look exactly like a bug in
whichever code you are checking -- and compares the harmonic frequencies and
the Hessian's Frobenius norm.

Why those two and not the matrix. The Hessian itself is not in mqc's JSON
output, and adding a 3N by 3N block to every Hessian run to make this script
easier would be the tail wagging the dog. The Frobenius norm is one scalar that
every element contributes to, and the frequencies are what the matrix is for --
between them a wrong element has nowhere to hide, while neither requires a new
output format.

**Frequencies are only compared where PySCF returns an all-real set.** A
harmonic analysis projects out six modes, and which six is well defined only at
a stationary point; away from one the two codes can partition the eighteen
differently, and PySCF hands back a complex frequency where mqc reports a small
real one. That is a disagreement about bookkeeping, not about the Hessian --
measured on the water dimer, where the matrices agree element by element to
2.6e-9 while one projected mode comes back as 765i from PySCF and 0.11 cm-1
from mqc. So the norm is checked always and the frequencies when they mean
something, rather than pinning a geometry that happens to work.

**On the tolerances.** These are looser than they look and still much tighter
than anything physical. The two codes agree element by element to between 3e-11
and 4e-7 relative, growing with system and basis size, and that residual has
been chased: it is not the coupled-perturbed convergence on either side
(ablated at 1.9e-10 for mqc and 1.2e-12 for PySCF), not integral screening on
either side (3.6e-11 and 1.4e-12), and not the reference energy, which agrees
to 1e-13. What is left is the reference orbital energies differing at 3e-10 and
the two integral libraries not being the same code. At 4e-7 relative on a
Hessian, a 4000 cm-1 mode moves by under a thousandth of a wavenumber.
"""
import json
import pathlib
import subprocess
import sys
import tempfile

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tools" / "cpu_validation"))
from gen_cpu_validation import bse_to_pyscf  # noqa: E402

# Angstrom, as a deck would give them.
WATER = [("O", (0.0, 0.0, 0.0)),
         ("H", (0.0, 0.757, 0.5862)),
         ("H", (0.0, -0.757, 0.5862))]
DIMER = [("O", (-0.0702, -0.0227, 0.0)),
         ("H", (0.8688, 0.1657, 0.0)),
         ("H", (-0.4917, 0.8332, 0.0)),
         ("O", (2.748, 0.1409, 0.0)),
         ("H", (3.0208, -0.3971, -0.7602)),
         ("H", (3.0208, -0.3971, 0.7602))]

CASES = [("water", WATER, ["sto-3g", "6-31g", "cc-pvdz"]),
         ("dimer", DIMER, ["sto-3g", "6-31g"])]

# A thousandth of a wavenumber on modes of a few thousand, and a part in a
# million on the norm. Both are far below the accuracy of the model itself and
# far above the disagreement measured above; a real error in the assembly moves
# them by orders of magnitude, not by factors.
FREQ_TOL = 1.0e-2      # cm-1
NORM_REL_TOL = 1.0e-6


def run_mqc(binary, atoms, basis, workdir):
    """One `driver: Hessian` deck through mqc; returns frequencies and the norm."""
    xyz = workdir / "geom.xyz"
    lines = [str(len(atoms)), ""]
    lines += [f"{s}  {c[0]:.10f} {c[1]:.10f} {c[2]:.10f}" for s, c in atoms]
    xyz.write_text("\n".join(lines) + "\n")

    deck = workdir / "case.json"
    deck.write_text(json.dumps({
        "schema": {"name": "case", "version": "1.0"},
        "model": {"method": "HF", "basis": basis},
        "driver": "Hessian",
        "molecules": [{"xyz": "geom.xyz",
                       "molecular_charge": 0,
                       "molecular_multiplicity": 1}],
    }))

    proc = subprocess.run([str(binary), "case.json"], cwd=workdir,
                          capture_output=True, text=True)
    out = workdir / "output_case.json"
    if proc.returncode != 0 or not out.exists():
        raise SystemExit(f"mqc failed on {basis}:\n{proc.stdout[-2000:]}")

    # The analytic path announces itself; without that line this compared the
    # finite-difference Hessian to PySCF and would have passed while testing
    # nothing that this script exists to test.
    if "computing the analytic Hessian" not in proc.stdout:
        raise SystemExit(f"{basis}: mqc fell back to finite differences, so there "
                         f"is no analytic Hessian here to check")

    data = json.load(open(out))
    block = data[next(iter(data))]
    return (np.array(block["vibrational_analysis"]["frequencies_cm1"]),
            block["hessian_frobenius_norm"])


def run_pyscf(atoms, basis):
    from pyscf import gto, scf
    from pyscf.hessian import thermo

    mol = gto.Mole()
    mol.atom = [[s, tuple(c)] for s, c in atoms]
    mol.unit = "Angstrom"
    mol.basis = {s: bse_to_pyscf(basis, s) for s, _ in atoms}
    mol.cart = False
    mol.verbose = 0
    mol.build()

    mf = scf.RHF(mol)
    mf.conv_tol = 1e-13
    mf.kernel()

    h = mf.Hessian().kernel()               # (natm, natm, 3, 3)
    natm = mol.natm
    n = 3 * natm
    flat = np.zeros((n, n))
    for i in range(natm):
        for j in range(natm):
            flat[3 * i:3 * i + 3, 3 * j:3 * j + 3] = h[i, j]

    info = thermo.harmonic_analysis(mol, h)
    return np.asarray(info["freq_wavenumber"]), float(np.linalg.norm(flat))


def main():
    binary = REPO / "build" / "mqc"
    if not binary.exists():
        print(f"  MISSING {binary} -- build mqc first")
        return 1

    failures = 0
    for name, atoms, bases in CASES:
        for basis in bases:
            with tempfile.TemporaryDirectory() as tmp:
                ours_freq, ours_norm = run_mqc(binary, atoms, basis, pathlib.Path(tmp))
            theirs_freq, theirs_norm = run_pyscf(atoms, basis)

            dn = abs(ours_norm - theirs_norm) / theirs_norm
            ok = dn <= NORM_REL_TOL
            note = ""

            if np.iscomplexobj(theirs_freq) and np.abs(theirs_freq.imag).max() > 0:
                note = "  (frequencies skipped: not a stationary point)"
                df = None
            else:
                # mqc reports every mode including the six it projects out;
                # PySCF returns the vibrations only. Compare the tail.
                real = np.real(theirs_freq)
                k = len(real)
                df = np.abs(np.sort(ours_freq)[-k:] - np.sort(real))
                ok = ok and df.max() <= FREQ_TOL

            head = f"  {'ok  ' if ok else 'FAIL'} {name:6s} {basis:10s} "
            if df is None:
                print(head + f"|dnorm|/norm {dn:8.2e} (<= {NORM_REL_TOL:.0e}){note}")
            else:
                print(head + f"{k} modes  max |dnu| {df.max():8.2e} cm-1 "
                             f"(<= {FREQ_TOL:.0e})  |dnorm|/norm {dn:8.2e} "
                             f"(<= {NORM_REL_TOL:.0e})")
            if not ok:
                failures += 1
                if df is not None and df.max() > FREQ_TOL:
                    worst = int(np.argmax(df))
                    print(f"       mode {worst}: {np.sort(ours_freq)[-k:][worst]:.4f} vs "
                          f"{np.sort(np.real(theirs_freq))[worst]:.4f} cm-1")

    print()
    if failures:
        print(f"[hessian] {failures} case(s) disagree with PySCF")
        return 1
    print("[hessian] the analytic RHF Hessian matches PySCF's on every case")
    return 0


if __name__ == "__main__":
    sys.exit(main())
