#!/usr/bin/env python3
"""Compare our multipole integrals against PySCF's.

    ./build/check_multipole && python3 validation/check_multipole.py

Three questions, and they fail independently:

1. **Are the integrals right, elementwise?** Against ``mol.intor('int1e_r')``,
   ``int1e_rr`` and ``int1e_rrr`` on the same basis and the same origin. This is
   the check that catches a wrong component order or a wrong block offset -- a
   transposed quadrupole still contracts to a believable number against a
   density, so nothing downstream would object.

2. **Does the origin actually reach libcint?** It travels through ``env``, whose
   slot constants are 0-based and are not converted by the Fortran interface, so
   an off-by-one puts it in a neighbouring slot and the moments come back about
   the wrong point. The shifted-origin case is what settles it: about (0,0,0) a
   dropped origin is indistinguishable from a delivered one, and for a *neutral*
   molecule even the dipole is invariant either way. The quadrupole is not.

3. **Is the dipole right?** ``-Tr(D r)`` plus the nuclear term, against PySCF's
   own. Worth stating that this number did not previously exist on the CPU path
   at all -- only the GPU backend computed a dipole -- so this is the first
   opportunity to check it against anything.
"""

import sys
import pathlib

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tools" / "cpu_validation"))
from gen_cpu_validation import bse_to_pyscf, molecule_form, CARTESIAN  # noqa: E402

#: **Do not add a Pople basis with sp shells to this comparison without permuting.**
#: PySCF orders a Cartesian 6-31G* oxygen as `1s 2s 3s 2px..2pz 3px..3pz d`, grouping
#: all s functions then all p. Our reader hands libcint the s and p of each shared-
#: exponent group adjacently, giving `1s 2s 2px..2pz 3s 3px..3pz d`, which is also
#: GAMESS's order. So an elementwise comparison of AO-indexed matrices against PySCF
#: is only valid for a basis without sp shells -- the cases here are all like that,
#: which is why they pass.
#:
#: This was found the hard way: using PySCF's overlap with our own localized orbitals
#: reported them as non-orthonormal (`c^T S c` = 0.97, 0.75, 0.74, 0.43) when
#: `check_localize` verifies orthonormality to 2e-15 against our own overlap. The
#: metric was in the wrong order, not the orbitals.

ATOMS = [("O", (0.0, 0.0, 0.0)),
         ("H", (0.0, 0.0, 0.9584)),
         ("H", (0.9268, 0.0, -0.2400))]

CASES = [("sto-3g", (0.0, 0.0, 0.0)),
         ("cc-pvdz", (0.0, 0.0, 0.0)),
         ("cc-pvdz", (1.3, -0.7, 2.1))]

# Relative, because the octopole's matrix elements are an order of magnitude
# larger than the dipole's and an absolute bound would be a different test for
# each order.
#
# The floor is not our arithmetic. It is the basis files: BSE publishes exponents
# and coefficients to ten significant figures while PySCF's own tables carry
# more, which is worth ~1e-9 relative on any integral and is documented in
# gen_cpu_validation as looking exactly like a bug in whichever code is checked.
# Both sides here read the *same* basis_sets/ JSON, so what is left is that
# rounding propagating differently through two contraction routines.
INTEGRAL_TOL = 1e-8
DIPOLE_TOL = 1e-7


def read_dump(path):
    tokens = path.read_text().split()
    basis = tokens[0]
    values = [float(x) for x in tokens[1:]]
    nao = int(values[0])
    at = 1
    origin = np.array(values[at:at + 3]); at += 3
    dipole = np.array(values[at:at + 3]); at += 3
    energy = values[at]; at += 1
    out = {}
    for name, ncomp in (("dip", 3), ("quad", 9), ("oct", 27)):
        n = nao * nao * ncomp
        # Fortran order: (nao, nao, ncomp)
        block = np.array(values[at:at + n]).reshape((nao, nao, ncomp), order="F")
        at += n
        # to PySCF's (ncomp, nao, nao)
        out[name] = np.transpose(block, (2, 0, 1))
    out["overlap"] = np.array(values[at:at + nao * nao]).reshape((nao, nao), order="F")
    return dict(basis=basis, nao=nao, origin=origin, dipole=dipole,
                energy=energy, **out)


def main():
    from pyscf import gto, scf

    failures = 0
    for index, (_basis, _origin) in enumerate(CASES, start=1):
        path = pathlib.Path(f"/tmp/mqc_multipole_{index}.txt")
        if not path.exists():
            print(f"  MISSING {path} -- run ./build/check_multipole first")
            failures += 1
            continue
        ours = read_dump(path)
        basis, origin = ours["basis"], tuple(ours["origin"])

        symbols = {a[0] for a in ATOMS}
        mol = gto.Mole()
        mol.atom = ATOMS
        mol.unit = "Angstrom"
        mol.basis = {s: bse_to_pyscf(basis, s) for s in symbols}
        mol.cart = molecule_form(basis, symbols) == CARTESIAN
        mol.verbose = 0
        mol.build()
        mol.set_common_orig(origin)

        line = f"  {basis:8s} origin ({origin[0]:4.1f},{origin[1]:5.1f},{origin[2]:4.1f})"
        ok = True
        for name, intor in (("dip", "int1e_r"), ("quad", "int1e_rr"),
                            ("oct", "int1e_rrr")):
            theirs = mol.intor(intor)
            scale = max(np.max(np.abs(theirs)), 1.0)
            worst = np.max(np.abs(ours[name] - theirs))/scale
            line += f"  {name} {worst:.1e}"
            ok = ok and worst < INTEGRAL_TOL

        mf = scf.RHF(mol)
        mf.conv_tol = 1e-11
        mf.kernel()
        their_dipole = mf.dip_moment(unit="AU", verbose=0)
        d = float(np.max(np.abs(ours["dipole"] - their_dipole))
                  / max(float(np.max(np.abs(their_dipole))), 1.0))
        line += f"  dipole {d:.1e}"
        ok = ok and d < DIPOLE_TOL

        print(("  ok  " if ok else "  FAIL") + line)
        if not ok:
            failures += 1

    # The origin has to *matter*. If a shifted origin leaves the quadrupole
    # unchanged, it never reached libcint -- and neither the elementwise check
    # above nor any dipole would have said so, because PySCF would be asked with
    # the same origin and a neutral molecule's dipole is invariant regardless.
    a = pathlib.Path("/tmp/mqc_multipole_2.txt")
    b = pathlib.Path("/tmp/mqc_multipole_3.txt")
    if a.exists() and b.exists():
        qa, qb = read_dump(a)["quad"], read_dump(b)["quad"]
        moved = float(np.max(np.abs(qa - qb)))
        if moved < 1.0:
            print(f"  FAIL the quadrupole barely moved when the origin did "
                  f"({moved:.1e}); the origin is not reaching libcint")
            failures += 1
        else:
            print(f"  ok   shifting the origin moves the quadrupole by {moved:.3f}, "
                  f"so it is being applied")

    print()
    if failures:
        print(f"[multipole] {failures} FAILURE(S)")
        return 1
    print("[multipole] integrals, origin handling and the dipole all match PySCF")
    return 0


if __name__ == "__main__":
    sys.exit(main())
