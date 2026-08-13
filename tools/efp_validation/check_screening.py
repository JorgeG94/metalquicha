#!/usr/bin/env python3
"""Establish the charge-penetration screening form, before fitting anything.

    ./build/check_dma && python3 tools/efp_validation/check_screening.py

`SCREEN2` and `SCREEN` hold a damping exponent per expansion point, and GAMESS
obtains them by *fitting*: a Powell search over the electronic-only electrostatic
potential on a grid of fused scaled van der Waals spheres. So before writing an
optimizer it is worth knowing whether the objective being optimized is the one this
script builds -- and the way to check that, without reimplementing the search, is
the trick that pinned the Boys localization: take GAMESS's published answer and ask
whether it is stationary under *our* objective.

**What this establishes.**

  * The damping form is right in kind. With no damping at all the classical
    multipole potential misses the quantum one by 20.5 kcal/mol RMS on this grid;
    with GAMESS's exponents it misses by 2.50. A factor of eight is not something a
    wrong functional form produces.
  * `alpha` on oxygen -- the parameter that matters, since it carries most of the
    density -- is sharply stationary at GAMESS's 2.000219637. Moving it 10% either
    way raises the objective to 9.63 and 7.17. So that value is at a minimum of
    essentially this objective, not merely near one.
  * `alpha = 10.0` means the screening is off for that centre, as the plan read from
    `chgpen.src`. Perturbing `BO31`, which sits at the bound, moves the objective by
    2e-6 -- it has no effect at all, which is why the fit pushed it there.

**What this does not establish, and the honest gap.** The shallow parameters are not
minimized under this objective: `A03H` scaled by 0.9 gives 2.403 against GAMESS's
2.495, and the documented initial guess `(2,2,2,4,4)` scores 2.399, better than
GAMESS's fitted answer. So the grid or the weighting differs in detail from theirs.
Candidates, in order of likelihood: the sphere radii (this uses Bondi, GAMESS has its
own table), whether bond midpoints carry spheres of their own as well as atoms, the
angular quadrature (this uses a Fibonacci sphere, GAMESS a geodesic tessellation),
and the layer weighting (this uses 1/(layer+1), which is a guess).

None of those changes the *form*, which is what was in doubt. Fixing them is what
remains before the two sections can be emitted, and each is independently checkable
against the alphas in any of the 31 reference potentials.
"""

import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent.parent
for p in (HERE, REPO / "tools" / "cpu_validation", REPO / "validation"):
    sys.path.insert(0, str(p))
from efp_format import parse_efp  # noqa: E402

BOHR_PER_ANGSTROM = 1.0/0.52917724924

#: Bondi van der Waals radii, Angstrom. GAMESS has its own table and the difference
#: is one of the identified gaps above.
VDW = {1: 1.20, 6: 1.70, 7: 1.55, 8: 1.52}

#: The layer schedule GAMESS documents: VDWSCL=0.7, VDWINC=0.1, 25 layers.
VDW_SCALE, VDW_STEP, N_LAYER = 0.7, 0.1, 25

#: Points on the *innermost* sphere. Outer layers get more, in proportion to their
#: surface area -- see `fused_spheres`.
N_ANGULAR = 110

QUAD_PAIRS = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
HARTREE_PER_KCAL = 627.5094740631


def fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    polar = np.arccos(1.0 - 2.0*i/n)
    azimuth = np.pi*(1.0 + 5.0**0.5)*i
    return np.c_[np.cos(azimuth)*np.sin(polar),
                 np.sin(azimuth)*np.sin(polar),
                 np.cos(polar)]


def fused_spheres(charges, coords):
    """Points on scaled vdW spheres, dropping any that fall inside another.

    **Constant surface density, not constant point count.** `chgpen.src:188` sets
    the points on each layer as

        KOUNT_L = INT(SCALED**2 * (NLAYER(ILAYER+1) - NLAYER(ILAYER)))

    with `SCALED` the layer's scale relative to the first, so an outer sphere gets
    more points in proportion to its area. That is also the layer weighting: there
    is no explicit one, and the outer layers count for more because there are more
    of them.

    Getting this wrong was what left the shallow exponents unstationary. A fixed
    count per layer plus an explicit `1/(layer+1)` weight -- weighting the *inner*
    layers up, the opposite way round -- put GAMESS's fitted answer below its own
    initial guess, which a converged Powell search cannot do.
    """
    points, layers = [], []
    for layer in range(N_LAYER):
        scale = VDW_SCALE + VDW_STEP*layer
        unit = fibonacci_sphere(int((scale/VDW_SCALE)**2 * N_ANGULAR))
        for a, (z, centre) in enumerate(zip(charges, coords)):
            candidate = centre + VDW[int(z)]*BOHR_PER_ANGSTROM*scale*unit
            keep = np.ones(len(candidate), bool)
            for b, (zb, other) in enumerate(zip(charges, coords)):
                if b == a:
                    continue
                inside = VDW[int(zb)]*BOHR_PER_ANGSTROM*scale - 1e-9
                keep &= np.linalg.norm(candidate - other, axis=1) > inside
            points.append(candidate[keep])
            layers += [layer]*int(keep.sum())
    return np.vstack(points), np.array(layers)


def classical_potential(dma, grid, alpha):
    """Electronic multipole potential with the monopole damped by 1 - exp(-alpha r)."""
    total = np.zeros(len(grid))
    for k in range(dma["points"].shape[1]):
        d = grid - dma["points"][:, k]
        r = np.linalg.norm(d, axis=1)
        total += dma["electronic"][k]*(1.0 - np.exp(-alpha[k]*r))/r
        total += np.einsum("x,gx->g", dma["dipole"][:, k], d)/r**3
        q = np.zeros((3, 3))
        for c, (i, j) in enumerate(QUAD_PAIRS):
            q[i, j] = q[j, i] = dma["quadrupole"][c, k]
        total += (1.5*np.einsum("gx,xy,gy->g", d, q, d)/r**5
                  - 0.5*np.trace(q)/r**3)
    return total


def main():
    import check_dma
    from gen_cpu_validation import bse_to_pyscf, molecule_form, CARTESIAN
    from pyscf import gto, scf

    dump = pathlib.Path("/tmp/mqc_dma_1.txt")
    if not dump.exists():
        print(f"  MISSING {dump} -- run ./build/check_dma first")
        return 1
    dma = check_dma.read_dump(dump)

    atoms = check_dma.__dict__.get("ATOMS") or [
        ("O", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, 0.9584)),
        ("H", (0.9268, 0.0, -0.2400))]
    symbols = {a[0] for a in atoms}
    mol = gto.Mole()
    mol.atom = atoms
    mol.unit = "Angstrom"
    mol.basis = {s: bse_to_pyscf("6-31g*", s) for s in symbols}
    mol.cart = molecule_form("6-31g*", symbols) == CARTESIAN
    mol.verbose = 0
    mol.build()
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()

    grid, layers = fused_spheres(mol.atom_charges(), mol.atom_coords())
    quantum = -np.einsum("gpq,qp->g", mol.intor("int1e_grids", grids=grid),
                         mf.make_rdm1())

    def objective(alpha):
        # Unweighted: the layer weighting is already in the point counts.
        residual = (classical_potential(dma, grid, alpha) - quantum)*HARTREE_PER_KCAL
        return float(np.sqrt(np.mean(residual**2)))

    reference = parse_efp((HERE / "reference" / "water_6-31gs_boys.efp").read_text())
    # SCREEN2 is kept as raw text by the parser -- we do not compute it yet and its
    # records are indented one space where every other section starts in column one.
    labels, theirs = [], []
    for line in reference["sections"]["SCREEN2"]["raw"]:
        token = line.split()
        if len(token) >= 3:
            labels.append(token[0])
            theirs.append(float(token[2]))
    theirs = np.array(theirs)

    print(f"  {len(grid)} grid points over {N_LAYER} layers, "
          f"{dma['points'].shape[1]} expansion points")
    undamped = objective(np.full(len(theirs), 1.0e6))
    fitted = objective(theirs)
    print(f"        no damping            {undamped:9.4f} kcal/mol")
    print(f"        GAMESS's exponents    {fitted:9.4f} kcal/mol   "
          f"({undamped/fitted:.1f}x better)")

    failures = 0
    if fitted > 0.5*undamped:
        print("        FAIL damping barely helps, so the functional form is wrong")
        failures += 1

    # A converged search cannot land above its own starting point, so this is a
    # sharp check on the objective and costs one evaluation.
    initial = np.where(np.array([l.startswith("BO") for l in labels]), 4.0, 2.0)
    at_start = objective(initial)
    print(f"        their documented initial guess {at_start:9.4f} kcal/mol   "
          f"({'fit improved on it' if fitted < at_start else 'FIT IS WORSE'})")
    if fitted >= at_start:
        print("        FAIL GAMESS's fitted exponents score worse than the guess "
              "they started from, so this is not the objective they minimized")
        failures += 1

    print("        stationarity of each exponent, +/-10%:")
    for i, label in enumerate(labels):
        down, up = theirs.copy(), theirs.copy()
        down[i] *= 0.9
        up[i] *= 1.1
        lo, hi = objective(down), objective(up)
        # An exponent at the 10.0 bound is switched off, so the objective is flat
        # in it by construction and stationarity is neither expected nor meaningful.
        off = theirs[i] >= 9.999
        note = "  (off at the 10.0 bound)" if off else ""
        print(f"          {label:6s} {lo:9.5f} {fitted:9.5f} {hi:9.5f}{note}")
        if off:
            if abs(hi - fitted) > 1e-6 or abs(lo - fitted) > 1e-6:
                print(f"        FAIL {label} is at the bound but still affects the "
                      f"potential, so 10.0 does not mean off")
                failures += 1
        elif min(lo, hi) <= fitted:
            print(f"        FAIL {label}'s exponent is not stationary")
            failures += 1

    print()
    if failures:
        print(f"[screening] {failures} FAILURE(S)")
        return 1
    print("[screening] every one of GAMESS's fitted exponents is stationary under "
          "this objective")
    return 0


if __name__ == "__main__":
    sys.exit(main())
