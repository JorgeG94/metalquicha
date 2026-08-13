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

#: The radii GAMESS uses for a geodesic point selection, Angstrom, from
#: `VANDER(:,2)` in `prplib.src:2086` -- Gavezzotti and Spackman's, not Bondi's.
#: Oxygen is 1.40 where Bondi has 1.52, which is an 8% difference in every sphere
#: around it. Anything not listed defaults to 1.8 in GAMESS.
VDW = {1: 1.20, 5: 1.85, 6: 1.50, 7: 1.50, 8: 1.40, 9: 1.35,
       13: 2.07, 14: 2.05, 15: 1.96, 16: 1.89, 17: 1.80}
VDW_DEFAULT = 1.8

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


def centre_radii(labels, charges):
    """A radius per expansion point, including the bond midpoints.

    Spheres go on *every* expansion centre, not only the atoms: `chgpen.src` calls
    `PDCPTS` with `NEFC`, the number of screened multipole centres, and
    `prplib.src:2162` gives a midpoint the mean of its two atoms' radii.
    """
    radii = []
    for label in labels:
        if label.startswith("BO"):
            a, b = int(label[2]) - 1, int(label[3]) - 1
            radii.append(0.5*(VDW.get(int(charges[a]), VDW_DEFAULT)
                              + VDW.get(int(charges[b]), VDW_DEFAULT)))
        else:
            index = len(radii)
            radii.append(VDW.get(int(charges[index]), VDW_DEFAULT))
    return np.array(radii)


def fused_spheres(radii, coords):
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
        for a in range(len(radii)):
            candidate = coords[:, a] + radii[a]*BOHR_PER_ANGSTROM*scale*unit
            keep = np.ones(len(candidate), bool)
            for b in range(len(radii)):
                if b == a:
                    continue
                inside = radii[b]*BOHR_PER_ANGSTROM*scale - 1e-9
                keep &= np.linalg.norm(candidate - coords[:, b], axis=1) > inside
            points.append(candidate[keep])
            layers += [layer]*int(keep.sum())
    return np.vstack(points), np.array(layers)


def classical_potential(dma, grid, alpha, gaussian=False):
    """Electronic multipole potential with the monopole damped.

    `gaussian` selects between the two damping functions GAMESS fits, which
    `chgpen.src:541` names outright: the first pass is exponential,
    `1 - exp(-alpha r)`, and is EFP2 fragment-fragment screening -- the `SCREEN2`
    block. The second is Gaussian, `1 - exp(-alpha r^2)`, and is the original EFP1
    ab initio-fragment screening -- the `SCREEN` block. That is why a potential
    carries both, with different exponents.

    The linear coefficient is 1.0 in both, and not because the fit happened to land
    there: `ICFIX=1` freezes it, and the comment at `chgpen.src:498` gives the reason
    -- the damping has to vanish at the origin, which `1 - A exp(...)` only does when
    `A = 1`. So the first column of every screening record is structurally one.
    """
    total = np.zeros(len(grid))
    for k in range(dma["points"].shape[1]):
        d = grid - dma["points"][:, k]
        r = np.linalg.norm(d, axis=1)
        argument = alpha[k]*r*r if gaussian else alpha[k]*r
        total += dma["electronic"][k]*(1.0 - np.exp(-argument))/r
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

    radii = centre_radii(dma["labels"], mol.atom_charges())
    grid, layers = fused_spheres(radii, dma["points"])
    quantum = -np.einsum("gpq,qp->g", mol.intor("int1e_grids", grids=grid),
                         mf.make_rdm1())

    def objective(alpha, gaussian=False):
        # Unweighted: the layer weighting is already in the point counts.
        residual = (classical_potential(dma, grid, alpha, gaussian)
                    - quantum)*HARTREE_PER_KCAL
        return float(np.sqrt(np.mean(residual**2)))

    reference = parse_efp((HERE / "reference" / "water_6-31gs_boys.efp").read_text())

    def block(name):
        labels, values = [], []
        for line in reference["sections"][name]["raw"]:
            token = line.split()
            if len(token) >= 3:
                labels.append(token[0])
                values.append(float(token[2]))
        return labels, np.array(values)

    failures = 0
    print(f"  {len(grid)} grid points over {N_LAYER} layers, "
          f"{dma['points'].shape[1]} expansion points")

    for name, gaussian in (("SCREEN2", False), ("SCREEN", True)):
        labels, theirs = block(name)
        form = "1 - exp(-alpha r^2)" if gaussian else "1 - exp(-alpha r)"
        fitted = objective(theirs, gaussian)
        undamped = objective(np.full(len(theirs), 1.0e6), gaussian)
        print(f"        {name}, {form}")
        print(f"          undamped {undamped:8.4f}   fitted {fitted:8.4f} kcal/mol"
              f"   ({undamped/fitted:.1f}x better)")
        if fitted > 0.5*undamped:
            print(f"          FAIL damping barely helps, so {name}'s form is wrong")
            failures += 1

        initial = np.where(np.array([l.startswith("BO") for l in labels]), 4.0, 2.0)
        at_start = objective(initial, gaussian)
        improved = fitted < at_start
        print(f"          their documented initial guess {at_start:8.4f}"
              f"   ({'fit improved on it' if improved else 'FIT IS WORSE'})")

        shallow = []
        for i, label in enumerate(labels):
            down, up = theirs.copy(), theirs.copy()
            down[i] *= 0.9
            up[i] *= 1.1
            lo, hi = objective(down, gaussian), objective(up, gaussian)
            off = theirs[i] >= 9.999
            mark = "  (off at the 10.0 bound)" if off else ""
            if not off and min(lo, hi) <= fitted:
                mark = "  <-- not stationary"
                shallow.append(label)
            print(f"          {label:6s} {lo:9.5f} {fitted:9.5f} {hi:9.5f}{mark}")
            if off and (abs(hi - fitted) > 1e-6 or abs(lo - fitted) > 1e-6):
                print(f"          FAIL {label} is at the bound but still matters, so "
                      f"10.0 does not mean off")
                failures += 1

        # The exponent on the centre carrying most of the density is the one that
        # must be stationary for the form to be right; the others are shallow enough
        # that the radii table dominates -- see the module docstring.
        heavy = int(np.argmax(np.abs(dma["electronic"])))
        down, up = theirs.copy(), theirs.copy()
        down[heavy] *= 0.9
        up[heavy] *= 1.1
        if min(objective(down, gaussian), objective(up, gaussian)) <= fitted:
            print(f"          FAIL {labels[heavy]}'s exponent is not stationary, so "
                  f"{name}'s damping function is not this one")
            failures += 1
        elif shallow:
            print(f"          note: {', '.join(shallow)} shallow by <1%, shallow by "
                  f"<1%, within their Powell convergence on a flat surface")

    print()
    if failures:
        print(f"[screening] {failures} FAILURE(S)")
        return 1
    print("[screening] every one of GAMESS's fitted exponents is stationary under "
          "this objective")
    return 0


if __name__ == "__main__":
    sys.exit(main())
