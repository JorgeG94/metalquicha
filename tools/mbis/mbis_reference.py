#!/usr/bin/env python3
"""MBIS charges from denspart, for test_mqc_mbis.

Minimal Basis Iterative Stockholder (Verstraelen et al., J. Chem. Theory
Comput. 12, 3894 (2016)): each atom's pro-density is a sum of exponential
shells, N_i a_i^3 / (8 pi) exp(-a_i r), one shell per row of the periodic
table it reaches, and their populations and exponents are fitted to the
molecule's own density.

The reference is denspart (github.com/theochem/denspart), the MBIS code of the
method's author, used as an external program only -- it is LGPL and nothing of
it is copied here. It reaches the answer a different way from mqc: it
minimises the extended Kullback-Leibler divergence with SciPy's trust-constr,
where mqc runs the MBIS fixed-point iteration, and the two meet at the same
stationary point. It also starts from its own table of atomic parameters, where
mqc fills the shells in aufbau order with Slater-screened exponents.

The start matters. The fixed point has more than one stationary point when a
shell is nearly empty: from an even split with exponents 2Z down to 2, LiF's
lithium collapsed its valence shell onto the core and came out at +1.000,
against the +0.940 both codes now reach. In Li+(H2O) the lithium charge is
flat to about 5e-3 across stationary points -- denspart +0.99984, mqc
+1.00407, and mqc from denspart's start +0.99909 -- so it is not a reference
case: there the method, not either code, leaves the number undetermined.

The molecular density is PySCF's, on PySCF's own 150x974 grid, with the orbital
basis read from this repository's JSON. Spin populations are not something
denspart computes; they are the converged pro-model's weights applied to the
spin density here.

    pip install git+https://github.com/theochem/grid.git   # PyPI qc-grid predates SciPy 1.15
    pip install git+https://github.com/theochem/denspart.git
    python3 tools/mbis/mbis_reference.py
"""
import pathlib
import sys

import numpy as np
from denspart.mbis import MBISProModel
from denspart.vh import optimize_reduce_pro_model
from grid.basegrid import Grid
from pyscf import dft, gto, scf

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "spherical_atom"))
from spherical_atom_reference import bse_to_pyscf  # noqa: E402

MOLECULAR_GRID = (150, 974)

# Bohr: the water of test_mqc_fukui and test_mqc_hirshfeld.
WATER = [("O", (0.0, 0.0, 0.0)),
         ("H", (0.0, 1.4308, 1.1078)),
         ("H", (0.0, -1.4308, 1.1078))]
# Bohr: lithium fluoride near its equilibrium bond, an ionic molecule whose
# lithium keeps a valence shell holding under a tenth of an electron.
LIF = [("Li", (0.0, 0.0, 0.0)),
       ("F", (0.0, 0.0, 2.955))]
# Bohr: hydrogen sulfide, for a three-shell atom.
H2S = [("S", (0.0, 0.0, 0.0)),
       ("H", (0.0, 1.8187, 1.7419)),
       ("H", (0.0, -1.8187, 1.7419))]


def molecule(atoms, basis, charge=0, spin=0):
    symbols = sorted({a[0] for a in atoms})
    return gto.M(atom=atoms, unit="Bohr", charge=charge, spin=spin,
                 basis={s: bse_to_pyscf(basis, gto.charge(s)) for s in symbols},
                 cart=False, verbose=0)


def converge(mol):
    mf = scf.RHF(mol) if mol.spin == 0 else scf.UHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.run()
    assert mf.converged
    dm = np.asarray(mf.make_rdm1())
    if dm.ndim == 2:
        return dm, None
    return dm[0] + dm[1], dm[0] - dm[1]


def mbis(mol, total, spin=None):
    grids = dft.Grids(mol)
    grids.atom_grid = MOLECULAR_GRID
    grids.build()
    ni = dft.numint.NumInt()
    ao = ni.eval_ao(mol, grids.coords)
    rho = ni.eval_rho(mol, ao, total)
    grid = Grid(grids.coords, grids.weights)
    model0 = MBISProModel.from_geometry(mol.atom_charges(), mol.atom_coords())
    model, _ = optimize_reduce_pro_model(model0, grid, rho, gtol=1e-12, maxiter=5000,
                                         density_cutoff=0.0)
    spins = None
    if spin is not None:
        spin_rho = ni.eval_rho(mol, ao, spin)
        pro = np.array([model.compute_proatom(a, grids.coords) for a in range(mol.natm)])
        spins = (pro / pro.sum(axis=0)) @ (grids.weights * spin_rho)
    return model, spins


def show(label, values):
    print("  %-10s" % label + " ".join("%17.12f" % v for v in values))


def report(name, model, spins=None, total=0.0):
    print(name)
    show("charges", model.charges)
    if spins is not None:
        show("spins", spins)
    print("  sum        %.3e" % (model.charges.sum() - total))
    for fn in model.fns:
        print("  atom %d  population %15.10f  exponent %15.10f"
              % (fn.iatom + 1, fn.population, fn.exponent))


def main():
    neutral = molecule(WATER, "cc-pvdz")
    model, _ = mbis(neutral, converge(neutral)[0])
    report("water_ccpvdz RHF", model)

    cation = molecule(WATER, "cc-pvdz", charge=1, spin=1)
    model, spins = mbis(cation, *converge(cation))
    report("water_ccpvdz_cation UHF", model, spins, 1.0)

    h2s = molecule(H2S, "cc-pvdz")
    model, _ = mbis(h2s, converge(h2s)[0])
    report("h2s_ccpvdz RHF", model)

    lif = molecule(LIF, "cc-pvdz")
    model, _ = mbis(lif, converge(lif)[0])
    report("lif_ccpvdz RHF", model)


if __name__ == "__main__":
    main()
