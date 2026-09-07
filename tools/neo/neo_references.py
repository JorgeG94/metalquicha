#!/usr/bin/env python3
"""Reference energies for the NEO validation cases, from PySCF-NEO.

Run inside a venv that has Yang Yang's PySCF fork installed
(https://github.com/theorychemyang/pyscf; this file was last run at commit
f9c0266). The mqc venv's PySCF has no `neo` module, which is why these
numbers are pasted into `tools/cpu_validation/gen_cpu_validation.py` as
NEO_CASES rather than computed there. Geometries are the generator's.

    python tools/neo/neo_references.py
"""
from pyscf import neo

CASES = [
    # label, atoms (Angstrom), electronic basis, quantum nuclei, xc, epc, grid level
    ("h2_hf", "H 0 0 0; H 0 0 0.7414", "ccpvdz", [0, 1], None, None, None),
    ("hcn_hf", "H 0 0 0; C 0 0 1.064; N 0 0 2.220", "ccpvdz", [0], None, None, None),
    ("hcn_b3lyp5", "H 0 0 0; C 0 0 1.064; N 0 0 2.220", "ccpvdz", [0], "HYB_GGA_XC_B3LYP5", None, 3),
    ("hcn_b3lyp5_epc17-2", "H 0 0 0; C 0 0 1.064; N 0 0 2.220", "ccpvdz", [0], "HYB_GGA_XC_B3LYP5", "17-2", 3),
]

for label, atoms, basis, quantum, xc, epc, level in CASES:
    mol = neo.Mole()
    mol.build(atom=atoms, basis=basis, quantum_nuc=quantum, nuc_basis="pb4d", verbose=0)
    if xc is None:
        mf = neo.HF(mol)
    else:
        mf = neo.KS(mol, xc=xc, epc=epc)
        mf.components["e"].grids.level = level
    mf.conv_tol = 1e-12
    energy = mf.scf()
    print(f"{label:22s} {energy:.12f}")
