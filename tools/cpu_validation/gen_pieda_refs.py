"""Independent FMO2 (exact ESP) + PIEDA reference numbers, in PySCF.

Ported from the pieda_l3 worktree's scratch script (`pieda_scratch/fmo2_pieda_ref.py`)
onto this repository's own basis reader (`bse_to_pyscf`, `gen_cpu_validation.py`),
so a basis name here means exactly what it means to `mqc` -- see the module
docstring of `gen_cpu_validation.py` for why that matters to the eighth decimal
on a Pople set.

Terms, per pair IJ, all Hartree, `E'` an internal energy (`E - Tr(D u)`):

    dE_IJ = E'_IJ - E'_I - E'_J + Tr[(D_IJ - D_I (+) D_J) u_IJ]   (FMO2 pair energy)
    Ees   = Tr[D_I V^J_nuc] + Tr[D_J V^I_nuc] + sum D^I D^J (mu nu|la si) + Enn(I,J)
    E'^HL = E[D_HL; h + u_IJ] - Tr[D_HL u_IJ],   D_HL = 2 C (C^T S C)^-1 C^T,
            C = [C_I, C_J]
    Eex   = E'^HL - E'_I - E'_J - Ees
    Ect+mix = dE_IJ - Ees - Eex

`C_I` is recovered from the converged monomer density alone, by a pivoted
Cholesky of `D_I / 2` -- exactly what `mqc_czt_pieda.f90`'s
`cholesky_occupied_orbitals` does with LAPACK `dpstrf`.

Usage, from the repository root (a configured build tree finds the basis
JSON automatically through `gen_cpu_validation.basis_file`):

    python3 tools/cpu_validation/gen_pieda_refs.py water3_cyclic.xyz 6-31g w3 '[[0,1,2],[3,4,5],[6,7,8]]'
    python3 tools/cpu_validation/gen_pieda_refs.py gly3_water_pair.xyz 6-31g glyw \\
        '[[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23],[24,25,26]]'

The xyz path is resolved against `validation/inputs/sample_inputs/` first, and
against the current directory otherwise -- both `run_validation.py` and this
script's own geometries can be named without a path.
"""
from __future__ import annotations

import itertools
import json
import pathlib
import sys

import numpy as np
from scipy.linalg import lapack

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from gen_cpu_validation import XYZ_DIR, bse_to_pyscf  # noqa: E402

from pyscf import gto, lib, scf  # noqa: E402

ANG = 1.0 / lib.param.BOHR


def read_xyz(path):
    candidate = XYZ_DIR / path
    if not candidate.exists():
        candidate = pathlib.Path(path)
    lines = candidate.read_text().splitlines()
    n = int(lines[0])
    atoms = []
    for line in lines[2:2 + n]:
        parts = line.split()
        atoms.append((parts[0], [float(x) for x in parts[1:4]]))
    return atoms


def make_mol(atoms, basis):
    symbols = sorted({a[0] for a in atoms})
    return gto.M(
        atom=[(s, np.array(x) * ANG) for s, x in atoms],
        unit="Bohr",
        basis={s: bse_to_pyscf(basis, s) for s in symbols},
        cart=False,
        verbose=0,
    )


def rhf(mol, hcore, dm0=None):
    """A closed-shell SCF pinned to `hcore`, converged tight enough that its
    noise never competes with the ~1e-9 gate this reference is checked at."""
    mf = scf.RHF(mol)
    mf.conv_tol = 1.0e-13
    mf.conv_tol_grad = 1.0e-9
    mf.max_cycle = 300
    mf.get_hcore = lambda *a: hcore
    mf.kernel(dm0=dm0)
    if not mf.converged:
        raise SystemExit("monomer or pair SCF did not converge")
    return mf


def cholesky_occupied(d):
    """`C_occ` from a converged closed-shell density, exactly as
    `mqc_czt_pieda.cholesky_occupied_orbitals` recovers it: a pivoted
    Cholesky of `D/2` (LAPACK `dpstrf`), no diagonalisation."""
    a = 0.5 * d.copy()
    n = a.shape[0]
    factor, piv, rank, _ = lapack.dpstrf(a, lower=1, tol=-1.0)
    l_factor = np.tril(factor)[:, :rank]
    c = np.zeros((n, rank))
    c[piv - 1, :] = l_factor
    return c, rank


def main(xyz, basis, tag, frags):
    atoms = read_xyz(xyz)
    full = make_mol(atoms, basis)
    ao_slices = full.aoslice_by_atom()

    def ao_of(ats):
        return np.concatenate([np.arange(ao_slices[a, 2], ao_slices[a, 3]) for a in ats])

    nao = full.nao
    kinetic = full.intor("int1e_kin")
    v_atom = []
    for a in range(full.natm):
        with full.with_rinv_origin(full.atom_coord(a)):
            v_atom.append(-full.atom_charge(a) * full.intor("int1e_rinv"))
    v_atom = np.array(v_atom)
    charges = full.atom_charges()
    coords = full.atom_coords()
    n_frag = len(frags)
    idx = [ao_of(f) for f in frags]
    fmols = [make_mol([atoms[a] for a in f], basis) for f in frags]
    density = [None] * n_frag
    e_internal = [0.0] * n_frag

    def coulomb_of(dms):
        dm = np.zeros((nao, nao))
        for k, d in dms:
            dm[np.ix_(idx[k], idx[k])] += d
        return scf.hf.get_jk(full, dm, hermi=1, with_k=False)[0]

    def field_on(group):
        outside = [k for k in range(n_frag) if k not in group]
        u = sum((v_atom[a] for k in outside for a in frags[k]), np.zeros((nao, nao)))
        dms = [(k, density[k]) for k in outside if density[k] is not None]
        if dms:
            u = u + coulomb_of(dms)
        return u

    def own_h(group):
        atoms_here = [a for k in group for a in frags[k]]
        return kinetic + sum(v_atom[a] for a in atoms_here)

    # Outer (monomer) self-consistency: the same loop `mqc_czt_fmo`'s
    # `calculate_monomers` runs, to the tolerance the deck asks for.
    previous = None
    for outer in range(200):
        updated = []
        for i in range(n_frag):
            ii = idx[i]
            u = field_on([i]) if outer > 0 else np.zeros((nao, nao))
            h = (own_h([i]) + u)[np.ix_(ii, ii)]
            mf = rhf(fmols[i], h, dm0=density[i])
            updated.append(mf.make_rdm1())
            e_internal[i] = mf.e_tot - np.sum(mf.make_rdm1() * u[np.ix_(ii, ii)])
        density = updated
        total = sum(e_internal)
        if previous is not None and abs(total - previous) < 1.0e-11:
            break
        previous = total
    print(f"[{tag}] outer iterations {outer}, sum E'_I = {sum(e_internal):.12f}")
    for i in range(n_frag):
        d = density[i]
        s_block = full.intor("int1e_ovlp")[np.ix_(idx[i], idx[i])]
        idem = np.abs(d @ s_block @ d - 2 * d).max()
        _, rank = cholesky_occupied(d)
        n_occ = round(np.trace(d @ s_block) / 2)
        print(f"[{tag}] monomer {i + 1}: E'={e_internal[i]:.12f} "
              f"|DSD-2D|={idem:.2e} rank={rank} n_occ={n_occ}")

    total = sum(e_internal)
    rows = []
    for i, j in itertools.combinations(range(n_frag), 2):
        group = [i, j]
        ij = np.concatenate([idx[i], idx[j]])
        u = field_on(group)[np.ix_(ij, ij)]
        h = own_h(group)[np.ix_(ij, ij)]
        dimer_mol = make_mol([atoms[a] for a in frags[i] + frags[j]], basis)
        n_i = len(idx[i])
        d_split = np.zeros((len(ij), len(ij)))
        d_split[:n_i, :n_i] = density[i]
        d_split[n_i:, n_i:] = density[j]
        mf = rhf(dimer_mol, h + u, dm0=d_split)
        d_ij = mf.make_rdm1()
        e_pair_internal = mf.e_tot - np.sum(d_ij * u)
        response = np.sum((d_ij - d_split) * u)
        delta_e = e_pair_internal - e_internal[i] - e_internal[j] + response
        total += delta_e

        d_i = np.zeros_like(d_split)
        d_i[:n_i, :n_i] = density[i]
        d_j = np.zeros_like(d_split)
        d_j[n_i:, n_i:] = density[j]
        v_j_on_i = sum(v_atom[a] for a in frags[i])[np.ix_(ij, ij)]
        v_i_on_j = sum(v_atom[a] for a in frags[j])[np.ix_(ij, ij)]
        j_of_j = scf.hf.get_jk(dimer_mol, d_j, hermi=1, with_k=False)[0]
        e_nn = sum(charges[a] * charges[b] / np.linalg.norm(coords[a] - coords[b])
                   for a in frags[i] for b in frags[j])
        ees = np.sum(d_i * v_i_on_j) + np.sum(d_j * v_j_on_i) + np.sum(d_i * j_of_j) + e_nn

        s_dimer = dimer_mol.intor("int1e_ovlp")
        c_i, rank_i = cholesky_occupied(density[i])
        c_j, rank_j = cholesky_occupied(density[j])
        c_union = np.zeros((len(ij), c_i.shape[1] + c_j.shape[1]))
        c_union[:n_i, :c_i.shape[1]] = c_i
        c_union[n_i:, c_i.shape[1]:] = c_j
        m_mat = c_union.T @ s_dimer @ c_union
        d_hl = 2 * c_union @ np.linalg.solve(m_mat, c_union.T)
        v_coulomb, v_exchange = scf.hf.get_jk(dimer_mol, d_hl, hermi=1)
        e_hl_total = np.sum(d_hl * (h + u)) + 0.5 * np.sum(d_hl * (v_coulomb - 0.5 * v_exchange))
        e_hl_total += dimer_mol.energy_nuc()
        e_hl_prime = e_hl_total - np.sum(d_hl * u)

        eex = e_hl_prime - e_internal[i] - e_internal[j] - ees
        ect_mix = delta_e - ees - eex
        rows.append((i + 1, j + 1, delta_e, ees, eex, ect_mix))
        print(f"[{tag}] pair {i + 1}-{j + 1}: dE={delta_e:.12f} Ees={ees:.12f} "
              f"Eex={eex:.12f} Ect+mix={ect_mix:.12f} "
              f"rank(C_I)={rank_i} rank(C_J)={rank_j} minEig(CtSC)={np.linalg.eigvalsh(m_mat).min():.4f}")
    print(f"[{tag}] FMO2 total = {total:.12f}")
    supermolecule = scf.RHF(full)
    supermolecule.conv_tol = 1.0e-12
    supermolecule.verbose = 0
    print(f"[{tag}] supermolecular RHF = {supermolecule.kernel():.12f}")
    return rows


if __name__ == "__main__":
    xyz_arg, basis_arg, tag_arg = sys.argv[1], sys.argv[2], sys.argv[3]
    frags_arg = json.loads(sys.argv[4])
    main(xyz_arg, basis_arg, tag_arg, frags_arg)
