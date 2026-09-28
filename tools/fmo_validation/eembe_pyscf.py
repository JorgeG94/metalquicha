#!/usr/bin/env python3
"""Independent PySCF reference for metalquicha's EE-MBE energy.

Reproduces the electrostatically-embedded many-body expansion (EE-MBE) that
``keywords.fragmentation.method: "ee-mbe"`` selects, at level 2, Hartree-Fock
and MP2, RI-MP2 and SCS-MP2 on top of the same embedded Hartree-Fock
reference.

Run from anywhere, e.g.::

    <python> tools/fmo_validation/eembe_pyscf.py --method hf
    <python> tools/fmo_validation/eembe_pyscf.py --method mp2 --all-electron
    <python> tools/fmo_validation/eembe_pyscf.py --method ri-mp2
    <python> tools/fmo_validation/eembe_pyscf.py --method scs-mp2

============================================================================
Semantics, established by reading the code (not by guessing)
============================================================================

``keywords.fragmentation.method: "ee-mbe"`` is read in ``mqc_driver.f90``
(``run_fragmented_calculation``, around the ``config%expansion_kind ==
"ee-mbe"`` branch). It builds an ``fmo_context_t`` with ``expansion = "mbe"``
and ``esp = "ptc"`` (point-charge embedding) -- the same machinery FMO uses
(``expansion = "fmo"``, ``esp = "exact"``), differing only in the field a
fragment sees and how the pieces are summed. The whole run happens in
``backends/cenzontle/fragments/mqc_czt_fmo.f90``.

**The outer (self-consistent charge, "SCC") loop** -- ``run_fmo2`` ->
``calculate_monomers``:

  1. Every monomer is solved once *in vacuo* (``solve_fragment(..., bare=.true.)``)
     to get a starting set of atomic charges.
  2. Then, for each outer pass: ``all_charges`` sums each atom's Mulliken
     charge from whichever fragment currently "owns" it (every atom but a
     detached one belongs to exactly one fragment; EE-MBE on this system cuts
     no bonds, so it is a plain per-atom charge). Every fragment is then
     re-solved **against that one fixed set of charges** -- a synchronous
     (Jacobi) update, not Gauss-Seidel: within a pass every fragment reads the
     *previous* pass's charges, never a value another fragment already
     produced this same pass (``calculate_monomers``'s own comment: solving
     into a copy and swapping afterwards is deliberate, "2lty's first pass
     came out 0.44 Hartree apart on one rank and on four" when it wasn't).
  3. Convergence is tested on the sum of every fragment's *internal* energy
     (``frag(:)%energy``, defined below) between passes:
     ``res%outer_change = abs(e_sum - e_prev)``, converged when
     ``res%outer_change < opts%outer_tol``. The default is
     ``keywords.fragmentation.outer_tol`` -> ``fmo_tolerance`` = ``1.0e-7``
     Hartree (`mqc_config_types.f90`), and the cap is
     ``fmo_max_outer`` = 50 passes. The water-trimer deck used here overrides
     neither.

  Charges: ``far_field`` (``keywords.fragmentation.far_field``, default
  ``"mulliken"`` -- `mqc_config_adapter.f90`) selects Mulliken over the
  converged fragment density, via ``fragment_charges`` ->
  ``mulliken_charges`` (`backends/cenzontle/analysis/mqc_czt_charges.f90`):
  ``q_A = Z_A - sum_{mu in A} (D S)_mu,mu``, with ``D`` the *total* (both
  spins) density.

**The embedding operator** -- ``embedding_operator`` in ``mqc_czt_fmo.f90``.
``effective_resppc`` forces the near/far cutoff to exactly 0 whenever
``esp == "ptc"``, so ``near_fragments`` (cutoff ``best <= resppc`` with
``resppc = 0``) never returns a fragment as "near": under EE-MBE *every*
fragment outside the group being solved is treated as a set of far-field
point charges, none through the exact two-electron (Coulomb) term. For a
group ``G`` (a monomer or an n-mer) and every atom ``A`` outside ``G``:

  ``u = -sum_A q_A * <mu| 1/|r-R_A| |nu>``

(``embedding_operator``'s own comment: "An electron carries charge -1, so a
positive charge lowers its energy: the operator is
``-sum_g w_g/|r - R_g|``"; ``esp_matrices``, in
`backends/cenzontle/integrals/mqc_czt_esp.f90`, returns the bare
``<mu|1/|r-R_g||nu>`` matrix with no charge folded in). This is added to the
one-electron Hamiltonian (``request%h_extra``) before the SCF runs
(`inner_scf`, `nmer_term`). No cut bonds exist on this system, so
``group_own_charge`` (``own_q``) and ``shared`` are always zero/false, and
``near`` is always empty -- both the AFO-specific and the "near fragment"
machinery are inert here; this script implements only the point-charge case.

**No nucleus-point-charge interaction energy is added anywhere.** ``u`` only
ever enters the *electronic* Hamiltonian; the SCF's own nuclear-repulsion
term (``mol.energy_nuc()`` in PySCF terms) covers only the group's own
nuclei. This differs from PySCF's ``pyscf.qmmm.mm_charge`` decorator, whose
docstring says plainly it adds "the interaction between the nuclei in QM
region and the MM charges" on top of the Hcore term -- so that decorator is
*not* used here; this script reproduces only the Hcore modification
(``QMMMSCF.get_hcore`` in ``pyscf/qmmm/itrf.py``, which is exactly
``h1e += einsum('kpq,k->pq', int1e_grids, -charges)``, the same sign
convention as ``u`` above) and leaves ``energy_nuc()`` untouched.

**Which energies enter the EE-MBE total** -- ``inner_scf`` /
``solve_fragment_method``
(`backends/cenzontle/fragments/mqc_czt_fragment_solver.f90`):
``outcome%reference`` is the converged SCF energy *including* ``Tr(D u)``
(because ``u`` was added straight into Hcore, this term falls out of the
ordinary SCF energy expression with no separate bookkeeping);
``outcome%energy = outcome%reference + outcome%correlation``;
``outcome%internal = outcome%energy - Tr(D_scf u)`` when there is a field.
``nmer_term`` picks, for ``opts%expansion == "mbe"`` (EE-MBE):
``e_internal = outcome%energy`` (the *full*, embedded, correlated energy) and
``e_resp = 0`` -- EE-MBE, unlike FMO, does **not** subtract a response term
``Tr((D-D_split) u)``; that response bookkeeping is FMO-only.  Likewise
``run_fmo2`` sets each monomer's contribution to ``frag(i)%energy_total``
(``outcome%energy``, the full embedded value) rather than ``frag(i)%energy``
(the internal one) exactly when ``opts%expansion == "mbe"``.

**Assembling the total** -- ``calculate_polymers`` + ``subtract_subsets``
(`backends/cenzontle/fragments/mqc_czt_subsets.f90`). Every subset of the
fragment list up to the truncation level is enumerated; each term's
"correction" starts as the group's own value (``energy_total`` for size 1,
``e_internal + e_resp`` for size >= 2) and is then reduced by every proper
subset's already-final correction:

  ``dE_S = f(S) - sum_{T subset S, T != S, T != {}} dE_T``

At level 2 this is exactly ``dE_{IJ} = E'_{IJ} - E'_I - E'_J``, and

  ``E_total = sum_I E'_I + sum_{pairs} dE_{IJ}``

which is what ``res%energy = res%monomer_sum + res%pair_sum`` computes. No
resdim/"separated pair" electrostatics-only shortcut applies here:
``keywords.fragmentation.resdim`` is refused outright for EE-MBE
(``mqc_config_adapter.f90``, ``check ... resdim``) and forced to 0 by
``mqc_driver.f90``, and ``separated_pairs`` returns "none separated" whenever
``opts%resdim <= 0``.

**Correlation (MP2 / RI-MP2 / SCS-MP2)** -- as
``mqc_czt_fragment_solver::solve_fragment_method`` / ``fragment_correlation``
add it:

  * The SCC loop (embedding charges) is run in Hartree-Fock only -- there is
    no MP2 density feeding the charges anywhere in this code, embedded or
    not.
  * Correlation is added *after* the embedded HF converges, on that SCF's own
    canonical orbitals and orbital energies (``fragment_correlation`` takes
    ``scf`` -- the already-converged ``rhf_result_t`` -- and calls
    ``run_czt_mp2``/``run_czt_ri_mp2`` on it): i.e. plain MP2 in the embedded
    HF molecular-orbital basis, exactly as ``pyscf.mp.MP2(mf)`` would give
    for the same ``mf``.
  * ``outcome%energy = outcome%reference + outcome%correlation`` --
    correlation is simply added to the embedded HF energy (``reference``,
    which already includes ``Tr(D_HF u)``); it does not re-enter the field or
    get its own ``Tr(D u)`` term. Combined with the ``expansion == "mbe"``
    rule above, a correlated monomer/dimer's contribution to the EE-MBE total
    is exactly ``E'_X = E_HF,embedded(X) + E_corr(X)``, and the total is
    assembled from those ``E'_X`` by the same
    ``sum_I E'_I + sum_pairs(E'_{IJ}-E'_I-E'_J)`` formula as Hartree-Fock.
  * Frozen core: ``method%n_frozen_core`` defaults to -1, meaning
    ``core_orbital_count(real_z)`` (`src/core/mqc_elements.f90`) -- the
    number of *filled shells below the valence one*, summed over the group's
    own real (non-ghost) atoms: 0 for H/He, 1 (the 1s) for Li..Ne, and so on.
    For an all-oxygen-and-hydrogen system that is 1 orbital per oxygen and 0
    per hydrogen -- i.e. each water's own 1s(O) -- matching a dimer's 2 and a
    monomer's 1. ``method%freeze_core`` defaults to ``.true.``
    (`mqc_cuest_iface.f90`).
  * RI-MP2's auxiliary basis is ``model.aux_basis`` when the deck names one,
    else ``cuest_scf_settings_t%aux_basis_set``'s own default,
    ``"def2-universal-jkfit"`` (`mqc_cuest_iface.f90`) -- the *same* default
    used for SCF density fitting; ``correlation_aux_basis`` in
    ``mqc_czt_bridge.f90`` warns (but does not refuse) that a JKFIT set is
    not a correlation-fitting (RIFIT) set, and uses it anyway when nothing
    else was named. That is the default this script uses too, read out of
    this repository's own ``basis_sets/def2-universal-jkfit.json`` file
    (never PySCF's bundled table -- see ``bse_to_pyscf`` below).
  * SCS-MP2's spin-component scale factors are
    ``cuest_scf_settings_t%scs_os = 1.2``, ``%scs_ss = 1/3``
    (`mqc_method_config.f90`); plain MP2 and RI-MP2 use 1.0/1.0.

**Basis form.** 6-31G has no shell above p, so the spherical-vs-Cartesian
question (`src/basis/mqc_json_basis_reader.f90`) never actually bites for
this deck's basis; every group is built the same way this repository's own
``tools/cpu_validation/gen_cpu_validation.py`` builds any PySCF reference,
through ``bse_to_pyscf``/``molecule_form``, which is imported (not
reimplemented) from that file below to guarantee the identical basis-set
JSON is read the identical way -- PySCF's own internal Pople tables round
Pople-basis coefficients differently in the eighth decimal, which would look
exactly like a bug in this script and is not one.

============================================================================
What is, and is not, checked against the binary
============================================================================

Only Hartree-Fock is checked against a copy of ``build/mqc`` run on
``validation/inputs/cpu/mqc/fmo/eembe_water3.json`` (the ``model.method:
"hf"``, ``keywords.fragmentation: {method: "ee-mbe", level: 2}`` deck this
task named), because MP2 is not yet wired into FMO/EE-MBE in the binary (see
above) -- there is nothing to run it against yet. The MP2/RI-MP2/SCS-MP2
numbers below are produced from the semantics above alone, so that a test
written once the Fortran side lands has a number to check against that was
derived independently of it.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools" / "cpu_validation"))
import gen_cpu_validation as gcv  # noqa: E402  (repo basis-JSON loader; read-only import)

from pyscf import gto, mp, scf  # noqa: E402

# ---------------------------------------------------------------------------
# The deck this script reproduces:
# validation/inputs/cpu/mqc/fmo/eembe_water3.json
# ---------------------------------------------------------------------------
GEOM_FILE = REPO / "validation" / "inputs" / "sample_inputs" / "w3.xyz"
BASIS = "6-31g"
FRAGMENTS = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]  # 0-based, file order (per the deck)
PAIRS = [(0, 1), (1, 2), (0, 2)]

# fmo_tolerance / fmo_max_outer defaults, mqc_config_types.f90
OUTER_TOL = 1.0e-7
MAX_OUTER = 50

# SCF is converged far tighter than mqc's own inner tolerances
# (fmo_scf_energy_tol=1e-9, fmo_scf_density_tol=1e-7) so this script's answer
# is the tight-convergence limit mqc's looser inner tolerance is aiming at.
SCF_CONV_TOL = 1.0e-12
SCF_MAX_CYCLE = 200

# cuest_scf_settings_t defaults, mqc_cuest_iface.f90 / mqc_method_config.f90
AUX_BASIS_DEFAULT = "def2-universal-jkfit"
SCS_OS = 1.2
SCS_SS = 1.0 / 3.0


def read_xyz(path: Path):
    lines = path.read_text().splitlines()
    n = int(lines[0].strip())
    atoms = []
    for line in lines[2:2 + n]:
        parts = line.split()
        sym = parts[0]
        x, y, z = (float(v) for v in parts[1:4])
        atoms.append((sym, x, y, z))
    if len(atoms) != n:
        raise ValueError(f"{path}: header says {n} atoms, found {len(atoms)}")
    return atoms


def core_orbital_count(atomic_numbers) -> int:
    """Mirror src/core/mqc_elements.f90::core_orbital_count exactly."""
    n_core = 0
    for z in atomic_numbers:
        if z <= 10:
            if z > 2:
                n_core += 1
        elif z <= 18:
            n_core += 5
        elif z <= 36:
            n_core += 9
        elif z <= 54:
            n_core += 18
        else:
            n_core += 27
    return n_core


def build_mol(atoms, indices, basis, verbose=0):
    """One group's Mole, atoms end to end in `indices` order (assemble_group)."""
    symbols = [atoms[i][0] for i in indices]
    coords = [atoms[i][1:4] for i in indices]
    mol = gto.Mole()
    mol.atom = list(zip(symbols, coords))
    mol.unit = "Angstrom"
    uniq = sorted(set(symbols))
    mol.basis = {s: gcv.bse_to_pyscf(basis, s) for s in uniq}
    mol.cart = gcv.molecule_form(basis, uniq) == gcv.CARTESIAN
    mol.charge = 0
    mol.spin = 0
    mol.verbose = verbose
    mol.build()
    return mol


def mulliken_charges(mol, dm) -> np.ndarray:
    """q_A = Z_A - sum_{mu in A} (D S)_mu,mu -- mqc_czt_charges.f90::mulliken_charges."""
    s = mol.intor("int1e_ovlp")
    pop_ao = np.einsum("pq,qp->p", dm, s)
    charges = np.array(mol.atom_charges(), dtype=float)
    for ia in range(mol.natm):
        p0, p1 = mol.aoslice_by_atom()[ia, 2:4]
        charges[ia] -= pop_ao[p0:p1].sum()
    return charges


def embedding_h1(mol, ext_coords_bohr: np.ndarray, ext_charges: np.ndarray):
    """u = -sum_A q_A <mu|1/|r-R_A||nu> -- mqc_czt_fmo.f90::embedding_operator.

    Same sign convention pyscf.qmmm's own get_hcore uses (see module
    docstring); unlike that decorator, no nucleus-charge term is added here.
    Returns None (no field) when there is nothing outside the group.
    """
    if ext_charges.size == 0:
        return None
    j3c = mol.intor("int1e_grids", hermi=1, grids=ext_coords_bohr)
    return -np.einsum("kpq,k->pq", j3c, ext_charges)


def run_embedded_rhf(mol, u: np.ndarray | None, dm0=None):
    """One embedded RHF, converged tight. u enters Hcore only; energy_nuc is untouched."""
    mf = scf.RHF(mol)
    mf.conv_tol = SCF_CONV_TOL
    mf.conv_tol_grad = 1.0e-10
    mf.max_cycle = SCF_MAX_CYCLE
    if u is not None:
        h1e0 = mf.get_hcore(mol)

        def get_hcore(molx=None, _h1e0=h1e0, _u=u):
            return _h1e0 + _u

        mf.get_hcore = get_hcore
    e_total = mf.kernel(dm0) if dm0 is not None else mf.kernel()
    if not mf.converged:
        raise RuntimeError("embedded RHF did not converge")
    dm = mf.make_rdm1()
    if u is not None:
        e_internal = e_total - np.einsum("pq,pq->", dm, u)
    else:
        e_internal = e_total
    return mf, e_total, e_internal, dm


def outside_indices(n_atoms: int, inside):
    inside_set = set(inside)
    return [a for a in range(n_atoms) if a not in inside_set]


def run_hf_scc(atoms, basis=BASIS, verbose=False):
    """The outer Jacobi SCC loop, then the three dimers. Returns everything downstream needs."""
    n_atoms = len(atoms)
    full_mol = build_mol(atoms, range(n_atoms), basis)
    coords_bohr = full_mol.atom_coords()
    frag_mols = [build_mol(atoms, idx, basis) for idx in FRAGMENTS]

    # Pass 0: every monomer in vacuo (solve_fragment(..., bare=.true.)).
    charges = np.zeros(n_atoms)
    frag = []
    for idx, mol in zip(FRAGMENTS, frag_mols):
        mf, e_tot, e_int, dm = run_embedded_rhf(mol, u=None)
        q = mulliken_charges(mol, dm)
        charges[idx] = q
        frag.append(dict(mf=mf, dm=dm, e_total=e_tot, e_internal=e_int))
    e_prev = sum(f["e_internal"] for f in frag)

    outer_trace = []
    for outer in range(1, MAX_OUTER + 1):
        q_all = charges.copy()  # the pass's fixed field: last pass's charges
        new_frag = []
        new_charges = charges.copy()
        for idx, mol in zip(FRAGMENTS, frag_mols):
            out = outside_indices(n_atoms, idx)
            u = embedding_h1(mol, coords_bohr[out], q_all[out])
            mf, e_tot, e_int, dm = run_embedded_rhf(mol, u, dm0=None)
            q = mulliken_charges(mol, dm)
            new_charges[idx] = q
            new_frag.append(dict(mf=mf, dm=dm, e_total=e_tot, e_internal=e_int))
        e_sum = sum(f["e_internal"] for f in new_frag)
        change = abs(e_sum - e_prev)
        outer_trace.append((outer, e_sum, change))
        if verbose:
            print(f"  outer {outer:2d}  internal sum {e_sum:20.12f}  change {change:.3e}")
        frag, charges = new_frag, new_charges
        if change < OUTER_TOL:
            break
        e_prev = e_sum
    else:
        raise RuntimeError(f"SCC loop did not converge in {MAX_OUTER} passes")

    monomer_hf_total = [f["e_total"] for f in frag]

    # Dimers, embedded in the converged charges.
    pair_mols = {}
    pair_hf_total = {}
    pair_mf = {}
    for (i, j) in PAIRS:
        idx = FRAGMENTS[i] + FRAGMENTS[j]
        mol = build_mol(atoms, idx, basis)
        out = outside_indices(n_atoms, idx)
        u = embedding_h1(mol, coords_bohr[out], charges[out])
        mf, e_tot, e_int, dm = run_embedded_rhf(mol, u)
        pair_mols[(i, j)] = mol
        pair_hf_total[(i, j)] = e_tot
        pair_mf[(i, j)] = mf

    monomer_sum = sum(monomer_hf_total)
    pair_delta = {p: pair_hf_total[p] - monomer_hf_total[p[0]] - monomer_hf_total[p[1]]
                  for p in PAIRS}
    pair_sum = sum(pair_delta.values())
    total = monomer_sum + pair_sum

    return dict(
        atoms=atoms, n_atoms=n_atoms, charges=charges,
        frag_mols=frag_mols, frag=frag, monomer_hf_total=monomer_hf_total,
        pair_mols=pair_mols, pair_mf=pair_mf, pair_hf_total=pair_hf_total,
        pair_delta=pair_delta, monomer_sum=monomer_sum, pair_sum=pair_sum,
        total=total, outer_trace=outer_trace,
    )


def correlation_energy(mf, real_z, scs_os, scs_ss, freeze_core, density_fitting, aux_basis):
    """fragment_correlation, mqc_czt_fragment_solver.f90: MP2/RI-MP2 on mf's own orbitals."""
    frozen = core_orbital_count(real_z) if freeze_core else 0
    if density_fitting:
        pt = mp.dfmp2.DFMP2(mf, frozen=frozen)
        symbols = sorted({mf.mol.atom_symbol(i) for i in range(mf.mol.natm)})
        auxbasis = {s: gcv.bse_to_pyscf(aux_basis, s) for s in symbols}
        pt.with_df = pt.with_df.__class__(mf.mol, auxbasis=auxbasis)
        pt.kernel()
        os_e, ss_e = pt.e_corr_os, pt.e_corr_ss
    else:
        pt = mp.MP2(mf, frozen=frozen)
        pt.kernel()
        os_e, ss_e = pt.e_corr_os, pt.e_corr_ss
    return scs_os * os_e + scs_ss * ss_e, frozen


def run_correlated(method, all_electron, basis=BASIS, aux_basis=AUX_BASIS_DEFAULT):
    scc = run_hf_scc(read_xyz(GEOM_FILE), basis)
    atoms = scc["atoms"]
    freeze_core = not all_electron
    density_fitting = method == "ri-mp2"
    if method == "scs-mp2":
        scs_os, scs_ss = SCS_OS, SCS_SS
    else:
        scs_os, scs_ss = 1.0, 1.0

    monomer_total = []
    monomer_corr = []
    for idx, mf in zip(FRAGMENTS, (f["mf"] for f in scc["frag"])):
        real_z = [atoms[i][0] for i in idx]
        real_z = [gto.charge(s) for s in real_z]
        corr, frozen = correlation_energy(mf, real_z, scs_os, scs_ss, freeze_core,
                                          density_fitting, aux_basis)
        monomer_corr.append(corr)
        monomer_total.append(scc["monomer_hf_total"][len(monomer_total)] + corr)

    pair_total = {}
    pair_corr = {}
    for p in PAIRS:
        idx = FRAGMENTS[p[0]] + FRAGMENTS[p[1]]
        real_z = [gto.charge(atoms[i][0]) for i in idx]
        corr, frozen = correlation_energy(scc["pair_mf"][p], real_z, scs_os, scs_ss,
                                          freeze_core, density_fitting, aux_basis)
        pair_corr[p] = corr
        pair_total[p] = scc["pair_hf_total"][p] + corr

    pair_delta = {p: pair_total[p] - monomer_total[p[0]] - monomer_total[p[1]] for p in PAIRS}
    monomer_sum = sum(monomer_total)
    pair_sum = sum(pair_delta.values())
    total = monomer_sum + pair_sum
    return dict(scc=scc, monomer_total=monomer_total, monomer_corr=monomer_corr,
                pair_total=pair_total, pair_corr=pair_corr, pair_delta=pair_delta,
                monomer_sum=monomer_sum, pair_sum=pair_sum, total=total)


def print_hf(scc):
    print("EE-MBE(2)/6-31G water trimer -- Hartree-Fock")
    print("SCC outer iterations:")
    for outer, e_sum, change in scc["outer_trace"]:
        print(f"  {outer:2d}  internal sum {e_sum:20.12f}  change {change:.3e}")
    print()
    for i, e in enumerate(scc["monomer_hf_total"], start=1):
        print(f"  monomer {i}  E'  {e:20.12f}")
    print(f"  monomer_sum      {scc['monomer_sum']:20.12f}")
    for p in PAIRS:
        i, j = p[0] + 1, p[1] + 1
        print(f"  pair {i}-{j}   dE  {scc['pair_delta'][p]:20.12f}")
    print(f"  pair_sum         {scc['pair_sum']:20.12f}")
    print(f"  TOTAL            {scc['total']:20.12f}")


def print_correlated(method, result):
    print(f"EE-MBE(2)/6-31G water trimer -- {method}"
          + (" (all-electron)" if result.get("all_electron") else " (frozen core)"))
    for i, (e, c) in enumerate(zip(result["monomer_total"], result["monomer_corr"]), start=1):
        print(f"  monomer {i}  E' = E_HF+corr  {e:20.12f}   (corr {c:20.12f})")
    print(f"  monomer_sum      {result['monomer_sum']:20.12f}")
    for p in PAIRS:
        i, j = p[0] + 1, p[1] + 1
        print(f"  pair {i}-{j}   dE  {result['pair_delta'][p]:20.12f}"
              f"   (E'_pair {result['pair_total'][p]:20.12f}, corr {result['pair_corr'][p]:20.12f})")
    print(f"  pair_sum         {result['pair_sum']:20.12f}")
    print(f"  TOTAL            {result['total']:20.12f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--method", choices=["hf", "mp2", "ri-mp2", "scs-mp2"], default="hf")
    ap.add_argument("--all-electron", action="store_true",
                    help="MP2 variants: no frozen core (default: freeze core_orbital_count)")
    ap.add_argument("--aux-basis", default=AUX_BASIS_DEFAULT,
                    help=f"RI-MP2 fitting basis (default: {AUX_BASIS_DEFAULT}, "
                         "cuest_scf_settings_t's own default)")
    args = ap.parse_args()

    if args.method == "hf":
        if args.all_electron:
            ap.error("--all-electron only applies to mp2/ri-mp2/scs-mp2")
        scc = run_hf_scc(read_xyz(GEOM_FILE))
        print_hf(scc)
    else:
        result = run_correlated(args.method, args.all_electron, aux_basis=args.aux_basis)
        result["all_electron"] = args.all_electron
        print_correlated(args.method, result)


if __name__ == "__main__":
    sys.exit(main())
