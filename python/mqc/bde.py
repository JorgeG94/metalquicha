"""Bond dissociation energies, from a parent molecule and the bonds to cut.

    import mqc
    from mqc import bde

    with mqc.session():
        result = bde.scan("methanol.xyz", bonds="X-H")
        print(result.table())          # weakest bond first
        row = result.bonds[0]
        print(row.bde_kcal, row.d0_kcal, row.dg_kcal)

For each bond ``(i, j)`` (0-based) the molecule is split there into the two
connected pieces, each piece is taken through a conformer search (off by default
for fragments), an optimization, a Hessian and a single point, and ::

    D_e      = E(A) + E(B) - E(AB)
    D0       = D_e + ZPE(A) + ZPE(B) - ZPE(AB)
    BDE(298) = H(A) + H(B) - H(AB)           the one usually quoted
    dG       = G(A) + G(B) - G(AB)

are reported in kcal/mol and kJ/mol. The parent is computed once and each
distinct fragment once, however many bonds produce it: the three equivalent C-H
bonds of a methyl group are one calculation. See `mqc_docs/source/bde.rst` for
the definitions, the radical treatment and what to distrust.

**The default is GFN2-xTB in the gas phase**, unlike `mqc.pka`, whose default is
ALPB water. Homolysis gives doublet radicals (``multiplicity`` follows from the
electron count) and a parent that is charged must say which fragment keeps the
charge. ``heterolytic=True`` gives a cation and an anion and needs
``cation=``. A bond in a ring is refused: cutting it leaves one molecule.

**Open-shell xTB.** The Fortran xTB path hands tblite the number of unpaired
electrons (``uhf = multiplicity - 1``, `mqc_method_xtb.f90`) and asks for one
spin channel (``new_wavefunction(..., nspin=1, ...)``) without adding the
spin-polarization interaction, so a doublet is an unpolarized calculation with
the extra electron in a singly occupied orbital. It runs, and radicals are
xTB's weak point; the result says so in ``warnings``. No ``<S^2>`` is reported
for xTB, and none of the ab initio backends' ``spin_squared`` reaches the JSON.

**Electronic entropy.** The Fortran thermochemistry has ``S_elec = R ln(2S+1)``
but no caller gives it the multiplicity, so it reports 0 for a radical. The
spin term is added here, for every species (`mqc._thermo.free_energy_correction`),
so ``dG`` and the free energies carry it. It does not enter ``D_e``, ``D0`` or
``BDE(298)``, which are enthalpies.

**An atom** (an H radical, a halogen) is not given a Hessian: it is one single
point plus ``H = E + 5/2 RT`` and ``S = S_trans + R ln(2S+1)``, with spin-orbit
coupling ignored.

Run it inside ``mqc.session()``, exactly as `mqc.pka`: the conformer and
optimization stages run the ``mqc`` executable (see `mqc._workflow`), the
Hessians and single points use the session.
"""

import dataclasses
import os

from . import _bde_math, _workflow
from ._bde_math import (  # noqa: F401  (re-exported: this is the public surface)
    BDEError,
    BDEResult,
    BondPlan,
    BondRow,
    FragmentPlan,
    SpeciesData,
    perceive_bonds,
    split_bond,
)
from ._thermo import Geometry

__all__ = [
    "Protocol",
    "Geometry",
    "BDEResult",
    "BondRow",
    "SpeciesData",
    "BDEError",
    "BondPlan",
    "FragmentPlan",
    "scan",
    "split",
    "perceive_bonds",
    "find_executable",
]

# The names the stages are looked up under at call time, so that this module's
# can be replaced (as `mqc.pka`'s are in its tests).
_run_mqc = _workflow.run_mqc
_find_executable = _workflow.find_executable


@dataclasses.dataclass
class Protocol(_workflow.StageProtocol):
    """What to run at each stage. **GFN2-xTB, gas phase** unless told otherwise.

    The four stage dicts are keyword arguments of `mqc.MBE`, as in
    `mqc.pka.Protocol`, and the defaults differ from it in the two places that
    matter: there is no solvent, and ``standard_state`` is False, so ``dG`` is at
    1 atm for every species. Putting ``"xtb": {"solvent": "water",
    "solvation_model": "alpb"}`` in the stages gives a solution-phase BDE, which
    is a different quantity (see the docs); set ``standard_state=True`` as well
    if the free energy should be at 1 M.

    ``parent_conformers`` and ``fragment_conformers`` switch the CREST search
    for the parent and for the fragments; it is on for the parent and **off for
    the fragments**, which are small and usually radicals CREST has no reason to
    treat better than the parent's own geometry. ``conformers=None`` switches it
    off everywhere. Everything else is as `mqc._workflow.StageProtocol`.
    """

    standard_state: bool = False
    workdir: str = "bde_work"
    parent_conformers: bool = True
    fragment_conformers: bool = False


def find_executable(protocol=None):
    """The ``mqc`` executable the conformer and optimization stages will run."""
    return _find_executable(protocol or Protocol())


def _as_geometry(parent):
    if isinstance(parent, Geometry):
        return parent
    if isinstance(parent, (str, os.PathLike)):
        return Geometry.from_xyz(parent)
    if isinstance(parent, (tuple, list)) and len(parent) == 2:
        return Geometry(*parent)
    raise TypeError(
        "the parent must be a Geometry, an .xyz path or (symbols, coords): its coordinates are "
        "needed to find the bonds, and an mqc.System cannot be read back"
    )


def split(
    parent,
    i,
    j,
    charge=0,
    multiplicity=1,
    heterolytic=False,
    charge_on=None,
    cation=None,
    multiplicities=None,
    tolerance=_bde_math.DEFAULT_TOLERANCE,
):
    """The two fragments of cutting bond ``(i, j)``, without computing anything.

    Returns a `BondPlan`: ``side_i`` holds atom ``i`` and ``side_j`` atom ``j``,
    each a `FragmentPlan` with the parent atoms it keeps, its charge,
    multiplicity and the key it is deduplicated on. Raises `BDEError` for a bond
    that is not there, a ring bond, or charges and spins that cannot be made
    consistent -- which `scan` does before it runs anything.
    """
    geometry = _as_geometry(parent)
    return _bde_math.plan_bonds(
        geometry,
        [(int(i), int(j))],
        charge,
        multiplicity,
        tolerance,
        heterolytic,
        charge_on,
        cation,
        multiplicities,
    )[0]


def _species_data(species, formula, key, evaluation, parent=False):
    """`SpeciesData` from an evaluation, including the vertical and relaxed energies."""
    records = evaluation.detail["conformers"]
    best = max(records, key=lambda r: r["weight"])
    vertical, relaxation = None, None
    if not parent and evaluation.e_vertical_hartree is not None:
        vertical = evaluation.e_vertical_hartree
        relaxation = (vertical - best["e_single_point_hartree"]) * _bde_math.HARTREE_TO_KCAL
    return SpeciesData(
        name=species.name,
        formula=formula,
        key=key,
        charge=species.charge,
        multiplicity=species.multiplicity,
        n_atoms=species.n_atoms,
        e_hartree=evaluation.e_hartree,
        zpe_kcal=evaluation.zpe_kcal,
        h_kcal=evaluation.h_kcal,
        g_kcal=evaluation.g_kcal,
        temperature_K=evaluation.temperature_K,
        n_imaginary=evaluation.detail["n_imaginary"],
        n_conformers=evaluation.detail["n_conformers"],
        e_vertical_hartree=vertical,
        relaxation_kcal=relaxation,
        detail=evaluation.detail,
    )


def _uses_xtb(protocol):
    return any(
        str((getattr(protocol, stage) or {}).get("method", "")).lower().startswith("gfn")
        for stage in ("optimize", "frequencies", "single_point")
    )


def _has_solvent(protocol):
    return any(
        (getattr(protocol, stage) or {}).get("xtb", {}).get("solvent")
        or (getattr(protocol, stage) or {}).get("solvent")
        for stage in ("conformers", "optimize", "frequencies", "single_point")
    )


def scan(
    parent,
    bonds="X-H",
    protocol=None,
    charge=0,
    multiplicity=1,
    *,
    heterolytic=False,
    charge_on=None,
    cation=None,
    multiplicities=None,
    tolerance=_bde_math.DEFAULT_TOLERANCE,
    prefix="bde",
    verbose=True,
):
    """Dissociation energies for some or all bonds of one molecule: a `BDEResult`.

    ``parent`` is a `Geometry`, an .xyz path or ``(symbols, coords)``.
    ``bonds`` is ``"X-H"`` (every bond to hydrogen), ``"all"`` or a list of
    ``(i, j)`` 0-based pairs. A ring bond picked by a keyword is skipped and
    reported; one named in a list is an error. Everything that can be wrong
    about the request -- a missing bond, a ring, a charge without ``charge_on``,
    an impossible spin -- is raised before any calculation starts.

    ``charge`` and ``multiplicity`` are the parent's. Homolytic fragments take
    their multiplicity from the electron count (doublets from a closed-shell
    parent); ``multiplicities=(m_i, m_j)`` overrides it, and is required when
    the parent is not a singlet. A charged parent needs ``charge_on="i"|"j"``,
    the side of the bond that keeps the charge: "i" is the fragment holding the
    *first* atom of the pair as written (the lower index for ``"X-H"`` and
    ``"all"``). ``heterolytic=True`` needs ``cation="i"|"j"`` in the same sense.

    The parent is computed once, then each distinct fragment once, starting from
    the parent's relaxed coordinates. The result is sorted by ``BDE(298)``.
    Sequential, rank 0, inside ``mqc.session()``.
    """
    protocol = protocol or Protocol()
    geometry = _as_geometry(parent)
    selected, skipped_rings = _bde_math.select_bonds(geometry, bonds, tolerance)
    if not selected:
        raise BDEError("no bonds to scan: the selection matched none that is not in a ring")
    plans = _bde_math.plan_bonds(
        geometry, selected, charge, multiplicity, tolerance, heterolytic, charge_on, cation, multiplicities
    )
    notes = []
    skipped = []
    for i, j in skipped_rings:
        skipped.append((i, j, f"{geometry.symbols[i]}-{geometry.symbols[j]} is in a ring: cutting it leaves one piece"))
    if skipped:
        notes.append(f"{len(skipped)} ring bond(s) were left out of the selection; see `skipped`")

    # Unique fragments, in the order the bonds first produce them. Two different
    # species that would share a (hash-stub) name get a counter.
    unique = {}
    names = {}
    for plan in plans:
        for side in (plan.side_i, plan.side_j):
            if side.key in unique:
                continue
            name = side.name
            while name in names and names[name] != side.key:
                name += "x"
            names[name] = side.key
            side.name = name
            unique[side.key] = side
    for plan in plans:  # names settled: point every side at its canonical one
        plan.side_i.name = unique[plan.side_i.key].name
        plan.side_j.name = unique[plan.side_j.key].name

    parent_formula = _bde_math.hill_formula(geometry.symbols)
    parent_species = _workflow.Species(f"{parent_formula}_parent", geometry, charge, multiplicity)

    def evaluate(species, use_conformers, vertical):
        return _workflow.evaluate_species(
            species,
            protocol,
            prefix,
            verbose,
            use_conformers=use_conformers,
            vertical=vertical,
            run=lambda *args: _run_mqc(*args),
            locate_executable=lambda p: _find_executable(p),
        )

    if verbose:
        print(f"parent {parent_formula}: {len(plans)} bond(s), {len(unique)} distinct fragment(s)", flush=True)
    parent_eval = evaluate(parent_species, protocol.parent_conformers, False)
    parent_data = _species_data(
        parent_species, parent_formula, _bde_math.species_key(geometry.symbols, perceive_bonds(geometry, tolerance), charge, multiplicity), parent_eval, parent=True
    )
    relaxed = parent_eval.lowest_geometry()
    if relaxed is not geometry and _bde_math.connectivity_changed(geometry, relaxed, tolerance):
        notes.append(
            "the relaxed parent has a different bond graph from the input at this tolerance; the fragments "
            "are cut along the INPUT connectivity and take the relaxed parent's coordinates"
        )

    species = {}
    vertical = protocol.optimize is not None
    for key, side in unique.items():
        fragment = _workflow.Species(
            side.name, _bde_math.fragment_geometry(relaxed, side.atoms), side.charge, side.multiplicity
        )
        evaluation = evaluate(fragment, protocol.fragment_conformers, vertical)
        species[side.name] = _species_data(fragment, _bde_math.hill_formula(side.symbols), key, evaluation)

    rows = []
    for plan in plans:
        a, b = species[plan.side_i.name], species[plan.side_j.name]
        d_e, d0, bde, dg = _bde_math.assemble_bond(parent_data, a, b)
        rows.append(
            BondRow(
                plan.i, plan.j, f"{geometry.symbols[plan.i]}-{geometry.symbols[plan.j]}",
                plan.side_i.atoms, plan.side_j.atoms, a.name, b.name,
                d_e, d0, bde, dg, parent_data.temperature_K, plan.heterolytic,
            )
        )

    all_species = [parent_data] + list(species.values())
    if any(s.multiplicity > 1 for s in all_species):
        if _uses_xtb(protocol):
            notes.append(
                "open-shell species were run as GFN-xTB with the unpaired electrons set and no spin "
                "polarization, which is xTB's weak point: expect radical stabilization to be off by "
                "several kcal/mol. No <S^2> is available to check for contamination."
            )
        else:
            notes.append("no <S^2> is read back for the open-shell species; check spin contamination in the output files")
    if any(s.n_atoms == 1 for s in species.values()):
        notes.append("single-atom fragments use H = E + 5/2 RT and S = S_trans + R ln(2S+1); spin-orbit coupling is ignored")
    if heterolytic and not _has_solvent(protocol):
        notes.append(
            "a gas-phase heterolytic dissociation separates charges with nothing to screen them, so it is "
            "hundreds of kcal/mol and not comparable with a homolytic one; use a solvent in the protocol"
        )
    if _has_solvent(protocol):
        notes.append(
            "a solvent is set: these are solution-phase dissociation enthalpies (the solvation model's "
            "free energy is in every energy, with no solvation enthalpy separated out), not gas-phase BDEs"
        )
    return BDEResult(
        parent_data, species, rows, parent_data.temperature_K, protocol=protocol.as_dict(), skipped=skipped, notes=notes
    )
