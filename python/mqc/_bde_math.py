"""The arithmetic and bookkeeping of a bond-dissociation-energy calculation, with no library in it.

Fragment generation (connectivity, splitting at a bond, charges and spins),
the key that decides two fragments are the same species, the assembly of
``D_e``, ``D0``, ``BDE(T)`` and ``dG`` from species energies, and the result
object. Nothing here imports `mqc` or loads `libmqc`, so `python/tests/test_bde.py`
runs without a build. The thermochemistry is `mqc._thermo`; the staged
evaluation of each species is `mqc._workflow`; `mqc.bde` joins them.

Units: energies in kcal/mol unless named ``_hartree`` or ``_kj``; geometries in
Angstrom; atom indices 0-based, as everywhere in the API.
"""

import dataclasses
import hashlib
import json
import math

try:
    from . import _thermo
except ImportError:  # loaded by path, as python/tests/test_bde.py does: no parent package
    import importlib.util as _util
    import os as _os
    import sys as _sys

    _spec = _util.spec_from_file_location(
        "mqc_thermo_under_bde_test", _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "_thermo.py")
    )
    _thermo = _util.module_from_spec(_spec)
    _sys.modules[_spec.name] = _thermo
    _spec.loader.exec_module(_thermo)

HARTREE_TO_KCAL = _thermo.HARTREE_TO_KCAL
KCAL_TO_KJ = _thermo.KCAL_TO_KJ

#: Slack on the radius sum, the value `System.perceive_bonds` and
#: `DEFAULT_BOND_TOLERANCE` in `mqc_bond_perception.f90` use.
DEFAULT_TOLERANCE = 1.2


class BDEError(_thermo.WorkflowError):
    """A bond-dissociation calculation that cannot be set up, or has no answer."""


# -- connectivity -----------------------------------------------------------

#: Covalent radii in Angstrom, Z = 1..96, copied from `CORDERO_RADII` in
#: `src/core/mqc_atomic_radii.f90` (Cordero et al., Dalton Trans. 2008, 2832),
#: which is what `element_covalent_radius` -- and so `System.perceive_bonds` --
#: uses. An element past curium has none and bonds to nothing, as there.
_CORDERO = (
    0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66,
    0.57, 0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05,
    1.02, 1.06, 2.03, 1.76, 1.7, 1.6, 1.53, 1.39,
    1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22, 1.2,
    1.19, 1.2, 1.2, 1.16, 2.2, 1.95, 1.9, 1.75,
    1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44,
    1.42, 1.39, 1.39, 1.38, 1.39, 1.4, 2.44, 2.15,
    2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96,
    1.94, 1.92, 1.92, 1.89, 1.9, 1.87, 1.87, 1.75,
    1.7, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32,
    1.45, 1.46, 1.48, 1.4, 1.5, 1.5, 2.6, 2.21,
    2.15, 2.06, 2.0, 1.96, 1.9, 1.87, 1.8, 1.69,
)


def covalent_radius(symbol):
    """Cordero covalent radius in Angstrom; 0 where none is tabulated."""
    z = _thermo.atomic_number(symbol)
    return _CORDERO[z - 1] if z <= len(_CORDERO) else 0.0


def perceive_bonds(geometry, tolerance=DEFAULT_TOLERANCE):
    """Bonded atom pairs ``(i, j)``, ``i < j``, by ``d < tolerance * (r_i + r_j)``.

    The rule of `mqc_bond_perception.f90`: Cordero radii, 1.2 by default, an
    element with no radius bonds to nothing. A heuristic: it has no bond
    orders and is wrong for hydrogen bonds a little shorter than 1.8 A and for
    long dative bonds.
    """
    radii = [covalent_radius(s) for s in geometry.symbols]
    bonds = []
    n = geometry.n_atoms
    for i in range(n):
        for j in range(i + 1, n):
            if radii[i] <= 0.0 or radii[j] <= 0.0:
                continue
            d = math.dist(geometry.coords[i], geometry.coords[j])
            if d < tolerance * (radii[i] + radii[j]):
                bonds.append((i, j))
    return bonds


def adjacency(n_atoms, bonds):
    """Neighbour lists, as ``[set, ...]``."""
    nbrs = [set() for _ in range(n_atoms)]
    for i, j in bonds:
        nbrs[i].add(j)
        nbrs[j].add(i)
    return nbrs


def _reach(nbrs, start, skip_edge=None):
    """Atoms reachable from ``start``, optionally with one bond removed."""
    seen = {start}
    stack = [start]
    while stack:
        a = stack.pop()
        for b in nbrs[a]:
            if skip_edge is not None and {a, b} == set(skip_edge):
                continue
            if b not in seen:
                seen.add(b)
                stack.append(b)
    return seen


def split_bond(geometry, i, j, tolerance=DEFAULT_TOLERANCE):
    """Cut the bond ``i-j``: ``(atoms_on_i_side, atoms_on_j_side)``, sorted index lists.

    Raises `BDEError` when ``i`` and ``j`` are not bonded at this tolerance,
    when the molecule is not one connected piece, and when the bond is in a
    ring -- removing it then leaves one molecule, a ring-opened diradical, not
    two fragments.
    """
    n = geometry.n_atoms
    i, j = int(i), int(j)
    if not (0 <= i < n and 0 <= j < n) or i == j:
        raise BDEError(f"bond ({i}, {j}): atom indices must be distinct and in 0..{n - 1}")
    nbrs = adjacency(n, perceive_bonds(geometry, tolerance))
    label = f"{geometry.symbols[i]}{i}-{geometry.symbols[j]}{j}"
    if j not in nbrs[i]:
        d = math.dist(geometry.coords[i], geometry.coords[j])
        raise BDEError(
            f"bond {label}: the atoms are {d:.3f} A apart and not bonded at tolerance {tolerance:g} "
            "(d < tolerance * sum of covalent radii). Check the indices, or raise `tolerance`."
        )
    if len(_reach(nbrs, 0)) != n:
        raise BDEError(
            "the parent is not one connected molecule at this tolerance; a bond dissociation "
            "energy needs a single molecule (give the fragments of a complex separately)"
        )
    side_i = _reach(nbrs, i, skip_edge=(i, j))
    if j in side_i:
        raise BDEError(
            f"bond {label} is in a ring: cutting it leaves one connected piece, so there are no "
            "two fragments to compute. A ring-opening energy is a different quantity; compute "
            "the open-chain diradical yourself and take the difference."
        )
    return sorted(side_i), sorted(set(range(n)) - side_i)


def select_bonds(geometry, spec="X-H", tolerance=DEFAULT_TOLERANCE):
    """Which bonds a scan covers: ``(bonds, skipped_ring_bonds)``.

    ``spec`` is ``"X-H"`` (every bond to a hydrogen), ``"all"``, or an explicit
    list of ``(i, j)``. A ring bond picked by a keyword is *skipped* and
    returned so that the caller can say so; one in an explicit list is an error,
    since the caller asked for it by name.
    """
    perceived = perceive_bonds(geometry, tolerance)
    nbrs = adjacency(geometry.n_atoms, perceived)

    def is_ring(i, j):
        return j in _reach(nbrs, i, skip_edge=(i, j))

    if isinstance(spec, str):
        key = spec.strip().lower().replace(" ", "")
        if key in ("x-h", "xh", "h"):
            chosen = [b for b in perceived if "H" in (geometry.symbols[b[0]], geometry.symbols[b[1]])]
        elif key == "all":
            chosen = list(perceived)
        else:
            raise ValueError(f"bonds must be 'X-H', 'all' or a list of (i, j) pairs, not {spec!r}")
        keep = [b for b in chosen if not is_ring(*b)]
        return keep, [b for b in chosen if is_ring(*b)]
    bonds = []
    for pair in spec:
        i, j = (int(pair[0]), int(pair[1]))
        # Kept in the order written: which side is "i" decides where a charge goes.
        if (i, j) not in bonds and (j, i) not in bonds:
            bonds.append((i, j))
    if not bonds:
        raise ValueError("no bonds to scan")
    for i, j in bonds:
        split_bond(geometry, i, j, tolerance)  # raises, with the reason, for a bad one
    return bonds, []


# -- charges and spins ------------------------------------------------------


def electron_count(symbols, charge):
    """Electrons of a species: sum of Z minus the charge."""
    return sum(_thermo.atomic_number(s) for s in symbols) - int(charge)


def _check_state(name, symbols, charge, multiplicity):
    n_el = electron_count(symbols, charge)
    if n_el < 0 or (multiplicity - 1) > n_el or (n_el - (multiplicity - 1)) % 2:
        raise BDEError(
            f"{name}: {n_el} electrons cannot have multiplicity {multiplicity} "
            f"(charge {charge:+d}). Give `multiplicities=` or `charge_on=` that agree."
        )


def _natural_multiplicity(symbols, charge):
    return 2 if electron_count(symbols, charge) % 2 else 1


def fragment_states(
    symbols_i,
    symbols_j,
    parent_charge=0,
    parent_multiplicity=1,
    heterolytic=False,
    charge_on=None,
    cation=None,
    multiplicities=None,
):
    """Charge and multiplicity of the two fragments: ``((q_i, m_i), (q_j, m_j))``.

    Homolysis shares the parent's charge as it is and each electron of the bond
    goes to one side, so the fragments are radicals (doublets) when the parent
    is a neutral closed shell. A **nonzero parent charge is never guessed**:
    ``charge_on`` (``"i"`` or ``"j"``, the side of the cut that keeps it) is
    required. ``heterolytic=True`` moves the bond pair to one side, which needs
    ``cation`` (``"i"`` or ``"j"``): that fragment is the cation (charge +1
    relative to the homolytic split), the other the anion (-1).

    The multiplicities follow from the electron count -- 2 for an odd number of
    electrons, else 1 -- which is the answer only for a closed-shell parent, so
    a parent with a nonzero spin needs ``multiplicities=(m_i, m_j)``. Whatever
    is given is checked against the electron count.
    """
    for name, side in (("charge_on", charge_on), ("cation", cation)):
        if side not in (None, "i", "j"):
            raise ValueError(f"{name} must be 'i' or 'j' (the side of the bond holding that atom), not {side!r}")
    parent_charge = int(parent_charge)
    if parent_charge:
        if charge_on is None:
            raise BDEError(
                f"the parent has charge {parent_charge:+d}; say which fragment keeps it with "
                "charge_on='i' (the side of the first atom of the bond) or 'j'"
            )
        q = (parent_charge, 0) if charge_on == "i" else (0, parent_charge)
    else:
        q = (0, 0)
    if heterolytic:
        if cation is None:
            raise BDEError("a heterolytic cleavage needs cation='i' or 'j': which fragment loses the electron pair")
        q = (q[0] + 1, q[1] - 1) if cation == "i" else (q[0] - 1, q[1] + 1)
    elif cation is not None:
        raise BDEError("cation= is for heterolytic=True; a homolytic cleavage has no cation")

    if multiplicities is None:
        if int(parent_multiplicity) != 1:
            raise BDEError(
                f"the parent has multiplicity {parent_multiplicity}, so how the spin divides between the "
                "fragments is a choice: give multiplicities=(m_i, m_j)"
            )
        m = (_natural_multiplicity(symbols_i, q[0]), _natural_multiplicity(symbols_j, q[1]))
        if not heterolytic and (m[0] == 1 or m[1] == 1):
            # One electron of the bond pair to each side makes both pieces odd. Two even pieces
            # means the pair went to one side: that is a heterolytic split under another name.
            raise BDEError(
                f"with charge_on={charge_on!r} the fragments have {electron_count(symbols_i, q[0])} and "
                f"{electron_count(symbols_j, q[1])} electrons, so both are closed shells and the bond "
                "pair has gone to one side: that is a heterolytic cleavage. Put the charge on the other "
                "fragment, pass heterolytic=True, or give multiplicities= if you mean something else."
            )
    else:
        m = (int(multiplicities[0]), int(multiplicities[1]))
    _check_state("fragment i", symbols_i, q[0], m[0])
    _check_state("fragment j", symbols_j, q[1], m[1])
    return (q[0], m[0]), (q[1], m[1])


# -- species identity -------------------------------------------------------


def hill_formula(symbols):
    """A Hill-order formula: C, H, then the rest alphabetically; plain alphabetical without C."""
    counts = {}
    for s in symbols:
        counts[s] = counts.get(s, 0) + 1
    order = sorted(counts)
    if "C" in counts:
        order = ["C"] + (["H"] if "H" in counts else []) + [s for s in order if s not in ("C", "H")]
    return "".join(f"{s}{counts[s] if counts[s] > 1 else ''}" for s in order)


def connectivity_hash(symbols, bonds):
    """A canonical hash of the bond graph, from colour refinement.

    Each atom starts as its element and is repeatedly relabelled by its own
    label and the sorted labels of its neighbours, until the number of distinct
    labels stops growing; the hash is of the sorted final labels. It does not
    depend on the order the atoms are listed in.

    **What it cannot see:** it is a graph invariant, not a canonical form, so
    two non-isomorphic graphs that colour refinement cannot tell apart share a
    hash (some regular graphs, such as the 3-regular graphs on six vertices
    with the same colours; none of those is a likely fragment of a small
    molecule, but the guarantee is absent). It carries no geometry, no
    stereochemistry, no bond orders and no conformer, so cis/trans isomers and
    enantiomers are one species here. The species key adds element counts,
    charge and multiplicity on top of it.
    """
    nbrs = adjacency(len(symbols), bonds)
    labels = [str(s) for s in symbols]
    n_classes = len(set(labels))
    for _ in range(len(symbols)):
        new = [
            hashlib.sha1((labels[a] + "|" + ",".join(sorted(labels[b] for b in nbrs[a]))).encode()).hexdigest()[:16]
            for a in range(len(symbols))
        ]
        labels = new
        classes = len(set(labels))
        if classes == n_classes:
            break
        n_classes = classes
    return hashlib.sha1("".join(sorted(labels)).encode()).hexdigest()[:12]


def species_key(symbols, bonds, charge, multiplicity):
    """The identity two fragments are deduplicated on: ``formula|q|m|hash``."""
    return f"{hill_formula(symbols)}|q{int(charge):+d}|m{int(multiplicity)}|{connectivity_hash(symbols, bonds)}"


def species_name(symbols, charge, multiplicity, key):
    """A readable, filename-safe name: formula, charge and spin, and a hash stub."""
    q = "" if not charge else (f"{abs(int(charge))}" if abs(int(charge)) > 1 else "") + ("p" if charge > 0 else "m")
    return f"{hill_formula(symbols)}{q}_m{int(multiplicity)}_{key.rsplit('|', 1)[1][:6]}"


@dataclasses.dataclass
class FragmentPlan:
    """One side of a cut: which parent atoms, and what it is."""

    atoms: list
    symbols: list
    charge: int
    multiplicity: int
    key: str
    name: str


@dataclasses.dataclass
class BondPlan:
    """A bond to cut and what it cuts into. ``side_i`` holds atom ``i``."""

    i: int
    j: int
    label: str
    side_i: FragmentPlan
    side_j: FragmentPlan
    heterolytic: bool = False


def plan_bonds(
    geometry,
    bonds,
    parent_charge=0,
    parent_multiplicity=1,
    tolerance=DEFAULT_TOLERANCE,
    heterolytic=False,
    charge_on=None,
    cation=None,
    multiplicities=None,
):
    """`BondPlan` for each bond: the split, the states, the keys. No calculation."""
    perceived = perceive_bonds(geometry, tolerance)
    _check_state("parent", geometry.symbols, int(parent_charge), int(parent_multiplicity))
    plans = []
    for i, j in bonds:
        atoms_i, atoms_j = split_bond(geometry, i, j, tolerance)
        sym_i = [geometry.symbols[a] for a in atoms_i]
        sym_j = [geometry.symbols[a] for a in atoms_j]
        (qi, mi), (qj, mj) = fragment_states(
            sym_i, sym_j, parent_charge, parent_multiplicity, heterolytic, charge_on, cation, multiplicities
        )
        sides = []
        for atoms, syms, q, m in ((atoms_i, sym_i, qi, mi), (atoms_j, sym_j, qj, mj)):
            local = {a: k for k, a in enumerate(atoms)}
            inner = [(local[a], local[b]) for a, b in perceived if a in local and b in local]
            key = species_key(syms, inner, q, m)
            sides.append(FragmentPlan(atoms, syms, q, m, key, species_name(syms, q, m, key)))
        label = f"{geometry.symbols[i]}-{geometry.symbols[j]}"
        plans.append(BondPlan(i, j, label, sides[0], sides[1], bool(heterolytic)))
    return plans


def fragment_geometry(geometry, atoms):
    """The sub-geometry of the given atoms, in order, as the parent has them."""
    return type(geometry)([geometry.symbols[a] for a in atoms], [geometry.coords[a] for a in atoms])


def connectivity_changed(reference, other, tolerance=DEFAULT_TOLERANCE):
    """Whether two geometries of one atom list have different bond graphs."""
    return set(perceive_bonds(reference, tolerance)) != set(perceive_bonds(other, tolerance))


# -- energies ---------------------------------------------------------------


@dataclasses.dataclass
class SpeciesData:
    """The numbers one species contributes: parent or fragment.

    ``e_hartree`` is the electronic energy at the single-point level,
    ``zpe_kcal`` the zero-point energy, ``h_kcal`` the enthalpy and ``g_kcal`` the
    free energy (electronic energy included, kcal/mol). ``e_vertical_hartree`` is
    the energy at the geometry the fragment had inside the parent, and
    ``relaxation_kcal`` its drop on relaxing (0 for an atom). Both are None for
    the parent and for a run that did not optimize.
    """

    name: str
    formula: str
    key: str
    charge: int
    multiplicity: int
    n_atoms: int
    e_hartree: float
    zpe_kcal: float
    h_kcal: float
    g_kcal: float
    temperature_K: float
    n_imaginary: int = 0
    n_conformers: int = 1
    e_vertical_hartree: float = None
    relaxation_kcal: float = None
    detail: dict = dataclasses.field(default_factory=dict)

    def to_dict(self):
        d = dataclasses.asdict(self)
        return d


@dataclasses.dataclass
class BondRow:
    """One bond's result. Energies are kcal/mol; ``*_kj`` properties give kJ/mol.

    ``bde_kcal`` is the enthalpy of dissociation at ``temperature_K`` -- the
    conventional ``BDE(298)`` -- ``d0_kcal`` the same at 0 K without thermal
    energy (``D_e`` plus the change in zero-point energy) and ``dg_kcal`` the
    free energy of dissociation, 1 atm standard state unless the protocol turned
    the 1 M term on.
    """

    i: int
    j: int
    label: str
    atoms_i: list
    atoms_j: list
    species_i: str
    species_j: str
    d_e_kcal: float
    d0_kcal: float
    bde_kcal: float
    dg_kcal: float
    temperature_K: float
    heterolytic: bool = False

    @property
    def d_e_kj(self):
        return self.d_e_kcal * KCAL_TO_KJ

    @property
    def d0_kj(self):
        return self.d0_kcal * KCAL_TO_KJ

    @property
    def bde_kj(self):
        return self.bde_kcal * KCAL_TO_KJ

    @property
    def dg_kj(self):
        return self.dg_kcal * KCAL_TO_KJ

    def to_dict(self):
        d = dataclasses.asdict(self)
        for key in ("d_e", "d0", "bde", "dg"):
            d[f"{key}_kj"] = getattr(self, f"{key}_kj")
        return d


def assemble_bond(parent, a, b):
    """``(D_e, D0, BDE(T), dG)`` in kcal/mol for ``parent -> a + b``.

    ``D_e = E(a) + E(b) - E(parent)``; ``D0 = D_e + dZPE``;
    ``BDE(T) = H(a) + H(b) - H(parent)``; ``dG = G(a) + G(b) - G(parent)``.
    """
    if not (parent.temperature_K == a.temperature_K == b.temperature_K):
        raise BDEError(
            "the parent and fragments were evaluated at different temperatures "
            f"({parent.temperature_K}, {a.temperature_K}, {b.temperature_K} K)"
        )
    d_e = (a.e_hartree + b.e_hartree - parent.e_hartree) * HARTREE_TO_KCAL
    d0 = d_e + (a.zpe_kcal + b.zpe_kcal - parent.zpe_kcal)
    bde = a.h_kcal + b.h_kcal - parent.h_kcal
    dg = a.g_kcal + b.g_kcal - parent.g_kcal
    return d_e, d0, bde, dg


class BDEResult:
    """Per-bond rows (weakest first) and per-species data.

    ``bonds`` are `BondRow`, sorted by ``bde_kcal``; ``parent`` and ``species``
    (``name -> SpeciesData``, each fragment once however many bonds share it)
    carry the energies behind them. ``skipped`` lists ``(i, j, reason)`` for bonds
    a keyword selection left out, and ``notes`` are warnings the scan itself
    raised. ``protocol`` is the protocol as plain data.
    """

    def __init__(self, parent, species, bonds, temperature, protocol=None, skipped=None, notes=None):
        self.parent = parent
        self.species = dict(species)
        self.bonds = sorted(bonds, key=lambda r: r.bde_kcal)
        self.temperature = float(temperature)
        self.protocol = protocol or {}
        self.skipped = list(skipped or [])
        self.notes = list(notes or [])

    def _all_species(self):
        return [self.parent] + list(self.species.values())

    @property
    def warnings(self):
        """What to distrust, in words."""
        out = list(self.notes)
        for sp in self._all_species():
            if sp.n_imaginary:
                out.append(
                    f"{sp.name}: {sp.n_imaginary} imaginary frequency(ies) across its conformers; "
                    "it is not a minimum, so its enthalpy and free energy are not a minimum's"
                )
            if sp.relaxation_kcal is not None and sp.relaxation_kcal < -0.1:
                out.append(
                    f"{sp.name}: the relaxed energy is {-sp.relaxation_kcal:.2f} kcal/mol ABOVE the vertical one; "
                    "the optimization or the conformer choice has gone somewhere worse"
                )
        for row in self.bonds:
            if row.bde_kcal < 0.0:
                out.append(f"bond {row.label} ({row.i}-{row.j}): BDE {row.bde_kcal:.1f} kcal/mol is negative")
        return out

    def table(self):
        """A text table of the bonds, weakest first, with the species behind it."""
        t = self.temperature
        p = self.parent
        lines = [
            f"Bond dissociation, parent {p.formula} (q {p.charge:+d}, m {p.multiplicity}), T = {t:.2f} K",
            f"energies in kcal/mol; BDE({t:.0f}) also in kJ/mol",
            "",
        ]
        head = f"{'bond':<8}{'atoms':<9}{'D_e':>8}{'D0':>8}{'BDE(' + format(t, '.0f') + ')':>10}{'dG':>8}{'kJ/mol':>9}  fragments"
        lines.append(head)
        lines.append("-" * len(head))
        for r in self.bonds:
            frag = f"{self.species[r.species_i].formula}{self._tag(r.species_i)} + {self.species[r.species_j].formula}{self._tag(r.species_j)}"
            lines.append(
                f"{r.label:<8}{f'{r.i}-{r.j}':<9}{r.d_e_kcal:>8.1f}{r.d0_kcal:>8.1f}{r.bde_kcal:>10.1f}"
                f"{r.dg_kcal:>8.1f}{r.bde_kj:>9.1f}  {frag}"
            )
        for i, j, reason in self.skipped:
            lines.append(f"{'skipped':<8}{f'{i}-{j}':<9}{reason}")
        notes = self.warnings
        if notes:
            lines.append("")
            lines.extend(f"warning: {w}" for w in notes)
        return "\n".join(lines)

    def _tag(self, name):
        sp = self.species[name]
        return ("" if not sp.charge else f"({sp.charge:+d})") + ("*" if sp.multiplicity > 1 else "")

    def to_dict(self):
        return {
            "temperature_K": self.temperature,
            "parent": self.parent.to_dict(),
            "species": {k: v.to_dict() for k, v in self.species.items()},
            "bonds": [r.to_dict() for r in self.bonds],
            "skipped": [{"i": i, "j": j, "reason": reason} for i, j, reason in self.skipped],
            "warnings": self.warnings,
            "protocol": self.protocol,
        }

    def to_json(self, path=None, indent=2):
        text = json.dumps(self.to_dict(), indent=indent, default=_json_default)
        if path is not None:
            with open(path, "w") as handle:
                handle.write(text + "\n")
        return text

    def __repr__(self):
        return f"<BDEResult {self.parent.formula}: {len(self.bonds)} bond(s)>"


def _json_default(obj):
    if hasattr(obj, "to_dict"):
        return obj.to_dict()
    return str(obj)
