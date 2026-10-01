"""Tests for the bond-dissociation workflow that need neither the library nor MPI nor a build.

    python3 python/tests/test_bde.py

Exits non-zero on the first failure; also collected by pytest if that is to hand.
As in `test_pka.py`, nothing here imports `mqc` the ordinary way (that loads
`libmqc`):

  * the arithmetic, `mqc/_bde_math.py` and `mqc/_thermo.py`, is loaded by path --
    neither imports anything from the package;
  * `mqc/bde.py` is imported under the stand-in package of `_standin.py`, whose
    `MBE` and `_check_label` are the real ones and whose physics is invented.

What this cannot test is anything on the other side of the C interface: that a
doublet GFN2 run converges, that the Hessian of the radical writes the keys read
here, that the optimizer leaves the file looked for. Those are checked by running
`python/examples/bde.py` against a real build.
"""

import importlib.util
import json
import math
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(HERE, "..", "mqc")
sys.path.insert(0, HERE)

from _standin import FakeResult, fake_package, fake_thermo  # noqa: E402


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(PKG, filename))
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


TH = _load("mqc_thermo_under_test", "_thermo.py")
B = _load("mqc_bde_math_under_test", "_bde_math.py")
T = 298.15
RT = TH.R_KCAL * T


def close(a, b, tol=1e-9, msg=""):
    assert abs(a - b) <= tol, f"{msg} {a!r} != {b!r} (|diff| {abs(a - b):.3e} > {tol:.1e})"


def raises(exc, fn, *args, match=None, **kwargs):
    try:
        fn(*args, **kwargs)
    except exc as error:
        if match is not None:
            assert match in str(error), f"{str(error)!r} does not mention {match!r}"
        return error
    raise AssertionError(f"{exc.__name__} not raised")


# -- molecules, Angstrom ----------------------------------------------------

#: C0 O1 H2 H3 H4 (on C) H5 (on O)
METHANOL = (
    ["C", "O", "H", "H", "H", "H"],
    [[-0.047, 0.664, 0.0], [-0.047, -0.758, 0.0], [-1.092, 0.969, 0.0],
     [0.439, 1.073, 0.890], [0.439, 1.073, -0.890], [0.858, -1.070, 0.0]],
)

ETHANE = (
    ["C", "C", "H", "H", "H", "H", "H", "H"],
    [[0.0, 0.0, 0.765], [0.0, 0.0, -0.765],
     [1.019, 0.0, 1.158], [-0.509, 0.882, 1.158], [-0.509, -0.882, 1.158],
     [-1.019, 0.0, -1.158], [0.509, -0.882, -1.158], [0.509, 0.882, -1.158]],
)

METHANE = (
    ["C", "H", "H", "H", "H"],
    [[0.0, 0.0, 0.0], [0.629, 0.629, 0.629], [-0.629, -0.629, 0.629],
     [-0.629, 0.629, -0.629], [0.629, -0.629, -0.629]],
)


def cyclohexane():
    """Chair, ring carbons first then a pair of hydrogens (axial, equatorial) on each."""
    symbols, coords = [], []
    for k in range(6):
        th = math.radians(60.0 * k)
        symbols.append("C")
        coords.append([1.45 * math.cos(th), 1.45 * math.sin(th), 0.25 if k % 2 == 0 else -0.25])
    for k in range(6):
        th = math.radians(60.0 * k)
        up = 1.0 if k % 2 == 0 else -1.0
        x, y, z = coords[k]
        symbols += ["H", "H"]
        coords.append([x, y, z + up * 1.09])
        coords.append([x + 1.03 * math.cos(th), y + 1.03 * math.sin(th), z - up * 0.36])
    return symbols, coords


def benzene():
    symbols, coords = [], []
    for k in range(6):
        th = math.radians(60.0 * k)
        symbols.append("C")
        coords.append([1.39 * math.cos(th), 1.39 * math.sin(th), 0.0])
    for k in range(6):
        th = math.radians(60.0 * k)
        symbols.append("H")
        coords.append([2.48 * math.cos(th), 2.48 * math.sin(th), 0.0])
    return symbols, coords


def geom(spec):
    return TH.Geometry(*spec)


# ---------------------------------------------------------------------------
#  Connectivity and fragments
# ---------------------------------------------------------------------------


def test_radii_and_connectivity_match_the_fortran_rule():
    # Cordero values, as in src/core/mqc_atomic_radii.f90, and the 1.2 tolerance.
    assert (B.covalent_radius("H"), B.covalent_radius("C"), B.covalent_radius("O")) == (0.31, 0.76, 0.66)
    assert B.covalent_radius("Cl") == 1.02 and len(B._CORDERO) == 96
    assert B.DEFAULT_TOLERANCE == 1.2
    assert B.perceive_bonds(geom(METHANOL)) == [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5)]
    # d < 1.2 (r_i + r_j): C-H at 1.09 A binds (limit 1.284) and stops at 1.30 A.
    near = TH.Geometry(["C", "H"], [[0, 0, 0], [0, 0, 1.28]])
    far = TH.Geometry(["C", "H"], [[0, 0, 0], [0, 0, 1.30]])
    assert B.perceive_bonds(near) == [(0, 1)] and B.perceive_bonds(far) == []
    assert len(B.perceive_bonds(geom(cyclohexane()))) == 18  # 6 C-C and 12 C-H, nothing across the ring
    assert len(B.perceive_bonds(geom(benzene()))) == 12


def test_split_methanol_oh_and_co():
    g = geom(METHANOL)
    atoms_c, atoms_h = B.split_bond(g, 1, 5)
    assert atoms_c == [0, 1, 2, 3, 4] and atoms_h == [5]
    atoms_c, atoms_o = B.split_bond(g, 0, 1)
    assert atoms_c == [0, 2, 3, 4] and atoms_o == [1, 5]

    oh = B.plan_bonds(g, [(1, 5)])[0]
    assert oh.label == "O-H"
    assert (oh.side_i.charge, oh.side_i.multiplicity) == (0, 2) and (oh.side_j.charge, oh.side_j.multiplicity) == (0, 2)
    assert oh.side_i.symbols == ["C", "O", "H", "H", "H"] and oh.side_j.symbols == ["H"]
    co = B.plan_bonds(g, [(0, 1)])[0]
    assert [s.symbols for s in (co.side_i, co.side_j)] == [["C", "H", "H", "H"], ["O", "H"]]
    assert all(s.multiplicity == 2 for s in (co.side_i, co.side_j))
    # Fragment geometry keeps the parent's coordinates for its atoms.
    sub = B.fragment_geometry(g, oh.side_i.atoms)
    assert sub.symbols == oh.side_i.symbols and sub.coords[1] == g.coords[1]


def test_split_ethane_cc_gives_two_identical_methyls():
    g = geom(ETHANE)
    plan = B.plan_bonds(g, [(0, 1)])[0]
    assert plan.side_i.atoms == [0, 2, 3, 4] and plan.side_j.atoms == [1, 5, 6, 7]
    assert plan.side_i.key == plan.side_j.key and plan.side_i.multiplicity == 2


def test_ring_bonds_are_refused():
    for spec, bond in ((cyclohexane(), (0, 1)), (benzene(), (0, 1)), (benzene(), (2, 3))):
        raises(B.BDEError, B.split_bond, geom(spec), *bond, match="ring")
    # ... while a C-H on the ring is an ordinary bond, and so is a substituent bond.
    assert B.split_bond(geom(benzene()), 0, 6) == ([0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11], [6])

    # A keyword selection skips ring bonds and says which; a list refuses them.
    cyc = geom(cyclohexane())
    bonds, skipped = B.select_bonds(cyc, "all")
    assert len(skipped) == 6 and len(bonds) == 12 and all(cyc.symbols[a] == "H" or cyc.symbols[b] == "H" for a, b in bonds)
    bonds, skipped = B.select_bonds(cyc, "X-H")
    assert len(bonds) == 12 and skipped == []
    raises(B.BDEError, B.select_bonds, cyc, [(0, 1)], match="ring")
    # Order as written is kept (it decides which side `charge_on` means), duplicates dropped.
    assert B.select_bonds(geom(METHANOL), [(5, 1), (1, 5), (0, 2)])[0] == [(5, 1), (0, 2)]
    raises(ValueError, B.select_bonds, geom(METHANOL), "weakest")


def test_bad_bonds_and_disconnected_parents():
    g = geom(METHANOL)
    raises(B.BDEError, B.split_bond, g, 2, 5, match="not bonded")  # two H atoms
    raises(B.BDEError, B.split_bond, g, 1, 1)
    raises(B.BDEError, B.split_bond, g, 0, 99)
    two = TH.Geometry(["O", "H", "H", "O", "H", "H"],
                      [[0, 0, 0.1], [0, 0.76, -0.47], [0, -0.76, -0.47],
                       [5, 0, 0.1], [5, 0.76, -0.47], [5, -0.76, -0.47]])
    raises(B.BDEError, B.split_bond, two, 0, 1, match="one connected molecule")


def test_charges_and_spins_are_explicit():
    c, n, nh2 = ["C", "H", "H", "H"], ["N", "H", "H", "H"], ["N", "H", "H"]
    # Neutral parent (CH3-NH2): radicals, from the electron count.
    assert B.fragment_states(c, nh2) == ((0, 2), (0, 2))
    # Charged parent (CH3-NH3+) below.
    # A charged parent never guesses where the charge sits.
    raises(B.BDEError, B.fragment_states, c, n, parent_charge=1, match="charge_on")
    assert B.fragment_states(c, n, parent_charge=1, charge_on="j") == ((0, 2), (1, 2))  # CH3. + NH3.+
    # CH3+ and NH3 are both closed shells: the bond pair went to one side, which is not homolysis.
    raises(B.BDEError, B.fragment_states, c, n, parent_charge=1, charge_on="i", match="heterolytic")
    raises(ValueError, B.fragment_states, c, nh2, charge_on="k")
    # Heterolytic: needs the cation named, gives closed shells (CH3+ + NH2-).
    raises(B.BDEError, B.fragment_states, c, nh2, heterolytic=True, match="cation")
    assert B.fragment_states(c, nh2, heterolytic=True, cation="i") == ((1, 1), (-1, 1))
    assert B.fragment_states(c, nh2, heterolytic=True, cation="j") == ((-1, 1), (1, 1))
    raises(B.BDEError, B.fragment_states, c, nh2, cation="i", match="homolytic")
    # A parent with spin: the split of the spin is the caller's.
    raises(B.BDEError, B.fragment_states, c, nh2, parent_multiplicity=2, match="multiplicities")
    assert B.fragment_states(c, nh2, parent_multiplicity=3, multiplicities=(2, 2)) == ((0, 2), (0, 2))
    # Given multiplicities are checked against the electron count.
    raises(B.BDEError, B.fragment_states, c, nh2, multiplicities=(1, 2), match="electrons")
    # Atoms: H. is a doublet, H+ has no electrons and is a singlet, H- is a singlet.
    assert B.fragment_states(["O", "H"], ["H"]) == ((0, 2), (0, 2))
    assert B.fragment_states(["O", "H"], ["H"], heterolytic=True, cation="j") == ((-1, 1), (1, 1))


def test_species_key_and_dedupe():
    # Methane's four C-H bonds: one pair of fragments.
    g = geom(METHANE)
    plans = B.plan_bonds(g, B.select_bonds(g, "X-H")[0])
    assert len(plans) == 4
    assert len({p.side_i.key for p in plans}) == 1 and len({p.side_j.key for p in plans}) == 1
    assert plans[0].side_i.key.startswith("CH3|q+0|m2|") and plans[0].side_j.key.startswith("H|q+0|m2|")
    # Methanol: the three methyl C-H bonds are one CH2OH fragment, the O-H is CH3O.
    m = geom(METHANOL)
    plans = B.plan_bonds(m, B.select_bonds(m, "X-H")[0])
    keys = {s.key for p in plans for s in (p.side_i, p.side_j)}
    assert len(keys) == 3, keys
    ch2oh, ch3o = plans[0].side_i.key, plans[3].side_i.key
    assert ch2oh != ch3o and ch2oh.split("|")[0] == ch3o.split("|")[0] == "CH3O"  # same formula: the hash tells them apart

    # The hash does not care about atom order ...
    sym = ["C", "O", "H", "H", "H", "H"]
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5)]
    perm = [3, 5, 0, 4, 2, 1]  # new index of old atom k
    sym2 = [None] * 6
    for old, new in enumerate(perm):
        sym2[new] = sym[old]
    bonds2 = [(perm[a], perm[b]) for a, b in bonds]
    assert B.connectivity_hash(sym, bonds) == B.connectivity_hash(sym2, bonds2)
    # ... and tells constitutional isomers apart: butane against isobutane (heavy atoms + H).
    def alkane(skeleton):
        n = len(skeleton) + 1
        symbols = ["C"] * n
        bonds = list(skeleton)
        valence = [0] * n
        for a, b in skeleton:
            valence[a] += 1
            valence[b] += 1
        for c in range(n):
            for _ in range(4 - valence[c]):
                symbols.append("H")
                bonds.append((c, len(symbols) - 1))
        return symbols, bonds
    butane = alkane([(0, 1), (1, 2), (2, 3)])
    isobutane = alkane([(0, 1), (0, 2), (0, 3)])
    assert sorted(butane[0]) == sorted(isobutane[0])
    assert B.connectivity_hash(*butane) != B.connectivity_hash(*isobutane)
    # charge and multiplicity are part of the identity, and names are filename-safe.
    k1 = B.species_key(["C", "H", "H", "H"], [(0, 1), (0, 2), (0, 3)], 0, 2)
    k2 = B.species_key(["C", "H", "H", "H"], [(0, 1), (0, 2), (0, 3)], 1, 1)
    assert k1 != k2 and B.species_name(["C", "H", "H", "H"], 1, 1, k2).startswith("CH3p_m1_")
    assert B.hill_formula(["H", "O", "H"]) == "H2O" and B.hill_formula(["Cl", "C", "H"]) == "CHCl"


# ---------------------------------------------------------------------------
#  Thermochemistry
# ---------------------------------------------------------------------------


def test_hydrogen_atom_thermochemistry():
    # H(298) - E = 5/2 RT: three halves from translation and RT from pV.
    corr = TH.atom_correction("H", 2, T, 1.0, standard_state=False)
    close(corr["h_corr"], 2.5 * RT, 1e-12)
    close(2.5 * RT, 1.4812, 1e-4)
    close(corr["e_trans"], 1.5 * RT, 1e-12)

    # Translational entropy, Sackur-Tetrode, written out again from the constants.
    h, kb, na_mass = 6.62607015e-34, 1.380649e-23, 1.008 * 1.66053906660e-27
    q_per_n = (2 * math.pi * na_mass * kb * T / h**2) ** 1.5 * kb * T / 101325.0
    s_trans = 1.98720425864 * (math.log(q_per_n) + 2.5)
    close(s_trans, 26.0146, 1e-3)
    close(corr["s_trans_cal"], s_trans, 1e-9)
    # Independent check: NIST S(H, 1 bar) = 114.716 J/(mol K) includes R ln 2; at 1 bar the
    # translational part is 26.0146 + R ln(1.01325) = 26.0407 cal = 108.95 J.
    close((s_trans + 1.98720425864 * math.log(1.01325)) * 4.184 + 8.314462618 * math.log(2.0), 114.716, 0.02)
    # With the spin degeneracy: R ln 2 = 1.377, total 27.392 cal/(mol K).
    close(corr["s_elec_cal"], 1.98720425864 * math.log(2.0), 1e-12)
    close(corr["s_elec_cal"], 1.3774, 1e-4)
    close(corr["s_trans_cal"] + corr["s_elec_cal"], 27.392, 1e-3)
    close(corr["ts_total"], T * (s_trans + corr["s_elec_cal"]) / 1000.0, 1e-9)
    close(corr["g_corr"], corr["h_corr"] - corr["ts_total"], 1e-12)
    assert corr["zpe"] == 0.0 and corr["n_imaginary"] == 0 and corr["s_vib_cal"] == 0.0
    # A heavier atom has more translational entropy, by 3/2 R ln(m'/m).
    cl = TH.atom_correction("Cl", 2, T, 1.0, standard_state=False)
    close(cl["s_trans_cal"] - corr["s_trans_cal"], 1.5 * 1.98720425864 * math.log(35.45 / 1.008), 1e-9)
    # Pressure and the standard state move G and not H.
    p2 = TH.atom_correction("H", 2, T, 2.0, standard_state=False)
    close(p2["g_corr"] - corr["g_corr"], RT * math.log(2.0), 1e-9)
    ss = TH.atom_correction("H", 2, T, 1.0, standard_state=True)
    close(ss["g_corr"] - corr["g_corr"], TH.standard_state_correction(T), 1e-12)
    raises(ValueError, TH.atom_correction, "Xx", 2)
    raises(ValueError, TH.electronic_entropy_cal, 0)


def test_enthalpy_correction_and_electronic_entropy_in_molecules():
    freqs = [0.0] * 6 + [500.0, 1000.0, 3000.0]
    thermo = fake_thermo()
    plain = TH.free_energy_correction(freqs, 3, thermo, standard_state=False)
    zpe = 0.5 * TH.R_KCAL * TH.CM1_TO_KELVIN * 4500.0
    e_vib = sum(TH.R_KCAL * TH.CM1_TO_KELVIN * nu / math.expm1(TH.CM1_TO_KELVIN * nu / T) for nu in (500.0, 1000.0, 3000.0))
    close(plain["h_corr"], zpe + e_vib + 2 * 0.0014164 * TH.HARTREE_TO_KCAL + RT, 1e-9)
    close(plain["g_corr"], plain["h_corr"] - plain["ts_total"], 1e-12)
    assert plain["multiplicity"] == 1 and plain["s_elec_cal"] == 0.0

    # The Fortran block says "multiplicity 1, no electronic entropy" for a radical too; the
    # multiplicity given here is what counts, and it moves G and not H.
    radical = TH.free_energy_correction(freqs, 3, thermo, standard_state=False, multiplicity=2)
    close(radical["h_corr"], plain["h_corr"], 1e-12)
    close(radical["g_corr"] - plain["g_corr"], -T * 1.98720425864 * math.log(2.0) / 1000.0, 1e-9)
    triplet = TH.free_energy_correction(freqs, 3, thermo, standard_state=False, multiplicity=3)
    close(triplet["s_elec_cal"], 1.98720425864 * math.log(3.0), 1e-12)
    # Left unset, it takes what the block reports, as before.
    thermo2 = fake_thermo(multiplicity_reported=2)
    thermo2["contributions"]["electronic"]["entropy_cal_mol_K"] = 1.3774
    close(TH.free_energy_correction(freqs, 3, thermo2, standard_state=False)["s_elec_cal"], 1.3774, 1e-12)


# ---------------------------------------------------------------------------
#  Assembly
# ---------------------------------------------------------------------------


def sd(name, e, zpe, h, g, mult=1, **kw):
    return B.SpeciesData(name, name, name, 0, mult, kw.pop("n_atoms", 3), e, zpe, h, g, T, **kw)


def test_assemble_d0_and_bde_from_toy_numbers():
    # parent: E -100.000 Eh, ZPE 30.0, H = E + 33.0, G = E + 20.0 (kcal/mol on top of E)
    # A: E -99.600, ZPE 25.0, H = E + 27.5, G = E + 14.0
    # B (an atom): E -0.300, ZPE 0, H = E + 1.5, G = E - 4.0
    hk = B.HARTREE_TO_KCAL
    parent = sd("P", -100.0, 30.0, -100.0 * hk + 33.0, -100.0 * hk + 20.0)
    a = sd("A", -99.6, 25.0, -99.6 * hk + 27.5, -99.6 * hk + 14.0)
    b = sd("B", -0.3, 0.0, -0.3 * hk + 1.5, -0.3 * hk - 4.0, n_atoms=1)
    d_e, d0, bde, dg = B.assemble_bond(parent, a, b)
    close(d_e, 0.1 * hk, 1e-6)  # -99.6 - 0.3 + 100 = 0.1 Eh
    close(d_e, 62.7509474, 1e-5)
    close(d0, d_e + (25.0 + 0.0 - 30.0), 1e-9)  # dZPE = -5.0
    close(d0, 57.7509474, 1e-5)
    close(bde, d_e + (27.5 + 1.5 - 33.0), 1e-6)  # dH_corr = -4.0
    close(bde, 62.7509474 - 4.0, 1e-5)
    close(dg, d_e + (14.0 - 4.0 - 20.0), 1e-6)
    close(dg, 52.7509474, 1e-5)
    # D0 < D_e (lost ZPE), and kJ/mol is 4.184 times it.
    row = B.BondRow(0, 1, "C-H", [0], [1], "A", "B", d_e, d0, bde, dg, T)
    close(row.bde_kj, bde * 4.184, 1e-12)
    assert row.to_dict()["d0_kj"] == row.d0_kj
    # Different temperatures cannot be combined.
    hot = B.SpeciesData("B", "B", "B", 0, 1, 1, -0.3, 0.0, 0.0, 0.0, 350.0)
    raises(B.BDEError, B.assemble_bond, parent, a, hot, match="temperature")


def test_result_table_json_and_warnings():
    hk = B.HARTREE_TO_KCAL
    parent = sd("CH4O_parent", -100.0, 30.0, -100.0 * hk + 33.0, -100.0 * hk + 20.0, n_atoms=6)
    a = sd("CH3O_m2_abc", -99.6, 25.0, -99.6 * hk + 27.5, -99.6 * hk + 14.0, mult=2, n_atoms=5,
           e_vertical_hartree=-99.597, relaxation_kcal=1.88, n_imaginary=1)
    b = sd("H_m2_def", -0.3, 0.0, -0.3 * hk + 1.5, -0.3 * hk - 4.0, mult=2, n_atoms=1)
    weak = sd("CH3O_m2_xyz", -99.7, 25.0, -99.7 * hk + 27.5, -99.7 * hk + 14.0, mult=2, n_atoms=5)

    def row(i, j, x):
        d_e, d0, bde, dg = B.assemble_bond(parent, x, b)
        return B.BondRow(i, j, "O-H", [i], [j], x.name, b.name, d_e, d0, bde, dg, T)

    result = B.BDEResult(parent, {a.name: a, b.name: b, weak.name: weak}, [row(1, 5, a), row(0, 2, weak)], T,
                         protocol={"qrrho": True}, skipped=[(0, 1, "ring")], notes=["a note"])
    assert [r.i for r in result.bonds] == [0, 1]  # sorted: the weaker bond first
    assert result.bonds[0].bde_kcal < result.bonds[1].bde_kcal
    text = result.table()
    assert "BDE(298)" in text and "CH4O" in text and "O-H" in text and "skipped" in text
    assert f"{result.bonds[0].bde_kcal:.1f}" in text and f"{result.bonds[0].bde_kj:.1f}" in text
    warnings = "\n".join(result.warnings)
    assert "a note" in warnings and "1 imaginary" in warnings and "CH3O_m2_abc" in warnings
    document = json.loads(result.to_json())
    assert document["temperature_K"] == T and len(document["bonds"]) == 2
    assert document["species"]["CH3O_m2_abc"]["relaxation_kcal"] == 1.88
    assert document["bonds"][0]["bde_kj"] == result.bonds[0].bde_kj
    assert document["parent"]["name"] == "CH4O_parent" and document["skipped"][0]["reason"] == "ring"
    path = os.path.join(tempfile.mkdtemp(), "r.json")
    try:
        result.to_json(path)
        assert json.load(open(path))["protocol"] == {"qrrho": True}
    finally:
        shutil.rmtree(os.path.dirname(path))
    # A negative bond energy is flagged.
    neg = B.BDEResult(parent, {a.name: a, b.name: b}, [B.BondRow(0, 1, "C-H", [0], [1], a.name, b.name, -1, -1, -2.0, -3, T)], T)
    assert any("negative" in w for w in neg.warnings)


# ---------------------------------------------------------------------------
#  Orchestration, against a stand-in for the library
# ---------------------------------------------------------------------------

#: Electronic energy of each species, picked out by what it is (see `physics`).
ENERGY = {"parent": -40.00, "CH3O": -39.55, "CH2OH": -39.56, "H": -0.40, "CH3": -39.8, "CH4": -40.5}


def _o_has_h(system):
    o = system.symbols.index("O")
    return any(
        s == "H" and math.dist(system.coords[o], c) < 1.2 for s, c in zip(system.symbols, system.coords)
    )


def physics(driver, system, label, n):
    n_atoms = system.n_atoms
    if "O" not in system.symbols:  # methane and its family
        e = {5: ENERGY["CH4"], 4: ENERGY["CH3"], 1: ENERGY["H"]}[n_atoms]
    elif n_atoms == 6:
        e = ENERGY["parent"]
    elif n_atoms == 5:
        e = ENERGY["CH2OH"] if _o_has_h(system) and system.symbols.count("H") == 3 and _c_has_two_h(system) else ENERGY["CH3O"]
    elif n_atoms == 1:
        e = ENERGY["H"]
    else:
        raise AssertionError(f"unexpected species of {n_atoms} atoms")
    if "_vert" in label:
        e += 0.003  # the geometry cut out of the parent is 1.88 kcal/mol above its minimum
    if driver == "hessian":
        assert n_atoms > 1, "a Hessian was asked for one atom"
        freqs = [0.0, 0.1, -0.1, 0.2, 0.3, -0.2] + [400.0 + 100.0 * i for i in range(3 * n_atoms - 6)]
        thermo = fake_thermo()
        thermo["total_energies_hartree"] = {"electronic": e}
        return FakeResult(e, label, freqs, thermo)
    return FakeResult(e, label)


def _c_has_two_h(system):
    c = system.symbols.index("C")
    return sum(1 for s, x in zip(system.symbols, system.coords) if s == "H" and math.dist(system.coords[c], x) < 1.3) == 2


def _expected_correction(n_atoms, multiplicity):
    if n_atoms == 1:
        return TH.atom_correction("H", multiplicity, T, 1.0, standard_state=False)
    freqs = [0.0, 0.1, -0.1, 0.2, 0.3, -0.2] + [400.0 + 100.0 * i for i in range(3 * n_atoms - 6)]
    return TH.free_energy_correction(freqs, n_atoms, fake_thermo(), standard_state=False, multiplicity=multiplicity)


def _fake_executable(calls, opt_energy):
    """Stand-in for `_run_mqc`: an optimizer that returns its input, and a CREST that finds one structure."""

    def run(executable, deck, directory, env):
        calls.append((deck, directory))
        with open(os.path.join(directory, deck)) as handle:
            document = json.load(handle)
        with open(os.path.join(directory, "start.xyz")) as handle:
            xyz = handle.read().split("\n", 2)
        n, body = xyz[0], xyz[2]
        if document["driver"] == "conformers":
            with open(os.path.join(directory, "crest_conformers.xyz"), "w") as out:
                out.write(f"{n}\n   {opt_energy:.8f}\n{body}")
        else:
            with open(os.path.join(directory, f"output_{deck[:-5]}_optimized.xyz"), "w") as out:
                out.write(f"{n}\nmetalquicha converged, E = {opt_energy:20.12f} Hartree\n{body}")

    return run


def test_scan_methanol_runs_each_species_once_and_assembles_bde():
    log = []
    workdir = tempfile.mkdtemp()
    try:
        with fake_package(log, physics, module="bde") as bde:
            protocol = bde.Protocol(conformers=None, optimize=None, workdir=workdir)
            result = bde.scan(METHANOL, "X-H", protocol, verbose=False)

            # One Hessian for the parent and each radical, no Hessian and one energy for H.
            runs = [(e["driver"], e["system"].n_atoms, e["system"].multiplicity) for e in log]
            assert runs == [("hessian", 6, 1), ("hessian", 5, 2), ("energy", 1, 2), ("hessian", 5, 2)], runs
            assert [e["system"].charges for e in log] == [[0]] * 4
            assert len(result.bonds) == 4 and len(result.species) == 3  # CH2OH, H, CH3O

            # Sorted: the C-H rows (CH2OH, lower energy fragment) come before the O-H.
            labels = [(r.label, r.i, r.j) for r in result.bonds]
            assert labels[-1] == ("O-H", 1, 5) and all(lab[0] == "C-H" for lab in labels[:-1])
            hk = TH.HARTREE_TO_KCAL

            parent = _expected_correction(6, 1)
            h_atom = _expected_correction(1, 2)
            rad = _expected_correction(5, 2)
            e_p = ENERGY["parent"]
            for r, e_frag in ((result.bonds[0], ENERGY["CH2OH"]), (result.bonds[-1], ENERGY["CH3O"])):
                d_e = (e_frag + ENERGY["H"] - e_p) * hk
                close(r.d_e_kcal, d_e, 1e-6)
                close(r.d0_kcal, d_e + rad["zpe"] + h_atom["zpe"] - parent["zpe"], 1e-6)
                close(r.bde_kcal, d_e + rad["h_corr"] + h_atom["h_corr"] - parent["h_corr"], 1e-6)
                close(r.dg_kcal, d_e + rad["g_corr"] + h_atom["g_corr"] - parent["g_corr"], 1e-6)
                close(r.bde_kj, r.bde_kcal * 4.184, 1e-9)
            # dG carries the electronic entropy of two doublets, BDE(298) does not: the radical's
            # free energy is lower by T R ln 2 than the closed-shell arithmetic would give.
            rad_closed = TH.free_energy_correction([0.0, 0.1, -0.1, 0.2, 0.3, -0.2] + [400.0 + 100.0 * i for i in range(9)], 5, fake_thermo(), standard_state=False)
            close(rad["g_corr"] - rad_closed["g_corr"], -RT * math.log(2.0), 1e-9)

            # The atom: analytic thermochemistry, no frequencies, and nothing imaginary.
            atom = next(s for s in result.species.values() if s.n_atoms == 1)
            close(atom.h_kcal, ENERGY["H"] * hk + 2.5 * RT, 1e-6)
            assert atom.n_imaginary == 0 and atom.detail["atom"] is True and atom.multiplicity == 2
            # The parent keeps its own, and every species is at 298.15 K.
            assert result.parent.name == "CH4O_parent" and result.temperature == T
            warnings = "\n".join(result.warnings)
            assert "open-shell" in warnings and "spin-orbit" in warnings and "S^2" in warnings
            text = result.table()
            assert "O-H" in text and "C-H" in text
            json.loads(result.to_json())
            assert result.protocol["standard_state"] is False  # a gas-phase default
    finally:
        shutil.rmtree(workdir)


def test_default_protocol_is_gas_phase_and_unsolvated():
    log = []
    with fake_package(log, physics, module="bde") as bde:
        p = bde.Protocol()
        for stage in (p.conformers, p.optimize, p.frequencies, p.single_point):
            assert stage == {"method": "gfn2", "verbosity": "error"}, stage  # no "xtb": {"solvent": ...}
        assert p.standard_state is False and p.fragment_conformers is False and p.parent_conformers is True
        assert p.workdir == "bde_work" and p.single_point_is_frequencies()
        assert p.conditions() == (298.15, 1.0)
        p = bde.Protocol(frequencies={"method": "gfn2", "hessian": {"temperature": 350.0, "pressure": 2.0}})
        assert p.conditions() == (350.0, 2.0)
        raises(ValueError, bde.Protocol, frequencies=None)
        raises(ValueError, bde.Protocol, optimize={"method": "gfn2", "driver": "optimize"})
    with fake_package(log, physics, module="pka") as pka:  # and the pKa default stays aqueous
        assert pka.Protocol().optimize["xtb"]["solvent"] == "water" and pka.Protocol().standard_state is True


def test_scan_stages_vertical_energy_and_conformers_only_for_the_parent():
    log, calls = [], []
    workdir = tempfile.mkdtemp()
    try:
        with fake_package(log, physics, module="bde") as bde:
            bde._run_mqc = _fake_executable(calls, -40.0)
            bde._find_executable = lambda protocol: "/fake/mqc"
            protocol = bde.Protocol(workdir=workdir)  # conformers on for the parent, off for fragments
            result = bde.scan(METHANE, "X-H", protocol, verbose=False)

            stages = [deck for deck, _ in calls]
            assert stages[:2] == ["conformers.json", "optimize.json"]  # the parent: search, then optimize
            assert stages.count("conformers.json") == 1  # fragments are not searched by default
            assert stages.count("optimize.json") == 2  # the parent and CH3; an atom needs none
            with open(os.path.join(calls[0][1], "conformers.json")) as handle:
                deck = json.load(handle)
            assert "xtb" not in deck["keywords"] and deck["model"]["method"] == "gfn2"  # gas phase
            assert deck["molecules"][0]["molecular_multiplicity"] == 1

            # Methane's four bonds are one calculation per species: parent, CH3, H.
            assert len(result.bonds) == 4 and len(result.species) == 2
            drivers = [(e["driver"], e["system"].n_atoms, e["label"]) for e in log]
            assert [d[:2] for d in drivers if "_vert" not in d[2]] == [("hessian", 5), ("hessian", 4), ("energy", 1)]
            vertical = [d for d in drivers if "_vert" in d[2]]
            assert len(vertical) == 1 and vertical[0][1] == 4  # the CH3 fragment; not the parent, not the atom
            ch3 = next(s for s in result.species.values() if s.n_atoms == 4)
            close(ch3.relaxation_kcal, 0.003 * TH.HARTREE_TO_KCAL, 1e-6)
            close(ch3.e_vertical_hartree - ch3.e_hartree, 0.003, 1e-12)
            h = next(s for s in result.species.values() if s.n_atoms == 1)
            assert h.relaxation_kcal in (None, 0.0) or abs(h.relaxation_kcal) < 1e-12
            assert result.parent.relaxation_kcal is None
            assert max(r.bde_kcal for r in result.bonds) - min(r.bde_kcal for r in result.bonds) < 1e-9
    finally:
        shutil.rmtree(workdir)


def test_scan_refuses_before_running_anything():
    log = []
    with fake_package(log, physics, module="bde") as bde:
        protocol = bde.Protocol(conformers=None, optimize=None)
        cyc = cyclohexane()
        raises(bde.BDEError, bde.scan, cyc, [(0, 1)], protocol, match="ring")
        raises(bde.BDEError, bde.scan, METHANOL, [(2, 5)], protocol, match="not bonded")
        raises(bde.BDEError, bde.scan, METHANOL, "X-H", protocol, charge=1, multiplicity=2, match="charge_on")
        raises(bde.BDEError, bde.scan, METHANOL, "X-H", protocol, charge=1, multiplicity=2, charge_on="i", match="multiplicities")
        raises(bde.BDEError, bde.scan, METHANOL, [(1, 5)], protocol, heterolytic=True, match="cation")
        raises(bde.BDEError, bde.scan, METHANOL, [(1, 5)], protocol, charge=1, multiplicity=1, match="electrons")
        raises(TypeError, bde.scan, object(), "X-H", protocol)
        raises(ValueError, bde.scan, METHANOL, "weakest", protocol)
        assert log == []
        # A selection that finds nothing is an error, not an empty table.
        raises(bde.BDEError, bde.scan, (["He", "He"], [[0, 0, 0], [0, 0, 3.0]]), "all", protocol)
        # A ring selected by keyword is skipped and reported, and never run.
        plan = bde.split(cyc, 0, 6)
        assert plan.side_j.symbols == ["H"] and plan.side_i.multiplicity == 2
        bonds, skipped = bde._bde_math.select_bonds(TH.Geometry(*cyc), "all")
        assert len(skipped) == 6 and len(bonds) == 12


def test_heterolytic_and_charged_plans():
    log = []
    with fake_package(log, physics, module="bde") as bde:
        plan = bde.split(METHANOL, 1, 5, heterolytic=True, cation="j")
        assert (plan.side_i.charge, plan.side_i.multiplicity) == (-1, 1)  # CH3O-
        assert (plan.side_j.charge, plan.side_j.multiplicity) == (1, 1)  # H+
        assert plan.heterolytic
        # An odd-electron "anion" of methanol is not a singlet, and the plan says so before anything runs.
        raises(bde.BDEError, bde.split, METHANOL, 1, 5, charge=-1, heterolytic=True, cation="j", match="electrons")
        raises(bde.BDEError, bde.split, METHANOL, 1, 5, charge=-1, multiplicity=2, heterolytic=True, cation="j",
               charge_on="i", match="multiplicities")


def test_one_atom_is_never_given_a_hessian_in_pka_either():
    # The shared evaluation sends any one-atom species down the atom path.
    log = []
    workdir = tempfile.mkdtemp()
    try:
        with fake_package(log, physics, module="pka") as pka:
            protocol = pka.Protocol(conformers=None, optimize=None, workdir=workdir)
            chloride = pka.Microstate("Cl-", pka.Geometry(["Cl"], [[0, 0, 0]]), charge=-1, n_protons=0)
            data = pka.evaluate_microstate(chloride, protocol, verbose=False)
            assert [e["driver"] for e in log] == ["energy"]
            close(data.free_energy, ENERGY["H"] * TH.HARTREE_TO_KCAL
                  + TH.atom_correction("Cl", 1, T, 1.0, standard_state=True)["g_corr"], 1e-6)
    finally:
        shutil.rmtree(workdir)


TESTS = [(name, fn) for name, fn in sorted(globals().items()) if name.startswith("test_") and callable(fn)]


def main():
    failures = 0
    for name, fn in TESTS:
        try:
            fn()
        except Exception as error:  # noqa: BLE001 -- report every failure, then exit non-zero
            failures += 1
            print(f"FAIL {name}: {type(error).__name__}: {error}")
        else:
            print(f"ok   {name}")
    print(f"\n{len(TESTS) - failures} passed, {failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
