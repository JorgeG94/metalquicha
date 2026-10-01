"""Tests for the pKa workflow that need neither the library nor MPI nor a build.

    python3 python/tests/test_pka.py

Exits non-zero on the first failure; also collected by pytest, if that is what
is to hand. Nothing here imports `mqc` the ordinary way, because that loads
`libmqc`. Two things are done instead:

  * the arithmetic, `mqc/_pka_math.py`, is loaded by path -- it imports nothing
    from the package, which is the reason it is a separate file;
  * `mqc/pka.py` is imported under the stand-in package of `_standin.py`, whose
    `MBE` and `_check_label` are the *real* ones, lifted out of `mqc/__init__.py`
    by parsing it, with only the calls that reach Fortran replaced. So the decks
    the workflow writes, the settings it passes and the labels it chooses are
    the real thing, and the physics is invented.

What this cannot test is anything on the other side of the C interface: that
the Hessian run really writes the keys read here, that CREST leaves the file
named here, that the optimizer's output is where it is looked for. Those are
read from the Fortran and the JSON writer, and are checked by running
`python/examples/pka.py` against a real build.
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

from _standin import fake_package  # noqa: E402  (the stand-in `mqc`, shared with test_bde.py)


def _load_math():
    spec = importlib.util.spec_from_file_location("mqc_pka_math_under_test", os.path.join(PKG, "_pka_math.py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


M = _load_math()
LN10 = math.log(10.0)
T = 298.15
RT = M.R_KCAL * T


def close(a, b, tol=1e-9, msg=""):
    assert abs(a - b) <= tol, f"{msg} {a!r} != {b!r} (|diff| {abs(a - b):.3e} > {tol:.1e})"


def raises(exc, fn, *args, **kwargs):
    try:
        fn(*args, **kwargs)
    except exc as error:
        return error
    raise AssertionError(f"{exc.__name__} not raised")


# ---------------------------------------------------------------------------
#  quasi-RRHO and the standard state
# ---------------------------------------------------------------------------


def test_qrrho_hand_values():
    # Worked by hand at 100 cm^-1 and 298.15 K, where the weight is exactly 1/2.
    #   u = 1.4387773538 * 100 / 298.15 = 0.482584
    #   S_RRHO / R = u/(e^u - 1) - ln(1 - e^-u) = 0.777817 + 0.960463 = 1.738280
    #   S_FR / R   = 1/2 + ln sqrt(8 pi^3 mu' kT / h^2), mu = h/(8 pi^2 nu),
    #                mu' = mu B/(mu + B), B = 1e-44 kg m^2  ->  1.436542
    close(M.qrrho_weight(100.0), 0.5, 1e-15)
    close(M.rrho_entropy_over_r(100.0, T), 1.738280, 2e-6)
    close(M._s_free_rotor_over_r(100.0, T), 1.436542, 2e-6)
    close(M.qrrho_entropy_over_r(100.0, T), 0.5 * (1.738280 + 1.436542), 2e-6)

    # The weight's other values: 1/(1 + (100/nu)^4).
    close(M.qrrho_weight(50.0), 1.0 / 17.0, 1e-15)
    close(M.qrrho_weight(200.0), 16.0 / 17.0, 1e-15)
    assert M.qrrho_weight(-5.0) == 0.0


def test_qrrho_free_rotor_independent_derivation():
    # The free-rotor entropy again, written out from the constants rather than
    # through the module's helper, at a frequency where it is not the limit.
    h, kb, c = 6.62607015e-34, 1.380649e-23, 299792458.0
    nu = 37.0
    mu = h / (8 * math.pi**2 * c * 100.0 * nu)
    mup = mu * 1e-44 / (mu + 1e-44)
    expected = 0.5 + 0.5 * math.log(8 * math.pi**3 * mup * kb * T / h**2)
    close(M._s_free_rotor_over_r(37.0, T), expected, 1e-12)


def test_qrrho_limits():
    # nu -> infinity: the oscillator, exactly what the Fortran RRHO sums.
    for nu in (1000.0, 2000.0, 3500.0):
        close(M.qrrho_entropy_over_r(nu, T), M.rrho_entropy_over_r(nu, T), 1e-4)
    close(M.qrrho_entropy_over_r(1.0e5, T), 0.0, 1e-9)

    # nu -> 0: the free rotor, which is *finite* where the oscillator diverges.
    assert M.rrho_entropy_over_r(1.0e-3, T) > 10.0
    h, kb = 6.62607015e-34, 1.380649e-23
    s_limit = 0.5 + 0.5 * math.log(8 * math.pi**3 * 1e-44 * kb * T / h**2)  # mu' -> B
    close(M.qrrho_entropy_over_r(1.0e-6, T), s_limit, 1e-3)
    s = [M.qrrho_entropy_over_r(nu, T) for nu in (1.0, 0.1, 0.01, 0.001)]
    assert all(b >= a - 1e-9 for a, b in zip(s, s[1:])) and s[-1] < s_limit, s  # rises to the cap, never past it

    raises(ValueError, M.qrrho_entropy_over_r, 0.0, T)


def test_standard_state_correction():
    # RT ln(24.465) at 298.15 K: the molar volume of an ideal gas at 1 atm, in
    # litres, against the 1 L of a molar solution.
    close(M.standard_state_correction(298.15), RT * math.log(24.4655), 2e-4)
    close(M.standard_state_correction(298.15), 1.89, 0.01)
    # Not a constant: it moves with T, and with the pressure the entropy used.
    close(M.standard_state_correction(350.0), M.R_KCAL * 350.0 * math.log(0.082057366 * 350.0), 1e-9)
    close(
        M.standard_state_correction(298.15, 2.0) - M.standard_state_correction(298.15, 1.0),
        -RT * math.log(2.0),
        1e-12,
    )


def test_default_proton_free_energy():
    # -6.28 + 1.89 - 265.9: the standard literature value, about -270.3.
    close(M.default_proton_free_energy(298.15), -270.29, 0.01)


# ---------------------------------------------------------------------------
#  Frequencies and the correction
# ---------------------------------------------------------------------------


WATER_FREQS = [-0.4, 0.1, 0.3, 1.2, -2.5, 4.0, 1650.0, 3800.0, 3900.0]


def test_split_frequencies():
    real, imag, tr_max = M.split_frequencies(WATER_FREQS, 3)
    assert real == [1650.0, 3800.0, 3900.0] and imag == [] and tr_max == 4.0

    # An imaginary vibration is not a translation or rotation, however small:
    # the six smallest magnitudes go, and a negative one among the rest is counted.
    freqs = [0.0, 0.1, -0.2, 0.3, 0.2, -0.1, -45.0, 1650.0, 3800.0]
    real, imag, tr_max = M.split_frequencies(freqs, 3)
    assert imag == [-45.0] and real == [1650.0, 3800.0] and tr_max == 0.3

    # Linear: five rotations/translations, so one more vibration.
    real, imag, _ = M.split_frequencies([0.0, 0.0, 0.1, -0.1, 0.2, 700.0, 2300.0, 2300.0, 4000.0][:9], 3, is_linear=True)
    assert len(real) == 4

    raises(M.PKaError, M.split_frequencies, [1.0] * 6, 3)  # not 3N: refuse, do not guess


def _thermo(temperature=T, pressure=1.0):
    return {
        "temperature_K": temperature,
        "pressure_atm": pressure,
        "is_linear": False,
        "contributions": {
            "translational": {"energy_hartree": 0.0014164, "entropy_cal_mol_K": 34.6},
            "rotational": {"energy_hartree": 0.0014164, "entropy_cal_mol_K": 10.5},
            "vibrational": {"energy_hartree": 0.0, "entropy_cal_mol_K": 0.0},
            "electronic": {"energy_hartree": 0.0, "entropy_cal_mol_K": 0.0},
        },
    }


def test_free_energy_correction_by_hand():
    # One real mode at 500 cm^-1 (three atoms, six zero modes) so every term can
    # be written out. The mode is above 100 cm^-1, w = 0.998402...
    freqs = [0.0] * 6 + [500.0, 1000.0, 3000.0]
    thermo = _thermo()
    corr = M.free_energy_correction(freqs, 3, thermo)
    zpe = 0.5 * M.R_KCAL * M.CM1_TO_KELVIN * (500.0 + 1000.0 + 3000.0)
    close(corr["zpe"], zpe, 1e-12)
    evib = sum(M.R_KCAL * M.CM1_TO_KELVIN * nu / math.expm1(M.CM1_TO_KELVIN * nu / T) for nu in (500.0, 1000.0, 3000.0))
    close(corr["e_vib"], evib, 1e-12)
    svib = M.R_CAL * sum(M.qrrho_entropy_over_r(nu, T) for nu in (500.0, 1000.0, 3000.0))
    close(corr["s_vib_cal"], svib, 1e-12)
    ts = T * (svib + 34.6 + 10.5) / 1000.0
    close(corr["ts_total"], ts, 1e-12)
    gas = zpe + evib + 2 * 0.0014164 * M.HARTREE_TO_KCAL + RT - ts
    close(corr["g_corr"], gas + M.standard_state_correction(T), 1e-9)
    assert corr["n_imaginary"] == 0 and corr["n_real"] == 3

    # Standard state off and plain RRHO: both are what they say.
    no_ss = M.free_energy_correction(freqs, 3, thermo, standard_state=False)
    close(corr["g_corr"] - no_ss["g_corr"], M.standard_state_correction(T), 1e-9)
    rrho = M.free_energy_correction(freqs, 3, thermo, qrrho=False)
    assert rrho["s_vib_cal"] > 0.0 and abs(rrho["g_corr"] - corr["g_corr"]) < 0.05


def test_soft_mode_costs_less_under_qrrho():
    # The reason for the correction: a 20 cm^-1 mode has a harmonic entropy far
    # above what a hindered rotor can have, so RRHO gives too low a free energy.
    freqs = [0.0] * 6 + [20.0, 1000.0, 3000.0]
    q = M.free_energy_correction(freqs, 3, _thermo())
    r = M.free_energy_correction(freqs, 3, _thermo(), qrrho=False)
    assert q["g_corr"] > r["g_corr"] + 0.3, (q["g_corr"], r["g_corr"])


def test_imaginary_modes_are_counted_never_dropped_silently():
    freqs = [0.0] * 6 + [-60.0, 1000.0, 3000.0]
    drop = M.free_energy_correction(freqs, 3, _thermo(), imaginary="drop")
    flip = M.free_energy_correction(freqs, 3, _thermo(), imaginary="flip")
    for corr in (drop, flip):
        assert corr["n_imaginary"] == 1 and corr["imaginary_cm1"] == [-60.0]
    assert drop["n_real"] == 2 and flip["g_corr"] != drop["g_corr"]
    raises(ValueError, M.free_energy_correction, freqs, 3, _thermo(), imaginary="ignore")


# ---------------------------------------------------------------------------
#  Boltzmann
# ---------------------------------------------------------------------------


def test_boltzmann_combination():
    # Two degenerate conformers: RT ln 2 lower.
    close(M.boltzmann_combine([10.0, 10.0], T), 10.0 - RT * math.log(2.0), 1e-12)
    # Far apart: the lowest.
    close(M.boltzmann_combine([0.0, 20.0], T), 0.0, 1e-8)
    # Weights normalise and order correctly.
    w = M.boltzmann_weights([0.0, 1.0, 2.0], T)
    close(sum(w), 1.0, 1e-15)
    assert w[0] > w[1] > w[2]
    close(w[0] / w[1], math.exp(1.0 / RT), 1e-9)
    # Absolute free energies are ~1e5 kcal/mol; the naive exp overflows.
    close(M.boltzmann_combine([-123456.0, -123456.0 + 0.5], T), -123456.0 - RT * math.log(1.0 + math.exp(-0.5 / RT)), 1e-8)
    raises(ValueError, M.boltzmann_combine, [], T)


# ---------------------------------------------------------------------------
#  The protonation model
# ---------------------------------------------------------------------------


def ms(name, q, n, g, cls=None):
    return M.MicrostateData(name, q, n, g, cls)


def model(states, **kwargs):
    return M.PKaModel(states, T, **kwargs)


def glycine(pka1, pka2, with_tautomer=False, shifts=None, ratio_kcal=None):
    """Cation / zwitterion / anion, built to have the pKas given (g_h = 0).

    With ``with_tautomer`` the neutral form is added at ``ratio_kcal`` above the
    zwitterion, which changes neither macroscopic pKa's *definition* nor, the
    pI -- the test below asks that it does not.
    """
    g_zw = 0.0
    g_cat = g_zw - pka1 * RT * LN10
    g_an = g_zw + pka2 * RT * LN10
    states = [
        ms("cation", +1, 2, g_cat, ("carboxylic", "ammonium")),
        ms("zwitterion", 0, 1, g_zw, ("ammonium",)),
        ms("anion", -1, 0, g_an, ()),
    ]
    if with_tautomer:
        states.append(ms("neutral", 0, 1, g_zw + ratio_kcal, ("carboxylic",)))
    return model(states, proton_free_energy=0.0, class_free_energy=shifts)


def test_two_independent_sites():
    # Two identical, independent sites of micro pKa p. Statistics move the
    # macroscopic constants: pKa1 = p - log10 2, pKa2 = p + log10 2.
    p = 6.3
    g2, g1, g0 = 0.0, p * RT * LN10, 2.0 * p * RT * LN10
    states = [
        ms("AH2", +2, 2, g2),
        ms("AH-a", +1, 1, g1),
        ms("AH-b", +1, 1, g1),
        ms("A2-", 0, 0, g0),
    ]
    mdl = model(states, proton_free_energy=0.0)
    macro = mdl.macro_pkas()
    assert [(d["from_n"], d["to_n"]) for d in macro] == [(2, 1), (1, 0)]
    close(macro[0]["pKa"], p - math.log10(2.0), 1e-9)
    close(macro[1]["pKa"], p + math.log10(2.0), 1e-9)
    for d in mdl.micro_pkas():
        close(d["pKa"], p, 1e-9)
    assert len(mdl.micro_pkas()) == 4  # AH2 -> a, b; a, b -> A2-


def test_populations_are_the_micro_pka():
    mdl = glycine(2.3, 9.6, with_tautomer=True, ratio_kcal=1.5)
    for d in mdl.micro_pkas():
        # At pH = pKa of s -> t the two are equally populated, by construction
        # of one set of populations -- the consistency the shifts must keep.
        pops = mdl.populations(d["pKa"])
        close(pops[d["acid"]], pops[d["base"]], 1e-9, d["acid"] + d["base"])
        close(sum(pops.values()), 1.0, 1e-12)


def test_glycine_pi_is_mean_of_pkas():
    mdl = glycine(2.34, 9.60)
    macro = mdl.macro_pkas()
    close(macro[0]["pKa"], 2.34, 1e-9)
    close(macro[1]["pKa"], 9.60, 1e-9)
    close(mdl.isoelectric_point(), 0.5 * (2.34 + 9.60), 1e-6)
    close(mdl.charge(mdl.isoelectric_point()), 0.0, 1e-6)

    # A neutral tautomer sharing the level with the zwitterion changes the
    # macroscopic pKas (they are partition sums) but pI stays their mean.
    mdl = glycine(2.34, 9.60, with_tautomer=True, ratio_kcal=0.8)
    macro = mdl.macro_pkas()
    assert macro[0]["pKa"] < 2.34 and macro[1]["pKa"] > 9.60 - 1.0
    close(mdl.isoelectric_point(), 0.5 * (macro[0]["pKa"] + macro[1]["pKa"]), 1e-6)


def test_class_shifts_keep_cycles_consistent():
    # With shifts on the carboxylic and ammonium sites, going cation -> anion
    # by either neutral form must cost the same: the shift belongs to the
    # protons a microstate carries, so it telescopes.
    shifts = {"carboxylic": 1.5, "ammonium": -2.0}
    mdl = glycine(2.0, 9.0, with_tautomer=True, ratio_kcal=1.0, shifts=shifts)
    by = {(d["acid"], d["base"]): d["pKa"] for d in mdl.micro_pkas()}
    via_neutral = by[("cation", "neutral")] + by[("neutral", "anion")]
    via_zwit = by[("cation", "zwitterion")] + by[("zwitterion", "anion")]
    close(via_neutral, via_zwit, 1e-9)
    # And the shift is the class of the proton that leaves.
    classes = {(d["acid"], d["base"]): d["site_class"] for d in mdl.micro_pkas()}
    assert classes[("cation", "zwitterion")] == "carboxylic"
    assert classes[("cation", "neutral")] == "ammonium"
    # populations still agree with the pKas under shifts.
    for d in mdl.micro_pkas():
        pops = mdl.populations(d["pKa"])
        close(pops[d["acid"]], pops[d["base"]], 1e-9)


def test_charge_is_monotonic_in_ph():
    cases = [
        glycine(2.34, 9.60),
        glycine(2.34, 9.60, with_tautomer=True, ratio_kcal=0.3),
        glycine(2.0, 9.0, with_tautomer=True, ratio_kcal=1.0, shifts={"carboxylic": 1.5, "ammonium": -2.0}),
    ]
    for mdl in cases:
        charges = [mdl.charge(0.05 * i) for i in range(0, 281)]
        assert all(b <= a + 1e-12 for a, b in zip(charges, charges[1:])), "charge increased with pH"
        assert charges[0] > 0.9 and charges[-1] < -0.9


def test_no_crossing_is_reported_not_an_edge():
    # Every microstate cationic: the charge is positive at every pH.
    cationic = model([ms("AH+", 2, 1, 0.0), ms("A", 1, 0, 3.0)], proton_free_energy=0.0)
    err = raises(M.NoIsoelectricPoint, cationic.isoelectric_point)
    assert "no sign change" in str(err)
    result = M.PKaResult([ms("AH+", 2, 1, 0.0), ms("A", 1, 0, 3.0)], T)
    assert result.pI is None and "no sign change" in result.pI_note

    # pKas that put the crossing below pH 0: the charge is already negative at 0.
    low = glycine(-3.0, -2.0)
    raises(M.NoIsoelectricPoint, low.isoelectric_point)
    high = glycine(15.0, 16.0)
    raises(M.NoIsoelectricPoint, high.isoelectric_point)


def test_charge_must_track_protons():
    raises(ValueError, model, [ms("a", 0, 1, 0.0), ms("b", 0, 0, 1.0)])
    raises(ValueError, model, [ms("a", 0, 1, 0.0), ms("a", -1, 0, 1.0)])  # duplicate name
    raises(ValueError, M.normalize_site_classes, ("x",), 2)
    assert M.normalize_site_classes("x", 2) == ("x", "x")
    assert M.normalize_site_classes(None, 1) == ("default",)


def test_calibration_recovers_known_shift():
    # Truth: shifted proton free energies per class. The "calculated" free
    # energies are what the uncalibrated method would give; the "experimental"
    # pKas are the truth model's micro pKas. Calibrating must recover the shifts.
    uncal = glycine(1.0, 8.0, with_tautomer=True, ratio_kcal=1.0)
    g_h = M.default_proton_free_energy(T)
    truth = {"carboxylic": g_h + 3.0, "ammonium": g_h - 2.0}
    free = {m.name: m.free_energy for m in uncal.microstates}
    true_model = M.PKaModel(
        [M.MicrostateData(m.name, m.charge, m.n_protons, m.free_energy, m.site_classes) for m in uncal.microstates],
        T,
        class_free_energy=truth,
    )
    pk = {(d["acid"], d["base"]): d["pKa"] for d in true_model.micro_pkas()}
    refs = [
        M.Reference("cation", "zwitterion", pk[("cation", "zwitterion")], "carboxylic"),
        M.Reference("zwitterion", "anion", pk[("zwitterion", "anion")], "ammonium"),
    ]
    cal = M.calibrate(free, refs, T)
    close(cal.class_free_energy["carboxylic"], truth["carboxylic"], 1e-9)
    close(cal.class_free_energy["ammonium"], truth["ammonium"], 1e-9)
    close(cal.shifts["carboxylic"], 3.0, 1e-9)
    assert all(abs(r) < 1e-9 for r in cal.residuals.values())

    # Through the result object, and the calibrated result reproduces the data.
    result = M.PKaResult(
        [M.MicrostateData(m.name, m.charge, m.n_protons, m.free_energy, m.site_classes) for m in uncal.microstates],
        T,
    ).with_calibration(refs)
    got = {(d["acid"], d["base"]): d["pKa"] for d in result.micro_pkas}
    close(got[("cation", "zwitterion")], refs[0].pka, 1e-9)
    close(got[("zwitterion", "anion")], refs[1].pka, 1e-9)
    assert result.calibration is not None and not any("uncalibrated" in w for w in result.warnings)


def test_calibration_is_a_mean_residual_and_a_relative_pka():
    free = {"A": 0.0, "B": 5.0, "C": 0.0, "D": 7.0}
    g_h = M.default_proton_free_energy(T)
    def p(a, b):
        return (free[b] + g_h - free[a]) / (RT * LN10)

    # Two references on one class, disagreeing by 0.4 units: the shift is their
    # mean residual and each is left with half the disagreement.
    refs = [M.Reference("A", "B", p("A", "B") + 1.0), M.Reference("C", "D", p("C", "D") + 1.4)]
    cal = M.calibrate(free, refs, T)
    close(cal.shifts["default"], 1.2 * RT * LN10, 1e-9)
    close(cal.residuals[("A", "B")], -0.2, 1e-9)
    close(cal.residuals[("C", "D")], 0.2, 1e-9)

    # One reference is a relative pKa: another acid is placed by its offset.
    one = M.calibrate(free, [M.Reference("A", "B", 4.76)], T)
    mdl = M.PKaModel(
        [ms("A", 0, 1, 0.0), ms("B", -1, 0, 5.0)], T, class_free_energy=one.class_free_energy
    )
    close(mdl.micro_pkas()[0]["pKa"], 4.76, 1e-9)
    raises(M.PKaError, M.calibrate, free, [M.Reference("A", "nope", 4.0)], T)
    raises(M.PKaError, M.calibrate, free, [], T)


# ---------------------------------------------------------------------------
#  Geometry and files
# ---------------------------------------------------------------------------

CREST_XYZ = """  3
      -74.96203654
O         -0.0000000000        0.0636623383        0.0000000001
H          0.7600000000       -0.5000000000        0.0000000000
H         -0.7600000000       -0.5000000000        0.0000000000
  3
      -74.96103654   some trailing words
O 0.0 0.0 0.0
H 0.0 0.9 0.0
H 0.9 0.0 0.0
"""


def test_geometry_and_xyz_parsing():
    ensemble = M.parse_xyz_ensemble(CREST_XYZ)
    assert len(ensemble) == 2
    close(ensemble[0][0], -74.96203654, 1e-12)
    close(ensemble[1][0], -74.96103654, 1e-12)
    assert ensemble[0][1].symbols == ["O", "H", "H"]
    again = M.parse_xyz_ensemble(ensemble[0][1].to_xyz("comment"))[0][1]
    assert again.symbols == ["O", "H", "H"]
    close(again.coords[1][0], 0.76, 1e-9)
    raises(ValueError, M.Geometry, ["O"], [[0.0, 0.0]])

    ok = "3\nmetalquicha converged, E =     -5.070544886125 Hartree\nO 0 0 0\nH 0 0 1\nH 0 1 0\n"
    geom, converged, energy = M.parse_optimized_xyz(ok)
    assert converged and geom.n_atoms == 3
    close(energy, -5.070544886125, 1e-12)
    bad = ok.replace("converged", "NOT CONVERGED")
    assert M.parse_optimized_xyz(bad)[1] is False


# ---------------------------------------------------------------------------
#  The result
# ---------------------------------------------------------------------------


def test_result_json_and_warnings():
    states = [
        M.MicrostateData("cation", 1, 2, -3.0, ("carboxylic", "ammonium"), {"n_conformers": 2, "n_imaginary": 1}),
        M.MicrostateData("zwitterion", 0, 1, 0.0, ("ammonium",), {"n_conformers": 1, "n_imaginary": 0}),
        M.MicrostateData("neutral", 0, 1, 1.0, ("carboxylic",)),
        M.MicrostateData("anion", -1, 0, 4.0),
    ]
    result = M.PKaResult(states, T, protocol={"qrrho": True})
    document = json.loads(result.to_json())
    assert {m["name"] for m in document["microstates"]} == {"cation", "zwitterion", "neutral", "anion"}
    assert len(document["micro_pkas"]) == 4  # cation -> {zwitterion, neutral} -> anion
    assert len(document["macro_pkas"]) == 2
    assert document["calibration"] is None
    text = "\n".join(result.warnings)
    assert "imaginary" in text and "within 1.00 kcal/mol" in text and "uncalibrated" in text
    # populations/charge delegate to the model.
    close(sum(result.populations(7.0).values()), 1.0, 1e-12)
    path = os.path.join(tempfile.mkdtemp(), "r.json")
    try:
        result.to_json(path)
        assert json.load(open(path))["temperature_K"] == T
    finally:
        shutil.rmtree(os.path.dirname(path))


# ---------------------------------------------------------------------------
#  The workflow, against a stand-in for the library
# ---------------------------------------------------------------------------


def _water(pka):
    return pka.Geometry(["O", "H", "H"], [[0.0, 0.0, 0.1], [0.0, 0.76, -0.47], [0.0, -0.76, -0.47]])


def _fake_executable(calls, ensemble_energies, opt_energies):
    """Stand-in for `_run_mqc`: writes what CREST and the optimizer would."""

    def run(executable, deck, directory, env):
        calls.append((executable, deck, directory, dict(env)))
        with open(os.path.join(directory, deck)) as handle:
            document = json.load(handle)
        assert os.path.exists(os.path.join(directory, "start.xyz"))
        if document["driver"] == "conformers":
            body = ""
            for e in ensemble_energies:
                body += f"3\n   {e:.8f}\nO 0 0 0.1\nH 0 0.76 -0.47\nH 0 -0.76 -0.47\n"
            with open(os.path.join(directory, "crest_conformers.xyz"), "w") as out:
                out.write(body)
        else:
            e = opt_energies.pop(0)
            name = f"output_{deck[:-5]}_optimized.xyz"
            with open(os.path.join(directory, name), "w") as out:
                out.write(f"3\nmetalquicha converged, E = {e:20.12f} Hartree\nO 0 0 0.1\nH 0 0.76 -0.47\nH 0 -0.76 -0.47\n")

    return run


def test_workflow_end_to_end_against_stand_in():
    log, calls = [], []
    workdir = tempfile.mkdtemp()
    try:
        with fake_package(log) as pka:
            # Third conformer is 5 kcal/mol up and outside the 3 kcal/mol window;
            # the second and third optimize to the same energy and are one.
            ensemble = [-75.0, -75.0 + 0.001, -75.0 + 0.008]
            pka._run_mqc = _fake_executable(calls, ensemble, [-75.1] * 4)
            pka._find_executable = lambda protocol: "/fake/mqc"

            protocol = pka.Protocol(workdir=workdir, energy_window_kcal=3.0, max_conformers=5)
            acid = pka.Microstate("acetic acid (neutral)", _water(pka), 0, 1)
            base = pka.Microstate("acetate", _water(pka), -1, 0)
            result = pka.run([acid, base], protocol, prefix="t")

            # Stages and decks.
            decks = [c for c in calls if c[1] == "conformers.json"]
            opts = [c for c in calls if c[1] == "optimize.json"]
            assert len(decks) == 2 and len(opts) == 4  # 2 conformers survive the window, per microstate
            assert calls[0][0] == "/fake/mqc"
            with open(os.path.join(decks[0][2], "conformers.json")) as handle:
                deck = json.load(handle)
            assert deck["driver"] == "conformers" and "fragmentation" not in deck["keywords"]
            assert deck["model"]["method"] == "gfn2"
            assert deck["keywords"]["xtb"] == {"solvent": "water", "solvation_model": "alpb"}
            assert deck["molecules"][0]["molecular_charge"] == 0 and deck["molecules"][0]["xyz"] == "start.xyz"
            with open(os.path.join(opts[2][2], "optimize.json")) as handle:
                assert json.load(handle)["molecules"][0]["molecular_charge"] == -1
            # No MPI launcher state leaks into the child.
            assert not any(k.startswith(("OMPI_", "PMIX_", "PMI_")) for c in calls for k in c[3])

            # In-session runs: Hessians, and no separate single points because the
            # default single point is the frequency calculation.
            drivers = [e["driver"] for e in log]
            assert set(drivers) == {"hessian"} and len(drivers) == 2 * 1  # duplicates collapse to one conformer each
            for entry in log:
                assert entry["system"].charges in ([0], [-1])
            assert [e["label"] for e in log] == ["t_acetic_acid__neutral__c0_freq", "t_acetate_c0_freq"]

            data = {m.name: m for m in result.microstates}
            assert set(data) == {"acetic acid (neutral)", "acetate"}
            for m in result.microstates:
                assert math.isfinite(m.free_energy) and m.detail["n_conformers"] == 1
                rec = m.detail["conformers"][0]
                close(rec["g_kcal_mol"], rec["e_single_point_hartree"] * M.HARTREE_TO_KCAL + rec["g_corr_kcal_mol"], 1e-6)
                close(rec["weight"], 1.0, 1e-12)
            assert len(result.macro_pkas) == 1 and result.pI is None
            json.loads(result.to_json())
            assert result.protocol["energy_window_kcal"] == 3.0

            # Calibrate afterwards, as the docs recommend, with a tuple reference.
            cal = result.with_calibration([("acetic acid (neutral)", "acetate", 4.76)])
            close(cal.micro_pkas[0]["pKa"], 4.76, 1e-9)
            close(result.micro_pkas[0]["pKa"], result.macro_pkas[0]["pKa"], 1e-9)  # one microstate per level
    finally:
        shutil.rmtree(workdir)


def test_workflow_separate_single_point_and_imaginary_surfaced():
    log, calls = [], []
    workdir = tempfile.mkdtemp()
    try:
        with fake_package(log) as pka:
            sp = {"method": "gfn1", "verbosity": "error"}  # a different level: must be its own run
            protocol = pka.Protocol(conformers=None, optimize=None, single_point=sp, workdir=workdir)
            assert not protocol.single_point_is_frequencies()
            a = pka.Microstate("imag-species", _water(pka), 0, 1)
            b = pka.Microstate("b", _water(pka), -1, 0)
            result = pka.run([a, b], protocol, prefix="p", verbose=False)
            assert calls == []  # both subprocess stages skipped
            assert [e["driver"] for e in log] == ["hessian", "energy", "hessian", "energy"]
            assert log[1]["kwargs"]["model"]["method"] == "gfn1"
            imag = {m.name: m.detail["n_imaginary"] for m in result.microstates}
            assert imag == {"imag-species": 1, "b": 0}
            assert any("imag-species: 1 imaginary frequency" in w for w in result.warnings)
    finally:
        shutil.rmtree(workdir)


def test_workflow_refuses_bad_input_before_running():
    log = []
    with fake_package(log) as pka:
        raises(ValueError, pka.Protocol, frequencies=None)
        raises(ValueError, pka.Protocol, single_point={"method": "gfn2", "driver": "energy"})
        raises(ValueError, pka.Protocol, max_conformers=0)
        a = pka.Microstate("a", _water(pka), 0, 1)
        raises(ValueError, pka.run, [])
        raises(ValueError, pka.run, [a, pka.Microstate("a", _water(pka), -1, 0)])
        raises(ValueError, pka.run, [a, pka.Microstate("a.", _water(pka), -1, 0), pka.Microstate("a_", _water(pka), -1, 0)])
        # charge / proton bookkeeping is checked before any calculation.
        raises(ValueError, pka.run, [a, pka.Microstate("b", _water(pka), 0, 0)])
        assert log == []
        # A System-only microstate cannot be searched.
        system = pka.Microstate("sys", pka.__dict__["Geometry"](["H"], [[0, 0, 0]]), 0, 0)
        fake = type("S", (), {"_handle": 1, "n_atoms": 1})()
        only = pka.Microstate("only", fake, 0, 0)
        raises(pka.PKaError, lambda: only.geometry)
        assert system.n_atoms == 1


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
