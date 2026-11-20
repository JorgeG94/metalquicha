"""Tests for the spectrum broadening, which needs neither the library nor a build.

    python3 python/tests/test_spectra.py

Exits non-zero on the first failure; also collected by pytest if that is to hand.
`mqc/spectra.py` imports nothing from the package, so it is loaded by path, as
the arithmetic modules of `test_pka.py` and `test_bde.py` are.

Every number checked here is one that can be derived by hand: a line shape's
area, the integral of a broadened band against the strength it was made from,
a published conversion factor, a cross section worked out on paper. What it
cannot test is that a real run writes the keys read here; the keys are the ones
`mqc_json_writer.f90` writes, and `python/examples/spectra.py` reads a real
output.
"""

import importlib.util
import json
import math
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(HERE, "..", "mqc")


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(PKG, filename))
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


spectra = _load("mqc_spectra_under_test", "spectra.py")


def close(a, b, rel=1.0e-9, abs_tol=0.0):
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_tol)


def area(x, y):
    """Trapezoid rule on a grid that need not be uniform."""
    return sum(
        0.5 * (y[k] + y[k + 1]) * (x[k + 1] - x[k]) for k in range(len(x) - 1)
    )


def vibrational_output(with_raman=True):
    """What a Hessian run on a bent triatomic writes, with invented numbers.

    All 3N = 9 modes are in the lists, as the Fortran writes them: six near
    zero for the translations and rotations, then the three vibrations.
    """
    block = {
        "n_modes": 9,
        "frequencies_cm1": [-0.03, -0.01, 0.0, 0.01, 0.02, 0.04, 1600.0, 3700.0, 3800.0],
        "ir_intensities_km_mol": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 70.0, 5.0, 50.0],
    }
    if with_raman:
        block["raman_activities_A4_amu"] = [0.3, 0.2, 0.1, 0.1, 0.2, 0.3, 4.0, 90.0, 40.0]
    return {"water": {"total_energy": -76.0, "vibrational_analysis": block}}


def excited_output():
    return {
        "benzene": {
            "excited_states": {
                "n_states": 3,
                "states": [
                    {"state": 1, "excitation_energy_ev": 4.0},
                    {
                        "state": 2,
                        "excitation_energy_ev": 5.0,
                        "oscillator_strength": 0.10,
                        "oscillator_strength_velocity": 0.08,
                    },
                    {
                        "state": 3,
                        "excitation_energy_ev": 6.2,
                        "oscillator_strength": 0.60,
                        "oscillator_strength_velocity": 0.55,
                    },
                ],
            }
        }
    }


# -- line shapes ------------------------------------------------------------


def test_lineshapes_have_unit_area_and_the_stated_width():
    grid = [k * 0.01 for k in range(-200000, 200001)]
    gauss = [spectra.gaussian(x, 0.0, 10.0) for x in grid]
    assert close(area(grid, gauss), 1.0, rel=1.0e-9)
    # A Lorentzian's tails are long: 2000 either side of a width of 10 still
    # leaves 2 * 5 / (pi * 2000) outside.
    lorentz = [spectra.lorentzian(x, 0.0, 10.0) for x in grid]
    assert close(area(grid, lorentz), 1.0 - 10.0 / (math.pi * 2000.0), rel=1.0e-6)
    for shape in (spectra.gaussian, spectra.lorentzian):
        assert close(shape(5.0, 0.0, 10.0), 0.5 * shape(0.0, 0.0, 10.0), rel=1.0e-12)
        assert close(shape(-5.0, 0.0, 10.0), 0.5 * shape(0.0, 0.0, 10.0), rel=1.0e-12)


# -- constants --------------------------------------------------------------


def test_conversion_factors_are_the_published_ones():
    assert close(spectra.EV_TO_CM1, 8065.543937, rel=1.0e-9)
    # A Lorentzian band of A km/mol and width w peaks at 27.648 A / w.
    peak = spectra.IR_KM_MOL_TO_EPSILON * 2.0 / math.pi
    assert close(peak, 27.648, rel=1.0e-4)
    # f = 4.32e-9 times the integrated molar absorptivity over wavenumber.
    assert close(spectra.F_PER_AREA, 4.319e-9, rel=1.0e-3)


# -- infrared ---------------------------------------------------------------


def test_ir_keeps_the_vibrations_and_conserves_each_band():
    result = spectra.ir(vibrational_output(), fwhm=10.0, lineshape="gaussian", step=0.5)
    assert [s[0] for s in result.sticks] == [1600.0, 3700.0, 3800.0]
    assert [s[1] for s in result.sticks] == [70.0, 5.0, 50.0]
    assert result.x[0] == 0.0 and result.x[-1] >= 3800.0 + 100.0
    # The area is ln10/100 of the intensity summed, by the definition of km/mol.
    total = area(result.x, result.y) * math.log(10.0) / 100.0
    assert close(total, 125.0, rel=1.0e-6), total
    tallest = result.peaks()[0]
    assert close(tallest[0], 1600.0, abs_tol=0.5)
    assert not result.warnings


def test_ir_scale_moves_the_lines_and_not_their_strengths():
    result = spectra.ir(vibrational_output(), scale=0.95)
    assert close(result.sticks[0][0], 1520.0)
    assert result.sticks[0][1] == 70.0


def test_imaginary_modes_are_left_out_and_reported():
    output = vibrational_output()
    output["water"]["vibrational_analysis"]["frequencies_cm1"][0] = -450.0
    result = spectra.ir(output)
    assert len(result.sticks) == 3
    assert len(result.warnings) == 1 and "450i" in result.warnings[0]


def test_missing_sections_are_refused_with_what_to_ask_for():
    for call, needle in (
        (lambda: spectra.ir({"x": {"total_energy": 1.0}}), "Hessian"),
        (lambda: spectra.raman(vibrational_output(with_raman=False)), "properties"),
        (lambda: spectra.uvvis(vibrational_output()), "excited_states"),
        (lambda: spectra.ir(vibrational_output(), lineshape="voigt"), "lorentzian"),
        (lambda: spectra.ir(vibrational_output(), fwhm=0.0), "positive"),
    ):
        try:
            call()
        except spectra.SpectrumError as problem:
            assert needle in str(problem), problem
        else:
            raise AssertionError("accepted what it should have refused: " + needle)


# -- Raman ------------------------------------------------------------------


def test_raman_cross_section_matches_a_hand_calculation():
    # Nitrogen-like: 2331 cm^-1, 20 A^4/amu, the 488 nm argon line, no
    # population factor. Each factor in CGS, multiplied out by hand:
    #   (2 pi)^4 / 45                     = 34.6347
    #   (20491.80 - 2331)^4               = 1.08779e17 cm^-4
    #   h / (8 pi^2 c nu)                 = 1.20090e-42 g cm^2
    #   S = 20e-32 / 1.66054e-24          = 1.20443e-7 cm^4/g
    value = spectra.raman_cross_section(2331.0, 20.0, 1.0e7 / 488.0, 0.0)
    assert close(value, 34.6347 * 1.08779e17 * 1.20090e-42 * 1.20443e-7, rel=2.0e-4)
    # Nitrogen at 488 nm is measured at about 5e-31 cm^2/sr.
    assert 1.0e-31 < value < 1.0e-30

    # The population factor is 1 / (1 - exp(-h c nu / k T)): negligible for a
    # stretch at room temperature, a factor of several for a soft mode.
    cold = spectra.raman_cross_section(100.0, 1.0, 1.0e7 / 532.0, 0.0)
    warm = spectra.raman_cross_section(100.0, 1.0, 1.0e7 / 532.0, 298.15)
    expected = 1.0 / (1.0 - math.exp(-1.4387768775 * 100.0 / 298.15))
    assert close(warm / cold, expected, rel=1.0e-8)
    assert expected > 2.5


def test_raman_weights_low_modes_up_against_the_activities():
    output = vibrational_output()
    by_activity = spectra.raman(output, activities=True)
    by_section = spectra.raman(output, laser_nm=532.0)
    assert [s[1] for s in by_activity.sticks] == [4.0, 90.0, 40.0]
    assert "activity" in by_activity.y_label and "cross section" in by_section.y_label
    # The bend's share of the stretch is larger as a cross section than as an
    # activity: 1/nu and a smaller Stokes shift both favour it.
    share_activity = by_activity.sticks[0][1] / by_activity.sticks[1][1]
    share_section = by_section.sticks[0][1] / by_section.sticks[1][1]
    assert share_section > 2.0 * share_activity
    # And the broadened curve integrates to the sum of the cross sections.
    total = area(by_section.x, by_section.y)
    assert close(total, sum(s[1] for s in by_section.sticks), rel=2.0e-2)


def test_raman_refuses_a_laser_below_the_mode():
    try:
        spectra.raman(vibrational_output(), laser_nm=5000.0)
    except spectra.SpectrumError as problem:
        assert "Stokes" in str(problem)
    else:
        raise AssertionError("a 2000 cm^-1 laser gave a line at 3800 cm^-1")


# -- UV-vis -----------------------------------------------------------------


def test_uvvis_band_area_gives_back_the_oscillator_strength():
    result = spectra.uvvis(excited_output(), fwhm_ev=0.3, unit="cm-1", points=20000)
    # The dark state is dropped and said to be; the bright ones are kept.
    assert [s[1] for s in result.sticks] == [0.10, 0.60]
    assert any("1 state(s) carry no oscillator strength" in w for w in result.warnings)
    total = area(result.x, result.y) * spectra.F_PER_AREA
    assert close(total, 0.70, rel=1.0e-4), total

    velocity = spectra.uvvis(excited_output(), gauge="velocity", unit="ev")
    assert [s[1] for s in velocity.sticks] == [0.08, 0.55]


def test_uvvis_axis_units_agree_with_one_another():
    in_ev = spectra.uvvis(excited_output(), unit="ev", points=4001)
    in_nm = spectra.uvvis(excited_output(), unit="nm", points=4001)
    assert close(in_nm.sticks[1][0], 1239.841984 / 6.2, rel=1.0e-8)
    assert in_nm.x[0] < in_nm.x[-1]
    # The same curve on either axis: its maximum is one height, at one energy.
    top_ev = max(zip(in_ev.y, in_ev.x))
    top_nm = max(zip(in_nm.y, in_nm.x))
    assert close(top_ev[0], top_nm[0], rel=1.0e-4)
    assert close(1239.841984 / top_nm[1], top_ev[1], rel=1.0e-3)
    # The band at 6.2 eV dominates, so that is where the maximum is.
    assert close(top_ev[1], 6.2, abs_tol=0.01)

    shifted = spectra.uvvis(excited_output(), unit="ev", shift_ev=-0.5)
    assert close(shifted.sticks[0][0], 4.5)


# -- files ------------------------------------------------------------------


def test_reads_a_file_and_writes_csv():
    folder = tempfile.mkdtemp()
    source = os.path.join(folder, "output_water.json")
    with open(source, "w") as handle:
        json.dump(vibrational_output(), handle)
    result = spectra.ir(source, xmin=1000.0, xmax=2000.0, step=10.0)
    curve = os.path.join(folder, "ir.csv")
    sticks = os.path.join(folder, "sticks.csv")
    result.to_csv(curve, sticks_path=sticks)
    with open(curve) as handle:
        lines = handle.read().splitlines()
    assert lines[0] == "# ir"
    assert lines[1] == "wavenumber / cm^-1,molar absorptivity / L mol^-1 cm^-1"
    assert len(lines) == 2 + 101
    first = lines[2].split(",")
    assert close(float(first[0]), 1000.0) and close(float(first[1]), result.y[0], rel=1e-7)
    with open(sticks) as handle:
        assert len(handle.read().splitlines()) == 2 + 3

    # The command line is the same thing with the arguments spelled out.
    code = spectra.main(["raman", source, "--csv", os.path.join(folder, "r.csv")])
    assert code == 0
    assert spectra.main(["uvvis", source]) == 1


def main():
    tests = [value for name, value in sorted(globals().items()) if name.startswith("test_")]
    for test in tests:
        test()
        print("ok  %s" % test.__name__)
    print("%d passed" % len(tests))


if __name__ == "__main__":
    main()
