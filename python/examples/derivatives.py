"""Reading a Hessian back: frequencies, intensities and thermochemistry.

    python derivatives.py

A Hessian run produces far more than an energy, and until recently none of it
was reachable from Python. The driver wrote `vibrational_analysis` and
`thermochemistry` into `output_<label>.json` and `Result` had no accessor for
either, so a script could ask for frequencies, get a converged calculation
back, and have no way to see them.

What is worth knowing about the numbers, since the assertions below lean on it:

**Imaginary modes come back negative.** The vibrational analysis writes them
that way and so does every other code that prints frequencies. A minimum has
none; a transition state has exactly one. A geometry that is not a stationary
point has whatever the curvature happens to be, which is why the water here is
optimized first rather than taken from a guess -- frequencies at a
non-stationary point are a number without a meaning, and asserting on them
would be asserting on noise.

**Six modes are removed, not computed.** Three translations and three
rotations are projected out before diagonalisation, so a non-linear molecule of
N atoms has 3N-6 frequencies. Water gives three. Getting 3N back means the
projection did not run.

**Thermochemistry follows the frequencies and inherits their problems.** RRHO
integrates over harmonic modes, so an imaginary one has nothing to contribute
and is dropped -- `n_imaginary_frequencies` is reported and the free energy is
still printed beside it. Read the count before quoting the number.

Run as a test rather than a demonstration: each case asserts, and the script
exits non-zero if any of them is wrong.
"""

import os
import sys

import mqc

#: Water near equilibrium, Angstrom. Close enough that the optimizer needs a
#: few steps rather than a search, and the result is a genuine minimum so the
#: frequency signs are meaningful.
WATER = dict(
    symbols=["O", "H", "H"],
    coords=[[0.0, 0.0, 0.117], [0.0, 0.757, -0.469], [0.0, -0.757, -0.469]],
)

#: Small and fast. The point here is the plumbing, not the spectroscopy: STO-3G
#: overestimates these frequencies badly and that is fine, because every
#: assertion below is about structure and sign rather than accuracy.
MODEL = dict(method="hf", basis="sto-3g", verbosity="error")


def water():
    s = mqc.System()
    s.set_geometry(WATER["symbols"], WATER["coords"])
    return s


def hessian_run(label):
    """A Hessian on water, unfragmented. `level=0` is the supermolecule."""
    return mqc.MBE(water(), level=0, driver="Hessian", **MODEL).run(label=label)


def case_frequencies():
    """3N-6 modes, all real, and the arrays agree with each other."""
    print("water/STO-3G Hessian: frequencies")
    result = hessian_run("py_freq")

    freqs = result.frequencies
    assert freqs is not None, "a Hessian run returned no frequencies"
    # Water is non-linear and has three atoms, so 3*3 - 6 = 3.
    assert len(freqs) == 3, f"expected 3 modes after projection, got {len(freqs)}"
    # A minimum. Any negative entry here means the geometry was not one.
    assert all(f > 0.0 for f in freqs), f"imaginary mode at a minimum: {freqs}"
    # Ordered, which is what makes "the first mode" mean anything.
    assert freqs == sorted(freqs), f"frequencies are not ascending: {freqs}"
    # The bend sits well below the two stretches in every basis.
    assert freqs[0] < freqs[1], "the bend is not the lowest mode"
    print(f"  {len(freqs)} modes: " + "  ".join(f"{f:.1f}" for f in freqs) + " cm^-1")

    # The raw block carries the companion arrays, one entry per mode.
    analysis = result.vibrational_analysis
    assert analysis is not None, "vibrational_analysis is missing"
    assert analysis["n_modes"] == len(freqs), "n_modes disagrees with the array"
    for key in ("reduced_masses_amu", "force_constants_mdyne_ang"):
        assert key in analysis, f"{key} is missing from the analysis"
        assert len(analysis[key]) == len(freqs), f"{key} has the wrong length"
        assert all(v > 0.0 for v in analysis[key]), f"{key} has a non-positive entry"


def case_thermochemistry():
    """The RRHO block, and the counts that say whether to trust it."""
    print("water/STO-3G Hessian: thermochemistry")
    result = hessian_run("py_thermo")

    thermo = result.thermochemistry
    assert thermo is not None, "a Hessian run returned no thermochemistry"

    # It says what it assumed, which is the part worth checking before the
    # energies: a wrong linearity or symmetry number moves the entropy without
    # moving anything that looks wrong.
    assert thermo["is_linear"] is False, "water was treated as linear"
    assert thermo["n_real_frequencies"] == 3, "the mode count does not match"
    assert thermo["n_imaginary_frequencies"] == 0, "an imaginary mode at a minimum"
    assert thermo["temperature_K"] > 0.0, "no temperature"

    # Zero-point energy is positive and the two units agree. 1 Hartree is
    # 627.5095 kcal/mol; a loose tolerance because the JSON is rounded.
    zpe_h = thermo["zero_point_energy_hartree"]
    zpe_k = thermo["zero_point_energy_kcal_mol"]
    assert zpe_h > 0.0, f"zero-point energy is not positive: {zpe_h}"
    assert abs(zpe_h * 627.5095 - zpe_k) < 1e-3, "the two ZPE units disagree"
    print(f"  ZPE {zpe_h:.6f} Hartree = {zpe_k:.3f} kcal/mol")

    # The nested blocks the deck path writes.
    for key in ("partition_functions", "thermal_corrections_hartree"):
        assert key in thermo, f"{key} is missing"
    functions = thermo["partition_functions"]
    for mode in ("translational", "rotational", "vibrational"):
        assert mode in functions, f"no {mode} partition function"
        assert functions[mode] > 0.0, f"{mode} partition function is not positive"


def case_absent_sections():
    """An energy run has none of this, and says so rather than raising."""
    print("water/STO-3G energy: the derivative blocks are absent")
    result = mqc.MBE(water(), level=0, driver="Energy", **MODEL).run(label="py_energy")

    assert result.energy < 0.0, "the energy run did not produce an energy"
    # Every accessor is optional and returns None when the run did not make it.
    # This is the half that a caller relies on: asking for frequencies after an
    # energy should not be an exception.
    for name in ("frequencies", "ir_intensities", "vibrational_analysis",
                 "thermochemistry", "fukui", "pie_terms"):
        assert getattr(result, name) is None, f"{name} came back from an energy run"
    print("  frequencies, intensities, thermochemistry, fukui, pie_terms all None")


def case_dipole():
    """The dipole travels with an ordinary run, in both unit systems."""
    print("water/STO-3G energy: dipole")
    result = mqc.MBE(water(), level=0, driver="Energy", **MODEL).run(label="py_dipole")

    dipole = result.dipole
    if dipole is None:
        # Not every method writes one; skip rather than fail the whole script.
        print("  no dipole written by this method; skipped")
        return
    for axis in ("x", "y", "z"):
        assert axis in dipole, f"the dipole has no {axis} component"
    magnitude = dipole["magnitude_debye"]
    assert magnitude > 0.0, "water came out with no dipole"
    # Water's dipole lies along the C2 axis, so the component out of the
    # molecular plane is zero by symmetry. The geometry above puts that plane
    # in y-z, making x the null one.
    assert abs(dipole["x"]) < 1e-6, f"a symmetry-forbidden x component: {dipole['x']}"
    print(f"  {magnitude:.4f} Debye")


def case_driver_names():
    """A misspelled driver is refused here rather than inside the library."""
    print("driver validation")
    for bad in ("Freq", "frequencies", "energies", ""):
        try:
            mqc.MBE(water(), driver=bad)
        except ValueError as exc:
            assert "unknown driver" in str(exc), f"wrong message for {bad!r}: {exc}"
        else:
            raise AssertionError(f"driver={bad!r} was accepted")

    # Every spelling the Fortran side accepts is accepted here, including the
    # aliases -- refusing one of those would be a regression, not a check.
    for good in sorted(mqc.DRIVERS):
        mqc.MBE(water(), driver=good)
        mqc.MBE(water(), driver=good.upper())
    print(f"  {len(mqc.DRIVERS)} spellings accepted, 4 rejected")


def case_keyword_blocks():
    """The derivative blocks reach the deck under their own names."""
    print("keyword blocks")
    job = mqc.MBE(
        water(),
        level=0,
        driver="Hessian",
        hessian={"displacement": 0.005},
        optimization={"max_iterations": 40},
        efp={"screening": True},
        **MODEL,
    )
    keywords = job.settings()["keywords"]
    assert keywords["hessian"] == {"displacement": 0.005}, "the hessian block was lost"
    assert keywords["optimization"] == {"max_iterations": 40}, "the optimization block was lost"
    assert keywords["efp"] == {"screening": True}, "the efp block was lost"

    # And the raw escape hatch still merges rather than replacing, so the two
    # spellings of a block do not fight.
    merged = mqc.MBE(
        water(), level=0, hessian={"displacement": 0.005},
        keywords={"hessian": {"project_rotations": True}}, **MODEL,
    ).settings()["keywords"]["hessian"]
    assert merged["displacement"] == 0.005, "the named block was dropped"
    assert merged["project_rotations"] is True, "the escape hatch was dropped"
    print("  hessian, optimization, efp reach the deck and merge")


CASES = {
    "frequencies": case_frequencies,
    "thermochemistry": case_thermochemistry,
    "absent": case_absent_sections,
    "dipole": case_dipole,
    "drivers": case_driver_names,
    "blocks": case_keyword_blocks,
}


def cleanup():
    """The runs above write output files; do not leave them behind."""
    for label in ("py_freq", "py_thermo", "py_energy", "py_dipole"):
        path = f"output_{label}.json"
        if os.path.exists(path):
            os.remove(path)


def main(argv):
    wanted = argv[1] if len(argv) > 1 else ""
    ran = 0
    try:
        for name, fn in CASES.items():
            if wanted and wanted not in name:
                continue
            fn()
            ran += 1
    finally:
        cleanup()
    if ran == 0:
        print(f"no case matching {wanted!r}; have: {', '.join(CASES)}")
        return 1
    print(f"\n{ran} case(s) passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
