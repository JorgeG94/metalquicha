"""Reading a Hessian back: frequencies, intensities and thermochemistry.

    python derivatives.py

A Hessian run produces far more than an energy, and until recently none of it
was reachable from Python. The driver wrote `vibrational_analysis` and
`thermochemistry` into `output_<label>.json` and `Result` had no accessor for
either, so a script could ask for frequencies, get a converged calculation
back, and have no way to see them.

What is worth knowing about the numbers, since the assertions below lean on it:

**All 3N modes come back, rigid-body ones included.** They are not projected
out: water returns nine frequencies, six of them within a hundredth of a
wavenumber of zero and three genuine vibrations. A caller wanting the
vibrations filters by magnitude, which is what the assertions below do. This
is worth knowing before comparing a mode count against another program, most
of which print 3N-6.

**Imaginary modes come back negative**, which is the convention every code
that prints frequencies uses. The catch is that the six rigid-body modes sit
at zero and their *sign* is numerical noise, so a perfectly good minimum
reports `n_imaginary_frequencies` of three. Water below is a minimum and does
exactly that. The count is about the near-zero modes, not about the chemistry;
what says a structure is a saddle point is a mode of real magnitude that came
out negative.

**Thermochemistry is built from the genuine vibrations only.** The zero-point
energy here is `0.5 * sum` over the three real modes and nothing else, which
the assertions check directly -- so the mode counts beside it can be confusing
while the energies are right.

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
    # The whole molecule as one monomer. The validator wants a partition
    # whether or not the expansion will use one, and `auto_monomers` refuses a
    # single covalently bonded molecule -- rightly, since where to cut it is a
    # chemical choice. At level 0 nothing is cut, so "all of it, once" is the
    # honest way to ask.
    s.set_monomers([[0, 1, 2]])
    return s


def hessian_run(label):
    """A Hessian on water, unfragmented. `level=0` is the supermolecule."""
    return mqc.MBE(water(), level=0, driver="Hessian", **MODEL).run(label=label)


#: Anything under this in cm^-1 is a rigid-body mode rather than a vibration.
#: Water's lowest genuine mode is the bend, two thousand wavenumbers up, so
#: nothing here is near the boundary.
RIGID_BODY_CUTOFF = 100.0

#: Wavenumbers per Hartree, for checking the zero-point energy against the
#: frequencies it was built from.
CM1_PER_HARTREE = 219474.6313632


def case_frequencies():
    """3N modes come back; three of them are the vibrations."""
    print("water/STO-3G Hessian: frequencies")
    result = hessian_run("py_freq")

    freqs = result.frequencies
    assert freqs is not None, "a Hessian run returned no frequencies"
    # All 3N, rigid-body modes included -- they are not projected out.
    assert len(freqs) == 9, f"expected 3N = 9 modes, got {len(freqs)}"
    # Ordered, which is what makes "the first mode" mean anything.
    assert freqs == sorted(freqs), f"frequencies are not ascending: {freqs}"

    vibrations = [f for f in freqs if abs(f) > RIGID_BODY_CUTOFF]
    rigid = [f for f in freqs if abs(f) <= RIGID_BODY_CUTOFF]
    # 3N-6 genuine vibrations for a non-linear three-atom molecule.
    assert len(vibrations) == 3, f"expected 3 vibrations, got {vibrations}"
    assert len(rigid) == 6, f"expected 6 rigid-body modes, got {rigid}"
    # The six sit at zero. Their sign is noise; their magnitude should not be.
    assert all(abs(f) < 1.0 for f in rigid), f"a rigid-body mode is not at zero: {rigid}"
    # A minimum: every genuine mode is positive.
    assert all(f > 0.0 for f in vibrations), f"imaginary vibration: {vibrations}"
    # The bend sits well below the two stretches in every basis.
    assert vibrations[0] < vibrations[1], "the bend is not the lowest vibration"
    print("  vibrations: " + "  ".join(f"{f:.1f}" for f in vibrations) + " cm^-1")

    # The raw block carries the companion arrays, one entry per mode.
    analysis = result.vibrational_analysis
    assert analysis is not None, "vibrational_analysis is missing"
    assert analysis["n_modes"] == len(freqs), "n_modes disagrees with the array"
    for key in ("reduced_masses_amu", "force_constants_mdyne_ang",
                "ir_intensities_km_mol"):
        assert key in analysis, f"{key} is missing from the analysis"
        assert len(analysis[key]) == len(freqs), f"{key} has the wrong length"
    # A reduced mass is a mass and is positive for every mode, rigid-body
    # included. A force constant is zero for those and positive for the rest,
    # which is the same statement as the frequencies in different units.
    assert all(m > 0.0 for m in analysis["reduced_masses_amu"]), "a non-positive mass"
    constants = [k for f, k in zip(freqs, analysis["force_constants_mdyne_ang"])
                 if abs(f) > RIGID_BODY_CUTOFF]
    assert all(k > 0.0 for k in constants), f"a vibration with no restoring force: {constants}"

    # And the list accessor is the same numbers as the block it comes from.
    assert result.ir_intensities == analysis["ir_intensities_km_mol"], \
        "ir_intensities disagrees with the block it reads"


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
    assert thermo["temperature_K"] > 0.0, "no temperature"
    # The counts are over all 3N modes, so the six rigid-body ones land in one
    # bucket or the other according to the sign of their numerical noise. Water
    # here is a minimum and still reports three imaginary. What the count can
    # say is that no *genuine* mode went negative, which bounds it at six.
    counted = thermo["n_real_frequencies"] + thermo["n_imaginary_frequencies"]
    assert counted == 9, f"the counts cover {counted} modes, not 3N"
    assert thermo["n_imaginary_frequencies"] <= 6, \
        "a genuine vibration came out imaginary at a minimum"

    # Zero-point energy is positive and the two units agree. 1 Hartree is
    # 627.5095 kcal/mol; a loose tolerance because the JSON is rounded.
    zpe_h = thermo["zero_point_energy_hartree"]
    zpe_k = thermo["zero_point_energy_kcal_mol"]
    assert zpe_h > 0.0, f"zero-point energy is not positive: {zpe_h}"
    assert abs(zpe_h * 627.5095 - zpe_k) < 1e-3, "the two ZPE units disagree"

    # And it is built from the genuine vibrations alone, which is the check
    # that matters given the mode counts above: half the sum of the three real
    # frequencies, converted from cm^-1. Getting the rigid-body modes into this
    # sum would not move it much, which is exactly why it needs asserting.
    vibrations = [f for f in result.frequencies if abs(f) > RIGID_BODY_CUTOFF]
    expected = 0.5 * sum(vibrations) / CM1_PER_HARTREE
    assert abs(zpe_h - expected) < 1e-6, \
        f"ZPE {zpe_h} is not half the sum of {vibrations}"
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
    # The session wraps everything, so a traceback closes it rather than
    # leaving the library holding a communicator.
    try:
        with mqc.session():
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
