"""pKa and isoelectric point, end to end: GFN2-xTB with ALPB water.

    python pka.py            conformer search, optimization, Hessian, single point
    python pka.py --quick    skip the conformer search (one structure each)

Acetic acid / acetate for a pKa, and glycine's four microstates for the pKas
and the isoelectric point. The geometries are hand-built and deliberately
rough -- the optimization stage is what makes them minima.

**What this asserts is identities and sanity, not agreement with experiment.**
Uncalibrated GFN2/ALPB pKas are expected to be several units off, which is the
reason `calibrate` exists; the experimental numbers are printed beside the
computed ones so the size of that is visible. What is asserted:

  * every free energy is finite and every microstate has a conformer;
  * where two microstates are the only ones at their level, the micro and
    macro pKa are the same number;
  * at pH = pKa the two microstates are equally populated, and the net charge
    never rises with pH;
  * a calibrated reference is reproduced exactly (that is what a
    calibration is) and a transferred shift gives glycine pKas that are
    ordered, with the isoelectric point between them.

Two of the four stages run the ``mqc`` executable, because the library refuses
the ``optimize`` and ``conformers`` drivers through a session. A build without
tblite, without the executable, or without CREST is reported and skipped rather
than failed -- except CREST, where the run falls back to no conformer search
and says so, because a build without it is still a build.

Run from a scratch directory: the Hessian and energy runs leave
``output_pka_*.json`` in the working directory and the executable stages work in
``pka_example_work/``.
"""

import math
import sys

import mqc
from mqc import pka

# -- hand-built geometries, Angstrom ----------------------------------------

ACETIC_ACID = (
    ["C", "C", "O", "O", "H", "H", "H", "H"],
    [[0.0000, 1.5100, 0.0000], [0.0000, 0.0000, 0.0000], [1.0392, -0.6000, 0.0000],
     [-1.1691, -0.6750, 0.0000], [-0.9509, -1.6201, 0.0000], [-1.0277, 1.8733, 0.0000],
     [0.5138, 1.8733, -0.8900], [0.5138, 1.8733, 0.8900]],
)

ACETATE = (
    ["C", "C", "O", "O", "H", "H", "H"],
    [[0.0000, 1.5400, 0.0000], [0.0000, 0.0000, 0.0000], [1.1176, -0.5818, 0.0000],
     [-1.1176, -0.5818, 0.0000], [-1.0277, 1.9033, 0.0000], [0.5138, 1.9033, -0.8900],
     [0.5138, 1.9033, 0.8900]],
)

#: N, CA, C, O, O, then the hydrogens: on CA, then on N, then on the carboxyl O.
GLYCINE_CATION = (  # +NH3-CH2-COOH
    ["N", "C", "C", "O", "O", "H", "H", "H", "H", "H", "H"],
    [[-1.2731, -0.7350, 0.0000], [0.0000, 0.0000, 0.0000], [1.1644, -0.9770, 0.0000],
     [2.3014, -0.5632, 0.0000], [0.9317, -2.2967, 0.0000], [1.7960, -2.7371, 0.0000],
     [0.0548, 0.6269, 0.8900], [0.0548, 0.6269, -0.8900], [-1.0848, -1.7477, 0.0000],
     [-1.8132, -0.4862, -0.8410], [-1.8132, -0.4862, 0.8410]],
)

GLYCINE_ZWITTERION = (  # +NH3-CH2-COO-
    ["N", "C", "C", "O", "O", "H", "H", "H", "H", "H"],
    [[-1.2731, -0.7350, 0.0000], [0.0000, 0.0000, 0.0000], [1.1644, -0.9770, 0.0000],
     [2.3484, -0.5461, 0.0000], [0.9456, -2.2179, 0.0000], [0.0548, 0.6269, 0.8900],
     [0.0548, 0.6269, -0.8900], [-1.0848, -1.7477, 0.0000], [-1.8132, -0.4862, -0.8410],
     [-1.8132, -0.4862, 0.8410]],
)

GLYCINE_NEUTRAL = (  # H2N-CH2-COOH, the tautomer the zwitterion is compared against
    ["N", "C", "C", "O", "O", "H", "H", "H", "H", "H"],
    [[-1.2731, -0.7350, 0.0000], [0.0000, 0.0000, 0.0000], [1.1644, -0.9770, 0.0000],
     [2.3014, -0.5632, 0.0000], [0.9317, -2.2967, 0.0000], [1.7960, -2.7371, 0.0000],
     [0.0548, 0.6269, 0.8900], [0.0548, 0.6269, -0.8900], [-1.6231, -0.9371, 0.9256],
     [-1.6231, -0.9371, -0.9256]],
)

GLYCINE_ANION = (  # H2N-CH2-COO-
    ["N", "C", "C", "O", "O", "H", "H", "H", "H"],
    [[-1.2731, -0.7350, 0.0000], [0.0000, 0.0000, 0.0000], [1.1644, -0.9770, 0.0000],
     [2.3484, -0.5461, 0.0000], [0.9456, -2.2179, 0.0000], [0.0548, 0.6269, 0.8900],
     [0.0548, 0.6269, -0.8900], [-1.6231, -0.9371, 0.9256], [-1.6231, -0.9371, -0.9256]],
)

EXPERIMENT = {"acetic": 4.76, "glycine_pKa1": 2.34, "glycine_pKa2": 9.60, "glycine_pI": 5.97}

#: What a build without a piece says when asked for it. Reduced builds are
#: builds, not regressions, so these are skipped and named; any other failure is one.
NO_TBLITE = ("tblite", "not built", "build with")
NO_EXECUTABLE = ("none was found",)
NO_CREST = ("cannot sample conformers", "MQC_ENABLE_CREST")


def geometry(spec):
    return pka.Geometry(*spec)


def microstates():
    acetic = [
        pka.Microstate("acetic acid", geometry(ACETIC_ACID), charge=0, n_protons=1,
                       site_class="carboxylic"),
        pka.Microstate("acetate", geometry(ACETATE), charge=-1, n_protons=0),
    ]
    glycine = [
        pka.Microstate("glycine cation", geometry(GLYCINE_CATION), charge=+1, n_protons=2,
                       site_class=("carboxylic", "ammonium")),
        pka.Microstate("glycine neutral", geometry(GLYCINE_NEUTRAL), charge=0, n_protons=1,
                       site_class="carboxylic"),
        pka.Microstate("glycine zwitterion", geometry(GLYCINE_ZWITTERION), charge=0,
                       n_protons=1, site_class="ammonium"),
        pka.Microstate("glycine anion", geometry(GLYCINE_ANION), charge=-1, n_protons=0),
    ]
    return acetic, glycine


def tblite_available():
    """One xTB single point that either returns or names the missing option."""
    system = mqc.System(symbols=["O", "H", "H"],
                        coords=[[0.0, 0.0, 0.1], [0.0, 0.77, -0.47], [0.0, -0.77, -0.47]])
    system.set_monomers([[0, 1, 2]])
    try:
        mqc.MBE(system, level=0, method="gfn2", verbosity="error").run(
            label="probe_pka_tblite", write_to_file=False)
    except mqc.MQCError as exc:
        return not any(hint in str(exc) for hint in NO_TBLITE)
    return True


def finite(x):
    return isinstance(x, float) and math.isfinite(x)


def run_set(states, protocol):
    """`pka.run`, retried without the conformer search if the build has no CREST."""
    try:
        return pka.run(states, protocol, prefix="pka")
    except pka.PKaError as exc:
        if protocol.conformers is None or not any(hint in str(exc) for hint in NO_CREST):
            raise
        print("  this build has no CREST: running without a conformer search")
        protocol.conformers = None
        return pka.run(states, protocol, prefix="pka")


def check_result_basics(result, expected_names):
    assert {m.name for m in result.microstates} == set(expected_names)
    for m in result.microstates:
        assert finite(m.free_energy), f"{m.name}: free energy {m.free_energy}"
        assert m.detail["n_conformers"] >= 1
        weights = [c["weight"] for c in m.detail["conformers"]]
        assert abs(sum(weights) - 1.0) < 1e-9, f"{m.name}: conformer weights {weights}"
    # Net charge never rises with pH, over the whole range.
    charges = [result.charge(0.25 * i) for i in range(57)]
    assert all(b <= a + 1e-9 for a, b in zip(charges, charges[1:])), "charge rose with pH"
    for warning in result.warnings:
        print("  note:", warning)


def acetic_acid(protocol):
    print("\nacetic acid / acetate")
    acid_base = microstates()[0]
    result = run_set(acid_base, protocol)
    check_result_basics(result, ["acetic acid", "acetate"])

    (micro,), (macro,) = result.micro_pkas, result.macro_pkas
    # One microstate per level: the micro and macro pKa are the same quantity.
    assert abs(micro["pKa"] - macro["pKa"]) < 1e-9
    pops = result.populations(micro["pKa"])
    assert abs(pops["acetic acid"] - 0.5) < 1e-9 and abs(pops["acetate"] - 0.5) < 1e-9
    assert result.pI is None, "a monoprotic acid has charges 0 and -1: no isoelectric point"
    print(f"  uncalibrated pKa {micro['pKa']:7.2f}   experiment {EXPERIMENT['acetic']:.2f}")

    # The calibration is what makes the number usable, and reproduces its own datum.
    ref = pka.Reference("acetic acid", "acetate", EXPERIMENT["acetic"], site_class="carboxylic")
    calibrated = result.with_calibration([ref])
    assert abs(calibrated.macro_pkas[0]["pKa"] - EXPERIMENT["acetic"]) < 1e-9
    shift = calibrated.calibration.shifts["carboxylic"]
    print(f"  carboxylic proton shift {shift:+.2f} kcal/mol "
          f"({shift / (rt_ln10(result)):+.2f} pKa units)")
    return result, calibrated


def rt_ln10(result):
    return pka._pka_math.thermal_rt(result.temperature) * pka._pka_math.LN10


def glycine(protocol, acetic_calibration):
    print("\nglycine")
    states = microstates()[1]
    result = run_set(states, protocol)
    check_result_basics(result, [m.name for m in states])

    macro = result.macro_pkas
    assert [(d["from_n"], d["to_n"]) for d in macro] == [(2, 1), (1, 0)]
    print(f"  uncalibrated pKa1 {macro[0]['pKa']:7.2f}   pKa2 {macro[1]['pKa']:7.2f}   "
          f"pI {result.pI if result.pI is None else round(result.pI, 2)}")
    print(f"  experiment    pKa1 {EXPERIMENT['glycine_pKa1']:7.2f}   "
          f"pKa2 {EXPERIMENT['glycine_pKa2']:7.2f}   pI {EXPERIMENT['glycine_pI']}")
    if result.pI is not None:
        assert macro[0]["pKa"] <= result.pI <= macro[1]["pKa"]

    # Every micro pKa is consistent with one set of populations.
    for d in result.micro_pkas:
        pops = result.populations(d["pKa"])
        assert abs(pops[d["acid"]] - pops[d["base"]]) < 1e-9
    # And the two routes from cation to anion cost the same, whichever neutral
    # form they pass through: the calibration shifts attach to sites, not paths.
    by = {(d["acid"], d["base"]): d["pKa"] for d in result.micro_pkas}
    via_neutral = by[("glycine cation", "glycine neutral")] + by[("glycine neutral", "glycine anion")]
    via_zwit = by[("glycine cation", "glycine zwitterion")] + by[("glycine zwitterion", "glycine anion")]
    assert abs(via_neutral - via_zwit) < 1e-9

    # Transfer the carboxylic shift fitted on acetic acid; fit the ammonium one
    # to glycine's own pKa2 (one reference is a relative pKa). pKa1 is then a
    # prediction relative to acetic acid, and the pI follows from both.
    carboxylic = acetic_calibration.calibration.class_free_energy["carboxylic"]
    ammonium_ref = pka.Reference("glycine zwitterion", "glycine anion",
                                 EXPERIMENT["glycine_pKa2"], site_class="ammonium")
    cal = pka.calibrate(result.free_energies, [ammonium_ref], result.temperature)
    final = pka.PKaResult(result.microstates, result.temperature,
                          class_free_energy={"carboxylic": carboxylic, **cal.class_free_energy},
                          calibration=cal)
    pka1, pka2 = final.macro_pkas[0]["pKa"], final.macro_pkas[1]["pKa"]
    assert finite(pka1) and finite(pka2) and pka1 < pka2, (pka1, pka2)
    assert final.pI is not None and pka1 < final.pI < pka2, (pka1, final.pI, pka2)
    print(f"  calibrated    pKa1 {pka1:7.2f}   pKa2 {pka2:7.2f}   pI {final.pI:.2f}")
    return final


def main(argv):
    quick = "--quick" in argv
    with mqc.session():
        if not tblite_available():
            print("skipped: this build has no tblite")
            return 0
        try:
            pka.find_executable()
        except pka.PKaError as exc:
            if any(hint in str(exc) for hint in NO_EXECUTABLE):
                print("skipped: the optimize and conformers stages need the mqc executable "
                      "(set MQC_EXECUTABLE)")
                return 0
            raise
        protocol = pka.Protocol(workdir="pka_example_work", max_conformers=3,
                                conformers=None if quick else pka.Protocol().conformers)
        _, acetic_calibrated = acetic_acid(protocol)
        final = glycine(protocol, acetic_calibrated)
        # The JSON is the artefact a script would keep.
        assert "isoelectric_point" in final.to_json()
    print("\npKa workflow: identities hold")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
