"""Infrared and UV-vis spectra of water, from two runs of the executable.

    MQC_EXECUTABLE=build/mqc python spectra.py

A Hessian run writes the frequencies and infrared intensities, an energy run
with `keywords.excited_states` writes the excitation energies and oscillator
strengths, and `mqc.spectra` turns either into a curve. The curves are written
beside the outputs as CSV, and as figures when matplotlib is installed.

Run as a test rather than a demonstration: it asserts what must hold of any
spectrum made this way, not agreement with experiment, and exits non-zero if
any of it fails.

  * every vibration the run reports is a stick, and the six modes near zero
    are not;
  * the area under the infrared curve is the summed intensity;
  * the area under the UV-vis curve is the summed oscillator strength.

With `properties.raman` on the Hessian deck the same output also carries Raman
activities, and `spectra.raman(output)` reads them the same way.
"""

import json
import math
import os
import subprocess
import sys
import tempfile

from mqc import spectra

WATER = {
    "symbols": ["O", "H", "H"],
    "geometry": [0.0, 0.0, 0.1173, 0.0, 0.7572, -0.4692, 0.0, -0.7572, -0.4692],
    "molecular_charge": 0,
    "molecular_multiplicity": 1,
}


def run(folder, name, driver, keywords):
    deck = {
        "schema": {"name": name, "version": "1.0"},
        "molecules": [WATER],
        "model": {"method": "hf", "basis": "6-31g"},
        "driver": driver,
        "keywords": keywords,
    }
    path = os.path.join(folder, name + ".json")
    with open(path, "w") as handle:
        json.dump(deck, handle)
    executable = os.environ.get("MQC_EXECUTABLE", "mqc")
    done = subprocess.run([executable, path], cwd=folder, capture_output=True, text=True)
    output = os.path.join(folder, "output_%s.json" % name)
    if done.returncode != 0 or not os.path.exists(output):
        sys.exit("the %s run failed:\n%s" % (name, (done.stdout + done.stderr)[-2000:]))
    return output


def area(x, y):
    return sum(0.5 * (y[k] + y[k + 1]) * (x[k + 1] - x[k]) for k in range(len(x) - 1))


def main():
    folder = tempfile.mkdtemp(prefix="mqc_spectra_")

    hessian = run(folder, "water_hessian", "Hessian", {})
    ir = spectra.ir(hessian, fwhm=20.0, lineshape="gaussian", scale=0.90)
    print("infrared, frequencies scaled by 0.90:")
    for position, intensity in ir.sticks:
        print("  %8.1f cm^-1   %8.2f km/mol" % (position, intensity))
    assert len(ir.sticks) == 3, "water has three vibrations, got %d" % len(ir.sticks)
    assert all(position > 1000.0 for position, _ in ir.sticks)
    summed = sum(intensity for _, intensity in ir.sticks)
    recovered = area(ir.x, ir.y) * math.log(10.0) / 100.0
    assert math.isclose(recovered, summed, rel_tol=1.0e-4), (recovered, summed)
    ir.to_csv(os.path.join(folder, "water_ir.csv"))

    excited = run(folder, "water_excited", "Energy", {"excited_states": {"n_states": 6}})
    uv = spectra.uvvis(excited, fwhm_ev=0.4, unit="cm-1", points=20000)
    print("UV-vis:")
    for position, strength in spectra.uvvis(excited, unit="nm").sticks:
        print("  %8.1f nm      f = %.4f" % (position, strength))
    summed = sum(strength for _, strength in uv.sticks)
    recovered = area(uv.x, uv.y) * spectra.F_PER_AREA
    assert summed > 0.0, "no state of water absorbs"
    assert math.isclose(recovered, summed, rel_tol=1.0e-3), (recovered, summed)
    uv.to_csv(os.path.join(folder, "water_uvvis.csv"))

    try:
        ir.plot(os.path.join(folder, "water_ir.png"))
        spectra.uvvis(excited, unit="nm").plot(os.path.join(folder, "water_uvvis.png"))
    except spectra.SpectrumError as problem:
        print("no figures: %s" % problem)
    for note in ir.warnings + uv.warnings:
        print("note: %s" % note)
    print("written to %s" % folder)


if __name__ == "__main__":
    main()
