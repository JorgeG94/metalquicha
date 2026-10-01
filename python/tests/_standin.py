"""A stand-in for the `mqc` package, for the workflow tests that need no build.

`fake_package` imports one workflow module (`mqc.pka`, `mqc.bde`) under a package
whose `MBE` and `_check_label` are the *real* ones, lifted out of
`mqc/__init__.py` by parsing it, with only the calls that reach Fortran replaced
by an invented physics. So the decks the workflow writes, the settings it passes
and the labels it chooses are the real thing; the numbers are not.

The thermochemistry block it returns mimics the Fortran one's known gap: it
reports ``spin_multiplicity: 1`` and a zero electronic entropy whatever the
system's multiplicity, because that is what `compute_thermochemistry` does when
no caller passes the multiplicity.
"""

import ast
import contextlib
import importlib
import json
import os
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(HERE, "..", "mqc")


def lift_from_init():
    """The real `MBE`, `_check_label`, `_merge` and `DRIVERS`, without `_ffi`."""
    with open(os.path.join(PKG, "__init__.py")) as handle:
        tree = ast.parse(handle.read())
    keep = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", "") == "DRIVERS" for t in node.targets):
            keep.append(node)
        elif isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in (
            "MBE",
            "_check_label",
            "_merge",
            "MQCError",
        ):
            keep.append(node)
    namespace = {"json": json, "os": os, "ctypes": None, "_ffi": None, "HARTREE_TO_EV": 27.211386245988}
    exec(compile(ast.Module(body=keep, type_ignores=[]), "mqc/__init__.py (lifted)", "exec"), namespace)
    return namespace


def fake_thermo(temperature=298.15, pressure=1.0, multiplicity_reported=1):
    return {
        "temperature_K": temperature,
        "pressure_atm": pressure,
        "is_linear": False,
        "spin_multiplicity": multiplicity_reported,
        "contributions": {
            "translational": {"energy_hartree": 0.0014164, "entropy_cal_mol_K": 34.6},
            "rotational": {"energy_hartree": 0.0014164, "entropy_cal_mol_K": 10.5},
            "vibrational": {"energy_hartree": 0.0, "entropy_cal_mol_K": 0.0},
            "electronic": {"energy_hartree": 0.0, "entropy_cal_mol_K": 0.0},
        },
    }


class FakeResult:
    def __init__(self, energy, label, freqs=None, thermo=None):
        self.energy = energy
        self.label = label
        self.frequencies = freqs
        self.thermochemistry = thermo


def default_physics(driver, system, label, n):
    """What `mqc.pka`'s tests have always used: a drifting energy, 3N-6 real modes."""
    energy = -50.0 - 0.001 * (n % 3)
    if driver == "hessian":
        n_atoms = system.n_atoms
        freqs = [0.0, 0.1, -0.1, 0.2, 0.3, -0.2] + [300.0 + 100.0 * i for i in range(3 * n_atoms - 6)]
        if "imag" in label:
            freqs[6] = -80.0
        thermo = fake_thermo()
        thermo["total_energies_hartree"] = {"electronic": energy}
        return FakeResult(energy, label, freqs, thermo)
    return FakeResult(energy, label)


@contextlib.contextmanager
def fake_package(log, physics=None, module="pka"):
    """`mqc` with the real settings and label check, and invented physics.

    ``physics(driver, system, label, n)`` returns a `FakeResult`; ``n`` counts the
    runs. Every ``mqc`` module is removed from ``sys.modules`` on the way in and
    restored on the way out, so one test's import does not leak into the next.
    """
    physics = physics or default_physics
    ns = lift_from_init()
    pkg = types.ModuleType("mqc")
    pkg.__path__ = [PKG]
    pkg._check_label = ns["_check_label"]
    pkg.MQCError = ns["MQCError"]
    real_mbe = ns["MBE"]
    counter = {"n": 0}

    class System:
        def __init__(self, symbols=None, coords=None, charge=0, multiplicity=1):
            self._handle = object()
            self.symbols, self.coords, self.charge, self.multiplicity = symbols, coords, charge, multiplicity
            self.n_atoms = len(symbols)

        def set_monomers(self, monomers, charges=None, multiplicities=None):
            self.monomers, self.charges, self.multiplicities = monomers, charges, multiplicities
            return self

    class MBE(real_mbe):
        def run(self, label="mqc", write_to_file=True):
            pkg._check_label(label)  # the real one: a bad label fails here as it would there
            counter["n"] += 1
            log.append({"label": label, "driver": self.driver, "kwargs": self.settings(), "system": self.system})
            return physics(str(self.driver).lower(), self.system, label, counter["n"])

    pkg.System, pkg.MBE, pkg.DRIVERS = System, MBE, ns["DRIVERS"]
    saved = {k: v for k, v in sys.modules.items() if k == "mqc" or k.startswith("mqc.")}
    for k in saved:
        del sys.modules[k]
    sys.modules["mqc"] = pkg
    try:
        loaded = importlib.import_module(f"mqc.{module}")
        setattr(pkg, module, loaded)
        yield loaded
    finally:
        for k in [k for k in sys.modules if k == "mqc" or k.startswith("mqc.")]:
            del sys.modules[k]
        sys.modules.update(saved)
