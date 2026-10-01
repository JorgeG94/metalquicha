"""pKa and isoelectric point from microstate free energies.

    import mqc
    from mqc import pka

    with mqc.session():
        acid = pka.Microstate("AcOH", pka.Geometry.from_xyz("acoh.xyz"), charge=0, n_protons=1)
        base = pka.Microstate("AcO-", pka.Geometry.from_xyz("aco.xyz"), charge=-1, n_protons=0)
        result = pka.run([acid, base])
        result = result.with_calibration([pka.Reference("AcOH", "AcO-", 4.76)])
        print(result.macro_pkas, result.pI)

Per microstate: a conformer search, an optimization of the survivors, a Hessian,
a single point, and one free energy ``G = E_sp + G_corr`` with the conformers
Boltzmann-combined. The microstates are *given*, not enumerated: which
protonation states and tautomers exist is a chemical judgement, and the free
energies are only as good as the list. See `mqc_docs/source/pka.rst` for the
protocol, the calibration and what to distrust.

**Run it inside ``mqc.session()``**, which only rank 0 leaves, so everything
here -- including the subprocesses below -- happens once. The calls are made one
after another; the Hessian and single point use every rank of the session.

**Two of the four stages do not go through the session.** The library refuses
``driver="optimize"`` and ``driver="conformers"`` through the C interface (both
drive ``run_calculation`` rather than being driven by it), so those two run the
``mqc`` executable on a deck written to ``Protocol.workdir``. That is a real
limitation and the plumbing for it is the least pleasant part of this module;
see the ``TODO(mqc)`` notes. The executable is ``Protocol.executable``, else
``$MQC_EXECUTABLE``, else ``mqc`` on the PATH, else ``build/mqc`` in the source
tree. MPI's launcher variables are removed from its environment: it is a
separate single-rank program, and CREST refuses to sample on more than one.

The formulas -- quasi-RRHO, the standard state, Boltzmann combination, the
protonation model, calibration, the result -- are in `mqc._pka_math`, which
imports nothing from the library and is what the tests exercise.
"""

import copy
import dataclasses
import glob
import os
import re
import shutil
import subprocess

from . import _pka_math
from ._pka_math import (  # noqa: F401  (re-exported: this is the public surface)
    Calibration,
    Geometry,
    MicrostateData,
    NoIsoelectricPoint,
    PKaError,
    PKaModel,
    PKaResult,
    Reference,
    boltzmann_combine,
    calibrate,
    default_proton_free_energy,
    free_energy_correction,
    qrrho_entropy_over_r,
    qrrho_weight,
    standard_state_correction,
)

__all__ = [
    "Microstate",
    "Protocol",
    "Geometry",
    "Reference",
    "PKaResult",
    "PKaModel",
    "PKaError",
    "NoIsoelectricPoint",
    "Calibration",
    "run",
    "evaluate_microstate",
    "find_executable",
    "calibrate",
    "default_proton_free_energy",
    "standard_state_correction",
    "free_energy_correction",
    "qrrho_entropy_over_r",
    "qrrho_weight",
]

#: GFN2-xTB with ALPB water. One dict per stage, so that a stage can be given a
#: different level of theory without the others noticing: swapping the
#: ``single_point`` for DFT is a change to one entry.
_XTB_WATER = {
    "method": "gfn2",
    "xtb": {"solvent": "water", "solvation_model": "alpb"},
    "verbosity": "error",
}

#: Keys of an `MBE` call that the workflow owns. A stage dict naming one would
#: be overridden or, worse, silently honoured in one stage and not another.
_RESERVED = ("system", "level", "driver")


def _default_stage():
    return copy.deepcopy(_XTB_WATER)


@dataclasses.dataclass
class Protocol:
    """What to run at each stage, and how to turn it into a free energy.

    The four stage dicts are keyword arguments of `mqc.MBE`; the workflow
    supplies the system, ``level=0`` and the driver. ``conformers`` and
    ``optimize`` may be None to skip the stage and use the structure as given.

    ``conformers`` is the *refinement* level CREST re-ranks its ensemble with;
    CREST samples with its own GFN2 and takes the solvent from ``xtb`` (only
    ALPB and GBSA with a named solvent can be passed to it).

    ``energy_window_kcal`` and ``max_conformers`` choose which of the ensemble
    is carried forward, by the refinement energy and before the free energies
    are known. A window narrower than the frequency-dependent part of the free
    energy can move is a choice to drop conformers that would have mattered.

    ``qrrho``, ``nu0``, ``alpha`` and ``b_av`` are Grimme's quasi-RRHO
    settings; ``standard_state`` adds the 1 atm -> 1 M term;
    ``imaginary`` is ``"drop"`` or ``"flip"`` (see `free_energy_correction`).
    """

    conformers: dict = dataclasses.field(default_factory=_default_stage)
    optimize: dict = dataclasses.field(default_factory=_default_stage)
    frequencies: dict = dataclasses.field(default_factory=_default_stage)
    single_point: dict = dataclasses.field(default_factory=_default_stage)
    energy_window_kcal: float = 3.0
    max_conformers: int = 5
    qrrho: bool = True
    nu0: float = _pka_math.QRRHO_NU0
    alpha: float = _pka_math.QRRHO_ALPHA
    b_av: float = _pka_math.QRRHO_B_AV
    standard_state: bool = True
    imaginary: str = "drop"
    allow_unconverged: bool = False
    workdir: str = "pka_work"
    executable: str = None
    subprocess_env: dict = None

    def __post_init__(self):
        for stage in ("conformers", "optimize", "frequencies", "single_point"):
            block = getattr(self, stage)
            if block is None:
                if stage in ("frequencies", "single_point"):
                    raise ValueError(f"the {stage} stage cannot be skipped")
                continue
            block = dict(block)
            for key in _RESERVED:
                if key in block:
                    raise ValueError(
                        f"{stage}: {key!r} is set by the workflow, not by the protocol"
                    )
            setattr(self, stage, block)
        if self.max_conformers < 1:
            raise ValueError("max_conformers must be at least 1")
        if self.imaginary not in ("drop", "flip"):
            raise ValueError("imaginary must be 'drop' or 'flip'")

    def as_dict(self):
        """The protocol as plain data, for the result to carry."""
        return dataclasses.asdict(self)

    def single_point_is_frequencies(self):
        """Whether the single point is the same calculation the Hessian already did."""
        return self.single_point == self.frequencies


class Microstate:
    """One protonation state or tautomer to be computed.

    ``system`` is a `Geometry` (preferred), an .xyz path, a ``(symbols,
    coords)`` pair, or an ``mqc.System``. A ``System`` cannot be read back --
    it has an atom count and no coordinates -- so it works only with a protocol
    that skips both the conformer search and the optimization.

    ``charge`` is the net charge and ``n_protons`` the number of *ionizable*
    protons it carries (acetic acid 1, acetate 0; glycine cation 2, both
    neutral forms 1, anion 0). Charge minus ``n_protons`` must be the same for
    every microstate of one molecule.

    ``site_class`` says what kind of site each of those protons is, for
    calibration: None (all one class), a string (every proton), or one label per
    proton. Tautomers at the same level differ in it -- the neutral glycine
    carries its proton on the carboxyl, the zwitterion on the amine -- and
    that is what lets the carboxyl and ammonium sites be calibrated apart.
    """

    def __init__(self, name, system, charge, n_protons, site_class=None, multiplicity=1):
        if not str(name).strip():
            raise ValueError("a microstate needs a name")
        self.name = str(name)
        self.charge = int(charge)
        self.n_protons = int(n_protons)
        self.multiplicity = int(multiplicity)
        self.site_classes = _pka_math.normalize_site_classes(site_class, self.n_protons)
        self._system = None
        self._geometry = None
        if isinstance(system, Geometry):
            self._geometry = system
        elif isinstance(system, (str, os.PathLike)):
            self._geometry = Geometry.from_xyz(system)
        elif isinstance(system, (tuple, list)) and len(system) == 2:
            self._geometry = Geometry(*system)
        elif hasattr(system, "_handle"):
            self._system = system
        else:
            raise TypeError(
                "system must be a Geometry, an .xyz path, (symbols, coords) or an mqc.System"
            )

    @property
    def geometry(self):
        if self._geometry is None:
            raise PKaError(
                f"{self.name}: this microstate was given as an mqc.System, which cannot "
                "be read back, so it cannot be searched or optimized. Give a Geometry, or "
                "set Protocol(conformers=None, optimize=None)."
            )
        return self._geometry

    @property
    def n_atoms(self):
        return self._geometry.n_atoms if self._geometry is not None else self._system.n_atoms

    def to_system(self, geometry=None):
        """An ``mqc.System`` for this microstate, optionally at another geometry."""
        if geometry is None and self._system is not None:
            return self._system
        from . import System

        geometry = geometry or self.geometry
        system = System(
            symbols=geometry.symbols,
            coords=geometry.coords,
            charge=self.charge,
            multiplicity=self.multiplicity,
        )
        # The whole molecule as one monomer, as every level-0 run here does, with
        # the charge on the monomer as well as on the system.
        system.set_monomers(
            [list(range(geometry.n_atoms))],
            charges=[self.charge],
            multiplicities=[self.multiplicity],
        )
        return system

    def __repr__(self):
        return f"<Microstate {self.name} q={self.charge:+d} n_H={self.n_protons}>"


# ---------------------------------------------------------------------------
#  The stages that run the executable
# ---------------------------------------------------------------------------


def _slug(name):
    """A name made safe for a label and a directory: no dots, slashes or spaces."""
    return re.sub(r"[^A-Za-z0-9_-]", "_", str(name))


def _find_executable(protocol):
    candidates = [protocol.executable, os.environ.get("MQC_EXECUTABLE"), shutil.which("mqc")]
    here = os.path.dirname(os.path.abspath(__file__))
    for up in ("..", "../..", "../../.."):
        candidates.append(os.path.join(here, up, "build", "mqc"))
    for path in candidates:
        if path and os.path.isfile(path) and os.access(path, os.X_OK):
            return os.path.abspath(path)
    raise PKaError(
        "the conformer and optimization stages run the mqc executable, and none was found. "
        "Set Protocol(executable=...) or MQC_EXECUTABLE, or skip those stages with "
        "Protocol(conformers=None, optimize=None)."
    )


def find_executable(protocol=None):
    """The ``mqc`` executable the conformer and optimization stages will run.

    Raises `PKaError` naming how to point at one when there is none -- which is
    also the way to ask in advance whether those stages can run at all.
    """
    return _find_executable(protocol or Protocol())


#: Variables an MPI launcher sets that would make a child `mqc` try to join the
#: parent's job rather than start as the single rank it must be.
_MPI_PREFIXES = ("OMPI_", "PMIX_", "PMI_", "HYDRA_", "I_MPI_", "MPICH_")


def _child_env(protocol):
    env = {k: v for k, v in os.environ.items() if not k.startswith(_MPI_PREFIXES)}
    if protocol.subprocess_env:
        env.update(protocol.subprocess_env)
    return env


def _deck(stage, driver, geometry, microstate, directory):
    """Write ``start.xyz`` and a deck for ``driver`` beside it; return the deck name.

    The deck is what `mqc.MBE` would send, less the fragmentation block (a
    deck for an optimization has none, and ``level=0`` is a spelling of "not
    fragmented" that only the in-process path needs) plus the molecule.
    """
    import json

    from . import MBE

    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, "start.xyz"), "w") as handle:
        handle.write(geometry.to_xyz(microstate.name))
    document = MBE(None, level=0, driver=driver, **stage).settings()
    document["keywords"].pop("fragmentation", None)
    document["molecules"] = [
        {
            "xyz": "start.xyz",
            "molecular_charge": microstate.charge,
            "molecular_multiplicity": microstate.multiplicity,
        }
    ]
    name = f"{driver}.json"
    with open(os.path.join(directory, name), "w") as handle:
        json.dump(document, handle, indent=2)
    return name


def _run_mqc(executable, deck, directory, env):
    """Run the executable on a deck; raise with its output if it fails."""
    done = subprocess.run(
        [executable, deck], cwd=directory, env=env, capture_output=True, text=True
    )
    if done.returncode != 0:
        tail = "\n".join((done.stdout + "\n" + done.stderr).strip().splitlines()[-15:])
        raise PKaError(f"mqc exited {done.returncode} on {directory}/{deck}:\n{tail}")


# TODO(mqc): the conformer ensemble reaches Python only as the files CREST leaves
# in the working directory (`crest_conformers.xyz`, energies on the comment line),
# and only from a deck run through the executable: `mqc_session` and the C API
# refuse the `conformers` driver, and `run_conformer_search` returns no ensemble
# to its caller. Same for `optimize`, whose result is `output_<deck>_optimized.xyz`.
# Both belong behind the C interface, returning energies and geometries.
def sample_conformers(microstate, protocol, directory, executable, env):
    """CREST ensemble for one microstate: ``[(energy_hartree, Geometry), ...]``, lowest first."""
    name = _deck(protocol.conformers, "conformers", microstate.geometry, microstate, directory)
    _run_mqc(executable, name, directory, env)
    path = os.path.join(directory, "crest_conformers.xyz")
    if not os.path.exists(path):
        raise PKaError(f"{microstate.name}: the conformer search left no {path}")
    with open(path) as handle:
        ensemble = _pka_math.parse_xyz_ensemble(handle.read())
    if not ensemble or any(e is None for e, _ in ensemble):
        raise PKaError(f"{path}: could not read an energy from every structure's comment line")
    ensemble.sort(key=lambda item: item[0])
    return ensemble


def optimize_structure(microstate, geometry, protocol, directory, executable, env):
    """Optimize one geometry: ``(Geometry, converged, energy_hartree)``."""
    name = _deck(protocol.optimize, "optimize", geometry, microstate, directory)
    _run_mqc(executable, name, directory, env)
    path = os.path.join(directory, f"output_{name[:-5]}_optimized.xyz")
    if not os.path.exists(path):
        # `optimized_xyz_path` derives the name from the output document's, and
        # a rename on the Fortran side should not read as a failed optimization.
        found = sorted(glob.glob(os.path.join(directory, "*_optimized.xyz")))
        if not found:
            raise PKaError(f"{microstate.name}: the optimization wrote no *_optimized.xyz in {directory}")
        path = found[0]
    with open(path) as handle:
        return _pka_math.parse_optimized_xyz(handle.read())


# ---------------------------------------------------------------------------
#  The stages that run in the session
# ---------------------------------------------------------------------------


def _frequency_stage(microstate, geometry, protocol, label):
    """Hessian, then the correction. Returns ``(correction, electronic_energy_hartree)``."""
    from . import MBE

    system = microstate.to_system(geometry)
    result = MBE(system, level=0, driver="hessian", **protocol.frequencies).run(label)
    freqs, thermo = result.frequencies, result.thermochemistry
    if freqs is None or thermo is None:
        raise PKaError(
            f"{label}: the Hessian run wrote no vibrational_analysis/thermochemistry block"
        )
    correction = free_energy_correction(
        freqs,
        microstate.n_atoms,
        thermo,
        nu0=protocol.nu0,
        alpha=protocol.alpha,
        b_av=protocol.b_av,
        standard_state=protocol.standard_state,
        imaginary=protocol.imaginary,
        qrrho=protocol.qrrho,
    )
    return correction, float(thermo["total_energies_hartree"]["electronic"])


def _single_point(microstate, geometry, protocol, label):
    from . import MBE

    system = microstate.to_system(geometry)
    return MBE(system, level=0, driver="energy", **protocol.single_point).run(label).energy


def evaluate_microstate(microstate, protocol=None, prefix="pka", verbose=True, _executable=None):
    """One microstate to one free energy: `MicrostateData` with per-conformer detail.

    Search, optimize, Hessian, single point, then Boltzmann-combine. Every
    conformer is kept in ``detail["conformers"]`` with its energies, correction,
    imaginary frequencies and final weight; ``detail["n_imaginary"]`` is the
    total across conformers. Imaginary frequencies are counted here and
    surfaced by `PKaResult.warnings`, not dropped.
    """
    protocol = protocol or Protocol()
    from . import _check_label

    _check_label(prefix)
    slug = _slug(microstate.name)
    root = os.path.join(protocol.workdir, slug)
    say = (lambda text: print(f"  [{microstate.name}] {text}", flush=True)) if verbose else (lambda t: None)

    need_exe = protocol.conformers is not None or protocol.optimize is not None
    executable = _executable or (_find_executable(protocol) if need_exe else None)
    env = _child_env(protocol)

    # 1. Structures to carry forward.
    if protocol.conformers is not None:
        say("conformer search")
        ensemble = sample_conformers(
            microstate, protocol, os.path.join(root, "conformers"), executable, env
        )
        lowest = ensemble[0][0]
        window = protocol.energy_window_kcal / _pka_math.HARTREE_TO_KCAL
        kept = [item for item in ensemble if item[0] - lowest <= window][: protocol.max_conformers]
        say(f"{len(ensemble)} conformers, {len(kept)} within {protocol.energy_window_kcal:g} kcal/mol")
    else:
        kept = [(None, microstate.geometry if microstate._geometry is not None else None)]

    # 2. Optimize them.
    structures = []
    for k, (e_conf, geom) in enumerate(kept):
        converged, e_opt = True, None
        if protocol.optimize is not None:
            say(f"optimizing conformer {k}")
            geom, converged, e_opt = optimize_structure(
                microstate, geom, protocol, os.path.join(root, f"opt_c{k}"), executable, env
            )
            if not converged and not protocol.allow_unconverged:
                raise PKaError(
                    f"{microstate.name}: the optimization of conformer {k} did not converge. "
                    "A Hessian there is not a minimum's; raise optimization steps or set "
                    "Protocol(allow_unconverged=True) to carry on and be warned."
                )
        structures.append({"k": k, "geometry": geom, "e_conf": e_conf, "e_opt": e_opt, "converged": converged})

    # Two conformers that optimized to the same minimum would be counted twice in
    # the Boltzmann sum, biasing the free energy down by RT ln 2.
    unique = []
    for s in structures:
        if s["e_opt"] is not None and any(
            u["e_opt"] is not None and abs(u["e_opt"] - s["e_opt"]) < 1.0e-6 for u in unique
        ):
            say(f"conformer {s['k']} optimized to the same energy as an earlier one; dropped")
            continue
        unique.append(s)

    # 3. Hessian and single point.
    records = []
    for s in unique:
        label = f"{prefix}_{slug}_c{s['k']}"
        _check_label(label + "_freq")
        say(f"frequencies, conformer {s['k']}")
        correction, e_freq = _frequency_stage(microstate, s["geometry"], protocol, label + "_freq")
        if protocol.single_point_is_frequencies():
            e_sp, reused = e_freq, True
        else:
            say(f"single point, conformer {s['k']}")
            e_sp, reused = _single_point(microstate, s["geometry"], protocol, label + "_sp"), False
        g = e_sp * _pka_math.HARTREE_TO_KCAL + correction["g_corr"]
        records.append(
            {
                "conformer": s["k"],
                "e_conformer_hartree": s["e_conf"],
                "e_optimized_hartree": s["e_opt"],
                "e_single_point_hartree": e_sp,
                "single_point_reused_frequency_energy": reused,
                "g_corr_kcal_mol": correction["g_corr"],
                "g_kcal_mol": g,
                "n_imaginary": correction["n_imaginary"],
                "imaginary_cm1": correction["imaginary_cm1"],
                "tr_max_cm1": correction["tr_max_cm1"],
                "optimization_converged": s["converged"],
                "temperature_K": correction["temperature_K"],
            }
        )

    temperature = records[0]["temperature_K"]
    weights = _pka_math.boltzmann_weights([r["g_kcal_mol"] for r in records], temperature)
    for r, w in zip(records, weights):
        r["weight"] = w
    g_micro = boltzmann_combine([r["g_kcal_mol"] for r in records], temperature)
    detail = {
        "n_conformers": len(records),
        "n_imaginary": sum(r["n_imaginary"] for r in records),
        "temperature_K": temperature,
        "conformers": records,
    }
    say(f"G = {g_micro:.3f} kcal/mol from {len(records)} conformer(s)")
    return MicrostateData(
        microstate.name, microstate.charge, microstate.n_protons, g_micro, microstate.site_classes, detail
    )


def run(microstates, protocol=None, references=None, prefix="pka", verbose=True):
    """Compute every microstate and return a `PKaResult`.

    ``references`` is an optional list of `Reference` (or ``(acid, base, pKa)``
    / ``(acid, base, pKa, site_class)`` tuples); with it the result is
    calibrated, without it the proton free energy is the literature value and
    the absolute pKas are not to be trusted. `PKaResult.with_calibration` does
    the same afterwards, which is the better order: the microstates are the
    expensive part and the reference pKas are a decision.

    Inside ``mqc.session()``. Sequential, rank 0.
    """
    protocol = protocol or Protocol()
    microstates = list(microstates)
    if not microstates:
        raise ValueError("no microstates")
    names = [m.name for m in microstates]
    if len(set(names)) != len(names) or len({_slug(n) for n in names}) != len(names):
        raise ValueError("microstate names must be unique, and unique once made filename-safe")
    # Charge and proton count are checked before anything is run: a bad list
    # found after an hour of CREST is the expensive way to learn it.
    _pka_math.PKaModel(
        [MicrostateData(m.name, m.charge, m.n_protons, 0.0, m.site_classes) for m in microstates]
    )
    data = []
    for m in microstates:
        if verbose:
            print(f"microstate {m.name}", flush=True)
        data.append(evaluate_microstate(m, protocol, prefix, verbose))
    temps = {d.detail["temperature_K"] for d in data}
    if len(temps) != 1:
        raise PKaError(f"the microstates were evaluated at different temperatures: {sorted(temps)}")
    result = PKaResult(data, temps.pop(), protocol=protocol.as_dict())
    if references:
        result = result.with_calibration(references)
    return result
