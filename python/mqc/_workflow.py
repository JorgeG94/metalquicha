"""The staged evaluation of one species, shared by `mqc.pka` and `mqc.bde`.

One species (a protonation microstate, a molecule, a radical fragment) goes
through a conformer search, an optimization of the survivors, a Hessian and a
single point, and comes out as an energy, an enthalpy and a free energy with the
conformers Boltzmann-combined. What a species *means* -- which proton it carries,
which bond it came from -- is the caller's; this module knows a name, a charge, a
multiplicity and a geometry.

**Two of the four stages do not go through the session.** The library refuses
``driver="optimize"`` and ``driver="conformers"`` through the C interface (both
drive ``run_calculation`` rather than being driven by it), so those two run the
``mqc`` executable on a deck written to ``Protocol.workdir``; see the
``TODO(mqc)`` notes. The executable is ``Protocol.executable``, else
``$MQC_EXECUTABLE``, else ``mqc`` on the PATH, else ``build/mqc`` in the source
tree. MPI's launcher variables are removed from its environment: it is a
separate single-rank program, and CREST refuses to sample on more than one.

A single atom skips all four: it has no conformers, nothing to optimize and no
vibrations, so it is one single point plus the analytic thermochemistry of
`mqc._thermo.atom_correction`.

The formulas are in `mqc._thermo`, which imports nothing from the library.
"""

import copy
import dataclasses
import glob
import json
import os
import re
import shutil
import subprocess

from . import _thermo
from ._thermo import Geometry, WorkflowError

#: Keys of an `MBE` call that the workflow owns. A stage dict naming one would
#: be overridden or, worse, silently honoured in one stage and not another.
_RESERVED = ("system", "level", "driver")

_STAGES = ("conformers", "optimize", "frequencies", "single_point")


def _gas_stage():
    return {"method": "gfn2", "verbosity": "error"}


@dataclasses.dataclass
class StageProtocol:
    """What to run at each stage, and how to turn it into a free energy.

    The four stage dicts are keyword arguments of `mqc.MBE`; the workflow
    supplies the system, ``level=0`` and the driver. ``conformers`` and
    ``optimize`` may be None to skip the stage and use the structure as given.
    A workflow subclasses this to choose its own defaults (`mqc.pka.Protocol`:
    GFN2 with ALPB water, `mqc.bde.Protocol`: GFN2 in the gas phase).

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

    The temperature and pressure of the thermochemistry are the ``hessian``
    block of the ``frequencies`` stage (``{"hessian": {"temperature": 310.0}}``),
    298.15 K and 1 atm when it names none.
    """

    conformers: dict = dataclasses.field(default_factory=_gas_stage)
    optimize: dict = dataclasses.field(default_factory=_gas_stage)
    frequencies: dict = dataclasses.field(default_factory=_gas_stage)
    single_point: dict = dataclasses.field(default_factory=_gas_stage)
    energy_window_kcal: float = 3.0
    max_conformers: int = 5
    qrrho: bool = True
    nu0: float = _thermo.QRRHO_NU0
    alpha: float = _thermo.QRRHO_ALPHA
    b_av: float = _thermo.QRRHO_B_AV
    standard_state: bool = True
    imaginary: str = "drop"
    allow_unconverged: bool = False
    workdir: str = "mqc_work"
    executable: str = None
    subprocess_env: dict = None

    def __post_init__(self):
        for stage in _STAGES:
            block = getattr(self, stage)
            if block is None:
                if stage in ("frequencies", "single_point"):
                    raise ValueError(f"the {stage} stage cannot be skipped")
                continue
            block = dict(block)
            for key in _RESERVED:
                if key in block:
                    raise ValueError(f"{stage}: {key!r} is set by the workflow, not by the protocol")
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

    def conditions(self):
        """``(temperature_K, pressure_atm)`` the thermochemistry is evaluated at."""
        block = (self.frequencies or {}).get("hessian") or {}
        return float(block.get("temperature", 298.15)), float(block.get("pressure", 1.0))


def default_stage(**overrides):
    """A fresh GFN2 gas-phase stage dict, with ``overrides`` merged in."""
    stage = _gas_stage()
    stage.update(copy.deepcopy(overrides))
    return stage


class Species:
    """One thing to be computed: a name, a charge, a multiplicity, a structure.

    ``system`` is a `Geometry` (preferred), an .xyz path, a ``(symbols,
    coords)`` pair, or an ``mqc.System``. A ``System`` cannot be read back --
    it has an atom count and no coordinates -- so it works only with a protocol
    that skips both the conformer search and the optimization, and not for a
    single atom (whose element is needed).
    """

    #: What a subclass calls itself in an error message.
    kind = "species"

    def __init__(self, name, system, charge=0, multiplicity=1):
        if not str(name).strip():
            raise ValueError(f"a {self.kind} needs a name")
        self.name = str(name)
        self.charge = int(charge)
        self.multiplicity = int(multiplicity)
        if self.multiplicity < 1:
            raise ValueError("multiplicity must be a positive integer")
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
            raise TypeError("system must be a Geometry, an .xyz path, (symbols, coords) or an mqc.System")

    @property
    def geometry(self):
        if self._geometry is None:
            raise WorkflowError(
                f"{self.name}: this {self.kind} was given as an mqc.System, which cannot "
                "be read back, so it cannot be searched or optimized. Give a Geometry, or "
                "set Protocol(conformers=None, optimize=None)."
            )
        return self._geometry

    @property
    def n_atoms(self):
        return self._geometry.n_atoms if self._geometry is not None else self._system.n_atoms

    def to_system(self, geometry=None):
        """An ``mqc.System`` for this species, optionally at another geometry."""
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
        return f"<{type(self).__name__} {self.name} q={self.charge:+d} m={self.multiplicity}>"


# ---------------------------------------------------------------------------
#  The stages that run the executable
# ---------------------------------------------------------------------------


def slug(name):
    """A name made safe for a label and a directory: no dots, slashes or spaces."""
    return re.sub(r"[^A-Za-z0-9_-]", "_", str(name))


def find_executable(protocol):
    """The ``mqc`` executable the conformer and optimization stages run.

    Raises `WorkflowError` naming how to point at one when there is none.
    """
    candidates = [protocol.executable, os.environ.get("MQC_EXECUTABLE"), shutil.which("mqc")]
    here = os.path.dirname(os.path.abspath(__file__))
    for up in ("..", "../..", "../../.."):
        candidates.append(os.path.join(here, up, "build", "mqc"))
    for path in candidates:
        if path and os.path.isfile(path) and os.access(path, os.X_OK):
            return os.path.abspath(path)
    raise WorkflowError(
        "the conformer and optimization stages run the mqc executable, and none was found. "
        "Set Protocol(executable=...) or MQC_EXECUTABLE, or skip those stages with "
        "Protocol(conformers=None, optimize=None)."
    )


#: Variables an MPI launcher sets that would make a child `mqc` try to join the
#: parent's job rather than start as the single rank it must be.
_MPI_PREFIXES = ("OMPI_", "PMIX_", "PMI_", "HYDRA_", "I_MPI_", "MPICH_")


def child_env(protocol):
    env = {k: v for k, v in os.environ.items() if not k.startswith(_MPI_PREFIXES)}
    if protocol.subprocess_env:
        env.update(protocol.subprocess_env)
    return env


def write_deck(stage, driver, geometry, species, directory):
    """Write ``start.xyz`` and a deck for ``driver`` beside it; return the deck name.

    The deck is what `mqc.MBE` would send, less the fragmentation block (a
    deck for an optimization has none, and ``level=0`` is a spelling of "not
    fragmented" that only the in-process path needs) plus the molecule.
    """
    from . import MBE

    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, "start.xyz"), "w") as handle:
        handle.write(geometry.to_xyz(species.name))
    document = MBE(None, level=0, driver=driver, **stage).settings()
    document["keywords"].pop("fragmentation", None)
    document["molecules"] = [
        {
            "xyz": "start.xyz",
            "molecular_charge": species.charge,
            "molecular_multiplicity": species.multiplicity,
        }
    ]
    name = f"{driver}.json"
    with open(os.path.join(directory, name), "w") as handle:
        json.dump(document, handle, indent=2)
    return name


def run_mqc(executable, deck, directory, env):
    """Run the executable on a deck; raise with its output if it fails."""
    done = subprocess.run([executable, deck], cwd=directory, env=env, capture_output=True, text=True)
    if done.returncode != 0:
        tail = "\n".join((done.stdout + "\n" + done.stderr).strip().splitlines()[-15:])
        raise WorkflowError(f"mqc exited {done.returncode} on {directory}/{deck}:\n{tail}")


# TODO(mqc): the conformer ensemble reaches Python only as the files CREST leaves
# in the working directory (`crest_conformers.xyz`, energies on the comment line),
# and only from a deck run through the executable: `mqc_session` and the C API
# refuse the `conformers` driver, and `run_conformer_search` returns no ensemble
# to its caller. Same for `optimize`, whose result is `output_<deck>_optimized.xyz`.
# Both belong behind the C interface, returning energies and geometries.
def sample_conformers(species, protocol, directory, executable, env, run=None):
    """CREST ensemble for one species: ``[(energy_hartree, Geometry), ...]``, lowest first."""
    run = run or run_mqc
    name = write_deck(protocol.conformers, "conformers", species.geometry, species, directory)
    run(executable, name, directory, env)
    path = os.path.join(directory, "crest_conformers.xyz")
    if not os.path.exists(path):
        raise WorkflowError(f"{species.name}: the conformer search left no {path}")
    with open(path) as handle:
        ensemble = _thermo.parse_xyz_ensemble(handle.read())
    if not ensemble or any(e is None for e, _ in ensemble):
        raise WorkflowError(f"{path}: could not read an energy from every structure's comment line")
    ensemble.sort(key=lambda item: item[0])
    return ensemble


def optimize_structure(species, geometry, protocol, directory, executable, env, run=None):
    """Optimize one geometry: ``(Geometry, converged, energy_hartree)``."""
    run = run or run_mqc
    name = write_deck(protocol.optimize, "optimize", geometry, species, directory)
    run(executable, name, directory, env)
    path = os.path.join(directory, f"output_{name[:-5]}_optimized.xyz")
    if not os.path.exists(path):
        # `optimized_xyz_path` derives the name from the output document's, and
        # a rename on the Fortran side should not read as a failed optimization.
        found = sorted(glob.glob(os.path.join(directory, "*_optimized.xyz")))
        if not found:
            raise WorkflowError(f"{species.name}: the optimization wrote no *_optimized.xyz in {directory}")
        path = found[0]
    with open(path) as handle:
        return _thermo.parse_optimized_xyz(handle.read())


# ---------------------------------------------------------------------------
#  The stages that run in the session
# ---------------------------------------------------------------------------


def frequency_stage(species, geometry, protocol, label):
    """Hessian, then the correction. Returns ``(correction, electronic_energy_hartree)``."""
    from . import MBE

    system = species.to_system(geometry)
    result = MBE(system, level=0, driver="hessian", **protocol.frequencies).run(label)
    freqs, thermo = result.frequencies, result.thermochemistry
    if freqs is None or thermo is None:
        raise WorkflowError(f"{label}: the Hessian run wrote no vibrational_analysis/thermochemistry block")
    correction = _thermo.free_energy_correction(
        freqs,
        species.n_atoms,
        thermo,
        nu0=protocol.nu0,
        alpha=protocol.alpha,
        b_av=protocol.b_av,
        standard_state=protocol.standard_state,
        imaginary=protocol.imaginary,
        qrrho=protocol.qrrho,
        multiplicity=species.multiplicity,
    )
    return correction, float(thermo["total_energies_hartree"]["electronic"])


def single_point(species, geometry, protocol, label):
    """Electronic energy in Hartree at the single-point level."""
    from . import MBE

    system = species.to_system(geometry)
    return MBE(system, level=0, driver="energy", **protocol.single_point).run(label).energy


# ---------------------------------------------------------------------------
#  One species
# ---------------------------------------------------------------------------


@dataclasses.dataclass
class SpeciesEvaluation:
    """What the stages produced for one species.

    ``e_hartree`` is the single-point electronic energy, ``zpe_kcal`` the
    zero-point energy, ``h_kcal`` the enthalpy ``E + H_corr`` and ``g_kcal`` the
    free energy ``E + G_corr`` in kcal/mol. With more than one conformer ``g_kcal``
    is their Boltzmann combination and the other three are averaged with the
    same weights, which neglects the (small) enthalpy of mixing.
    ``e_vertical_hartree`` is the single-point energy at the geometry the species
    was *given*, before any search or optimization -- None unless asked for.
    ``detail`` is JSON-safe; ``geometries`` is not, and holds the geometry of each
    conformer in the order of ``detail["conformers"]``.
    """

    name: str
    e_hartree: float
    zpe_kcal: float
    h_kcal: float
    g_kcal: float
    temperature_K: float
    detail: dict
    geometries: list
    e_vertical_hartree: float = None

    def lowest_geometry(self):
        """The geometry of the conformer with the largest weight, or None."""
        weights = [c["weight"] for c in self.detail["conformers"]]
        return self.geometries[weights.index(max(weights))]


def _evaluate_atom(species, protocol, prefix, say, run_single_point):
    geom = species.geometry
    if geom.n_atoms != 1:
        raise WorkflowError(f"{species.name}: internal error, not an atom")
    from . import _check_label

    temperature, pressure = protocol.conditions()
    label = f"{prefix}_{slug(species.name)}_atom"
    _check_label(label)
    say("single atom: no conformers, optimization or Hessian; analytic thermochemistry")
    e_sp = run_single_point(species, geom, protocol, label)
    corr = _thermo.atom_correction(
        geom.symbols[0], species.multiplicity, temperature, pressure, protocol.standard_state
    )
    e_kcal = e_sp * _thermo.HARTREE_TO_KCAL
    record = {
        "conformer": 0,
        "e_conformer_hartree": None,
        "e_optimized_hartree": None,
        "e_single_point_hartree": e_sp,
        "single_point_reused_frequency_energy": False,
        "zpe_kcal_mol": 0.0,
        "h_corr_kcal_mol": corr["h_corr"],
        "h_kcal_mol": e_kcal + corr["h_corr"],
        "g_corr_kcal_mol": corr["g_corr"],
        "g_kcal_mol": e_kcal + corr["g_corr"],
        "n_imaginary": 0,
        "imaginary_cm1": [],
        "tr_max_cm1": 0.0,
        "optimization_converged": True,
        "temperature_K": corr["temperature_K"],
        "weight": 1.0,
    }
    detail = {
        "n_conformers": 1,
        "n_imaginary": 0,
        "temperature_K": corr["temperature_K"],
        "atom": True,
        "conformers": [record],
    }
    say(f"G = {record['g_kcal_mol']:.3f} kcal/mol (atom)")
    return SpeciesEvaluation(
        species.name, e_sp, 0.0, record["h_kcal_mol"], record["g_kcal_mol"], corr["temperature_K"], detail, [geom], e_sp
    )


def evaluate_species(
    species,
    protocol,
    prefix,
    verbose=True,
    *,
    use_conformers=True,
    vertical=False,
    executable=None,
    run=None,
    locate_executable=None,
):
    """One species to one free energy, with per-conformer detail.

    Search, optimize, Hessian, single point, then Boltzmann-combine. Every
    conformer is kept in ``detail["conformers"]`` with its energies, correction,
    imaginary frequencies and final weight; ``detail["n_imaginary"]`` is the
    total across conformers. Imaginary frequencies are counted here and
    surfaced by the caller's warnings, not dropped.

    ``use_conformers=False`` skips the search for this species whatever the
    protocol says. ``vertical=True`` adds one single point at the geometry the
    species was given, for the relaxation energy. ``run`` and
    ``locate_executable`` replace `run_mqc` and `find_executable`; a caller
    passes them to let its own module's names be patched.
    """
    run = run or run_mqc
    locate_executable = locate_executable or find_executable
    from . import _check_label

    _check_label(prefix)
    name = slug(species.name)
    root = os.path.join(protocol.workdir, name)
    say = (lambda text: print(f"  [{species.name}] {text}", flush=True)) if verbose else (lambda t: None)

    if species._geometry is not None and species._geometry.n_atoms == 1:
        return _evaluate_atom(species, protocol, prefix, say, single_point)

    search = use_conformers and protocol.conformers is not None
    need_exe = search or protocol.optimize is not None
    executable = executable or (locate_executable(protocol) if need_exe else None)
    env = child_env(protocol)

    start = species._geometry
    e_vertical = None
    if vertical and start is not None:
        label = f"{prefix}_{name}_vert"
        _check_label(label)
        say("single point at the starting geometry")
        e_vertical = single_point(species, start, protocol, label)

    # 1. Structures to carry forward.
    if search:
        say("conformer search")
        ensemble = sample_conformers(species, protocol, os.path.join(root, "conformers"), executable, env, run)
        lowest = ensemble[0][0]
        window = protocol.energy_window_kcal / _thermo.HARTREE_TO_KCAL
        kept = [item for item in ensemble if item[0] - lowest <= window][: protocol.max_conformers]
        say(f"{len(ensemble)} conformers, {len(kept)} within {protocol.energy_window_kcal:g} kcal/mol")
    else:
        kept = [(None, start)]

    # 2. Optimize them.
    structures = []
    for k, (e_conf, geom) in enumerate(kept):
        converged, e_opt = True, None
        if protocol.optimize is not None:
            say(f"optimizing conformer {k}")
            geom, converged, e_opt = optimize_structure(
                species, geom, protocol, os.path.join(root, f"opt_c{k}"), executable, env, run
            )
            if not converged and not protocol.allow_unconverged:
                raise WorkflowError(
                    f"{species.name}: the optimization of conformer {k} did not converge. "
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
        label = f"{prefix}_{name}_c{s['k']}"
        _check_label(label + "_freq")
        say(f"frequencies, conformer {s['k']}")
        correction, e_freq = frequency_stage(species, s["geometry"], protocol, label + "_freq")
        if protocol.single_point_is_frequencies():
            e_sp, reused = e_freq, True
        else:
            say(f"single point, conformer {s['k']}")
            e_sp, reused = single_point(species, s["geometry"], protocol, label + "_sp"), False
        e_kcal = e_sp * _thermo.HARTREE_TO_KCAL
        records.append(
            {
                "conformer": s["k"],
                "e_conformer_hartree": s["e_conf"],
                "e_optimized_hartree": s["e_opt"],
                "e_single_point_hartree": e_sp,
                "single_point_reused_frequency_energy": reused,
                "zpe_kcal_mol": correction["zpe"],
                "h_corr_kcal_mol": correction["h_corr"],
                "h_kcal_mol": e_kcal + correction["h_corr"],
                "g_corr_kcal_mol": correction["g_corr"],
                "g_kcal_mol": e_kcal + correction["g_corr"],
                "n_imaginary": correction["n_imaginary"],
                "imaginary_cm1": correction["imaginary_cm1"],
                "tr_max_cm1": correction["tr_max_cm1"],
                "optimization_converged": s["converged"],
                "temperature_K": correction["temperature_K"],
            }
        )

    temperature = records[0]["temperature_K"]
    weights = _thermo.boltzmann_weights([r["g_kcal_mol"] for r in records], temperature)
    for r, w in zip(records, weights):
        r["weight"] = w
    g_micro = _thermo.boltzmann_combine([r["g_kcal_mol"] for r in records], temperature)
    average = lambda key: sum(w * r[key] for w, r in zip(weights, records))  # noqa: E731
    detail = {
        "n_conformers": len(records),
        "n_imaginary": sum(r["n_imaginary"] for r in records),
        "temperature_K": temperature,
        "conformers": records,
    }
    say(f"G = {g_micro:.3f} kcal/mol from {len(records)} conformer(s)")
    return SpeciesEvaluation(
        species.name,
        average("e_single_point_hartree"),
        average("zpe_kcal_mol"),
        average("h_kcal_mol"),
        g_micro,
        temperature,
        detail,
        [s["geometry"] for s in unique],
        e_vertical,
    )
