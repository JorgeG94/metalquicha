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

The staged evaluation of a species is `mqc._workflow`, shared with `mqc.bde`.
The formulas -- quasi-RRHO, the standard state, Boltzmann combination
(`mqc._thermo`, shared) and the protonation model, calibration and the result
(`mqc._pka_math`) -- import nothing from the library and are what the tests
exercise.
"""

import copy
import dataclasses

from . import _pka_math, _workflow
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

# The names the stages had in this module before they were shared. `run` and
# `evaluate_microstate` look `_run_mqc` and `_find_executable` up here at call
# time, so a caller (or a test) that replaces them in this module is honoured.
_slug = _workflow.slug
_run_mqc = _workflow.run_mqc
_find_executable = _workflow.find_executable
_child_env = _workflow.child_env
_deck = _workflow.write_deck
sample_conformers = _workflow.sample_conformers
optimize_structure = _workflow.optimize_structure

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


def _default_stage():
    return copy.deepcopy(_XTB_WATER)


@dataclasses.dataclass
class Protocol(_workflow.StageProtocol):
    """What to run at each stage, and how to turn it into a free energy.

    The four stage dicts are keyword arguments of `mqc.MBE`; the workflow
    supplies the system, ``level=0`` and the driver. ``conformers`` and
    ``optimize`` may be None to skip the stage and use the structure as given.
    **The default is GFN2 with ALPB water at every stage, and the 1 atm -> 1 M
    standard-state term on** (`mqc.bde.Protocol` is the gas-phase counterpart).

    ``conformers`` is the *refinement* level CREST re-ranks its ensemble with;
    CREST samples with its own GFN2 and takes the solvent from ``xtb`` (only
    ALPB and GBSA with a named solvent can be passed to it).

    ``energy_window_kcal`` and ``max_conformers`` choose which of the ensemble
    is carried forward, by the refinement energy and before the free energies
    are known. A window narrower than the frequency-dependent part of the free
    energy can move is a choice to drop conformers that would have mattered.

    ``qrrho``, ``nu0``, ``alpha`` and ``b_av`` are Grimme's quasi-RRHO
    settings; ``standard_state`` adds the 1 atm to 1 M term;
    ``imaginary`` is ``"drop"`` or ``"flip"`` (see `free_energy_correction`).
    """

    conformers: dict = dataclasses.field(default_factory=_default_stage)
    optimize: dict = dataclasses.field(default_factory=_default_stage)
    frequencies: dict = dataclasses.field(default_factory=_default_stage)
    single_point: dict = dataclasses.field(default_factory=_default_stage)
    workdir: str = "pka_work"


class Microstate(_workflow.Species):
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

    kind = "microstate"

    def __init__(self, name, system, charge, n_protons, site_class=None, multiplicity=1):
        super().__init__(name, system, charge, multiplicity)
        self.n_protons = int(n_protons)
        self.site_classes = _pka_math.normalize_site_classes(site_class, self.n_protons)

    def __repr__(self):
        return f"<Microstate {self.name} q={self.charge:+d} n_H={self.n_protons}>"


def find_executable(protocol=None):
    """The ``mqc`` executable the conformer and optimization stages will run.

    Raises `PKaError` naming how to point at one when there is none -- which is
    also the way to ask in advance whether those stages can run at all.
    """
    return _find_executable(protocol or Protocol())


def evaluate_microstate(microstate, protocol=None, prefix="pka", verbose=True, _executable=None):
    """One microstate to one free energy: `MicrostateData` with per-conformer detail.

    Search, optimize, Hessian, single point, then Boltzmann-combine. Every
    conformer is kept in ``detail["conformers"]`` with its energies, correction,
    imaginary frequencies and final weight; ``detail["n_imaginary"]`` is the
    total across conformers. Imaginary frequencies are counted here and
    surfaced by `PKaResult.warnings`, not dropped. The stages themselves are in
    `mqc._workflow`.
    """
    protocol = protocol or Protocol()
    evaluation = _workflow.evaluate_species(
        microstate,
        protocol,
        prefix,
        verbose,
        executable=_executable,
        run=lambda *args: _run_mqc(*args),
        locate_executable=lambda p: _find_executable(p),
    )
    return MicrostateData(
        microstate.name,
        microstate.charge,
        microstate.n_protons,
        evaluation.g_kcal,
        microstate.site_classes,
        evaluation.detail,
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
