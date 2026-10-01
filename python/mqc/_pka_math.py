"""The arithmetic of a pKa / isoelectric-point calculation, with no library in it.

Everything in `mqc.pka` that is a formula lives here, and nothing here imports
`mqc` or loads `libmqc`. That split exists so the arithmetic can be tested --
against hand values, analytic limits and a model whose answer is known -- on a
machine with no Fortran compiler, and so that a number which looks wrong can be
checked without a session, MPI or a build.

The method-independent part -- quasi-RRHO, the standard state, Boltzmann
combination, `Geometry` and the .xyz readers -- is in `mqc._thermo`, shared with
`mqc.bde`, and is re-exported here under its old names. What is left is the
protonation model: the proton, calibration, and the result.

Units are `mqc._thermo`'s: free energies in kcal/mol, frequencies in cm^-1.
"""

import json
import math

try:
    from . import _thermo
except ImportError:  # loaded by path, as python/tests/test_pka.py does: no parent package
    import importlib.util as _util
    import os as _os
    import sys as _sys

    _spec = _util.spec_from_file_location(
        "mqc_thermo_under_test", _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "_thermo.py")
    )
    _thermo = _util.module_from_spec(_spec)
    _sys.modules[_spec.name] = _thermo
    _spec.loader.exec_module(_thermo)

# Re-exported: these are what `mqc.pka` and the tests have always found here.
HARTREE_TO_KCAL = _thermo.HARTREE_TO_KCAL
R_KCAL = _thermo.R_KCAL
R_CAL = _thermo.R_CAL
KB_HARTREE = _thermo.KB_HARTREE
CM1_TO_KELVIN = _thermo.CM1_TO_KELVIN
LN10 = _thermo.LN10
_H = _thermo._H
_KB = _thermo._KB
_C = _thermo._C
_R_L_ATM = _thermo._R_L_ATM
QRRHO_NU0 = _thermo.QRRHO_NU0
QRRHO_ALPHA = _thermo.QRRHO_ALPHA
QRRHO_B_AV = _thermo.QRRHO_B_AV
thermal_rt = _thermo.thermal_rt
standard_state_correction = _thermo.standard_state_correction
_s_rrho_over_r = _thermo._s_rrho_over_r
_s_free_rotor_over_r = _thermo._s_free_rotor_over_r
qrrho_weight = _thermo.qrrho_weight
qrrho_entropy_over_r = _thermo.qrrho_entropy_over_r
rrho_entropy_over_r = _thermo.rrho_entropy_over_r
split_frequencies = _thermo.split_frequencies
free_energy_correction = _thermo.free_energy_correction
logsumexp = _thermo.logsumexp
boltzmann_weights = _thermo.boltzmann_weights
boltzmann_combine = _thermo.boltzmann_combine
Geometry = _thermo.Geometry
parse_xyz_ensemble = _thermo.parse_xyz_ensemble
parse_optimized_xyz = _thermo.parse_optimized_xyz

#: Literature proton, 1 atm gas -> 1 M aqueous, for the uncalibrated default.
#: G_gas(H+) = H - TS = 5/2 RT - T S_trans, which is -6.28 kcal/mol at 298.15 K
#: (Sackur-Tetrode with the proton's mass); delta G_solv is Tissandier et al.
#: Both are quoted for a *1 atm* gas and a *1 M* solution, which is why the
#: standard-state term is added when they are combined.
G_GAS_PROTON = -6.28
DG_SOLV_PROTON = -265.9

DEFAULT_CLASS = "default"

#: The workflows' shared error. `PKaError` is another name for it, so that what
#: the shared arithmetic raises is caught by the name this module has always had.
PKaError = _thermo.WorkflowError


class NoIsoelectricPoint(PKaError):
    """The net charge never changes sign on the pH interval asked for."""


# ---------------------------------------------------------------------------
#  The protonation model
# ---------------------------------------------------------------------------


class MicrostateData:
    """What the thermodynamics needs to know about one microstate.

    ``free_energy`` is in kcal/mol and is already the Boltzmann combination of
    the conformers. ``site_classes`` is one label per ionizable proton the
    microstate carries, which is what lets a calibration shift attach to a
    *site* and so stay path independent -- see `PKaModel`.
    """

    def __init__(self, name, charge, n_protons, free_energy, site_classes=None, detail=None):
        self.name = str(name)
        self.charge = int(charge)
        self.n_protons = int(n_protons)
        self.free_energy = float(free_energy)
        self.site_classes = normalize_site_classes(site_classes, self.n_protons)
        self.detail = detail or {}


def normalize_site_classes(site_class, n_protons):
    """One class label per carried proton.

    ``None`` is the single class ``"default"``; a string applies to every
    proton the microstate carries; a sequence must have one entry each.
    """
    if site_class is None:
        return (DEFAULT_CLASS,) * n_protons
    if isinstance(site_class, str):
        return (site_class,) * n_protons
    classes = tuple(str(c) for c in site_class)
    if len(classes) != n_protons:
        raise ValueError(
            f"{len(classes)} site classes for {n_protons} protons: give one class per "
            "ionizable proton the microstate carries, or a single string for all of them"
        )
    return classes


def default_proton_free_energy(temperature=298.15):
    """The uncalibrated G(H+) for aqueous solution at 1 M, kcal/mol.

    ``G_gas(H+, 1 atm) + RT ln(24.46) + dG_solv(H+) = -6.28 + 1.89 - 265.9``,
    about -270.3 kcal/mol. It is consistent with solute free energies that
    carry the same standard-state term, and at GFN2 level it is expected to be
    several pKa units wrong -- the proton is a number the electronic-structure
    method does not get to choose, and the anion's solvation is the weak point.
    That is what `calibrate` is for.
    """
    return G_GAS_PROTON + standard_state_correction(temperature) + DG_SOLV_PROTON


class PKaModel:
    """Microstate free energies, the proton, and everything derived from them.

    A microstate ``s`` at proton chemical potential ``mu`` has the grand
    potential ::

        Omega_s = G_s - n_s * mu - sum_{c in sites(s)} delta_c

    where ``mu(pH) = G_H+ - RT ln(10) pH`` and ``delta_c = G_H+(c) - G_H+`` is
    the calibration shift of site class ``c``. Populations are
    ``exp(-Omega_s/RT)``, normalised. The shift is attached to the *protons a
    microstate carries* and not to a transition, so the free energy of every
    thermodynamic cycle is the same whichever way round it is walked: the
    micro pKa of ``s -> t`` is ``(G_t - G_s + G_H+(c)) / (RT ln 10)`` with ``c``
    the class of the proton that leaves, and the macro pKas and the charge curve
    are built from the same ``Omega`` and so cannot disagree with them.

    Slope is fixed at one (shift-only). A fitted slope would make the pKa of a
    transition depend on how many protons the microstate carries.

    Net charge has to change one-to-one with the proton count, because the
    populations are over protonation states of one molecule; that is checked.
    """

    def __init__(self, microstates, temperature=298.15, proton_free_energy=None, class_free_energy=None):
        self.microstates = list(microstates)
        if not self.microstates:
            raise ValueError("a model needs at least one microstate")
        names = [m.name for m in self.microstates]
        if len(set(names)) != len(names):
            raise ValueError("microstate names must be unique")
        offsets = {m.charge - m.n_protons for m in self.microstates}
        if len(offsets) != 1:
            raise ValueError(
                "charge - n_protons must be the same for every microstate: a microstate "
                "that differs from another by a proton must differ by one unit of charge"
            )
        self.temperature = float(temperature)
        self.g_h = (
            default_proton_free_energy(self.temperature)
            if proton_free_energy is None
            else float(proton_free_energy)
        )
        self.class_g_h = dict(class_free_energy or {})
        self._by_name = {m.name: m for m in self.microstates}

    # -- the proton --------------------------------------------------------

    @property
    def rt(self):
        return thermal_rt(self.temperature)

    def proton_free_energy(self, site_class=DEFAULT_CLASS):
        """G(H+) in kcal/mol for one site class; the model's own if it has no entry."""
        return self.class_g_h.get(site_class, self.g_h)

    def _shift(self, site_class):
        return self.proton_free_energy(site_class) - self.g_h

    def _epsilon(self, m):
        return sum(self._shift(c) for c in m.site_classes)

    def effective_free_energy(self, name):
        """``G_s - sum delta_c`` -- what the microstate is worth once the shifts are in."""
        m = self._by_name[name]
        return m.free_energy - self._epsilon(m)

    # -- pKas --------------------------------------------------------------

    def micro_pkas(self):
        """Microscopic pKas between every pair differing by exactly one proton.

        A list of dicts ``{acid, base, pKa, site_class}``; ``site_class`` is
        the class of the proton that leaves. When the acid carries several
        protons of different classes and the base differs from it by one of
        them, the class is the one the base is missing; an ambiguous pair (the
        acid has a class the base does not account for) is reported with the
        class that differs, found by multiset difference, and ``None`` when
        it is not a single class.
        """
        out = []
        for a in self.microstates:
            for b in self.microstates:
                if a.n_protons - b.n_protons != 1:
                    continue
                lost = _multiset_difference(a.site_classes, b.site_classes)
                cls = lost[0] if len(lost) == 1 else None
                eff_a, eff_b = self.effective_free_energy(a.name), self.effective_free_energy(b.name)
                pka = (eff_b - eff_a + self.g_h) / (self.rt * LN10)
                out.append({"acid": a.name, "base": b.name, "pKa": pka, "site_class": cls})
        return out

    def level_free_energy(self, n_protons):
        """``-RT ln sum_s exp(-G_s^eff / RT)`` over the microstates at one level."""
        members = [m for m in self.microstates if m.n_protons == n_protons]
        if not members:
            raise KeyError(n_protons)
        return -self.rt * logsumexp(
            [-self.effective_free_energy(m.name) / self.rt for m in members]
        )

    def macro_pkas(self):
        """Macroscopic pKas from partition sums over each protonation level.

        ``pKa(n -> n-1)`` is the pH at which the two levels are equally
        populated, ``(A_{n-1} + G_H+ - A_n) / (RT ln 10)`` with ``A`` the level
        free energy from `level_free_energy`. Returned highest level first, as
        dicts ``{from_n, to_n, pKa}``. Levels with no microstate are skipped
        rather than interpolated, so a gap in the supplied ``n_protons`` is a
        missing pKa and not a wrong one.
        """
        levels = sorted({m.n_protons for m in self.microstates}, reverse=True)
        out = []
        for hi, lo in zip(levels, levels[1:]):
            if hi - lo != 1:
                continue
            pka = (self.level_free_energy(lo) + self.g_h - self.level_free_energy(hi)) / (
                self.rt * LN10
            )
            out.append({"from_n": hi, "to_n": lo, "pKa": pka})
        return out

    # -- populations -------------------------------------------------------

    def mu(self, ph):
        """Proton chemical potential at a pH, kcal/mol."""
        return self.g_h - self.rt * LN10 * ph

    def populations(self, ph):
        """``{name: fraction}`` at a pH, from ``exp(-Omega_s/RT)``."""
        mu = self.mu(ph)
        logs = [
            -(m.free_energy - m.n_protons * mu - self._epsilon(m)) / self.rt
            for m in self.microstates
        ]
        norm = logsumexp(logs)
        return {m.name: math.exp(v - norm) for m, v in zip(self.microstates, logs)}

    def charge(self, ph):
        """Mean net charge at a pH."""
        pops = self.populations(ph)
        return sum(m.charge * pops[m.name] for m in self.microstates)

    def isoelectric_point(self, lo=0.0, hi=14.0, tol=1.0e-9):
        """The pH at which the mean net charge is zero, by bisection.

        Raises `NoIsoelectricPoint` when the charge does not change sign on
        ``[lo, hi]`` -- a molecule whose every microstate is cationic, or whose
        calibrated pKas put the crossing outside the interval -- rather than
        returning an end of the interval, which would read as an answer.
        """
        f_lo, f_hi = self.charge(lo), self.charge(hi)
        if not (f_lo > 0.0 > f_hi):
            raise NoIsoelectricPoint(
                f"net charge is {f_lo:+.3f} at pH {lo:g} and {f_hi:+.3f} at pH {hi:g}: "
                "no sign change, so there is no isoelectric point on this interval"
            )
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            if self.charge(mid) > 0.0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)


def _multiset_difference(a, b):
    remaining = list(b)
    lost = []
    for item in a:
        if item in remaining:
            remaining.remove(item)
        else:
            lost.append(item)
    return lost


# ---------------------------------------------------------------------------
#  Calibration
# ---------------------------------------------------------------------------


class Reference:
    """One calibration datum: an acid, its conjugate base, and a measured pKa.

    ``acid`` and ``base`` are microstate names. ``site_class`` is the class of
    the proton that leaves, and the one whose shift this datum informs.
    """

    def __init__(self, acid, base, pka, site_class=DEFAULT_CLASS):
        self.acid = acid
        self.base = base
        self.pka = float(pka)
        self.site_class = site_class


class Calibration:
    """The fitted proton free energies and what they left over."""

    def __init__(self, class_free_energy, shifts, residuals, counts):
        self.class_free_energy = class_free_energy  #: class -> G(H+), kcal/mol
        self.shifts = shifts  #: class -> shift from the uncalibrated G(H+), kcal/mol
        self.residuals = residuals  #: (acid, base) -> exp - calc after the fit, pKa units
        self.counts = counts  #: class -> number of references

    def to_dict(self):
        return {
            "class_free_energy_kcal_mol": dict(self.class_free_energy),
            "shifts_kcal_mol": dict(self.shifts),
            "residuals_pKa": {f"{a}->{b}": r for (a, b), r in self.residuals.items()},
            "n_references": dict(self.counts),
        }


def calibrate(free_energies, references, temperature=298.15, proton_free_energy=None):
    """Fit one proton free energy per site class from reference pKas.

    ``free_energies`` maps microstate name to its free energy in kcal/mol (a
    `PKaResult` has them). For each reference the uncalibrated pKa is
    ``(G_base + G_H+ - G_acid) / (RT ln 10)``; the shift for a class is the
    **mean residual** of its references, converted to kcal/mol, and the class's
    G(H+) is the default plus that shift. Shift-only -- the slope stays one.

    One reference per class makes this a *relative* pKa: the shift is whatever
    makes that one compound right, and every other pKa on the same class is then
    its offset from it, which is the quantity the method is better at than the
    absolute.

    Returns a `Calibration`; hand ``calibration.class_free_energy`` to
    `PKaModel` (``class_free_energy=``) or use `PKaResult.with_calibration`.
    """
    g_h = (
        default_proton_free_energy(temperature) if proton_free_energy is None else proton_free_energy
    )
    rt_ln10 = thermal_rt(temperature) * LN10
    by_class = {}
    for ref in references:
        try:
            g_acid, g_base = free_energies[ref.acid], free_energies[ref.base]
        except KeyError as exc:
            raise PKaError(f"reference names a microstate that was not computed: {exc}") from None
        calc = (g_base + g_h - g_acid) / rt_ln10
        by_class.setdefault(ref.site_class, []).append((ref, ref.pka - calc))
    if not by_class:
        raise PKaError("calibration needs at least one reference")
    shifts, class_g, counts, residuals = {}, {}, {}, {}
    for cls, items in by_class.items():
        mean = sum(r for _, r in items) / len(items)
        shifts[cls] = mean * rt_ln10
        class_g[cls] = g_h + shifts[cls]
        counts[cls] = len(items)
        for ref, r in items:
            residuals[(ref.acid, ref.base)] = r - mean
    return Calibration(class_g, shifts, residuals, counts)


# ---------------------------------------------------------------------------
#  The result
# ---------------------------------------------------------------------------


class PKaResult:
    """Microstate free energies and what follows from them.

    ``microstates`` is a list of `MicrostateData` whose ``detail`` carries the
    per-conformer record the workflow kept (``conformers``, ``n_conformers``,
    ``n_imaginary``). The thermodynamics is a `PKaModel`; this adds the pieces a
    reader of a finished calculation wants -- the numbers, the warnings, JSON --
    and `with_calibration`, which returns a new result rather than editing this
    one, so the uncalibrated numbers stay available beside the calibrated.
    """

    #: Two microstates of the same protonation level closer than this are
    #: flagged: GFN2 does not order tautomers and zwitterions reliably, so a
    #: difference this small is not a ranking.
    CLOSE_KCAL = 3.0

    def __init__(self, microstates, temperature=298.15, class_free_energy=None, calibration=None, protocol=None):
        self.microstates = list(microstates)
        self.temperature = float(temperature)
        self.calibration = calibration
        self.protocol = protocol
        self.model = PKaModel(self.microstates, self.temperature, class_free_energy=class_free_energy)

    def with_calibration(self, references):
        """A new result whose proton free energies are fitted to ``references``.

        `Reference` objects, or ``(acid, base, pKa)`` / ``(acid, base, pKa,
        site_class)`` tuples.
        """
        references = [r if isinstance(r, Reference) else Reference(*r) for r in references]
        cal = calibrate(self.free_energies, references, self.temperature)
        return PKaResult(
            self.microstates, self.temperature, cal.class_free_energy, cal, self.protocol
        )

    @property
    def free_energies(self):
        """``{name: G}`` in kcal/mol, the Boltzmann combination over conformers."""
        return {m.name: m.free_energy for m in self.microstates}

    @property
    def micro_pkas(self):
        return self.model.micro_pkas()

    @property
    def macro_pkas(self):
        return self.model.macro_pkas()

    @property
    def pI(self):
        """The isoelectric point, or None; see `pI_note` for why it is None."""
        try:
            return self.model.isoelectric_point()
        except NoIsoelectricPoint:
            return None

    @property
    def pI_note(self):
        """Empty when there is an isoelectric point, the reason when there is not."""
        try:
            self.model.isoelectric_point()
        except NoIsoelectricPoint as exc:
            return str(exc)
        return ""

    def populations(self, ph):
        return self.model.populations(ph)

    def charge(self, ph):
        return self.model.charge(ph)

    @property
    def warnings(self):
        """Things a reader should look at before trusting a number, as strings."""
        out = []
        for m in self.microstates:
            n_imag = m.detail.get("n_imaginary", 0)
            if n_imag:
                out.append(
                    f"{m.name}: {n_imag} imaginary frequenc{'y' if n_imag == 1 else 'ies'} "
                    "among its conformers"
                )
        by_level = {}
        for m in self.microstates:
            by_level.setdefault(m.n_protons, []).append(m)
        for n, group in sorted(by_level.items()):
            group = sorted(group, key=lambda m: m.free_energy)
            for first, second in zip(group, group[1:]):
                gap = second.free_energy - first.free_energy
                if gap < self.CLOSE_KCAL:
                    out.append(
                        f"{first.name} and {second.name} (both {n} protons) are within "
                        f"{gap:.2f} kcal/mol: the ordering is not reliable at this level"
                    )
        if self.calibration is None:
            out.append(
                "uncalibrated: the proton free energy is a literature value and absolute "
                "pKas are expected to be several units off at xTB level"
            )
        return out

    def to_dict(self):
        pi = self.pI
        return {
            "temperature_K": self.temperature,
            "proton_free_energy_kcal_mol": self.model.g_h,
            "class_free_energy_kcal_mol": dict(self.model.class_g_h),
            "calibration": self.calibration.to_dict() if self.calibration else None,
            "microstates": [
                {
                    "name": m.name,
                    "charge": m.charge,
                    "n_protons": m.n_protons,
                    "site_classes": list(m.site_classes),
                    "free_energy_kcal_mol": m.free_energy,
                    **m.detail,
                }
                for m in self.microstates
            ],
            "micro_pkas": self.micro_pkas,
            "macro_pkas": self.macro_pkas,
            "isoelectric_point": pi,
            "isoelectric_point_note": self.pI_note,
            "warnings": self.warnings,
            "protocol": self.protocol,
        }

    def to_json(self, path=None, indent=2):
        """The result as JSON text; written to ``path`` too when one is given."""
        text = json.dumps(self.to_dict(), indent=indent, default=_json_default)
        if path is not None:
            with open(path, "w") as handle:
                handle.write(text + "\n")
        return text

    def __repr__(self):
        return f"<PKaResult {len(self.microstates)} microstates, pI={self.pI}>"


def _json_default(obj):
    if hasattr(obj, "to_dict"):
        return obj.to_dict()
    raise TypeError(f"{type(obj).__name__} is not JSON serializable")
