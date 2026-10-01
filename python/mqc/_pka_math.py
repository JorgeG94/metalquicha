"""The arithmetic of a pKa / isoelectric-point calculation, with no library in it.

Everything in `mqc.pka` that is a formula lives here, and nothing here imports
`mqc` or loads `libmqc`. That split exists so the arithmetic can be tested --
against hand values, analytic limits and a model whose answer is known -- on a
machine with no Fortran compiler, and so that a number which looks wrong can be
checked without a session, MPI or a build.

Units, throughout: free energies and corrections in **kcal/mol**, entropies in
**cal/(mol K)** where they come from the Fortran thermochemistry block (which
writes them that way) and in units of R where they are computed here, frequencies
in cm^-1, temperature in K. Hartree appears only where an energy is read from or
handed back to a calculation.
"""

import json
import math

# -- constants --------------------------------------------------------------
#
# Values the Fortran side already has (`mqc_physical_constants.F90`) are
# repeated here digit for digit rather than re-derived, so the two thermochemistry
# paths differ by method and not by a constant. The SI ones are CODATA 2018 and
# are used only by the free-rotor entropy, which has no Fortran counterpart.

HARTREE_TO_KCAL = 627.5094740631
R_KCAL = 1.98720425864e-3  #: kcal/(mol K)
R_CAL = 1.98720425864  #: cal/(mol K)
KB_HARTREE = 3.1668115634556e-6  #: Hartree/K
CM1_TO_KELVIN = 1.4387773538277  #: hc/k in K per cm^-1
LN10 = math.log(10.0)

_H = 6.62607015e-34  #: J s
_KB = 1.380649e-23  #: J/K
_C = 299792458.0  #: m/s
_R_L_ATM = 0.082057366080960  #: L atm/(mol K)

#: Literature proton, 1 atm gas -> 1 M aqueous, for the uncalibrated default.
#: G_gas(H+) = H - TS = 5/2 RT - T S_trans, which is -6.28 kcal/mol at 298.15 K
#: (Sackur-Tetrode with the proton's mass); delta G_solv is Tissandier et al.
#: Both are quoted for a *1 atm* gas and a *1 M* solution, which is why the
#: standard-state term is added when they are combined.
G_GAS_PROTON = -6.28
DG_SOLV_PROTON = -265.9

#: Grimme's quasi-RRHO parameters (Chem. Eur. J. 2012, 18, 9955).
QRRHO_NU0 = 100.0  #: cm^-1, the switch between oscillator and rotor
QRRHO_ALPHA = 4.0
QRRHO_B_AV = 1.0e-44  #: kg m^2, the average moment of inertia of a free rotor

DEFAULT_CLASS = "default"


class PKaError(RuntimeError):
    """The workflow could not produce a number it would stand behind."""


class NoIsoelectricPoint(PKaError):
    """The net charge never changes sign on the pH interval asked for."""


# ---------------------------------------------------------------------------
#  Thermochemistry
# ---------------------------------------------------------------------------


def thermal_rt(temperature):
    """RT in kcal/mol."""
    return R_KCAL * temperature


def standard_state_correction(temperature=298.15, pressure_atm=1.0):
    """The 1 atm gas to 1 M solution standard-state change, kcal/mol.

    ``RT ln(R T / P)`` with R in L atm/(mol K) is the work of compressing the
    ideal-gas molar volume (24.46 L at 298.15 K and 1 atm) to one litre. It is
    added to a solute's free energy because the Fortran translational entropy
    is evaluated at ``pressure_atm`` while an aqueous pKa is for 1 mol/L.
    Computed from T so that a calculation at another temperature is not
    silently given the 298.15 K constant.
    """
    return thermal_rt(temperature) * math.log(_R_L_ATM * temperature / pressure_atm)


def _s_rrho_over_r(nu, temperature):
    """Harmonic-oscillator entropy of one mode, in units of R."""
    u = CM1_TO_KELVIN * nu / temperature
    if u > 700.0:
        return 0.0
    em1 = math.expm1(u)
    return u / em1 - math.log1p(-math.exp(-u))


def _s_free_rotor_over_r(nu, temperature, b_av=QRRHO_B_AV):
    """Free-rotor entropy of one mode, in units of R (Grimme's eq. 4 and 5).

    The mode is treated as a rotor whose moment of inertia is the one the
    oscillator would have, mu = h / (8 pi^2 nu), capped against the average
    moment of inertia ``b_av`` of a free rotor: mu' = mu b_av / (mu + b_av).
    """
    nu_hz = _C * 100.0 * nu
    mu = _H / (8.0 * math.pi**2 * nu_hz)
    mu_prime = mu * b_av / (mu + b_av)
    return 0.5 + math.log(math.sqrt(8.0 * math.pi**3 * mu_prime * _KB * temperature / _H**2))


# TODO(mqc): the quasi-RRHO entropy, the free-rotor term, and the 1 atm -> 1 M
# standard-state correction belong in `mqc_thermochemistry.f90`, next to the RRHO
# vibrational entropy they replace, so that the deck path and `thermochemistry`
# in the JSON report the same free energy as this module. They are here because
# a Python-side change needs no build; until they move, `thermal_corrections_hartree
# .to_gibbs` in the JSON is plain RRHO and does not agree with `g_corr` below.


def qrrho_weight(nu, nu0=QRRHO_NU0, alpha=QRRHO_ALPHA):
    """Grimme's damping function, w = 1 / (1 + (nu0/nu)^alpha).

    One at high frequency, where a mode is an oscillator, and zero as it
    softens into a rotor; one half at ``nu0``.
    """
    if nu <= 0.0:
        return 0.0
    return 1.0 / (1.0 + (nu0 / nu) ** alpha)


def qrrho_entropy_over_r(nu, temperature, nu0=QRRHO_NU0, alpha=QRRHO_ALPHA, b_av=QRRHO_B_AV):
    """Quasi-RRHO entropy of one mode, in units of R.

    ``w S_RRHO + (1 - w) S_FR``. Only the entropy is interpolated, as in the
    paper; the enthalpy stays harmonic.
    """
    if nu <= 0.0:
        raise ValueError("qrrho_entropy_over_r takes a real, positive frequency")
    w = qrrho_weight(nu, nu0, alpha)
    return w * _s_rrho_over_r(nu, temperature) + (1.0 - w) * _s_free_rotor_over_r(
        nu, temperature, b_av
    )


def rrho_entropy_over_r(nu, temperature):
    """Plain harmonic-oscillator entropy of one mode, in units of R."""
    return _s_rrho_over_r(nu, temperature)


def split_frequencies(frequencies, n_atoms, is_linear=False):
    """Separate a 3N list into vibrations and the modes that are not.

    The Fortran list has the translations and rotations *in it*, at or near
    zero, and its ``n_imaginary_frequencies`` counts every negative entry,
    including a -0.3 cm^-1 rotational residual. Neither is useful here, so the
    ``3N - 6`` (``3N - 5`` linear, ``3N - 3`` for one atom) modes of smallest
    magnitude are taken as the translations and rotations and the rest are the
    vibrations; of those, a negative one is imaginary.

    Returns ``(real, imaginary, tr_max)``: the positive vibrations, the
    imaginary ones as negative numbers, and the largest magnitude among the
    modes dropped as translation or rotation -- the number to look at if a
    soft real mode might have been swallowed by them.
    """
    freqs = [float(f) for f in frequencies]
    if len(freqs) != 3 * n_atoms:
        raise PKaError(
            f"{len(freqs)} frequencies for {n_atoms} atoms: the vibrational analysis "
            f"is expected to return 3N = {3 * n_atoms}, translations and rotations included"
        )
    n_tr = 3 if n_atoms == 1 else (5 if is_linear else 6)
    order = sorted(range(len(freqs)), key=lambda i: abs(freqs[i]))
    dropped = set(order[:n_tr])
    tr_max = max((abs(freqs[i]) for i in dropped), default=0.0)
    vib = [f for i, f in enumerate(freqs) if i not in dropped]
    real = [f for f in vib if f > 0.0]
    imaginary = [f for f in vib if f < 0.0]
    return real, imaginary, tr_max


def free_energy_correction(
    frequencies,
    n_atoms,
    thermo,
    nu0=QRRHO_NU0,
    alpha=QRRHO_ALPHA,
    b_av=QRRHO_B_AV,
    standard_state=True,
    imaginary="drop",
    qrrho=True,
):
    """G_corr for one structure from its frequencies and thermochemistry block.

    ``thermo`` is the dict ``Result.thermochemistry`` returns. From it come the
    conditions (``temperature_K``, ``pressure_atm``), whether the molecule is
    linear, and the translational, rotational and electronic contributions,
    which are used as they are: it is only the vibrational part that is
    recomputed, because that is where the rigid-rotor-harmonic-oscillator
    approximation is wrong for the soft modes of a flexible molecule.

    The correction is ``ZPE + E_vib + E_trans + E_rot + RT - T S`` with the
    entropy ``S_trans + S_rot + S_elec + S_vib`` and ``S_vib`` quasi-RRHO, plus
    the 1 atm -> 1 M standard-state term when ``standard_state``. It is added to
    an electronic energy to give a free energy; it contains none.

    ``imaginary`` is ``"drop"`` (an imaginary mode contributes nothing, which is
    what the Fortran thermochemistry does) or ``"flip"`` (taken at its absolute
    value, the usual repair for the small imaginary mode of a loosely
    converged structure). Either way they are counted and returned, never
    discarded without a trace. ``qrrho=False`` gives plain RRHO, for comparison.

    Returns a dict, all energies in kcal/mol: ``g_corr`` (the one to use) and
    the pieces ``zpe``, ``e_vib``, ``e_trans``, ``e_rot``, ``rt``,
    ``ts_total``, ``s_vib_cal`` (cal/mol/K), ``standard_state``, together with
    ``temperature_K``, ``n_real``, ``n_imaginary``, ``imaginary_cm1`` and
    ``tr_max_cm1``.
    """
    if imaginary not in ("drop", "flip"):
        raise ValueError("imaginary must be 'drop' or 'flip'")
    temperature = float(thermo["temperature_K"])
    pressure = float(thermo.get("pressure_atm", 1.0))
    contributions = thermo["contributions"]
    is_linear = bool(thermo.get("is_linear", False))

    real, imag, tr_max = split_frequencies(frequencies, n_atoms, is_linear)
    modes = list(real)
    if imaginary == "flip":
        modes += [abs(f) for f in imag]

    zpe = 0.0
    e_vib = 0.0
    s_vib_over_r = 0.0
    for nu in modes:
        theta = CM1_TO_KELVIN * nu
        zpe += 0.5 * R_KCAL * theta
        u = theta / temperature
        if u < 700.0:
            e_vib += R_KCAL * theta / math.expm1(u)
        if qrrho:
            s_vib_over_r += qrrho_entropy_over_r(nu, temperature, nu0, alpha, b_av)
        else:
            s_vib_over_r += rrho_entropy_over_r(nu, temperature)
    s_vib_cal = R_CAL * s_vib_over_r

    e_trans = contributions["translational"]["energy_hartree"] * HARTREE_TO_KCAL
    e_rot = contributions["rotational"]["energy_hartree"] * HARTREE_TO_KCAL
    s_other_cal = (
        contributions["translational"]["entropy_cal_mol_K"]
        + contributions["rotational"]["entropy_cal_mol_K"]
        + contributions["electronic"]["entropy_cal_mol_K"]
    )
    rt = thermal_rt(temperature)
    ts_total = temperature * (s_vib_cal + s_other_cal) / 1000.0
    ss = standard_state_correction(temperature, pressure) if standard_state else 0.0
    g_corr = zpe + e_vib + e_trans + e_rot + rt - ts_total + ss
    return {
        "g_corr": g_corr,
        "zpe": zpe,
        "e_vib": e_vib,
        "e_trans": e_trans,
        "e_rot": e_rot,
        "rt": rt,
        "ts_total": ts_total,
        "s_vib_cal": s_vib_cal,
        "standard_state": ss,
        "temperature_K": temperature,
        "pressure_atm": pressure,
        "n_real": len(real),
        "n_imaginary": len(imag),
        "imaginary_cm1": imag,
        "tr_max_cm1": tr_max,
    }


# ---------------------------------------------------------------------------
#  Combining conformers and microstates
# ---------------------------------------------------------------------------


def logsumexp(values):
    """log(sum(exp(v))), stable for large magnitudes."""
    values = list(values)
    top = max(values)
    if top == float("-inf"):
        return top
    return top + math.log(sum(math.exp(v - top) for v in values))


def boltzmann_weights(free_energies, temperature):
    """Normalised weights exp(-G/RT) / sum, for free energies in kcal/mol."""
    rt = thermal_rt(temperature)
    logs = [-g / rt for g in free_energies]
    norm = logsumexp(logs)
    return [math.exp(v - norm) for v in logs]


def boltzmann_combine(free_energies, temperature):
    """One free energy for a set of conformers, -RT ln sum exp(-G_i/RT)."""
    if not free_energies:
        raise ValueError("no conformers to combine")
    rt = thermal_rt(temperature)
    return -rt * logsumexp([-g / rt for g in free_energies])


# ---------------------------------------------------------------------------
#  Geometry
# ---------------------------------------------------------------------------


class Geometry:
    """Element symbols and Cartesian coordinates in Angstrom.

    The library's ``System`` cannot be read back -- it has an atom count and
    no coordinates -- and the conformer and optimization stages need to write
    the structure out and read a different one back, so the workflow keeps its
    own copy. Plain data, so a `Microstate` can be built and checked without a
    session.
    """

    def __init__(self, symbols, coords):
        self.symbols = [str(s).strip().capitalize() for s in symbols]
        self.coords = [[float(x) for x in xyz] for xyz in coords]
        if len(self.symbols) != len(self.coords) or any(len(c) != 3 for c in self.coords):
            raise ValueError("symbols and coords must be N symbols and N rows of three numbers")
        if not self.symbols:
            raise ValueError("a geometry needs at least one atom")

    @property
    def n_atoms(self):
        return len(self.symbols)

    @classmethod
    def from_xyz(cls, path):
        with open(path) as handle:
            geoms = parse_xyz_ensemble(handle.read())
        if len(geoms) != 1:
            raise ValueError(f"{path} holds {len(geoms)} structures, expected one")
        return geoms[0][1]

    def to_xyz(self, comment=""):
        lines = [str(self.n_atoms), comment]
        for sym, (x, y, z) in zip(self.symbols, self.coords):
            lines.append(f"{sym:<3s} {x:18.10f} {y:18.10f} {z:18.10f}")
        return "\n".join(lines) + "\n"

    def __repr__(self):
        return f"<Geometry {self.n_atoms} atoms>"


def parse_xyz_ensemble(text):
    """Read a (multi-)structure .xyz into ``[(energy_or_None, Geometry), ...]``.

    The energy is the first number on the comment line, which is where CREST
    puts the absolute energy (Hartree) and where mqc's own optimized-geometry
    file puts nothing parseable at the start -- see `parse_optimized_xyz`.
    """
    lines = text.splitlines()
    out = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        body = lines[i + 2 : i + 2 + n]
        if len(body) != n:
            raise ValueError("truncated .xyz structure")
        symbols, coords = [], []
        for line in body:
            parts = line.split()
            symbols.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
        energy = None
        for token in comment.split():
            try:
                energy = float(token)
                break
            except ValueError:
                continue
        out.append((energy, Geometry(symbols, coords)))
        i += 2 + n
    return out


def parse_optimized_xyz(text):
    """Read mqc's ``output_<label>_optimized.xyz``.

    Returns ``(Geometry, converged, energy_hartree)``. The comment line is
    ``metalquicha converged, E = <f20.12> Hartree`` or ``metalquicha NOT
    CONVERGED, E = ...`` (`write_optimized_xyz` in the geometry optimizer); the
    file is written whichever it was, so the word has to be read.
    """
    lines = text.splitlines()
    comment = lines[1] if len(lines) > 1 else ""
    structures = parse_xyz_ensemble(text)
    if len(structures) != 1:
        raise ValueError("an optimized-geometry file holds one structure")
    converged = "NOT CONVERGED" not in comment.upper() and "CONVERGED" in comment.upper()
    energy = None
    if "E =" in comment:
        try:
            energy = float(comment.split("E =", 1)[1].split()[0])
        except (ValueError, IndexError):
            energy = None
    return structures[0][1], converged, energy


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
