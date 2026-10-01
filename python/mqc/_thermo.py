"""The thermochemistry every workflow shares, with no library in it.

Quasi-RRHO, the standard state, the conformer Boltzmann sum, the one-atom
thermochemistry, and the `Geometry` and .xyz parsing the staged workflows
(`mqc.pka`, `mqc.bde`) pass around. Nothing here imports `mqc` or loads
`libmqc`, so all of it is tested -- against hand values and analytic limits --
on a machine with no Fortran compiler.

Units, throughout: free energies and corrections in **kcal/mol**, entropies in
**cal/(mol K)** where they come from the Fortran thermochemistry block (which
writes them that way) and in units of R where they are computed here, frequencies
in cm^-1, temperature in K. Hartree appears only where an energy is read from or
handed back to a calculation.
"""

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

#: Grimme's quasi-RRHO parameters (Chem. Eur. J. 2012, 18, 9955).
QRRHO_NU0 = 100.0  #: cm^-1, the switch between oscillator and rotor
QRRHO_ALPHA = 4.0
QRRHO_B_AV = 1.0e-44  #: kg m^2, the average moment of inertia of a free rotor

KCAL_TO_KJ = 4.184  #: thermochemical calorie
_AMU_KG = 1.66053906660e-27  #: `AMU_TO_KG`, `mqc_physical_constants.F90`
_ATM_PA = 101325.0  #: `ATM_TO_PA`



class WorkflowError(RuntimeError):
    """A staged workflow could not produce a number it would stand behind."""


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
        raise WorkflowError(
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
    multiplicity=None,
):
    """G_corr and H_corr for one structure from its frequencies and thermochemistry block.

    ``thermo`` is the dict ``Result.thermochemistry`` returns. From it come the
    conditions (``temperature_K``, ``pressure_atm``), whether the molecule is
    linear, and the translational, rotational and electronic contributions,
    which are used as they are: it is only the vibrational part that is
    recomputed, because that is where the rigid-rotor-harmonic-oscillator
    approximation is wrong for the soft modes of a flexible molecule.

    The enthalpy correction is ``H_corr = ZPE + E_vib + E_trans + E_rot + RT``
    and the free-energy correction is ``G_corr = H_corr - T S`` with the
    entropy ``S_trans + S_rot + S_elec + S_vib`` and ``S_vib`` quasi-RRHO, plus
    the 1 atm -> 1 M standard-state term when ``standard_state``. Both are
    added to an electronic energy; they contain none. The enthalpy is the
    harmonic one -- quasi-RRHO interpolates the entropy only.

    **The electronic entropy is computed here, as ``R ln(multiplicity)``, and
    the Fortran block's own value is not used.** The Fortran thermochemistry
    has the term but no caller passes it the multiplicity (every
    ``compute_thermochemistry`` call leaves ``spin_multiplicity`` at its
    default of 1), so for a radical its ``S_elec`` is 0 and its
    ``spin_multiplicity`` reads 1. ``multiplicity=None`` takes whatever the
    block reports, which is the old behaviour and is right for a closed shell.

    ``imaginary`` is ``"drop"`` (an imaginary mode contributes nothing, which is
    what the Fortran thermochemistry does) or ``"flip"`` (taken at its absolute
    value, the usual repair for the small imaginary mode of a loosely
    converged structure). Either way they are counted and returned, never
    discarded without a trace. ``qrrho=False`` gives plain RRHO, for comparison.

    Returns a dict, all energies in kcal/mol: ``g_corr`` (the one to use) and
    the pieces ``zpe``, ``e_vib``, ``e_trans``, ``e_rot``, ``rt``,
    ``ts_total``, ``s_vib_cal`` (cal/mol/K), ``standard_state``, together with
    ``temperature_K``, ``n_real``, ``n_imaginary``, ``imaginary_cm1`` and
    ``tr_max_cm1``; and ``h_corr`` and ``s_elec_cal``.
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
    if multiplicity is None:
        s_elec_cal = contributions["electronic"]["entropy_cal_mol_K"]
        multiplicity = int(thermo.get("spin_multiplicity", 1))
    else:
        s_elec_cal = electronic_entropy_cal(multiplicity)
    s_other_cal = (
        contributions["translational"]["entropy_cal_mol_K"]
        + contributions["rotational"]["entropy_cal_mol_K"]
        + s_elec_cal
    )
    rt = thermal_rt(temperature)
    ts_total = temperature * (s_vib_cal + s_other_cal) / 1000.0
    ss = standard_state_correction(temperature, pressure) if standard_state else 0.0
    h_corr = zpe + e_vib + e_trans + e_rot + rt
    g_corr = h_corr - ts_total + ss
    return {
        "g_corr": g_corr,
        "h_corr": h_corr,
        "s_elec_cal": s_elec_cal,
        "multiplicity": int(multiplicity),
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
#  Elements, and one atom
# ---------------------------------------------------------------------------

#: H through Og, the order of `mqc_elements.f90`.
ELEMENTS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
    "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I "
    "Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt "
    "Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr "
    "Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()

#: Standard atomic masses in amu, copied from `element_masses` in
#: `src/core/mqc_elements.f90` so that an atom's translational entropy is
#: evaluated with the mass the molecules' is.
_MASSES = (
    1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999,
    18.998, 20.18, 22.99, 24.305, 26.982, 28.085, 30.974, 32.06,
    35.45, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996,
    54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.63,
    74.922, 78.971, 79.904, 83.798, 85.468, 87.62, 88.906, 91.224,
    92.906, 95.95, 98.0, 101.07, 102.91, 106.42, 107.87, 112.41,
    114.82, 118.71, 121.76, 127.6, 126.9, 131.29, 132.91, 137.33,
    138.91, 140.12, 140.91, 144.24, 145.0, 150.36, 151.96, 157.25,
    158.93, 162.5, 164.93, 167.26, 168.93, 173.05, 174.97, 178.49,
    180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59,
    204.38, 207.2, 208.98, 209.0, 210.0, 222.0, 223.0, 226.0,
    227.0, 232.04, 231.04, 238.03, 237.0, 244.0, 243.0, 247.0,
    247.0, 251.0, 252.0, 257.0, 258.0, 259.0, 262.0, 267.0,
    268.0, 271.0, 272.0, 270.0, 276.0, 281.0, 280.0, 285.0,
    284.0, 289.0, 288.0, 293.0, 294.0, 294.0,
)


def atomic_number(symbol):
    """Z for an element symbol, case-insensitive; ValueError for anything else."""
    try:
        return ELEMENTS.index(str(symbol).strip().capitalize()) + 1
    except ValueError:
        raise ValueError(f"unknown element symbol {symbol!r}") from None


def atomic_mass(symbol):
    """Standard atomic mass in amu."""
    return _MASSES[atomic_number(symbol) - 1]


def electronic_entropy_cal(multiplicity):
    """``R ln(2S+1)`` in cal/(mol K): the entropy of the spin degeneracy.

    A single electronic state at the given multiplicity, nothing thermally
    populated above it. For an atom that ignores the spin-orbit multiplet
    (a halogen's ``2P`` is really ``2P3/2`` with a degeneracy of four, and the
    ``2P1/2`` level a few kcal/mol above it), which is the same level of
    approximation as the 2S+1 itself.
    """
    if int(multiplicity) < 1:
        raise ValueError("multiplicity must be a positive integer")
    return R_CAL * math.log(int(multiplicity))


def translational_entropy_cal(mass_amu, temperature=298.15, pressure_atm=1.0):
    """Sackur-Tetrode translational entropy, cal/(mol K), at ``pressure_atm``.

    ``S = R [ ln( (2 pi m kT / h^2)^(3/2) kT / P ) + 5/2 ]``. For hydrogen at
    298.15 K and 1 atm this is 26.015 cal/(mol K) (NIST's 114.7 J/(mol K) for
    H at 1 bar is that plus ``R ln 2`` and the 1 atm -> 1 bar term).
    """
    m = mass_amu * _AMU_KG
    q_over_n = (2.0 * math.pi * m * _KB * temperature / _H**2) ** 1.5 * _KB * temperature / (
        pressure_atm * _ATM_PA
    )
    return R_CAL * (math.log(q_over_n) + 2.5)


def atom_correction(symbol, multiplicity, temperature=298.15, pressure_atm=1.0, standard_state=True):
    """The thermal correction for a single atom, in the keys of `free_energy_correction`.

    An atom has no vibration and no rotation, so there is no Hessian to ask
    for: ``E = 3/2 RT``, ``H = E + RT = 5/2 RT`` and ``S`` is the Sackur-Tetrode
    translational entropy plus ``R ln(multiplicity)``. Spin-orbit coupling is
    ignored. Everything else in the dict is zero or empty, so a one-atom
    species goes through the same arithmetic as a molecule.
    """
    rt = thermal_rt(temperature)
    s_trans = translational_entropy_cal(atomic_mass(symbol), temperature, pressure_atm)
    s_elec = electronic_entropy_cal(multiplicity)
    ts = temperature * (s_trans + s_elec) / 1000.0
    ss = standard_state_correction(temperature, pressure_atm) if standard_state else 0.0
    h_corr = 2.5 * rt
    return {
        "g_corr": h_corr - ts + ss,
        "h_corr": h_corr,
        "s_elec_cal": s_elec,
        "multiplicity": int(multiplicity),
        "zpe": 0.0,
        "e_vib": 0.0,
        "e_trans": 1.5 * rt,
        "e_rot": 0.0,
        "rt": rt,
        "ts_total": ts,
        "s_vib_cal": 0.0,
        "s_trans_cal": s_trans,
        "standard_state": ss,
        "temperature_K": float(temperature),
        "pressure_atm": float(pressure_atm),
        "n_real": 0,
        "n_imaginary": 0,
        "imaginary_cm1": [],
        "tr_max_cm1": 0.0,
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
