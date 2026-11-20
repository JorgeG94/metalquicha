"""Plottable infrared, Raman and UV-vis spectra from an output file.

    from mqc import spectra

    ir = spectra.ir("output_water.json", fwhm=20.0, scale=0.96)
    ir.to_csv("water_ir.csv")
    ir.plot("water_ir.png")                # needs matplotlib

    raman = spectra.raman("output_water.json", laser_nm=532.0)
    uv = spectra.uvvis("output_benzene.json", fwhm_ev=0.4, unit="nm")
    h1 = spectra.nmr("output_ethanol.json", "H", reference="output_tms.json")

or, with no Python written at all,

    python3 python/mqc/spectra.py ir output_water.json --csv water_ir.csv

A calculation writes *sticks*: a position and a strength per vibration or per
excited state. A spectrum is those sticks each given a line shape and added up,
with the strength converted to what an instrument reports. This module does
that and nothing else. It runs no calculation, imports nothing from the rest of
the package and needs only the standard library, so it works on a machine with
no build. Two things ask for more, and say so when it is missing: `plot` needs
matplotlib, and the exact simulation of NMR multiplets needs numpy.

What is read, from the JSON a run writes:

  ========  ==========================================  ======================
  spectrum  keys                                        from
  ========  ==========================================  ======================
  ir        ``frequencies_cm1``,                        ``driver: "Hessian"``
            ``ir_intensities_km_mol``
  raman     ``frequencies_cm1``,                        the same, with
            ``raman_activities_A4_amu``                 ``properties.raman``
  uvvis     ``excitation_energy_ev``,                   ``keywords.excited_states``
            ``oscillator_strength``
  nmr       ``isotropic`` and ``element`` per atom      ``properties.nmr``
  ========  ==========================================  ======================

Units of what comes back:

  * **ir**: molar absorptivity, L mol^-1 cm^-1, against cm^-1. The area under
    a band is its intensity: ``epsilon = 100 A / ln 10`` times a line shape of
    unit area, with ``A`` in km/mol.
  * **raman**: differential cross section per unit shift,
    cm^2 sr^-1 per cm^-1, against the Stokes shift in cm^-1. An *activity* is
    not an intensity: the cross section carries the laser frequency to the
    fourth power and the thermal population of the mode, see `raman`.
  * **uvvis**: molar absorptivity, L mol^-1 cm^-1, against eV, nm or cm^-1.
    The area under a band over wavenumber is ``f / 4.32e-9``.
  * **nmr**: nuclei per ppm against the chemical shift in ppm, so the area
    under a line is the number of nuclei in it. Shifts are against a
    reference shielding that has to be supplied, and with the couplings of a
    `properties.spin_spin_coupling` run the lines are split into their
    multiplets; see `nmr`.

The line width is not computed by anything here. It is a parameter standing for
everything the calculation left out -- lifetime, solvent, temperature, the
vibrational structure of an electronic band -- and the defaults (20 cm^-1 for a
vibration, 0.4 eV for an electronic band, 0.02 ppm for a nucleus) are
conventions, not predictions.
"""

import json
import math
import sys

__all__ = ["Spectrum", "SpectrumError", "ir", "raman", "uvvis", "nmr", "load"]


# -- constants --------------------------------------------------------------
#
# CODATA 2018, SI. Used for the Raman cross section and the oscillator-strength
# conversion, both of which are easier to get right in SI and convert once.

_H = 6.62607015e-34  #: J s
_KB = 1.380649e-23  #: J/K
_C = 299792458.0  #: m/s
_NA = 6.02214076e23  #: 1/mol
_E = 1.602176634e-19  #: C
_ME = 9.1093837015e-31  #: kg
_EPS0 = 8.8541878128e-12  #: F/m
_AMU_KG = 1.66053906660e-27

EV_TO_CM1 = _E / (_H * _C * 100.0)  #: 8065.54...
CM1_TO_KELVIN = _H * _C * 100.0 / _KB  #: hc/k, K per cm^-1

#: epsilon [L mol^-1 cm^-1] = IR_KM_MOL_TO_EPSILON * A [km/mol] * shape [1/cm^-1]
IR_KM_MOL_TO_EPSILON = 100.0 / math.log(10.0)

#: f = F_PER_AREA * integral of epsilon [L mol^-1 cm^-1] over wavenumber [cm^-1].
#: In SI the factor is 4 eps0 m_e c^2 ln10 / (N_A e^2), for epsilon in m^2/mol
#: over m^-1; one L mol^-1 cm^-1 cm^-1 is 10 m/mol.
F_PER_AREA = 4.0 * _EPS0 * _ME * _C**2 * math.log(10.0) / (_NA * _E**2) * 10.0

LINESHAPES = ("lorentzian", "gaussian")


class SpectrumError(ValueError):
    """The output file does not hold what the spectrum asked for needs."""


# ---------------------------------------------------------------------------
#  Line shapes
# ---------------------------------------------------------------------------


def lorentzian(x, centre, fwhm):
    """A Lorentzian of unit area and full width `fwhm` at half maximum."""
    half = 0.5 * fwhm
    return half / (math.pi * ((x - centre) ** 2 + half * half))


def gaussian(x, centre, fwhm):
    """A Gaussian of unit area and full width `fwhm` at half maximum."""
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    arg = (x - centre) / sigma
    return math.exp(-0.5 * arg * arg) / (sigma * math.sqrt(2.0 * math.pi))


def _shape(name):
    if name == "lorentzian":
        return lorentzian
    if name == "gaussian":
        return gaussian
    raise SpectrumError(
        "lineshape is %r; accepted: %s" % (name, ", ".join(LINESHAPES))
    )


def _grid(lo, hi, step):
    if not step > 0.0:
        raise SpectrumError("step must be positive, got %r" % (step,))
    if not hi > lo:
        raise SpectrumError("the range is empty: %r to %r" % (lo, hi))
    count = int(math.floor((hi - lo) / step + 1.0e-9)) + 1
    return [lo + k * step for k in range(count)]


def _broaden(positions, heights, grid, fwhm, lineshape):
    """Sum over sticks of height times a unit-area line at each grid point."""
    if not fwhm > 0.0:
        raise SpectrumError("the line width must be positive, got %r" % (fwhm,))
    shape = _shape(lineshape)
    return [
        sum(h * shape(x, p, fwhm) for p, h in zip(positions, heights)) for x in grid
    ]


# ---------------------------------------------------------------------------
#  The result
# ---------------------------------------------------------------------------


class Spectrum:
    """A broadened spectrum and the sticks it was made from.

    `x` and `y` are the curve. `sticks` is a list of ``(position, strength)``
    in the units of `stick_label`, positions on the same axis as `x`: these
    are the calculation's own numbers after any scaling, and what to quote.
    `warnings` lists anything dropped or assumed on the way.
    """

    def __init__(self, kind, x, y, sticks, x_label, y_label, stick_label, warnings):
        self.kind = kind
        self.x = x
        self.y = y
        self.sticks = sticks
        self.x_label = x_label
        self.y_label = y_label
        self.stick_label = stick_label
        self.warnings = warnings

    def __repr__(self):
        return "Spectrum(%s, %d sticks, %d points, %g to %g)" % (
            self.kind,
            len(self.sticks),
            len(self.x),
            self.x[0],
            self.x[-1],
        )

    def peaks(self):
        """Local maxima of the curve, as ``(x, y)``, tallest first."""
        found = []
        for k in range(1, len(self.y) - 1):
            if self.y[k] > self.y[k - 1] and self.y[k] >= self.y[k + 1]:
                found.append((self.x[k], self.y[k]))
        return sorted(found, key=lambda point: -point[1])

    def to_csv(self, path, sticks_path=None):
        """Write the curve as two columns, and optionally the sticks beside it."""
        with open(path, "w") as handle:
            handle.write("# %s\n" % self.kind)
            handle.write("%s,%s\n" % (_csv(self.x_label), _csv(self.y_label)))
            for x, y in zip(self.x, self.y):
                handle.write("%.6f,%.8e\n" % (x, y))
        if sticks_path is not None:
            with open(sticks_path, "w") as handle:
                handle.write("# %s, sticks\n" % self.kind)
                handle.write("%s,%s\n" % (_csv(self.x_label), _csv(self.stick_label)))
                for x, y in self.sticks:
                    handle.write("%.6f,%.8e\n" % (x, y))

    def plot(self, path=None, ax=None, sticks=True, peaks_down=None, **line_kwargs):
        """Draw the curve with matplotlib, with the sticks on a second axis.

        Returns the axes. `path` saves the figure; `ax` draws into axes of
        the caller's. A vibrational spectrum is drawn with the wavenumber
        falling to the right and an NMR one with the shift falling to the
        right, as the instruments draw them, and an infrared
        one with its bands hanging from the top, as a transmittance spectrum
        shows them. `peaks_down` overrides that either way; the quantity
        plotted is the absorptivity whichever way up it is drawn.
        """
        try:
            import matplotlib

            if ax is None and path is not None:
                matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            raise SpectrumError(
                "plot needs matplotlib, which is not installed. to_csv writes the "
                "same curve for any other plotting program."
            )
        if peaks_down is None:
            peaks_down = self.kind == "ir"
        figure = None
        if ax is None:
            figure, ax = plt.subplots(figsize=(7.0, 4.0))
        ax.plot(self.x, self.y, **line_kwargs)
        ax.set_xlabel(self.x_label)
        ax.set_ylabel(self.y_label)
        ax.set_ylim(bottom=0.0)
        if sticks and self.sticks:
            twin = ax.twinx()
            twin.vlines(
                [s[0] for s in self.sticks],
                0.0,
                [s[1] for s in self.sticks],
                colors="0.6",
                linewidth=1.0,
            )
            twin.set_ylabel(self.stick_label)
            twin.set_ylim(bottom=0.0)
            if peaks_down:
                twin.invert_yaxis()
            ax.set_zorder(twin.get_zorder() + 1)
            ax.patch.set_visible(False)
        ax.set_xlim(self.x[0], self.x[-1])
        if self.kind in ("ir", "raman", "nmr"):
            ax.invert_xaxis()
        if peaks_down:
            ax.invert_yaxis()
        if figure is not None:
            figure.tight_layout()
            if path is not None:
                figure.savefig(path, dpi=200)
                plt.close(figure)
        return ax


def _csv(label):
    return label.replace(",", ";")


# ---------------------------------------------------------------------------
#  Reading the output
# ---------------------------------------------------------------------------


def load(source):
    """The output as a dict, from a path or from a dict already loaded."""
    if isinstance(source, dict):
        return source
    with open(source) as handle:
        return json.load(handle)


def _find(node, key):
    """The first object named `key` anywhere in the output, depth first.

    An unfragmented run holds its sections one level down, under the deck's
    name; looking for the key rather than assuming the path keeps this
    working when that nesting is something else.
    """
    if isinstance(node, dict):
        if key in node and isinstance(node[key], dict):
            return node[key]
        for value in node.values():
            found = _find(value, key)
            if found is not None:
                return found
    elif isinstance(node, list):
        for value in node:
            found = _find(value, key)
            if found is not None:
                return found
    return None


def _vibrations(source, strength_key, what, hint, scale, min_frequency):
    """Scaled frequencies and strengths of the real vibrations, and warnings."""
    block = _find(load(source), "vibrational_analysis")
    if block is None:
        raise SpectrumError(
            'the output has no "vibrational_analysis" section. %s spectra need a '
            'run with "driver": "Hessian".' % what
        )
    if strength_key not in block:
        raise SpectrumError(
            'the vibrational analysis has no "%s". %s' % (strength_key, hint)
        )
    frequencies = block["frequencies_cm1"]
    strengths = block[strength_key]
    if len(frequencies) != len(strengths):
        raise SpectrumError(
            '"frequencies_cm1" and "%s" differ in length' % strength_key
        )
    if not scale > 0.0:
        raise SpectrumError("the frequency scale must be positive, got %r" % (scale,))

    # The list holds all 3N modes: the translations and rotations near zero,
    # and any imaginary mode as a negative number. Neither is a line.
    warnings = []
    imaginary = [f for f in frequencies if f < -min_frequency]
    if imaginary:
        warnings.append(
            "%d imaginary mode(s) left out (%s cm^-1): the geometry is not a "
            "minimum, and the spectrum is of a structure that is not one."
            % (len(imaginary), ", ".join("%.0fi" % -f for f in imaginary))
        )
    kept = [
        (f * scale, s) for f, s in zip(frequencies, strengths) if f >= min_frequency
    ]
    if not kept:
        raise SpectrumError(
            "no vibration above %g cm^-1 in the output" % min_frequency
        )
    return [k[0] for k in kept], [k[1] for k in kept], warnings


def _vibrational_grid(positions, xmin, xmax, step, fwhm):
    lo = 0.0 if xmin is None else xmin
    hi = max(positions) + 10.0 * fwhm + 100.0 if xmax is None else xmax
    if xmax is None:
        hi = math.ceil(hi / 100.0) * 100.0
    return _grid(lo, hi, step)


# ---------------------------------------------------------------------------
#  The three spectra
# ---------------------------------------------------------------------------


def ir(
    source,
    fwhm=20.0,
    scale=1.0,
    lineshape="lorentzian",
    xmin=None,
    xmax=None,
    step=1.0,
    min_frequency=10.0,
):
    """The infrared absorption spectrum of a Hessian run.

    `fwhm` is the full width at half maximum in cm^-1. `scale` multiplies
    every harmonic frequency, which is how the systematic overestimate of a
    method is usually taken out; the intensities are left alone. Modes below
    `min_frequency` -- the six near zero, and imaginary ones -- are left out.
    """
    positions, intensities, warnings = _vibrations(
        source,
        "ir_intensities_km_mol",
        "Infrared",
        "The method's Hessian path did not produce dipole derivatives.",
        scale,
        min_frequency,
    )
    grid = _vibrational_grid(positions, xmin, xmax, step, fwhm)
    heights = [IR_KM_MOL_TO_EPSILON * a for a in intensities]
    curve = _broaden(positions, heights, grid, fwhm, lineshape)
    return Spectrum(
        "ir",
        grid,
        curve,
        list(zip(positions, intensities)),
        "wavenumber / cm^-1",
        "molar absorptivity / L mol^-1 cm^-1",
        "intensity / km mol^-1",
        warnings,
    )


def raman_cross_section(frequency, activity, laser_cm1, temperature):
    """The Stokes differential cross section of one mode, in cm^2 sr^-1.

    For the usual arrangement -- light scattered at a right angle to a
    linearly polarized beam and to its polarization, both polarizations of
    the scattered light collected -- which is what the activity
    ``45 a'^2 + 7 gamma'^2`` is the invariant of:

        d sigma / d Omega = (2 pi)^4 / 45 * (nu0 - nu)^4
                            * h / (8 pi^2 c nu) * S / (1 - exp(-h c nu / k T))

    with the wavenumbers in cm^-1 and `activity` in A^4/amu. The fourth power
    is of the *scattered* frequency, and the last factor is the thermal
    population of the mode, which lifts the low-frequency lines.
    `temperature` zero or negative leaves the population factor out.
    """
    if not laser_cm1 > frequency:
        raise SpectrumError(
            "the laser at %.0f cm^-1 is below the mode at %.0f cm^-1: there is "
            "no Stokes line" % (laser_cm1, frequency)
        )
    # CGS, where the formula is usually written: the polarizability is a
    # volume, so the activity is cm^4/g.
    activity_cgs = activity * 1.0e-32 / (_AMU_KG * 1.0e3)
    h_cgs = _H * 1.0e7
    c_cgs = _C * 1.0e2
    amplitude = h_cgs / (8.0 * math.pi**2 * c_cgs * frequency)
    population = 1.0
    if temperature > 0.0:
        population = 1.0 / (1.0 - math.exp(-CM1_TO_KELVIN * frequency / temperature))
    return (
        (2.0 * math.pi) ** 4
        / 45.0
        * (laser_cm1 - frequency) ** 4
        * amplitude
        * activity_cgs
        * population
    )


def raman(
    source,
    laser_nm=532.0,
    temperature=298.15,
    fwhm=20.0,
    scale=1.0,
    lineshape="lorentzian",
    xmin=None,
    xmax=None,
    step=1.0,
    min_frequency=10.0,
    activities=False,
):
    """The Stokes Raman spectrum of a Hessian run that asked for `properties.raman`.

    The activities in the output are converted to cross sections at the laser
    wavelength `laser_nm` and the `temperature` in kelvin by
    `raman_cross_section`, then broadened. With ``activities=True`` that
    conversion is skipped and the activities themselves are broadened, which
    is what many published "calculated Raman spectra" are and is not what a
    spectrometer records: it overstates the high-frequency stretches against
    everything below them.

    The polarizability behind the activities is the static one, so the
    spectrum is an off-resonance one whatever `laser_nm` says: the laser
    enters through the fourth-power factor alone.
    """
    positions, activity, warnings = _vibrations(
        source,
        "raman_activities_A4_amu",
        "Raman",
        'Ask for them with "properties": {"raman": {}} on the Hessian run.',
        scale,
        min_frequency,
    )
    if not laser_nm > 0.0:
        raise SpectrumError("laser_nm must be positive, got %r" % (laser_nm,))
    grid = _vibrational_grid(positions, xmin, xmax, step, fwhm)
    if activities:
        heights = list(activity)
        y_label = "Raman activity / A^4 amu^-1 per cm^-1"
        stick_label = "activity / A^4 amu^-1"
    else:
        laser_cm1 = 1.0e7 / laser_nm
        heights = [
            raman_cross_section(f, s, laser_cm1, temperature)
            for f, s in zip(positions, activity)
        ]
        y_label = "cross section / cm^2 sr^-1 per cm^-1"
        stick_label = "cross section / cm^2 sr^-1"
    curve = _broaden(positions, heights, grid, fwhm, lineshape)
    return Spectrum(
        "raman",
        grid,
        curve,
        list(zip(positions, heights)),
        "Raman shift / cm^-1",
        y_label,
        stick_label,
        warnings,
    )


UVVIS_UNITS = ("ev", "nm", "cm-1")


def uvvis(
    source,
    fwhm_ev=0.4,
    lineshape="gaussian",
    unit="nm",
    gauge="length",
    shift_ev=0.0,
    xmin=None,
    xmax=None,
    points=2000,
):
    """The electronic absorption spectrum of a run with `keywords.excited_states`.

    Each state is broadened in *energy*, with full width `fwhm_ev`, whatever
    `unit` the axis is reported in: a band of constant width in eV is not one
    of constant width in nm. `unit` is ``"nm"``, ``"ev"`` or ``"cm-1"``, and
    `xmin`/`xmax` are in it. `gauge` picks the oscillator strength:
    ``"length"`` or ``"velocity"``. `shift_ev` is added to every excitation
    energy, for the rigid shift a method's systematic error is often
    corrected by.

    States with no oscillator strength in the output -- triplets out of a
    closed shell -- carry none here either and are left out.
    """
    if unit not in UVVIS_UNITS:
        raise SpectrumError(
            "unit is %r; accepted: %s" % (unit, ", ".join(UVVIS_UNITS))
        )
    if gauge not in ("length", "velocity"):
        raise SpectrumError("gauge is %r; accepted: length, velocity" % (gauge,))
    if points < 2:
        raise SpectrumError("points must be at least 2, got %r" % (points,))
    block = _find(load(source), "excited_states")
    if block is None or "states" not in block:
        raise SpectrumError(
            'the output has no "excited_states" section. A UV-vis spectrum needs '
            '"keywords": {"excited_states": {...}} on an energy run.'
        )
    key = "oscillator_strength" if gauge == "length" else "oscillator_strength_velocity"
    warnings = []
    energies, strengths = [], []
    dark = 0
    for state in block["states"]:
        if key not in state:
            dark += 1
            continue
        energies.append(state["excitation_energy_ev"] + shift_ev)
        strengths.append(state[key])
    if dark:
        warnings.append(
            "%d state(s) carry no oscillator strength in the output and are left "
            "out (spin-forbidden from this reference)." % dark
        )
    if not energies:
        raise SpectrumError(
            'no state in the output has "%s": nothing absorbs' % key
        )
    if min(energies) <= 0.0:
        raise SpectrumError("an excitation energy is not positive after the shift")
    warnings.append(
        "%d state(s) were computed, the highest at %.2f eV: absorption above that "
        "is missing, not absent." % (len(energies), max(energies))
    )

    # The grid is uniform in the unit asked for; the broadening is in energy.
    pad = 3.0 * fwhm_ev
    lo_ev = max(min(energies) - pad, 0.05)
    hi_ev = max(energies) + pad
    if unit == "ev":
        lo = lo_ev if xmin is None else xmin
        hi = hi_ev if xmax is None else xmax
    elif unit == "cm-1":
        lo = lo_ev * EV_TO_CM1 if xmin is None else xmin
        hi = hi_ev * EV_TO_CM1 if xmax is None else xmax
    else:
        lo = _ev_to_nm(hi_ev) if xmin is None else xmin
        hi = _ev_to_nm(lo_ev) if xmax is None else xmax
    if not hi > lo or not lo > 0.0:
        raise SpectrumError("the range is empty or not positive: %r to %r" % (lo, hi))
    grid = [lo + k * (hi - lo) / (points - 1) for k in range(points)]
    if unit == "ev":
        grid_ev = grid
    elif unit == "cm-1":
        grid_ev = [x / EV_TO_CM1 for x in grid]
    else:
        grid_ev = [_ev_to_nm(x) for x in grid]

    # The area under a band over wavenumber is f / F_PER_AREA, so a line of
    # unit area in eV is divided by the eV-to-wavenumber factor.
    heights = [f / F_PER_AREA / EV_TO_CM1 for f in strengths]
    curve = _broaden(energies, heights, grid_ev, fwhm_ev, lineshape)

    if unit == "ev":
        positions = energies
        x_label = "energy / eV"
    elif unit == "cm-1":
        positions = [e * EV_TO_CM1 for e in energies]
        x_label = "wavenumber / cm^-1"
    else:
        positions = [_ev_to_nm(e) for e in energies]
        x_label = "wavelength / nm"
    return Spectrum(
        "uvvis",
        grid,
        curve,
        list(zip(positions, strengths)),
        x_label,
        "molar absorptivity / L mol^-1 cm^-1",
        "oscillator strength",
        warnings,
    )


def _ev_to_nm(value):
    """eV to nm, which is its own inverse."""
    return 1.0e7 / (value * EV_TO_CM1)


def _shieldings(source, nucleus):
    """`(atom, isotropic shielding)` for every atom of element `nucleus`."""
    block = _find(load(source), "nmr_shielding")
    if block is None or "atoms" not in block:
        raise SpectrumError(
            'the output has no "nmr_shielding" section. An NMR spectrum needs '
            '"properties": {"nmr": {}} on the run.'
        )
    atoms = block["atoms"]
    if atoms and "element" not in atoms[0]:
        raise SpectrumError(
            "the shieldings in this output carry no element, so the nuclei "
            "cannot be picked out. Rerun with a current build."
        )
    wanted = nucleus.strip().lower()
    return [
        (entry.get("atom", index + 1), entry["isotropic"])
        for index, entry in enumerate(atoms)
        if entry["element"].strip().lower() == wanted
    ]


def nmr(
    source,
    nucleus="H",
    reference=None,
    groups=None,
    couplings=None,
    frequency_mhz=400.0,
    hetero=None,
    method="exact",
    fwhm=0.02,
    fwhm_hz=None,
    lineshape="lorentzian",
    xmin=None,
    xmax=None,
    points=None,
):
    """The NMR spectrum of one element: shifts, and multiplets when given couplings.

    A calculation gives shieldings; a spectrum is drawn in shifts,

        delta = sigma_reference - sigma

    so `reference` is required. It is either the reference compound's
    shielding in ppm, or the output of an NMR run on the reference compound
    -- tetramethylsilane for H, C and Si -- in which case the shieldings of
    its `nucleus` atoms are averaged. **The reference has to come from the
    same method and basis**: the errors of a shielding largely cancel in the
    difference, and only then.

    `groups` is a list of lists of atom numbers, counted from one as the
    output counts them, whose shieldings are averaged before anything else.
    That is how the three hydrogens of a rotating methyl group become one
    line: the calculation is of one frozen geometry and gives them three
    shifts, the experiment sees their average. Atoms in no group stand alone.

    Without `couplings` every line is a singlet: the positions and the
    areas are a spectrum's, its multiplets are not. `fwhm` is in ppm.

    `couplings` is the output of a run with `properties.spin_spin_coupling`
    on the same molecule with the same atom order -- it may be the same file
    as `source`, or one from another method, which is usual since couplings
    want a functional. The lines are then split. A splitting is in hertz and
    a shift in ppm, so the picture depends on the spectrometer:
    `frequency_mhz` is the frequency of the observed nucleus, 400 for
    hydrogen on a "400 MHz" instrument and 100.6 for carbon on the same one.
    The line width becomes `fwhm_hz`, one hertz unless given.

    `method="exact"` diagonalises the spin Hamiltonian of the observed
    nuclei, which gives the leaning of multiplets towards each other and
    every other second-order effect, needs numpy, and is limited to
    `MAX_EXACT_SPINS` nuclei. `method="first_order"` splits each line by its
    neighbours' counts, needs nothing, and is right only where shifts are far
    apart compared with the couplings. Either way a set in `groups` is
    treated as equivalent nuclei: one shift, and one coupling, the average,
    to each nucleus outside it.

    `hetero` lists other elements whose nuclei also split the lines, such as
    ``["F"]``. Only hydrogen, fluorine and phosphorus qualify: spin one half
    and the whole of their element. Couplings to carbon are left out, as 13C
    is one carbon in ninety.
    """
    if reference is None:
        raise SpectrumError(
            "nmr needs a reference: the reference compound's shielding in ppm, "
            "or the output of an NMR run on it at the same method and basis."
        )
    shieldings = _shieldings(source, nucleus)
    if not shieldings:
        raise SpectrumError("the output has no %s nucleus" % nucleus)
    warnings = []
    if isinstance(reference, (int, float)):
        sigma_ref = float(reference)
    else:
        in_reference = _shieldings(reference, nucleus)
        if not in_reference:
            raise SpectrumError("the reference output has no %s nucleus" % nucleus)
        values = [value for _, value in in_reference]
        sigma_ref = sum(values) / len(values)
        spread = max(values) - min(values)
        if spread > 0.05:
            warnings.append(
                "the %d %s shieldings of the reference differ by %.2f ppm and "
                "were averaged to %.3f ppm." % (len(values), nucleus, spread, sigma_ref)
            )

    by_atom = dict(shieldings)
    group_of = {}  # atom -> index of the set it is averaged with
    grouped = set()
    for index, group in enumerate(groups or []):
        missing = [atom for atom in group if atom not in by_atom]
        if missing or not group:
            raise SpectrumError(
                "group %r names atoms that are not %s nuclei of this output: %r"
                % (list(group), nucleus, missing)
            )
        if grouped.intersection(group):
            raise SpectrumError("an atom is in more than one group: %r" % (list(group),))
        grouped.update(group)
        mean = sum(by_atom[atom] for atom in group) / len(group)
        for atom in group:
            by_atom[atom] = mean
            group_of[atom] = index
    for atom, _ in shieldings:
        group_of.setdefault(atom, len(group_of) + len(groups or []) + atom)
    atoms = [atom for atom, _ in shieldings]
    shifts = {atom: sigma_ref - by_atom[atom] for atom in atoms}

    if couplings is None:
        lines = [(shifts[atom], 1.0) for atom in atoms]
        if fwhm_hz is not None:
            fwhm = fwhm_hz / _positive(frequency_mhz, "frequency_mhz")
        warnings.append(
            "no spin-spin couplings: every line is a singlet, and nuclei that "
            "exchange are separate lines unless grouped."
        )
        x_label = "%s chemical shift / ppm" % nucleus
    else:
        _positive(frequency_mhz, "frequency_mhz")
        if fwhm_hz is None:
            fwhm_hz = 1.0
        fwhm = fwhm_hz / frequency_mhz
        lines = _coupled_lines(
            atoms, shifts, group_of, couplings, nucleus, hetero, frequency_mhz, method, warnings
        )
        x_label = "%s chemical shift / ppm, at %g MHz" % (nucleus, frequency_mhz)

    # Lines that agree to well inside a line width are one line, so the
    # sticks read as the spectrum does.
    lines.sort()
    merged = []
    for shift, count in lines:
        if merged and abs(shift - merged[-1][0]) < 0.02 * fwhm:
            total = merged[-1][1] + count
            merged[-1] = ((merged[-1][0] * merged[-1][1] + shift * count) / total, total)
        else:
            merged.append((shift, count))
    positions = [line[0] for line in merged]
    counts = [line[1] for line in merged]

    pad = max(25.0 * fwhm, 0.1 * (max(positions) - min(positions)))
    lo = min(positions) - pad if xmin is None else xmin
    hi = max(positions) + pad if xmax is None else xmax
    if not hi > lo:
        raise SpectrumError("the range is empty: %r to %r" % (lo, hi))
    if points is None:
        # Enough to draw the narrowest line: eight points across its width.
        points = min(max(4000, int(8.0 * (hi - lo) / fwhm) + 1), 400000)
    if points < 2:
        raise SpectrumError("points must be at least 2, got %r" % (points,))
    grid = [lo + k * (hi - lo) / (points - 1) for k in range(points)]
    curve = _broaden(positions, counts, grid, fwhm, lineshape)
    return Spectrum(
        "nmr",
        grid,
        curve,
        list(zip(positions, counts)),
        x_label,
        "nuclei per ppm",
        "nuclei",
        warnings,
    )


def _positive(value, name):
    if not value > 0.0:
        raise SpectrumError("%s must be positive, got %r" % (name, value))
    return value


#: Elements whose observed isotope has spin one half and is all, or nearly
#: all, of the element: the ones that split a neighbour's line in an ordinary
#: spectrum. Carbon is absent because 13C is one nucleus in ninety.
SPIN_HALF_ABUNDANT = ("H", "F", "P")
MAX_EXACT_SPINS = 12
MAX_HETERO_PARTNERS = 8


def _coupling_matrix(couplings):
    """`{(atom, atom): J in Hz}` and `{atom: element}` from an output."""
    block = _find(load(couplings), "spin_spin_coupling")
    if block is None or "j_hz" not in block:
        raise SpectrumError(
            'the couplings output has no "spin_spin_coupling" section. Ask for it '
            'with "properties": {"spin_spin_coupling": {}} on the run.'
        )
    numbers = [entry["atom"] for entry in block["atoms"]]
    elements = {entry["atom"]: entry["element"].strip() for entry in block["atoms"]}
    table = {}
    for a, row in zip(numbers, block["j_hz"]):
        for b, value in zip(numbers, row):
            table[(a, b)] = value
    return table, elements


def _coupled_lines(atoms, shifts, group_of, couplings, nucleus, hetero, frequency, method, warnings):
    """`(shift in ppm, intensity)` of the coupled spectrum of `atoms`."""
    table, elements = _coupling_matrix(couplings)
    missing = [atom for atom in atoms if atom not in elements]
    if missing or any(elements[a].lower() != nucleus.strip().lower() for a in atoms):
        raise SpectrumError(
            "the couplings output does not describe the same molecule: its atoms "
            "and their elements differ from the shieldings'."
        )
    if method not in ("exact", "first_order"):
        raise SpectrumError("method is %r; accepted: exact, first_order" % (method,))

    # Nuclei averaged together are magnetically equivalent by assumption: one
    # shift, and one coupling to each nucleus outside the set. The couplings
    # inside the set then do not show and are left out.
    members = {}
    for atom in atoms:
        members.setdefault(group_of[atom], []).append(atom)

    def mean_j(set_a, set_b):
        return sum(table[(a, b)] for a in set_a for b in set_b) / (len(set_a) * len(set_b))

    sets = list(members.values())
    j_sets = [
        [0.0 if i == k else mean_j(sets[i], sets[k]) for k in range(len(sets))]
        for i in range(len(sets))
    ]
    shift_sets = [shifts[group[0]] for group in sets]

    # Nuclei of another element that are all spin one half. Their own
    # frequencies are megahertz away, so each only adds or subtracts half its
    # coupling to a nucleus's frequency, with equal odds.
    partners = []
    for element in hetero or []:
        name = element.strip().capitalize()
        if name.lower() == nucleus.strip().lower():
            continue
        if name not in SPIN_HALF_ABUNDANT:
            raise SpectrumError(
                "hetero element %r: only spin-one-half nuclei that make up their "
                "element are handled (%s)." % (element, ", ".join(SPIN_HALF_ABUNDANT))
            )
        partners += [atom for atom, e in elements.items() if e.capitalize() == name]
    j_hetero = [[mean_j(group, [partner]) for partner in partners] for group in sets]

    lines = []
    if method == "exact":
        if len(atoms) > MAX_EXACT_SPINS:
            raise SpectrumError(
                "%d coupled %s nuclei is %d spin states, too many to diagonalise "
                'here (the limit is %d nuclei). Use method="first_order".'
                % (len(atoms), nucleus, 2 ** len(atoms), MAX_EXACT_SPINS)
            )
        if len(partners) > MAX_HETERO_PARTNERS:
            raise SpectrumError(
                "%d coupled nuclei of other elements is %d simulations, too many "
                'here (the limit is %d). Use method="first_order".'
                % (len(partners), 2 ** len(partners), MAX_HETERO_PARTNERS)
            )
        spin_set = [index for index, group in enumerate(sets) for _ in group]
        n = len(spin_set)
        j_spins = [[j_sets[spin_set[a]][spin_set[b]] for b in range(n)] for a in range(n)]
        for state in range(2 ** len(partners)):
            offsets = [
                sum(
                    j_hetero[index][k] * (0.5 if state >> k & 1 else -0.5)
                    for k in range(len(partners))
                )
                for index in range(len(sets))
            ]
            larmor = [shift_sets[q] * frequency + offsets[q] for q in spin_set]
            for freq, weight in _exact_lines(larmor, j_spins):
                lines.append((freq / frequency, weight / 2 ** len(partners)))
    else:
        for index, group in enumerate(sets):
            pattern = [(shift_sets[index] * frequency, float(len(group)))]
            for other, partner in enumerate(sets):
                if other != index:
                    pattern = _split(pattern, j_sets[index][other], len(partner))
            for k in range(len(partners)):
                pattern = _split(pattern, j_hetero[index][k], 1)
            lines += [(freq / frequency, weight) for freq, weight in pattern]
        close = [
            (a, b)
            for a in range(len(sets))
            for b in range(a)
            if abs(j_sets[a][b]) > 0.1
            and abs(shift_sets[a] - shift_sets[b]) * frequency < 10.0 * abs(j_sets[a][b])
        ]
        if close:
            warnings.append(
                "first-order splitting was used where %d pair(s) of shifts are "
                "within ten couplings of each other: the real multiplets lean "
                'towards each other there, which method="exact" shows.' % len(close)
            )
    warnings.append(
        "couplings to 13C and other dilute or quadrupolar nuclei are left out, "
        "as they are from an ordinary spectrum."
    )
    return lines


def _split(lines, coupling, partners):
    """Each line split by `partners` equivalent spin-one-half nuclei."""
    if partners == 0 or coupling == 0.0:
        return lines
    out = []
    for freq, weight in lines:
        for up in range(partners + 1):
            share = math.comb(partners, up) / 2.0**partners
            out.append((freq + coupling * (up - 0.5 * partners), weight * share))
    return out


def _exact_lines(frequencies, couplings):
    """Transitions of coupled spin-one-half nuclei, as `(Hz, intensity)`.

    The Hamiltonian ``sum_i nu_i Iz_i + sum_{i<j} J_ij I_i . I_j`` conserves
    the total projection, so it is diagonalised one projection at a time and
    the lines are the transitions between neighbouring projections, each as
    strong as the square of the lowering operator's element between its two
    states. Normalised so that the intensities sum to the number of nuclei.
    """
    try:
        import numpy as np
    except ImportError:
        raise SpectrumError(
            "the exact spin simulation needs numpy, which is not installed. "
            'method="first_order" needs nothing.'
        )
    n = len(frequencies)
    blocks = {}
    for state in range(2**n):
        blocks.setdefault(bin(state).count("1"), []).append(state)

    def m(state, i):
        return 0.5 if state >> i & 1 else -0.5

    solved = {}
    for up, states in blocks.items():
        index = {state: k for k, state in enumerate(states)}
        ham = np.zeros((len(states), len(states)))
        for state, k in index.items():
            energy = sum(frequencies[i] * m(state, i) for i in range(n))
            for i in range(n):
                for j in range(i):
                    energy += couplings[i][j] * m(state, i) * m(state, j)
                    if (state >> i & 1) != (state >> j & 1):
                        flipped = state ^ (1 << i) ^ (1 << j)
                        ham[index[flipped], k] += 0.5 * couplings[i][j]
            ham[k, k] = energy
        solved[up] = (index,) + tuple(np.linalg.eigh(ham))

    lines = []
    for up in range(1, n + 1):
        upper_index, upper_e, upper_v = solved[up]
        lower_index, lower_e, lower_v = solved[up - 1]
        lowering = np.zeros((len(lower_index), len(upper_index)))
        for state, k in upper_index.items():
            for i in range(n):
                if state >> i & 1:
                    lowering[lower_index[state ^ (1 << i)], k] = 1.0
        strength = (lower_v.T @ lowering @ upper_v) ** 2
        for a in range(len(lower_e)):
            for b in range(len(upper_e)):
                if strength[a, b] > 1.0e-9:
                    lines.append((float(upper_e[b] - lower_e[a]), float(strength[a, b])))
    total = sum(weight for _, weight in lines)
    # Combination lines, where several nuclei flip at once, are allowed and
    # for coupled nuclei far apart in shift are a part in a billion: dropped.
    return [
        (freq, weight * n / total) for freq, weight in lines if weight * n / total > 1.0e-6
    ]


# ---------------------------------------------------------------------------
#  Command line
# ---------------------------------------------------------------------------


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(
        prog="spectra.py",
        description="Broaden the lines in a metalquicha output into a spectrum.",
    )
    parser.add_argument("kind", choices=("ir", "raman", "uvvis", "nmr"))
    parser.add_argument("output", help="the output_<name>.json of the run")
    parser.add_argument("--csv", help="write the curve here")
    parser.add_argument("--sticks", help="write the sticks here (needs --csv)")
    parser.add_argument("--plot", help="write a figure here (needs matplotlib)")
    parser.add_argument(
        "--peaks-up", action="store_true", help="draw infrared bands upwards"
    )
    parser.add_argument(
        "--fwhm",
        type=float,
        help="cm^-1 for ir/raman, eV for uvvis, ppm for nmr (Hz with --couplings)",
    )
    parser.add_argument("--nucleus", default="H", help="element of an nmr spectrum")
    parser.add_argument(
        "--reference",
        help="nmr: the reference shielding in ppm, or the reference's output file",
    )
    parser.add_argument(
        "--couplings", help="nmr: output with spin_spin_coupling; splits the lines"
    )
    parser.add_argument(
        "--frequency-mhz", type=float, default=400.0, help="nmr: of the observed nucleus"
    )
    parser.add_argument("--first-order", action="store_true", help="nmr: no numpy needed")
    parser.add_argument(
        "--hetero", action="append", default=[], help="nmr: another coupled element"
    )
    parser.add_argument(
        "--group",
        action="append",
        default=[],
        help="nmr: atoms averaged into one line, e.g. 3,4,5 (repeatable)",
    )
    parser.add_argument("--lineshape", choices=LINESHAPES)
    parser.add_argument("--scale", type=float, default=1.0, help="frequency scale")
    parser.add_argument("--laser-nm", type=float, default=532.0)
    parser.add_argument("--temperature", type=float, default=298.15)
    parser.add_argument("--activities", action="store_true")
    parser.add_argument("--unit", choices=UVVIS_UNITS, default="nm")
    parser.add_argument("--gauge", choices=("length", "velocity"), default="length")
    parser.add_argument("--shift-ev", type=float, default=0.0)
    parser.add_argument("--xmin", type=float)
    parser.add_argument("--xmax", type=float)
    args = parser.parse_args(argv)

    common = {"xmin": args.xmin, "xmax": args.xmax}
    if args.lineshape:
        common["lineshape"] = args.lineshape
    try:
        if args.kind == "nmr":
            reference = args.reference
            try:
                reference = float(reference)
            except (TypeError, ValueError):
                pass
            groups = [[int(atom) for atom in group.split(",")] for group in args.group]
            if args.fwhm:
                common["fwhm_hz" if args.couplings else "fwhm"] = args.fwhm
            result = nmr(
                args.output,
                args.nucleus,
                reference=reference,
                groups=groups,
                couplings=args.couplings,
                frequency_mhz=args.frequency_mhz,
                hetero=args.hetero,
                method="first_order" if args.first_order else "exact",
                **common
            )
        elif args.kind == "uvvis":
            if args.fwhm:
                common["fwhm_ev"] = args.fwhm
            result = uvvis(
                args.output,
                unit=args.unit,
                gauge=args.gauge,
                shift_ev=args.shift_ev,
                **common
            )
        else:
            if args.fwhm:
                common["fwhm"] = args.fwhm
            if args.kind == "ir":
                result = ir(args.output, scale=args.scale, **common)
            else:
                result = raman(
                    args.output,
                    scale=args.scale,
                    laser_nm=args.laser_nm,
                    temperature=args.temperature,
                    activities=args.activities,
                    **common
                )
        if args.csv:
            result.to_csv(args.csv, sticks_path=args.sticks)
        if args.plot:
            result.plot(args.plot, peaks_down=False if args.peaks_up else None)
    except (SpectrumError, OSError) as problem:
        sys.stderr.write("spectra: %s\n" % problem)
        return 1

    for warning in result.warnings:
        sys.stderr.write("note: %s\n" % warning)
    print("%-28s %s" % (result.x_label, result.stick_label))
    for position, strength in result.sticks:
        print("%-28.2f %.6g" % (position, strength))
    return 0


if __name__ == "__main__":
    sys.exit(main())
