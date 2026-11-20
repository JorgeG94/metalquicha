================
Plotting spectra
================

A calculation writes *sticks*: a position and a strength for each vibration or
excited state. ``mqc.spectra`` turns the sticks in an output file into a curve
that can be laid over a measured spectrum -- infrared, Raman, UV-vis or NMR --
with the strength converted to what an instrument reports.

It runs no calculation and changes nothing in the Fortran. It reads the JSON a
run wrote, imports nothing from the rest of the package and needs only the
Python standard library, so it works on a laptop with no build of the program.
Only ``plot`` needs matplotlib; ``to_csv`` writes the same curve for any other
plotting program.

.. code-block:: python

   from mqc import spectra

   ir = spectra.ir("output_water.json", fwhm=20.0, scale=0.96)
   ir.to_csv("water_ir.csv")
   ir.plot("water_ir.png")

   raman = spectra.raman("output_water.json", laser_nm=532.0)
   uv = spectra.uvvis("output_benzene.json", fwhm_ev=0.4, unit="nm")
   h1 = spectra.nmr("output_ethanol.json", "H", reference="output_tms.json")

or from a shell, with no Python written:

.. code-block:: bash

   python3 python/mqc/spectra.py ir output_water.json --csv water_ir.csv --plot water_ir.png
   python3 python/mqc/spectra.py raman output_water.json --laser-nm 785 --csv water_raman.csv
   python3 python/mqc/spectra.py uvvis output_benzene.json --unit ev --fwhm 0.3 --csv uv.csv
   python3 python/mqc/spectra.py nmr output_ethanol.json --nucleus H \
       --reference output_tms.json --group 5,6,7 --group 8,9 --plot ethanol_1h.png

The command prints the sticks and writes the curve where ``--csv`` says.
``python/examples/spectra.py`` runs a Hessian and an excited-state deck on water
and makes both spectra from them.

What each spectrum needs
------------------------

.. list-table::
   :header-rows: 1
   :widths: 12 44 44

   * - Spectrum
     - The run
     - Keys read
   * - ``ir``
     - ``"driver": "Hessian"``
     - ``frequencies_cm1``, ``ir_intensities_km_mol``
   * - ``raman``
     - the same, with ``"properties": {"raman": {}}``; see
       :doc:`vibrational_analysis`
     - ``frequencies_cm1``, ``raman_activities_A4_amu``
   * - ``uvvis``
     - ``"driver": "Energy"`` with ``keywords.excited_states``; see
       :doc:`excited_states`
     - ``excitation_energy_ev``, ``oscillator_strength``
   * - ``nmr``
     - ``"properties": {"nmr": {}}``, for the molecule and for the reference
       compound; see :doc:`nmr_shielding`
     - ``isotropic`` and ``element`` of each atom

An output without the section a spectrum needs is refused with the name of
what to ask for.

What comes back
---------------

Each function returns a ``Spectrum``:

- ``x``, ``y``: the curve.
- ``sticks``: ``(position, strength)`` per line, on the same axis as ``x``.
  These are the calculation's own numbers after any scaling, and what to quote.
- ``x_label``, ``y_label``, ``stick_label``: the units, as text.
- ``warnings``: anything left out or assumed, such as an imaginary mode.
- ``peaks()``: the maxima of the curve, tallest first.
- ``to_csv(path, sticks_path=None)`` and ``plot(path=None, ax=None)``.

``plot`` draws the curve with the sticks on a second axis, and a vibrational
spectrum with the wavenumber falling to the right, as the instruments draw it.
**An infrared spectrum is drawn with its bands hanging from the top**, the way
a transmittance spectrum shows them; what is plotted is still the
absorptivity, on an inverted axis. ``peaks_down=False`` (``--peaks-up`` on the
command line) draws it the other way, and ``peaks_down=True`` inverts any
other spectrum.
Pass ``ax`` to draw into a figure of your own, for instance to overlay an
experiment.

Infrared
--------

.. code-block:: python

   spectra.ir(source, fwhm=20.0, scale=1.0, lineshape="lorentzian",
              xmin=None, xmax=None, step=1.0, min_frequency=10.0)

The curve is the molar absorptivity in L mol\ :sup:`-1` cm\ :sup:`-1` against
wavenumber. An intensity :math:`A` in km/mol is the area of its band,

.. math::

   \varepsilon(\tilde\nu) = \frac{100}{\ln 10}\,\sum_i A_i\, g(\tilde\nu - \tilde\nu_i)

with :math:`g` a line shape of unit area, so a Lorentzian band of width
:math:`w` peaks at :math:`27.65\,A/w`.

``scale`` multiplies every frequency and leaves the intensities alone. Harmonic
frequencies are systematically too high -- by about ten per cent at
Hartree-Fock and three to five at a hybrid functional -- and a single factor
per method and basis is the usual correction.

The output lists all :math:`3N` modes. The six near zero are the translations
and rotations and are left out, as is anything below ``min_frequency``. **An
imaginary mode is left out too and reported in** ``warnings``: the geometry is
then not a minimum, and the spectrum is of a structure that is not one.

Raman
-----

.. code-block:: python

   spectra.raman(source, laser_nm=532.0, temperature=298.15, fwhm=20.0,
                 scale=1.0, lineshape="lorentzian", activities=False, ...)

**An activity is not an intensity.** What the run writes,
:math:`S_i = 45 a_i'^2 + 7\gamma_i'^2`, is a property of the molecule alone.
What a spectrometer records also depends on the laser and on how many
molecules are in each vibrational level, and for light scattered at a right
angle to a polarized beam it is the differential cross section

.. math::

   \frac{d\sigma_i}{d\Omega} = \frac{(2\pi)^4}{45}\,
      (\tilde\nu_0 - \tilde\nu_i)^4\,
      \frac{h}{8\pi^2 c\,\tilde\nu_i}\,
      \frac{S_i}{1 - e^{-hc\tilde\nu_i/kT}}

The curve is that, in cm\ :sup:`2` sr\ :sup:`-1` per cm\ :sup:`-1`, against
the Stokes shift. The :math:`1/\tilde\nu_i` and the thermal factor both lift
the low-frequency lines against the stretches, so a plot of the activities
themselves -- which ``activities=True`` gives, and which is what many
published "calculated Raman spectra" are -- overstates the C-H and O-H region.

**The spectrum is an off-resonance one whatever the laser.** The activities
come from the static polarizability, and the laser wavelength enters through
the fourth-power factor only. Near an electronic absorption the real spectrum
is a resonance Raman one and looks nothing like this.

UV-vis
------

.. code-block:: python

   spectra.uvvis(source, fwhm_ev=0.4, lineshape="gaussian", unit="nm",
                 gauge="length", shift_ev=0.0, xmin=None, xmax=None, points=2000)

The curve is the molar absorptivity against wavelength, energy
(``unit="ev"``) or wavenumber (``unit="cm-1"``). An oscillator strength is the
area of its band over wavenumber,

.. math::

   f = 4.32\times10^{-9} \int \varepsilon(\tilde\nu)\, d\tilde\nu

**The broadening is done in energy whatever the axis**: a band of constant
width in eV is not one of constant width in nm, and broadening in wavelength
makes the short-wavelength bands too narrow. A Gaussian is the default because
the width of an electronic band is mostly unresolved vibrational structure,
which a Lorentzian's long tails misrepresent.

Three things to keep in mind when reading the result:

- **It ends where the states end.** Only the states asked for are there, and
  ``warnings`` says where the highest one is; absorption above it is missing,
  not absent. Ask for more states than the region of interest needs.
- **Dark states are left out.** A triplet out of a closed shell has no
  oscillator strength in the output and contributes nothing.
- ``gauge="velocity"`` uses the velocity-gauge strengths. The two gauges agree
  in a complete basis, so the difference between the curves is a measure of
  the basis set.

``shift_ev`` adds a constant to every excitation energy, for the rigid shift a
method's systematic error is often corrected by.

NMR
---

.. code-block:: python

   spectra.nmr(source, nucleus="H", reference=None, groups=None, fwhm=0.02,
               lineshape="lorentzian", xmin=None, xmax=None, points=4000)

A calculation gives shieldings and a spectrum is drawn in shifts,

.. math::

   \delta = \sigma_\text{ref} - \sigma

so ``reference`` is required: either the reference compound's shielding in
ppm, or the output of an NMR run on the reference compound, whose ``nucleus``
shieldings are then averaged. For H, C and Si that compound is
tetramethylsilane. **Run the reference at the same method and basis as the
molecule.** A shielding's error is tens of ppm and mostly the same for every
nucleus of an element, so it cancels in the difference, and only when both
numbers carry it.

The curve is in nuclei per ppm against the shift, so the area under a line is
the number of nuclei in it, and the axis runs downfield to the left.

``groups`` averages nuclei that exchange. The calculation is of one frozen
geometry, in which the three hydrogens of a methyl group are in three
different places and have three shifts; the experiment sees one line at their
average, because the group rotates. List the atom numbers of each such set,
counted from one as the output counts them:

.. code-block:: python

   spectra.nmr("output_ethanol.json", "H", reference="output_tms.json",
               groups=[[5, 6, 7], [8, 9]])

Atoms in no group stand alone, and nuclei whose shifts coincide are reported as
one line with their count.

For ethanol at Hartree-Fock/cc-pVDZ against tetramethylsilane at the same
level this gives hydrogen lines at 1.24 ppm (three nuclei), 3.17 (two) and
0.25 (the hydroxyl, gas phase), and carbon lines at 17.5 and 54.1 ppm. The
measured values in solution are 1.2 and 3.7, and 18 and 58.

Two things this is not:

- **There are no spin-spin couplings, so every line is a singlet.** The
  positions and the areas are a spectrum's; the triplet and quartet of an
  ethyl group are not there. The couplings are a separate calculation that is
  not implemented.
- **The shieldings are closed-shell Hartree-Fock**, which is all
  ``properties.nmr`` computes. A hydroxyl or amine hydrogen is also moved by
  ppm by hydrogen bonding in solution, which a gas-phase molecule knows
  nothing about.

The line width
--------------

Nothing here computes a line width. It stands for everything the calculation
left out -- the lifetime, the solvent, the temperature, the vibrational
structure of an electronic band -- and the defaults, 20 cm\ :sup:`-1` for a
vibration, 0.4 eV for an electronic band and 0.02 ppm for a nucleus, are
conventions. The area under
each band does not depend on it; the peak height does, inversely.

What is not here
----------------

- **NMR multiplets.** There are no spin-spin couplings; see `NMR`_.
- **Anharmonicity.** Overtones, combination bands and Fermi resonances are
  absent; the frequencies are harmonic and ``scale`` is the only correction.
- **Vibronic structure** of electronic bands, emission, and circular dichroism.
- **Fragmented runs.** The vibrational analysis of an unfragmented calculation
  is what is read.
