=============================
pKa and the isoelectric point
=============================

``mqc.pka`` computes acid dissociation constants and the isoelectric point of a
small molecule from the free energies of its protonation states. It is a Python
workflow on top of the interface described in :doc:`python_interface`: every
calculation it makes is an ordinary ``mqc.MBE`` run, and nothing in the Fortran
changed for it.

It is built for molecules of thirty to forty atoms with several ionizable
sites, at GFN2-xTB with ALPB water at every stage. The one stage that is
expected to move to something better is the single point, and the design keeps
that a change to one dictionary; see `Swapping in a DFT single point`_.

.. code-block:: python

   import mqc
   from mqc import pka

   with mqc.session():
       acid = pka.Microstate("AcOH", pka.Geometry.from_xyz("acoh.xyz"),
                             charge=0, n_protons=1, site_class="carboxylic")
       base = pka.Microstate("AcO-", pka.Geometry.from_xyz("aco.xyz"),
                             charge=-1, n_protons=0)

       result = pka.run([acid, base])
       print(result.macro_pkas)            # uncalibrated: read the caveats first

       calibrated = result.with_calibration(
           [pka.Reference("AcOH", "AcO-", 4.76, site_class="carboxylic")])
       print(calibrated.to_json())

``python/examples/pka.py`` runs this end to end on acetic acid and on glycine,
whose pI it computes.

What it does, and what it does not
==================================

**The microstates are given, not enumerated.** You supply each protonation state
and each tautomer worth considering, as a geometry and a charge. Deciding that
glycine has a zwitterion, or that a given imidazole can be protonated on either
nitrogen, is a chemical judgement the workflow does not make, and a missing
microstate is a missing term in a partition sum rather than an error. Listing
them is the part of the calculation that is yours.

**It reports free energies, and then a calibrated pKa.** An uncalibrated pKa
from a semiempirical method is a number in the right unit and not much else; the
section below on calibration is the reason the workflow has the shape it has.

The protocol
============

Each microstate goes through four stages, each with its own level of theory:

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Stage
     - What it does
   * - ``conformers``
     - CREST samples conformers of the structure given, in solvent. The
       survivors inside an energy window (``energy_window_kcal``, default 3)
       and at most ``max_conformers`` (default 5) go forward.
   * - ``optimize``
     - Each goes to a minimum.
   * - ``frequencies``
     - A Hessian at each minimum, for the thermal correction.
   * - ``single_point``
     - The electronic energy at each minimum.

The free energy of a conformer is

.. math::

   G = E_\mathrm{sp} + G_\mathrm{corr}, \qquad
   G_\mathrm{corr} = \mathrm{ZPE} + E_\mathrm{vib} + E_\mathrm{trans}
     + E_\mathrm{rot} + RT - T S + RT\ln 24.46

and the conformers of a microstate are combined as
:math:`G_s = -RT \ln \sum_i e^{-G_i/RT}`. Every conformer is kept in the
result, with its energies, its correction, its weight and its imaginary
frequencies.

Each stage is a dictionary of ``mqc.MBE`` keyword arguments, and the default for
all four is GFN2 with ALPB water:

.. code-block:: python

   protocol = pka.Protocol(
       conformers={"method": "gfn2", "xtb": {"solvent": "water", "solvation_model": "alpb"}},
       energy_window_kcal=3.0,
       max_conformers=5,
   )

``conformers=None`` or ``optimize=None`` skips the stage and uses the structure
as given.

Two of the stages run the executable
------------------------------------

The library refuses ``driver="optimize"`` and ``driver="conformers"`` through a
session or the C interface: both drive ``run_calculation`` rather than being
driven by it. So those two stages write a deck under ``Protocol.workdir`` and run
the ``mqc`` executable on it, once per structure, and read back what it leaves:
``crest_conformers.xyz`` for the ensemble (the absolute energy on each comment
line is the refinement level's) and ``output_optimize_optimized.xyz`` for the
optimized geometry. The executable is ``Protocol.executable``, else
``$MQC_EXECUTABLE``, else ``mqc`` on the ``PATH``, else ``build/mqc``.

This is a limitation rather than a design, and it has consequences. The
executable is a separate single-rank program, which is what CREST requires, so
those stages do not use the ranks of your job; the MPI launcher's environment
variables are removed from the child's environment so it starts as a singleton.
``OMP_NUM_THREADS`` is passed through and CREST samples on that many threads.
The Hessians and single points do use the whole session.

Thermochemistry
===============

The Hessian run already writes a thermochemistry block, and the workflow reads
the translational, rotational and electronic terms from it as they are. It
recomputes the vibrational part, for three reasons.

**Quasi-RRHO entropy.** A flexible molecule has soft modes, and a harmonic
oscillator's entropy diverges as its frequency goes to zero while a real hindered
rotor's does not. Grimme's interpolation
[Chem. Eur. J. 18, 9955 (2012)] weights each mode between the oscillator and a
free rotor,

.. math::

   S = w(\nu)\,S_\mathrm{RRHO} + \bigl(1 - w(\nu)\bigr)\,S_\mathrm{FR}, \qquad
   w(\nu) = \frac{1}{1 + (\nu_0/\nu)^{\alpha}},

with :math:`\nu_0 = 100\ \mathrm{cm^{-1}}`, :math:`\alpha = 4`, and the free
rotor's moment of inertia capped at the average
:math:`B_\mathrm{av} = 10^{-44}\ \mathrm{kg\,m^2}`. Only the entropy is
interpolated; the enthalpy stays harmonic. A mode at 100 cm\ :sup:`-1` counts
half oscillator and half rotor, a mode at 500 cm\ :sup:`-1` is an oscillator to
a part in six hundred, and ``Protocol(qrrho=False)`` gives plain RRHO to compare
against. The thermal correction in ``thermal_corrections_hartree.to_gibbs`` in
the JSON output is plain RRHO and will not agree with the ``g_corr`` here; that is
expected. Moving the quasi-RRHO treatment into ``mqc_thermochemistry.f90`` is
noted there and in the Python source as a ``TODO(mqc)``.

**Translations and rotations are removed by count.** The vibrational analysis
returns :math:`3N` modes with the translations and rotations at or near zero,
and its imaginary-frequency count includes a ``-0.3`` cm\ :sup:`-1` rotational
residual. The workflow drops the six (five for a linear molecule) modes of
smallest magnitude as translations and rotations, and counts a negative one among
the rest as imaginary. The largest magnitude among the dropped modes is kept as
``tr_max_cm1``, which is the number to look at if you suspect a real soft mode was
absorbed by them.

**Imaginary frequencies are counted, never dropped silently.** Each conformer
carries ``n_imaginary`` and the values; the result's ``warnings`` list names every
microstate that has any. By default an imaginary mode contributes nothing, as in
the Fortran thermochemistry; ``Protocol(imaginary="flip")`` takes its absolute
value instead, the usual repair for the small imaginary mode of a loosely
converged structure. A large one means the structure is not a minimum and the
free energy is not to be used.

The standard state
------------------

The translational entropy is evaluated for a gas at 1 atm, and a pKa is for 1
mol/L. Compressing the ideal-gas molar volume (24.46 L at 298.15 K) to one litre
costs

.. math:: RT \ln\!\left(\frac{RT}{P}\right) = RT \ln 24.46 \approx 1.89\ \mathrm{kcal/mol}

and it is added to every solute. It is computed from the temperature and pressure
in the thermochemistry block, so a run at another temperature is not given the
298.15 K constant. ``Protocol(standard_state=False)`` leaves it out.

The term is the same for every microstate, so it cancels in every pKa difference
and is absorbed by the proton parameter; what it protects is the uncalibrated
absolute cycle, which uses literature values quoted for exactly this standard
state. GFN2's ALPB and the ``use_shift`` key of ``keywords.xtb`` may carry a
reference-state shift of their own; if both are on the solute is shifted twice, by
a constant, with the same harmlessness for relative pKas.

The proton
==========

A pKa is a difference of free energies that contains the free energy of the
proton in solution, and **no electronic-structure method here computes that
number**. It enters as a parameter, one per *site class*:

.. math:: \mathrm{p}K_a(s \to t) = \frac{G_t - G_s + G_{\mathrm{H}^+}(c)}{RT \ln 10}

where :math:`c` is the class of the proton that leaves.

Uncalibrated
------------

Without references the default is the literature cycle,
:math:`G_\mathrm{gas}(\mathrm{H^+}) = -6.28`, :math:`\Delta G_\mathrm{solv}(\mathrm{H^+}) = -265.9`
and the 1 atm to 1 M term, :math:`-270.3\ \mathrm{kcal/mol}` at 298.15 K. At xTB
level **expect it to be several pKa units off**, and anions worst: the
solvation of a localized negative charge is the weak point of ALPB, and the error
goes into :math:`G_t` directly. One pKa unit is :math:`RT\ln 10 = 1.36\ \mathrm{kcal/mol}`
at 298.15 K, which is a smaller error than any of the contributions above is
reliably accurate to. The result says so in its ``warnings``.

Calibrated
----------

``calibrate`` fits the shift of each site class from reference
``(acid, base, experimental pKa)`` triples: the shift is the *mean residual* of the
references on that class, converted to kcal/mol, and the slope stays one. A
fitted slope would make a transition's pKa depend on how many protons the
microstate carries, and break the consistency below.

One reference is a relative pKa: the shift is whatever makes that compound right,
and every other pKa on the class is then its offset from it, which is the quantity
the method is better at than the absolute. The example transfers the carboxylic
shift fitted on acetic acid to glycine's carboxyl group.

.. code-block:: python

   refs = [pka.Reference("AcOH", "AcO-", 4.76, site_class="carboxylic"),
           pka.Reference("PhOH", "PhO-", 9.99, site_class="phenol")]
   calibrated = result.with_calibration(refs)

``with_calibration`` returns a new result, so the uncalibrated numbers stay
available beside the calibrated ones.

Keeping the cycles consistent
-----------------------------

The shift belongs to the *protons a microstate carries*, not to a transition.
``site_class`` is one label per ionizable proton (a string applies to all of them),
and a microstate's grand potential at proton chemical potential :math:`\mu` is

.. math:: \Omega_s = G_s - n_s \mu - \sum_{c \in \mathrm{sites}(s)} \bigl(G_{\mathrm{H}^+}(c) - G_{\mathrm{H}^+}\bigr)

Micro pKas, macro pKas, populations and the charge curve are all built from the
same :math:`\Omega_s`, so they cannot disagree: going from the glycine cation to
the anion through the neutral tautomer costs exactly what it costs through the
zwitterion, and at ``pH = pKa(s -> t)`` the two microstates are equally
populated. That is checked in the test suite and again in the example. A tautomer
pair is where the labels matter: the neutral glycine carries its proton on the
carboxyl (``"carboxylic"``) and the zwitterion on the amine (``"ammonium"``).

Micro to macro, and the isoelectric point
=========================================

``micro_pkas`` lists the microscopic constant between every pair of microstates
that differ by one proton. ``macro_pkas`` gives the macroscopic constant between
protonation level :math:`n` and :math:`n-1` from partition sums over every
microstate at each level,

.. math::

   \mathrm{p}K_a(n \to n{-}1) = \frac{A_{n-1} + G_{\mathrm{H}^+} - A_n}{RT \ln 10},
   \qquad A_n = -RT \ln \sum_{s:\,n_s = n} e^{-G_s^\mathrm{eff}/RT}

which is the pH at which the two levels are equally populated. Two independent
identical sites of microscopic pKa :math:`p` give macroscopic constants
:math:`p - \log_{10} 2` and :math:`p + \log_{10} 2`, which the tests check.
Charge minus ``n_protons`` must be the same for every microstate, and levels
that are not consecutive give no macroscopic constant rather than an
interpolated one.

``populations(pH)`` and ``charge(pH)`` evaluate

.. math:: w_s(\mathrm{pH}) \propto e^{-\Omega_s(\mathrm{pH})/RT}, \qquad
          \mu(\mathrm{pH}) = G_{\mathrm{H}^+} - RT\ln 10\; \mathrm{pH}

with log-sum-exp, since the absolute free energies are of order :math:`10^5`
kcal/mol. The isoelectric point is the root of ``charge(pH)`` on [0, 14] by
bisection. **When the net charge does not change sign on that interval
there is no pI and the result says so**: ``pI`` is ``None``, ``pI_note`` gives the
reason, and ``isoelectric_point()`` raises ``NoIsoelectricPoint``. A monoprotic
acid is the usual case (charges 0 and -1); an uncalibrated molecule whose pKas
fall off the end of the scale is the other, and returning an edge of the interval
there would read as an answer. For a three-state cation, zwitterion and anion
the pI is exactly the mean of the two macroscopic pKas.

Caveats
=======

**GFN2 does not order tautomers and zwitterions reliably.** Whether glycine's
zwitterion is below its neutral form in implicit water is a balance of a few
kcal/mol, and a semiempirical method's error on it is as large. The result's
``warnings`` flags every pair of microstates at the same protonation level within
3 kcal/mol, the range in which the ordering is not a ranking. A populations table
that depends on that ordering inherits its uncertainty.

**Anions and ALPB.** Carboxylates, phenolates and deprotonated amines are where
implicit solvation does worst, because the first solvation shell is a set of
hydrogen bonds the continuum does not have. Calibrate per site class, and do not
transfer a shift between classes that solvate differently.

**Imaginary modes.** A Hessian is only a minimum's if the optimization converged;
the workflow refuses an unconverged optimization unless
``Protocol(allow_unconverged=True)``, and reports every imaginary frequency it
sees.

**Conformers are selected before the free energies exist.** The window applies to
the conformer search's refinement energy, not to :math:`G`, and two conformers that
optimize to the same minimum are merged (they would otherwise be counted twice,
lowering :math:`G` by :math:`RT \ln 2`). Enantiomeric conformers are not
merged.

**The symmetry number is whatever the thermochemistry block used**, which is 1
unless the deck says otherwise; the workflow does not detect it. Rotational
symmetry changes a free energy by :math:`RT\ln\sigma`, a few tenths of a kcal/mol
for a methyl group, and is the same on both sides of a pKa when the group is
unchanged.

**Closed shells only in practice.** ``Microstate`` takes a ``multiplicity`` and the
electronic entropy term uses it, but nothing here has been exercised on an open
shell.

Swapping in a DFT single point
==============================

The single point is the only stage whose result needs to be better than xTB, and
the structures and frequencies stay at GFN2. Replace one entry:

.. code-block:: python

   protocol = pka.Protocol(
       single_point={
           "method": "dft", "functional": "pbe0", "basis": "def2-svp",
           "pcm": {"method": "iefpcm", "dielectric": 78.3553},
       },
   )

The free energy becomes :math:`E_\mathrm{DFT} + G_\mathrm{corr}^\mathrm{GFN2}`, where
the correction is the frequency stage's and so carries no DFT contribution. When
the single point is not the same calculation as the Hessian, the workflow runs it
as its own energy calculation at each conformer's geometry; when the two are equal
(the default) the Hessian run's electronic energy is reused rather than computed a
second time.

Two things change with it. The solvation energy is now the continuum model's, and
``keywords.pcm`` and ``keywords.xtb`` have different conventions about the
reference state (see :doc:`continuum_solvation`), so re-calibrate rather than
carrying shifts over from a GFN2 run. And the electronic energy now contains a
solvation free energy from a different model than the one the structures were
relaxed in, which is normal practice and worth knowing about when a conformer's
ranking moves.

Reading the result
==================

``PKaResult`` carries, for each microstate, the Boltzmann-combined free energy, the
number of conformers and the imaginary-frequency count, with the per-conformer
record under ``detail["conformers"]``. ``micro_pkas`` and ``macro_pkas`` are lists
of dictionaries, ``pI`` is a number or ``None``, and ``to_json()`` returns all of it
(and writes it, given a path). Labels of the runs the workflow makes are
``<prefix>_<microstate>_c<k>_freq`` and ``..._sp``, with anything not a letter, digit,
underscore or hyphen in the microstate's name replaced by an underscore, and their
``output_*.json`` files are left in the working directory as the record of what was
computed.

Run it inside ``mqc.session()`` and only once: after the session is entered only
rank 0 continues, and the calls are made one after another.
