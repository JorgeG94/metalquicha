Continuum solvation (PCM)
=========================

A molecule in solution sits in a polarizable medium: its charge distribution
polarizes the solvent, and the polarization acts back on the molecule. A
polarizable continuum model (PCM) replaces the solvent with a dielectric
continuum outside a molecule-shaped cavity, puts apparent charges on the cavity
surface, solves for them self-consistently with the SCF density, and adds their
potential to the Fock matrix and their interaction energy -- with its
variational factor of one half -- to the total.

The keyword reference lives with the other keywords in :doc:`input_files`
(``keywords.pcm``). This page is what the model is, which backend does what,
why the cavity convention matters more than anything else in it, and how each
implementation was validated.

Why you would use it
--------------------

The sharpest case is an ion that does not exist in the gas phase. A 51-atom
dianion (C\ :sub:`19`\ H\ :sub:`16`\ N\ :sub:`8`\ O\ :sub:`8`, charge -2) at
PBE0/cc-pVDZ is electronically *unbound* in vacuum: six of its 126 occupied
orbitals come out with positive orbital energy and the HOMO sits at +0.045
Hartree. That is not a convergence problem to push through -- adding diffuse
functions lowers the energy by another 99 kcal/mol, which is the signature of
electrons escaping toward the continuum limit, not of a basis converging. Every
gas-phase property of that wave function is an artifact of the basis being too
small to let the electrons leave.

With a water cavity (IEF-PCM, :math:`\varepsilon` = 78.3553) the same molecule
is entirely ordinary:

.. list-table::
   :header-rows: 1

   * - Quantity
     - Gas phase
     - IEF-PCM water
   * - HOMO (Hartree)
     - +0.045
     - -0.203
   * - Occupied orbitals above zero
     - 6
     - 0
   * - HOMO-LUMO gap (eV)
     - 0.4
     - 3.39
   * - Solvation energy (kcal/mol)
     - --
     - about -200

A deck for it, on the CPU backend:

.. code-block:: json

   {
     "molecules": [
       { "xyz": "mol.xyz", "molecular_charge": -2, "molecular_multiplicity": 1 }
     ],
     "model": { "method": "dft", "functional": "pbe0", "basis": "cc-pvdz" },
     "keywords": {
       "scf": { "maxiter": 60, "tolerance": 1e-08 },
       "dft": { "grid_level": 3 },
       "pcm": { "method": "iefpcm", "dielectric": 78.3553, "angular_points": 302 }
     },
     "driver": "Energy"
   }

The same gas-phase SCF oscillates for as long as it is allowed to run; the
solvated one converges in 13 iterations, and the continuum adds almost nothing
to the iteration cost -- the surface-charge equations are factorized once per
geometry and each iteration solves them directly.

Two backends, two implementations
---------------------------------

**CPU (libcint) backend.** Implements the smooth switching/Gaussian (SWIG)
discretization of Lange and Herbert [J. Chem. Phys. 133, 244111 (2010)],
following PySCF's ``pyscf.solvent.pcm`` term for term: Lebedev points per
atomic sphere, the quintic switching function with buried points discarded,
per-point Gaussian charges with the fitted exponents of
[J. Chem. Phys. 122, 194110 (2005)], and the same S and D matrices. Two models
are solved, chosen by ``keywords.pcm.method``:

- ``cpcm`` -- conductor-like; the charge equation is scaled by
  :math:`(\varepsilon-1)/\varepsilon`.
- ``iefpcm`` -- the integral equation formalism, whose D-matrix terms make it
  exact for ideal cavities; internal scale :math:`(\varepsilon-1)/(\varepsilon+1)`.

Restricted and unrestricted, Hartree-Fock and any non-double-hybrid
functional, with or without density fitting. **Energies only** -- see the
refusals below.

**cuEST (GPU) backend.** The cavity and the charge solve belong to the cuEST
library (an iSWIG surface with a preconditioned conjugate-gradient solve); mqc
supplies the radii, the dielectric and the exponents. One fixed model -- a deck
asking for ``iefpcm`` there is refused rather than quietly given the other
model's charges. The GPU path has the PCM nuclear gradient; the CPU path does
not.

Both backends build their cavity from the same radii table
(``mqc_pcm_radii``: Bondi, filled in from Mantina, scaled by
``radii_scale`` = 1.2), so a deck moved between them describes the same cavity.

The cavity decides the answer
-----------------------------

The solvation energy scales roughly as :math:`1/R` in the cavity radii, so the
radii are worth more than everything else in the model combined -- and codes
disagree about them *by default*. mqc uses Bondi's van der Waals radii
(hydrogen 1.20 Angstrom) scaled by 1.2. PySCF and gpu4pyscf default to a
"modified Bondi" table whose only change is hydrogen at **1.1** Angstrom, times
the same 1.2.

That one radius is not a detail. On the dianion above (16 hydrogens), mqc
against gpu4pyscf's *default* cavity differs by about 5 mHartree in total
energy and 0.04 eV in the gap -- numbers that look exactly like a bug in
whichever code is being checked. It is not a bug; it is two defensible
hydrogen radii. To make PySCF or gpu4pyscf build **mqc's** cavity:

.. code-block:: python

   from pyscf.data import radii
   cm.radii_table = 1.2 * radii.VDW    # Bondi, H = 1.20 A -- matches mqc
   cm.lebedev_order = 29               # 302 points, matches angular_points: 302

With the cavity matched, the two implementations agree to the ninth decimal
place (below). Without it, no amount of grid tightening will close the gap.

Validation
----------

**CPU backend, against PySCF** (same model, same dielectric, same cavity radii,
same Lebedev order; basis exponents read from ``basis_sets/`` on both sides,
because PySCF's internal tables carry more digits than BSE writes and that
alone is 2.5e-8 on water):

.. list-table::
   :header-rows: 1

   * - Case
     - Model
     - Agreement (Hartree)
   * - Water, RHF/STO-3G
     - C-PCM
     - 3e-12
   * - Hydroxide anion, RHF/6-31G
     - IEF-PCM
     - 5e-12
   * - Water, PBE0/cc-pVDZ
     - IEF-PCM
     - 3e-8 (grid quadrature, not PCM)
   * - H\ :sub:`2`\ O\ :sup:`+`, UHF/6-31G
     - IEF-PCM
     - 9e-9
   * - The 51-atom dianion, PBE0/cc-pVDZ
     - IEF-PCM
     - 1.2e-9 vs gpu4pyscf, matched cavity

The first two run in the CPU validation suite
(``validation/validation_tests_cpu.json``, at the manifest's global 1e-9). The
physics that needs no external reference is pinned in
``test/test_mqc_pcm.f90``: an isolated sphere's tessellation integrates to
exactly :math:`4\pi R^2`, buried points are accounted for one by one, and a
central charge in a spherical cavity reproduces Born's closed-form energy --
where IEF-PCM must land on the *true* :math:`(\varepsilon-1)/\varepsilon` Born
energy despite its internal :math:`(\varepsilon-1)/(\varepsilon+1)` factor,
which is the D-matrix terms doing their job and the first thing to break if
they are misscaled.

**cuEST backend, against PySCF's C-PCM**, on water in water at PBE0/def2-SVP
and :math:`\varepsilon` = 78.39, at the defaults:

.. list-table::
   :header-rows: 1

   * - Quantity
     - metalquicha (cuEST)
     - PySCF C-PCM
   * - dielectric energy
     - -10.2104 mHartree
     - -10.3625 mHartree
   * - solvation energy
     - -5.804 kcal/mol
     - -5.904 kcal/mol
   * - solvated dipole
     - 2.21754 D
     - 2.2162 D

The dipole is the sharp one: agreeing to 0.0013 D means the potential reaching
the Fock matrix is right, not merely that an energy came out plausible. The
residual 1.5% on the dielectric energy is the difference between two
independently discretised cavities and is not expected to close -- cuEST's
surface construction is its own, where the CPU backend's is PySCF's by design,
which is why the CPU numbers above agree to twelve digits and these agree to
three.

Two limits of the cuEST path were checked as well, and both are worth knowing
because they need no external reference. At :math:`\varepsilon` = 1.0001 the
dielectric energy falls to -7e-7 Hartree and the total returns to the
gas-phase energy, agreeing with PySCF's to 5.5 microhartree -- the continuum
contributes nothing when it should contribute nothing. And the surface
resolution matters more than it looks: at 110 angular points the dielectric
energy was 21% short, which is why ``angular_points`` defaults to 302.

Fukui indices in solution
-------------------------

A Fukui analysis is three SCFs -- the neutral and the two ions -- and all three
are solved in the same continuum. Ask for ``keywords.pcm`` and
``properties.fukui`` together and the ions inherit the neutral's cavity; there
is nothing further to switch on.

.. code-block:: json

   {
     "model": {"method": "hf", "basis": "6-31g"},
     "keywords": {
       "scf": {"maxiter": 100, "tolerance": 1e-12},
       "pcm": {"method": "cpcm", "dielectric": 78.3553, "angular_points": 302}
     },
     "driver": "Energy",
     "properties": {"fukui": {"population": "chelpg"}}
   }

**This used to be refused, and the reason it was is worth keeping in mind.**
The ions were solved in the gas phase while the neutral carried the continuum,
which makes the ionisation potential a difference across two different physics.
That is not a small error: on water in 6-31G the solvent moves the IP by
0.131 Hartree and the electron affinity by 0.114 Hartree -- 3.6 and 3.1
electronvolts, the size of the quantities themselves rather than a correction
to them. Both move the way they should, the cation and the anion each being
stabilised by the dielectric.

**The result is the adiabatic quantity, not the vertical one.** All three states
are solved with equilibrium solvation, so the solvent is allowed to relax around
each charge state. A vertical ionisation potential -- what a photoelectron
experiment measures, and what Fukui indices are conventionally taken to be --
would freeze the slow orientational polarisation at the neutral's value and let
only the fast electronic part respond, through an optical dielectric that
``pcm_config_t`` does not carry. For a redox potential or a solution-phase
reactivity index the adiabatic number is the one wanted; for a vertical
spectroscopic quantity it is not, and the difference in a polar solvent is a
few tenths of an electronvolt.

The cavity is built once from the geometry and reused by all three states,
since the nuclei do not move between them.

What is refused, and why
------------------------

These are deliberate refusals, not oversights. Each would otherwise return a
converged, plausible number computed in a quietly mixed phase, with nothing in
the output to say so.

- **The PCM nuclear gradient (CPU backend).** The energy's derivative carries
  the cavity moving with the nuclei and the charges responding, and neither
  term is built. An analytic gradient without them would be missing the
  solvent's pull on every atom. Run the energy on the CPU, or the gradient on
  cuEST, which has it.
- **Correlation on a solvated reference** -- MP2, coupled cluster, and the
  double hybrids, whose perturbative term is an MP2 by another spelling. The
  orbitals would carry the continuum while the correlation treatment ignored
  it; whether that is the PTE scheme or an inconsistency should be a decision,
  not an accident.
- **The bonding energy decomposition**, which rebuilds atom energies from
  gas-phase operators and would drop the dielectric term without saying so.
- **The analytic Hessian** quietly stands down rather than refusing: a solvated
  run falls through to the finite-difference path, which is *correct* for the
  continuum -- each displaced energy rebuilds its own cavity.

The xTB path's CPCM (``keywords.xtb``) is tblite's own model with its own
cavity and shares nothing with ``keywords.pcm``; see :doc:`input_files`.
