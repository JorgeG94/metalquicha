Effective Fragment Potentials (MAKEFP)
======================================

The ``MakeFP`` driver builds an **effective fragment potential** for a molecule and
writes it as a ``.efp`` file. It computes no interaction energy of its own: the
output is a set of parameters describing one fragment, to be used afterwards --
either by another program or by a later fragment calculation here.

This is GAMESS's ``RUNTYP=MAKEFP``, and the file it writes is in GAMESS's format,
so a potential produced here can be handed straight to GAMESS.

Running one
-----------

.. code-block:: json

   {
     "schema": {"name": "water_makefp", "version": "1.0"},
     "molecules": [{
       "xyz": "water.xyz",
       "molecular_charge": 0,
       "molecular_multiplicity": 1
     }],
     "model": {"method": "hf", "basis": "6-31g*"},
     "driver": "MakeFP"
   }

.. code-block:: bash

   ./mqc water_makefp.json

The potential is written beside the deck, named after it -- ``water_makefp.json``
leaves ``water_makefp.efp`` -- in the same way the JSON output is placed. ``MakeEFP``
is accepted as a spelling of the driver as well.

Only rank zero does the work. It is one SCF and a set of response solves, all of
them threaded inside the integral backend, so distributing it would mean every rank
recomputing the same potential to overwrite the same path. Run it under ``mpirun``
if you like; the extra ranks return immediately.

Example decks live in ``validation/inputs/cpu/mqc/makefp/``.

What the file contains
----------------------

All seventeen sections GAMESS's reader recognises:

* **Electrostatics** -- distributed multipoles through the octupole, at the atoms
  and at bond midpoints, with charge-penetration screening fitted in both the
  exponential and Gaussian damping forms
* **Polarization** -- static and dynamic (imaginary-frequency) dipole
  polarizability tensors, distributed over localized orbitals
* **Dispersion** -- the ``E6``, ``E7`` and ``E8`` coefficients, from the
  dipole--dipole, dipole--quadrupole and quadrupole--quadrupole responses at twelve
  Casimir--Polder frequencies
* **Exchange repulsion and charge transfer** -- the localized orbitals, their Fock
  matrix, the canonical orbitals and the basis set itself. These terms are not
  computed here at all; they are formed at use time from overlaps between fragments,
  so what the file carries is the data to form them from

Basis set requirements
----------------------

The orbital basis must be **Cartesian**, with shells no higher than ``f``. GAMESS
writes and reads projection data in the Cartesian basis, and the ordering and
normalization of each shell had to be matched against its output shell by shell.
``s``, ``p``, ``d`` and ``f`` are mapped; ``g`` and above are refused with a clear
message rather than written wrongly.

In practice this means a Pople set -- ``6-31g*``, ``6-31g(2df,p)`` -- and not a
Dunning one. A spherical basis is refused, and Dunning sets are spherical here
because that is how they are defined.

Density fitting
---------------

The dynamic response is what a potential costs: building its Hessian exactly takes
one Fock build per occupied-virtual pair. Fitting the integrals replaces that with
two matrix products, which at cc-pVQZ is 268 seconds against 0.074.

It reuses the SCF's own keywords, meaning the same thing here as there:

.. code-block:: json

   "model": {
     "method": "hf",
     "basis": "6-31g*",
     "aux_basis": "def2-universal-jkfit"
   },
   "keywords": {"scf": {"density_fitting": true}}

The auxiliary basis is read in whatever angular form the orbital basis is in, so a
spherical fitting set can serve a Cartesian orbital basis. What that costs, measured
against the exact potential from the same deck: the multipoles, the projection data,
the Fock matrix and both screening fits come out bit-identical, the dynamic
polarizabilities differ by about ``1e-04`` relative, and GAMESS's own energies from
the two files differ by ``5e-08`` Hartree on the dispersion total -- about
``3e-05`` kcal/mol. Nothing else moves, none of the other terms depending on the
response.

Fitting is an approximation and is therefore requested rather than inferred.

Tuning the response
-------------------

``keywords.efp`` holds what the stages after the SCF read. The SCF itself is
configured by ``keywords.scf`` and reads it here exactly as on any other driver;
there is deliberately no second spelling of a tolerance in this group.

.. code-block:: json

   "keywords": {
     "efp": {
       "dynamic_tolerance": 1.0e-7,
       "dynamic_maxiter": 200,
       "response": "auto",
       "vdw_scale": 0.7
     }
   }

.. list-table::
   :header-rows: 1
   :widths: 24 16 60

   * - Key
     - Default
     - Meaning
   * - ``dynamic_tolerance``
     - ``1.0e-7``
     - Residual the frequency-dependent response solve converges to, relative to
       its right-hand side
   * - ``dynamic_maxiter``
     - ``200``
     - Iterations that solve gets per system before it declines to converge
   * - ``response``
     - ``"auto"``
     - ``"auto"``, ``"dense"`` or ``"matrix_free"``; how that solve is carried out
   * - ``vdw_scale``
     - ``0.7``
     - Innermost layer of the charge-penetration screening grid, as a fraction of a
       van der Waals radius. GAMESS's ``VDWSCL``

**Only the first of these changes the answer.** The dynamic polarizabilities, and
every dispersion coefficient taken from them, are converged to
``dynamic_tolerance`` and no further, so loosening it moves the numbers in the
file. The EFMO literature runs the equivalent solves at ``5e-5`` and reports about
``3e-6`` Hartree of error in a total interaction energy for it -- Sattasathuchana
et al., JCTC **20**, 2445 (2024) -- but that measurement is of their accumulation,
not this one, which is the reason the default here is still ``1e-7``. The other
three change only how the answer is reached.

``response`` picks between two routes to the same equations. The dense one builds
the ``(A+B)`` and ``(A-B)`` operators and factorises them, which costs one Fock
build per occupied-virtual pair and then makes every frequency and perturbation
free; the matrix-free one never forms them, and its cost is the iteration count
instead. ``"auto"`` decides on size alone: dense unless three ``n_ov^2`` matrices
would be too large to hold. That threshold is a judgement about memory rather than
about wall clock, so the two forcing values exist to let a molecule be run both
ways and timed. On water in ``6-31g`` the two routes agree to ``1.8e-07`` on every
dynamic block and leave every other section of the file bit-identical, which is
what a comparison of the two should look like.

A forced ``"matrix_free"`` declines the auxiliary basis along with the operator:
there is no Hessian for it to fit. The run says so rather than fitting nothing
quietly.

Agreement with GAMESS
---------------------

For a water dimer built from one potential, against GAMESS's own MAKEFP output for
the same molecule and basis, comparing the energies GAMESS itself computes from each
file:

.. list-table::
   :header-rows: 1
   :widths: 40 30

   * - Term
     - Difference (Hartree)
   * - Electrostatics
     - ``1e-10``
   * - Polarization
     - ``1e-09``
   * - Exchange repulsion
     - exact
   * - Dispersion ``E6``
     - ``2e-09``
   * - Dispersion ``E7``
     - ``4e-10``
   * - Dispersion ``E8``
     - ``2e-07``
   * - Charge transfer
     - exact

The two terms that agree exactly are the ones the file only carries data for rather
than computing, so GAMESS forms them itself from identical inputs.

One known difference is the charge-penetration screening fit. Ours and GAMESS's
disagree by about 0.7 kcal/mol on a water hexamer, and the disagreement is genuine
rather than a bug in the transfer: our fit finds a real minimum of the objective
where GAMESS's sits on the bound of its own search range.
