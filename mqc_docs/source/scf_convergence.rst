================
SCF Convergence
================

An SCF that does not converge is not one problem. It is at least three, and they
want different answers: the iteration is oscillating between two densities, the
iteration is crawling because the gap is small, or the state being asked for does
not exist. Only the first two are convergence problems. This page covers the aids
Metalquicha has for them, and how to tell the third case apart before spending a
week on it.

What is available
-----------------

DIIS is on by default and does most of the work. Beyond it there are two
energy-based accelerators, EDIIS and ADIIS, which open an SCF that DIIS opens
badly and then hand back to it; and there is level shifting. Both are described
below. There is no damping, no Fermi smearing and no second-order fallback; see
:doc:`capabilities` for the standing list.

.. _level-shifting:

Level shifting
--------------

A level shift adds a constant to the virtual block of the Fock matrix before each
diagonalisation, leaving the occupied block alone. The orbital energies it
diagonalises are the true ones for the occupied orbitals and the true ones plus
``level_shift`` for the virtuals, so the gap the next density is built through is
wider than the real one by exactly that amount.

That is the whole idea. A density built through an artificially wide gap moves
less from one iteration to the next, because the occupied-virtual mixing that
carries it is damped by the larger denominator. An SCF that was swinging past its
solution and back stops swinging.

.. code-block:: json

   "keywords": {
     "scf": {
       "level_shift": 0.5
     }
   }

The value is in Hartree and the default is ``0.0``, meaning off. It applies to
Hartree-Fock and Kohn-Sham alike, restricted and unrestricted, density-fitted and
conventional, and to the reference SCF underneath a CASSCF or ORMAS run.

A negative value is refused rather than clamped to zero. A negative shift lowers
the virtuals *toward* the occupied set, which narrows the gap and drives exactly
the oscillation the shift was asked to damp -- anyone typing one has the sign
backwards, and silently using zero would hide that.

Choosing a value
^^^^^^^^^^^^^^^^

Start at ``0.2`` and go up. Values between ``0.2`` and ``1.0`` Hartree cover most
cases. Larger shifts damp harder and cost more iterations; a shift of several
Hartree will usually converge something, slowly, and is worth trying once before
concluding that a system is hopeless.

There is no cost to a shift beyond iteration count -- two matrix products per
cycle against matrices already in hand -- so the tradeoff is entirely between
iterations saved far from the solution and iterations spent near it.

When it will not help
^^^^^^^^^^^^^^^^^^^^^

A level shift widens the gap the iteration *uses*. It cannot widen the gap the
molecule *has*, and it cannot conjure a bound state.

The case that looks most like a convergence failure and is not: an anion whose
extra electron is not bound at the chosen basis and level of theory. The highest
occupied orbital comes out at a positive energy, often with another orbital a
fraction of a millihartree away, and the SCF wanders between near-degenerate
states that are all equally unbound. No amount of shifting fixes this, because
there is nothing to converge *to* -- the density that minimises the energy wants
to put the electron at infinity, and a finite basis merely prevents it from
saying so.

The tell is the orbital energies, not the iteration count. If the HOMO is
positive and the HOMO-LUMO gap is a hundredth of an electronvolt, stop tuning the
solver. The usual fix is physical rather than numerical: put the system in a
continuum, where the solvent reaction field stabilises the excess charge and the
state becomes bound. See :doc:`continuum_solvation`.

What the shift does not change
------------------------------

Nothing that leaves the SCF carries the shift.

The shift is tapered off before convergence -- it is applied only while the
density is still moving by more than a hundred times the density tolerance, so
with the default ``1e-6`` it is gone by the time the RMS change reaches ``1e-4``
-- and convergence is not declared on a shifted iteration even if the tolerances
are met there. The orbitals and orbital energies handed back therefore belong to
the unshifted Fock operator.

This matters more than it looks. Those eigenvalues are read back downstream as
MP2 and coupled-cluster denominators, as the weights of the energy-weighted
density in every analytic gradient, and as the occupied energies and response
poles of an EFP fragment potential. A shift left in would move all of them by an
amount nothing downstream could recognise as a shift, and every one of those
numbers would be quietly wrong rather than visibly wrong.

The unit tests pin this: at a shift of 0.5 Hartree the RHF energy, the full
orbital energy spectrum, an MP2 correlation energy built on the result, and the
UHF energy all agree with the unshifted run to converged precision. The shift
changes the path, not the answer.

Where it is applied
-------------------

Two details of the implementation are worth knowing, because both are places a
level shift is commonly got wrong.

**After DIIS, not before.** The DIIS error vector is built from the unshifted
Fock matrix, and the extrapolation happens first; the shift is added to the
extrapolated matrix on its way to the diagonaliser. Shifting first would mean the
vectors DIIS stores are no longer a subspace of Fock matrices, so its
extrapolation is no longer extrapolating the thing it is meant to. It would still
converge to something, which is what makes that ordering expensive to find.

**Through a projector, not an orbital rotation.** Adding a constant to the
virtual block requires the virtual projector, which completeness supplies without
ever forming the virtual orbitals: since
:math:`C_o C_o^T + C_v C_v^T = S^{-1}`, the shift operator is
:math:`S C_v C_v^T S = S - \tfrac{1}{2} S D S` for a closed shell, and
:math:`S - S D_\sigma S` per spin for an open one. That is two matrix products
against the density already in hand.

Backend support
---------------

Level shifting is implemented on the CPU path. The cuEST GPU backend accepts the
keyword but does not currently apply it, so a GPU run with ``level_shift`` set
converges as though it were absent.

.. _accelerators:

Accelerators
------------

``keywords.scf.accelerator`` chooses which accelerator opens the SCF:
``diis`` (the default), ``ediis`` or ``adiis``. A name outside those three is
refused rather than ignored.

DIIS extrapolates from the error vectors of previous iterations, and it is very
good once those iterations are near enough to the answer to be informative. Far
from convergence they are not, and the extrapolation can be worse than the
iteration it replaces -- the classic symptom being an SCF that oscillates
between two densities rather than settling into either. EDIIS and ADIIS build
their step from the energies and densities instead, minimising a model that is
bounded below, so they cannot make that particular mistake.

They are also slower per iteration and less accurate near convergence, which is
why they do not simply replace DIIS.

The handover
^^^^^^^^^^^^

**Naming an accelerator asks for a different opening, not a different endgame.**
EDIIS and ADIIS run only while the commutator :math:`FDS - SDF` is above
``ACCEL_SWITCH``, currently ``1e-2``; below that the SCF hands over to DIIS for
the rest of the run. This is deliberate: the energy-based schemes earn their
keep exactly where DIIS is unreliable and lose to it everywhere else.

Two consequences are worth knowing before you conclude the setting is broken.

**Identical iteration counts are the expected result on an easy case.** If the
commutator starts below ``1e-2``, or drops below it in the first cycle, the
handover fires immediately and the run is a DIIS run. An N\ :sub:`2` CAS(6,6)
reference converges in 11 iterations under ``diis`` and 11 under ``ediis``; only
the trajectory differs, in the eighth decimal. **Do not use the iteration count
to check that the keyword took effect.** Two things do say so: an invalid name
is refused, and a non-default accelerator prints a line naming itself when the
SCF starts.

**EDIIS can stall just above the switch.** Because the threshold is a fixed
number and the stall is a property of the system, an SCF can settle at a
commutator slightly the wrong side of it and stay there, running the slower
algorithm precisely where it has stopped helping. Water / 6-31G / RHF, with
``accelerator: ediis``:

.. code-block:: text

    iter          energy        dD          commutator
       3   -75.983865252943   1.334E-03      1.975E-02
       4   -75.984014411771   2.563E-12      1.191E-02
       5   -75.984014411771   3.520E-04      1.191E-02
       6   -75.984037428504   1.059E-03      9.313E-03

It sits at ``1.191E-02`` against a switch of ``1e-2``, with the density
momentarily not moving at all, and only then falls through and grinds down.
The same molecule takes 8 iterations under DIIS and 12 under EDIIS. The symptom
a user sees is "EDIIS is slower"; the cause is that it is stuck just the wrong
side of a threshold. How long it lasts varies with the system -- runs of two and
of seven iterations have both been observed.

``ACCEL_SWITCH`` is a compile-time parameter. There is no keyword for it, so
when this happens the available answer is to use ``diis``, not to move the
threshold.

When to reach for one
^^^^^^^^^^^^^^^^^^^^^

Use ``ediis`` or ``adiis`` on an SCF that oscillates from a bad starting point --
a stretched bond, a transition-metal complex, a poor guess -- and DIIS on
anything that is merely slow. An SCF that crawls because the gap is small is not
helped by either; that is what :ref:`level-shifting` is for. The two are
independent and can be set together.

Backend support
^^^^^^^^^^^^^^^

Implemented on the CPU path, for Hartree-Fock, DFT, and the reference SCF of a
CASSCF or CASCI.

**The cuEST GPU backend refuses anything but** ``diis``. Its extrapolation is
device-resident and Pulay-only; no energy-based scheme is implemented there. A
deck naming ``ediis`` or ``adiis`` for a GPU run is refused by name rather than
answered with a DIIS run that says nothing about it, and a misspelled name is
refused by the same branch. Use the CPU backend if you want one of them.
