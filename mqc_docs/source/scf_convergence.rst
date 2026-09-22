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
badly and then hand back to it; there is level shifting; and there is a
second-order finish, which lets DIIS open and then converges the rest by
trust-region Newton on the orbital rotations. All three are described below.
There is no damping and no Fermi smearing; see :doc:`capabilities` for the
standing list.

They answer different questions, and the three-way split at the top of this page
is how to choose:

- **Oscillating** -- reach for :ref:`accelerators` (``ediis`` or ``adiis``), or
  for the second-order finish.
- **Crawling, with a small gap** -- reach for :ref:`level-shifting`. Neither
  accelerator helps.
- **Converged to the wrong thing** -- that is not a convergence problem at all,
  and it is the case :ref:`second-order-convergence` exists for.

.. _when-is-it-converged:

When is it converged
--------------------

Two things have to be true, and they are pyscf's two:

.. math::

   |E_n - E_{n-1}| < \texttt{tolerance}
   \quad\text{and}\quad
   \max_{\mu\nu} |(FDS - SDF)_{\mu\nu}| < \texttt{gradient\_tolerance}

The energy stopped moving, **and** the Fock matrix commutes with the density.

The iteration table shows both, under ``dE`` and ``diis``. It does **not** show a
density change: pyscf computes one, calls it ``norm_ddm`` and prints it without
testing it, and a number on the convergence table that is not part of the
convergence test reads as though it were. The commutator column is the one a
stalled SCF has to be read from -- see the EDIIS example below, where the density
change goes to 1e-11 while the commutator sits at 1.2e-2. The density change is
still computed, because the level shift tapers off it.

.. list-table::
   :header-rows: 1
   :widths: 28 18 54

   * - Key
     - Default
     - What it bounds
   * - ``keywords.scf.tolerance``
     - ``1e-9``
     - The change in energy between iterations
   * - ``keywords.scf.gradient_tolerance``
     - ``sqrt(tolerance)``
     - The commutator, as a max element. pyscf's ``conv_tol_grad``
   * - ``keywords.scf.density_tolerance``
     - ``1e-6``
     - **Nothing, for convergence.** Still sets where the level shift tapers off
   * - ``keywords.scf.convergence_metric``
     - ``standard``
     - Which measure decides. ``tolerance`` is read in its units

Choosing the measure
^^^^^^^^^^^^^^^^^^^^

``convergence_metric`` names the quantity ``tolerance`` bounds:

- ``standard`` -- the default, and what this program has always done: the
  energy **and** the commutator, the latter at ``sqrt(tolerance)`` or at
  ``gradient_tolerance`` when that is given.
- ``commutator`` -- the commutator alone, with ``tolerance`` read as its bound
  rather than an energy. ``diis`` and ``gradient`` are accepted spellings of
  the same thing: :math:`FDS - SDF` in the orthogonal basis is both what DIIS
  extrapolates against and the orbital gradient, pyscf's ``norm_gorb``.
- ``energy`` -- the change in energy alone. **The weakest of the three**, and
  see the warning below before choosing it.
- ``density`` -- the RMS change in the density matrix alone.

.. code-block:: json

   "keywords": { "scf": { "convergence_metric": "commutator", "tolerance": 1e-6 } }

Note the units change with the measure. A commutator of ``1e-6`` is about as
converged as an energy of ``1e-12``, because the energy's error falls as the
*square* of the commutator -- so a number carried across from a deck written
against ``standard`` will be far tighter than intended, and one written the
other way far looser.

.. warning::

   **The energy metric can stop early, and silently.** An accelerator that
   *interpolates* rather than extrapolates can hold ``dE`` still while the
   commutator stays large: EDIIS on water/6-31G sits at ``dE`` of 5.7e-13 with
   a commutator of 1.1e-2, and an energy-only test stops there -- 5.8e-5
   hartree from the answer. Selecting ``energy`` alongside ``ediis`` or
   ``adiis`` prints a warning for this reason. Prefer ``commutator`` if you
   want one measure, or ``standard`` if you want the safe default.

**Why the commutator has to be in it.** ``dE`` and ``dD`` say the iteration
stopped moving. They do not say it stopped at a stationary point, and an SCF can
hold both small while :math:`FDS - SDF` is nowhere near zero. Any scheme that
*interpolates* rather than extrapolates will do it, because a stalled
interpolation moves the density hardly at all. EDIIS on water/6-31G pins at a
commutator of 1.2e-2 for seven iterations while ``dE`` sits at 1e-12 and ``dD``
at 1e-11; without the commutator the SCF stops there, 5.9e-5 Hartree from the
answer. On water/def2-SVP the same failure has been seen at 0.11 Hartree. That
stall is a property of the scheme rather than an accident -- see
:ref:`accelerators`, which is where the schemes that cause it are described.

.. _gradient-tolerance:

Setting ``gradient_tolerance`` yourself
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Leave it alone if you want an energy. **Set it if you want a density.**

Near convergence the energy is variational, so its error falls as the *square*
of the commutator; the density's error falls only linearly. The default
``sqrt(tolerance)`` is the right scaling for the first and about three orders too
loose for the second. Anything computed from the converged density rather than
from its energy inherits that -- multipoles, atomic charges, a dipole, and any
finite-difference of one of those.

.. code-block:: json

   "keywords": { "scf": { "gradient_tolerance": 1e-8 } }

**It cannot be had by tightening** ``tolerance`` **instead.** Reaching a
commutator of 1e-8 through ``sqrt`` would need ``tolerance`` of 1e-16, which is
below what a molecular energy resolves in double precision, and the SCF then
never converges at all. The two demands are separate and take separate numbers,
which is why pyscf exposes ``conv_tol_grad`` rather than only deriving it. The
paths inside the program that consume a density -- the EFP fragment potentials,
the charge partitioning -- set it for themselves and say so where they do.

**There is a floor.** The commutator cannot be resolved below the scatter that
unordered OpenMP reduction merges put into it, and that scatter is a band whose
top moves with thread count: on AlH3/6-31G it plateaus at 9.8e-14 on one thread
and 3.1e-10 on sixteen. A threshold inside that band is sampling noise rather
than measuring convergence -- one at 1e-10 made the same binary converge at 1 and
8 threads and fail at 2 and 4. Larger systems on more threads should be expected
to floor higher.

Backends measure it differently
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CPU path takes the largest element of :math:`FDS - SDF`. The cuEST backend
takes a Frobenius norm over the whole matrix. **The same number does not mean the
same thing to both**, so ``gradient_tolerance`` is not portable between them at
present; on cuEST, leaving it unset keeps that backend's historical behaviour,
where the threshold was read from ``density_tolerance``. Harmonising the two
measures needs a GPU to verify on.

**One consequence for GPU runs, worth knowing before it surprises you.** The
energy leg of the test now defaults to ``1e-9`` on both backends, where it was
``1e-6``; the commutator leg on cuEST still falls back to ``density_tolerance``,
which has not moved. So a GPU run naming no tolerances at all converges *tighter
than it used to*, and takes more iterations to do it. ``gradient_tolerance``
cannot loosen that, because it moves only the second leg -- ``tolerance`` is the
knob. The asymmetry is deliberate: giving the GPU the CPU's derived default would
change what it converges to, on a measure that is not the CPU's, and that cannot
be verified without the hardware.

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

``keywords.scf.accelerator`` chooses the SCF's convergence scheme: ``diis``
(the default), ``ediis``, ``adiis`` or ``soscf``. A name outside those four is
refused rather than ignored.

The first three are the subject of this section. ``soscf`` is a different kind
of thing and has its own section below -- where ``ediis`` and ``adiis`` change
how the SCF *opens* and hand back to DIIS, ``soscf`` leaves the opening to DIIS
and changes how it *finishes*. See :ref:`second-order-convergence`.

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

.. _second-order-convergence:

Second-order SCF
----------------

Everything above is a first-order method. DIIS, EDIIS, ADIIS and the level shift
all work on a sequence of Fock matrices and densities; none of them knows
anything about the *curvature* of the energy. That is usually enough, and it is
very cheap.

A second-order SCF does know. It parametrises the orbitals as :math:`C \to C
\exp(\kappa)` with :math:`\kappa` antisymmetric, takes the gradient and the
Hessian of the energy with respect to those rotations, and steps by Newton's
method inside a trust region. Two things follow from that, and the second is
the one that makes it worth the cost.

**It converges in fewer iterations.** Newton's method converges quadratically
near a solution where DIIS does not, so an SCF that oscillates or stalls under
DIIS often simply stops doing so.

**It can tell a minimum from a saddle point, and DIIS cannot.** This is the real
reason. A first-order SCF stops wherever the gradient vanishes, and a saddle
point has a vanishing gradient -- so a converged DIIS run that has landed on one
reports an ordinary energy, ordinary orbital energies, and nothing at all to say
a lower solution exists a short rotation away. Everything built on top of it --
a correlation energy, a gradient, a frequency -- then describes the wrong
reference. The second-order SCF refuses to converge on negative curvature and
follows it out instead.

N\ :sub:`2` at 1.6 A in 6-31G is the case to keep in mind. DIIS converges from
the core, GWH and SAD guesses alike onto a saddle. The second-order run reaches
a genuine minimum between 0.020 and 0.242 hartree lower -- a difference no
tolerance would have caught, because both runs converged.

Asking for it
^^^^^^^^^^^^^

Two spellings, one feature. Either turns it on, and neither turns the other off:

.. code-block:: json

   "keywords": { "scf": { "accelerator": "soscf" } }

.. code-block:: json

   "keywords": { "scf": { "second_order": true, "soscf_start": 1e-2 } }

The accelerator spelling exists because this is a choice of convergence scheme,
made where ``diis`` and ``ediis`` are made. ``second-order`` and
``second_order`` are accepted as the accelerator name too, and case does not
matter.

.. list-table::
   :header-rows: 1
   :widths: 30 16 54

   * - Key
     - Default
     - What it does
   * - ``keywords.scf.accelerator``
     - ``diis``
     - ``soscf`` selects the second-order finish
   * - ``keywords.scf.second_order``
     - ``false``
     - The same request, as a flag
   * - ``keywords.scf.soscf_start``
     - ``1e-2``
     - The commutator :math:`\max|FDS - SDF|` at which DIIS hands over.
       Governs **both** spellings
   * - ``keywords.scf.stability``
     - ``false``
     - A separate feature. See below

DIIS always opens, whichever spelling you used
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Naming** ``soscf`` **does not start the SCF second order, and nothing can.**

A Newton step is a step on a quadratic model of the energy, and that model is
the energy only near the point it was built at. An initial guess is not near the
solution. A Newton step from one is long, arbitrary in direction and routinely
uphill, and a second-order SCF started from a guess diverges -- which is why
every one in practice is preceded by something else. Here that something else is
the ordinary DIIS SCF, run until :math:`\max|FDS - SDF|` falls below
``soscf_start``.

So ``soscf_start`` is not a setting the accelerator spelling lets you skip. It
governs both routes identically, and it is what makes the method robust rather
than a curiosity:

- **Raising it** hands Newton a worse starting point, and the backtracking
  search then spends Fock builds recovering from long steps.
- **Lowering it** spends DIIS iterations doing work the Newton steps would have
  done in fewer.

The iteration at which the handover happened is logged, and the second-order
phase prints its own table -- energy, energy change, commutator, gradient, trust
radius, curvature, and the Fock builds that row cost.

What it costs
^^^^^^^^^^^^^

**Count Fock builds, not iterations.** One second-order iteration costs one Fock
build per trial step, accepted or rejected, plus one per Hessian-vector product
inside the Newton solve -- up to twenty of those, though the residual test
usually stops well short. On seven well-behaved closed shells all converged to a
commutator of 1e-9, the second-order path took 5-6 iterations and 7-8 energy
builds against DIIS's 9-14 and 10-15, but spent 13-18 Hessian-vector products on
top. That is *fewer iterations and roughly twice the integral passes*.

On a well-behaved molecule, then, this is the wrong choice and DIIS is the right
one. It earns its cost in exactly two places: where DIIS oscillates or stalls,
and where the answer has to be a minimum rather than merely a stationary point.

:ref:`second-order-scf` in the input reference carries the step-by-step
mechanics -- semicanonicalisation, the Krylov solve, the trust-region
backtracking -- and the full measurement table.

What it will not do
^^^^^^^^^^^^^^^^^^^

Each of these is refused by name, before the SCF runs, rather than approximated
or silently demoted to DIIS.

- **An unrestricted reference.** The rotations parametrised here are the
  closed-shell ones; the open-shell space is larger and is not implemented.
- **A continuum solvent.** The orbital-rotation Hessian carries no response of
  the surface charges, so the step would be taken on the curvature of a
  different energy than the one being minimised. See
  :doc:`continuum_solvation`.
- **A Fock projector** -- frozen orbitals, as the AFO bond treatment uses. The
  rotations the projector forbids are not excluded from the Newton step's
  parameter space, so the step would break the constraint.
- **The GPU backend.** There is no orbital-rotation Hessian under cuEST. A GPU
  deck naming ``soscf``, or setting ``second_order``, is refused rather than
  answered with a DIIS run that says nothing about it.
- **The MCSCF reference SCF, the Fukui ions, and the EFP fragment potentials.**
  The first has its own second-order orbital optimiser; the second is open
  shell; the third is not implemented.

A restricted Kohn-Sham reference **is** supported -- the exchange-correlation
kernel enters the Hessian through the same operator the analytic Hessian's
coupled-perturbed solve uses.

One asymmetry is worth knowing rather than discovering: with a density-fitted
reference the Hessian is built from exact integrals, so it is the curvature of a
slightly different surface than the one being minimised. That degrades the
convergence *rate*, not the answer -- the gradient and the energy a step is
accepted on are the fitted ones, so the stationary point reached is the fitted
SCF's.

It is not the stability analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``keywords.scf.stability`` is a **separate feature**, and the two are easy to
conflate because they diagonalise the same matrix.

``stability`` runs *after* a converged SCF and answers a question: it finds the
lowest eigenvalue of the electronic Hessian and reports whether the solution is
a minimum. It changes nothing. The second-order SCF uses the same Hessian
*during* the iteration, to take steps with.

They also differ in what their curvature numbers are worth, and the difference
matters. The stability analysis converges an eigenvector, so its number is the
smallest curvature. The second-order SCF reads its curvature off a small Krylov
subspace, and a value from a subspace is always an **upper bound** on the true
smallest eigenvalue. So a negative one is proof of a saddle, and a positive one
is not proof of a minimum. The second-order SCF can therefore refuse to converge
at a saddle, and it cannot certify a minimum. ``stability`` is what certifies a
minimum.

Running both is the belt-and-braces combination, and they are independent
keywords:

.. code-block:: json

   "keywords": {
     "scf": {
       "accelerator": "soscf",
       "stability": true
     }
   }

See :ref:`scf-stability` for what the analysis reports and what its verdict does
and does not cover.

When it still will not converge
-------------------------------

In rough order of what to try, and the first step is the one most often skipped.

**1. Find out which problem you have.** Read the iteration table, not the exit
status. An energy swinging up and down is oscillation; an energy falling
steadily by ever-smaller amounts with a commutator that will not follow is a
small gap; a run that converges quickly and cleanly is not a convergence problem
at all and belongs at step 6.

**2. Start somewhere better.** A convergence problem is often a guess problem,
and it is far cheaper to fix there. ``keywords.scf.guess`` defaults to ``auto``,
which is ``sad`` on the CPU path. ``basis_set_projection`` converges a small
basis first and projects the density up a ladder, which is the answer when the
large basis is itself the difficulty -- diffuse functions especially -- or when
the SCF converges to the wrong state. It is not free, and it is the wrong tool
for an SCF that merely crawls. See :doc:`scf_guess`.

**3. Match the aid to the symptom.** ``ediis`` or ``adiis`` for oscillation,
``level_shift`` of 0.2 to 1.0 for a small gap. They are independent and can be
set together. Widening ``diis_size`` to 12-20 is the first thing to try on an
SCF that converges monotonically but slowly, where a level shift would only make
it slower.

**4. Turn off what might be hiding the problem.** ``incremental_fock: false``
forces a full Fock build every iteration. A run that then converges was being
broken by accumulated increments, which is worth knowing. ``diis: false`` is a
diagnostic in the same spirit rather than a setting.

**5. Go second order.** ``accelerator: "soscf"``. This is also the point at
which to suspect that the difficulty is real curvature rather than a bad path
through it.

**6. Ask whether there is anything to converge to.** Two cases, and no solver
setting reaches either.

*An unbound anion.* The extra electron is not bound at this basis and level of
theory: the HOMO comes out at a positive energy, often with another orbital a
fraction of a millihartree away, and the SCF wanders between near-degenerate
states that are all equally unbound. The tell is the orbital energies, not the
iteration count. The fix is physical -- put the system in a continuum, where the
reaction field stabilises the excess charge and the state becomes bound. See
:doc:`continuum_solvation`.

*A converged answer that is wrong.* The SCF stopped at a saddle point. Nothing
in the iteration table says so. Set ``stability: true`` and read the verdict; if
it is a saddle, ``accelerator: "soscf"`` will follow the negative curvature out
to the real minimum.

**7. Accept a non-converged result deliberately, if that is the right call.**
``allow_crap_scf: true`` keeps the last iterate rather than failing. In a
fragmented run this is often the only way to finish at all -- a handful of
fragments out of millions will not converge, and stopping on the first wastes
the other million. The fragments that failed are named in the output, so the run
can be followed up rather than trusted.
