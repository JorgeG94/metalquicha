==================
Analytic Hessians
==================

A Hessian can be had two ways. The semi-numerical path differences analytic
gradients: ``6N`` of them for ``N`` atoms, each a converged SCF. The analytic
path differentiates the energy expression twice and runs one SCF. On glycine
trimer -- 24 atoms, 139 basis functions, PBE/6-31G, grid 3, twenty threads --
the semi-numerical Hessian takes **1962 seconds** and 144 gradient
evaluations. That is the number the analytic path exists to beat, and it is
why this is worth having even though the semi-numerical one is correct.

What is covered
===============

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Family
     - Analytic?
     - Examples
   * - Hartree-Fock
     - yes
     - ``hf``
   * - LDA
     - yes
     - ``svwn``, ``lda``
   * - GGA
     - yes
     - ``pbe``, ``blyp``
   * - Hybrid GGA
     - yes
     - ``b3lyp``, ``pbe0``
   * - meta-GGA
     - yes
     - ``tpss``, ``m06-l``
   * - Range-separated hybrid
     - yes
     - ``wb97x``, ``cam-b3lyp``
   * - VV10 non-local correlation
     - **no**
     - ``b97m-v``, ``wb97m-v``
   * - Unrestricted (UKS)
     - **no**
     - any open shell
   * - Density fitted
     - **no**
     - ``keywords.scf.density_fitting``

Anything not covered falls back to the semi-numerical path rather than
failing, so a deck that asks for a Hessian gets one either way. The difference
is how long it takes.

Every covered functional is checked against ``pyscf.hessian.rks`` with
``grid_response = False``, water at STO-3G, on every element of the Hessian:

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 25

   * - Functional
     - Worst element
     - Functional
     - Worst element
   * - ``lda_x``
     - 1.83e-08
     - ``pbe0``
     - 1.74e-08
   * - ``pbe``
     - 1.77e-08
     - ``tpss``
     - 1.93e-08
   * - ``b3lyp``
     - 1.46e-08
     - ``m06-l``
     - 6.56e-08
   * - ``cam-b3lyp``
     - 1.71e-08
     - ``wb97x``
     - 1.27e-08

against Hessian elements of order one.

How it is assembled
===================

Three pieces, and they only mean anything added together:

**The nuclear repulsion**, which is arithmetic.

**The explicit term** -- every second derivative of the integrals and of the
exchange-correlation energy, contracted against the density the SCF converged
to and held fixed. Exact exchange is scaled by the functional's mixing
fraction, so a pure functional has none and a hybrid has its own.

**The response** -- what the Hessian gains from letting the density relax. This
is the coupled-perturbed problem, and for a Kohn-Sham reference three things in
it are easy to get wrong and were each wrong once here:

* the operator needs the exchange-correlation kernel, not only Coulomb and
  scaled exchange;
* the derivative Fock matrix needs the exchange-correlation potential's own
  nuclear derivative, which the integral routines cannot supply;
* the energy-weighted density must be built from the **Kohn-Sham** Fock
  matrix, not a Hartree-Fock one.

Each produces a Hessian that is symmetric, translationally invariant and wrong,
which is why they are named here.

A range-separated functional splits its exchange over an ``erf`` kernel, so
every one of those places needs a *second* exchange pass at the screened
:math:`\omega` on top of the full-range one. No new integrals are involved:
libfint switches kernels through an environment slot rather than a separate
entry point, so an attenuated pass is the same loop over a modified copy of the
environment. That last word is the trap. Two routines here built the attenuated
copy and then handed the integral the *unattenuated* one, so the long-range
pass quietly returned full-range integrals scaled by the long-range
coefficient. The result stayed symmetric and translationally invariant --
full-range integrals are perfectly good integrals of the wrong operator -- and
sat 2.1 from PySCF for wB97X. Neither bug was findable from the Hessian's
shape; both fell out immediately once it was differenced against a gradient
that had been validated first.

.. _hessian-grid-response:

The grid does not respond
=========================

The quadrature is atom-centred, so a rigorous second derivative carries the
grid points moving with their owner and the partition weights being reset by
every nucleus. Both are omitted, which is what ``pyscf.hessian.rks`` does by
default -- its ``grid_response`` attribute is ``False`` -- so the two are
comparable term for term rather than approximately.

The omission is visible in translational invariance, which a Hessian should
satisfy exactly. On water at STO-3G, grid level 3, PySCF's own Hessian violates
it by 5.4e-4 for LDA, 5.7e-4 for PBE, 1.2e-3 for TPSS and 6.7e-3 for M06-L.
Ours tracks those. It shrinks as the grid is refined -- 5.4e-4 at level 3 and
2.6e-6 at level 9 for LDA -- because an exact quadrature does not depend on
where its points are.

For frequencies this is small. If it matters for what you are doing, there are
two ways to get it back. Refining the grid shrinks it, as the table above
shows. Or take the semi-numerical path, which does **not** carry this omission:
it differences the physical ``xc_gradient``, and that gradient includes both
the point-motion and the partition-weight terms explicitly. That is precisely
why ``xc_gradient_fixed_grid`` had to be written to test this one -- the two
gradients disagree by exactly the omitted term, so differencing the physical
one is not the instrument for checking a fixed-grid Hessian, and is the better
instrument for computing a responsive one.

.. _hessian-missing:

What is missing, and why
========================

Each of these is refused rather than approximated. A Hessian that silently
drops a term is worse than one you cannot have: the frequencies come out
plausible.

**VV10 non-local correlation** (``b97m-v``, ``wb97m-v``). Not implemented. The
*gradient* is -- and it turned out not to need new machinery, because the pair
sum is spent producing ``vrho`` and ``vsigma`` and what follows is the ordinary
GGA contraction. The second derivative is not the same story: it is a
double-grid object. PySCF implements it in three parts behind
dedicated kernels. On water the missing term is worth about 1.4e-3 in the
Hessian -- small, but the functional exists for dispersion-bound systems and
nothing says it stays small there.

**Unrestricted references.** The response machinery is restricted-only. An
open-shell Kohn-Sham Hessian needs both spin channels through every term and a
spin-resolved kernel. Mechanical rather than deep, and not done.

**Density fitting.** The Hessian has no fitted path at all, so everything is
conventional four-index. This one is blocked on integrals rather than algebra:
only the first derivatives ``int3c2e_ip1`` and ``int3c2e_ip2`` exist, and a
fitted Hessian wants ``int3c2e_ipip1``, ``ipvip1``, ``ip1ip2`` and
``int2c2e_ipip1``, none of which are among the procedures the integral library
exposes. Doing it conventionally first is the right order anyway: fitting
trades accuracy for speed, and without an exact Hessian to difference against
there is no way to say whether a fitted one is converged or merely fast.

**PCM.** A continuum's contribution to a second derivative is not implemented,
and its interaction with the response terms has not been worked out.

How a deck reaches it
=====================

It does not ask. ``driver: Hessian`` requests the analytic Hessian, the backend
takes it where the table above says yes, and everything else falls back to
central differences of analytic gradients without saying so as an error --
declining is not failing, because a Hessian has a second way to get one that is
correct and merely slower. A run that took the analytic path logs
``computing the analytic Hessian`` at verbose level, which is the only way to
tell from the outside.

There is no keyword to force either path. A deck that wants the fallback can
have it by asking for something the analytic path declines -- density fitting,
say -- but that is a side effect rather than a control, and if a real need for
one appears it should be a keyword rather than a trick.

Accuracy on larger bases
========================

The 1e-8 figures above are water at STO-3G, and they do not hold everywhere.
Against PySCF, a functional agrees to 1e-9 relative at STO-3G, 8e-6 at 6-31G
and 2e-5 at cc-pVDZ, while Hartree-Fock stays at 1e-8 on all three. The
residual is the two codes' quadrature rather than the derivatives: it moves
with the grid level and vanishes at STO-3G, where the density is compact enough
that where the points sit stops mattering.

In frequency terms the worst of it is 0.16 cm-1 on a mode at 3894 cm-1. That is
below anything the underlying model resolves, but it is not the agreement the
STO-3G numbers suggest, and it has not been chased further. ``check_hessian.py``
in the validation suite carries both figures, one tolerance per case, rather
than loosening every case to the worst one.
