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
     - yes
     - ``b97m-v``, ``wb97m-v``
   * - MP2, all-electron
     - yes
     - ``mp2`` with ``keywords.correlation.freeze_core: false``
   * - MP2, frozen core
     - yes
     - the default -- ``freeze_core`` is on unless the deck turns it off
   * - Double hybrid, all-electron
     - yes
     - ``b2plyp``, ``b2gp-plyp``, ``mpw2plyp`` with
       ``keywords.correlation.freeze_core: false``
   * - Double hybrid, frozen core
     - **no**
     - the default -- ``freeze_core`` is on unless the deck turns it off
   * - Spin-scaled or fitted MP2
     - **no**
     - ``scs-mp2``, ``sos-mp2``, ``ri-mp2``
   * - Coupled cluster
     - **no**
     - ``ccsd``, ``ccsd(t)``
   * - Unrestricted (UKS)
     - **no**
     - any open shell
   * - Density fitted
     - **no**
     - ``keywords.scf.density_fitting``

Anything not covered falls back to the semi-numerical path rather than
failing, so a deck that asks for a Hessian gets one either way. The difference
is how long it takes.

Note what the MP2 rows imply together: a plain ``mp2`` Hessian deck freezes
its core by default and now takes the analytic path with it -- the frozen
count flows into the correlation assembly, whose core rotations carry it.
A run that took the analytic path logs
``computing the analytic MP2 Hessian``.

The double-hybrid rows say the opposite of the MP2 ones about a frozen core,
and the asymmetry is real rather than an oversight: the frozen blocks the MP2
Hessian uses are built and correct, but they have been validated over a
Hartree-Fock reference only, never against a Kohn-Sham operator carrying its
kernel in the response. A frozen-core double-hybrid deck therefore falls back
rather than reusing them on trust. Its *gradient* refuses outright for the same
reason, which means such a deck does not reach the Hessian at all.

A double hybrid that took the analytic path logs ``computing the analytic
Hessian`` followed by ``computing the perturbative term's Hessian``. Seeing
only the first line on a ``b2plyp`` deck would mean the Hessian excluded what
the energy included -- which is what this used to do, returning the underlying
hybrid's Hessian with a Frobenius norm of 2.146683 against a true 2.120638 on
water at STO-3G.

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
   * - ``b97m-v``
     - 2.11e-05
     -
     -

against Hessian elements of order one.

``b97m-v`` sits apart, and the gap is understood rather than tolerated: mqc
holds the NLC grid fixed like every other grid, while PySCF's
``_get_vnlc_deriv1`` hard-codes the NLC grid response even under
``grid_response = False``. The difference between the two Hessians is exactly
that term, and it behaves like one -- 2.11e-05 with the NLC grid at level 1
(the default), 3.56e-06 at level 2 and 5.98e-07 at level 3, on both codes'
matched grids. A missing derivative would not move when the grid did; this
shrinks the way quadrature response must.

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

An MP2 Hessian adds a fourth piece on top of those three: the **correlation
block**, itself three groups -- the second-derivative integrals contracted
against the unperturbed relaxed densities, the orbital response through the
skeleton Lagrangian, and the 2n+1 density response, which is one perturbed
Z-vector solve per nuclear coordinate on top of the coupled-perturbed
equations the reference already solves. The second-derivative two-electron
integrals dominate the cost and are generated once, feeding the correlation
block and the reference skeleton from the same sweep. The assembled block is
checked element by element against pycc's analytic MP2 Hessian
(water/6-31G, symmetric and asymmetric geometries, agreement at ~5e-12
all-electron and ~6e-12 frozen-core) and
transitively against a seven-point finite difference of our own correlation
gradient; the symmetry and translational sum rules it must earn, and the
identity between its reference skeleton and the standalone Hartree-Fock
Hessian, are pinned in ``test_mqc_mp2_hessian_assembly``.

A frozen core changes the orbital response, not just the index ranges: the
core--active block of the orbital rotation takes the canonical Brillouin
divide with its coupling terms, the skeleton carriers gain the closed-form
pair-rotation augmentation, and the perturbed core--active rotation is a
Sylvester relation whose off-diagonal Fock-derivative couplings the naive
quotient rule drops -- an error pycc measured at ~7e-7 on this molecule and
basis. All three are in the assembly and gated against pycc's frozen-core
reference intermediates, so the default ``freeze_core`` deck is served
analytically.

A ``-V`` functional's non-local correlation contributes in the same three
places -- an explicit second derivative, a Fock-derivative term, and a kernel
in the response operator -- each on the NLC grid, each a pair sum rather than
a per-point expression. The response kernel is applied to every perturbation
in a batch at once, because its pair sweep costs the same whether it carries
one trial density or a dozen, and a per-perturbation application would
multiply the only expensive part by :math:`3N` on every iteration of the
solve.

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

The NLC grid of a ``-V`` functional is held fixed the same way, and here the
two codes part company: PySCF's ``_get_vnlc_deriv1`` includes the NLC grid
response unconditionally, with a comment warning that omitting it can shift
the coupled-perturbed solution enough for ~1e-3 Hessian error on some systems
(it names H2O2 and H2CO). mqc omits it for the same reason it omits the
semilocal one -- a fixed-grid second derivative is exactly differentiable
against a fixed-grid gradient, which is what every VV10 piece here was
validated by. On water the whole disagreement is 2e-5 and dies with grid
refinement, as the table above records; on a system where it does not, refine
``keywords.dft.nlc_grid_level``.

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

**Spin-scaled MP2.** The same reason its gradient is refused: a scaled MP2 is
a different energy, and its derivatives are not the unscaled ones rescaled --
the amplitudes enter the relaxed density and the response equations unscaled.

**RI-MP2.** The fitted correlation Hessian does not exist yet, and the fitted
reference is excluded for the same reason as everywhere else in this table.

**Double hybrids that are not GGA-based, all-electron and unscreened.** Four
cases fall back, and they divide into two kinds. A meta-GGA or a
range-separated double hybrid is missing something nameable: the kernel's third
derivative in ``tau`` for the first, and an operator at a screened omega for the
second, where the correlation block currently scales exchange by a single
``exx_fraction``. Both of those *object* -- ``xc_kernel2_apply`` and
``xc_kernel_deriv`` refuse a meta-GGA rather than returning a partial kernel.

A VV10 double hybrid is the case worth naming separately, because it would not
object. ``ks_hessian`` carries the non-local term, so the reference block is
complete; it is the perturbative term's kernel derivatives that do not, and
they would silently return the semilocal part alone. The guard in the bridge
stands in for a refusal the kernel routines do not make.

A frozen core is the fourth, and is covered above.

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
and 2e-5 at cc-pVDZ, while Hartree-Fock stays at 1e-8 on all three.

**That residual is the reference's, not the quadrature's.** psi4's analytic
restricted Hartree-Fock Hessian is exactly symmetric and sits 2.3e-9 from ours
on water/cc-pVDZ, while PySCF's sits 4.0e-8 from both -- so on the one case a
third code can arbitrate, the disagreement is PySCF's. Its size is visible
without a third code: a same-atom 3x3 block must be symmetric, PySCF mirrors
its *inter*-atom blocks by construction (those agree to 2e-16) and computes the
same-atom ones directly, and those come back asymmetric by 9.3e-9 for
Hartree-Fock and 1.06e-5 for PBE/cc-pVDZ. Ours are symmetric to 2e-10, and the
disagreement tracks that asymmetry case by case.
``pyscf.hessian.rhf.solve_mo1`` also calls ``cphf.solve`` without passing
``tol``, so its coupled-perturbed threshold is the library default and no
attribute reaches it -- ``max_cycle``, ``conv_tol_cpscf`` and ``mol.precision``
all leave the result bit-identical.

The quadrature reading was the first one here too, and what rules it out is
that the disagreement does not move with the grid: 1.7e-5, 8.5e-6 and 1.6e-5
for PBE/cc-pVDZ at levels 3, 5 and 7, while the energy holds at 1e-12
throughout. What *is* quadrature is the analytic Hessian itself, which is not
grid-converged until level 7 -- level 3 is 5.7e-3 out on a meta-GGA -- but both
codes carry that equally and it cancels in the comparison.

In frequency terms the worst of it is 0.16 cm-1 on a mode at 3894 cm-1, below
anything the underlying model resolves.

There is no third opinion for the density-functional cases: psi4 has no
analytic Kohn-Sham Hessian and silently differences gradients instead, which
puts the grid response back in and lands 7.1e-4 away -- a useful confirmation
that the omission is real and shared, and no use at all as an arbiter. So those
cases on cc-pVDZ are held to 1e-4, still three orders inside what these bounds
exist to catch: the three assembly errors found while this was written were
0.87, 0.26 and 0.12.

Validated by
============

``validation/inputs/cpu/mqc/hessian/`` -- every element of every Hessian
against ``pyscf.hessian.rhf`` or ``pyscf.hessian.rks`` with
``grid_response = False``, carried in ``validation_tests_cpu.json`` and
generated by ``tools/cpu_validation``. A manifest case per functional rather
than the standalone ``check_hessian.py``, so a regression is caught by the same
sweep as everything else instead of by remembering to run a script.

The meta-GGAs on d functions are validated there rather than by differencing,
and that is deliberate. The tau channel is the most grid-sensitive of the
three, so on a level-3 grid a central difference of the exchange-correlation
gradient does not converge at all: swept from ``h = 1e-2`` down to ``6.25e-4``
it returns 0.589, 0.792, 0.740, 0.768, 0.748 for TPSS and 4.04, -0.30, 0.67,
0.53, 0.63 for M06-L. A second difference of the *energy* scatters the same
way, so it is the function that is not smooth and not either derivative. At
levels 7 and 9 the sweep converges, and converges to the analytic Hessian.
