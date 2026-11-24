==================
Analytic Hessians
==================

A Hessian can be had three ways. The semi-numerical path differences analytic
gradients: ``6N`` of them for ``N`` atoms, each a converged SCF. The analytic
path differentiates the energy expression twice and runs one SCF. The numerical
path differences energies alone, ``1 + 18 N**2`` of them, and needs neither (see
`How a deck reaches it`_). On glycine
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
   * - RI-MP2, all-electron or frozen core
     - yes
     - ``ri-mp2`` over an exact reference (``model.aux_basis`` named,
       ``keywords.scf.density_fitting`` off)
   * - Double hybrid, named auxiliary basis
     - yes
     - a double hybrid above with ``model.aux_basis`` set; all-electron
   * - Spin-scaled MP2
     - **no**
     - ``scs-mp2``, ``sos-mp2``
   * - Coupled cluster
     - **no**
     - ``ccsd``, ``ccsd(t)``
   * - Unrestricted (UKS)
     - **no**
     - any open shell
   * - Density-fitted Hartree-Fock and Kohn-Sham
     - yes
     - ``hf`` or ``dft`` with ``keywords.scf.density_fitting``, provided the
       auxiliary metric inverts without dropping a mode. LDA, GGA, meta-GGA,
       global hybrids and VV10
   * - Density-fitted range-separated hybrid
     - **no**
     - ``wb97x``, ``cam-b3lyp`` and the like with
       ``keywords.scf.density_fitting``
   * - RI-MP2 over a density-fitted reference
     - yes
     - ``ri-mp2`` with ``keywords.scf.density_fitting``, frozen core allowed
   * - Double hybrid with both halves fitted
     - yes
     - ``b2plyp`` and the like with ``model.aux_basis`` and
       ``keywords.scf.density_fitting``, under the conditions the unfitted
       double hybrid has
   * - Conventional MP2 over a density-fitted reference
     - **no**
     - ``mp2`` with ``keywords.scf.density_fitting``

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

An ``ri-mp2`` deck logs ``computing the analytic RI-MP2 Hessian`` instead. Its
correlation block (``mqc_czt_ri_mp2_hessian``) differentiates the fitted energy
the deck reports: the reference keeps exact integrals and only the correlation
functional contracts the fitted tensor, with no four-index array anywhere. A
double hybrid whose deck names an auxiliary basis takes the same block, with
the functional's kernel in every response. It needs an auxiliary metric that
drops no mode; one that does is declined at run time, with a warning, and the
Hessian is taken semi-numerically. The analytic gradient of these decks solves
its Z-vector with a fitted operator, which sits about 1e-7 from the
exact-reference derivative, so the Hessian and a difference of that gradient
agree to that and not further.

The :math:`3N` Z-vectors of the density response are solved together: each
conjugate-gradient iteration is one Fock-type build for a chunk of up to
``keywords.hessian.response_batch`` of them, the machinery the reference's own
coupled-perturbed solve uses, instead of one build per right-hand side. Every
right-hand side is converged on its own, to a 2-norm residual no larger than
``response_tolerance`` times its norm. On water trimer in cc-pVDZ this took that stage
from 45 s to 9 s on four threads, and the run from 64 s to 30 s.

The fitted block adds up what it will hold at once before it allocates -- a
few stacks of :math:`n_{mo}^2 \times 3N`, the three-index integrals, three
:math:`o^2 v^2` amplitude arrays and one atom's batch of derivative integrals --
and stops with an error stating the bytes needed if that exceeds the memory
budget (``system.memory_gb``, else four fifths of what the machine reports
available). It does not fall back by itself: ``"semi-numerical"`` is the
answer for a problem that does not fit.

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

**Orbitals removed for linear dependence.** The response equations take the
orbitals to span the basis. Where ``linear_dependence_threshold`` projected
functions out they do not, and every analytic Hessian here then returns a
converged, symmetric matrix that is wrong: 0.13 hartree/bohr\ :sup:`2` on HCN
in cc-pVDZ with four functions removed. Such a run is declined and its
gradient differenced instead.

**PCM.** A continuum's contribution to a second derivative is not implemented,
and its interaction with the response terms has not been worked out.

**Spin-scaled MP2.** The same reason its gradient is refused: a scaled MP2 is
a different energy, and its derivatives are not the unscaled ones rescaled --
the amplitudes enter the relaxed density and the response equations unscaled.

**A fitted reference.** A fitted energy is not the exact one with fitted
integrals put in, and its Hessian is the derivative of the expression the SCF
minimised. With :math:`V = (\mu\nu|P)`, :math:`M = (P|Q)` and
:math:`c = M^{-1}V`, every perturbation contributes
:math:`U^x = V^x - M^x c`, and the fixed-density second derivative is the
second-derivative integrals against two fixed densities plus a product of
first derivatives, :math:`\sum \Gamma\, U^x M^{-1} U^y`. That product is
where the fit's own response to a moved auxiliary basis lives. A Hessian that
carries the auxiliary basis's response only to first order is still symmetric
and still translationally invariant, and on water is 1e-4 of an element away
from the derivative of the fitted gradient.

Hartree-Fock and Kohn-Sham take this route. Nothing in it is built per perturbation at the
size of the three-index integrals: the exchange-type contractions go through
the occupied orbitals, each atom touches only the derivative integrals whose
differentiated function it owns, and the product term is accumulated a batch
of auxiliary functions at a time. The response solver applies its operator in
the molecular basis, :math:`N^2 N_{occ} N_{aux}` per trial vector, and the
relaxed mean field is one more application of it. On a water tetramer in
cc-pVDZ with cc-pVTZ-JKFIT, on four threads, the whole calculation takes 3.4 s
against 13.8 s for the exact analytic Hessian and about a minute for differencing the
fitted gradient. The integrals are still held whole, about four arrays of
:math:`N^2 N_{aux}` at the peak.

It declines when the auxiliary metric drops a mode, because
:math:`-M^{-1} M^x M^{-1}` is the derivative of the inverse only when none is
dropped.

Under a functional the fitted terms carry its exact-exchange fraction -- none
for a pure functional, the mixing fraction for a hybrid -- and the
exchange-correlation second derivatives, potential derivative and kernel go
exactly where the unfitted Kohn-Sham Hessian puts them; none of the three
knows how Coulomb and exchange were built. The kernel is the one part of the
response operator that needs the trial densities on the grid, so they are
built for it alone. That grid work is most of a Kohn-Sham Hessian either way:
on the same water tetramer a fitted PBE Hessian takes 25 s against 36 s
unfitted, of which the kernel inside the response solve is 16 s and 22 s.

The reference is PySCF's fitted Kohn-Sham Hessian with
``auxbasis_response = 2`` and agrees to 1e-8 or better, on bases with d
functions too -- tighter than the unfitted comparison, where PySCF's own
same-atom blocks limit it to 1e-5. VV10 has no such reference, because PySCF
carries the non-local grid's response and this code holds that grid fixed; the
fitted Hessian sits the same 4.7e-4 from PySCF's as the unfitted one does, and
the two offsets agree to 2e-8.

**Range-separated hybrids over a fitted reference.** Declined. The long-range
exchange has its own attenuated fit in the SCF and no derivative of it here.

**RI-MP2 over a fitted reference.** The fitted correlation block uses the
reference's operator in five places: the operator itself on a general density,
its first derivative at fixed orbitals, the two coupled-perturbed solves, and
the second derivative of the relaxed correlation density against the
reference's. Over an exact reference each is a pass over four-index
integrals. Over a fitted one each comes from the fit instead -- the same
three-index machinery as the fitted Hartree-Fock Hessian, generalised from the
reference's own density to any symmetric one -- and the block is added to the
fitted Hartree-Fock Hessian rather than the exact one. One auxiliary basis
serves both fits, because ``model.aux_basis`` is the only place a deck names
one.

This is the case where the analytic Hessian and the differenced gradient are
derivatives of exactly the same energy, with no operator standing in for
another: the gradient's Z-vector is solved against the fitted reference's own
operator. A four-point stencil of that gradient agrees to 1e-8. Over an exact
reference the same comparison stops near 1e-5, because there the gradient's
fitted Z-vector operator is an approximation to the exact one.

It is also the faster of the two. On a cyclic water trimer in cc-pVDZ with
cc-pVDZ-RIFIT and a frozen core, on four threads, the whole calculation takes
5.2 s over a fitted reference against 29 s over an exact one. Every stage that
called a four-index Fock build -- the two coupled-perturbed solves, the Fock
responses of the rotations, and the operator derivatives inside the perturbed
ladder -- drops by a factor of four to twenty, and what is left is the
correlation ladder itself.

**A double hybrid over a fitted reference.** The same five branches with the
functional's exchange fraction, and the kernel terms where they were: the
kernel's share of the fold, its skeleton derivative and the potential's second
derivative are all evaluated on the grid and none of them reads a two-electron
integral. The coupled-perturbed solves take the fitted operator and the kernel
together, as the fitted Kohn-Sham Hessian's do. The reference Hessian under
the block is the fitted Kohn-Sham one. A seven-point stencil of the fitted
gradient on a fixed grid agrees with the block to 8e-10.

**The nuclear responses are solved once.** The reference's Hessian and the
fitted correlation block need the same first-order orbitals, from the same
equations at the same tolerance, so the block takes the reference's and skips
its own coupled-perturbed solve; the log says ``nuclear responses taken from
the reference's Hessian``. That holds for RI-MP2 and for a fitted double
hybrid, over an exact or a fitted reference. The conventional MP2 block
solves at a tighter tolerance than the deck's and keeps its own.

**A double hybrid's kernel terms take every perturbation in one pass over the
grid**, and are threaded over grid blocks. On 1,3,5-trifluorobenzene, B2PLYP
in 6-31G on four threads, the perturbed ladder's first stage went from 195 s
to 47 s and the whole Hessian from 343 s to 151 s, the Hessian itself
unchanged to 1e-11. The conventional (not fitted) double-hybrid block still
takes one perturbation at a time.

**Conventional MP2 over a fitted reference.** Still excluded. The conventional
block contracts four-index integrals throughout.

**Double hybrids that are not GGA-based, all-electron and unscreened.** Read
this entry against the table above rather than instead of it: a meta-GGA, a
range-separated hybrid and a VV10 functional all take the analytic path as
*functionals*, and only their **double hybrids** fall back. The line between
the two is the order of the functional derivative each needs, not the method.
A Kohn-Sham Hessian needs the second, :math:`f_{xc}`, which
``xc_kernel_apply`` provides at every rung and ``vv10_kernel_apply`` extends to
the non-local term. A double hybrid's perturbative block needs the *third*,
:math:`g_{xc}`, plus that kernel's own nuclear derivative -- and those two
routines, ``xc_kernel2_apply`` and ``xc_kernel_deriv``, stop at the GGA rung.

Four cases fall back, and they divide into two kinds. A meta-GGA or a
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

``keywords.hessian.method`` chooses among three routes. Its default,
``"analytic"``, is the behaviour described above. There is no way to force the
analytic route where the table declines a method.

The analytic route is a single-rank route. An unfragmented Hessian launched on
several MPI ranks distributes its displacements across them instead (gradients,
or energies for ``"numerical"``), and the log says so. A fragmented run on
several ranks does the same, fanning each fragment or PIE term out into its
gradient displacements unless ``method`` is ``"numerical"``. Every one of those
displacements is taken at ``keywords.hessian.displacement``, the same step the
single-rank and unfragmented runs use, so a run gives the same Hessian on one
rank as on many.

``"semi-numerical"`` forces the fallback. The backend skips every analytic
Hessian, Hartree-Fock, Kohn-Sham, double hybrid and MP2 alike, and takes
central differences of analytic gradients, :math:`6N + 1` of them for
:math:`N` atoms; the run logs ``analytic Hessian skipped`` once, where it would
otherwise have computed the analytic one. It exists chiefly for the MP2
Hessian, whose memory grows as :math:`n_{mo}^4 \times 3N` and can be too large
where the :math:`6N` gradients are not. It needs the method's analytic
gradient, so a method with none (coupled cluster, a continuum solvent on the
CPU backend, a spin-scaled MP2) is refused by the gradient's own message, which
then points at ``"numerical"``.

``"numerical"`` takes no gradient and no analytic Hessian at all. The diagonal
is :math:`(E(+h_i) - 2E(0) + E(-h_i))/h^2` and each off-diagonal element the
four-point difference

.. math::

   H_{ij} = \frac{E(+h_i,+h_j) - E(+h_i,-h_j) - E(-h_i,+h_j) + E(-h_i,-h_j)}
                 {4h^2}

with :math:`h` the ``displacement``. Both are accurate to :math:`O(h^2)`, and
the matrix is symmetric by construction. That is :math:`1 + 2(3N) +
4 \cdot \tfrac{3N(3N-1)}{2} = 1 + 18N^2` energy evaluations: 163 for water,
289 for four atoms, against :math:`6N + 1` gradients for the semi-numerical
route. Because the second difference divides the energies by :math:`h^2`, their
noise is amplified by :math:`10^4` at the default step, so converge the SCF
tightly. On water at STO-3G, with ``keywords.scf.tolerance`` at ``1e-12`` and
the default step, the frequencies from this route agreed with the analytic ones
to 0.10 cm\ :sup:`-1` and the semi-numerical route's to 0.06 cm\ :sup:`-1`;
looser convergence costs accuracy here first. The same energies give a
central-difference gradient for nothing, which is reported beside the Hessian.
It is the only route for a method whose gradient does not exist, and on several
MPI ranks the energy evaluations are spread across them like the
semi-numerical route's gradients. A fragmented run hands each fragment (or PIE
term) to one rank, which differences it serially, instead of flattening the
displacements into the queue.

IR intensities need the dipole derivative. The analytic route supplies them for
Hartree-Fock and Kohn-Sham, and for a double hybrid through its relaxed
density. The semi-numerical and numerical routes difference the reference
density's dipole at the single displacements, which is right for Hartree-Fock,
Kohn-Sham and xTB; for a double hybrid it is not the relaxed density's, so they
leave the dipole derivatives out rather than report the wrong ones. An MP2 or
coupled cluster Hessian has none on any route.

What a deck can tune
====================

The coupled-perturbed solve behind the analytic Hessian takes three settings,
all under ``keywords.hessian``, and ``method`` chooses the route:

.. code-block:: json

   "keywords": {
     "hessian": {
       "response_batch": 12,
       "response_tolerance": 1e-9,
       "response_max_iter": 50,
       "method": "analytic"
     }
   }

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Key
     - Default
     - What it does
   * - ``response_batch``
     - 12
     - Perturbation densities carried through one pass over the integrals.
       The response equations are ``n_ov`` applications of one operator, and
       batching lets a single integral pass serve that many right-hand sides
       instead of one. Raising it trades memory for passes: every density in
       a batch is resident at once, so the ceiling is what fits rather than
       what is fastest in principle.
   * - ``response_tolerance``
     - 1e-9
     - Convergence of the coupled-perturbed solve.
   * - ``response_max_iter``
     - 50
     - Iteration cap for that solve. A solve that hits it is reported rather
       than quietly returning its last iterate.

Two things worth knowing before reaching for these.

``keywords.efp`` has its own ``response_batch``, a different key on a different
group driving the EFP polarization response. Setting one and watching the other
is the easy mistake.

And ``response_batch`` is a **performance dial, not an accuracy one** -- with
one exception worth naming. Changing it reorders floating-point accumulation,
which moves a double-hybrid derivative by around 5e-10. That is the noise floor
those quantities already carry rather than anything the batch introduces, so
varying this is a way to *measure* that scatter when chasing a small
discrepancy, and not a way to remove it.

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
