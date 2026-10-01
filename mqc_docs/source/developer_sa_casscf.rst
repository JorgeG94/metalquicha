.. _developer_sa_casscf:

===============================================
SA-CASSCF Gradients: Conventions and Equations
===============================================

This page fixes the equations for state-averaged CASSCF energies and
analytic gradients in the conventions this code already uses, and names the
existing routine each term reuses and the routines that are new. The
reference implementation checked against is PySCF 2.14's
``pyscf/grad/sacasscf.py`` and ``pyscf/mcscf/newton_casscf.py``.

Conventions this code already has
==================================

One- and two-particle density matrices
---------------------------------------

``active_space_rdms`` (``src/methods/ci/mqc_rdm.f90``) builds

.. math::

   D_{pq} = \langle \Psi | E_{pq} | \Psi \rangle, \qquad
   d_{pqrs} = \langle \Psi | E_{pq} E_{rs} - \delta_{qr} E_{ps} | \Psi \rangle

with :math:`E_{pq}` summed over both spins -- **spin-traced**, not spin-orbital
-- and in **chemist ordering**: :math:`d_{pqrs}` contracts against
:math:`(pq|rs)`, not the physicist's :math:`\langle pq|rs\rangle`. The module
docstring states the convention is PySCF's ``make_rdm12``, which is also
Helgaker's and GAMESS's, so nothing here is a private choice.

There is **no extra factor of 1/2 folded into** :math:`d_{pqrs}` itself; the
1/2 appears only in the energy contraction. ``rdm_energy`` (same file) rebuilds
the active-space energy as

.. math::

   E = \sum_{pq} h_{pq} D_{pq} + \frac12 \sum_{pqrs} (pq|rs)\, d_{pqrs}

and is used as a self-check against the CI eigenvalue, which catches exactly
the two failure modes the module docstring calls out: a transposed index or a
missing factor of two, either of which still gives a plausible-looking energy.

``run_czt_casci`` (``mqc_czt_casci.f90``) calls ``active_space_rdms`` **once**,
on ``result%ci_vector`` -- the lowest root only -- regardless of how many roots
``n_roots`` asked the Davidson to converge. ``davidson_lowest`` returns every
root's vector in ``davidson_result_t%vectors(n_alpha, n_beta, n_roots)`` and
every root's energy in ``%values(n_roots)``, both already computed; only the
1- and 2-RDM of root 1 are ever built from them.  ``casci_result_t%energies``
carries every root's total energy (core + active), ascending, but
``casci_result_t%dm1``/``%dm2`` belong to root 1 alone. State averaging needs
every root's density, so building the others is new work (a loop over
``result%vectors(:, :, i)`` through ``active_space_rdms``, not a new routine),
covered under the file-level plan below.

Orbital rotation parameter :math:`\kappa`
------------------------------------------

Orbitals move as :math:`C \to C \exp(\kappa)` with :math:`\kappa` antisymmetric
(``mqc_czt_mcscf.f90`` module header). The sign is fixed by finite differences
in ``test_mqc_mcscf.f90``, not asserted -- both :math:`C\exp(\kappa)` and
:math:`C\exp(-\kappa)` appear in the literature as "the" parametrization, and
the two give gradients that are the same expression with opposite sign.  The
resulting orbital gradient is

.. math::

   g_{pq} = \frac{\partial E}{\partial \kappa_{pq}} = 2\,(F_{qp} - F_{pq})

built by ``orbital_gradient`` from the generalised Fock matrix :math:`F`
(``mcscf_fock_t%general``, built by ``generalized_fock``). :math:`g_{pq}` is
antisymmetric by construction, which is what makes it the gradient of an
antisymmetric parametrization at all.

**Non-redundant blocks.** ``is_redundant`` (same file) marks a rotation
:math:`(p,q)` redundant when :math:`p` and :math:`q` fall in the same orbital
class (inactive/active/virtual) -- for a *complete* active space; an
occupation-restricted (ORMAS) space additionally splits the active-active case
by subspace, through ``subspace_of``. The three non-redundant blocks for a
CAS are exactly **inactive-active, inactive-virtual and active-virtual**;
``rotation_parameters`` enumerates them as a flat `(rows, cols)` list with
:math:`p > q` only, since :math:`\kappa_{pq} = -\kappa_{qp}` is one variable.

**The rows of** :math:`F` **depend on the orbital's class**:

.. math::

   F_{in} = 2\,(FI_{ni} + FA_{ni}) \quad (n\ \text{inactive}), \qquad
   F_{tn} = \sum_u D_{tu} FI_{nu} + \sum_{uvw} d_{tuvw}\,(nu|vw) \quad (t\ \text{active}), \qquad
   F_{an} = 0 \quad (a\ \text{virtual})

with :math:`FI` the inactive (closed-shell) Fock and :math:`FA` the mean field
of the active density alone, both built by one call each to
``build_fock_direct`` inside ``generalized_fock``. **Both** :math:`FI` and
:math:`FA` **and the accumulation over** :math:`D_{tu}`, :math:`d_{tuvw}` are
**exactly linear in** ``dm1``/``dm2``: every term is one of the density
matrices multiplying integrals that do not depend on the density. This is
stated by the module ("nothing here can have a Hessian that disagrees with its
gradient") but is worth making explicit for state averaging, since it is the
fact phase 1 leans on (below).

``one_index_fock`` differentiates :math:`F` along one rotation
:math:`\kappa`: it is *the same construction*, with each coefficient matrix in
turn replaced by :math:`C\kappa`, so

.. math::

   (H\kappa)_{pq} = 2\,(\tilde F_{qp} - \tilde F_{pq})

with :math:`\tilde F` built by ``one_index_fock`` exactly as ``F`` is built by
``generalized_fock``. ``orbital_hessian`` reads the full non-redundant Hessian
out of this by calling ``one_index_fock`` once per parameter (one unit
:math:`\kappa` per column) and symmetrising. **`one_index_fock` is also
exactly linear in `dm1`/`dm2`**, for the same reason: the two densities enter
only as multiplicative weights on integrals or one-index-transformed
integrals, never multiplied by each other.

``czt_mcscf_gradient``'s two-electron split
---------------------------------------------

``czt_mcscf_gradient`` (``mqc_czt_mcscf_gradient.f90``) is explicitly a
**single-state, no-Z-vector** gradient: valid only because a fully optimised
CASSCF is stationary in both orbitals and CI, so every response term is
identically zero and what remains is a contraction of differentiated integrals
against the densities already in hand.

The two-electron term is split by how many active indices it carries, to keep
the cost at :math:`n_{ao}^4` (once) rather than :math:`n_{ao}^4` per active
index combination:

.. math::

   E_{2e} = \frac12\, D \cdot (J - K/2)(D) \;+\; \frac12 \sum_{tuvw} \delta d_{tuvw}\,(tu|vw)

with :math:`D = D_{\text{inactive}} + D_{\text{active}}` the **total**
one-particle density in the AO basis, and

.. math::

   \delta d_{tuvw} = d_{tuvw} - D_{tu} D_{vw} + \tfrac12 D_{tv} D_{uw}

the active two-particle density with its own mean-field part subtracted
(``cumulant_two_particle_density``). This is an **algebraic identity** true
for any consistent :math:`(D, d)` pair, not a physicality assumption -- it
holds equally for a relaxed effective density, which is what phase 4 needs
(below). The first term goes through the same ``two_electron_deriv`` contraction
a closed-shell SCF gradient uses; only the second needs the general four-index
machinery (``active_two_electron_gradient``, ``gamma_block``), built from
:math:`n_{active}^4` rather than from the basis. At :math:`n_{active}=0`,
:math:`\delta d` is empty and the whole module reduces to the closed-shell SCF
gradient -- the check the module's own docstring names.

The Pulay term contracts the overlap derivative against the **generalised
Fock** :math:`\tfrac12(F + F^T)` (called ``weighted`` in the code), which for a
converged, single-state MCSCF is the energy-weighted density.  For SA-CASSCF
this is exactly the piece that changes: which :math:`(D, d)` (and hence which
:math:`F`) feeds this contraction is what phase 4's relaxed density answers.

``n_roots`` in CASCI
---------------------

``run_czt_casci``'s ``n_roots`` argument already asks ``davidson_lowest`` for
that many roots and returns every one of their energies and vectors (see
above); nothing about **converging several roots** is new. What is new is
*acting* on more than the first: no state averaging, no transition density,
and no Z-vector exist anywhere in this tree yet. ``run_czt_casscf``
(``mqc_czt_mcscf.f90``) always calls ``solve_ci`` with an implicit
``n_roots = 1`` (the optional argument is never passed), and takes the ground
state's density matrices alone into ``generalized_fock``/``orbital_gradient``/
``orbital_hessian`` each macro-iteration.

Spin
----

The CI (``mqc_ci.f90``, after PySCF's ``direct_spin1``) works in the
determinant basis at fixed :math:`M_S` and imposes no :math:`S`. Its lowest
roots can be of any spin, and that matters here: for CAS(2,2) at planar
C2H4/6-31G* the second root is the triplet, and at twisted C2H4 the triplet is
the lowest root. A singlet state average needs the CI restricted to singlets.
When :math:`n_\alpha = n_\beta`, a CI vector :math:`c(i_\alpha, i_\beta)`
of even :math:`S` is symmetric under :math:`i_\alpha \leftrightarrow
i_\beta` and one of odd :math:`S` is antisymmetric. Symmetrising every
Davidson guess and correction vector therefore excludes the triplets exactly,
which is what PySCF's ``direct_spin0`` does. A spin penalty
(``fix_spin_``) is the approximate alternative; on twisted C2H4 it stalled
PySCF's orbital gradient near :math:`10^{-6}`. The CI Lagrange multipliers
:math:`\bar c_{I,J}` live in the same symmetric subspace, so the phase 3
projector also symmetrises. Every root's :math:`\langle S^2\rangle` is
reported so that a spin mix-up is visible.

The SA-CASSCF Lagrangian and Z-vector, in these conventions
=============================================================

Orbitals are shared across states, so write the active-space integrals and
generalised-Fock machinery as functions of :math:`\kappa` alone (at fixed
reference orbitals :math:`C`), and each state's CI vector :math:`c_J` as an
eigenvector of the active-space Hamiltonian :math:`H(\kappa)` built at those
orbitals. Weights :math:`w_J \ge 0`, :math:`\sum_J w_J = 1`, are fixed (not
optimised).

.. math::

   E_J(\kappa, c_J) = \langle c_J | H(\kappa) | c_J \rangle, \qquad
   E_{SA}(\kappa, \{c_J\}) = \sum_J w_J\, E_J(\kappa, c_J)

**The key simplification this code gets for free.** Because
``generalized_fock``/``orbital_gradient``/``orbital_hessian`` are exactly
linear in :math:`(D, d)`,

.. math::

   \frac{\partial E_{SA}}{\partial \kappa}
   = \sum_J w_J \frac{\partial E_J}{\partial \kappa}
   = g\big(D_{SA}, d_{SA}\big), \qquad
   D_{SA} = \sum_J w_J D_J,\quad d_{SA} = \sum_J w_J d_J

i.e. the **SA orbital gradient and orbital Hessian are `orbital_gradient` and
`orbital_hessian` called once, at the weighted-summed densities** -- not a sum
of :math:`N` separate gradients/Hessians. This is what makes phase 1 (the SA
*energy*) a small change: build :math:`D_{SA}, d_{SA}` from every root's
``active_space_rdms`` and hand them to the orbital optimiser exactly as today's
single-root :math:`D, d` are.

The Lagrangian for one root's energy :math:`E_I`
--------------------------------------------------

The orbitals are stationary for :math:`E_{SA}`, not for :math:`E_I` alone, so
:math:`E_I` needs a Lagrangian carrying that constraint and the CI
eigenvalue constraints for every state:

.. math::

   L_I = E_I(\kappa, c_I)
       + \bar\kappa_I \cdot \left.\frac{\partial E_{SA}}{\partial \kappa}\right|_{\kappa=0}
       + \sum_J w_J\, \bar c_{I,J} \cdot \big(H(\kappa{=}0) - E_J\big) c_J

:math:`\bar\kappa_I` (an antisymmetric matrix, one non-redundant parameter per
entry of ``rotation_parameters``) and :math:`\bar c_{I,J}` (one vector per
state :math:`J`, same shape as :math:`c_J`) are Lagrange multipliers,
determined by making :math:`L_I` stationary in :math:`\kappa` and every
:math:`c_J` at the converged point. That stationarity condition is the
**Z-vector / CP-MCSCF equation**:

.. math::

   H_{SA}\,\begin{bmatrix}\bar\kappa_I \\ \{\bar c_{I,J}\}\end{bmatrix}
   = -\begin{bmatrix}\partial E_I/\partial\kappa \\ \{\partial E_I/\partial c_J\}\end{bmatrix}

:math:`H_{SA}` is the Hessian of :math:`E_{SA}` (not of :math:`E_I`) in the
joint orbital+CI space, and **does not depend on which root's RHS is being
solved** -- the fusion argument the plan is built on.

The right-hand side, precisely
--------------------------------

**CI part is zero.** :math:`c_I` is an exact eigenvector of :math:`H(\kappa=0)`
(the Davidson converged it to whatever tolerance ``run_czt_casci`` was given),
so :math:`\partial E_I/\partial c_J` vanishes identically for :math:`J = I`
(the usual variational-eigenvector argument) and is not defined/needed for
:math:`J \ne I` (:math:`E_I` does not depend on another state's coefficients
at all). PySCF's own code computes this generically through
``newton_casscf.gen_g_hop`` on a *single-state* CASCI object pinned to
:math:`c_I`, and that CI block is zero to numerical precision for a converged
CI -- confirmed by reading ``get_wfn_response`` in ``sacasscf.py``, which
builds the joint gradient of a plain (non-averaged) CASSCF at
:math:`(\kappa=0, c_I)` and only ever *inserts* it at state :math:`I`'s slot,
leaving every other state's CI RHS at zero. So in this code's language:

.. math::

   \text{RHS}_{\text{orb}} = \left.\frac{\partial E_I}{\partial\kappa}\right|_{c_I}
   = \texttt{orbital\_gradient}\big(\texttt{generalized\_fock}(D_I, d_I)\big), \qquad
   \text{RHS}_{\text{CI},J} = 0 \ \ \forall J

i.e. the RHS orbital block is **exactly today's single-root orbital gradient**,
evaluated at root :math:`I`'s own densities but at the shared (SA-optimised)
orbitals -- one call to routines that already exist, no new physics. Only
:math:`H_{SA}` (left-hand side) is new.

The SA Hessian blocks
-----------------------

**Orbital-orbital.** By the same linearity argument as the gradient, this
block is ``orbital_hessian`` called once at :math:`(D_{SA}, d_{SA})` --
weighted SA densities, not a sum of per-root Hessians.

**CI-CI, block** :math:`J`. Acting on a trial vector :math:`\bar c_J`, this
block is :math:`w_J \cdot 2(H(\kappa{=}0) - E_J)\,\bar c_J` -- the ordinary CI
Hamiltonian shifted by that state's own energy, scaled by its weight, and
**block-diagonal across states before projection** (a state's CI space only
talks to itself through this block; the Hessian couples states only through
the orbital-CI block below and through the redundancy projector).

**Orbital-CI coupling.** The mixed second derivative
:math:`\partial^2 E_{SA}/\partial\kappa\,\partial c_J`, weighted by :math:`w_J`
and summed over :math:`J`. This is the one block with no existing mqc
counterpart at all: it is what ``one_index_fock`` becomes when *one* of the
transformed quantities is a transition density between :math:`c_J` and a CI
trial vector rather than a one-index-rotated density -- i.e. it needs
transition RDMs (phase 2) composed with the ``one_index_fock``-style
machinery (phase 3).

**Redundancy projection.** With equal weights, :math:`E_{SA}` is the trace of
:math:`H` over :math:`\mathrm{span}\{c_1,\dots,c_N\}`, so rotations that mix
the averaged states into each other leave it unchanged: they are redundant and
are projected out. With unequal weights they are not redundant. A rotation by
:math:`\theta` between eigenstates :math:`J` and :math:`K` changes
:math:`E_{SA}` by :math:`(w_J - w_K)(E_K - E_J)\,\theta^2`, so the direction has
non-zero curvature and has to stay in the Hessian. PySCF's gradient code does
not handle that case: ``Gradients.__init__`` in ``pyscf/grad/sacasscf.py``
raises ``NotImplementedError`` when :math:`\max(w) - \min(w) > 10^{-8}`. The
gradient here is likewise **equal weights only**; the energy (phase 1) takes
any weights.

The projector PySCF actually applies (``project_Aop`` in the same file) removes,
from the CI part of a Hessian-vector product for state :math:`i`, its overlap
with **every** :math:`c_j` in the SA space sharing its spin sector, not only
:math:`c_i`:

.. math::

   (Ax_{ci})_i \;\leftarrow\; (Ax_{ci})_i - \sum_{j:\ \text{spin}(j)=\text{spin}(i)} \langle (Ax_{ci})_i, c_j\rangle\, c_j

i.e. **project onto the complement of** :math:`\mathrm{span}\{c_1,\dots,c_N\}`
**restricted to the matching spin/spatial-symmetry sector**, applied
separately to every :math:`\bar c_J`'s residual during the iterative solve.
This is the precise form to implement in phase 3; "orthogonal to :math:`c_J`
alone" (only the same-index state) is not what PySCF does and would leave a
near-null direction in the Hessian whenever two averaged states are close in
energy.

Relaxed density assembly for :math:`dE_I/dR`
================================================

PySCF's own total gradient (``lagrange.Gradients.kernel``) is literally a sum
of two pieces, computed by two independent calls:

.. math::

   \frac{dE_I}{dR} = \underbrace{\frac{dE_I}{dR}\bigg|_{\kappa,\,c_I\ \text{fixed}}}_{\texttt{get\_ham\_response}}
   \;+\; \underbrace{\bar\kappa_I \cdot \frac{d}{dR}\frac{\partial E_{SA}}{\partial\kappa}
   + \sum_J w_J\, \bar c_{I,J}\cdot\frac{d}{dR}(H-E_J)c_J}_{\texttt{get\_LdotJnuc}}

The first term is **exactly today's `czt_mcscf_gradient`, called on root**
:math:`I`\ **'s own** :math:`(D_I, d_I)` at the shared SA orbitals -- no new
code. The second term ("the Lagrange response") is what is new, and it
decomposes into an orbital-response piece and a CI-response piece that both
land on the *same* contraction machinery ``czt_mcscf_gradient`` already has,
fed a different effective density:

Orbital-response piece (from :math:`\bar\kappa_I`)
-----------------------------------------------------

:math:`\bar\kappa_I` one-index-transforms the **SA** densities exactly the way
``one_index_fock`` already builds its :math:`d_{\text{inactive}}`,
:math:`d_{\text{active}}` intermediates from a trial :math:`\kappa`:

.. math::

   \tilde D^{\bar\kappa_I} = \big[\bar\kappa_I, D_{SA}\big]\text{-type one-index transform}, \qquad
   \tilde d^{\bar\kappa_I} = \text{the four-term one-index transform of } d_{SA}

(PySCF's ``Lorb_dot_dgorb_dx`` builds exactly this pair, calling them
``dm1L``/the implicit two-particle piece, from ``Lorb`` -- the code's
:math:`\bar\kappa_I`.) The resulting effective one/two-particle density is fed
to the derivative-integral contraction the same way :math:`(D, d)` are today.

CI-response piece (from :math:`\bar c_{I,J}`)
------------------------------------------------

The symmetrised **transition** density between the Lagrange multiplier and the
state it is attached to, summed over the SA space:

.. math::

   \tilde D^{ci} = \sum_J w_J\, \big(\langle \bar c_{I,J} | E_{pq} | c_J\rangle + \langle c_J | E_{pq} | \bar c_{I,J}\rangle\big), \qquad
   \tilde d^{ci} \text{ likewise from the transition 2-RDM}

matching PySCF's ``Lci_dot_dgci_dx``, which forms ``trans_rdm12(Lci, ci)``
symmetrised. This needs the transition-RDM routine phase 2 adds beside
``active_space_rdms`` (bra :math:`\ne` ket).

Overlap / energy-weighted-density term
-----------------------------------------

The Pulay term contracts the overlap derivative against the generalised Fock
of *whichever* density the energy expression is built from. Because
``generalized_fock`` is linear in :math:`(D,d)`, phase 4 does **not** need
PySCF's three-way split of ``dme0`` into a root, an orbital-response and a
CI-response generalised Fock (three separate ``get_jk`` calls in PySCF): one
call to a `generalized_fock`-shaped contraction on the **total relaxed**
density

.. math::

   D_I^{\text{relaxed}} = D_I + \tilde D^{\bar\kappa_I} + \tilde D^{ci}, \qquad
   d_I^{\text{relaxed}} = d_I + \tilde d^{\bar\kappa_I} + \tilde d^{ci}

gives the same energy-weighted density in one pass, and the cumulant split
(``cumulant_two_particle_density``) applies unchanged, since it is an algebraic
identity in whatever :math:`(D,d)` pair it is given (stated above). This is
the concrete target for phase 4: build one relaxed :math:`(D_I, d_I)` pair per
root and hand it to a generalised
``czt_mcscf_gradient``-shaped routine, rather than porting PySCF's
term-by-term decomposition.

Checks
========

- :math:`n_{states}=1` must give :math:`\bar\kappa_I = 0` and every
  :math:`\bar c_{I,J}=0` **exactly**: with one state the RHS orbital block is
  the only nonzero piece of the RHS, and it already equals
  :math:`\partial E_{SA}/\partial\kappa` (since :math:`E_{SA}=E_1`), so the
  Z-vector equation reads :math:`H_{SA}\,\bar\kappa_I = -\partial
  E_{SA}/\partial\kappa`, which is solved at :math:`\bar\kappa_I=0` because a
  converged CASSCF has :math:`\partial E_{SA}/\partial\kappa = 0` already. The
  relaxed density then collapses to :math:`D_1, d_1` and phase 4 must
  reproduce ``czt_mcscf_gradient`` **bit for bit** on this path.
- :math:`\sum_I w_I\, dE_I/dR = dE_{SA}/dR`, and the right-hand side needs no
  response at all (:math:`E_{SA}` is stationary in both :math:`\kappa` and
  every :math:`c_J` by construction) -- a pure sanity identity with no
  tolerance-fitting escape hatch.

Reference systems
-----------------

``tools/sa_casscf/pyscf_ref.py`` produces the PySCF numbers.

- ``c2h4_twisted.xyz`` is twisted by 90 degrees **and** has its second CH2
  group pyramidalised. At the unpyramidalised D2d geometry, SA-2-CAS(2,2)
  breaks the symmetry: the zwitterionic S1 puts its charge on one carbon, so
  there are two mirror-image solutions with equal energies and gradients that
  swap C1 and C2. An implementation can converge to either one, so a correct
  gradient can fail a comparison. Pyramidalising one end leaves one lowest
  solution.
- PySCF's SA-CASSCF orbital gradient stalls near :math:`10^{-7}` even with the
  Newton solver. A root energy is not stationary in the orbitals, so that
  residual leaves about :math:`10^{-6}` Hartree/Bohr of scatter in a central
  difference with :math:`h = 10^{-3}`. PySCF's analytic gradient is the
  sharper reference, and finite differences of this code's own energies,
  which converge further, arbitrate below :math:`10^{-6}`.

The fused-kernel data flow
=============================

What already exists, and can likely be reused as-is
-------------------------------------------------------

**The multi-density Fock build the plan calls "new" already exists.**
``build_fock_direct_many`` (``backends/cenzontle/integrals/mqc_czt_direct.f90``)
builds :math:`F = H + J - K/2` for an arbitrary stack of densities,
``densities(n_ao, n_ao, n_set)``, over **one** pass of the shell-quartet loop,
with Schwarz screening shared across the whole batch (so a batch is bit-for-bit
the same as the densities run one at a time). It already supports an
antisymmetric-density mode (for the antisymmetric one-index-transformed
densities phase 3/4/5 need), a ``k_scale``/``j_scale`` split, and range
separation. It is already used by ``mqc_czt_response_product.f90`` and
``mqc_czt_cphf.F90`` (the RHF CPHF solver) to batch several perturbations'
Fock builds together. **Phase 5's "multi-density Fock build variant" should
start by trying `build_fock_direct_many` directly**, and only write something
new if the SA Hessian's density shapes (inactive+active decomposition,
occupied-column-only potentials as in `transformed_potential`) do not fit its
`(n_ao, n_ao, n_set)` contract.

**The block-CG-with-many-RHS pattern already exists**, in
``mqc_czt_cphf.F90``'s ``cphf_solve``: it takes ``perturbations(n_ao, n_ao,
n_perturbations)``, solves every right-hand side together through
``response_operator``/``response_product`` (which is what calls
``build_fock_direct_many``), and returns one response `U_ai` per perturbation.
This is architecturally the template for the SA block PCG (phase 3/5) -- same
"one operator apply builds every RHS's Fock at once" idea -- but it is *not*
reusable directly: CPHF solves only the occupied-virtual MO block of a
single-determinant reference, with no CI coupling and no inactive-active or
active-virtual distinction. The plan's own phrase, "a pattern to borrow, not
reuse," is confirmed by reading it.

What does not exist and needs a stacked variant
----------------------------------------------------

``one_electron_deriv``/``two_electron_deriv``/``iprinv_deriv_at``
(``backends/cenzontle/derivatives/mqc_czt_gradient.f90``) take **one**
density (or none, for the raw derivative-integral matrices) and return **one**
gradient contribution. Concretely:

- ``one_electron_deriv(mol, matrix, which)`` returns the bare derivative
  integral matrix ``(n_ao, n_ao, 3)`` -- it takes no density at all, so it
  needs no stacking; it is reused unchanged for every root, and only the
  *contraction* against it (done by hand in ``czt_mcscf_gradient``, not inside
  this routine) needs to run once per root or be batched.
- ``two_electron_deriv(mol, density, vhf, error, ...)`` **does** take a single
  ``density(:,:)`` and return a single ``vhf(:,:,:)``. This is the routine a
  fused multi-root pass needs a stacked twin of: `two_electron_deriv_many(mol,
  densities(n_ao,n_ao,n_set), vhfs(n_ao,n_ao,3,n_set), ...)`, sharing the
  shell-quartet loop and the screening the same way
  ``build_fock_direct_many`` shares it for the ordinary (non-derivative) case.
  No such routine exists yet.
- The active two-electron gradient path (``active_two_electron_gradient``,
  ``gamma_block``, both in ``mqc_czt_mcscf_gradient.f90``) is per-root already
  in a different sense: it is driven by ``ddm2``, which is :math:`n_{active}^4`
  and cheap to hold for every root simultaneously, so the natural fusion there
  is stacking the **first-index AO block** dimension across roots inside
  ``gamma_block``'s existing shell-blocked loop, not touching
  ``two_electron_mp2_terms`` itself (it already takes one ``gamma_blk`` per
  call; a stacked call would pass a wider one and expect a wider gradient
  return, which is the actual interface change needed there).

File-level plan, phase by phase
==================================

Phase 1 -- SA-CASSCF energy (this branch, after this document)
-----------------------------------------------------------------

**Confirmed exactly as the plan states**: ``mqc_json_schema.f90``'s
``mcscf_keys()`` allow-lists neither ``n_states`` nor a weights key, and its
own docstring says why -- "`mcscf_config_t` carries fields for state averaging
... and none of that is implemented." The fields already exist, at every layer
that is *not* the JSON path:

- ``mcscf_config_t`` in ``src/methods/dispatch/mqc_method_config.f90`` already
  has ``n_states`` (default 1) and ``state_weights`` (unallocated).
- ``mcscf_options_t`` in ``src/methods/dispatch/mqc_method_mcscf.f90`` already
  has the same two fields, with a comment stating plainly they are "not
  reachable from a deck."
- ``mqc_method_factory.F90``'s ``configure_mcscf`` already copies
  ``config%mcscf%n_states``/``state_weights`` into ``m%options``.

**What is actually missing**, traced end to end:

1. The **flat** ``mqc_config_t`` in ``src/io/mqc_config_types.f90`` (what the
   JSON reader fills in directly) has no ``mcscf_n_states``/
   ``mcscf_state_weights`` fields at all -- add them, defaulting to ``1`` and
   unallocated, next to the other ``mcscf_*`` fields.
2. ``mcscf_keys()`` in ``mqc_json_schema.f90`` must allow ``"n_states"`` and a
   weights key (``"weights"``, an array) -- and its docstring's justification
   for refusing them needs to be replaced, not just its allow-list.
3. ``mqc_json_config_reader.f90`` needs an ``optional_int`` for
   ``keywords.mcscf.n_states`` and a new array reader for
   ``keywords.mcscf.weights`` (there is no existing ``optional_real_array``
   helper for a top-level array at time of writing; check
   ``read_ormas_partition`` in the same file for the pattern a JSON array of
   numbers is read with, since ORMAS's ``subspaces`` is the nearest existing
   example of an array key, even though its elements are integers).
4. ``mqc_config_adapter.f90`` needs to copy the two new flat fields into
   ``driver_config%method_config%mcscf%n_states``/``%state_weights`` --
   currently **not copied at all**, which is why ``config%mcscf%n_states``
   always reads back its default of 1 regardless of what a deck might someday
   say.
5. ``src/methods/dispatch/mqc_method_mcscf.f90``'s ``mcscf_run`` builds
   ``settings%mcscf`` (a ``cuest_scf_settings_t``'s ``mcscf_config_t``,
   ``mqc_cuest_iface.f90``) field by field from ``this%options``, and **does
   not copy** ``n_states``/``state_weights`` into it -- add that copy.
6. ``run_czt_mcscf`` in ``backends/cenzontle/mqc_czt_bridge.f90`` must refuse
   ``settings%mcscf%n_states > 1`` combined with
   ``allocated(settings%mcscf%ormas_subspaces)`` by name (SA + ORMAS is out of
   scope; there is no transition-density machinery for a restricted space and
   no reason to build it before the CAS case works), following the existing
   refusal style in that routine (``result%error%set(ERROR_VALIDATION, ...)``,
   ``result%has_error = .true.``, early ``return`` -- see the PCM/charges/
   bond-order refusals a few lines above the CASSCF dispatch, or
   ``mcscf_gradient_into``'s CASCI-gradient refusal).
7. **Singlet selection** (see Spin, above): when the target multiplicity is 1
   and :math:`n_\alpha = n_\beta`, the Davidson symmetrises its guess and
   correction vectors under alpha/beta transposition, behind an option so
   that single-state CASSCF and CASCI results do not move. Every root's
   :math:`\langle S^2\rangle` goes to the output.
8. ``run_czt_casscf`` (``mqc_czt_mcscf.f90``) needs new optional
   ``n_states``/``weights`` arguments. Inside the macro loop, replace
   ``call solve_ci(..., ci, ...)`` (implicit one root) with a request for
   ``n_states`` roots, then build :math:`D_{SA}, d_{SA}` as
   :math:`\sum_J w_J D_J, \sum_J w_J d_J` by calling ``active_space_rdms`` once
   per root on ``ci%vectors(:, :, J)`` (the restricted-space ORMAS path is
   refused before this point, so no equivalent change is needed in
   ``ormas_density_matrices`` for phase 1). Feed :math:`(D_{SA}, d_{SA})`
   into ``generalized_fock``/``orbital_gradient``/``orbital_hessian`` exactly
   where ``dm1``/``dm2`` are used today -- this is the linearity payoff above:
   no change to any of those three routines. ``casscf_result_t`` needs an
   ``energies(n_states)`` field (mirroring ``casci_result_t``) so every root's
   energy reaches the bridge; ``%dm1``/``%dm2`` should probably become the SA
   densities (what the orbital optimiser actually used), with a note that a
   per-root density is not carried here yet (phase 4 territory).
9. JSON output: every root's energy plus :math:`E_{SA}`, following the
   ``excitation_energies``-style array pattern already in
   ``mqc_json_output_types.f90``/``mqc_json_writer.f90`` for excited states.

Phase 2 -- transition RDMs
------------------------------

New routine beside ``active_space_rdms`` in ``mqc_rdm.f90`` (proposed:
``transition_space_rdms(bra, ket, alpha, beta, tdm1, tdm2, error)``), same
spin-traced, chemist-ordered convention, ``bra /= ket`` in general. Reuses the
same ``apply_excitations``/``excitations_block`` machinery
``active_space_rdms`` already calls, contracting the bra's excited vectors
against the ket's rather than a vector against itself. Gate: ``bra = ket``
reproduces ``active_space_rdms`` bit for bit; orthogonal states give trace
zero; compare against PySCF's ``trans_rdm12``.

Phase 3 -- matrix-free SA Hessian-vector product
----------------------------------------------------

New module (proposed ``mqc_czt_sa_hessian.f90``), building the HVP described
above: orbital-orbital block via ``one_index_fock`` at SA densities (reuse
unchanged), orbital-CI coupling via ``one_index_fock``-style transforms fed a
transition density from phase 2, CI-CI block via the existing
``sigma_operator_t``/Davidson machinery in ``mqc_ci``/``mqc_davidson.f90``
(the :math:`(H-E_J)` action is exactly what a sigma build already computes),
and the redundancy projector exactly as read out of PySCF's ``project_Aop``
above. Gate: the orbital block matches ``orbital_hessian`` (at SA densities) on
random vectors; the full HVP matches a finite difference of the SA gradient
(orbital+CI) to about :math:`10^{-7}`.

Phase 4 -- single-root gradient via Z-vector
------------------------------------------------

New routine (proposed ``czt_sa_casscf_gradient`` beside
``czt_mcscf_gradient`` in a new ``mqc_czt_mcscf_gradient`` sibling or the same
file): one PCG solve of :math:`H_{SA}\,x = -\text{RHS}_I` (one RHS), assembling
the relaxed :math:`(D_I, d_I)` per the section above, and calling a
generalised ``czt_mcscf_gradient``-shaped contraction on it. Gate: vs PySCF
per root and vs finite differences; :math:`n_{states}=1` bit-identical to
``czt_mcscf_gradient``; the weighted-sum identity; C2H4 planar and twisted.

Phase 5 -- fused multi-root
-------------------------------

Block PCG (all :math:`N` RHS at once, reusing phase 3's HVP applied to a block
of trial vectors), the multi-density Fock build (try
``build_fock_direct_many`` first, per above), and the stacked
derivative-integral contraction (``two_electron_deriv_many``, new, per above).
``keywords.mcscf.gradient_roots`` (``"all"`` default under SA, or an explicit
list) added the same way as the phase 1 keys. Gate: per-root agreement with
phase 4 to :math:`\le 10^{-10}`; repeated timings, :math:`N=1..4` on PSB3,
reported as cost(N)/cost(1).

Phase 6 -- driver/JSON, docs, example deck
-----------------------------------------------

All-root gradients and pairwise differences in JSON output; a user-facing doc
page (``mqc_docs/source/sa_casscf.rst`` or folded into an existing MCSCF page);
an example deck, C2H4/6-31G*/SA-2-CAS(2,2); a ``run_validation.py`` entry if
cheap enough for the default manifest.

Open questions
==============

- **Unequal weights.** The gradient is equal weights only, as in PySCF.
  Unequal weights would need the in-space rotations kept in the Hessian, not
  projected (see Redundancy projection).
- **Spin selection.** Singlets by alpha/beta symmetrisation of the CI vectors,
  as in PySCF's ``direct_spin0``. Other target spins are not handled.
