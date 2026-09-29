======================
State-Averaged CASSCF
======================

A state-averaged CASSCF (SA-CASSCF) optimises one set of orbitals against
several CI roots at once, rather than the lowest one alone. It is what a
near-degeneracy or a conical intersection needs: optimising the orbitals for
one root only can make them a poor description of a nearby state, or lose
track of which root is which as a geometry changes. Asked for with
``keywords.mcscf`` on an ordinary ``"method": "casscf"`` deck -- state
averaging is a property of how many roots the orbitals are optimised
against, not a different method.

Keys
====

.. code-block:: json

   "mcscf": {
     "n_active_electrons": 2,
     "n_active_orbitals": 2,
     "n_states": 2,
     "weights": [0.5, 0.5],
     "gradient_roots": "all"
   }

- ``n_states`` (default: 1, an ordinary CASSCF): how many roots the orbitals
  are optimised against.
- ``weights`` (default: equal, ``1/n_states``): one non-negative weight per
  state, summing to one. Any weights are accepted for the energy. A
  **gradient needs equal weights** (below).
- ``gradient_roots`` (default: every state): on a ``"driver": "Gradient"``
  deck only, which roots to differentiate -- the string ``"all"`` or an
  explicit list of 1-based root indices, e.g. ``[1, 3]``. Refused, rather
  than silently ignored, when given with ``n_states = 1`` (there is only one
  root) or off a Gradient driver (it would never be read).

**Singlets only, for now.** State averaging needs the reference to be a
singlet -- equal active alpha and beta electrons -- and is refused
otherwise. The determinant CI imposes no total spin by itself, so its
lowest roots can be of any spin (two electrons in two orbitals already has a
triplet among its four lowest determinants), and an unrestricted state
average would mix whatever spins the Davidson happened to converge to. With
a singlet reference the Davidson is restricted to the symmetric
alpha/beta subspace whenever ``n_states > 1``, which excludes every
non-singlet root exactly rather than approximately. Every state's
``<S^2>`` reaches the output so a spin mix-up would be visible if the
restriction were ever wrong.

Refused rather than silently ignored
=====================================

- ``n_states > 1`` together with ``ormas``: there is no transition-density
  machinery for a restricted active space.
- ``n_states > 1`` on a CASCI (``optimize_orbitals: false``): averaging is a
  property of what the orbitals are optimised against, and a CASCI never
  moves them.
- ``n_states > 1`` on a ``Hessian`` driver: there is no analytic CASSCF
  Hessian at all yet, single-state or averaged.
- Unequal weights on a ``Gradient`` driver: a rotation between two
  unequally-weighted averaged states is not redundant and needs a
  curvature term this code does not carry, exactly as in PySCF's own
  SA-CASSCF gradient code.
- ``gradient_roots`` with ``n_states = 1``, or off a Gradient driver: it
  would be parsed and then never read.
- State averaging under fragmentation (``keywords.fragmentation``): nothing
  yet carries a fragment's per-root energies, spins or gradients through the
  fragment machinery or MPI packing, so a fragmented deck with
  ``n_states > 1`` is refused rather than silently reporting only the
  averaged total.
- CASPT2/NEVPT2 corrections are not implemented, and no keyword accepts
  them, state-averaged or not.

What the gradient is
=====================

The orbitals are optimised for :math:`E_{SA} = \sum_J w_J E_J`, not for any
one root, so a single root's energy is not stationary in the orbitals and
needs its own Lagrangian (the Z-vector / CP-MCSCF equation) to
differentiate correctly. Every requested root shares one Hessian state and
one block preconditioned-CG solve, so asking for several roots together
costs less than computing them one at a time. See :doc:`developer_sa_casscf`
for the equations.

**The top-level ``gradient`` is** :math:`dE_{SA}/dR` **, not any one root's
own gradient.** It equals the weight-averaged sum of every root's own
gradient exactly (the SA energy needs no response to differentiate, unlike
any single root), and it is what an SA geometry optimisation should follow.
A root's own gradient is under ``mcscf_states`` in the output, below.

Output
======

Every root's energy, ``<S^2>`` and weight reach the JSON output under an
``mcscf_states`` section, alongside ``E_SA``. On a Gradient driver, every
root named in ``gradient_roots`` additionally carries its own gradient and
gradient norm, and a ``gradient_differences`` array holds :math:`g_i - g_j`
for every pair of requested roots -- what a conical-intersection search or a
surface-hopping trajectory needs to find where two surfaces come close:

.. code-block:: json

   "mcscf_states": {
     "n_states": 2,
     "e_sa_hartree": -77.864049844180,
     "states": [
       {
         "state": 1, "energy_hartree": -77.920166075660, "s2": 0.0,
         "weight": 0.5,
         "gradient_norm": 0.159,
         "gradient_units": "hartree/bohr",
         "gradient": [[...], "..."]
       },
       {
         "state": 2, "energy_hartree": -77.807933612710, "s2": 0.0,
         "weight": 0.5,
         "gradient_norm": 0.068,
         "gradient_units": "hartree/bohr",
         "gradient": [[...], "..."]
       }
     ],
     "gradient_differences": [
       {"states": [1, 2], "gradient_norm": 0.137, "gradient": [[...], "..."]}
     ]
   }

``total_energy`` at the top level of the same output is ``E_SA``, not root
1's own energy -- said here rather than left to be assumed, since it is the
one fact a consumer reading the numbers cannot see on its own.

At ``info`` log level, every requested root's energy, gradient norm and the
Z-vector solver's CG iteration count print as the run goes; at ``verbose``,
each root's full gradient table follows.

Example deck
============

``tools/sa_casscf/c2h4_twisted_sa2_grad_6-31gs.json``: twisted ethylene
(twisted 90 degrees, with the second CH2 group pyramidalised to leave a
single lowest solution -- see :doc:`developer_sa_casscf` for why), SA-2 over
a CAS(2,2), 6-31G*, both roots' gradients on a Gradient driver:

.. code-block:: json

   {
     "molecules": [{"xyz": "c2h4_twisted.xyz", "molecular_charge": 0,
                    "molecular_multiplicity": 1}],
     "model": {"method": "casscf", "basis": "6-31g*"},
     "keywords": {
       "scf": {"maxiter": 200, "tolerance": 1e-12},
       "mcscf": {
         "n_active_electrons": 2, "n_active_orbitals": 2,
         "max_macro_iter": 400, "orbital_convergence": 1e-10,
         "n_states": 2, "weights": [0.5, 0.5], "gradient_roots": "all"
       }
     },
     "driver": "Gradient"
   }

``tools/sa_casscf/pyscf_ref.py`` produces the independent PySCF reference
this deck's numbers are checked against.

Nonadiabatic couplings
======================

A Gradient run can also return the derivative coupling
:math:`d_{IJ} = \langle \Psi_I|\partial/\partial R\,\Psi_J\rangle` between
pairs of averaged roots. List the pairs, 1-based:

.. code-block:: json

   "mcscf": {"n_states": 2, "gradient_roots": "all", "nac_pairs": [[1, 2]]}

Each pair adds an entry to ``mcscf_states.nonadiabatic_couplings``:

- ``coupling``: :math:`d_{IJ}`, per atom, in 1/Bohr.
- ``interstate_coupling``: :math:`h_{IJ} = (E_J - E_I)\,d_{IJ} =
  \langle I|\partial H/\partial R|J\rangle`, in Hartree/Bohr. PySCF's
  ``mult_ediff=True`` scales by :math:`E_I - E_J` instead, so it has the
  opposite sign.
- ``csf_term``: the part of :math:`h_{IJ}` from the determinant (CSF)
  derivative. It is included in ``coupling`` and ``interstate_coupling``,
  and corresponds to PySCF's ``use_etfs=False``. It is not translationally
  invariant on its own.
- ``energy_difference_hartree``: :math:`E_J - E_I`.

The overall sign of a coupling depends on the arbitrary phases of the CI
vectors, so it can differ from another code's by a factor of -1 per pair.
The same refusals as for gradients apply: equal weights only, not under
fragmentation. See :doc:`developer_sa_casscf`, "Nonadiabatic couplings",
for the equations. ``tools/sa_casscf/c2h4_twisted_sa2_nac_6-31gs.json`` is an
example deck.

Against PySCF 2.14's ``pyscf.nac.sacasscf``, up to that sign, largest
difference over all components of :math:`d_{IJ}` and :math:`h_{IJ}`, with and
without the CSF term:

=====================================================  =====================
system                                                 max difference
=====================================================  =====================
C2H4 planar, 6-31G*, SA-2-CAS(2,2)                     :math:`1\times10^{-8}`
C2H4 twisted, 6-31G*, SA-2-CAS(2,2)                    :math:`2\times10^{-7}`
H2O (one bond stretched), 6-31G, SA-3-CAS(4,4), pairs  :math:`9\times10^{-7}`
=====================================================  =====================

A symmetric molecule whose averaged states break its symmetry has two
mirror-image SA solutions, and a coupling computed at either is correct. This
happens for H2O SA-3 at the equilibrium geometry and for C2H4 at the
symmetric 90 degree twist. Two codes, or two runs, can land on different ones.

Accuracy against PySCF
=======================

Both roots' gradients, against PySCF 2.14's ``sacasscf`` Lagrange-multiplier
gradient (fed this repository's own basis JSON, not PySCF's internal
tables): LiH/STO-3G and C2H4/6-31G* planar and twisted agree to
:math:`3\times10^{-10}` to :math:`5\times10^{-9}` Hartree/Bohr, well inside
PySCF's own iterative-solver floor. A finite difference of this code's own
root energies (central, tight CASSCF convergence) agrees with the analytic
gradient to :math:`3\times10^{-8}` to :math:`3\times10^{-6}` on the twisted
geometry, where PySCF's own analytic gradient is the sharper of the two
references.
