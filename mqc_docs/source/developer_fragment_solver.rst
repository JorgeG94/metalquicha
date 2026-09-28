========================================
The fragment solver: design and progress
========================================

**Status: design agreed (see Decisions at the end); phase 1 built.**
The survey below describes the code as it is now; the sections after it
describe the planned change. They will be rewritten as each phase lands.

Plain MBE and GMBE already run any method: each fragment is a
``physical_fragment_t`` handed to ``qc_method_t%calc_energy``, built by
``create_method``. FMO, EE-MBE and EFMO do not work that way. They solve their
fragments inside the CPU backend, so they can add an embedding operator,
freeze orbitals at a cut and start from a guess density, and every one of
those solves calls ``run_czt_rhf`` directly. The rule this design enforces:

* A fragmentation scheme names no method. It hands a fragment to **one
  fragment-solver interface** and gets back an energy and a density.
* Each method states what it can do under a scheme through **one capability
  query**.
* An unsupported combination is refused in **one place**, by name, before any
  work starts.

Survey: where RHF is hard-coded today
=====================================

Scheme call sites
-----------------

``backends/cenzontle/fragments/mqc_czt_fmo.f90`` runs both FMO and EE-MBE.
EE-MBE is ``esp = "ptc"`` and ``expansion = "mbe"``.

* ``inner_scf``: the monomer SCF. ``run_czt_rhf`` has four branches, depending
  on whether ``h_extra=u`` and ``projector=proj`` are present. It uses
  ``scf%energy``, ``scf%density`` and ``scf%converged``. The internal energy
  is ``scf%energy - sum(scf%density*u)``. A fragment that does not converge is
  refused regardless of ``allow_crap_scf``.
* ``nmer_term``: dimers and n-mers. Four branches, with guess
  ``SCF_GUESS_PROJ`` from ``d_split`` (D_I ⊕ D_J). A second attempt uses
  ``level_shift = max(deck, 0.5)``; SOSCF is not the fallback because it
  refuses a projector. It computes ``e_internal = E - Tr(D u)`` and
  ``e_resp = Tr((D - d_split) u)``. Under ``expansion = "mbe"`` it subtracts
  neither.
* ``pieda_pair_term``: ``cholesky_occupied_orbitals(density, nelec/2)``, then
  ``hl_prime_energy``, which evaluates the RHF energy functional of D_HL.
* ``pieda_dispersion_term``: ``dispersion_apply(kind, "hf", ...)``.
* ``build_fragments``: refuses an odd electron count, because the method is
  closed-shell.
* ``es_dimer_energy``, ``embedding_operator``, ``cut_embedding``,
  ``local_coulomb`` and ``fragment_charges`` use only the fragment densities,
  so they do not depend on the method once a density is given.

``backends/cenzontle/fragments/mqc_czt_efmo.f90``:

* ``nmer_energy`` and ``cut_nmer_energy``: ``run_czt_rhf`` on the near
  n-mers. There is no embedding, by the design of EFMO; a projector is used
  under AFO.
* ``fragment_correlation``: ``run_czt_mp2`` or ``run_czt_ri_mp2`` on a
  reference's orbitals.
* ``fragment_potential`` calls ``make_efp_potential``, whose own
  ``run_czt_rhf`` supplies the monomer energy ``E_I^0``.

These stay Hartree-Fock by construction, and do not go through the solver:

* ``bond_lmo_set`` and ``bond_hybrid`` in ``mqc_czt_afo.f90``. The model
  system only supplies the frozen orbitals.
* ``make_efp_potential`` in ``mqc_czt_efp_potential.f90``. MAKEFP is a
  Hartree-Fock construction.

Refusals and gaps in ``src``
----------------------------

* ``fmo_method_refusal`` in ``mqc_many_body_expansion.f90``, called from
  ``run_fragmented_calculation``, refuses every method except HF under FMO and
  EE-MBE. This is the stopgap that this work removes.
* ``run_efmo_energy`` in ``mqc_driver.f90`` refuses the following, each by
  name: a non-Energy driver, SCS/SOS-MP2, any method other than HF, MP2 or
  RI-MP2, and ``unrestricted``. ``run_efmo`` in the backend refuses MP2 across
  a cut.
* **Gap:** an FMO or EE-MBE deck with ``driver: "Gradient"`` is not refused.
  ``fmo_run_serial`` and ``fmo_run_distributed`` ignore ``calc_type`` and
  report an energy.
* **Gap:** ``model.unrestricted`` is silently ignored under FMO.

The SCF settings drift
----------------------

The duplication of ``hf_options_t`` and ``dft_options_t`` was closed on
2026-08-29 (``bb95a4867a``, ``05b9631e20``, ``e90f30c12b``). Both now extend
``scf_options_t``, which extends ``scf_numerics_t``, and ``apply_scf_settings``
is the one copy into ``cuest_scf_settings_t``.

What is still duplicated sits in ``mqc_driver.f90``: four hand copies of
``method_config%scf`` into a ``scf_numerics_t``. They are
``expansion%scf_drive`` (FMO and EE-MBE), ``efmo_scf``, ``makefp_scf`` and
``neo_scf``, and each copies a different subset:

* The FMO copy has no ``second_order`` or ``soscf_start``.
* The MAKEFP and NEO copies have no ``convergence_metric`` and no
  ``allow_crap_scf``.
* EFMO's tolerances are constants in ``run_efmo_energy``.

Method entry points
-------------------

* ``run_czt_rhf`` is the only entry that takes ``h_extra``, ``projector`` and
  ``guess_density``. Restricted Kohn-Sham is the same routine with ``xc=``
  present, and range-separated hybrids are handled inside ``xc_context_t``.
* ``run_czt_uhf`` takes none of the three, so under this design there is no
  unrestricted fragment solve.
* ``run_czt_mp2``, ``run_czt_ri_mp2``, ``run_czt_rccsd`` and ``run_czt_ccsd``
  take a converged reference's orbitals and orbital energies. They rebuild
  only the two-electron integrals and assume canonical orbitals. An SCF solved
  with ``h_extra`` therefore carries its embedding into them through the
  orbital energies, with nothing else needed. With a projector, the frozen
  virtuals sit at the shift (1e3 Eh) and would still be correlated. Nothing
  excludes them yet.
* D3 and D4 are added in the method layer (``dft_run`` via
  ``dispersion_apply``), not in the backend.
* ``run_czt_hf``, the unfragmented dispatcher, passes neither ``h_extra`` nor
  ``projector``. It returns no density.

The design
==========

Where it lives
--------------

``backends/cenzontle/fragments/mqc_czt_fragment_solver.f90`` (phase 1b, built).
It has to be in the backend because the embedding operator, the projector and
the molecule are all backend objects built over the fragment's own basis. FMO,
EE-MBE and EFMO call it; nothing in those three calls ``run_czt_rhf`` any
more.

The method arrives as the ``cuest_scf_settings_t`` that an unfragmented run of
the same deck would build. ``hf_run`` and ``dft_run`` build that settings
object today; the part that does it is factored out into a routine the driver
can call without running anything. ``run_czt_fmo`` and ``run_czt_efmo`` take
it in place of ``basis_name`` and ``scf_drive``, and their stubs change with
them. A fragment is then solved with the same settings as the whole molecule
would be, including the functional, grid, density fitting, MP2 variant and
frozen core. ``keywords.fragmentation.fmo_scf_*`` still override the
fragment tolerances, as they do now.

The same change closes the drift: ``scf_numerics_t`` is the parent component
of ``scf_options_t``, so the four hand copies become one conversion in one
place. That place is used by MAKEFP and NEO as well.

The interface
-------------

.. code-block:: fortran

   type :: fragment_request_t
      real(dp), allocatable :: h_extra(:, :)           ! embedding u over mol's AOs
      type(fock_projector_t), allocatable :: projector  ! allocated only when it holds frozen orbitals
      integer, allocatable :: guess                      ! SCF_GUESS_*; unallocated = backend default
      real(dp), allocatable :: guess_density(:, :)      ! total density, e.g. D_I (+) D_J
      type(scf_numerics_t) :: drive                      ! how the SCF is driven
      integer :: max_iter                                ! fragment tolerances (fmo_scf_*)
      real(dp) :: energy_tol, density_tol
      real(dp), allocatable :: grad_tol
      logical :: retry_level_shift                       ! one retry at >= 0.5 Eh if unconverged
      logical :: verbose
   end type

   type :: fragment_outcome_t
      real(dp) :: energy        ! E: reference + correlation, including Tr(D u)
      real(dp) :: internal      ! E' = E - Tr(D_esp u); energy when there is no u
      real(dp) :: reference     ! the SCF energy
      real(dp) :: correlation   ! zero for HF (and DFT)
      real(dp), allocatable :: density(:, :)  ! D_esp: the reference SCF's total density
      type(rhf_result_t) :: scf                ! the reference SCF itself
      logical :: converged
   end type

   subroutine solve_fragment_method(method, mol, nelec, real_z, request, outcome, error, &
                                    reference, aux, label)

``method`` is the ``cuest_scf_settings_t``. The caller builds ``mol``, plus
``aux`` for RI-MP2. ``real_z`` holds the atomic numbers of the real atoms, for
the frozen-core count.

``reference`` is an optional, already converged ``rhf_result_t``. EFMO passes
MAKEFP's own SCF, so that the monomer is not solved twice; the solver then only
adds correlation. Inside the routine:

* **HF:** ``run_czt_rhf`` with exactly the arguments the call sites pass
  today.
* **DFT:** the same call with an ``xc_context_t`` built from ``settings``.
  Dispersion is handled according to the answer to question 2 below.
* **MP2, SCS/SOS-MP2 and RI-MP2:** HF as above, then ``run_czt_mp2`` or
  ``run_czt_ri_mp2`` on the embedded orbitals, scaled by ``scs_ss`` and
  ``scs_os``. The frozen core is counted from ``real_z``, so ghosts do not
  count. This is built, and EFMO uses it.
* **CC:** later, the same shape as MP2.

The level-shifted n-mer retry and E' are computed in the solver. An
unconverged SCF is reported in ``outcome%converged``, and each caller keeps its
own refusal of it. FMO still refuses regardless of ``allow_crap_scf``; EFMO
honours it.

Which density defines the ESP
-----------------------------

**HF and DFT.** Both are variational in D, and ``u`` enters linearly through
``h``. So D_esp is the SCF density, ``E' = E - Tr(D u)``, and
``Tr(ΔD u)`` uses the n-mer's SCF density against ``d_split``. For HF this is
what the code does today.

**MP2 family.** This follows GAMESS FMO-MP2. D_esp is the embedded HF
density, and every ESP, every set of charges and ``Tr(ΔD u)`` come from HF
densities. Correlation is added per fragment and per n-mer on top of the
embedded HF:

* ``E' = E'_HF + E_corr``.
* ``ΔE_IJ`` gains ``E_corr(IJ) - E_corr(I) - E_corr(J)``.

The relaxed MP2 density would be more rigorous, but it is not what GAMESS
does, it would need a Z-vector per fragment, and it would break comparability
with GAMESS references.

The capability query and the one refusal site
---------------------------------------------

``src/methods/dispatch/mqc_fragment_capabilities.f90`` (phase 1a, built). It
is in ``src`` because the refusal must happen in the driver, before any backend
work:

.. code-block:: fortran

   type :: fragment_needs_t          ! what the deck asks
      logical :: cut                 ! bond_breaking = "afo"
      logical :: unrestricted        ! model.unrestricted
      integer :: calc_type
      logical :: pieda
   end type

   type :: fragment_capabilities_t   ! what the method can do under the scheme
      logical :: runs, cut, unrestricted, gradient, pieda
   end type

   pure function fragment_capabilities(scheme, method_config) result(cap)
   function fragment_refusal(scheme, method_config, needs) result(message)

Every method that runs at all accepts an embedding operator, so embedding is
not a capability. ``pieda`` says whether the method has its own PIEDA terms.
Which schemes offer PIEDA is still decided by ``check_pieda_support``.

``fragment_refusal`` is the single refusal site. ``run_fragmented_calculation``
calls it for FMO and EE-MBE, and ``run_efmo_energy`` calls it for EFMO. It
replaced ``fmo_method_refusal`` and the method refusals in
``run_efmo_energy``. The backend's MP2-across-a-cut check in ``run_efmo`` is
now only a guard against a caller that skipped it. The message names the
method, the scheme and the key to change.

Refusals that depend on the data rather than the method stay where the data
is. An odd electron count after a cut is one example.

.. list-table:: Capabilities as planned, by phase
   :header-rows: 1

   * - method
     - FMO / EE-MBE
     - with a cut (AFO)
     - EFMO
     - PIEDA
   * - HF
     - yes (phase 1)
     - yes
     - yes
     - yes
   * - DFT, including RSH, D3/D4
     - yes (phase 2)
     - yes (phase 2)
     - no: MAKEFP is HF
     - phase 4
   * - MP2, SCS/SOS, RI-MP2
     - yes (phase 3)
     - no: frozen virtuals would be correlated
     - yes, as now (SCS/SOS still refused)
     - phase 4
   * - CC
     - later
     - no
     - later
     - no
   * - MCSCF, xTB, SAPT, EFP
     - no
     - no
     - no
     - no

Unrestricted is refused everywhere, because ``run_czt_uhf`` has no embedding
or projector. A gradient is refused everywhere. Both gaps above are closed.

Phases and gates
================

0. **This page.** Approved before any code is written.
1. **The interface, with HF through it.** Phase 1a, done: the capability
   query and refusal site, and ``scf_numerics_from_config`` in place of the
   driver's four copies, bit-identical on the FMO, EFMO, NEO and MAKEFP decks.
   Phase 1b, done: ``mqc_czt_fragment_solver`` is the only place FMO, EE-MBE
   and EFMO solve a fragment or n-mer (the AFO model system and MAKEFP stay
   Hartree-Fock by construction), and ``run_czt_fmo``/``run_czt_efmo`` take
   the ``cuest_scf_settings_t`` that ``method_backend_settings`` builds.

   Gate: every existing FMO, EFMO, EE-MBE, AFO and PIEDA test passes
   unchanged, and totals are bit-identical at one thread on the FMO, EFMO,
   NEO and MAKEFP decks in ``validation/inputs/cpu/mqc/``, against the parent
   commit rebuilt in the same tree. Two builds of the *same* source can differ
   in the last bit of EFMO's far-pair EFP terms, because link-time
   optimisation's partitioning can inline them differently. A reference
   built at another time is therefore not a valid comparison.
2. **DFT under FMO**, including range-separated hybrids and D3/D4.

   Gate: two fragments at full level reproduce the supermolecule, the same
   identity the covalent EFMO tests use. An FMO2 water trimer matches an
   independent reference.
3. **The MP2 family** as correlation on embedded HF.

   Gate: the same full-level identity, and an FMO2 reference from GAMESS.
4. **PIEDA.**

   * For MP2, Edi is ``Ec(IJ) - Ec(I) - Ec(J)``.
   * For DFT, ``E'^HL`` takes ``E_XC[D_HL]``, which costs one extra
     quadrature per pair.
   * ``pieda_dispersion`` stays an option.
5. **Later:** CC, and anything else the capability query allows.

Decisions
=========

1. **ESP density for correlated methods: HF**, as GAMESS does. Every ESP,
   set of charges and ``Tr(ΔD u)`` comes from the embedded HF density.
   Correlation is added per fragment and per n-mer.
2. **D3/D4 under FMO-DFT: match GAMESS.** Whether that is once for the whole
   system or per fragment and n-mer is to be checked against GAMESS before
   phase 2. Phase 1 does not depend on it.
3. **EFMO with DFT stays refused.** MAKEFP is HF, so EFMO allows HF and the
   MP2 family only. It still routes through the solver.
