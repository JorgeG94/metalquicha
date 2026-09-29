========================================
The fragment solver: design and progress
========================================

**Status: design agreed (see Decisions at the end); phases 1 to 4 built, dispersion included.**
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

The AFO model system goes through the solver too:

* ``bond_lmo_set`` and ``bond_hybrid`` in ``mqc_czt_afo.f90`` solve the model
  system through ``solve_fragment_method`` at ``afo_options_t%method``, which
  ``build_afo_context`` fills from the deck's method with the correlation
  switched off. The model only supplies the frozen orbitals, so under Kohn-Sham
  they are Kohn-Sham orbitals at the deck's functional, and under MP2 they are
  the Hartree-Fock reference's. An unallocated ``method`` is Hartree-Fock and
  reaches ``run_czt_rhf`` with the arguments it always had.

``make_efp_potential`` in ``mqc_czt_efp_potential.f90`` stays Hartree-Fock by
construction and does not go through the solver. MAKEFP is a Hartree-Fock
construction.

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
  ``dispersion_apply``) for an unfragmented run. Under FMO and EE-MBE they are
  added in the fragment solver, per fragment and n-mer (Decision 2).
* ``run_czt_hf``, the unfragmented dispatcher, passes neither ``h_extra`` nor
  ``projector``. It returns no density.

The design
==========

Where it lives
--------------

``backends/cenzontle/fragments/mqc_czt_fragment_solver.f90`` (phases 1b and 2,
built).
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
      character(len=16) :: dispersion                    ! "none", "d3bj" or "d4"
      real(dp), allocatable :: dispersion_xyz(:, :)     ! the real atoms, Bohr
      real(dp) :: dispersion_charge                      ! only D4 reads it
   end type

   type :: fragment_outcome_t
      real(dp) :: energy        ! E: reference + correlation + dispersion, including Tr(D u)
      real(dp) :: internal      ! E' = E - Tr(D_esp u); energy when there is no u
      real(dp) :: reference     ! the SCF energy
      real(dp) :: correlation   ! zero for HF (and DFT)
      real(dp) :: dispersion    ! inside energy and internal; zero unless asked for
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
* **DFT:** built. The same call with an ``xc_context_t`` built from
  ``settings``, restricted only. ``fragment_refusal`` refuses a double
  hybrid, which it identifies by name through ``xc_spec_from_name``, so its
  PT2 part is never silently dropped. The solver has a guard for the same
  case. Dispersion, when the request asks for it, is added here (Decision 2).
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

``src/methods/dispatch/mqc_fragment_capabilities.f90`` (phases 1a and 2,
built). It is in ``src`` because the refusal must happen in the driver, before
any backend work:

.. code-block:: fortran

   type :: fragment_needs_t          ! what the deck asks
      logical :: cut                 ! bond_breaking = "afo"
      logical :: unrestricted        ! model.unrestricted
      integer :: calc_type
      logical :: pieda
      logical :: pieda_dispersion    ! keywords.fragmentation.pieda_dispersion
      logical :: dispersion          ! keywords.dft.dispersion
   end type

   type :: fragment_capabilities_t   ! what the method can do under the scheme
      logical :: runs, cut, unrestricted, gradient, pieda, pieda_dispersion, dispersion
   end type

   pure function fragment_capabilities(scheme, method_config) result(cap)
   function fragment_refusal(scheme, method_config, needs) result(message)

Every method that runs at all accepts an embedding operator, so embedding is
not a capability. ``pieda`` says whether the method has its own PIEDA terms.
Which schemes offer PIEDA is still decided by ``check_pieda_support``.
``pieda_dispersion`` says whether PIEDA's empirical ``Edi`` can be added beside
the method's terms: it can for Hartree-Fock and Kohn-Sham, and not for the MP2
family, whose ``Edi`` is its correlation. ``dispersion`` is true for Kohn-Sham
under FMO and EE-MBE and false for every other method and for EFMO -- see
Decision 2. ``fragment_refusal`` also refuses ``dispersion`` together with
``pieda_dispersion``, which would count the pair dispersion twice.

``fragment_refusal`` is the single refusal site. ``run_fragmented_calculation``
calls it for FMO and EE-MBE, and ``run_efmo_energy`` calls it for EFMO. It
replaced ``fmo_method_refusal`` and the method refusals in
``run_efmo_energy``. The backend's MP2-across-a-cut check in ``run_efmo`` is
now only a guard against a caller that skipped it. The message names the
method, the scheme and the key to change.

Refusals that depend on the data rather than the method stay where the data
is. An odd electron count after a cut is one example.

.. list-table:: Capabilities, by phase
   :header-rows: 1

   * - method
     - FMO / EE-MBE
     - with a cut (AFO)
     - EFMO
     - PIEDA
     - D3/D4
   * - HF
     - yes (phase 1)
     - yes
     - yes
     - yes
     - refused: no functional for the damping (Decision 2)
   * - DFT, including RSH
     - yes (phase 2)
     - yes (phase 2)
     - no: MAKEFP is HF
     - yes (phase 4), FMO only
     - yes, per fragment and n-mer (Decision 2)
   * - DFT, double hybrid
     - no: PT2 fraction not added
     - no
     - no: MAKEFP is HF
     - no
     - refused
   * - MP2, SCS/SOS, RI-MP2
     - yes (phase 3)
     - no: frozen virtuals would be correlated
     - yes, as before (SCS/SOS still refused there)
     - yes (phase 4), FMO only; Edi is the correlation interaction
     - n/a
   * - CC
     - later
     - no
     - later
     - no
     - n/a
   * - MCSCF, xTB, SAPT, EFP
     - no
     - no
     - no
     - no
     - n/a

Unrestricted is refused everywhere, because ``run_czt_uhf`` has no embedding
or projector. A gradient is refused everywhere. Both gaps above are closed.

Phases and gates
================

0. **This page.** Approved before any code is written.
1. **The interface, with HF through it.** Phase 1a, done: the capability
   query and refusal site, and ``scf_numerics_from_config`` in place of the
   driver's four copies, bit-identical on the FMO, EFMO, NEO and MAKEFP decks.
   Phase 1b, done: ``mqc_czt_fragment_solver`` is the only place FMO, EE-MBE
   and EFMO solve a fragment or n-mer (MAKEFP stays Hartree-Fock by
   construction; the AFO model system follows the deck's functional, see
   Decisions), and ``run_czt_fmo``/``run_czt_efmo`` take
   the ``cuest_scf_settings_t`` that ``method_backend_settings`` builds.

   Gate: every existing FMO, EFMO, EE-MBE, AFO and PIEDA test passes
   unchanged, and totals are bit-identical at one thread on the FMO, EFMO,
   NEO and MAKEFP decks in ``validation/inputs/cpu/mqc/``, against the parent
   commit rebuilt in the same tree. Two builds of the *same* source can differ
   in the last bit of EFMO's far-pair EFP terms, because link-time
   optimisation's partitioning can inline them differently. A reference
   built at another time is therefore not a valid comparison.
2. **DFT under FMO**, including range-separated hybrids. Built, restricted
   only, for FMO and EE-MBE, with or without a detached bond; D3/D4 is
   added per fragment and n-mer (Decision 2). A double hybrid is refused by
   ``fragment_refusal``.

   Gate: two fragments at full level reproduce the supermolecule, the same
   identity the covalent EFMO tests use -- ``test/test_mqc_fmo_dft.f90``, on
   a GGA, a global hybrid and a range-separated hybrid, under FMO and
   EE-MBE, with and without a detached bond, and a three-fragment FMO3 water
   trimer at full level.

   A full-level identity cannot see the embedding. When the fragment count
   equals the level, every embedded term cancels out of the total. What
   checks embedded Kohn-Sham is the Hellmann-Feynman test in
   ``test/test_mqc_czt_fragment_solver.f90``: with ``h + lambda u``, dE/dlambda
   equals Tr(D u) to about 1e-10 for HF, PBE and CAM-B3LYP.

   End to end through the driver, the FMO3 water-trimer deck with
   ``"dft"``/``"pbe"`` gives the unfragmented PBE energy to 5e-12 Eh.

   Still owed: an FMO2 comparison against an independent reference, which
   needs GAMESS.
3. **The MP2 family** as correlation on embedded HF. Built, for FMO and
   EE-MBE, without a detached bond (Decision 1: the ESP, the charges and
   ``Tr(dD u)`` stay the embedded Hartree-Fock's; correlation is added per
   fragment and per n-mer, as GAMESS FMO-MP2 does). Restricted MP2, SCS-MP2,
   SOS-MP2 and RI-MP2 all run; a double hybrid stays refused as
   before, and a separated pair (``resdim``) gets no pair-level correlation.

   The outer (monomer) self-consistent-charge loop still solves plain
   Hartree-Fock on every pass -- adding MP2 to it would be wasted work and
   would perturb the convergence test on the sum of monomer energies. One
   more monomer pass, with correlation, runs once the loop has settled, under
   the field it converged to and started from each monomer's own converged
   density (``SCF_GUESS_PROJ``), so it typically finishes in a few
   iterations. An n-mer's correlation is added by ``solve_fragment_method``
   in the same call as its reference, since an n-mer is solved once and not
   iterated.

   Gate: the same full-level identity as phase 2 (water dimer and trimer, FMO
   and EE-MBE, MP2, SCS-MP2 and RI-MP2, frozen core on and off) --
   ``test/test_mqc_fmo_mp2.f90``.

   The embedding is checked independently. ``tools/fmo_validation/eembe_pyscf.py``
   is a PySCF reimplementation of EE-MBE, fed this repository's basis JSON.
   On the water trimer at level 2 it reproduces the Hartree-Fock total to
   3e-11 Eh, and MP2, RI-MP2 and SCS-MP2 to 9e-9 Eh. The MP2 agreement is
   asserted in ``test_mqc_fmo_mp2``.

   Still owed: an FMO2-MP2 reference, which needs GAMESS. FMO's exact ESP is
   not in the PySCF replica.
4. **PIEDA.** Built, for FMO only (``check_pieda_support`` still refuses
   EE-MBE and EFMO).

   * **MP2 family.** As GAMESS's PIEDA/MP2: ``Ees``, ``Eex`` and ``Ect+mix``
     are the Hartree-Fock terms from the embedded Hartree-Fock densities and
     the monomers' Hartree-Fock internal energies (``fragment_t%energy -
     fragment_t%correlation``). ``Edi = Ec(IJ) - Ec(I) - Ec(J)`` is inside
     ``dE_IJ``, and ``Ect+mix = dE_IJ - Ees - Eex - Edi``. The pair's own
     correlation is carried out of ``nmer_term`` before ``subtract_subsets``.
     ``fmo_result_t%edi_in_energy`` (JSON ``edi_in_energy``) says that ``Edi``
     is in ``dE_IJ``, so the printed ``Total`` is not ``dE_IJ + Edi``.
   * **Kohn-Sham.** ``E'^HL`` is the Kohn-Sham functional at ``D_HL``,
     ``E_XC[D_HL]``, its exact-exchange fraction and its range separation,
     which costs one extra quadrature per pair. It is ``density_energy`` in
     ``mqc_czt_rhf``, the ``assemble_fock`` the SCF's own energy comes from,
     with the ``xc_context_t`` that ``fragment_xc_context`` builds on the
     pair's molecule from the same settings the solver uses. A detached bond
     enters only through the construction of ``D_HL`` (``pieda_hl``).
   * ``pieda_dispersion`` stays an option at Hartree-Fock and under Kohn-Sham,
     with the functional's damping parameters, and is refused with the MP2
     family (``fragment_refusal``, ``pieda_dispersion`` need): its ``Edi``
     already holds the dispersion.

   Gate: ``test/test_mqc_fmo_pieda_methods.f90``. The four terms close every
   pair energy for PBE, PBE0, CAM-B3LYP, MP2, SCS-MP2 and RI-MP2; an MP2 run's
   ``Ees`` and ``Eex`` equal Hartree-Fock's to 1e-10; the ``Edi`` of MP2, RI-MP2
   and SCS-MP2 equals separately computed correlation energies; the union
   energy at one fragment's density is that fragment's Kohn-Sham energy for
   PBE, B3LYP and CAM-B3LYP; Kohn-Sham ``Eex`` vanishes with distance and
   repels at contact; and PBE decomposes the pairs next to a detached bond in
   both ``pieda_hl`` modes.
5. **Later:** CC, and anything else the capability query allows.

Decisions
=========

1. **ESP density for correlated methods: HF**, as GAMESS does. Every ESP,
   set of charges and ``Tr(ΔD u)`` comes from the embedded HF density.
   Correlation is added per fragment and per n-mer.
2. **D3/D4 under FMO-DFT: match GAMESS.** GAMESS applies D3 per fragment
   and per n-mer (``dftdis.src``, ``DFTDSM``/``DFTDSMI``; confirmed from a
   run), not once for the whole system. **Done.** ``fmo_options_t%dispersion``
   (``"d3bj"`` or ``"d4"``, from ``keywords.dft.dispersion``) reaches every
   ``solve_fragment_method`` call through ``fragment_request_t``, which carries
   the group's real atoms (never a ghost centre or a split nucleus), their
   coordinates and the sum of the members' declared charges. The solver adds
   the correction once to ``outcome%energy`` and ``outcome%internal`` -- it does
   not depend on the density, so ``E'`` changes by it alone -- and reports it
   in ``outcome%dispersion``. A monomer's is recomputed on each pass of the
   monomer loop and replaces the last, so it is counted once. A separated pair
   runs no SCF, so ``calculate_polymers`` adds ``E_D(IJ) - E_D(I) - E_D(J)``
   to its term directly. PIEDA reports the same increment as ``edi``, inside
   the pair energy; ``Eex`` takes the monomers' dispersion back out of their
   ``E'``. Gate: ``test/test_mqc_fmo_dispersion.f90``, and the PBE-D3(BJ) case
   of ``check_fmo_mpi`` for the distributed path.
   Whether GAMESS also adds it to a separated dimer is still to be confirmed.
3. **EFMO with DFT stays refused.** MAKEFP is HF, so EFMO allows HF and the
   MP2 family only. It still routes through the solver.

4. **The AFO model system under Kohn-Sham: solve it at the deck's
   functional**, as GAMESS does. GAMESS solves a cut bond's model system with
   the deck's functional under FMO-DFT; its log shows ``FINAL R-PBE ENERGY`` for
   the model. Done: ``bond_lmo_set`` and ``bond_hybrid`` go through
   ``solve_fragment_method`` with the functional, grid and settings of the
   fragments, so cut FMO-DFT freezes Kohn-Sham-derived orbitals. The model's own
   convergence is unchanged. Under MP2 the model is the Hartree-Fock reference,
   because it supplies orbitals and no energy. Hartree-Fock results are
   bit-identical to before. Gate: ``test_mqc_afo_orbital`` (the PBE orbital set
   against an independent PBE solve, and unlike the Hartree-Fock one) and
   ``propane_in_three_fragments_pbe_freezes_a_kohn_sham_orbital`` in
   ``test_mqc_fmo_dft``. No GAMESS reference for cut FMO-DFT exists yet: its
   own PBE model SCF did not converge on Gly3 plus water in STO-3G.

Open, for review
================

* None at present.
