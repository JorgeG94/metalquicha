!! The effective fragment molecular orbital energy
module mqc_czt_efmo
   !! EFMO (Sattasathuchana, Xu, Bertoni, Kim, Leang, Pham, Gordon, JCTC 20,
   !! 2445 (2024), eq 6; Steinmann, Fedorov, Jensen, JPCA 114, 8705 (2010)):
   !!
   !!     E = sum over near groups S, |S| <= n, of ( dE_S^0 - dE_S^pol )
   !!       + sum_{I<J, R_IJ >  R_cut} ( E_IJ^Coul + E_IJ^disp + E_IJ^ExRep + E_IJ^CT )
   !!       + E_pol^total
   !!
   !! **EFMO to any many-body order.** `dE_S^0` is the many-body difference of
   !! the *in-vacuo* energies of `S` and its subsets and `dE_S^pol` is the same
   !! difference applied to the induction energy of the same group's potentials,
   !! both from [[mqc_czt_subsets]]. At `n = 2` that is eq 6 of the paper
   !! written out -- `dE_I^0 = E_I^0`, `dE_IJ^0 = E_IJ^0 - E_I^0 - E_J^0`,
   !! `dE_I^pol = 0` and `dE_IJ^pol = E_IJ^pol` -- and the level-two total is
   !! unchanged to the last bit. At `n = 3` the group terms are the usual
   !! three-body forms, `E_IJK^pol - E_IJ^pol - E_IK^pol - E_JK^pol` being the
   !! non-additive part of the trimer's induction.
   !!
   !! **Nothing in the expansion is two-body specific**, which is why this
   !! generalises at all: no group's Hamiltonian depends on its environment, so
   !! the many-body differences telescope exactly. GAMESS stops at two because
   !! it enumerates the levels by hand, not because the method does.
   !!
   !! **A group is near only if EVERY pair inside it is near.** The far half of
   !! the energy is pairwise by construction -- the four effective-fragment
   !! terms are two-body and the induction is already all orders in
   !! `E_pol^total` -- so there is no far n-body term, and a group holding a far
   !! pair would count that pair twice. `efmo_near_subsets` implements exactly
   !! that, and the criterion is inherited by subsets, which is what makes the
   !! filtered group list a valid one to difference over.
   !!
   !! **At level = N with `R_cut` huge the whole polarization correction
   !! cancels**: the induction series telescopes to `E_pol^total`, the last term
   !! of the energy, and the in-vacuo series telescopes to the supersystem's own
   !! energy. So EFMO at full level *is* the unfragmented calculation, which is
   !! the sharpest available check on the subset-induction bookkeeping and is
   !! asserted in `test_mqc_czt_efmo`.
   !!
   !! **Nothing here is self-consistent across fragments.** `E_I^0` and `E_IJ^0`
   !! are *in vacuo* energies -- no embedding field, no monomer loop -- which is
   !! what separates EFMO from FMO and what lets a diffuse basis work: there is
   !! no neighbouring point charge for a diffuse function to fall onto. The
   !! coupling between fragments is carried entirely by the effective fragment
   !! potentials, one per fragment, each built by MAKEFP from the same SCF that
   !! produced `E_I^0`.
   !!
   !! **`E_S^pol` is subtracted from every near group, and it is not small.**
   !! A quantum dimer already contains the mutual induction of its two
   !! fragments, and `E_pol^total` -- the induction solved over every fragment
   !! at once -- contains it too, so one copy has to go. What survives,
   !! `E_pol^total - sum_IJ E_IJ^pol`, is the many-body part of the induction:
   !! measured on three waters at four Angstrom it is 44 per cent of the total,
   !! because the energy is quadratic in the field and the square of a sum keeps
   !! cross terms no pair has. It is a term of the method, not a residue.
   !!
   !! **Where each number comes from**, all of it Phase 1 work:
   !!
   !! * `E_I^0` is `pot%scf_energy`, the SCF `make_efp_potential` runs anyway.
   !! * the potential becomes an `efp_fragment_t` through `potential_to_fragment`
   !!   rather than through a written `.efp` file.
   !! * the QM/EFP split is `efmo_split_pairs` on the vdW-scaled `R_IJ` of eq 2.
   !! * the far pairs go through `efp_pair_terms`, which builds a two-fragment
   !!   system per pair with the charge-penetration screening on. **The
   !!   system-wide `electrostatic_energy` is deliberately never called here**:
   !!   it would put every fragment's points in one set, which is the same sum
   !!   only when *every* pair is a far pair, and there is no mask to exclude
   !!   the near ones.
   !! * `E_pol^total` is `polarization_energy` over every fragment at once,
   !!   which is the call `efp_interaction_energy` makes internally.
   !!
   !! **Cartesian throughout.** `make_efp_potential` forces 6d/10f, because a
   !! `.efp` is read by GAMESS, so `E_I^0` is a Cartesian energy. The dimer SCFs
   !! here are built the same way: with a spherical dimer the difference
   !! `E_IJ^0 - E_I^0 - E_J^0` would be taken between two different models and
   !! would not be an interaction energy at all. Below d functions the two forms
   !! coincide and the choice is invisible.
   !!
   !! **No embedding, quantum or effective, anywhere.** Every group's SCF runs
   !! in vacuo: the neighbouring fragments' potentials are not in its
   !! Hamiltonian. That is deliberate and is what the whole expansion above
   !! rests on. Mutually polarized QM/EFP embedding -- the fragment SCF solved
   !! in the field of the other potentials, with the induction inside the SCF --
   !! is a **separate method for later**, not a refinement of this one: it makes
   !! a group's energy depend on its environment, so the many-body differences
   !! no longer telescope and the level = N identity above stops holding.
   !!
   !! **Closed shell, whole molecules.** Covalent fragments are Phase 5; a
   !! partition that cuts a bond is refused here rather than capped, since a
   !! cap's multipoles would act on the partner across the cut.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_calculation_defaults, only: DEFAULT_VDW_SCALE, DEFAULT_DYNAMIC_TOL, &
                                       DEFAULT_DYNAMIC_MAXITER, DEFAULT_RESPONSE_BATCH, &
                                       EFP_RESPONSE_AUTO
   use mqc_scf_types, only: scf_numerics_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_atomic_guess, only: build_restricted_guess
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_efp_potential, only: efp_potential_t, make_efp_potential
   use mqc_czt_efp_read, only: efp_fragment_t
   use mqc_czt_efp_convert, only: potential_to_fragment
   use mqc_czt_efp_energy, only: efp_pair_energy_t, efp_pair_terms, &
                                 subset_polarization_energy
   use mqc_czt_efp_interaction, only: efp_system_t, build_efp_system, polarization_energy
   use mqc_czt_efmo_pairs, only: efmo_split_pairs, efmo_near_subsets
   use mqc_czt_subsets, only: subtract_subsets, n_choose
   use mqc_czt_mp2, only: mp2_result_t, run_czt_mp2, run_czt_ri_mp2
   use mqc_elements, only: core_orbital_count
   use mqc_program_limits, only: EFMO_CORR_NONE, EFMO_CORR_MP2, EFMO_CORR_RI_MP2
   use mqc_czt_efp_serialize, only: EFP_HEADER_INTS, fragment_header, &
                                    fragment_buffer_sizes, fragment_pack, fragment_unpack
   use pic_mpi_lib, only: comm_t, allreduce, MPI_SUM
   use mqc_timing, only: timing_report_t
   implicit none
   private

   public :: EFMO_CORR_NONE, EFMO_CORR_MP2, EFMO_CORR_RI_MP2
   public :: efmo_options_t
   public :: efmo_pair_t
   public :: efmo_result_t
   public :: run_efmo

   type :: efmo_options_t
      !! What to run, and how hard
      character(len=64) :: basis = ""
         !! **Empty on purpose, and refused rather than defaulted.** This field
         !! used to start at "6-31g", which no run ever saw: every caller
         !! overwrites it from the deck, and a deck that omits `model.basis`
         !! gets "sto-3g" from `mqc_method_config`. So the initialiser named a
         !! basis nothing was ever computed in, which is worse than no default
         !! at all -- a plumbing bug that lost the deck's basis would have
         !! silently produced 6-31G numbers.
      real(dp) :: rcut = 2.0_dp
         !! `R_cut` of eq 2, **unitless**: each interatomic distance is divided
         !! by the two van der Waals radii, so 1 is contact. A pair at or inside
         !! it is a quantum dimer, a pair beyond it is four EFP terms. At or
         !! below zero every pair is effective, which is EFP with in-vacuo
         !! monomers; huge, every pair is quantum, which is FMO2 in vacuo plus
         !! the many-body induction. Both limits run, and both are tested.
      integer :: level = 2
         !! Truncate the many-body expansion of the near groups here.
         !! `keywords.fragmentation.level`, the same key MBE and FMO read.
         !!
         !! Two is EFMO as published and as GAMESS runs it. One is the fragment
         !! sum alone -- every near pair then contributes nothing, which is a
         !! legitimate truncation and a poor one. Above two the near groups are
         !! trimers and beyond, each an in-vacuo SCF and a subset induction; at
         !! the fragment count the whole expansion is exact, and with `rcut`
         !! huge it reproduces the unfragmented energy.
         !!
         !! **The cost is the binomial.** There are C(N, n) groups of size n
         !! before the near criterion thins them, so level three on twenty
         !! fragments is up to 1140 SCFs against 190 for level two. No level is
         !! refused; the count is reported before any of them is computed.
      logical :: charge_transfer = .true.
         !! Include `E_IJ^CT` in the far pairs. GAMESS's EFMO has it; the 2012
         !! method left it out, so it is switchable rather than assumed.
      real(dp) :: induction_damping = 0.0_dp
         !! Tang-Toennies-like damping of the induction field, `a` in
         !! `1 - exp(-a R^2)(1 + a R^2)`. **Zero is off**, which is what every
         !! reference pinned before Phase 4 was computed with; GAMESS's EFMO
         !! runs 0.6 for a cluster of whole molecules and 0.1 where a fragment
         !! was cut across a bond, and 0.6 is what closes our induction onto
         !! its numbers. Applied to `E_IJ^pol` and to `E_pol^total` alike --
         !! the two are the same solver on different systems and damping one
         !! without the other would leave the difference in the many-body
         !! remainder.
      integer :: correlation = EFMO_CORR_NONE
         !! `model.method`: `hf` leaves this alone, `mp2` and `ri-mp2` set it.
         !!
         !! **The whole of what a correlated EFMO is.** Eq 6 says nothing about
         !! the level of theory: `E_I^0` and `E_IJ^0` are in-vacuo energies of
         !! whatever model, so running MP2 after each of those SCFs turns the
         !! fragment sum and every near-dimer correction correlated and leaves
         !! the effective-fragment half exactly as it was. The far pairs and
         !! the induction come from the potentials, which are built from the
         !! Hartree-Fock density either way -- that is the method, not an
         !! approximation taken here: MAKEFP is a Hartree-Fock construction.
      character(len=64) :: corr_aux_basis = ""
         !! `model.aux_basis`, the fitting set `EFMO_CORR_RI_MP2` uses. Distinct
         !! from `aux_basis` above, which fits MAKEFP's response Hessian: one is
         !! a correlation-fitting set and the other a Coulomb-fitting one, and
         !! a run may want both or neither.
      logical :: freeze_core = .true.
      integer :: n_frozen_core = -1
         !! `keywords.correlation`. A negative count is derived per fragment
         !! from its elements, which makes the dimer's core the sum of its two
         !! monomers' -- so `E_IJ - E_I - E_J` differences the same set of
         !! correlated orbitals on both sides.
      character(len=32) :: guess = "auto"
         !! Initial guess for every SCF here, monomer and dimer alike.
      character(len=64) :: aux_basis = ""
         !! Fit the MAKEFP response Hessian against this basis. Empty is exact.
      type(scf_numerics_t) :: scf
         !! How every SCF is driven -- the accelerator, DIIS subspace, level
         !! shift, linear-dependence threshold and incremental Fock switch. Its
         !! `max_iter`, `energy_tol` and `density_tol` are not read: the four
         !! fields below are, and they are passed positionally so they win.
      integer :: scf_max_iter = 200
      real(dp) :: scf_energy_tol = 1.0e-10_dp
      real(dp) :: scf_density_tol = 1.0e-8_dp
      real(dp) :: scf_grad_tol = 1.0e-8_dp
         !! MAKEFP's own defaults, deliberately tighter than a whole-system
         !! run's. The dimer SCF has to be converged to the same place the
         !! monomer one is, because their difference is the interaction energy
         !! and is four orders smaller than either.
      real(dp) :: vdw_scale = DEFAULT_VDW_SCALE
      logical :: quadrupole_blocks = .true.
      real(dp) :: dynamic_tolerance = DEFAULT_DYNAMIC_TOL
      integer :: dynamic_maxiter = DEFAULT_DYNAMIC_MAXITER
      logical :: allow_crap_response = .false.
      integer :: response = EFP_RESPONSE_AUTO
      integer :: response_batch = DEFAULT_RESPONSE_BATCH
         !! `keywords.efp`, forwarded whole to `make_efp_potential`. Nothing
         !! here is read by this module; it configures the stages of MAKEFP
         !! after the SCF.
      logical :: verbose = .false.
         !! Let MAKEFP report its own stages. The EFMO table is written at info
         !! level either way.
   end type efmo_options_t

   type :: efmo_pair_t
      !! One fragment pair, in whichever of the two lists it landed
      integer :: i = 0, j = 0
      real(dp) :: r = 0.0_dp
         !! `R_IJ`, the vdW-scaled closest approach. Unitless.
      logical :: qm = .false.
         !! True: a dimer SCF ran and `e_dimer`/`e_pair_pol` are filled. False:
         !! the four EFP terms below are.
      real(dp) :: e_dimer = 0.0_dp             !! `E_IJ^0`, in vacuo
      real(dp) :: e_pair_pol = 0.0_dp          !! `E_IJ^pol`, subtracted
      real(dp) :: electrostatics = 0.0_dp
      real(dp) :: dispersion = 0.0_dp
      real(dp) :: exchange_repulsion = 0.0_dp
      real(dp) :: charge_transfer = 0.0_dp
   end type efmo_pair_t

   type :: efmo_result_t
      !! The total, and the six sums it is made of
      !!
      !! `energy` is exactly
      !! `monomer_sum + nmer_correction - induction_correction + far_electrostatics
      !! + far_dispersion + far_exchange_repulsion + far_charge_transfer
      !! + polarization_total`, with `induction_correction` held positive and
      !! subtracted, since that is how eq 6 writes it.
      real(dp) :: energy = 0.0_dp
      real(dp) :: monomer_sum = 0.0_dp
         !! `sum_I E_I^0`
      real(dp) :: nmer_correction = 0.0_dp
         !! `sum over near groups with |S| >= 2 of dE_S^0`. At level two that is
         !! `sum (E_IJ^0 - E_I^0 - E_J^0)` over the quantum dimers exactly.
      real(dp) :: induction_correction = 0.0_dp
         !! `sum over near groups of dE_S^pol`, **subtracted** from the total.
         !! At level two that is `sum E_IJ^pol` over the quantum dimers; above
         !! it, the pair sum plus the non-additive remainders of the larger
         !! groups. Reported with the sign it has as a sum, not the sign it
         !! enters with, so it can be compared against another code's directly.
      real(dp) :: far_electrostatics = 0.0_dp
      real(dp) :: far_dispersion = 0.0_dp
      real(dp) :: far_exchange_repulsion = 0.0_dp
      real(dp) :: far_charge_transfer = 0.0_dp
         !! The four EFP terms, summed over the effective dimers
      real(dp) :: polarization_total = 0.0_dp
         !! `E_pol^total`, induction over every fragment at once
      real(dp), allocatable :: level_vacuum(:)
         !! `sum over |S| = m of dE_S^0`, for m = 1 to the level run. Slot one
         !! is `monomer_sum`; the rest add up to `nmer_correction`. Per level
         !! rather than one number because the whole question a level-three run
         !! answers is how much the three-body term was worth.
      real(dp), allocatable :: level_induction(:)
         !! `sum over |S| = m of dE_S^pol`. Slot one is zero -- one fragment
         !! induces nothing -- and the rest add up to `induction_correction`.
      integer, allocatable :: level_count(:)
         !! How many near groups of each size there were
      real(dp) :: monomer_correlation = 0.0_dp
      real(dp) :: nmer_correlation = 0.0_dp
         !! How much of `monomer_sum` and `nmer_correction` above is
         !! correlation rather than Hartree-Fock. Reported, not summed: the two
         !! sums already hold them, because a correlated `E_I^0` *is* the
         !! monomer energy of eq 6. Zero on a Hartree-Fock run.
      real(dp), allocatable :: monomer_energy(:)   !! `E_I^0`
      type(efmo_pair_t), allocatable :: pairs(:)
         !! Every pair, quantum ones first, each carrying its `R_IJ`
      integer :: n_qm_pairs = 0
      integer :: n_efp_pairs = 0
      integer :: n_qm_groups = 0
         !! Near groups of two or more fragments -- the SCFs the near half
         !! cost. Equal to `n_qm_pairs` at level two.
   end type efmo_result_t

contains

   subroutine run_efmo(atomic_numbers, symbols, coordinates, owner, fragment_charges, &
                       opts, res, error, comm)
      !! One EFMO energy, from a system already partitioned into fragments
      !!
      !! `owner(i)` is the fragment of atom `i`, numbered from one with no gaps
      !! -- the same partition `run_fmo2` takes. Coordinates are Bohr.
      !! `fragment_charges(k)` is fragment `k`'s net charge; a charged fragment
      !! needs nothing extra, its monomer SCF and every dimer holding it just
      !! carry the charge.
      integer, intent(in) :: atomic_numbers(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: coordinates(:, :)     !! (3, n_atoms), Bohr
      integer, intent(in) :: owner(:)
      integer, intent(in) :: fragment_charges(:)    !! (n_fragments)
      type(efmo_options_t), intent(in) :: opts
      type(efmo_result_t), intent(out) :: res
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm
         !! Spread the monomers and the quantum groups over these ranks. Every
         !! rank runs this same routine on the same geometry and comes out with
         !! the same total; nothing is gathered to a leader.
         !!
         !! **What is distributed is what costs.** A monomer is a MAKEFP -- an
         !! SCF, a localization and twelve frequency-dependent response solves
         !! -- against one SCF for a group, so the monomer loop is what the
         !! balance is struck on, exactly as `run_fmo2` balances on its
         !! fragments. The far pairs and the induction stay replicated: they are
         !! milliseconds beside a potential, and replicating them means every
         !! rank reaches the same total without a second reduction.

      type(efp_fragment_t), allocatable :: frags(:)
      type(efp_pair_energy_t), allocatable :: far(:)
      real(dp), allocatable :: shifts(:, :)
      integer, allocatable :: qm_pairs(:, :), efp_pairs(:, :)
      integer, allocatable :: terms(:, :), term_size(:)
      real(dp), allocatable :: separation(:, :)
      integer, allocatable :: count_of(:)
      real(dp), allocatable :: monomer_corr(:)
      integer :: n_terms
      type(timing_report_t) :: clk
         !! Where an EFMO run's wall time goes, stage by stage. The paper's Fig
         !! 7 makes the same split and the claim it supports -- that the
         !! effective-fragment half is free beside the quantum half -- is one a
         !! run should be able to check on its own system rather than take.
      integer :: n_atoms, n_frag, k, p

      n_atoms = size(atomic_numbers)
      if (size(owner) /= n_atoms .or. size(coordinates, 2) /= n_atoms &
          .or. size(symbols) /= n_atoms) then
         call error%set(ERROR_VALIDATION, "efmo: the owner list, the symbols and the "// &
                        "coordinates must cover every atom")
         return
      end if
      if (n_atoms < 1) then
         call error%set(ERROR_VALIDATION, "efmo: there are no atoms to fragment")
         return
      end if
      if (len_trim(opts%basis) == 0) then
         call error%set(ERROR_VALIDATION, "efmo: no orbital basis was named. Every "// &
                        "caller sets it from the deck, so an empty one is a plumbing "// &
                        "fault rather than a request for a default -- and guessing a "// &
                        "basis here would return plausible numbers for a basis "// &
                        "nobody asked for.")
         return
      end if
      if (minval(owner) < 1) then
         call error%set(ERROR_VALIDATION, "efmo: every atom must belong to a fragment "// &
                        "numbered from one")
         return
      end if
      n_frag = maxval(owner)
      if (size(fragment_charges) /= n_frag) then
         call error%set(ERROR_VALIDATION, "efmo: the system has "//to_char(n_frag)// &
                        " fragments but "//to_char(size(fragment_charges))//" charges")
         return
      end if

      call fragment_counts(owner, n_frag, count_of, error)
      if (error%has_error()) return

      allocate (res%monomer_energy(n_frag), source=0.0_dp)
      allocate (monomer_corr(n_frag), source=0.0_dp)
      allocate (frags(n_frag), shifts(3, n_frag))
      ! Every potential is built at the geometry it is used at, so no fragment
      ! is placed and no rigid transform is looked for. That is why
      ! `place_fragment` never appears here.
      shifts = 0.0_dp

      call clk%start()
      call clk%begin("monomers (MAKEFP)")
      call build_potentials(atomic_numbers, symbols, coordinates, owner, count_of, &
                            fragment_charges, opts, frags, res%monomer_energy, &
                            monomer_corr, error, comm)
      call clk%lap("monomers (MAKEFP)")
      if (error%has_error()) return
      res%monomer_sum = sum(res%monomer_energy)
      res%monomer_correlation = sum(monomer_corr)

      ! One matrix of separations decides both halves: which pairs are far, and
      ! which groups are near enough to be solved quantum mechanically.
      call efmo_split_pairs(owner, atomic_numbers, coordinates, opts%rcut, &
                            qm_pairs, efp_pairs, error, r=separation)
      if (error%has_error()) return
      res%n_qm_pairs = size(qm_pairs, 2)
      res%n_efp_pairs = size(efp_pairs, 2)
      allocate (res%pairs(res%n_qm_pairs + res%n_efp_pairs))

      call efmo_near_subsets(separation, opts%rcut, opts%level, terms, term_size, &
                             n_terms, error)
      if (error%has_error()) return
      res%n_qm_groups = count(term_size(1:n_terms) >= 2)
      if (is_leader(comm)) call announce_cost(res, opts, n_frag)

      call clk%begin("quantum groups")
      call quantum_subsets(atomic_numbers, symbols, coordinates, owner, &
                           fragment_charges, terms, term_size, n_terms, separation, &
                           frags, shifts, opts, monomer_corr, res, error, comm)
      call clk%lap("quantum groups")
      if (error%has_error()) return

      ! The far half, and it stays pairwise however high the level goes: the
      ! effective-fragment electrostatics, exchange repulsion, dispersion and
      ! charge transfer are all two-body terms and the induction is already all
      ! orders in `E_pol^total`, so there is no far n-body term to add.
      ! `efp_pair_terms` takes the pair list directly, so the near pairs
      ! contribute nothing here -- which is the point, since their
      ! electrostatics, exchange and dispersion are inside their group's SCF.
      call clk%begin("effective-fragment pairs")
      far = efp_pair_terms(frags, shifts, efp_pairs, error, &
                           charge_transfer_on=opts%charge_transfer)
      if (error%has_error()) return
      do k = 1, res%n_efp_pairs
         p = res%n_qm_pairs + k
         res%pairs(p)%i = efp_pairs(1, k)
         res%pairs(p)%j = efp_pairs(2, k)
         res%pairs(p)%qm = .false.
         res%pairs(p)%r = separation(efp_pairs(1, k), efp_pairs(2, k))
         res%pairs(p)%electrostatics = far(k)%electrostatics
         res%pairs(p)%dispersion = far(k)%dispersion
         res%pairs(p)%exchange_repulsion = far(k)%exchange_repulsion
         res%pairs(p)%charge_transfer = far(k)%charge_transfer
      end do
      res%far_electrostatics = sum(far%electrostatics)
      res%far_dispersion = sum(far%dispersion)
      res%far_exchange_repulsion = sum(far%exchange_repulsion)
      res%far_charge_transfer = sum(far%charge_transfer)
      call clk%lap("effective-fragment pairs")

      call clk%begin("induction over all fragments")
      call total_polarization(frags, shifts, opts%induction_damping, &
                              res%polarization_total, error)
      call clk%lap("induction over all fragments")
      call clk%finish()
      if (error%has_error()) return

      res%energy = res%monomer_sum + res%nmer_correction - res%induction_correction &
                   + res%far_electrostatics + res%far_dispersion &
                   + res%far_exchange_repulsion + res%far_charge_transfer &
                   + res%polarization_total

      if (is_leader(comm)) then
         call report(res, opts)
         call clk%report("EFMO")
      end if
   end subroutine run_efmo

   function spread_over(comm) result(many)
      !! Whether there is more than one rank to spread over
      type(comm_t), intent(in), optional :: comm
      logical :: many

      many = .false.
      if (.not. present(comm)) return
      many = comm%size() > 1
   end function spread_over

   function mine(task, comm) result(owned)
      !! Whether this rank owns a task, round robin
      !!
      !! The same rule `run_fmo2` uses. Round robin rather than blocked because
      !! the tasks here are near enough equal in cost -- one potential is one
      !! potential -- and a contiguous split would leave the last rank short
      !! whenever the count is not a multiple of the rank count.
      integer, intent(in) :: task
      type(comm_t), intent(in), optional :: comm
      logical :: owned

      owned = .true.
      if (.not. spread_over(comm)) return
      owned = mod(task - 1, comm%size()) == comm%rank()
   end function mine

   function is_leader(comm) result(leads)
      !! Whether this rank writes the run's shared log lines
      type(comm_t), intent(in), optional :: comm
      logical :: leads

      leads = .true.
      if (.not. spread_over(comm)) return
      leads = comm%rank() == 0
   end function is_leader

   subroutine share_failure(error, comm)
      !! Make an error on any rank an error on every rank
      !!
      !! Without this the first rank to fail returns while the others block in
      !! the next reduction, and the run hangs rather than stops. The message
      !! is not shared -- it exists only where it was raised -- so a rank that
      !! was fine is told which of its neighbours was not and to look there.
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm

      integer :: status(1)

      if (.not. spread_over(comm)) return
      status = 0
      if (error%has_error()) status = 1
      call allreduce(comm, status, 1, MPI_SUM)
      if (status(1) == 0 .or. error%has_error()) return
      call error%set(ERROR_VALIDATION, "efmo: "//to_char(status(1))//" of "// &
                     to_char(comm%size())//" ranks failed on the fragments they "// &
                     "own. The message is on those ranks; this one had no error "// &
                     "of its own and stops so that the run does not hang in the "// &
                     "next reduction.")
   end subroutine share_failure

   subroutine exchange_potentials(frags, opts, error, comm)
      !! Give every rank every fragment's potential
      !!
      !! A rank built a third of the potentials and needs all of them: the far
      !! pairs and the one induction over every fragment are not decomposable
      !! by owner. So each potential is flattened into an integer and a real
      !! buffer that is zero where this rank computed nothing, the buffers are
      !! summed across ranks, and the fragments this rank does not own are
      !! rebuilt from the sum. Bit for bit what the owning rank had, which is
      !! what makes one rank and four agree to the last digit rather than to
      !! the eight decimals a written `.efp` would carry.
      !!
      !! Two reductions, because the *layout* has to be known before the
      !! contents can be: a fixed-size header of counts and flags first, then a
      !! body whose length every rank works out from the reduced headers.
      type(efp_fragment_t), intent(inout) :: frags(:)
      type(efmo_options_t), intent(in) :: opts
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm

      integer, allocatable :: headers(:), ibuf(:)
      integer, allocatable :: ni(:), nr(:), ioff(:), roff(:)
      real(dp), allocatable :: rbuf(:)
      integer :: n_frag, k, total_i, total_r

      if (.not. spread_over(comm)) return
      n_frag = size(frags)

      ! One flat vector rather than a matrix: `allreduce` takes rank-one
      ! buffers, and a header is a fixed stride, so fragment `k`'s slice is
      ! arithmetic.
      allocate (headers(EFP_HEADER_INTS*n_frag), source=0)
      do k = 1, n_frag
         if (mine(k, comm)) then
            call fragment_header(frags(k), headers(head_at(k):head_at(k) + EFP_HEADER_INTS - 1))
         end if
      end do
      call allreduce(comm, headers, size(headers), MPI_SUM)

      allocate (ni(n_frag), nr(n_frag), ioff(n_frag), roff(n_frag))
      total_i = 0
      total_r = 0
      do k = 1, n_frag
         call fragment_buffer_sizes(headers(head_at(k):head_at(k) + EFP_HEADER_INTS - 1), &
                                    ni(k), nr(k))
         ioff(k) = total_i
         roff(k) = total_r
         total_i = total_i + ni(k)
         total_r = total_r + nr(k)
      end do
      if (total_r == 0) return

      allocate (ibuf(max(total_i, 1)), source=0)
      allocate (rbuf(max(total_r, 1)), source=0.0_dp)
      do k = 1, n_frag
         if (.not. mine(k, comm)) cycle
         call fragment_pack(frags(k), headers(head_at(k):head_at(k) + EFP_HEADER_INTS - 1), &
                            ibuf(ioff(k) + 1:ioff(k) + ni(k)), &
                            rbuf(roff(k) + 1:roff(k) + nr(k)), error)
         if (error%has_error()) exit
      end do
      call share_failure(error, comm)
      if (error%has_error()) return

      if (total_i > 0) call allreduce(comm, ibuf, total_i, MPI_SUM)
      call allreduce(comm, rbuf, total_r, MPI_SUM)

      do k = 1, n_frag
         if (mine(k, comm)) cycle
         call fragment_unpack(headers(head_at(k):head_at(k) + EFP_HEADER_INTS - 1), &
                              ibuf(ioff(k) + 1:ioff(k) + ni(k)), &
                              rbuf(roff(k) + 1:roff(k) + nr(k)), frags(k), error)
         if (error%has_error()) return
      end do
      if (opts%verbose .and. is_leader(comm)) then
         call logger%verbose("  efmo: shared "//to_char(n_frag)//" potentials over "// &
                             to_char(comm%size())//" ranks, "//to_char(total_r)// &
                             " reals")
      end if
   end subroutine exchange_potentials

   pure function head_at(k) result(at)
      !! Where fragment `k`'s header starts in the flat header vector
      integer, intent(in) :: k
      integer :: at

      at = (k - 1)*EFP_HEADER_INTS + 1
   end function head_at

   subroutine fragment_counts(owner, n_frag, count_of, error)
      !! How many atoms each fragment holds
      !!
      !! The atoms of a fragment need not be contiguous in the deck, so there is
      !! no offset to keep and `gather` builds the index list instead. What is
      !! checked here is that the numbering has no gap, which would otherwise
      !! reach `make_efp_potential` as a molecule with no atoms.
      integer, intent(in) :: owner(:)
      integer, intent(in) :: n_frag
      integer, allocatable, intent(out) :: count_of(:)
      type(error_t), intent(inout) :: error

      integer :: i

      allocate (count_of(n_frag))
      count_of = 0
      do i = 1, size(owner)
         count_of(owner(i)) = count_of(owner(i)) + 1
      end do
      if (any(count_of == 0)) then
         call error%set(ERROR_VALIDATION, "efmo: fragment numbering has a gap -- some "// &
                        "fragment between one and "//to_char(n_frag)//" holds no atoms")
      end if
   end subroutine fragment_counts

   pure function gather(owner, k) result(idx)
      !! The system indices of fragment `k`'s atoms, in deck order
      integer, intent(in) :: owner(:)
      integer, intent(in) :: k
      integer, allocatable :: idx(:)

      integer :: i, n

      n = count(owner == k)
      allocate (idx(n))
      n = 0
      do i = 1, size(owner)
         if (owner(i) == k) then
            n = n + 1
            idx(n) = i
         end if
      end do
   end function gather

   subroutine build_potentials(z, symbols, xyz, owner, count_of, charges, opts, &
                               frags, monomer_energy, correlation, error, comm)
      !! One MAKEFP per fragment, and `E_I^0` off the same SCF
      !!
      !! **The whole cost of an EFMO run is here.** A potential is an SCF, a
      !! localization and twelve frequency-dependent response solves, against
      !! one SCF for a dimer, so the monomer loop dominates and is what Phase 4
      !! distributes.
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: owner(:), count_of(:), charges(:)
      type(efmo_options_t), intent(in) :: opts
      type(efp_fragment_t), intent(out) :: frags(:)
      real(dp), intent(out) :: monomer_energy(:)
      real(dp), intent(out) :: correlation(:)
         !! How much of each `monomer_energy` is correlation. Reported, and
         !! differenced out of each dimer's so the dimer correction can be
         !! reported the same way.
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm

      type(efp_potential_t) :: pot
      type(rhf_result_t) :: scf
      integer, allocatable :: idx(:)
      character(len=:), allocatable :: aux
      real(dp) :: e_corr
      integer :: k

      correlation = 0.0_dp
      monomer_energy = 0.0_dp
      do k = 1, size(count_of)
         ! Round robin, so a rank does every `size()`th potential. What is left
         ! at zero here is filled by the reduction below rather than left out.
         if (.not. mine(k, comm)) cycle
         idx = gather(owner, k)
         call logger%verbose("  efmo: fragment "//to_char(k)//" of "// &
                             to_char(size(count_of))//", "//to_char(size(idx))//" atoms")
         ! One call whether or not an auxiliary basis was named: an absent
         ! optional passed on as an actual argument arrives absent.
         if (len_trim(opts%aux_basis) > 0) then
            aux = trim(opts%aux_basis)
            call make_efp_potential(z(idx), symbols(idx), xyz(:, idx), trim(opts%basis), &
                                    "FRAG"//to_char(k), pot, error, charge=charges(k), &
                                    verbose=opts%verbose, aux_basis=aux, &
                                    guess=trim(opts%guess), &
                                    energy_tol=opts%scf_energy_tol, &
                                    density_tol=opts%scf_density_tol, &
                                    grad_tol_in=opts%scf_grad_tol, scf_in=opts%scf, &
                                    max_iter_in=opts%scf_max_iter, &
                                    vdwscl=opts%vdw_scale, &
                                    quadrupole_blocks=opts%quadrupole_blocks, &
                                    dynamic_tol=opts%dynamic_tolerance, &
                                    dynamic_maxiter=opts%dynamic_maxiter, &
                                    response=opts%response, &
                                    allow_crap_response=opts%allow_crap_response, &
                                    response_batch=opts%response_batch, scf_out=scf)
         else
            call make_efp_potential(z(idx), symbols(idx), xyz(:, idx), trim(opts%basis), &
                                    "FRAG"//to_char(k), pot, error, charge=charges(k), &
                                    verbose=opts%verbose, guess=trim(opts%guess), &
                                    energy_tol=opts%scf_energy_tol, &
                                    density_tol=opts%scf_density_tol, &
                                    grad_tol_in=opts%scf_grad_tol, scf_in=opts%scf, &
                                    max_iter_in=opts%scf_max_iter, &
                                    vdwscl=opts%vdw_scale, &
                                    quadrupole_blocks=opts%quadrupole_blocks, &
                                    dynamic_tol=opts%dynamic_tolerance, &
                                    dynamic_maxiter=opts%dynamic_maxiter, &
                                    response=opts%response, &
                                    allow_crap_response=opts%allow_crap_response, &
                                    response_batch=opts%response_batch, scf_out=scf)
         end if
         if (error%has_error()) return

         ! `E_I^0` *is* that SCF. Running a second one here would be a second
         ! determinant nothing compares against.
         monomer_energy(k) = pot%scf_energy
         call potential_to_fragment(pot, frags(k), error)
         call pot%destroy()
         if (error%has_error()) return

         ! And the correlation on those same orbitals, when the deck asked for
         ! it. **On the SCF MAKEFP already ran**, which is why `scf_out` exists:
         ! a second SCF here would converge to the same place and cost as much
         ! as everything after it.
         if (opts%correlation /= EFMO_CORR_NONE) then
            call fragment_correlation(z(idx), symbols(idx), xyz(:, idx), &
                                      sum(z(idx)) - charges(k), scf, opts, &
                                      e_corr, error)
            if (error%has_error()) return
            monomer_energy(k) = monomer_energy(k) + e_corr
            correlation(k) = e_corr
         end if
      end do
      call share_failure(error, comm)
      if (error%has_error()) return

      if (spread_over(comm)) then
         call allreduce(comm, monomer_energy, size(monomer_energy), MPI_SUM)
         call allreduce(comm, correlation, size(correlation), MPI_SUM)
      end if
      call exchange_potentials(frags, opts, error, comm)
   end subroutine build_potentials

   subroutine fragment_correlation(z, symbols, xyz, nelec, scf, opts, energy, error)
      !! One fragment's MP2 correlation energy, on orbitals already converged
      !!
      !! The molecule is rebuilt rather than passed in, because the two callers
      !! get theirs from different places -- `make_efp_potential` builds and
      !! destroys its own -- and rebuilding it is a basis-set lookup against an
      !! SCF. **Cartesian**, the same as every SCF in this module and for the
      !! same reason: these orbitals came from a Cartesian molecule and an MO
      !! transform against a spherical one would be nonsense rather than a
      !! small error.
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: nelec
      type(rhf_result_t), intent(in) :: scf
      type(efmo_options_t), intent(in) :: opts
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: mol, aux
      type(mp2_result_t) :: mp2
      integer :: frozen

      energy = 0.0_dp
      if (opts%correlation == EFMO_CORR_NONE) return

      frozen = opts%n_frozen_core
      if (frozen < 0) frozen = core_orbital_count(z)
      if (.not. opts%freeze_core) frozen = 0

      call build_czt_molecule(z, symbols, xyz, trim(opts%basis), mol, error, &
                              force_cartesian=.true.)
      if (error%has_error()) return

      if (opts%correlation == EFMO_CORR_RI_MP2) then
         if (len_trim(opts%corr_aux_basis) == 0) then
            call error%set(ERROR_VALIDATION, "efmo: a fitted correlation needs an "// &
                           "auxiliary basis. Set model.aux_basis, or ask for 'mp2' "// &
                           "rather than 'ri-mp2'.")
            call mol%destroy()
            return
         end if
         ! The fitting set in the orbital basis's angular form, as everywhere
         ! else here: libcint builds all three centres of a fitting integral in
         ! one form, and the orbital basis is Cartesian because MAKEFP's is.
         call build_czt_molecule(z, symbols, xyz, trim(opts%corr_aux_basis), aux, &
                                 error, force_cartesian=.true.)
         if (error%has_error()) then
            call mol%destroy()
            return
         end if
         call run_czt_ri_mp2(mol, aux, scf%orbitals, scf%orbital_energies, &
                             nelec/2, scf%energy, mp2, error, n_frozen=frozen)
         call aux%destroy()
      else
         call run_czt_mp2(mol, scf%orbitals, scf%orbital_energies, &
                          nelec/2, scf%energy, mp2, error, n_frozen=frozen)
      end if
      call mol%destroy()
      if (error%has_error()) return

      ! `same_spin + opposite_spin`, unscaled. Spin-component scaling is not
      ! offered: SCS-MP2 fragment energies would be a different method and the
      ! deck names it separately, so a run that asked for it is refused above
      ! rather than silently given plain MP2.
      energy = mp2%same_spin + mp2%opposite_spin
   end subroutine fragment_correlation

   subroutine announce_cost(res, opts, n_frag)
      !! What the near half is about to cost, before any of it is paid
      !!
      !! **The binomial is the whole story**, so it is said out loud rather than
      !! left in a docstring: `C(N, n)` groups of size `n` before the near
      !! criterion thins them, each one an in-vacuo SCF. A level a user picked
      !! without doing that arithmetic is better met with a warning than with a
      !! refusal -- the run may be exactly what was wanted -- and better with a
      !! warning than with silence, because the cost is superlinear in a number
      !! typed as a single digit.
      type(efmo_result_t), intent(in) :: res
      type(efmo_options_t), intent(in) :: opts
      integer, intent(in) :: n_frag

      integer :: level, m, unscreened

      level = min(opts%level, n_frag)
      call logger%info("  efmo: "//to_char(res%n_qm_groups)//" quantum groups up to "// &
                       "level "//to_char(level)//", "//to_char(res%n_efp_pairs)// &
                       " effective-fragment pairs")
      if (level < 3) return

      unscreened = 0
      do m = 2, level
         unscreened = unscreened + n_choose(n_frag, m)
      end do
      call logger%warning("efmo at level "//to_char(level)//" on "//to_char(n_frag)// &
                          " fragments enumerates up to "//to_char(unscreened)// &
                          " groups, of which the cutoff kept "// &
                          to_char(res%n_qm_groups)//". The count is the binomial "// &
                          "C(N, n) and each group is an SCF, so a level raised by "// &
                          "one is not a small change.")
   end subroutine announce_cost

   subroutine quantum_subsets(z, symbols, xyz, owner, charges, terms, term_size, &
                              n_terms, separation, frags, shifts, opts, monomer_corr, &
                              res, error, comm)
      !! Every near group: one in-vacuo SCF, its subset induction, and the
      !! many-body difference of both
      !!
      !! **Two expansions, one operator.** `dE_S^0` and `dE_S^pol` are the same
      !! difference from [[mqc_czt_subsets]] applied to two different quantities
      !! -- the group's in-vacuo energy and the induction energy of the group's
      !! own potentials -- and the near sum of the EFMO energy is the first minus
      !! the second. Applying it to the induction as well is the whole of the
      !! generalization: at level two `dE_IJ^pol` is `E_IJ^pol` and this reduces
      !! to eq 6 bit for bit, and at level N with a huge cutoff the induction
      !! series telescopes to `E_pol^total` and cancels it.
      !!
      !! Groups of one are not computed here: their energy is the monomer energy
      !! `build_potentials` already has and their induction is zero, and they are
      !! in the list only because the difference operator needs them to be.
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: owner(:), charges(:)
      integer, intent(in) :: terms(:, :), term_size(:)
         !! The near groups, smallest first, from `efmo_near_subsets`
      integer, intent(in) :: n_terms
      real(dp), intent(in) :: separation(:, :)
         !! `R_IJ` for every pair, for the reported pair table
      type(efp_fragment_t), intent(in) :: frags(:)
      real(dp), intent(in) :: shifts(:, :)
      type(efmo_options_t), intent(in) :: opts
      real(dp), intent(in) :: monomer_corr(:)
         !! Each `E_I^0`'s correlation part, differenced the same way so the
         !! correlation inside the group corrections can be reported beside them.
      type(efmo_result_t), intent(inout) :: res
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm

      integer, allocatable :: idx(:), members(:)
      real(dp), allocatable :: raw_vac(:), raw_pol(:), raw_corr(:)
      real(dp), allocatable :: group_vac(:), group_pol(:), group_corr(:)
      integer :: t, m, a, task, n_groups, level, pair

      n_groups = count(term_size(1:n_terms) >= 2)
      level = 1
      if (n_terms > 0) level = maxval(term_size(1:n_terms))
      allocate (group_vac(max(n_groups, 1)), group_pol(max(n_groups, 1)), &
                group_corr(max(n_groups, 1)), source=0.0_dp)

      task = 0
      do t = 1, n_terms
         m = term_size(t)
         if (m < 2) cycle
         task = task + 1

         ! Round robin over the flat group list, which at level two is the pair
         ! list in the order it was built -- so a level-two run distributes
         ! exactly as it did before there were larger groups to distribute.
         if (.not. mine(task, comm)) cycle
         members = terms(1:m, t)
         idx = group_atoms(owner, members)
         call logger%verbose("  efmo: group "//members_text(members)//", "// &
                             to_char(size(idx))//" atoms")
         call nmer_energy(z(idx), symbols(idx), xyz(:, idx), sum(charges(members)), &
                          opts, group_vac(task), group_corr(task), error)
         if (error%has_error()) exit

         ! The same induction solver on the group's fragments alone -- the same
         ! screening, the same static field rank, the same tolerance as the
         ! total below. A group solved any other way would leave a residue in
         ! the induction expansion that looks like non-additive induction.
         group_pol(task) = subset_polarization_energy(frags, shifts, members, error, &
                                                      damping=opts%induction_damping)
         if (error%has_error()) exit
      end do
      call share_failure(error, comm)
      if (error%has_error()) return

      ! Each rank filled only its own groups and left the rest at zero, so a sum
      ! gathers them. Summed *before* they are differenced and accumulated, not
      ! after, so the order the terms are added in does not depend on the rank
      ! count and one rank against four is an identity rather than a rounding
      ! question.
      if (spread_over(comm)) then
         call allreduce(comm, group_vac, size(group_vac), MPI_SUM)
         call allreduce(comm, group_pol, size(group_pol), MPI_SUM)
         call allreduce(comm, group_corr, size(group_corr), MPI_SUM)
      end if

      ! The three quantities on every group, singles included: the monomer
      ! energies the potentials came with, no induction at all, and the monomer
      ! correlation. Differencing needs them in the same list as the groups.
      allocate (raw_vac(max(n_terms, 1)), raw_pol(max(n_terms, 1)), &
                raw_corr(max(n_terms, 1)), source=0.0_dp)
      task = 0
      pair = 0
      do t = 1, n_terms
         m = term_size(t)
         if (m == 1) then
            a = terms(1, t)
            raw_vac(t) = res%monomer_energy(a)
            raw_corr(t) = monomer_corr(a)
            cycle
         end if
         task = task + 1
         raw_vac(t) = group_vac(task)
         raw_pol(t) = group_pol(task)
         raw_corr(t) = group_corr(task)
         if (m /= 2) cycle
         ! The pair table, filled on every rank from reduced numbers so it is
         ! identical everywhere. Pairs come out of the enumeration in the order
         ! `efmo_split_pairs` built its quantum list in, which is the order
         ! `res%pairs` reserves its first slots in.
         pair = pair + 1
         res%pairs(pair)%i = terms(1, t)
         res%pairs(pair)%j = terms(2, t)
         res%pairs(pair)%qm = .true.
         res%pairs(pair)%r = separation(terms(1, t), terms(2, t))
         res%pairs(pair)%e_dimer = group_vac(task)
         res%pairs(pair)%e_pair_pol = group_pol(task)
      end do

      call subtract_subsets(terms, term_size, n_terms, raw_vac)
      call subtract_subsets(terms, term_size, n_terms, raw_pol)
      call subtract_subsets(terms, term_size, n_terms, raw_corr)

      allocate (res%level_vacuum(level), res%level_induction(level), source=0.0_dp)
      allocate (res%level_count(level), source=0)
      do t = 1, n_terms
         m = term_size(t)
         res%level_count(m) = res%level_count(m) + 1
         if (m == 1) cycle
         res%level_vacuum(m) = res%level_vacuum(m) + raw_vac(t)
         res%level_induction(m) = res%level_induction(m) + raw_pol(t)
      end do
      ! Slot one is the fragment sum, taken from `monomer_sum` rather than
      ! re-added here: `dE_I^0` *is* `E_I^0`, and one sum of the same numbers is
      ! enough.
      res%level_vacuum(1) = res%monomer_sum
      do m = 2, level
         res%nmer_correction = res%nmer_correction + res%level_vacuum(m)
         res%induction_correction = res%induction_correction + res%level_induction(m)
      end do
      do t = 1, n_terms
         if (term_size(t) < 2) cycle
         res%nmer_correlation = res%nmer_correlation + raw_corr(t)
      end do
   end subroutine quantum_subsets

   pure function group_atoms(owner, members) result(idx)
      !! Every atom of a group of fragments, fragment by fragment in group order
      integer, intent(in) :: owner(:)
      integer, intent(in) :: members(:)
      integer, allocatable :: idx(:)

      integer :: k

      ! Fragment by fragment in group order, which is the order the group's
      ! basis functions come out in and the order a dimer was assembled in
      ! before there were larger groups.
      idx = gather(owner, members(1))
      do k = 2, size(members)
         idx = [idx, gather(owner, members(k))]
      end do
   end function group_atoms

   function members_text(members) result(text)
      !! A group's fragments as `1-2-5`, for a log line
      integer, intent(in) :: members(:)
      character(len=:), allocatable :: text

      integer :: k

      text = to_char(members(1))
      do k = 2, size(members)
         text = text//"-"//to_char(members(k))
      end do
   end function members_text

   subroutine nmer_energy(z, symbols, xyz, charge, opts, energy, correlation, error)
      !! One group's restricted Hartree-Fock energy, in vacuo
      !!
      !! Two fragments or twenty: the group arrives as one atom list and this is
      !! an ordinary closed-shell SCF on it. Cartesian, to match the monomer
      !! SCFs `make_efp_potential` ran; see the module header.
      !!
      !! **No embedding of any kind**, and that is a property of the method and
      !! not a simplification: the fragments outside the group are not in this
      !! Hamiltonian, their interaction with it being carried by the
      !! effective-fragment pair terms and the one total induction. It is also
      !! what makes the many-body differences telescope, so an embedded variant
      !! would be a different method rather than a better version of this one.
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: symbols(:)
      real(dp), intent(in) :: xyz(:, :)
      integer, intent(in) :: charge
      type(efmo_options_t), intent(in) :: opts
      real(dp), intent(out) :: energy
      real(dp), intent(out) :: correlation
         !! The correlation part of `energy`, which already holds it. Zero on a
         !! Hartree-Fock run.
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: guess_density(:, :)
      integer :: guess_kind, nelec

      energy = 0.0_dp
      correlation = 0.0_dp
      nelec = sum(z) - charge
      if (nelec < 2 .or. mod(nelec, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "efmo: a group with "//to_char(nelec)// &
                        " electrons is not closed-shell. EFMO is restricted "// &
                        "Hartree-Fock for now, so the fragment charges have to "// &
                        "leave every fragment and every group with an even count.")
         return
      end if

      call build_czt_molecule(z, symbols, xyz, trim(opts%basis), mol, error, &
                              force_cartesian=.true.)
      if (error%has_error()) return
      call build_restricted_guess(mol, trim(opts%guess), guess_kind, guess_density, error)
      if (error%has_error()) then
         call mol%destroy()
         return
      end if

      call run_czt_rhf(mol, nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                       opts%scf_density_tol, opts%verbose, scf, error, &
                       guess=guess_kind, guess_density=guess_density, &
                       grad_tol=opts%scf_grad_tol, scf=opts%scf)
      call mol%destroy()
      if (error%has_error()) return
      if (.not. scf%converged .and. .not. opts%scf%allow_crap_scf) then
         call error%set(ERROR_VALIDATION, "efmo: a group SCF did not converge, so the "// &
                        "many-body correction it feeds is not trustworthy. Set "// &
                        "keywords.scf.allow_crap_scf to finish anyway.")
         return
      end if
      energy = scf%energy

      ! The same correlation step the monomers got, on the group's own
      ! orbitals. `dE_S^0` is then a difference of energies of one model, which
      ! is the only reading under which it is an interaction energy. The frozen
      ! core is counted from this group's elements, so it is the sum of its
      ! fragments' cores and both sides of the difference freeze the same set.
      call fragment_correlation(z, symbols, xyz, nelec, scf, opts, correlation, error)
      if (error%has_error()) return
      energy = energy + correlation
   end subroutine nmer_energy

   subroutine total_polarization(frags, shifts, damping, energy, error)
      !! `E_pol^total`: induction over every fragment at once
      !!
      !! The same call `efp_interaction_energy` makes internally, on the same
      !! system, so the pair terms subtracted from the near dimers cancel
      !! against exactly what is here.
      type(efp_fragment_t), intent(in) :: frags(:)
      real(dp), intent(in) :: shifts(:, :)
      real(dp), intent(in) :: damping
         !! The same number every `E_IJ^pol` was solved with; see
         !! `efmo_options_t%induction_damping`.
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error

      type(efp_system_t) :: system

      energy = 0.0_dp
      if (size(frags) < 2) return
      call build_efp_system(frags, shifts, system, error)
      if (error%has_error()) return
      energy = polarization_energy(system, frags, error, damping=damping)
      call system%destroy()
   end subroutine total_polarization

   subroutine report(res, opts)
      !! The breakdown, at info level
      type(efmo_result_t), intent(in) :: res
      type(efmo_options_t), intent(in) :: opts

      character(len=160) :: line
      integer :: k

      call logger%info("============================================================")
      call logger%info("  EFMO, level "//to_char(min(opts%level, &
                                                     size(res%monomer_energy)))//", R_cut = "//to_char(opts%rcut)// &
                       " (unitless)")
      if (opts%induction_damping > 0.0_dp) then
         call logger%info("  induction field damped, a = "// &
                          to_char(opts%induction_damping))
      end if
      call logger%info("------------------------------------------------------------")
      call logger%info("  fragment           E_I^0 / Hartree")
      do k = 1, size(res%monomer_energy)
         write (line, "(A,I6,F26.10)") "  ", k, res%monomer_energy(k)
         call logger%info(trim(line))
      end do
      call logger%info("------------------------------------------------------------")
      call logger%info("  pair      R_IJ  class      contribution / Hartree")
      do k = 1, size(res%pairs)
         if (res%pairs(k)%qm) then
            write (line, "(A,I4,A,I4,F8.3,A,F22.10)") "  ", res%pairs(k)%i, "-", &
               res%pairs(k)%j, res%pairs(k)%r, "  QM  ", &
               res%pairs(k)%e_dimer - res%monomer_energy(res%pairs(k)%i) &
               - res%monomer_energy(res%pairs(k)%j) - res%pairs(k)%e_pair_pol
         else
            write (line, "(A,I4,A,I4,F8.3,A,F22.10)") "  ", res%pairs(k)%i, "-", &
               res%pairs(k)%j, res%pairs(k)%r, "  EFP ", &
               res%pairs(k)%electrostatics + res%pairs(k)%dispersion &
               + res%pairs(k)%exchange_repulsion + res%pairs(k)%charge_transfer
         end if
         call logger%info(trim(line))
      end do
      if (size(res%level_vacuum) >= 2) then
         ! Per level rather than one lumped correction: what a level-three run
         ! is *for* is the size of its three-body term, and a single number
         ! cannot say whether the expansion is converging.
         call logger%info("------------------------------------------------------------")
         call logger%info("  level  groups        sum dE_S^0          sum dE_S^pol")
         do k = 1, size(res%level_vacuum)
            write (line, "(A,I5,I8,F20.10,F22.10)") "  ", k, res%level_count(k), &
               res%level_vacuum(k), res%level_induction(k)
            call logger%info(trim(line))
         end do
      end if
      call logger%info("------------------------------------------------------------")
      call logger%info("  monomers            sum E_I^0        "//to_char(res%monomer_sum))
      call logger%info("  QM groups           sum dE_S^0       "//to_char(res%nmer_correction))
      call logger%info("  QM groups           - sum dE_S^pol   "//to_char(-res%induction_correction))
      call logger%info("  EFP dimers          Coulomb          "//to_char(res%far_electrostatics))
      call logger%info("  EFP dimers          dispersion       "//to_char(res%far_dispersion))
      call logger%info("  EFP dimers          exchange rep.    "//to_char(res%far_exchange_repulsion))
      call logger%info("  EFP dimers          charge transfer  "//to_char(res%far_charge_transfer))
      call logger%info("  all fragments       E_pol^total      "//to_char(res%polarization_total))
      if (opts%correlation /= EFMO_CORR_NONE) then
         ! Inside the two sums above, not beside them: a correlated `E_I^0` is
         ! the monomer energy of eq 6 and not a term added to it.
         call logger%info("    of which correlation, monomers  "// &
                          to_char(res%monomer_correlation))
         call logger%info("    of which correlation, QM groups "// &
                          to_char(res%nmer_correlation))
      end if
      call logger%info("------------------------------------------------------------")
      call logger%info("  QM groups "//to_char(res%n_qm_groups)//" (of which pairs "// &
                       to_char(res%n_qm_pairs)//"), EFP dimers "// &
                       to_char(res%n_efp_pairs))
      call logger%info("  EFMO total energy   "//to_char(res%energy)//" Hartree")
      call logger%info("============================================================")
   end subroutine report

end module mqc_czt_efmo
