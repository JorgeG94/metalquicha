!! Bonding analysis on a converged wave function
module mqc_libcint_bonding
   !! One entry point that takes a converged set of molecular orbitals and
   !! produces the quasi-atomic bonding picture: atomic populations and charges,
   !! then every orbital pair worth reporting with its bond order, its kinetic
   !! bond order and what kind of orbital each end is.
   !!
   !! The pipeline itself lives in `mqc_aambs`, `mqc_libcint_quao` and
   !! `mqc_libcint_quao_report`. What is here is the assembly: which spaces to
   !! build, in what order, and what to hand to the report.
   !!
   !! **This is a Hartree-Fock analysis unless a correlated density is passed.**
   !! `quasi_atomic_orbitals` accepts one, through `occupations`; without it
   !! what is analysed is the reference determinant.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_io, only: to_char
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_aambs, only: aambs_dimensions, aambs_dimensions_t, aambs_element_counts
   use mqc_libcint_integrals, only: libcint_molecule_t, mixed_basis_overlap
   use mqc_libcint_quao, only: build_aambs_molecule, mo_aambs_overlap, &
                               valence_virtual_orbitals, vvo_result_t, &
                               aambs_atom_ranges, quasi_atomic_orbitals, &
                               quao_result_t, orient_quasi_atomic_orbitals, &
                               kinetic_bond_orders
   use mqc_physical_constants, only: HARTREE_TO_KCALMOL
   use mqc_timing, only: timing_report_t
   use mqc_libcint_atomic_guess, only: free_atom_energies
   use mqc_libcint_casci, only: active_space_integrals, run_libcint_casci, casci_result_t
   use mqc_rdm, only: active_space_rdms, rdm_energy
   use mqc_determinants, only: link_table_t, build_link_table
   use mqc_ci_transform, only: string_rotation_matrix, transform_ci_vector, &
                               rotation_is_orthogonal, orthogonality_defect, &
                               scatter_restricted_ci
   use mqc_ormas_space, only: ormas_space_t
   use mqc_libcint_ieda, only: kinetic_decomposition, print_kinetic_decomposition, &
                               nuclear_attraction_per_atom, quao_nuclear_attraction, &
                               nuclear_decomposition, print_nuclear_decomposition, &
                               combine_quao_sets, quao_interference, kinetic_total, &
                               quao_eris, two_electron_decomposition, &
                               print_two_electron_decomposition, &
                               nuclear_repulsion_pairs, print_total_decomposition, &
                               print_interatomic_split, screened_nucleus_split, &
                               project_no_sharing, &
                               active_cumulant, quao_projection, transform_cumulant
   use mqc_libcint_quao_report, only: quao_labels_t, label_quasi_atomic_orbitals, &
                                      print_quao_report
   implicit none
   private

   public :: run_quao_analysis
   public :: BONDING_NONE, BONDING_GMS_QUAO
   public :: NO_SHARING_CI_TRANSFORM, NO_SHARING_CI_RESOLVE
   public :: bonding_analysis_kind
   public :: bonding_analysis_name
   public :: no_sharing_ci_kind
   public :: valence_wavefunction_t

   integer, parameter :: BONDING_NONE = 0
   integer, parameter :: BONDING_GMS_QUAO = 1
      !! The Ruedenberg quasi-atomic analysis as GAMESS implements it. Named for
      !! the reference implementation rather than the papers, because the labels
      !! and the thresholds are GAMESS's and only the underlying quantities are
      !! the papers'.

   real(dp), parameter :: SUM_RULE_TOL = 1.0e-6_dp
      !! The population analysis must account for every electron. Loose next to
      !! the others because it is a sum over the whole molecule, so it carries
      !! the rounding of every term in it.

   real(dp), parameter :: SPAN_TOL = 1.0e-6_dp
      !! How far the valence virtual space may fall short of spanning what the
      !! minimal basis asks for before the analysis is refused.

   real(dp), parameter :: BALANCE_TOL = 1.0e-8_dp
      !! A decomposition has to add back up to the thing it decomposed. Tighter
      !! than the sum rule above: these compare a handful of terms rather than a
      !! sum over the molecule.

   real(dp), parameter :: OCCUPATION_FLOOR = 1.0e-14_dp
      !! Occupations below this contribute nothing and are skipped, which also
      !! keeps them out of divisions.

   integer, parameter :: NO_SHARING_CI_TRANSFORM = 0
      !! Carry a CI vector into the quasi-atomic basis rather than solving
      !! there: the one the calculation already converged if it is over the
      !! full valence space, and otherwise one solved here in the molecular
      !! orbital basis. The default.
   integer, parameter :: NO_SHARING_CI_RESOLVE = 1
      !! Solve it again, in the quasi-atomic basis.

   type :: valence_wavefunction_t
      !! A multiconfigurational wave function the calculation already converged
      !!
      !! Offered to the no-sharing analysis so that it need not converge its own
      !! when the two would be the same wave function. Whether they are is
      !! decided by `inherit_valence_ci`, not by the caller: the analysis is
      !! defined on the *full valence* space.
      !!
      !! `ci` and `orbitals` belong together -- the coefficients are over
      !! determinants of those orbitals, in that order.
      real(dp), allocatable :: ci(:, :)
         !! (n_alpha_strings, n_beta_strings). Left unallocated by a restricted
         !! space, which carries `ci_flat` and `ormas` instead.
      real(dp), allocatable :: ci_flat(:)
         !! The coefficients of an occupation-restricted space, whose
         !! determinants are not a rectangle. Written out over the complete
         !! space before use -- see `scatter_restricted_ci`: the rotation
         !! carries such a wave function *out* of its own space, so the complete
         !! one is the only place the answer fits.
      type(ormas_space_t) :: ormas
         !! The partition `ci_flat` is addressed by
      real(dp), allocatable :: orbitals(:, :)   !! (n_ao, n_active)
      integer :: n_inactive = 0
      integer :: n_active = 0
      integer :: n_alpha = 0
      integer :: n_beta = 0
      real(dp) :: energy = 0.0_dp
         !! What the calculation reported for it. The transformed vector has to
         !! reproduce this against the quasi-atomic Hamiltonian, which is what
         !! makes accepting the offer safe.
   end type valence_wavefunction_t

contains

   pure function bonding_analysis_kind(name) result(kind)
      !! Parse a deck's `properties.bonding_analysis` value
      character(len=*), intent(in) :: name
      integer :: kind

      select case (trim(adjustl(name)))
      case ("", "none")
         kind = BONDING_NONE
      case ("gms_quao", "quao")
         kind = BONDING_GMS_QUAO
      case default
         kind = -1
      end select
   end function bonding_analysis_kind

   pure function bonding_analysis_name(kind) result(name)
      !! The spelling a deck would use, for messages
      integer, intent(in) :: kind
      character(len=:), allocatable :: name

      select case (kind)
      case (BONDING_GMS_QUAO)
         name = "gms_quao"
      case default
         name = "none"
      end select
   end function bonding_analysis_name

   pure function no_sharing_ci_kind(name) result(kind)
      !! Parse a deck's `properties.bonding_analysis.no_sharing_ci` value
      !!
      !! Returns `-1` for anything else, which the caller refuses. The schema
      !! catches a misspelling before a deck reaches here, but a caller through
      !! the C or Python interface never passes through the schema.
      character(len=*), intent(in) :: name
      integer :: kind

      select case (trim(adjustl(name)))
      case ("", "transform")
         kind = NO_SHARING_CI_TRANSFORM
      case ("resolve")
         kind = NO_SHARING_CI_RESOLVE
      case default
         kind = -1
      end select
   end function no_sharing_ci_kind

   subroutine run_quao_analysis(mol, atomic_numbers, element_symbols, coordinates, &
                                orbitals, n_electrons, error, verbose, threshold, &
                                occupations, active_orbitals, active_dm1, active_dm2, &
                                reference_energy, energy_decomposition, no_sharing, &
                                no_sharing_ci, valence_wavefunction, &
                                restrict_localization, atom_energy, &
                                free_atom_energy, pair_energy, pair_classical, &
                                formation_energy)
      !! The quasi-atomic bonding analysis, start to finish
      ! TODO(mqc): the dummy arguments are interrupted by local declarations --
      ! `aambs` through `labels` sit between `occupations` and `active_orbitals`
      ! -- so the interface cannot be read off the top of the declaration block.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)     !! (3, natm), Bohr
      real(dp), intent(in) :: orbitals(:, :)        !! Converged, (n_ao, n_mo)
      integer, intent(in) :: n_electrons
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose
      real(dp), intent(in), optional :: threshold
         !! kcal/mol. Pairs whose kinetic bond order is weaker are counted and
         !! not printed.
      real(dp), intent(in), optional :: occupations(:)
         !! (n_mo) occupation of each orbital, which must then be a *natural*
         !! orbital -- the density has to be diagonal in the basis given for
         !! these numbers to be its whole content.
         !!
         !! Absent means the reference determinant: two electrons in each of the
         !! first `n_occupied` orbitals. Present is what a correlated wave
         !! function supplies, and is the only difference between analysing one
         !! and analysing the other.

      type(libcint_molecule_t) :: aambs
      type(aambs_dimensions_t) :: dims
      type(vvo_result_t) :: vvo
      type(quao_result_t) :: quao
      type(quao_labels_t) :: labels
      real(dp), intent(in), optional :: active_orbitals(:, :)
         !! (n_ao, n_active), the orbitals `active_dm2` is expressed in. These
         !! are the optimised active orbitals, not the natural ones the rest of
         !! the analysis runs on: the two-particle density belongs to the basis
         !! it was computed in, and is projected from there.
      real(dp), intent(in), optional :: active_dm1(:, :)
      real(dp), intent(in), optional :: active_dm2(:, :, :, :)
         !! Present together with the two above, or not at all. Their presence
         !! is what makes the energy decomposition correlated rather than a
         !! single-determinant one.
      real(dp), intent(in), optional :: reference_energy
         !! The energy the calculation reported, to check the decomposition
         !! against. Required for a correlated wave function and optional
         !! otherwise: the reference built here from the orbitals is a
         !! determinant expression, so it is the energy of the wave function
         !! only when the wave function is a determinant. Passing the CI energy
         !! is what keeps the correlated case checked rather than trusted.
      logical, intent(in), optional :: energy_decomposition
         !! Resolve the energy onto atoms and atom pairs. Off by default on
         !! cost: the two-electron term needs the dense `n_ao^4` integral array,
         !! where the bonding tables need only one-electron integrals.
      logical, intent(in), optional :: no_sharing
         !! Run the no-sharing analysis, which needs a full valence CI expanded
         !! over the quasi-atomic orbitals and projects it. Off by default
         !! because that CI is factorial in the valence shell -- ethane is
         !! eleven million determinants and benzene is out of reach.
      character(len=*), intent(in), optional :: no_sharing_ci
         !! How that expansion is obtained: `"transform"`, the default, or
         !! `"resolve"`. See `no_sharing_analysis`.
      logical, intent(in), optional :: restrict_localization
         !! Confine the localization to the wave function's
         !! occupation-restricted subspaces, so no rotation mixes two of them.
         !! Off by default: the constraint costs atomic character, and buys only
         !! the ability to keep a restricted wave function in its own space,
         !! which matters when writing it out over the complete one will not fit.
      type(valence_wavefunction_t), intent(in), optional :: valence_wavefunction
         !! A converged multiconfigurational wave function to use instead of
         !! solving one, if it happens to be over the full valence space. Its
         !! suitability is checked here rather than asserted by the caller, and
         !! an unsuitable one costs nothing but a line saying so.
      real(dp), intent(out), optional, allocatable :: atom_energy(:)
      real(dp), intent(out), optional, allocatable :: free_atom_energy(:)
      real(dp), intent(out), optional, allocatable :: pair_energy(:, :)
      real(dp), intent(out), optional, allocatable :: pair_classical(:, :)
      real(dp), intent(out), optional :: formation_energy
         !! The decomposition, for a caller that wants the numbers rather than
         !! the tables. Left unallocated when the energy decomposition did not
         !! run, which is what distinguishes "not asked for" from "came out
         !! zero".
         !!
         !! `pair_energy` and `pair_classical` carry the full pair energy in
         !! both (A,B) and (B,A), as everything in the decomposition does, so
         !! the total is `sum(atom_energy) + 0.5*sum(pair_energy)`.

      real(dp), allocatable :: s_mbs(:, :), mixed(:, :), projection(:, :)
      real(dp), allocatable :: valence_internal(:, :), kinetic(:, :)
      real(dp), allocatable :: kbo(:, :), interference(:, :), populations(:)
      real(dp), allocatable :: valence_density(:, :)
      real(dp), allocatable :: kin_intra(:), kin_inter(:, :)
      real(dp), allocatable :: v_atom_ao(:, :, :), v_atom_quao(:, :, :)
      real(dp), allocatable :: nuc_intra(:), nuc_inter(:, :)
      real(dp), allocatable :: full_interference(:, :), h_core(:, :), density_ao(:, :)
      type(quao_result_t) :: core_quao, full_quao
      real(dp) :: one_electron, reference, two_electron, total_energy
      real(dp), allocatable :: eri_ao(:, :, :, :), eri_quao(:, :, :, :)
      real(dp), allocatable :: two_intra(:), two_inter(:, :)
      real(dp), allocatable :: rep_inter(:, :), tot_intra(:), tot_inter(:, :)
      real(dp), allocatable :: s_ao(:, :), u_active(:, :)
      real(dp), allocatable :: cumulant(:, :, :, :), cumulant_quao(:, :, :, :)
      logical :: correlated, want_energy, adopted, starved, restrict
      type(timing_report_t) :: clk
      integer :: ci_route
      integer, allocatable :: restriction(:)
      real(dp) :: span_deficit, formation
      real(dp), allocatable :: free_energy(:), adaptation(:)
      real(dp), allocatable :: nuc_coulomb(:, :), two_coulomb(:, :), classical(:, :)
      real(dp), allocatable :: nuc_core(:), nuc_val(:)
      integer, allocatable :: core_off(:), core_n(:), val_off(:), val_n(:)
      integer, allocatable :: order(:)
      character(len=160) :: line
      character(len=8) :: label
      integer :: natm, iatom, i, core, valence
      logical :: loud

      if (error%has_error()) return
      ! An `intent(out)` scalar is undefined until it is written, and this one is
      ! written only where the decomposition runs. Defined here so a caller that
      ! reads it unconditionally reads a zero rather than a stack value.
      if (present(formation_energy)) formation_energy = 0.0_dp
      loud = .true.
      if (present(verbose)) loud = verbose
      natm = size(atomic_numbers)

      ! The dimensions first: they are pure counting and they refuse the cases
      ! this cannot handle -- an element past xenon, an odd electron count --
      ! before any integral is built.
      call clk%start()
      call aambs_dimensions(atomic_numbers, n_electrons, dims, error)
      if (error%has_error()) return

      call build_aambs_molecule(atomic_numbers, element_symbols, coordinates, aambs, error)
      if (error%has_error()) return
      call aambs%overlap(s_mbs)
      call mixed_basis_overlap(mol, aambs, mixed, error)
      if (error%has_error()) return

      call mo_aambs_overlap(orbitals, mixed, s_mbs, projection, error)
      call valence_virtual_orbitals(orbitals, projection, dims, vvo, error)
      call aambs_atom_ranges(atomic_numbers, aambs, core_off, core_n, val_off, &
                             val_n, error)
      if (error%has_error()) return

      ! The valence-internal space: occupied orbitals that are not core,
      ! followed by the valence-virtual ones.
      allocate (valence_internal(size(orbitals, 1), dims%n_valocc + vvo%n_vvo))
      valence_internal(:, 1:dims%n_valocc) = &
         orbitals(:, dims%n_core + 1:dims%n_occupied)
      valence_internal(:, dims%n_valocc + 1:) = vvo%orbitals

      ! **Unless the calculation already has a valence space, in which case use
      ! it.** The extraction above is the recipe for a wave function that does
      ! not span the valence shell, recovering the missing part by an SVD
      ! projection of the minimal-basis orbitals on the virtual ones.
      !
      ! A full valence MCSCF has no missing part: its active space *is* the
      ! valence space, and re-deriving one here would produce a different
      ! subspace of the same dimension, so the analysis would decompose a wave
      ! function the calculation never computed. GAMESS takes the same branch at
      ! `vvos.src:540`.
      call adopt_valence_space(valence_wavefunction, dims, vvo%n_vvo, n_electrons, &
                               valence_internal, loud, adopted)

      ! An occupation-restricted partition constrains the localization, and can
      ! only do so against the orbitals it was defined over, so it travels only
      ! when those orbitals were adopted above. No rotation may then cross a
      ! subspace, which is what keeps a restricted wave function invariant under
      ! the transformation.
      restrict = .false.
      if (present(restrict_localization)) restrict = restrict_localization
      if (restrict .and. adopted .and. present(valence_wavefunction)) then
         if (valence_wavefunction%ormas%n_subspaces > 1) then
            restriction = valence_wavefunction%ormas%first_orbital
            if (loud) then
               write (line, "(a,i0,a)") "    localization                confined "// &
                  "to ", size(restriction), " occupation-restricted subspaces"
               call logger%info(trim(line))
            end if
         end if
      end if

      ! The density in the valence-internal basis. For a reference determinant it
      ! is two on the occupied diagonal and zero elsewhere, which is what
      ! `quasi_atomic_orbitals` assumes when nothing is passed. For a correlated
      ! one it is neither diagonal nor idempotent, and has to be projected:
      !
      !     D_vi = W^T diag(n) W ,   W = C^T S V
      !
      ! with `V` the valence-internal orbitals and `C` the natural orbitals.
      ! **The projection is not optional even though the density is diagonal in
      ! `C`.** The valence-virtual orbitals are combinations of virtuals with
      ! *different* occupations, so the density restricted to this space has
      ! off-diagonal elements even when the full one does not. Reading
      ! occupations straight off the diagonal would silently drop them.
      if (present(occupations)) then
         call project_density(mol, orbitals, occupations, valence_internal, &
                              valence_density, error)
         if (error%has_error()) return
         if (allocated(restriction)) then
            call quasi_atomic_orbitals(atomic_numbers, valence_internal, dims%n_valocc, &
                                       mixed, val_off, val_n, quao, error, &
                                       valence_density=valence_density, &
                                       subspaces=restriction, starved=starved)
            if (error%has_error()) return
            if (starved) then
               ! The partition cannot give every atom an orbital. Not a fault:
               ! these subspaces are usually grouped by orbital energy, so one
               ! strong atom can win every slot. GAMESS stops here; this does
               ! not have to, because the unconstrained construction plus
               ! writing the wave function out over the complete space reaches
               ! the same answer.
               call logger%warning("  the occupation-restricted subspaces cannot "// &
                                   "give every atom a quasi-atomic orbital, so the "// &
                                   "localization is left unconstrained and the wave "// &
                                   "function is written out over the complete space "// &
                                   "instead")
               deallocate (restriction)
               call quasi_atomic_orbitals(atomic_numbers, valence_internal, &
                                          dims%n_valocc, mixed, val_off, val_n, quao, &
                                          error, valence_density=valence_density)
            end if
         else
            call quasi_atomic_orbitals(atomic_numbers, valence_internal, dims%n_valocc, &
                                       mixed, val_off, val_n, quao, error, &
                                       valence_density=valence_density)
         end if
      else
         call quasi_atomic_orbitals(atomic_numbers, valence_internal, dims%n_valocc, &
                                    mixed, val_off, val_n, quao, error)
      end if
      call clk%lap("quasi-atomic orbitals")
      call orient_quasi_atomic_orbitals(quao, error)
      if (error%has_error()) return
      call clk%lap("orientation")

      call mol%kinetic(kinetic)
      call kinetic_bond_orders(quao, kinetic, kbo, error, interference)
      call label_quasi_atomic_orbitals(quao, mol, aambs, mixed, s_mbs, &
                                       atomic_numbers, coordinates, labels, error)
      if (error%has_error()) return

      if (loud) then
         ! Populations first. The core electrons are added back because the
         ! construction only ever sees the valence.
         allocate (populations(natm))
         populations = 0.0_dp
         do i = 1, quao%n_quao
            populations(quao%atom_of(i)) = populations(quao%atom_of(i)) &
                                           + quao%population_bond_order(i, i)
         end do
         do iatom = 1, natm
            call aambs_element_counts(atomic_numbers(iatom), core, valence, error)
            if (error%has_error()) return
            populations(iatom) = populations(iatom) + 2.0_dp*real(core, dp)
         end do

         call logger%info("")
         call logger%info("  ==== quasi-atomic bonding analysis "// &
                          "================================")
         call logger%info("")
         call logger%info("  atoms")
         call logger%info("     atom      population      charge")
         do iatom = 1, natm
            ! Into a fixed-width buffer rather than straight through an `i0`:
            ! the atom number is variable width, so writing it inline shifts
            ! every column right the moment a molecule has ten atoms.
            write (label, "(a,i0)") trim(adjustl(element_symbols(iatom)))//" ", iatom
            write (line, "(4x,a8,f13.6,f12.6)") label, populations(iatom), &
               real(atomic_numbers(iatom), dp) - populations(iatom)
            call logger%info(trim(line))
         end do
         write (line, "(4x,a8,f13.6)") "total   ", sum(populations)
         call logger%info(trim(line))

         ! **The populations need not sum to the electron count, and when they
         ! do not it is a statement about the analysis rather than an error.**
         ! The quasi-atomic orbitals span the occupied valence orbitals plus the
         ! valence-virtual ones, and nothing else. A reference determinant puts
         ! every valence electron inside that space, so the sum is exact. A
         ! correlated density does not: the valence-virtual orbitals are chosen
         ! to look like free-atom orbitals, not to be natural orbitals, so
         ! whatever occupation sits outside their span is not recovered.
         !
         ! It is small, and it is the measure of how much of the correlated wave
         ! function this analysis is describing. Printed rather than hidden.
         if (abs(sum(populations) - real(n_electrons, dp)) > SUM_RULE_TOL) then
            write (line, "(a,f10.6,a,i0,a)") &
               "    outside the valence space ", &
               real(n_electrons, dp) - sum(populations), " of ", n_electrons, &
               " electrons"
            call logger%info(trim(line))
         end if

         call print_quao_report(.true., quao, labels, interference, element_symbols, &
                                dims%n_core, kinetic_bond_order=kbo, &
                                threshold=threshold)

         ! The same interference matrix regrouped by atom rather than by orbital
         ! pair. Unscaled, so unlike the kinetic bond orders above these are
         ! energies and they add up to one -- the valence kinetic energy, since
         ! the core is not in this basis.
      end if

      ! Outside the `loud` block that used to hold it. It is read below as
      ! `loud .and. want_energy`, and Fortran does not promise to stop at a
      ! false first operand -- so a caller passing `verbose = .false.` read it
      ! undefined. No caller passes `verbose` today, which is the only reason
      ! this has never been seen.
      want_energy = .false.
      if (present(energy_decomposition)) want_energy = energy_decomposition
      if (present(no_sharing)) then
         ! Asking for the no-sharing analysis is asking for the decomposition
         ! it is part of; refusing on a technicality would be unhelpful.
         if (no_sharing) want_energy = .true.
      end if
      ! TODO(mqc): `want_energy` is assigned only inside the `if (loud)` block
      ! above and read in `if (loud .and. want_energy)` below, so a caller
      ! passing `verbose = .false.` reads it undefined and silently gets no
      ! energy decomposition however `energy_decomposition` was set.

      ci_route = NO_SHARING_CI_TRANSFORM
      if (present(no_sharing_ci)) then
         if (len_trim(no_sharing_ci) > 0) ci_route = no_sharing_ci_kind(no_sharing_ci)
      end if
      if (ci_route < 0) then
         call error%set(ERROR_VALIDATION, "the no-sharing CI route is '"// &
                        trim(adjustl(no_sharing_ci))//"'. It is 'transform' or "// &
                        "'resolve'.")
         return
      end if

      if (loud .and. want_energy) then
         ! The decomposition needs the core, which the bonding analysis above
         ! does not: it carries most of the kinetic energy and most of the
         ! nuclear attraction, and without it the totals below would sum to a
         ! number nothing else in the calculation knows. Built by the same
         ! construction, handed the core orbitals and the core minimal-basis
         ! ranges. Left unoriented deliberately -- see `quao_interference`.
         if (dims%n_core > 0) then
            call quasi_atomic_orbitals(atomic_numbers, orbitals(:, 1:dims%n_core), &
                                       dims%n_core, mixed, core_off, core_n, &
                                       core_quao, error)
            if (error%has_error()) return
            call combine_quao_sets(core_quao, quao, full_quao, error)
         else
            full_quao = quao
         end if
         if (error%has_error()) return

         call quao_interference(full_quao, kinetic, full_interference, error)
         if (error%has_error()) return
         call kinetic_decomposition(full_quao, full_interference, natm, kin_intra, &
                                    kin_inter, error)
         if (error%has_error()) return
         call print_kinetic_decomposition(kin_intra, kin_inter, element_symbols)

         ! The other half of the one-electron energy. Split by which nucleus is
         ! attracting, so a term carries three atomic labels rather than two and
         ! an atom's density sitting in a foreign nuclear field becomes an
         ! interaction between that pair rather than a property of one atom.
         call nuclear_attraction_per_atom(mol, atomic_numbers, coordinates, v_atom_ao, error)
         if (error%has_error()) return
         call quao_nuclear_attraction(full_quao, v_atom_ao, v_atom_quao, error)
         if (error%has_error()) return
         call nuclear_decomposition(full_quao, full_quao%population_bond_order, &
                                    v_atom_quao, natm, nuc_intra, nuc_inter, error, &
                                    coulomb=nuc_coulomb)
         if (error%has_error()) return
         call print_nuclear_decomposition(nuc_intra, nuc_inter, element_symbols)

         ! A valence electron does not see the bare nucleus but one screened by
         ! the core, so its attraction to its own nucleus is divided in the
         ! ratio of the two charges. Pure relabelling: the shares must add back
         ! to the term they came from, which is asserted here because a split
         ! that quietly lost part of it would still print two plausible columns.
         call screened_nucleus_split(full_quao, full_quao%population_bond_order, &
                                     v_atom_quao, dims%n_core, atomic_numbers, natm, &
                                     nuc_core, nuc_val, error)
         if (error%has_error()) return
         if (maxval(abs(nuc_core + nuc_val - nuc_intra)) > BALANCE_TOL) then
            call error%set(ERROR_VALIDATION, "the screened-nucleus split does not sum "// &
                           "back to the own-nucleus attraction it divides, so the two "// &
                           "charge fractions do not add to one.")
            return
         end if
         call logger%info("     own nucleus, by charge      to core     to valence")
         do iatom = 1, natm
            write (label, "(a,i0)") trim(adjustl(element_symbols(iatom)))//" ", iatom
            write (line, "(4x,a8,18x,2f14.3)") label, &
               nuc_core(iatom)*1000.0_dp, nuc_val(iatom)*1000.0_dp
            call logger%info(trim(line))
         end do

         ! The two-electron energy, where a term carries four atomic labels
         ! rather than three and both a bra and a ket distribution can be spread
         ! across two atoms at once.
         call mol%eris(eri_ao)
         call quao_eris(full_quao, eri_ao, eri_quao, error)
         if (error%has_error()) return

         ! What correlation adds. The determinant expression below is evaluated
         ! over the whole quasi-atomic basis whatever the wave function; the
         ! cumulant is the part no determinant can produce, it is zero unless
         ! all four of its indices are active, and it is carried here from the
         ! orbitals it was computed in rather than from the natural ones.
         correlated = present(active_dm2) .and. present(active_dm1) .and. &
                      present(active_orbitals)
         if (correlated) then
            call mol%overlap(s_ao)
            call quao_projection(full_quao, s_ao, active_orbitals, u_active, &
                                 span_deficit, error)
            if (error%has_error()) return
            ! The active orbitals are not quite inside the span -- they were
            ! optimised to lower an energy, the valence-virtual ones to look
            ! atomic. Reported for the same reason the population shortfall is.
            if (span_deficit > SPAN_TOL) then
               write (line, "(a,es10.2)") &
                  "    active orbitals outside the quasi-atomic span ", span_deficit
               call logger%info(trim(line))
            end if
            call active_cumulant(active_dm1, active_dm2, cumulant)
            call transform_cumulant(u_active, cumulant, cumulant_quao, error)
            if (error%has_error()) return
            call two_electron_decomposition(full_quao, full_quao%population_bond_order, &
                                            eri_quao, natm, two_intra, two_inter, &
                                            two_electron, error, cumulant=cumulant_quao, &
                                            coulomb=two_coulomb)
         else
            call two_electron_decomposition(full_quao, full_quao%population_bond_order, &
                                            eri_quao, natm, two_intra, two_inter, &
                                            two_electron, error, coulomb=two_coulomb)
         end if
         if (error%has_error()) return
         call print_two_electron_decomposition(two_intra, two_inter, element_symbols)

         ! **The check the whole thing exists to pass**, further down: the energy
         ! the decomposition arrives at against the energy of the same wave
         ! function computed in the atomic-orbital basis, with nothing
         ! quasi-atomic anywhere in it. For a reference determinant the two are
         ! equal to rounding; for a correlated density the difference is the
         ! occupation living outside the valence-virtual span, which is a
         ! statement about the analysis rather than an error and so is printed.
         !
         ! The nuclear repulsion needs no decomposing -- it is a sum over pairs
         ! already -- so adding it in the same convention resolves the whole
         ! energy into atoms and pairs with nothing left outside.
         call nuclear_repulsion_pairs(atomic_numbers, coordinates, rep_inter, error)
         if (error%has_error()) return
         allocate (tot_intra(natm), tot_inter(natm, natm))
         tot_intra = kin_intra + nuc_intra + two_intra
         tot_inter = kin_inter + nuc_inter + two_inter + rep_inter
         call print_total_decomposition(tot_intra, tot_inter, element_symbols)

         ! Everything an electrostatic model could account for: each atom's
         ! density in the other's nuclear field, the repulsion of their two
         ! charge clouds, and the repulsion of their nuclei. The rest is
         ! interference, and the kinetic contribution is entirely interference
         ! by construction -- two orbitals on different atoms is what the word
         ! means.
         allocate (classical(natm, natm))
         classical = nuc_coulomb + two_coulomb + rep_inter
         call print_interatomic_split(tot_inter, classical, element_symbols)

         ! The atoms as they would be on their own. Subtracting them turns a
         ! table of large numbers into the energy of formation: what it costs to
         ! prepare each atom in the shape the molecule needs, against what the
         ! pairs give back.
         call free_atom_energies(mol, free_energy, error)
         if (error%has_error()) return
         allocate (adaptation(natm))
         adaptation = tot_intra - free_energy
         formation = sum(adaptation) + 0.5_dp*sum(tot_inter)
         call print_formation(adaptation, free_energy, tot_intra, tot_inter, &
                              element_symbols, formation)

         ! Handed out here rather than at the end of the routine, so a caller
         ! that asked for the decomposition gets it even if the no-sharing CI it
         ! also asked for turns out to be unaffordable.
         if (present(atom_energy)) atom_energy = tot_intra
         if (present(free_atom_energy)) free_atom_energy = free_energy
         if (present(pair_energy)) pair_energy = tot_inter
         if (present(pair_classical)) pair_classical = classical
         if (present(formation_energy)) formation_energy = formation

         ! What the molecule would be if the atoms shared electrons without ever
         ! lending them. Last, because it is a separate calculation rather than
         ! a rearrangement of this one.
         if (present(no_sharing)) then
            if (no_sharing) then
               call clk%lap("energy decomposition")
               call no_sharing_analysis(mol, full_quao, quao, dims%n_core, &
                                        atomic_numbers, natm, &
                                        orbitals(:, 1:dims%n_core), valence_internal, &
                                        ci_route, error, valence_wavefunction)
               if (error%has_error()) return
               call clk%lap("no-sharing wave function")
            end if
         end if

         ! Subtracting a constant from every atomic term cannot change what the
         ! pieces sum to, so this must equal the total less the free atoms. It
         ! catches a free-atom energy landing on the wrong atom, which nothing
         ! else here would notice.
         if (abs(formation - (kinetic_total(tot_intra, tot_inter) - sum(free_energy))) &
             > BALANCE_TOL) then
            call error%set(ERROR_VALIDATION, "the energy of formation does not match "// &
                           "the total less the free atoms, so the atomic references "// &
                           "were not subtracted from the atoms they belong to.")
            return
         end if

         one_electron = kinetic_total(kin_intra, kin_inter) + &
                        kinetic_total(nuc_intra, nuc_inter)
         call mol%core_hamiltonian(h_core)
         call one_electron_reference(orbitals, dims%n_occupied, occupations, &
                                     h_core, density_ao, reference)
         total_energy = one_electron + two_electron + mol%nuclear_repulsion()
         ! The reference built from the orbitals is a determinant expression and
         ! describes a correlated wave function not at all, so for one of those
         ! the energy the calculation reported is what the decomposition is held
         ! against instead. Either way it is a number reached without any of
         ! this arithmetic, which is the whole point of comparing.
         if (correlated .and. .not. present(reference_energy)) then
            call error%set(ERROR_VALIDATION, "a correlated energy decomposition was "// &
                           "asked for without the energy to check it against. The "// &
                           "reference built here assumes a single determinant and "// &
                           "would disagree by the whole correlation energy.")
            return
         end if
         if (present(reference_energy)) then
            reference = reference_energy
         else
            reference = reference + two_electron_reference(density_ao, eri_ao) + &
                        mol%nuclear_repulsion()
         end if

         call logger%info("")
         write (line, "(4x,a30,f16.6,a)") "one-electron", one_electron, " hartree"
         call logger%info(trim(line))
         write (line, "(4x,a30,f16.6,a)") "two-electron", two_electron, " hartree"
         call logger%info(trim(line))
         write (line, "(4x,a30,f16.6,a)") "nuclear repulsion", &
            mol%nuclear_repulsion(), " hartree"
         call logger%info(trim(line))
         write (line, "(4x,a30,f16.6,a)") "total", total_energy, " hartree"
         call logger%info(trim(line))
         if (present(reference_energy)) then
            write (line, "(4x,a30,f16.6,a)") "reported by the calculation", reference, &
               " hartree"
         else
            write (line, "(4x,a30,f16.6,a)") "from the orbitals directly", reference, &
               " hartree"
         end if
         call logger%info(trim(line))
         if (abs(kinetic_total(tot_intra, tot_inter) - total_energy) > BALANCE_TOL) then
            call error%set(ERROR_VALIDATION, "the atom and atom-pair table does not "// &
                           "sum to the total energy, so the four contributions were "// &
                           "not accumulated consistently.")
            return
         end if
         if (abs(total_energy - reference) > BALANCE_TOL) then
            write (line, "(4x,a30,f16.6,a)") "outside the quasi-atomic span", &
               reference - total_energy, " hartree"
            call logger%info(trim(line))
         end if

         ! Two diagnostics saying whether to believe the tables above. The
         ! atomic character is how much of each orbital stayed on its own atom,
         ! one when the orbitals are exactly free-atom ones. The gap is the
         ! separation the valence-virtual selection cut through; **anything
         ! narrow means the valence space was not cleanly separable and the
         ! whole analysis is on sand.**
         call logger%info("")
         write (line, "(a,f8.4,a,i0,a)") "    atomic character   ", &
            quao%atomic_character, "   (", quao%sweeps, " refinement sweeps)"
         call logger%info(trim(line))
         write (line, "(a,f8.4,a,f8.4)") "    valence gap        ", &
            vvo%smallest_retained, "   against rejected ", vvo%largest_rejected
         call logger%info(trim(line))

         ! Per orbital as well as averaged, because the average cannot tell every
         ! orbital being mediocre from a few being poor.
         if (allocated(quao%character_of)) then
            call logger%info("")
            call logger%info("    atomic character per orbital, worst first")
            call sort_ascending(quao%character_of, order)
            do i = 1, min(size(order), 12)
               write (line, "(a,i4,a,a,i0,a,f9.4)") "      orbital ", order(i), &
                  "  on ", trim(adjustl(element_symbols(quao%atom_of(order(i))))), &
                  quao%atom_of(order(i)), "   ", quao%character_of(order(i))
               call logger%info(trim(line))
            end do
            if (size(order) > 12) then
               write (line, "(a,i0,a)") "      (", size(order) - 12, " better ones not shown)"
               call logger%info(trim(line))
            end if
            deallocate (order)
         end if
         call logger%info("")
         deallocate (populations)
      end if

      call clk%finish()
      call clk%report("quasi-atomic analysis", loud)

      call aambs%destroy()
      deallocate (s_mbs, mixed, projection, valence_internal, kinetic, kbo)
      deallocate (interference, core_off, core_n, val_off, val_n)
   end subroutine run_quao_analysis

   subroutine project_density(mol, natural, occupations, valence_internal, &
                              density, error)
      !! `D_vi = W^T diag(n) W` with `W = C^T S V`
      !!
      !! The one-particle density, expressed in the valence-internal basis. `C`
      !! must be natural orbitals -- the density diagonal in that basis with the
      !! given occupations -- because that is the assumption the whole
      !! expression rests on.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: natural(:, :)           !! (n_ao, n_mo)
      real(dp), intent(in) :: occupations(:)          !! (n_mo)
      real(dp), intent(in) :: valence_internal(:, :)  !! (n_ao, n_val)
      real(dp), allocatable, intent(out) :: density(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: overlap(:, :), work(:, :), weights(:, :), scaled(:, :)
      integer :: n_ao, n_mo, n_val, i

      if (error%has_error()) return
      n_ao = size(natural, 1)
      n_mo = size(natural, 2)
      n_val = size(valence_internal, 2)
      if (size(occupations) /= n_mo) then
         call error%set(ERROR_VALIDATION, "there are "//to_char(size(occupations))// &
                        " occupations for "//to_char(n_mo)//" orbitals.")
         return
      end if

      call mol%overlap(overlap)
      allocate (work(n_ao, n_val), weights(n_mo, n_val))
      call pic_gemm(overlap, valence_internal, work)
      call pic_gemm(natural, work, weights, transa="T")

      allocate (scaled(n_mo, n_val), density(n_val, n_val))
      do i = 1, n_mo
         scaled(i, :) = occupations(i)*weights(i, :)
      end do
      call pic_gemm(weights, scaled, density, transa="T")

      deallocate (overlap, work, weights, scaled)
   end subroutine project_density

   subroutine one_electron_reference(orbitals, n_occupied, occupations, h, &
                                     density_ao, energy)
      !! The one-electron energy straight from the converged orbitals
      !!
      !! Built without touching anything quasi-atomic, so that comparing it
      !! against the decomposition is a check rather than a restatement.
      !! `sum(gamma * h)` rather than a trace with a transpose, since both
      !! matrices are symmetric.
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_occupied
      real(dp), intent(in), optional :: occupations(:)
      real(dp), intent(in) :: h(:, :)
      real(dp), allocatable, intent(out) :: density_ao(:, :)
      real(dp), intent(out) :: energy

      integer :: n_ao, i, n_use

      n_ao = size(orbitals, 1)
      allocate (density_ao(n_ao, n_ao))
      density_ao = 0.0_dp

      if (present(occupations)) then
         n_use = min(size(occupations), size(orbitals, 2))
         do i = 1, n_use
            if (abs(occupations(i)) < OCCUPATION_FLOOR) cycle
            density_ao = density_ao + occupations(i)* &
                         spread(orbitals(:, i), 2, n_ao)*spread(orbitals(:, i), 1, n_ao)
         end do
      else
         do i = 1, n_occupied
            density_ao = density_ao + 2.0_dp* &
                         spread(orbitals(:, i), 2, n_ao)*spread(orbitals(:, i), 1, n_ao)
         end do
      end if

      energy = sum(density_ao*h)
   end subroutine one_electron_reference

   pure function two_electron_reference(density_ao, eri_ao) result(energy)
      !! The two-electron energy straight from the atomic-orbital density
      !!
      !! The same `Gamma = gamma gamma - (1/2) gamma gamma` as the decomposition
      !! uses, contracted in the basis the integrals were computed in, and built
      !! independently so that comparing the two is a check.
      real(dp), intent(in) :: density_ao(:, :)
      real(dp), intent(in) :: eri_ao(:, :, :, :)
      real(dp) :: energy

      integer :: n, mu, nu, la, si

      n = size(density_ao, 1)
      energy = 0.0_dp
      do mu = 1, n
         do nu = 1, n
            do la = 1, n
               do si = 1, n
                  energy = energy + 0.5_dp*(density_ao(mu, nu)*density_ao(la, si) &
                                            - 0.5_dp*density_ao(mu, si)*density_ao(la, nu)) &
                           *eri_ao(mu, nu, la, si)
               end do
            end do
         end do
      end do
   end function two_electron_reference

   subroutine print_formation(adaptation, free_energy, intra, inter, &
                              element_symbols, formation)
      !! What it costs to prepare the atoms, against what the pairs give back
      !!
      !! **The adaptation energy is positive and that is not a bug.** An atom in
      !! a molecule is deformed, which costs energy against the free atom;
      !! binding comes from the interatomic terms being more negative than the
      !! adaptation is positive.
      real(dp), intent(in) :: adaptation(:), free_energy(:), intra(:), inter(:, :)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: formation

      character(len=160) :: line
      character(len=16) :: label
      integer :: natm, a

      natm = size(adaptation)
      call logger%info("")
      call logger%info("  energy of formation")
      call logger%info("     atom          in molecule      free atom     adaptation")
      do a = 1, natm
         write (label, "(a,i0)") trim(adjustl(element_symbols(a)))//" ", a
         write (line, "(4x,a10,3f15.6)") label, intra(a), free_energy(a), adaptation(a)
         call logger%info(trim(line))
      end do
      call logger%info("")
      write (line, "(4x,a24,f16.6,a)") "adaptation total", sum(adaptation), " hartree"
      call logger%info(trim(line))
      write (line, "(4x,a24,f16.6,a)") "interatomic total", 0.5_dp*sum(inter), " hartree"
      call logger%info(trim(line))
      write (line, "(4x,a24,f16.6,a,f12.3,a)") "energy of formation", formation, &
         " hartree", formation*HARTREE_TO_KCALMOL, " kcal/mol"
      call logger%info(trim(line))
   end subroutine print_formation

   subroutine adopt_valence_space(offered, dims, n_vvo, n_electrons, &
                                  valence_internal, loud, adopted)
      !! Use the calculation's own valence space, when it has one
      !!
      !! Overwrites `valence_internal` with the offered active orbitals if that
      !! active space is the full valence shell, and leaves it alone otherwise.
      !!
      !! The test is on counts, and counts are all that can be tested: the
      !! valence dimension is a property of the *atoms*, so an active space of
      !! that size holding that many electrons on that many inactive orbitals is
      !! the full valence shell or is a coincidence. GAMESS tests the same thing
      !! at `locsvd.src:3676`.
      !!
      !! A coincidence is not silent: an active space of the right size that is
      !! not the valence shell shows up in the printed `atomic character`.
      type(valence_wavefunction_t), intent(in), optional :: offered
      type(aambs_dimensions_t), intent(in) :: dims
      integer, intent(in) :: n_vvo, n_electrons
      real(dp), intent(inout) :: valence_internal(:, :)
      logical, intent(in) :: loud
      logical, intent(out) :: adopted
         !! Whether the offer was taken. More than a log line: an
         !! occupation-restricted partition is expressed in the active orbitals,
         !! so it means nothing against a valence space derived here instead.

      character(len=160) :: line
      integer :: n_val

      adopted = .false.
      if (.not. present(offered)) return
      if (.not. allocated(offered%orbitals)) return

      n_val = dims%n_valocc + n_vvo
      if (offered%n_active /= n_val) return
      if (offered%n_inactive /= dims%n_core) return
      if (offered%n_alpha + offered%n_beta /= n_electrons - 2*dims%n_core) return
      if (size(offered%orbitals, 2) /= n_val) return
      if (size(offered%orbitals, 1) /= size(valence_internal, 1)) return

      valence_internal = offered%orbitals
      adopted = .true.
      if (loud) then
         write (line, "(a,i0,a,i0,a)") "    valence space               the "// &
            "calculation's own, CAS(", offered%n_alpha + offered%n_beta, ",", &
            n_val, ")"
         call logger%info(trim(line))
      end if
   end subroutine adopt_valence_space

   subroutine inherit_valence_ci(mol, offered, valence_quao, valence_internal, &
                                 n_core, n_active, na, nb, to_quao, declined)
      !! Whether a converged wave function is this one, and the rotation to it
      !!
      !! Returns `to_quao` allocated -- the rotation from the orbitals the
      !! offered CI vector is expanded in to the quasi-atomic ones -- when the
      !! offer can be taken up, and `declined` saying why when it cannot.
      !!
      !! **Nothing here is taken on trust.** A deck's active space is whatever
      !! the deck asked for, and the no-sharing analysis is defined on the full
      !! valence space. Four things have to hold, and the third is the one that
      !! cannot be read off a dimension:
      !!
      !!   1. The vector is a rectangle. A restricted space carries `ci_flat`
      !!      instead, and its determinants are not closed under a rotation of
      !!      the active orbitals in any case.
      !!   2. The counts agree -- same inactive, active and electrons.
      !!   3. The active orbitals **span the same space** as the valence-internal
      !!      set. A CASSCF has moved its orbitals, so equal dimensions prove
      !!      nothing; what proves it is that `<active|S|valence>` comes out
      !!      orthogonal, which two orthonormal bases of the same space give and
      !!      two bases of different spaces do not.
      !!   4. Closed shell, since the projection downstream is.
      !!
      !! Even then the energy is checked afterwards, in the caller: this routine
      !! establishes only that the offer is plausible.
      type(libcint_molecule_t), intent(in) :: mol
      type(valence_wavefunction_t), intent(in) :: offered
      type(quao_result_t), intent(in) :: valence_quao
      real(dp), intent(in) :: valence_internal(:, :)
      integer, intent(in) :: n_core, n_active, na, nb
      real(dp), allocatable, intent(out) :: to_quao(:, :)
      character(len=:), allocatable, intent(out) :: declined

      real(dp), allocatable :: s_ao(:, :), work(:, :), to_valence(:, :)
      integer :: n_val

      declined = ""
      n_val = size(valence_internal, 2)

      if (.not. allocated(offered%ci) .and. .not. allocated(offered%ci_flat)) then
         declined = "it carries no coefficients"
         return
      end if
      if (offered%n_active /= n_active .or. offered%n_inactive /= n_core) then
         declined = "it has "//to_char(offered%n_active)//" active orbitals in "// &
                    to_char(offered%n_inactive)//" inactive, against "// &
                    to_char(n_active)//" in "//to_char(n_core)
         return
      end if
      if (offered%n_alpha /= na .or. offered%n_beta /= nb) then
         declined = "it holds "//to_char(offered%n_alpha)//" alpha and "// &
                    to_char(offered%n_beta)//" beta active electrons, against "// &
                    to_char(na)//" and "//to_char(nb)
         return
      end if
      if (.not. allocated(offered%orbitals)) then
         declined = "the orbitals its coefficients are expanded in did not come with it"
         return
      end if
      if (size(offered%orbitals, 2) /= n_active) then
         declined = "its orbital set is not the size its active space claims"
         return
      end if

      ! <active | S | valence-internal>. Orthogonal exactly when the two
      ! orthonormal sets span the same space, which is the question.
      call mol%overlap(s_ao)
      allocate (work(size(s_ao, 1), n_val), to_valence(n_active, n_val))
      call pic_gemm(s_ao, valence_internal, work)
      call pic_gemm(offered%orbitals, work, to_valence, transa="T")
      deallocate (work, s_ao)

      if (.not. rotation_is_orthogonal(to_valence)) then
         declined = "its active orbitals span a different space from the valence one, "// &
                    "by "//to_char(orthogonality_defect(to_valence))
         return
      end if

      ! Compose: active -> valence-internal -> quasi-atomic.
      allocate (to_quao(n_active, valence_quao%n_quao))
      call pic_gemm(to_valence, valence_quao%to_valence_internal, to_quao)
      deallocate (to_valence)
   end subroutine inherit_valence_ci

   pure subroutine sort_ascending(values, order)
      !! Index permutation putting `values` in ascending order
      !!
      !! Ascending because the interesting end is the low one: a single orbital
      !! that failed to become atomic says more than a dozen that succeeded.
      real(dp), intent(in) :: values(:)
      integer, allocatable, intent(out) :: order(:)

      integer :: n, i, j, pick, hold

      n = size(values)
      allocate (order(n))
      do i = 1, n
         order(i) = i
      end do
      do i = 1, n - 1
         pick = i
         do j = i + 1, n
            if (values(order(j)) < values(order(pick))) pick = j
         end do
         hold = order(i)
         order(i) = order(pick)
         order(pick) = hold
      end do
   end subroutine sort_ascending

   subroutine no_sharing_analysis(mol, full_quao, valence_quao, n_core, &
                                  atomic_numbers, natm, core_orbitals, &
                                  valence_internal, ci_route, error, offered)
      !! The no-sharing wave function, and what charge transfer is worth
      !!
      !! Three steps, of which only the first is expensive. A full valence CI is
      !! needed **expanded over the quasi-atomic orbitals** -- "how many
      !! electrons are on this atom" is a question only an atomic basis can
      !! answer. Its coefficients are then struck out wherever an atom is not
      !! neutral, and the energy of what remains is rebuilt from the projected
      !! density matrices.
      !!
      !! **There are two ways to reach that expansion and they are not equally
      !! good.** A full valence CI is invariant under rotation of its active
      !! orbitals, so the wave function can be had either by solving in the
      !! quasi-atomic basis directly (`NO_SHARING_CI_RESOLVE`) or by solving in
      !! the molecular orbital basis and carrying the vector across with the
      !! orbital transformation (`NO_SHARING_CI_TRANSFORM`, the default).
      !!
      !! Prefer the transform. The two matrices are the same to an orthogonal
      !! similarity, but the Davidson's starting vector and diagonal
      !! preconditioner both assume a dominant configuration the quasi-atomic
      !! basis does not have. The re-solve is kept as an independent route to
      !! the same number, which is the only thing that would catch the
      !! transformation being wrong.
      !!
      !! The difference between the two energies is the charge-transfer
      !! stabilisation. `E(Psi-0)` must come out **above** `E(Psi)`, since a
      !! projection cannot lower a variational energy, and that is asserted.
      type(libcint_molecule_t), intent(in) :: mol
      type(quao_result_t), intent(in) :: full_quao, valence_quao
      integer, intent(in) :: n_core, natm
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: core_orbitals(:, :)
         !! (n_ao, n_core) the molecular orbitals the core quasi-atomic ones were
         !! built from. They span the same space, so the inactive energy and the
         !! mean field are the same either way; they are here so the transform
         !! route can solve in a basis where the reference determinant is a
         !! determinant of the basis.
      real(dp), intent(in) :: valence_internal(:, :)
         !! (n_ao, n_val) the valence-internal orbitals, occupied valence
         !! followed by valence-virtual. `valence_quao%to_valence_internal` is
         !! the transformation from these to the quasi-atomic ones.
      integer, intent(in) :: ci_route
         !! `NO_SHARING_CI_TRANSFORM` or `NO_SHARING_CI_RESOLVE`
      type(error_t), intent(inout) :: error
      type(valence_wavefunction_t), intent(in), optional :: offered
         !! A wave function the calculation already converged. Taken up only if
         !! it is over this same full valence space, and ignored on the resolve
         !! route, whose whole purpose is to arrive independently.

      real(dp), allocatable :: h_eff(:, :), eri_act(:, :, :, :)
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :), ci(:, :)
      real(dp), allocatable :: mo_basis(:, :), string_rotation(:, :)
      type(timing_report_t) :: ns
      real(dp), allocatable :: to_quao(:, :), spread_ci(:, :)
      type(casci_result_t) :: cas
      type(link_table_t) :: alpha, beta
      integer, allocatable :: neutral(:)
      character(len=160) :: line
      character(len=:), allocatable :: declined
      real(dp) :: core_energy, e_psi, e_zero, recovered, e_check
      integer :: n_active, n_valence_electrons, na, nb, iatom, core_orbitals_count
      integer :: valence_orbitals, n_kept, n_total, n_ao
      logical :: inherited

      n_active = valence_quao%n_quao

      ! Neutral means each atom holding the valence electrons it owns.
      allocate (neutral(natm))
      n_valence_electrons = 0
      do iatom = 1, natm
         call aambs_element_counts(atomic_numbers(iatom), core_orbitals_count, &
                                   valence_orbitals, error)
         if (error%has_error()) return
         neutral(iatom) = atomic_numbers(iatom) - 2*core_orbitals_count
         n_valence_electrons = n_valence_electrons + neutral(iatom)
      end do
      if (mod(n_valence_electrons, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "the valence shell holds an odd number of "// &
                        "electrons, which this closed-shell analysis cannot describe.")
         return
      end if
      na = n_valence_electrons/2
      nb = na

      call logger%info("")
      call logger%info("  no-sharing analysis")
      write (line, "(a,i0,a,i0,a)") "     full valence CAS(", n_valence_electrons, &
         ",", n_active, ") over the quasi-atomic orbitals"
      call logger%info(trim(line))

      ! The active-space Hamiltonian in the quasi-atomic basis. Needed on both
      ! routes -- the projected wave function's energy is built against it --
      ! and on the resolve route it is what the CI is solved in as well.
      call ns%start()
      call active_space_integrals(mol, full_quao%orbitals, n_core, n_active, &
                                  h_eff, eri_act, core_energy, error)
      if (error%has_error()) return
      call ns%lap("active space integrals")

      inherited = .false.
      if (ci_route == NO_SHARING_CI_RESOLVE) then
         call logger%info("     CI solved in the quasi-atomic basis")
         call run_libcint_casci(mol, full_quao%orbitals, n_core, n_active, na, nb, &
                                cas, error, verbose=.false.)
         if (error%has_error()) return
         e_psi = cas%energy
         n_total = size(cas%ci_vector)
         ci = cas%ci_vector
      else
         ! Is there already a wave function over this space? A full valence
         ! CASSCF has converged the very thing this analysis is about to
         ! converge again, differing only in which orbitals it is expanded in --
         ! which is a rotation, not a calculation.
         inherited = .false.
         if (present(offered)) then
            call inherit_valence_ci(mol, offered, valence_quao, valence_internal, &
                                    n_core, n_active, na, nb, to_quao, declined)
            inherited = allocated(to_quao)
            if (.not. inherited) then
               ! A warning and not a note: what follows is a decomposition of a
               ! *different* wave function from the one the calculation
               ! reported -- the full valence one this analysis is defined on,
               ! rather than whatever active space the deck asked for. GAMESS
               ! refuses outright at `quao_eda4.src:145`; this one builds its
               ! own wave function, so the result is still meaningful.
               call logger%warning("  the converged wave function is not over the "// &
                                   "full valence space ("//declined//"), so the "// &
                                   "no-sharing analysis solves its own and decomposes "// &
                                   "that one instead")
            end if
         end if

         if (inherited) then
            e_psi = offered%energy
            if (allocated(offered%ci_flat)) then
               ! An occupation-restricted wave function, written out over the
               ! complete space first. It has to be: a restricted space is not
               ! closed under rotation of its active orbitals, so the rotated
               ! wave function has amplitude on determinants the restriction
               ! excluded and there is nowhere in the restricted space to put
               ! them.
               call scatter_restricted_ci(offered%ormas, offered%ci_flat, spread_ci, &
                                          error)
               if (error%has_error()) return
               n_total = size(spread_ci)
               write (line, "(a,i0,a,i0,a)") "     restricted wave function written "// &
                  "out over the complete space (", size(offered%ci_flat), " of ", &
                  n_total, " determinants carried amplitude)"
               call logger%info(trim(line))
            else
               n_total = size(offered%ci)
            end if
         else
            ! Solve where the reference determinant dominates, then re-expand.
            n_ao = size(valence_internal, 1)
            allocate (mo_basis(n_ao, n_core + n_active))
            if (n_core > 0) mo_basis(:, 1:n_core) = core_orbitals
            mo_basis(:, n_core + 1:) = valence_internal
            call run_libcint_casci(mol, mo_basis, n_core, n_active, na, nb, &
                                   cas, error, verbose=.false.)
            if (error%has_error()) return
            e_psi = cas%energy
            n_total = size(cas%ci_vector)
            to_quao = valence_quao%to_valence_internal
         end if

         call string_rotation_matrix(to_quao, na, string_rotation, error)
         if (error%has_error()) return
         call ns%lap("string rotation")
         if (inherited .and. allocated(spread_ci)) then
            call transform_ci_vector(string_rotation, string_rotation, spread_ci, &
                                     ci, error)
         else if (inherited) then
            call transform_ci_vector(string_rotation, string_rotation, offered%ci, &
                                     ci, error)
         else
            call transform_ci_vector(string_rotation, string_rotation, cas%ci_vector, &
                                     ci, error)
         end if
         if (error%has_error()) return
         call ns%lap("transform")
         if (inherited) then
            write (line, "(a,i0,a,i0,a)") "     the converged wave function "// &
               "transformed, no CI solved here (", size(string_rotation, 1), " x ", &
               size(string_rotation, 1), " string rotation)"
         else
            write (line, "(a,i0,a,i0,a)") "     CI solved in the molecular orbital "// &
               "basis and transformed (", size(string_rotation, 1), " x ", &
               size(string_rotation, 1), " string rotation)"
         end if
         call logger%info(trim(line))

         ! The transformed vector is the same state, so its energy against the
         ! quasi-atomic Hamiltonian is the energy it was solved at. This is the
         ! check the whole route rests on: it tests the string minors, their
         ! signs and the orbital transformation together, and nothing else here
         ! would notice any of them being wrong.
         call build_link_table(n_active, na, alpha, error)
         call build_link_table(n_active, nb, beta, error)
         if (error%has_error()) return
         call active_space_rdms(ci, alpha, beta, dm1, dm2, error)
         if (error%has_error()) return
         e_check = rdm_energy(h_eff, eri_act, dm1, dm2) + core_energy
         call ns%lap("invariance check (density matrices)")
         call alpha%destroy()
         call beta%destroy()
         deallocate (dm1, dm2)
         if (abs(e_check - e_psi) > 1.0e-9_dp) then
            ! Concatenated rather than written into `line`, which is 160
            ! characters and would overflow the record on the one path where
            ! being readable matters most.
            call error%set(ERROR_VALIDATION, "the transformed CI vector has a "// &
                           "different energy from the one it was solved at, by "// &
                           to_char(e_check - e_psi)//" hartree. A complete active "// &
                           "space is invariant under rotation of its active "// &
                           "orbitals, so the transformation is wrong.")
            return
         end if
         deallocate (string_rotation, to_quao)
      end if

      if (present(offered) .and. .not. inherited) then
         if (offered%energy /= 0.0_dp) then
            write (line, "(a,f18.12,a,f18.12)") "     decomposing E = ", e_psi, &
               " rather than the calculation's ", offered%energy
            call logger%warning(trim(line))
         end if
      end if

      ! What the route cost. The two solve the same matrix to an orthogonal
      ! similarity, so a difference here is the starting vector and the
      ! preconditioner and nothing else. An inherited wave function has no
      ! iterations of its own to report.
      if (.not. inherited) then
         write (line, "(a,i0,a,i0,a)") "     CI iterations        ", cas%iterations, &
            " (", cas%sigma_products, " sigma products)"
         call logger%info(trim(line))
      end if

      call project_no_sharing(valence_quao%atom_of, natm, neutral, na, nb, ci, &
                              recovered, n_kept, error)
      if (error%has_error()) return
      call ns%lap("neutral projection")

      ! The energy of the projection, from its own density matrices against the
      ! same active-space Hamiltonian the CI used.
      call build_link_table(n_active, na, alpha, error)
      call build_link_table(n_active, nb, beta, error)
      if (error%has_error()) return
      call active_space_rdms(ci, alpha, beta, dm1, dm2, error)
      if (error%has_error()) return
      e_zero = rdm_energy(h_eff, eri_act, dm1, dm2) + core_energy
      call ns%lap("projected density matrices")
      call ns%finish()
      call ns%report("no-sharing", .true.)

      write (line, "(a,i0,a,i0,a,f8.2,a)") "     neutral determinants ", n_kept, &
         " of ", n_total, ", holding ", 100.0_dp*recovered, "% of the squared norm"
      call logger%info(trim(line))
      ! Printed to twelve figures because it is checkable: a full valence CI is
      ! invariant under rotation of its active orbitals, so this must equal the
      ! same CAS run over the ordinary molecular orbitals.
      write (line, "(4x,a24,f18.12,a)") "E(Psi)", e_psi, " hartree"
      call logger%info(trim(line))
      write (line, "(4x,a24,f18.12,a)") "E(Psi-0)", e_zero, " hartree"
      call logger%info(trim(line))
      write (line, "(4x,a24,f16.6,a,f12.3,a)") "charge transfer", e_psi - e_zero, &
         " hartree", (e_psi - e_zero)*HARTREE_TO_KCALMOL, " kcal/mol"
      call logger%info(trim(line))

      call alpha%destroy()
      call beta%destroy()

      ! A projection cannot lower a variational energy. If it appears to, the
      ! coefficients being struck out were not this wave function's.
      if (e_zero < e_psi - 1.0e-10_dp) then
         call error%set(ERROR_VALIDATION, "the no-sharing wave function came out below "// &
                        "the one it was projected from, which a projection cannot do. "// &
                        "The determinants struck out did not belong to this space.")
         return
      end if
   end subroutine no_sharing_analysis

end module mqc_libcint_bonding
