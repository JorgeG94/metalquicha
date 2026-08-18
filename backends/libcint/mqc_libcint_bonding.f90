!! Bonding analysis on a converged wave function
module mqc_libcint_bonding
   !! One entry point that takes a converged set of molecular orbitals and
   !! produces the quasi-atomic bonding picture: atomic populations and charges,
   !! then every orbital pair worth reporting with its bond order, its kinetic
   !! bond order and what kind of orbital each end is.
   !!
   !! The pipeline itself lives in `mqc_aambs`, `mqc_libcint_quao` and
   !! `mqc_libcint_quao_report`, each validated on its own. What is here is only
   !! the assembly -- which spaces to build, in what order, and what to hand to
   !! the report -- so that a deck can ask for the analysis without a caller
   !! having to know that the valence-internal space is the occupied valence
   !! orbitals followed by the valence-virtual ones.
   !!
   !! **This is a Hartree-Fock analysis for now.** The construction takes a
   !! one-particle density matrix and does not care where it came from, and
   !! `quasi_atomic_orbitals` accepts one explicitly for exactly that reason.
   !! But nothing here supplies a correlated density yet, so what is analysed is
   !! the reference determinant.
   use pic_types, only: dp
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
   use mqc_libcint_quao_report, only: quao_labels_t, label_quasi_atomic_orbitals, &
                                      print_quao_report
   implicit none
   private

   public :: run_quao_analysis
   public :: BONDING_NONE, BONDING_GMS_QUAO
   public :: bonding_analysis_kind
   public :: bonding_analysis_name

   integer, parameter :: BONDING_NONE = 0
   integer, parameter :: BONDING_GMS_QUAO = 1
      !! The Ruedenberg quasi-atomic analysis as GAMESS implements it. Named for
      !! the reference implementation rather than the papers because the labels
      !! and the thresholds are GAMESS's, and only the underlying quantities are
      !! the papers'.

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

   subroutine run_quao_analysis(mol, atomic_numbers, element_symbols, coordinates, &
                                orbitals, n_electrons, error, verbose)
      !! The quasi-atomic bonding analysis, start to finish
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)     !! (3, natm), Bohr
      real(dp), intent(in) :: orbitals(:, :)        !! Converged, (n_ao, n_mo)
      integer, intent(in) :: n_electrons
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose

      type(libcint_molecule_t) :: aambs
      type(aambs_dimensions_t) :: dims
      type(vvo_result_t) :: vvo
      type(quao_result_t) :: quao
      type(quao_labels_t) :: labels
      real(dp), allocatable :: s_mbs(:, :), mixed(:, :), projection(:, :)
      real(dp), allocatable :: valence_internal(:, :), kinetic(:, :)
      real(dp), allocatable :: kbo(:, :), interference(:, :), populations(:)
      integer, allocatable :: core_off(:), core_n(:), val_off(:), val_n(:)
      character(len=160) :: line
      integer :: natm, iatom, i, core, valence
      logical :: loud

      if (error%has_error()) return
      loud = .true.
      if (present(verbose)) loud = verbose
      natm = size(atomic_numbers)

      ! The dimensions first: they are pure counting and they refuse the cases
      ! this cannot handle -- an element past xenon, an odd electron count --
      ! before any integral is built.
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

      call quasi_atomic_orbitals(atomic_numbers, valence_internal, dims%n_valocc, &
                                 mixed, val_off, val_n, quao, error)
      call orient_quasi_atomic_orbitals(quao, error)
      if (error%has_error()) return

      call mol%kinetic(kinetic)
      call kinetic_bond_orders(quao, kinetic, kbo, error, interference)
      call label_quasi_atomic_orbitals(quao, mol, aambs, mixed, s_mbs, &
                                       atomic_numbers, coordinates, labels, error)
      if (error%has_error()) return

      if (loud) then
         ! Populations first: the headline of the analysis, and the number a
         ! reader checks against chemical intuition before looking at anything
         ! else. The core electrons are added back because the construction only
         ! ever sees the valence.
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
         call logger%info("  quasi-atomic bonding analysis")
         call logger%info("")
         call logger%info("    atom      population        charge")
         do iatom = 1, natm
            write (line, "(a,i4,1x,a2,2f16.7)") "    ", iatom, &
               adjustl(element_symbols(iatom)), populations(iatom), &
               real(atomic_numbers(iatom), dp) - populations(iatom)
            call logger%info(trim(line))
         end do
         write (line, "(a,f16.7)") "    total                 ", sum(populations)
         call logger%info(trim(line))

         write (line, "(a,f10.6,a,i0,a)") "    atomic character  ", &
            quao%atomic_character, " after ", quao%sweeps, " refinement sweeps"
         call logger%info("")
         call logger%info(trim(line))
         write (line, "(a,f10.6,a,f10.6)") "    valence-virtual gap: retained ", &
            vvo%smallest_retained, " against rejected ", vvo%largest_rejected
         call logger%info(trim(line))

         call print_quao_report(.true., quao, labels, interference, element_symbols, &
                                dims%n_core)
         deallocate (populations)
      end if

      call aambs%destroy()
      deallocate (s_mbs, mixed, projection, valence_internal, kinetic, kbo)
      deallocate (interference, core_off, core_n, val_off, val_n)
   end subroutine run_quao_analysis

end module mqc_libcint_bonding
