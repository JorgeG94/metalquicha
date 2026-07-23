!! Fragment-level entry point for cuEST-backed SCF calculations
module mqc_cuest_driver
   !! One place where a `physical_fragment_t` becomes a converged SCF result.
   !!
   !! Hartree-Fock and Kohn-Sham differ only in whether a functional is named,
   !! so both methods funnel through here rather than each carrying its own
   !! copy of the basis loading, context acquisition and teardown.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cgto, only: molecular_basis_type
   use mqc_elements, only: element_number_to_symbol
   use mqc_basis_utils, only: find_basis_file, BASIS_FORMAT_GBS, BASIS_FORMAT_GAMESS
   use mqc_basis_reader, only: build_molecular_basis
   use mqc_basis_file_reader, only: basis_file_t, open_basis_file
   use mqc_gbs_reader, only: build_molecular_basis_gbs
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_cuest_context, only: cuest_context_t, get_cuest_context
   use mqc_cuest_integrals, only: cuest_system_t
   use mqc_cuest_functionals, only: functional_name_to_id
   use mqc_cuest_scf, only: scf_result_t, run_rhf_scf, run_uks_scf, spin_occupations
   use mqc_cuest_gradient, only: compute_scf_gradient
   implicit none
   private

   public :: cuest_scf_settings_t  !! Everything a cuEST SCF needs to run
   public :: run_cuest_scf         !! Fragment -> converged result

   type :: cuest_scf_settings_t
      !! Method-independent description of one cuEST SCF calculation
      character(len=32) :: basis_set = 'sto-3g'
         !! Orbital basis set name
      character(len=32) :: aux_basis_set = 'def2-universal-jkfit'
         !! Auxiliary (JKFIT) basis. Required: cuEST fits J and K always.
      character(len=32) :: functional = ''
         !! Exchange-correlation functional; empty means Hartree-Fock
      logical :: spherical = .true.
         !! Pure (spherical) vs Cartesian angular functions
      logical :: verbose = .false.
         !! Print the SCF iteration table
      integer :: device_rank = 0
         !! Node-local MPI rank. Decides which GPU this rank binds to; leaving
         !! it at zero in a multi-rank run puts every rank on device 0.
      logical :: unrestricted = .false.
         !! Force UHF/UKS even for a closed shell. With equal alpha and beta
         !! occupations and the same guess for both spins, this must reproduce
         !! the restricted energy exactly, which is how the unrestricted
         !! factors are checked.

      integer :: max_iter = 100
      real(dp) :: energy_tol = 1.0e-8_dp
      real(dp) :: density_tol = 1.0e-6_dp
      logical :: use_diis = .true.
      integer :: diis_size = 8

      integer :: radial_points = 75    !! XC grid radial points per atom
      integer :: angular_points = 302  !! XC grid Lebedev order
   end type cuest_scf_settings_t

contains

   subroutine run_cuest_scf(settings, fragment, result, want_gradient)
      !! Run a closed-shell cuEST SCF for one fragment
      !!
      !! With `want_gradient`, the analytic nuclear gradient is evaluated from
      !! the converged density before the cuEST objects are torn down -- they
      !! describe this geometry and this basis, so reusing them is both
      !! cheaper and safer than rebuilding.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient

      type(cuest_context_t), pointer :: context
      type(cuest_system_t) :: system
      type(scf_result_t) :: scf
      type(molecular_basis_type) :: orbital_basis, auxiliary_basis
      type(error_t) :: error
      real(dp), allocatable :: gradient(:, :)
      character(len=8), allocatable :: element_symbols(:)
      integer :: iatom, functional_id, n_alpha, n_beta
      logical :: need_gradient, unrestricted, occupations_ok

      need_gradient = .false.
      if (present(want_gradient)) need_gradient = want_gradient

      ! ---- which functional, if any ----------------------------------------
      if (len_trim(settings%functional) == 0) then
         functional_id = -1   ! pure Hartree-Fock, no XC plan
      else
         call functional_name_to_id(settings%functional, functional_id, error)
         if (error%has_error()) then
            call record_failure(result, error)
            return
         end if
      end if

      ! ---- basis sets -------------------------------------------------------
      allocate (element_symbols(fragment%n_atoms))
      do iatom = 1, fragment%n_atoms
         element_symbols(iatom) = element_number_to_symbol(fragment%element_numbers(iatom))
      end do

      call load_basis(settings%basis_set, element_symbols, orbital_basis, error)
      if (error%has_error()) then
         call record_failure(result, error)
         return
      end if

      call load_basis(settings%aux_basis_set, element_symbols, auxiliary_basis, error)
      if (error%has_error()) then
         call record_failure(result, error)
         return
      end if

      ! ---- per-rank handle and device scratch -------------------------------
      !
      ! One handle and one set of scratch pools per rank, reused by every
      ! fragment this rank evaluates. The node-local rank spreads ranks across
      ! the GPUs on the node; it is ignored after the first call, since the
      ! context is created once and then shared.
      call get_cuest_context(settings%device_rank, context, error)
      if (error%has_error()) then
         call record_failure(result, error)
         return
      end if

      ! ---- restricted or unrestricted? --------------------------------------
      !
      ! An odd electron count forces unrestricted regardless of what the input
      ! claims, so a mis-stated multiplicity fails as a validation error rather
      ! than silently halving an odd number of electrons.
      call spin_occupations(fragment%nelec, fragment%multiplicity, n_alpha, n_beta, occupations_ok)
      if (.not. occupations_ok) then
         call error%set(ERROR_VALIDATION, "Electron count and multiplicity are inconsistent")
         call record_failure(result, error)
         return
      end if
      unrestricted = (fragment%multiplicity /= 1) .or. (mod(fragment%nelec, 2) /= 0) &
                     .or. settings%unrestricted

      ! ---- build, solve, tear down ------------------------------------------
      if (unrestricted) then
         call system%create(context, fragment%element_numbers, fragment%coordinates, &
                            orbital_basis, auxiliary_basis, settings%spherical, &
                            n_alpha, functional_id, settings%radial_points, &
                            settings%angular_points, error, n_occ_beta=n_beta)
      else
         call system%create(context, fragment%element_numbers, fragment%coordinates, &
                            orbital_basis, auxiliary_basis, settings%spherical, &
                            fragment%nelec/2, functional_id, settings%radial_points, &
                            settings%angular_points, error)
      end if

      if (.not. error%has_error()) then
         if (unrestricted) then
            call run_uks_scf(system, fragment%element_numbers, fragment%coordinates, &
                             fragment%nelec, fragment%multiplicity, settings%max_iter, &
                             settings%energy_tol, settings%density_tol, settings%use_diis, &
                             settings%diis_size, settings%verbose, scf, error)
         else
            call run_rhf_scf(system, fragment%element_numbers, fragment%coordinates, &
                             fragment%nelec, settings%max_iter, settings%energy_tol, &
                             settings%density_tol, settings%use_diis, settings%diis_size, &
                             settings%verbose, scf, error)
         end if
      end if

      if (need_gradient .and. .not. error%has_error()) then
         call compute_scf_gradient(system, fragment%element_numbers, fragment%coordinates, &
                                   scf, gradient, error)
      end if

      call system%destroy()
      call orbital_basis%destroy()
      call auxiliary_basis%destroy()

      if (error%has_error()) then
         call record_failure(result, error)
         return
      end if

      result%energy%scf = scf%total_energy
      result%has_energy = .true.

      if (need_gradient) then
         call move_alloc(gradient, result%gradient)
         result%has_gradient = .true.
      end if
   end subroutine run_cuest_scf

   subroutine load_basis(basis_name, element_symbols, mol_basis, error)
      !! Locate and parse a basis set for the atoms of a fragment
      !!
      !! `.gbs` (Gaussian94) is preferred and `.txt` (GAMESS $DATA) accepted;
      !! `find_basis_file` decides which by what is actually on disk.
      character(len=*), intent(in) :: basis_name
      character(len=*), intent(in) :: element_symbols(:)
      type(molecular_basis_type), intent(out) :: mol_basis
      type(error_t), intent(out) :: error

      character(len=:), allocatable :: basis_path
      type(basis_file_t) :: basis_file
      integer :: basis_format

      if (len_trim(basis_name) == 0) then
         call error%set(ERROR_VALIDATION, "No basis set specified")
         return
      end if

      call find_basis_file(basis_name, basis_path, error, basis_format)
      if (error%has_error()) return

      select case (basis_format)
      case (BASIS_FORMAT_GBS)
         call build_molecular_basis_gbs(basis_path, element_symbols, mol_basis, error)
      case (BASIS_FORMAT_GAMESS)
         call open_basis_file(basis_file, basis_path, error)
         if (error%has_error()) return
         call build_molecular_basis(basis_file%data_section, element_symbols, mol_basis, error)
      case default
         call error%set(ERROR_VALIDATION, "Unrecognized basis file format for "//trim(basis_name))
      end select
   end subroutine load_basis

   subroutine record_failure(result, error)
      !! Mark a calculation as failed, carrying the diagnostic with it
      type(calculation_result_t), intent(inout) :: result
      type(error_t), intent(in) :: error

      result%error = error
      result%has_error = .true.
      result%has_energy = .false.
   end subroutine record_failure

end module mqc_cuest_driver
