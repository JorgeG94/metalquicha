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
   use mqc_basis_utils, only: find_basis_file
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_cuest_context, only: cuest_context_t, get_cuest_context
   use mqc_cuest_integrals, only: cuest_system_t
   use mqc_cuest_functionals, only: functional_name_to_id
   use mqc_cuest_scf, only: scf_result_t, run_rhf_scf, run_uks_scf, spin_occupations, &
                            SCF_GUESS_CORE, SCF_GUESS_GWH, SCF_GUESS_SAC
   use mqc_cuest_atomic_guess, only: build_sac_guess
   use mqc_cuest_gradient, only: compute_scf_gradient
   use mqc_cuest_iface, only: cuest_scf_settings_t
   implicit none
   private

   public :: run_cuest_scf         !! Fragment -> converged result

! cuest_scf_settings_t now lives in mqc_cuest_iface (src/), so the method
   ! files can reach it without pulling the backend into the fpm build.

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
      real(dp), allocatable :: guess_alpha(:, :), guess_beta(:, :)
      type(molecular_basis_type), allocatable :: atom_bases(:), atom_aux_bases(:)
      integer :: guess_type, n_guess_columns
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

      call load_basis(settings%basis_set, element_symbols, orbital_basis, error, "orbital")
      if (error%has_error()) then
         call record_failure(result, error)
         return
      end if

      call load_basis(settings%aux_basis_set, element_symbols, auxiliary_basis, error, "auxiliary")
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

      ! ---- which initial guess? ---------------------------------------------
      select case (trim(settings%guess))
      case ("core")
         guess_type = SCF_GUESS_CORE
      case ("sac")
         guess_type = SCF_GUESS_SAC
      case default
         guess_type = SCF_GUESS_GWH
      end select

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

      ! ---- superposed atomic guess, if asked for ----------------------------
      !
      ! Built before the molecular system so its width is known: the guess
      ! carries the summed atomic occupations, usually more columns than the
      ! molecule has electrons, and the device pools must be sized for that up
      ! front rather than grown later under a live system.
      n_guess_columns = 0
      if (guess_type == SCF_GUESS_SAC) then
         call build_atom_bases(settings, element_symbols, atom_bases, atom_aux_bases, error)
         if (.not. error%has_error()) then
            call build_sac_guess(context, fragment%element_numbers, orbital_basis, &
                                 auxiliary_basis, settings%spherical, functional_id, &
                                 settings%radial_points, settings%angular_points, &
                                 atom_bases, atom_aux_bases, guess_alpha, guess_beta, error)
         end if
         if (error%has_error()) then
            call record_failure(result, error)
            return
         end if
         n_guess_columns = size(guess_alpha, 2) + size(guess_beta, 2)
      end if

      ! ---- build, solve, tear down ------------------------------------------
      if (unrestricted) then
         call system%create(context, fragment%element_numbers, fragment%coordinates, &
                            orbital_basis, auxiliary_basis, settings%spherical, &
                            n_alpha, functional_id, settings%radial_points, &
                            settings%angular_points, error, n_occ_beta=n_beta, &
                            n_guess_columns=n_guess_columns)
      else
         call system%create(context, fragment%element_numbers, fragment%coordinates, &
                            orbital_basis, auxiliary_basis, settings%spherical, &
                            fragment%nelec/2, functional_id, settings%radial_points, &
                            settings%angular_points, error, &
                            n_guess_columns=n_guess_columns)
      end if

      if (.not. error%has_error()) then
         if (unrestricted) then
            if (guess_type == SCF_GUESS_SAC) then
               call run_uks_scf(system, fragment%element_numbers, fragment%coordinates, &
                                fragment%nelec, fragment%multiplicity, settings%max_iter, &
                                settings%energy_tol, settings%density_tol, settings%use_diis, &
                                settings%diis_size, settings%verbose, scf, error, &
                                guess=guess_type, guess_alpha=guess_alpha, guess_beta=guess_beta)
            else
               call run_uks_scf(system, fragment%element_numbers, fragment%coordinates, &
                                fragment%nelec, fragment%multiplicity, settings%max_iter, &
                                settings%energy_tol, settings%density_tol, settings%use_diis, &
                                settings%diis_size, settings%verbose, scf, error, guess=guess_type)
            end if
         else
            if (guess_type == SCF_GUESS_SAC) then
               call run_rhf_scf(system, fragment%element_numbers, fragment%coordinates, &
                                fragment%nelec, settings%max_iter, settings%energy_tol, &
                                settings%density_tol, settings%use_diis, settings%diis_size, &
                                settings%verbose, scf, error, guess=guess_type, &
                                guess_alpha=guess_alpha, guess_beta=guess_beta)
            else
               call run_rhf_scf(system, fragment%element_numbers, fragment%coordinates, &
                                fragment%nelec, settings%max_iter, settings%energy_tol, &
                                settings%density_tol, settings%use_diis, settings%diis_size, &
                                settings%verbose, scf, error, guess=guess_type)
            end if
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

      ! The dipole is what IR intensities are built from: the Hessian path
      ! collects it at every displacement and differences it.
      if (scf%has_dipole) then
         result%dipole = scf%dipole
         result%has_dipole = .true.
      end if

      if (need_gradient) then
         call move_alloc(gradient, result%gradient)
         result%has_gradient = .true.
      end if
   end subroutine run_cuest_scf

   subroutine build_atom_bases(settings, element_symbols, atom_bases, atom_aux_bases, error)
      !! One-atom orbital and auxiliary bases, for the free-atom guess runs
      type(cuest_scf_settings_t), intent(in) :: settings
      character(len=*), intent(in) :: element_symbols(:)
      type(molecular_basis_type), allocatable, intent(out) :: atom_bases(:), atom_aux_bases(:)
      type(error_t), intent(inout) :: error

      integer :: iatom, n_atoms

      n_atoms = size(element_symbols)
      allocate (atom_bases(n_atoms), atom_aux_bases(n_atoms))
      do iatom = 1, n_atoms
         call load_basis(settings%basis_set, element_symbols(iatom:iatom), &
                         atom_bases(iatom), error)
         if (error%has_error()) return
         call load_basis(settings%aux_basis_set, element_symbols(iatom:iatom), &
                         atom_aux_bases(iatom), error)
         if (error%has_error()) return
      end do
   end subroutine build_atom_bases

   subroutine load_basis(basis_name, element_symbols, mol_basis, error, role)
      !! Locate and parse a basis set for the atoms of a fragment
      !!
      !! Basis Set Exchange JSON is the only format; `find_basis_file` walks
      !! the search path and reports where it looked if nothing turns up.
      character(len=*), intent(in) :: basis_name
      character(len=*), intent(in) :: element_symbols(:)
      type(molecular_basis_type), intent(out) :: mol_basis
      type(error_t), intent(out) :: error
      character(len=*), intent(in), optional :: role  !! "orbital" or "auxiliary", for diagnostics

      character(len=:), allocatable :: basis_path

      if (len_trim(basis_name) == 0) then
         call error%set(ERROR_VALIDATION, "No basis set specified")
         return
      end if

      call find_basis_file(basis_name, basis_path, error)
      if (error%has_error()) return

      call build_molecular_basis_json(basis_path, element_symbols, mol_basis, error)
      if (error%has_error()) return

      call check_basis_covers_atoms(basis_name, basis_path, element_symbols, mol_basis, error, role)
   end subroutine load_basis

   subroutine check_basis_covers_atoms(basis_name, basis_path, element_symbols, mol_basis, error, role)
      !! Fail unless every atom came back with at least one shell
      !!
      !! The readers already refuse an element that is absent from the file,
      !! but each of the three has its own notion of "absent" and a block that
      !! parses to nothing would slip through all of them. An atom with no
      !! basis functions is not a small error: cuEST would build a system whose
      !! AO space simply omits that centre and converge an SCF for a different
      !! molecule than the one that was asked for. Refuse it here, once, on the
      !! only path any of the readers can reach the backend by.
      character(len=*), intent(in) :: basis_name, basis_path
      character(len=*), intent(in) :: element_symbols(:)
      type(molecular_basis_type), intent(in) :: mol_basis
      type(error_t), intent(inout) :: error
      character(len=*), intent(in), optional :: role

      character(len=:), allocatable :: what
      integer :: iatom

      what = "basis set"
      if (present(role)) what = trim(role)//" basis set"

      if (mol_basis%nelements /= size(element_symbols)) then
         call error%set(ERROR_VALIDATION, "The "//what//" "//trim(basis_name)// &
                        " produced entries for a different number of atoms than the fragment has")
         return
      end if

      do iatom = 1, size(element_symbols)
         if (mol_basis%elements(iatom)%nshells <= 0) then
            call error%set(ERROR_VALIDATION, "No "//what//" defined for element "// &
                           trim(adjustl(element_symbols(iatom)))//" (atom "//index_text(iatom)// &
                           ") in "//trim(basis_path)//". Add the element to the basis file or "// &
                           "choose a basis that covers it.")
            return
         end if
      end do
   end subroutine check_basis_covers_atoms

   pure function index_text(value) result(text)
      !! An integer as a trimmed string, for error messages
      integer, intent(in) :: value
      character(len=:), allocatable :: text
      character(len=12) :: buffer

      write (buffer, "(I0)") value
      text = trim(buffer)
   end function index_text

   subroutine record_failure(result, error)
      !! Mark a calculation as failed, carrying the diagnostic with it
      type(calculation_result_t), intent(inout) :: result
      type(error_t), intent(in) :: error

      result%error = error
      result%has_error = .true.
      result%has_energy = .false.
   end subroutine record_failure

end module mqc_cuest_driver
