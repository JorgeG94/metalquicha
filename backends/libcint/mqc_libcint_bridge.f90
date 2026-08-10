!! One fragment into the CPU SCF and back out
module mqc_libcint_bridge
   !! Turns a `physical_fragment_t` into a libcint molecule, runs RHF over it,
   !! and fills in the result the rest of the program reads.
   !!
   !! The mirror of `mqc_cuest_bridge`, and the same shape on purpose: one
   !! entry point taking the settings the method already assembled, so which
   !! backend runs is a choice made in one place rather than a difference
   !! spread through the method layer.
   !!
   !! Caps are ordinary atoms here. A capped fragment arrives with its cap
   !! hydrogens already in `element_numbers` and `coordinates`, and nothing
   !! about the integrals cares which atoms were in the original molecule --
   !! the redistribution of forces back onto the heavy atoms happens later and
   !! elsewhere.
   use pic_types, only: dp
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t, SCF_CONVERGED, SCF_NOT_CONVERGED
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_number_to_symbol
   use mqc_program_limits, only: MAX_ELEMENT_SYMBOL_LEN
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   implicit none
   private

   public :: run_libcint_hf
   public :: libcint_backend_available

contains

   pure function libcint_backend_available() result(available)
      !! Whether this build can run an SCF on the CPU
      logical :: available
      available = .true.
   end function libcint_backend_available

   subroutine run_libcint_hf(settings, fragment, result, want_gradient)
      !! Closed-shell HF for one fragment, on the CPU
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: error
      character(len=MAX_ELEMENT_SYMBOL_LEN), allocatable :: symbols(:)
      integer :: iatom, diis_size

      if (present(want_gradient)) then
         if (want_gradient) then
            ! libcint has the derivative entry points -- ip1 and ipovlp -- so
            ! this is a gap rather than a wall, but an untested gradient is
            ! worse than an absent one.
            call result%error%set(ERROR_VALIDATION, "the CPU backend has no gradients "// &
                                  "yet; run an energy, or build with cuEST")
            result%has_error = .true.
            return
         end if
      end if

      if (settings%unrestricted) then
         call result%error%set(ERROR_VALIDATION, "the CPU backend is closed-shell only; "// &
                               "this fragment asks for an unrestricted calculation")
         result%has_error = .true.
         return
      end if

      allocate (symbols(fragment%n_atoms))
      do iatom = 1, fragment%n_atoms
         symbols(iatom) = element_number_to_symbol(fragment%element_numbers(iatom))
      end do

      call build_libcint_molecule(fragment%element_numbers, symbols, &
                                  fragment%coordinates, trim(settings%basis_set), &
                                  mol, error)
      if (error%has_error()) then
         call result%error%set(ERROR_VALIDATION, error%get_message())
         result%has_error = .true.
         return
      end if

      diis_size = settings%diis_size
      if (.not. settings%use_diis) diis_size = 0

      call run_libcint_rhf(mol, fragment%nelec, settings%max_iter, settings%energy_tol, &
                           settings%density_tol, settings%verbose, scf, error, &
                           diis_vectors=diis_size)
      if (error%has_error()) then
         call result%error%set(ERROR_VALIDATION, error%get_message())
         result%has_error = .true.
         call mol%destroy()
         return
      end if

      ! Same policy as every other backend: a non-converged SCF is refused
      ! unless the run asked to keep it, because its energy has the right
      ! magnitude and nothing downstream can tell.
      if (scf%converged) then
         result%scf_status = SCF_CONVERGED
      else
         result%scf_status = SCF_NOT_CONVERGED
      end if
      result%scf_iterations = scf%iterations
      if (.not. scf%converged .and. .not. settings%allow_crap_scf) then
         call result%error%set(ERROR_VALIDATION, "SCF not converged in "// &
                               to_text(scf%iterations)//" cycles")
         result%has_error = .true.
         call mol%destroy()
         return
      end if

      result%energy%scf = scf%energy
      result%has_energy = .true.

      ! The frontier pair, from the orbital energies the SCF already produced.
      if (allocated(scf%orbital_energies)) then
         if (scf%n_occupied >= 1 .and. scf%n_occupied < size(scf%orbital_energies)) then
            result%homo = scf%orbital_energies(scf%n_occupied)
            result%lumo = scf%orbital_energies(scf%n_occupied + 1)
            result%has_orbitals = .true.
         end if
      end if

      call mol%destroy()
   end subroutine run_libcint_hf

   pure function to_text(value) result(out)
      integer, intent(in) :: value
      character(len=:), allocatable :: out
      character(len=16) :: buffer
      write (buffer, "(i0)") value
      out = trim(adjustl(buffer))
   end function to_text

end module mqc_libcint_bridge
