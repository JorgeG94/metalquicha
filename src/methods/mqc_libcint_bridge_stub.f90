!! Stand-in for the CPU SCF when the build has no libcint
module mqc_libcint_bridge
   !! Same name and same entry points as the real bridge, declining.
   !!
   !! A stub rather than a preprocessor guard at the call site, matching
   !! mqc_cuest_bridge_stub next door: the method layer reads the same either
   !! way, and `libcint_backend_available` is how it asks in advance so the
   !! refusal can name the build option.
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_error, only: ERROR_VALIDATION
   use mqc_cuest_iface, only: cuest_scf_settings_t
   implicit none
   private

   public :: run_libcint_hf
   public :: run_libcint_makefp
   public :: run_libcint_efp
   public :: libcint_backend_available

contains

   pure function libcint_backend_available() result(available)
      !! Whether this build can run an SCF on the CPU
      logical :: available
      available = .false.
   end function libcint_backend_available

   subroutine run_libcint_efp(potentials, fragment_sizes, fragment_atoms, &
                              coordinates, terms, error)
      !! No-op stand-in: an EFP interaction energy needs the CPU backend
      !!
      !! Every piece of it does -- the multipole and polarizability machinery, the
      !! integrals two fragments' basis sets share, and the potential reader itself
      !! all live behind `MQC_ENABLE_LIBCINT`.
      use pic_types, only: dp
      use mqc_program_limits, only: N_EFP_TERMS
      use mqc_error, only: error_t
      character(len=*), intent(in) :: potentials(:)
      integer, intent(in) :: fragment_sizes(:)
      integer, intent(in) :: fragment_atoms(:, :)
      real(dp), intent(in) :: coordinates(:, :)
      real(dp), intent(out) :: terms(N_EFP_TERMS)
      type(error_t), intent(inout) :: error

      terms = 0.0_dp
      call error%set(ERROR_VALIDATION, &
                     "EFP needs the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_LIBCINT=ON")
      if (len_trim(potentials(1)) < 0) return
      if (size(fragment_sizes) < 0 .or. size(fragment_atoms) < 0) return
      if (size(coordinates) < 0) return
   end subroutine run_libcint_efp

   subroutine run_libcint_makefp(atomic_numbers, element_symbols, coordinates, &
                                 basis_name, name, path, error, charge, verbose, &
                                 aux_basis)
      !! No-op stand-in: an effective fragment potential needs the CPU backend
      use pic_types, only: dp
      use mqc_error, only: error_t
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      character(len=*), intent(in) :: basis_name, name, path
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: charge   !! Net charge; a fragment may be an ion
      logical, intent(in), optional :: verbose
      character(len=*), intent(in), optional :: aux_basis

      call error%set(ERROR_VALIDATION, &
                     "MAKEFP needs the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_LIBCINT=ON")
      if (size(atomic_numbers) < 0 .or. size(coordinates) < 0) return
      if (len_trim(element_symbols(1)) < 0) return
      if (len_trim(basis_name)*len_trim(name)*len_trim(path) < 0) return
      if (present(charge) .or. present(verbose)) return
   end subroutine run_libcint_makefp

   subroutine run_libcint_hf(settings, fragment, result, want_gradient)
      !! No-op stand-in: report the missing backend, compute nothing
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient

      call result%error%set(ERROR_VALIDATION, &
                            "This calculation needs an integral backend; build with "// &
                            "-DMQC_ENABLE_LIBCINT=ON for the CPU one, or "// &
                            "-DMQC_ENABLE_CUEST=ON for the GPU one")
      result%has_error = .true.
      result%has_energy = .false.
      if (len_trim(settings%basis_set) == 0 .or. fragment%n_atoms < 0) return
      if (present(want_gradient)) return
   end subroutine run_libcint_hf

end module mqc_libcint_bridge
