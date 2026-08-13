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
   public :: libcint_backend_available

contains

   pure function libcint_backend_available() result(available)
      !! Whether this build can run an SCF on the CPU
      logical :: available
      available = .false.
   end function libcint_backend_available

   subroutine run_libcint_makefp(atomic_numbers, element_symbols, coordinates, &
                                 basis_name, name, path, error, verbose)
      !! No-op stand-in: an effective fragment potential needs the CPU backend
      use pic_types, only: dp
      use mqc_error, only: error_t
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      character(len=*), intent(in) :: basis_name, name, path
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose

      call error%set(ERROR_VALIDATION, &
                     "MAKEFP needs the CPU integral backend; build with "// &
                     "-DMQC_ENABLE_LIBCINT=ON")
      if (size(atomic_numbers) < 0 .or. size(coordinates) < 0) return
      if (len_trim(element_symbols(1)) < 0) return
      if (len_trim(basis_name)*len_trim(name)*len_trim(path) < 0) return
      if (present(verbose)) return
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
