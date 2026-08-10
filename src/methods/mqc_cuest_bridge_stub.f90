!! Stub cuEST bridge, for builds without the cuEST backend
module mqc_cuest_bridge
   !! The fpm build, and any CMake build with MQC_ENABLE_CUEST=OFF, compiles
   !! this. It reports that no integral backend is present rather than
   !! computing anything. The real implementation, compiled only by CMake with
   !! the backend, lives in `backends/cuest/backend/mqc_cuest_bridge.f90` and
   !! provides the same module and the same `run_cuest_scf` interface.
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_error, only: ERROR_VALIDATION
   implicit none
   private

   public :: run_cuest_scf

contains

   subroutine run_cuest_scf(settings, fragment, result, want_gradient)
      !! No-op stand-in: report the missing backend, compute nothing
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient

      call result%error%set(ERROR_VALIDATION, &
                            "This calculation needs the cuEST integral backend; "// &
                            "build with CMake and -DMQC_ENABLE_CUEST=ON")
      result%has_error = .true.
      result%has_energy = .false.
   end subroutine run_cuest_scf

end module mqc_cuest_bridge
