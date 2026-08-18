!! Real cuEST bridge, compiled only by CMake with the backend
module mqc_cuest_bridge
   !! Forwards to the cuEST driver. The stub form of this module, in
   !! `src/methods/mqc_cuest_bridge_stub.f90`, is what the fpm build compiles;
   !! CMake compiles this one instead when MQC_ENABLE_CUEST is ON, so the
   !! method files see a working `run_cuest_scf` without any preprocessor
   !! conditionals of their own.
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_cuest_driver, only: driver_run_cuest_scf => run_cuest_scf
   implicit none
   private

   public :: run_cuest_scf
   public :: cuest_backend_available

contains

   pure function cuest_backend_available() result(available)
      !! .true. -- this build has cuEST
      !!
      !! The stub next door answers .false., so the question is a property of
      !! the binary that linked rather than of the flags somebody believes they
      !! configured with.
      logical :: available

      available = .true.
   end function cuest_backend_available

   subroutine run_cuest_scf(settings, fragment, result, want_gradient)
      !! Hand the calculation to the cuEST driver
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient

      call driver_run_cuest_scf(settings, fragment, result, want_gradient)
   end subroutine run_cuest_scf

end module mqc_cuest_bridge
