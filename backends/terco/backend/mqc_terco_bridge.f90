!! Real terco bridge, compiled only by CMake with the backend
module mqc_terco_bridge
   !! Forwards to the terco driver. The stub form of this module, in
   !! `src/methods/stubs/mqc_terco_bridge_stub.f90`, is what every other build
   !! compiles; CMake compiles this one instead when MQC_ENABLE_TERCO is ON, so
   !! the method files see a working `run_terco_scf` without any preprocessor
   !! conditionals of their own.
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_terco_driver, only: driver_run_terco_scf => run_terco_scf
   implicit none
   private

   public :: run_terco_scf
   public :: terco_backend_available

contains

   pure function terco_backend_available() result(available)
      !! .true. -- this build has terco
      !!
      !! The stub next door answers .false., so the question is a property of
      !! the binary that linked rather than of the flags somebody believes they
      !! configured with.
      logical :: available

      available = .true.
   end function terco_backend_available

   subroutine run_terco_scf(settings, fragment, result, want_gradient)
      !! Hand the calculation to the terco driver
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient

      call driver_run_terco_scf(settings, fragment, result, want_gradient)
   end subroutine run_terco_scf

end module mqc_terco_bridge
