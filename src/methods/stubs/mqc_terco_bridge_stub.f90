!! Stub terco bridge, for builds without the terco backend
module mqc_terco_bridge
   !! The fpm build, and any CMake build with MQC_ENABLE_TERCO=OFF, compiles
   !! this. It reports that the backend was asked for and is not present rather
   !! than computing anything. The real implementation, compiled only by CMake
   !! with the backend, lives in `backends/terco/backend/mqc_terco_bridge.f90`
   !! and provides the same module and the same `run_terco_scf` interface.
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_error, only: ERROR_VALIDATION
   implicit none
   private

   public :: run_terco_scf
   public :: terco_backend_available

contains

   pure function terco_backend_available() result(available)
      !! .false. -- this build has no terco
      !!
      !! Asked before a calculation starts, so a deck that asks for terco on a
      !! build without it is refused while it is still a deck rather than after
      !! a fragment has been set up and handed to a bridge that computes
      !! nothing.
      logical :: available

      available = .false.
   end function terco_backend_available

   subroutine run_terco_scf(settings, fragment, result, want_gradient)
      !! No-op stand-in: report the missing backend, compute nothing
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient

      call result%error%set(ERROR_VALIDATION, &
                            "This calculation needs the terco integral backend; "// &
                            "build with CMake and -DMQC_ENABLE_TERCO=ON")
      result%has_error = .true.
      result%has_energy = .false.
   end subroutine run_terco_scf

end module mqc_terco_bridge
