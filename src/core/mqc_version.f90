!! Version information for metalquicha
module mqc_version
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: MQC_VERSION_STR, print_version

   character(len=*), parameter :: MQC_VERSION_STR = "0.2.0"
      !! Kept in step with `VERSION` in the top-level CMakeLists by hand.
      !! Nothing checks that they agree, and they had already drifted once.

contains

   subroutine print_version()
      !! Version, and which optional backends this binary actually has
      !!
      !! The feature list is here so that something outside the program can
      !! find out what this build can do without trying it and reading the
      !! failure. `run_validation.py` gates on it: a deck needing an optional
      !! backend is skipped rather than failed on a build configured without
      !! one, which is the difference between "not tested here" and "broken".
      !!
      !! Asked of the backends themselves rather than of preprocessor symbols.
      !! Each of these resolves to a real module or to its stub at link time,
      !! so the answer is a property of the binary that exists rather than of
      !! the flags somebody believes they configured with.
      use mqc_libcint_bridge, only: libcint_backend_available, xc_available, ecp_backend_available
      use mqc_method_factory, only: method_backend_built
      use mqc_method_types, only: METHOD_TYPE_GFN2
      use mqc_cuest_bridge, only: cuest_backend_available
      use mqc_dlfind_bridge, only: dlfind_available

      call logger%info("metalquicha version "//MQC_VERSION_STR)
      ! tblite has no bridge of its own to ask, so the question goes to the
      ! factory, which is preprocessed and knows whether the xTB method type
      ! has a branch to reach.
      call logger%info("features: libcint="//available(libcint_backend_available())// &
                       " libxc="//available(xc_available())// &
                       " tblite="//available(method_backend_built(METHOD_TYPE_GFN2))// &
                       " cuest="//available(cuest_backend_available())// &
                       " dlfind="//available(dlfind_available())// &
                       " ecp="//available(ecp_backend_available()))
   end subroutine print_version

   pure function available(is_available) result(text)
      !! "yes" or "no", in the shape a script can split on
      logical, intent(in) :: is_available
      character(len=:), allocatable :: text

      if (is_available) then
         text = "yes"
      else
         text = "no"
      end if
   end function available

end module mqc_version
