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
      call logger%info("metalquicha version "//MQC_VERSION_STR)
   end subroutine print_version

end module mqc_version
