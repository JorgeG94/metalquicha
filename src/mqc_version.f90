!! Version information for metalquicha
module mqc_version
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: MQC_VERSION_STR, print_version

   character(len=*), parameter :: MQC_VERSION_STR = "0.1.0"

contains

   subroutine print_version()
      call logger%info("metalquicha version "//MQC_VERSION_STR)
   end subroutine print_version

end module mqc_version
