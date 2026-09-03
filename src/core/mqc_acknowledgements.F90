!! Credit to the upstream libraries a run actually used
module mqc_acknowledgements
   !! Prints a banner naming the integral backend a calculation is about to run
   !! on, and thanking the people who wrote it.
   use pic_logger, only: logger => global_logger
   use mqc_method_types, only: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2
   use mqc_cuest_iface, only: parse_backend_name, method_runs_on_cuest, &
                              BACKEND_CUEST, BACKEND_CZT
   use mqc_cuest_bridge, only: cuest_backend_available
   use mqc_error, only: error_t
   implicit none
   private

   public :: print_acknowledgement   !! Banner for the backend this run will use
   public :: acknowledged_backend    !! Which backend that is; exposed for testing
   public :: ACK_NONE, ACK_CUEST, ACK_LIBCINT, ACK_TBLITE

   integer, parameter :: ACK_NONE = 0     !! Nothing to credit (unknown method)
   integer, parameter :: ACK_CUEST = 1
   integer, parameter :: ACK_LIBCINT = 2
   integer, parameter :: ACK_TBLITE = 3

contains

   subroutine print_acknowledgement(method_type, backend)
      !! Name the backend this run will use, and thank the people behind it
      !!
      !! Call once, from one rank. At `info` level, so `system.logger.level`
      !! can quiet it along with everything else.
      integer, intent(in) :: method_type  !! From `parse_method_string`
      character(len=*), intent(in) :: backend  !! The deck's `backend` request

      select case (acknowledged_backend(method_type, backend))
      case (ACK_CUEST)
         call banner("Integrals and SCF provided through NVIDIA's cuEST")
      case (ACK_LIBCINT)
#ifdef MQC_WITH_LIBFINT
         ! libfint is a port rather than a new derivation, so the credit is
         ! still Qiming Sun's; naming only the port would take it away from him,
         ! and naming only libcint would credit a library this build did not
         ! link. Both, in the order they matter.
         call banner("Integrals provided through libfint, a Fortran port of "// &
                     "Qiming Sun's LibCint")
#else
         call banner("Integrals provided through Qiming Sun's LibCint")
#endif
      case (ACK_TBLITE)
         call banner("Semi-empirical methods provided through the Grimme group's tblite")
      case default
         ! Nothing ran that has a backend to credit. Silence beats a guess.
      end select
   end subroutine print_acknowledgement

   function acknowledged_backend(method_type, backend) result(which)
      !! Which library will actually do the work, as one of ACK_*
      integer, intent(in) :: method_type
      character(len=*), intent(in) :: backend
      integer :: which

      integer :: kind
      type(error_t) :: parse_error

      if (method_type == METHOD_TYPE_GFN1 .or. method_type == METHOD_TYPE_GFN2) then
         which = ACK_TBLITE
         return
      end if

      call parse_backend_name(backend, kind, parse_error)
      if (parse_error%has_error()) then
         ! Unreachable through a deck -- the reader refuses an unknown backend
         ! long before this -- so say nothing rather than credit a guess.
         which = ACK_NONE
         return
      end if

      select case (kind)
      case (BACKEND_CUEST)
         which = ACK_CUEST
      case (BACKEND_CZT)
         which = ACK_LIBCINT
      case default
         if (cuest_backend_available() .and. method_runs_on_cuest(method_type)) then
            which = ACK_CUEST
         else
            which = ACK_LIBCINT
         end if
      end select
   end function acknowledged_backend

   subroutine banner(text)
      !! One line in a box sized to it
      !!
      !! Sized from the text rather than to a fixed width, so a reworded line
      !! cannot leave the box crooked and nobody has to count characters.
      character(len=*), intent(in) :: text

      character(len=:), allocatable :: rule

      rule = repeat("═", len(text) + 2)
      call logger%info(" ")
      call logger%info("    ╔"//rule//"╗")
      call logger%info("    ║ "//text//" ║")
      call logger%info("    ╚"//rule//"╝")
      call logger%info(" ")
   end subroutine banner

end module mqc_acknowledgements
