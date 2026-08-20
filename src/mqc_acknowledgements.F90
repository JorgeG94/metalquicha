!! Credit to the upstream libraries a run actually used
module mqc_acknowledgements
   !! Prints a banner naming the integral backend a calculation is about to run
   !! on, and thanking the people who wrote it.
   !!
   !! Metalquicha computes almost no integrals of its own. Every number it
   !! reports comes out of somebody else's library -- cuEST on the GPU, libfint
   !! or libcint on the CPU, tblite for the semi-empirical methods -- and a code
   !! that
   !! prints its own logo at startup and never names them has its priorities
   !! backwards. This says which one ran, in the output the user keeps.
   !!
   !! Which backend is named is *derived*, never declared: it is resolved by
   !! `acknowledged_backend` below from the same three facts the dispatch in
   !! `mqc_method_hf`/`mqc_method_dft` uses -- the method, the requested
   !! backend, and whether this build has cuEST. Thanking a library that did not
   !! run would be worse than thanking none.
   !!
   !! `ACK_LIBCINT` names the CPU integrals slot rather than a library: which of
   !! libfint and libcint filled it is a build-time choice, so the banner text
   !! is the one thing here decided by the preprocessor rather than derived at
   !! run time. `MQC_WITH_LIBCINT` cannot answer it -- both branches define it.
   use pic_logger, only: logger => global_logger
   use mqc_method_types, only: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2
   use mqc_cuest_iface, only: parse_backend_name, method_runs_on_cuest, &
                              BACKEND_CUEST, BACKEND_LIBCINT
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
      !!
      !! This mirrors the dispatch in `mqc_method_hf` and `mqc_method_dft`, and
      !! has to keep mirroring it: an explicit request is honoured, and `auto`
      !! resolves to cuEST exactly when this build has cuEST *and* the method is
      !! one cuEST implements. Anything else is the CPU path.
      !!
      !! One known divergence, on a cuEST build only. `hf_run`'s `auto` branch
      !! is a bare `#ifdef MQC_WITH_CUEST` with no check that cuEST implements
      !! the method, so an `auto` MP2 or CCSD is handed to `run_cuest_scf` there
      !! -- while the explicit `cuest` branch right above it refuses exactly
      !! that combination. This function follows the explicit branch, which is
      !! the behaviour the refusal describes. When `hf_run` grows the same guard
      !! on `auto`, the two agree and this note can go.
      !!
      !! The semi-empirical methods never reach either -- they are tblite's own
      !! and take no `backend` -- so they are answered before the request is
      !! even looked at.
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
      case (BACKEND_LIBCINT)
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
