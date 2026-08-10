module test_mqc_capi_session
   !! Tests the C entry points for the MPI session that a caller can reach
   !! without one running.
   !!
   !! `begin`/`end` drive the MPI lifecycle and never return on a worker, so
   !! they belong to the MPI validation harness rather than a serial unit test.
   !! What is left is the part a script queries before -- or instead of --
   !! starting a session: the state accessors, which must answer with defined
   !! sentinels rather than reading uninitialised MPI state. Without a session
   !! that means inactive, rank -1, size 0, not root, and no error yet.
   use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_capi_session, only: mqc_session_rank, mqc_session_size, &
                               mqc_session_is_root, mqc_session_active, &
                               mqc_session_last_error
   implicit none
   private
   public :: collect_mqc_capi_session_tests

contains

   subroutine collect_mqc_capi_session_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("no_session_is_inactive", test_inactive), &
                  new_unittest("no_session_is_not_root", test_not_root), &
                  new_unittest("no_session_rank_and_size", test_rank_size), &
                  new_unittest("no_error_yet", test_no_error) &
                  ]
   end subroutine collect_mqc_capi_session_tests

   subroutine test_inactive(error)
      !! MPI is not up until begin is called
      type(error_type), allocatable, intent(out) :: error

      call check(error, int(mqc_session_active()), 0, "no session is inactive")
   end subroutine test_inactive

   subroutine test_not_root(error)
      !! is_root reads rank 0; the sentinel rank is not it
      type(error_type), allocatable, intent(out) :: error

      call check(error, int(mqc_session_is_root()), 0, "no session is not root")
   end subroutine test_not_root

   subroutine test_rank_size(error)
      !! The sentinels are -1 and 0, not whatever was on the stack
      type(error_type), allocatable, intent(out) :: error

      call check(error, int(mqc_session_rank()), -1, "no session -> rank -1")
      if (allocated(error)) return
      call check(error, int(mqc_session_size()), 0, "no session -> size 0")
   end subroutine test_rank_size

   subroutine test_no_error(error)
      !! Nothing has failed, so the message is empty
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: BUFLEN = 512
      character(kind=c_char) :: buffer(BUFLEN)

      ! A short buffer is a no-op that must not overrun; a zero one likewise.
      buffer = 'x'
      call mqc_session_last_error(0_c_int, buffer)

      buffer = 'x'
      call mqc_session_last_error(int(BUFLEN, c_int), buffer)
      call check(error, buffer(1) == c_null_char, &
                 "an untouched session has no error message")
   end subroutine test_no_error

end module test_mqc_capi_session

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_capi_session, only: collect_mqc_capi_session_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_capi_session", collect_mqc_capi_session_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if

end program tester
