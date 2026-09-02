module test_mqc_cli_parser
   !! Tests for the command line parser.
   !!
   !! `mqc_cli_parser` re-exports `normalize_basis_name` and `find_basis_file`
   !! from `mqc_basis_utils`. Those are exercised in `test_mqc_basis_utils`;
   !! duplicating them here only meant two copies to update whenever the
   !! naming rules move, which is exactly what went stale last time.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_cli_parser, only: cli_args_type
   implicit none
   private
   public :: collect_mqc_cli_parser_tests

contains

   !! Collect all exported unit tests
   subroutine collect_mqc_cli_parser_tests(testsuite)
      !! Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("cli_args_destroy", test_cli_args_destroy) &
                  ]
   end subroutine collect_mqc_cli_parser_tests

   subroutine test_cli_args_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(cli_args_type) :: args

      args%xyz_file = "test.xyz"
      args%basis_name = "6-31G"

      call args%destroy()

      call check(error,.not. allocated(args%xyz_file), &
                 "xyz_file should be deallocated")
      if (allocated(error)) return

      call check(error,.not. allocated(args%basis_name), &
                 "basis_name should be deallocated")
   end subroutine test_cli_args_destroy

end module test_mqc_cli_parser

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_cli_parser, only: collect_mqc_cli_parser_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_cli_parser", collect_mqc_cli_parser_tests) &
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
