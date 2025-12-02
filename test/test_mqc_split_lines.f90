module test_mqc_split_lines_module
   !! Test module for split_lines subroutine
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_xyz_reader, only: split_lines
   implicit none
   private

   public :: collect_mqc_split_lines

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_split_lines(testsuite)
      !! Collect tests for split_lines functionality
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("test_simple_split", test_simple_split), &
                  new_unittest("test_empty_lines", test_empty_lines), &
                  new_unittest("test_debug_case", test_debug_case) &
                  ]
   end subroutine collect_mqc_split_lines

   subroutine test_simple_split(error)
      !! Test basic line splitting with LF endings
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: test_input = "1"//achar(10)//"comment"//achar(10)//"H 0 0 0"
      character(len=:), allocatable :: lines(:)
      integer :: nlines

      call split_lines(test_input, lines, nlines)

      call check(error, nlines, 3, "Expected 3 lines")
      if (allocated(error)) return

      call check(error, trim(lines(1)), "1", "First line should be '1'")
      if (allocated(error)) return

      call check(error, trim(lines(2)), "comment", "Second line should be 'comment'")
      if (allocated(error)) return

      call check(error, trim(lines(3)), "H 0 0 0", "Third line should be 'H 0 0 0'")
   end subroutine test_simple_split

   subroutine test_empty_lines(error)
      !! Test handling of empty lines
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: test_input = "1"//achar(10)//achar(10)//"H 0 0 0"
      character(len=:), allocatable :: lines(:)
      integer :: nlines

      call split_lines(test_input, lines, nlines)

      call check(error, nlines, 3, "Expected 3 lines with empty middle line")
      if (allocated(error)) return

      call check(error, trim(lines(1)), "1", "First line should be '1'")
      if (allocated(error)) return

      call check(error, trim(lines(2)), "", "Second line should be empty")
      if (allocated(error)) return

      call check(error, trim(lines(3)), "H 0 0 0", "Third line should be 'H 0 0 0'")
   end subroutine test_empty_lines

   subroutine test_debug_case(error)
      !! Test the specific case that was failing in the original tests
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: test_input = "1"//achar(10)//""//achar(10)//"H 0 0 0"
      character(len=:), allocatable :: lines(:)
      integer :: nlines
      integer :: i

      ! Add debug output
      write (*, '("DEBUG: Testing failing case input length: ", I0)') len(test_input)
      write (*, '("DEBUG: Input bytes:")')
      do i = 1, min(len(test_input), 10)
         write (*, '("  [", I0, "] = ", I0)') i, iachar(test_input(i:i))
      end do

      call split_lines(test_input, lines, nlines)

      write (*, '("DEBUG: Result nlines = ", I0)') nlines
      do i = 1, nlines
         write (*, '("DEBUG: Line ", I0, " = ''", A, "'' (len=", I0, ")")') i, lines(i), len_trim(lines(i))
      end do

      call check(error, nlines, 3, "Expected 3 lines")
      if (allocated(error)) return

      call check(error, trim(lines(1)), "1", "First line should be '1'")
      if (allocated(error)) return

      call check(error, trim(lines(2)), "", "Second line should be empty")
      if (allocated(error)) return

      call check(error, trim(lines(3)), "H 0 0 0", "Third line should be 'H 0 0 0'")
   end subroutine test_debug_case

end module test_mqc_split_lines_module

program tester
   !! Test program for split_lines subroutine using testdrive framework
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_split_lines_module
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_split_lines", collect_mqc_split_lines) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program tester
