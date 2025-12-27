module test_mqc_xyz_reader_detailed
   !! Detailed unit tests for XYZ reader components
   !!
   !! Tests individual scenarios and edge cases to isolate Intel compiler issues
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_geometry, only: geometry_type
   use mqc_xyz_reader, only: read_xyz_string, read_xyz_file
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_xyz_reader_detailed_tests

contains

   subroutine collect_mqc_xyz_reader_detailed_tests(testsuite)
      !! Collect all detailed unit tests for XYZ reader
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("test_minimal_xyz", test_minimal_xyz), &
                  new_unittest("test_whitespace_handling", test_whitespace_handling), &
                  new_unittest("test_debug_line_by_line", test_debug_line_by_line), &
                  new_unittest("test_direct_parsing", test_direct_parsing), &
                  new_unittest("test_string_construction", test_string_construction), &
                  new_unittest("test_newline_variants", test_newline_variants), &
                  new_unittest("test_intel_specific_cases", test_intel_specific_cases) &
                  ]
   end subroutine collect_mqc_xyz_reader_detailed_tests

   subroutine test_minimal_xyz(error)
      !! Test the most minimal possible XYZ file
      type(error_type), allocatable, intent(out) :: error

      type(geometry_type) :: geom
      integer :: stat
      character(len=:), allocatable :: errmsg
      character(len=*), parameter :: minimal_xyz = "1"//new_line('a')//""//new_line('a')//"H 0 0 0"

      write (*, *) "DEBUG: Testing minimal XYZ: '", minimal_xyz, "'"
      write (*, *) "DEBUG: Length of minimal_xyz: ", len(minimal_xyz)

      call read_xyz_string(minimal_xyz, geom, stat, errmsg)

      write (*, *) "DEBUG: stat = ", stat
      if (stat /= 0) then
         write (*, *) "DEBUG: errmsg = ", errmsg
      else
         write (*, *) "DEBUG: Success! geom%natoms = ", geom%natoms
         if (allocated(geom%elements)) then
            write (*, *) "DEBUG: elements allocated, size = ", size(geom%elements)
         else
            write (*, *) "DEBUG: elements not allocated"
         end if
         if (allocated(geom%coords)) then
            write (*, *) "DEBUG: coords allocated, size = ", size(geom%coords)
         else
            write (*, *) "DEBUG: coords not allocated"
         end if
      end if

      if (stat == 0) then
         call check(error, stat, 0, "Minimal XYZ should parse")
      else
         call check(error, stat, 0, "Minimal XYZ should parse: "//errmsg)
      end if
      if (allocated(error)) return

      if (stat == 0) then
         call check(error, geom%natoms, 1, "Should have 1 atom")
         if (allocated(error)) return

         write (*, *) "DEBUG: About to call geom%destroy()"
         call geom%destroy()
         write (*, *) "DEBUG: geom%destroy() completed successfully"
      end if

   end subroutine test_minimal_xyz

   subroutine test_whitespace_handling(error)
      !! Test various whitespace scenarios that might confuse Intel compiler
      type(error_type), allocatable, intent(out) :: error

      type(geometry_type) :: geom
      integer :: stat
      character(len=:), allocatable :: errmsg
      character(len=*), parameter :: spaced_xyz = " 1 "//new_line('a')// &
                                     "  Test  "//new_line('a')// &
                                     "  H   0.0   0.0   0.0  "

      call read_xyz_string(spaced_xyz, geom, stat, errmsg)

      if (stat == 0) then
         call check(error, stat, 0, "Whitespace XYZ should parse")
      else
         call check(error, stat, 0, "Whitespace XYZ should parse: "//errmsg)
      end if
      if (allocated(error)) return

      call check(error, geom%natoms, 1, "Should have 1 atom")
      if (allocated(error)) return

      call geom%destroy()
   end subroutine test_whitespace_handling

   subroutine test_debug_line_by_line(error)
      !! Test with detailed debugging of each step
      type(error_type), allocatable, intent(out) :: error

      type(geometry_type) :: geom
      integer :: stat
      character(len=:), allocatable :: errmsg
      character(len=*), parameter :: debug_xyz = "2"//new_line('a')// &
                                     "Debug test"//new_line('a')// &
                                     "H 0.0 0.0 0.0"//new_line('a')// &
                                     "He 1.0 0.0 0.0"

      write (*, *) "DEBUG: Input string length: ", len(debug_xyz)
      write (*, *) "DEBUG: Input string repr: '", debug_xyz, "'"

      ! Check that we can manually parse parts
      block
         integer :: test_natoms, test_stat
         character(len=100) :: first_line

         first_line = "2"
         read (first_line, *, iostat=test_stat) test_natoms
         write (*, *) "DEBUG: Manual parse of '2' gave stat=", test_stat, " natoms=", test_natoms
      end block

      call read_xyz_string(debug_xyz, geom, stat, errmsg)

      write (*, *) "DEBUG: Final result: stat=", stat
      if (stat /= 0) write (*, *) "DEBUG: Final error: ", errmsg

      if (stat == 0) then
         call check(error, stat, 0, "Debug XYZ should parse")
      else
         call check(error, stat, 0, "Debug XYZ should parse: "//errmsg)
      end if
      if (allocated(error)) return

      call check(error, geom%natoms, 2, "Should have 2 atoms")
      if (allocated(error)) return

      call geom%destroy()
   end subroutine test_debug_line_by_line

   subroutine test_direct_parsing(error)
      !! Test direct parsing operations that might fail with Intel
      type(error_type), allocatable, intent(out) :: error

      integer :: natoms, stat
      integer, parameter :: num_tests = 5
      character(len=100) :: test_strings(num_tests)
      integer :: i

      test_strings(1) = "1"
      test_strings(2) = "3"
      test_strings(3) = " 2 "
      test_strings(4) = "5   "
      test_strings(5) = "   10"

      do i = 1, 5
         read (test_strings(i), *, iostat=stat) natoms
         write (*, *) "DEBUG: Reading '", trim(test_strings(i)), "' gave stat=", stat, " value=", natoms

         call check(error, stat, 0, "Should read integer from: '"//trim(test_strings(i))//"'")
         if (allocated(error)) return
      end do

   end subroutine test_direct_parsing

   subroutine test_string_construction(error)
      !! Test how strings are constructed with new_line
      type(error_type), allocatable, intent(out) :: error

      character(len=:), allocatable :: test_str
      integer :: i

      test_str = "1"//new_line('a')//"test"//new_line('a')//"H 0 0 0"

      write (*, *) "DEBUG: Constructed string length: ", len(test_str)
      write (*, *) "DEBUG: String bytes:"
      do i = 1, len(test_str)
         write (*, '(A,I0,A,I0,A,A)') "  [", i, "] = ", ichar(test_str(i:i)), " '", test_str(i:i), "'"
      end do

      call check(error, len(test_str) > 5, "String should have reasonable length")

   end subroutine test_string_construction

   subroutine test_newline_variants(error)
      !! Test different newline character variants
      type(error_type), allocatable, intent(out) :: error

      type(geometry_type) :: geom
      integer :: stat
      character(len=:), allocatable :: errmsg
      character(len=*), parameter :: lf_xyz = "1"//achar(10)//"comment"//achar(10)//"H 0 0 0"
      character(len=*), parameter :: crlf_xyz = "1"//achar(13)//achar(10)//"comment"//achar(13)//achar(10)//"H 0 0 0"

      ! Test LF only
      call read_xyz_string(lf_xyz, geom, stat, errmsg)
      if (stat == 0) then
         call check(error, stat, 0, "LF-only XYZ should parse")
      else
         call check(error, stat, 0, "LF-only XYZ should parse: "//errmsg)
      end if
      if (allocated(error)) return
      call geom%destroy()

      ! Test CRLF
      call read_xyz_string(crlf_xyz, geom, stat, errmsg)
      if (stat == 0) then
         call check(error, stat, 0, "CRLF XYZ should parse")
      else
         call check(error, stat, 0, "CRLF XYZ should parse: "//errmsg)
      end if
      if (allocated(error)) return
      call geom%destroy()

   end subroutine test_newline_variants

   subroutine test_intel_specific_cases(error)
      !! Test cases that specifically might cause Intel compiler issues
      type(error_type), allocatable, intent(out) :: error

      type(geometry_type) :: geom
      integer :: stat
      character(len=:), allocatable :: errmsg

      ! Test the exact case from the failing test
      character(len=*), parameter :: water_xyz = &
                                     "3"//new_line('a')// &
                                     "Water molecule"//new_line('a')// &
                                     "O    0.000000    0.000000    0.119262"//new_line('a')// &
                                     "H    0.000000    0.763239   -0.477047"//new_line('a')// &
                                     "H    0.000000   -0.763239   -0.477047"

      write (*, *) "DEBUG: Testing original failing case"
      write (*, *) "DEBUG: String length: ", len(water_xyz)

      call read_xyz_string(water_xyz, geom, stat, errmsg)

      write (*, *) "DEBUG: Result stat=", stat
      if (stat /= 0) write (*, *) "DEBUG: Error message: ", errmsg

      if (stat == 0) then
         call check(error, stat, 0, "Original water test should work")
      else
         call check(error, stat, 0, "Original water test should work: "//errmsg)
      end if
      if (allocated(error)) return

      call check(error, geom%natoms, 3, "Water should have 3 atoms")
      if (allocated(error)) return

      call geom%destroy()
   end subroutine test_intel_specific_cases

end module test_mqc_xyz_reader_detailed

program tester_mqc_xyz_reader_detailed
   !! Test runner for detailed XYZ reader tests
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_xyz_reader_detailed, only: collect_mqc_xyz_reader_detailed_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_xyz_reader_detailed", collect_mqc_xyz_reader_detailed_tests) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_xyz_reader_detailed
