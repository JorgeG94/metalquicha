module test_mqc_xyz_reader
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_geometry, only: geometry_type
   use mqc_xyz_reader, only: read_xyz_string, read_xyz_file
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_xyz_reader_tests

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_xyz_reader_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("read_water_molecule", test_read_water_molecule), &
                  new_unittest("read_single_atom", test_read_single_atom), &
                  new_unittest("empty_comment_line", test_empty_comment_line), &
                  new_unittest("error_insufficient_lines", test_error_insufficient_lines), &
                  new_unittest("error_invalid_atom_count", test_error_invalid_atom_count), &
                  new_unittest("error_malformed_coordinates", test_error_malformed_coordinates) &
                  ]
   end subroutine collect_mqc_xyz_reader_tests

   subroutine test_read_water_molecule(error)
      type(error_type), allocatable, intent(out) :: error
      type(geometry_type) :: geom
      integer :: stat
      character(len=:), allocatable :: errmsg
      character(len=*), parameter :: test_xyz = &
                                     "3"//new_line('a')// &
                                     "Water molecule"//new_line('a')// &
                                     "O    0.000000    0.000000    0.119262"//new_line('a')// &
                                     "H    0.000000    0.763239   -0.477047"//new_line('a')// &
                                     "H    0.000000   -0.763239   -0.477047"

      call read_xyz_string(test_xyz, geom, stat, errmsg)

      call check(error, stat, 0, "read_xyz_string should succeed")
      if (allocated(error)) return

      call check(error, geom%natoms, 3, "Water molecule should have 3 atoms")
      if (allocated(error)) return

      call check(error, trim(geom%comment), "Water molecule", "Comment should match")
      if (allocated(error)) return

      call check(error, trim(geom%elements(1)), "O", "First element should be O")
      if (allocated(error)) return

      call check(error, abs(geom%coords(3, 1) - 0.119262_dp) < 1.0e-6_dp, &
                 "Z-coordinate of oxygen should match")
      if (allocated(error)) return

      call geom%destroy()
   end subroutine test_read_water_molecule

   subroutine test_read_single_atom(error)
      type(error_type), allocatable, intent(out) :: error
      type(geometry_type) :: geom
      integer :: stat
      character(len=:), allocatable :: errmsg

      call read_xyz_string("1"//new_line('a')// &
                           "Single carbon atom"//new_line('a')// &
                           "C  1.0  2.0  3.0", geom, stat, errmsg)

      call check(error, stat, 0, "read_xyz_string should succeed for single atom")
      if (allocated(error)) return

      call check(error, geom%natoms, 1, "Should have 1 atom")
      if (allocated(error)) return

      call check(error, trim(geom%elements(1)), "C", "Element should be C")
      if (allocated(error)) return

      call check(error, abs(geom%coords(1, 1) - 1.0_dp) < 1.0e-6_dp, &
                 "X-coordinate should be 1.0")
      if (allocated(error)) return

      call geom%destroy()
   end subroutine test_read_single_atom

   subroutine test_empty_comment_line(error)
      type(error_type), allocatable, intent(out) :: error
      type(geometry_type) :: geom
      integer :: stat
      character(len=:), allocatable :: errmsg

      call read_xyz_string("2"//new_line('a')// &
                           new_line('a')// &
                           "He  0.0  0.0  0.0"//new_line('a')// &
                           "Ne  5.0  0.0  0.0", geom, stat, errmsg)

      call check(error, stat, 0, "read_xyz_string should handle empty comment")
      if (allocated(error)) return

      call check(error, geom%natoms, 2, "Should have 2 atoms")
      if (allocated(error)) return

      call geom%destroy()
   end subroutine test_empty_comment_line

   subroutine test_error_insufficient_lines(error)
      type(error_type), allocatable, intent(out) :: error
      type(geometry_type) :: geom
      integer :: stat
      character(len=:), allocatable :: errmsg

      call read_xyz_string("3"//new_line('a')// &
                           "Should fail"//new_line('a')// &
                           "H  0.0  0.0  0.0", geom, stat, errmsg)

      call check(error, stat /= 0, "Should detect insufficient lines")
      if (allocated(error)) return
   end subroutine test_error_insufficient_lines

   subroutine test_error_invalid_atom_count(error)
      type(error_type), allocatable, intent(out) :: error
      type(geometry_type) :: geom
      integer :: stat
      character(len=:), allocatable :: errmsg

      call read_xyz_string("not_a_number"//new_line('a')// &
                           "Comment", geom, stat, errmsg)

      call check(error, stat /= 0, "Should detect invalid atom count")
      if (allocated(error)) return
   end subroutine test_error_invalid_atom_count

   subroutine test_error_malformed_coordinates(error)
      type(error_type), allocatable, intent(out) :: error
      type(geometry_type) :: geom
      integer :: stat
      character(len=:), allocatable :: errmsg

      call read_xyz_string("1"//new_line('a')// &
                           "Test"//new_line('a')// &
                           "C  1.0  invalid  3.0", geom, stat, errmsg)

      call check(error, stat /= 0, "Should detect malformed coordinates")
      if (allocated(error)) return
   end subroutine test_error_malformed_coordinates

end module test_mqc_xyz_reader

program tester_mqc_xyz_reader
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_xyz_reader, only: collect_mqc_xyz_reader_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_xyz_reader", collect_mqc_xyz_reader_tests) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_xyz_reader
