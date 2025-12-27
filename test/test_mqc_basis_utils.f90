module test_mqc_basis_utils
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_basis_utils, only: normalize_basis_name, find_basis_file
   implicit none
   private
   public :: collect_mqc_basis_utils_tests

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_basis_utils_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("normalize_basis_stars", test_normalize_stars), &
                  new_unittest("normalize_basis_plus", test_normalize_plus), &
                  new_unittest("normalize_basis_parens", test_normalize_parens), &
                  new_unittest("normalize_basis_complex", test_normalize_complex), &
                  new_unittest("normalize_basis_unchanged", test_normalize_unchanged), &
                  new_unittest("find_basis_file_exists", test_find_basis_exists), &
                  new_unittest("find_basis_file_not_found", test_find_basis_not_found) &
                  ]
   end subroutine collect_mqc_basis_utils_tests

   subroutine test_normalize_stars(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: result

      ! Single star
      result = normalize_basis_name("6-31G*")
      call check(error, trim(result), "6-31Gs", "6-31G* should become 6-31Gs")
      if (allocated(error)) return

      ! Double star
      result = normalize_basis_name("6-31G**")
      call check(error, trim(result), "6-31Gss", "6-31G** should become 6-31Gss")
   end subroutine test_normalize_stars

   subroutine test_normalize_plus(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: result

      ! Single plus
      result = normalize_basis_name("6-31+G")
      call check(error, trim(result), "6-31pG", "6-31+G should become 6-31pG")
      if (allocated(error)) return

      ! Double plus
      result = normalize_basis_name("6-311++G")
      call check(error, trim(result), "6-311ppG", "6-311++G should become 6-311ppG")
   end subroutine test_normalize_plus

   subroutine test_normalize_parens(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: result

      ! Single parentheses
      result = normalize_basis_name("6-31G(d)")
      call check(error, trim(result), "6-31Gd", "6-31G(d) should become 6-31Gd")
      if (allocated(error)) return

      ! Multiple in parentheses
      result = normalize_basis_name("6-311G(d,p)")
      call check(error, trim(result), "6-311Gdp", "6-311G(d,p) should become 6-311Gdp")
      if (allocated(error)) return

      ! Triple in parentheses
      result = normalize_basis_name("6-311G(2d,2p)")
      call check(error, trim(result), "6-311G2d2p", "6-311G(2d,2p) should become 6-311G2d2p")
   end subroutine test_normalize_parens

   subroutine test_normalize_complex(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: result

      ! Plus and star combined
      result = normalize_basis_name("6-31+G*")
      call check(error, trim(result), "6-31pGs", "6-31+G* should become 6-31pGs")
      if (allocated(error)) return

      ! Double plus and double star
      result = normalize_basis_name("6-311++G**")
      call check(error, trim(result), "6-311ppGss", "6-311++G** should become 6-311ppGss")
      if (allocated(error)) return

      ! Plus with parentheses
      result = normalize_basis_name("6-31+G(d)")
      call check(error, trim(result), "6-31pGd", "6-31+G(d) should become 6-31pGd")
   end subroutine test_normalize_complex

   subroutine test_normalize_unchanged(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: result

      ! Simple names should remain unchanged
      result = normalize_basis_name("6-31G")
      call check(error, trim(result), "6-31G", "6-31G should remain unchanged")
      if (allocated(error)) return

      result = normalize_basis_name("STO-3G")
      call check(error, trim(result), "STO-3G", "STO-3G should remain unchanged")
      if (allocated(error)) return

      result = normalize_basis_name("cc-pVDZ")
      call check(error, trim(result), "cc-pVDZ", "cc-pVDZ should remain unchanged")
   end subroutine test_normalize_unchanged

   subroutine test_find_basis_exists(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: filename
      integer :: stat
      character(len=:), allocatable :: errmsg
      logical :: file_exists

      ! Create basis_sets directory if it doesn't exist
      call execute_command_line("mkdir -p basis_sets", exitstat=stat)

      ! Create a test basis file in basis_sets directory
      call create_test_basis_file("basis_sets/test_basis.txt")

      ! Try to find it
      call find_basis_file("test_basis", filename, stat, errmsg)

      if (stat /= 0) then
         call check(error, .false., "Should find test_basis.txt: "//errmsg)
         call delete_file("basis_sets/test_basis.txt")
         return
      end if

      ! Check that the returned filename actually exists
      inquire (file=trim(filename), exist=file_exists)
      call check(error, file_exists, &
                 "Returned filename should exist: "//trim(filename))
      if (allocated(error)) then
         call delete_file("basis_sets/test_basis.txt")
         return
      end if

      ! Check that the filename is correct
      call check(error, trim(filename), "basis_sets/test_basis.txt", &
                 "Filename should be basis_sets/test_basis.txt")

      call delete_file("basis_sets/test_basis.txt")
   end subroutine test_find_basis_exists

   subroutine test_find_basis_not_found(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: filename
      integer :: stat
      character(len=:), allocatable :: errmsg

      ! Try to find non-existent basis
      call find_basis_file("NONEXISTENT-BASIS", filename, stat, errmsg)

      call check(error, stat /= 0, "Should fail to find non-existent basis")
   end subroutine test_find_basis_not_found

   ! Helper subroutines

   subroutine create_test_basis_file(filename)
      character(len=*), intent(in) :: filename
      integer :: unit
      open (newunit=unit, file=filename, status="replace")
      write (unit, '(a)') "$DATA"
      write (unit, '(a)') "$END"
      close (unit)
   end subroutine create_test_basis_file

   subroutine delete_file(filename)
      character(len=*), intent(in) :: filename
      integer :: unit, stat
      open (newunit=unit, file=filename, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete")
   end subroutine delete_file

end module test_mqc_basis_utils

program tester
   use iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_basis_utils, only: collect_mqc_basis_utils_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_basis_utils", collect_mqc_basis_utils_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program tester
