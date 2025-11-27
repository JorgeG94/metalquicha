module test_mqc_basis_reader
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_basis_reader, only: classify_line, parse_element_basis, &
                               build_molecular_basis, ang_mom_int_to_char
   use mqc_cgto, only: atomic_basis_type, molecular_basis_type
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_basis_reader_tests

   character(len=*), parameter :: test_basis = &
                                  "$DATA"//new_line('a')// &
                                  ""//new_line('a')// &
                                  "HYDROGEN"//new_line('a')// &
                                  "S 2"//new_line('a')// &
                                  "1 1.0 2.0"//new_line('a')// &
                                  "2 1.0 2.0"//new_line('a')// &
                                  ""//new_line('a')// &
                                  "L 2"//new_line('a')// &
                                  "1 1.0 2.0 3.0"//new_line('a')// &
                                  "2 4.0 5.0 6.0"//new_line('a')// &
                                  ""//new_line('a')// &
                                  "$END"

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_basis_reader_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("classify_data_start", test_classify_data_start), &
                  new_unittest("classify_data_end", test_classify_data_end), &
                  new_unittest("classify_element", test_classify_element), &
                  new_unittest("classify_shell", test_classify_shell), &
                  new_unittest("classify_empty", test_classify_empty), &
                  new_unittest("parse_hydrogen_basis", test_parse_hydrogen_basis), &
                  new_unittest("build_h2_molecular_basis", test_build_h2_molecular_basis), &
                  new_unittest("ang_mom_conversions", test_ang_mom_conversions) &
                  ]
   end subroutine collect_mqc_basis_reader_tests

   subroutine test_classify_data_start(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: line_type

      line_type = classify_line("$DATA")
      call check(error, line_type, 0, "Should classify $DATA as LINE_UNKNOWN (0)")
   end subroutine test_classify_data_start

   subroutine test_classify_data_end(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: line_type

      line_type = classify_line("$END")
      call check(error, line_type, 0, "Should classify $END as LINE_UNKNOWN (0)")
   end subroutine test_classify_data_end

   subroutine test_classify_element(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: line_type

      line_type = classify_line("HYDROGEN")
      call check(error, line_type, 1, "Should classify element name as LINE_ATOM (1)")
      if (allocated(error)) return

      line_type = classify_line("CARBON")
      call check(error, line_type, 1, "Should classify CARBON as LINE_ATOM (1)")
   end subroutine test_classify_element

   subroutine test_classify_shell(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: line_type

      line_type = classify_line("S 2")
      call check(error, line_type, 2, "Should classify shell definition as LINE_SHELL (2)")
      if (allocated(error)) return

      line_type = classify_line("L 3")
      call check(error, line_type, 2, "Should classify L shell as LINE_SHELL (2)")
      if (allocated(error)) return

      line_type = classify_line("P 1")
      call check(error, line_type, 2, "Should classify P shell as LINE_SHELL (2)")
   end subroutine test_classify_shell

   subroutine test_classify_empty(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: line_type

      line_type = classify_line("")
      call check(error, line_type, 0, "Should classify empty line as type 0")
      if (allocated(error)) return

      line_type = classify_line("   ")
      call check(error, line_type, 0, "Should classify whitespace as type 0")
   end subroutine test_classify_empty

   subroutine test_parse_hydrogen_basis(error)
      type(error_type), allocatable, intent(out) :: error
      type(atomic_basis_type) :: h_basis
      integer :: stat
      character(len=:), allocatable :: errmsg

      call parse_element_basis(test_basis, "HYDROGEN", h_basis, stat, errmsg)

      call check(error, stat, 0, "Should successfully parse hydrogen basis")
      if (allocated(error)) return

      call check(error, trim(h_basis%element), "HYDROGEN", "Element name should match")
      if (allocated(error)) return

      ! L shells are split into S+P, so: 1 S shell + 1 L shell (split into 2) = 3 total
      call check(error, h_basis%nshells >= 2, "Hydrogen should have at least 2 shells in test basis")
      if (allocated(error)) return

      ! Check first shell (S shell)
      call check(error, h_basis%shells(1)%ang_mom, 0, "First shell should be S (ang_mom=0)")
      if (allocated(error)) return

      call check(error, h_basis%shells(1)%nfunc, 2, "S shell should have 2 primitives")
      if (allocated(error)) return

      call h_basis%destroy()
   end subroutine test_parse_hydrogen_basis

   subroutine test_build_h2_molecular_basis(error)
      type(error_type), allocatable, intent(out) :: error
      type(molecular_basis_type) :: h2_basis
      integer :: stat, nbf
      character(len=:), allocatable :: errmsg
      character(len=*), dimension(2), parameter :: h2_atoms = ["HYDROGEN", "HYDROGEN"]

      call build_molecular_basis(test_basis, h2_atoms, h2_basis, stat, errmsg)

      call check(error, stat, 0, "Should successfully build H2 molecular basis")
      if (allocated(error)) return

      call check(error, h2_basis%nelements, 2, "H2 should have 2 atoms")
      if (allocated(error)) return

      nbf = h2_basis%num_basis_functions()
      call check(error, nbf > 0, "Should have positive number of basis functions")
      if (allocated(error)) return

      ! Both atoms should have the same basis
      call check(error, h2_basis%elements(1)%nshells, h2_basis%elements(2)%nshells, &
                 "Both H atoms should have same number of shells")
      if (allocated(error)) return

      call h2_basis%destroy()
   end subroutine test_build_h2_molecular_basis

   subroutine test_ang_mom_conversions(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=1) :: ang_char

      ang_char = ang_mom_int_to_char(0)
      call check(error, ang_char, 'S', "ang_mom 0 should be S")
      if (allocated(error)) return

      ang_char = ang_mom_int_to_char(1)
      call check(error, ang_char, 'P', "ang_mom 1 should be P")
      if (allocated(error)) return

      ang_char = ang_mom_int_to_char(2)
      call check(error, ang_char, 'D', "ang_mom 2 should be D")
      if (allocated(error)) return

      ang_char = ang_mom_int_to_char(3)
      call check(error, ang_char, 'F', "ang_mom 3 should be F")
      if (allocated(error)) return

      ! -1 is used internally for L shells but ang_mom_int_to_char returns '?' for unknown
      ang_char = ang_mom_int_to_char(-1)
      call check(error, ang_char, '?', "ang_mom -1 should return '?' (unknown)")
      if (allocated(error)) return
   end subroutine test_ang_mom_conversions

end module test_mqc_basis_reader

program tester_mqc_basis_reader
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_basis_reader, only: collect_mqc_basis_reader_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_basis_reader", collect_mqc_basis_reader_tests) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_basis_reader
