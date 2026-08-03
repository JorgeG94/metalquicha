module test_mqc_gbs_reader
   !! Unit tests for the Gaussian94 (.gbs) basis reader and shell normalization
   !!
   !! These run on any machine -- they touch no GPU and no cuEST -- so the
   !! basis half of the cuEST integration can be verified without queue time.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_gbs_reader, only: read_gbs_element, build_molecular_basis_gbs
   use mqc_basis_normalization, only: normalized_coefficients
   use mqc_cgto, only: atomic_basis_type, molecular_basis_type
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_gbs_reader_tests

   character(len=*), parameter :: TEST_FILE = "test_gbs_scratch.gbs"

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_gbs_reader_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("parse_hydrogen", test_parse_hydrogen), &
                  new_unittest("parse_second_element", test_parse_second_element), &
                  new_unittest("missing_element_errors", test_missing_element_errors), &
                  new_unittest("build_molecular_basis", test_build_molecular_basis), &
                  new_unittest("normalization_s_shell", test_normalization_s_shell), &
                  new_unittest("normalization_p_shell", test_normalization_p_shell) &
                  ]
   end subroutine collect_mqc_gbs_reader_tests

   subroutine write_test_file()
      !! Write a two-element Gaussian94 file mirroring the def2-SVP layout
      integer :: unit

      open (newunit=unit, file=TEST_FILE, status='replace', action='write')
      write (unit, '(A)') "spherical"
      write (unit, '(A)') ""
      write (unit, '(A)') "! a comment line that must be ignored"
      write (unit, '(A)') "****"
      write (unit, '(A)') "H     0"
      write (unit, '(A)') "S    3   1.00"
      write (unit, '(A)') "     13.0107010              0.19682158D-01"
      write (unit, '(A)') "      1.9622572              0.13796524"
      write (unit, '(A)') "      0.44453796             0.47831935"
      write (unit, '(A)') "S    1   1.00"
      write (unit, '(A)') "      0.12194962             1.0000000"
      write (unit, '(A)') "P    1   1.00"
      write (unit, '(A)') "      0.8000000              1.0000000"
      write (unit, '(A)') "****"
      write (unit, '(A)') "O     0"
      write (unit, '(A)') "S    2   1.00"
      write (unit, '(A)') "   2266.1767785             -0.53431809926D-02"
      write (unit, '(A)') "    340.87010191            -0.39890039230D-01"
      write (unit, '(A)') "P    1   1.00"
      write (unit, '(A)') "      0.2700058226            1.0000000"
      write (unit, '(A)') "****"
      close (unit)
   end subroutine write_test_file

   subroutine remove_test_file()
      !! Delete the scratch basis file
      integer :: unit
      logical :: exists

      inquire (file=TEST_FILE, exist=exists)
      if (.not. exists) return
      open (newunit=unit, file=TEST_FILE, status='old')
      close (unit, status='delete')
   end subroutine remove_test_file

   subroutine test_parse_hydrogen(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(atomic_basis_type) :: atom_basis
      type(error_t) :: parse_error

      call write_test_file()
      call read_gbs_element(TEST_FILE, "H", atom_basis, parse_error)
      call remove_test_file()

      call check(error,.not. parse_error%has_error(), "parsing H should succeed")
      if (allocated(error)) return

      call check(error, atom_basis%nshells, 3, "H should have 3 shells (2S + 1P)")
      if (allocated(error)) return

      call check(error, atom_basis%shells(1)%ang_mom, 0, "shell 1 should be S")
      if (allocated(error)) return
      call check(error, atom_basis%shells(3)%ang_mom, 1, "shell 3 should be P")
      if (allocated(error)) return

      call check(error, atom_basis%shells(1)%nfunc, 3, "first S shell has 3 primitives")
      if (allocated(error)) return

      ! The D exponent form must be read natively by list-directed input
      call check(error, atom_basis%shells(1)%coefficients(1), 0.019682158_dp, &
                 "D-exponent coefficient parsed", thr=1.0e-12_dp)
      if (allocated(error)) return

      call check(error, atom_basis%shells(1)%exponents(1), 13.0107010_dp, &
                 "first exponent parsed", thr=1.0e-10_dp)
   end subroutine test_parse_hydrogen

   subroutine test_parse_second_element(error)
      !! An element after the first must be found, not shadowed by it
      type(error_type), allocatable, intent(out) :: error
      type(atomic_basis_type) :: atom_basis
      type(error_t) :: parse_error

      call write_test_file()
      call read_gbs_element(TEST_FILE, "O", atom_basis, parse_error)
      call remove_test_file()

      call check(error,.not. parse_error%has_error(), "parsing O should succeed")
      if (allocated(error)) return

      call check(error, atom_basis%nshells, 2, "O should have 2 shells")
      if (allocated(error)) return

      call check(error, atom_basis%shells(1)%nfunc, 2, "O first shell has 2 primitives")
      if (allocated(error)) return

      call check(error, atom_basis%shells(1)%coefficients(2), -0.039890039230_dp, &
                 "negative D-exponent coefficient parsed", thr=1.0e-12_dp)
   end subroutine test_parse_second_element

   subroutine test_missing_element_errors(error)
      !! An element not in the file must be an error, not a silent empty basis
      type(error_type), allocatable, intent(out) :: error
      type(atomic_basis_type) :: atom_basis
      type(error_t) :: parse_error

      call write_test_file()
      call read_gbs_element(TEST_FILE, "Fe", atom_basis, parse_error)
      call remove_test_file()

      call check(error, parse_error%has_error(), "missing element must report an error")
   end subroutine test_missing_element_errors

   subroutine test_build_molecular_basis(error)
      !! Every atom gets its element's shells, in geometry order
      type(error_type), allocatable, intent(out) :: error
      type(molecular_basis_type) :: mol_basis
      type(error_t) :: parse_error
      character(len=2) :: symbols(3)

      symbols = ["O ", "H ", "H "]

      call write_test_file()
      call build_molecular_basis_gbs(TEST_FILE, symbols, mol_basis, parse_error)
      call remove_test_file()

      call check(error,.not. parse_error%has_error(), "building water basis should succeed")
      if (allocated(error)) return

      call check(error, mol_basis%nelements, 3, "one basis entry per atom")
      if (allocated(error)) return

      call check(error, mol_basis%elements(1)%nshells, 2, "atom 1 is O with 2 shells")
      if (allocated(error)) return
      call check(error, mol_basis%elements(2)%nshells, 3, "atom 2 is H with 3 shells")
      if (allocated(error)) return
      call check(error, mol_basis%elements(3)%nshells, 3, "atom 3 is H with 3 shells")
      if (allocated(error)) return

      ! num_basis_functions is Cartesian: O has S(1) + P(3) = 4
      call check(error, mol_basis%elements(1)%num_basis_functions(), 4, &
                 "O contributes 4 Cartesian functions")
   end subroutine test_build_molecular_basis

   subroutine test_normalization_s_shell(error)
      !! A single uncontracted s primitive normalizes to the known closed form
      !!
      !! For one primitive with coefficient 1, the normalized coefficient is
      !! (2a/pi)^(3/4), the standard 1s Gaussian normalization constant.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: norm_error
      real(dp) :: exponents(1), coefficients(1), normalized(1), expected
      real(dp), parameter :: PI = 3.14159265358979323846_dp

      exponents = [1.25_dp]
      coefficients = [1.0_dp]

      call normalized_coefficients(0, exponents, coefficients, 1.0_dp, normalized, norm_error)

      call check(error,.not. norm_error%has_error(), "s-shell normalization should succeed")
      if (allocated(error)) return

      expected = (2.0_dp*exponents(1)/PI)**0.75_dp
      call check(error, normalized(1), expected, &
                 "single s primitive matches (2a/pi)^(3/4)", thr=1.0e-12_dp)
   end subroutine test_normalization_s_shell

   subroutine test_normalization_p_shell(error)
      !! A single uncontracted p primitive normalizes to its closed form
      !!
      !! Expected: 2 * a^(1/2) * (2a/pi)^(3/4), the standard 2p constant.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: norm_error
      real(dp) :: exponents(1), coefficients(1), normalized(1), expected
      real(dp), parameter :: PI = 3.14159265358979323846_dp

      exponents = [0.8_dp]
      coefficients = [1.0_dp]

      call normalized_coefficients(1, exponents, coefficients, 1.0_dp, normalized, norm_error)

      call check(error,.not. norm_error%has_error(), "p-shell normalization should succeed")
      if (allocated(error)) return

      expected = 2.0_dp*sqrt(exponents(1))*(2.0_dp*exponents(1)/PI)**0.75_dp
      call check(error, normalized(1), expected, &
                 "single p primitive matches 2 sqrt(a) (2a/pi)^(3/4)", thr=1.0e-12_dp)
   end subroutine test_normalization_p_shell

end module test_mqc_gbs_reader

program tester_mqc_gbs_reader
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_gbs_reader, only: collect_mqc_gbs_reader_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_gbs_reader", collect_mqc_gbs_reader_tests) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_gbs_reader
