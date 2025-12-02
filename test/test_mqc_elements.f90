module test_mqc_elements
   !! Unit tests for the periodic table elements module
   !!
   !! Tests element symbol/number conversions, atomic masses,
   !! and error handling for the complete periodic table.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_elements, only: element_symbol_to_number, element_number_to_symbol, element_mass
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_elements_tests  !! Test suite collection function

contains

   subroutine collect_mqc_elements_tests(testsuite)
      !! Collect all exported unit tests for elements module
      !!
      !! Assembles complete test suite covering symbol/number conversions,
      !! atomic masses, case sensitivity, and error conditions.
      type(unittest_type), allocatable, intent(out) :: testsuite(:)  !! Test collection array

      testsuite = [ &
                  new_unittest("symbol_to_number_hydrogen", test_symbol_to_number_h), &
                  new_unittest("symbol_to_number_carbon", test_symbol_to_number_c), &
                  new_unittest("symbol_to_number_case_insensitive", test_symbol_case_insensitive), &
                  new_unittest("symbol_to_number_invalid", test_symbol_invalid), &
                  new_unittest("number_to_symbol_hydrogen", test_number_to_symbol_h), &
                  new_unittest("number_to_symbol_carbon", test_number_to_symbol_c), &
                  new_unittest("number_to_symbol_oganesson", test_number_to_symbol_og), &
                  new_unittest("number_to_symbol_invalid", test_number_to_symbol_invalid), &
                  new_unittest("element_mass_hydrogen", test_mass_hydrogen), &
                  new_unittest("element_mass_carbon", test_mass_carbon), &
                  new_unittest("element_mass_gold", test_mass_gold), &
                  new_unittest("element_mass_invalid", test_mass_invalid), &
                  new_unittest("round_trip_conversion", test_round_trip) &
                  ]
   end subroutine collect_mqc_elements_tests

   subroutine test_symbol_to_number_h(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: z

      z = element_symbol_to_number("H")
      call check(error, z, 1, "Hydrogen should have atomic number 1")
   end subroutine test_symbol_to_number_h

   subroutine test_symbol_to_number_c(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: z

      z = element_symbol_to_number("C")
      call check(error, z, 6, "Carbon should have atomic number 6")
   end subroutine test_symbol_to_number_c

   subroutine test_symbol_case_insensitive(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: z1, z2, z3

      z1 = element_symbol_to_number("He")
      z2 = element_symbol_to_number("HE")
      z3 = element_symbol_to_number("he")

      call check(error, z1, 2, "He should give atomic number 2")
      if (allocated(error)) return

      call check(error, z2, 2, "HE should give atomic number 2")
      if (allocated(error)) return

      call check(error, z3, 2, "he should give atomic number 2")
      if (allocated(error)) return
   end subroutine test_symbol_case_insensitive

   subroutine test_symbol_invalid(error)
      type(error_type), allocatable, intent(out) :: error
      integer :: z

      z = element_symbol_to_number("Xx")
      call check(error, z, 0, "Invalid symbol should return 0")
   end subroutine test_symbol_invalid

   subroutine test_number_to_symbol_h(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=2) :: sym

      sym = element_number_to_symbol(1)
      call check(error, trim(sym), "H", "Atomic number 1 should be H")
   end subroutine test_number_to_symbol_h

   subroutine test_number_to_symbol_c(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=2) :: sym

      sym = element_number_to_symbol(6)
      call check(error, trim(sym), "C", "Atomic number 6 should be C")
   end subroutine test_number_to_symbol_c

   subroutine test_number_to_symbol_og(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=2) :: sym

      sym = element_number_to_symbol(118)
      call check(error, trim(sym), "Og", "Atomic number 118 should be Oganesson")
   end subroutine test_number_to_symbol_og

   subroutine test_number_to_symbol_invalid(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=2) :: sym

      sym = element_number_to_symbol(0)
      call check(error, trim(sym), "Xx", "Invalid atomic number should return Xx")
      if (allocated(error)) return

      sym = element_number_to_symbol(200)
      call check(error, trim(sym), "Xx", "Out of range atomic number should return Xx")
   end subroutine test_number_to_symbol_invalid

   subroutine test_mass_hydrogen(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: mass

      mass = element_mass(1)
      call check(error, abs(mass - 1.008_dp) < 1.0e-3_dp, &
                 "Hydrogen mass should be approximately 1.008 amu")
   end subroutine test_mass_hydrogen

   subroutine test_mass_carbon(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: mass

      mass = element_mass(6)
      call check(error, abs(mass - 12.011_dp) < 1.0e-3_dp, &
                 "Carbon mass should be approximately 12.011 amu")
   end subroutine test_mass_carbon

   subroutine test_mass_gold(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: mass

      mass = element_mass(79)
      call check(error, abs(mass - 196.97_dp) < 1.0e-2_dp, &
                 "Gold mass should be approximately 196.97 amu")
   end subroutine test_mass_gold

   subroutine test_mass_invalid(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: mass

      mass = element_mass(0)
      call check(error, abs(mass) < 1.0e-10_dp, &
                 "Invalid atomic number should return 0 mass")
      if (allocated(error)) return

      mass = element_mass(200)
      call check(error, abs(mass) < 1.0e-10_dp, &
                 "Out of range atomic number should return 0 mass")
   end subroutine test_mass_invalid

   subroutine test_round_trip(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=2) :: sym
      integer :: z

      ! Test round trip: symbol -> number -> symbol
      z = element_symbol_to_number("N")
      sym = element_number_to_symbol(z)
      call check(error, trim(sym), "N", "Round trip conversion should preserve symbol")
      if (allocated(error)) return

      ! Test another element
      z = element_symbol_to_number("Fe")
      sym = element_number_to_symbol(z)
      call check(error, trim(sym), "Fe", "Round trip for Fe should work")
   end subroutine test_round_trip

end module test_mqc_elements

program tester_mqc_elements
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_elements, only: collect_mqc_elements_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_elements", collect_mqc_elements_tests) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_elements
