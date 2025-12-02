module test_mqc_frag_utils
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_frag_utils
   use pic_types, only: default_int, dp
   implicit none
   private
   public :: collect_mqc_frag_utils_tests

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_frag_utils_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("binomial_base_cases", test_binomial_base_cases), &
                  new_unittest("binomial_small_values", test_binomial_small_values), &
                  new_unittest("binomial_larger_values", test_binomial_larger_values), &
                  new_unittest("binomial_edge_cases", test_binomial_edge_cases), &
                  new_unittest("get_nfrags_single_level", test_get_nfrags_single), &
                  new_unittest("get_nfrags_multiple_levels", test_get_nfrags_multiple), &
                  new_unittest("get_nfrags_all_levels", test_get_nfrags_all), &
                  new_unittest("create_monomer_list_small", test_create_monomer_list_small), &
                  new_unittest("create_monomer_list_larger", test_create_monomer_list_larger), &
                  new_unittest("generate_fragment_list_dimers", test_generate_dimers), &
                  new_unittest("generate_fragment_list_trimers", test_generate_trimers), &
                  new_unittest("generate_fragment_list_mixed", test_generate_mixed) &
                  ]
   end subroutine collect_mqc_frag_utils_tests

   subroutine test_binomial_base_cases(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int) :: result

      ! C(n, 0) = 1 for any n
      result = binomial(5_default_int, 0_default_int)
      call check(error, result, 1_default_int, "C(5,0) should be 1")
      if (allocated(error)) return

      ! C(n, n) = 1 for any n
      result = binomial(5_default_int, 5_default_int)
      call check(error, result, 1_default_int, "C(5,5) should be 1")
      if (allocated(error)) return

      ! C(n, r) = 0 when r > n
      result = binomial(3_default_int, 5_default_int)
      call check(error, result, 0_default_int, "C(3,5) should be 0")
   end subroutine test_binomial_base_cases

   subroutine test_binomial_small_values(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int) :: result

      ! C(4, 2) = 6
      result = binomial(4_default_int, 2_default_int)
      call check(error, result, 6_default_int, "C(4,2) should be 6")
      if (allocated(error)) return

      ! C(5, 2) = 10
      result = binomial(5_default_int, 2_default_int)
      call check(error, result, 10_default_int, "C(5,2) should be 10")
      if (allocated(error)) return

      ! C(6, 3) = 20
      result = binomial(6_default_int, 3_default_int)
      call check(error, result, 20_default_int, "C(6,3) should be 20")
   end subroutine test_binomial_small_values

   subroutine test_binomial_larger_values(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int) :: result

      ! C(10, 3) = 120
      result = binomial(10_default_int, 3_default_int)
      call check(error, result, 120_default_int, "C(10,3) should be 120")
      if (allocated(error)) return

      ! C(10, 5) = 252
      result = binomial(10_default_int, 5_default_int)
      call check(error, result, 252_default_int, "C(10,5) should be 252")
      if (allocated(error)) return

      ! C(7, 4) = 35
      result = binomial(7_default_int, 4_default_int)
      call check(error, result, 35_default_int, "C(7,4) should be 35")
   end subroutine test_binomial_larger_values

   subroutine test_binomial_edge_cases(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int) :: result

      ! C(1, 1) = 1
      result = binomial(1_default_int, 1_default_int)
      call check(error, result, 1_default_int, "C(1,1) should be 1")
      if (allocated(error)) return

      ! C(10, 1) = 10
      result = binomial(10_default_int, 1_default_int)
      call check(error, result, 10_default_int, "C(10,1) should be 10")
   end subroutine test_binomial_edge_cases

   subroutine test_get_nfrags_single(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int) :: result

      ! 5 monomers, up to dimers (level 2)
      ! Should have C(5,1) + C(5,2) = 5 + 10 = 15
      result = get_nfrags(5_default_int, 2_default_int)
      call check(error, result, 15_default_int, "5 monomers up to dimers should be 15")
      if (allocated(error)) return

      ! 4 monomers, up to monomers only (level 1)
      ! Should have C(4,1) = 4
      result = get_nfrags(4_default_int, 1_default_int)
      call check(error, result, 4_default_int, "4 monomers up to level 1 should be 4")
   end subroutine test_get_nfrags_single

   subroutine test_get_nfrags_multiple(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int) :: result

      ! 5 monomers, up to trimers (level 3)
      ! C(5,1) + C(5,2) + C(5,3) = 5 + 10 + 10 = 25
      result = get_nfrags(5_default_int, 3_default_int)
      call check(error, result, 25_default_int, "5 monomers up to trimers should be 25")
      if (allocated(error)) return

      ! 6 monomers, up to dimers (level 2)
      ! C(6,1) + C(6,2) = 6 + 15 = 21
      result = get_nfrags(6_default_int, 2_default_int)
      call check(error, result, 21_default_int, "6 monomers up to dimers should be 21")
   end subroutine test_get_nfrags_multiple

   subroutine test_get_nfrags_all(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int) :: result

      ! 4 monomers, all levels (level 4)
      ! C(4,1) + C(4,2) + C(4,3) + C(4,4) = 4 + 6 + 4 + 1 = 15
      result = get_nfrags(4_default_int, 4_default_int)
      call check(error, result, 15_default_int, "4 monomers all levels should be 15 (2^4-1)")
      if (allocated(error)) return

      ! 3 monomers, all levels (level 3)
      ! C(3,1) + C(3,2) + C(3,3) = 3 + 3 + 1 = 7
      result = get_nfrags(3_default_int, 3_default_int)
      call check(error, result, 7_default_int, "3 monomers all levels should be 7 (2^3-1)")
   end subroutine test_get_nfrags_all

   subroutine test_create_monomer_list_small(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), allocatable :: monomers(:)
      integer(default_int) :: i

      ! Create a list of 3 monomers
      allocate (monomers(3))
      call create_monomer_list(monomers)

      call check(error, monomers(1), 1_default_int, "First monomer should be 1")
      if (allocated(error)) then
         deallocate (monomers)
         return
      end if

      call check(error, monomers(2), 2_default_int, "Second monomer should be 2")
      if (allocated(error)) then
         deallocate (monomers)
         return
      end if

      call check(error, monomers(3), 3_default_int, "Third monomer should be 3")

      deallocate (monomers)
   end subroutine test_create_monomer_list_small

   subroutine test_create_monomer_list_larger(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), allocatable :: monomers(:)
      integer(default_int) :: i
      logical :: all_correct

      ! Create a list of 10 monomers
      allocate (monomers(10))
      call create_monomer_list(monomers)

      all_correct = .true.
      do i = 1, 10
         if (monomers(i) /= i) then
            all_correct = .false.
            exit
         end if
      end do

      call check(error, all_correct, "All monomers should be numbered correctly from 1 to 10")

      deallocate (monomers)
   end subroutine test_create_monomer_list_larger

   subroutine test_generate_dimers(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), allocatable :: monomers(:)
      integer(default_int), allocatable :: polymers(:, :)
      integer(default_int) :: count, expected_count
      integer(default_int) :: max_level, n_monomers

      n_monomers = 4_default_int
      max_level = 2_default_int

      ! Create monomer list
      allocate (monomers(n_monomers))
      call create_monomer_list(monomers)

      ! Calculate expected fragments (only dimers, not monomers)
      expected_count = binomial(n_monomers, 2_default_int)

      ! Allocate output array
      allocate (polymers(expected_count, max_level))
      polymers = 0
      count = 0

      ! Generate fragments
      call generate_fragment_list(monomers, max_level, polymers, count)

      ! Check that we got the right number of dimers
      call check(error, count, expected_count, "Should generate 6 dimers from 4 monomers")
      if (allocated(error)) then
         deallocate (monomers, polymers)
         return
      end if

      ! Check first dimer is [1, 2]
      call check(error, polymers(1, 1), 1_default_int, "First dimer first element should be 1")
      if (allocated(error)) then
         deallocate (monomers, polymers)
         return
      end if

      call check(error, polymers(1, 2), 2_default_int, "First dimer second element should be 2")
      if (allocated(error)) then
         deallocate (monomers, polymers)
         return
      end if

      ! Check last dimer is [3, 4]
      call check(error, polymers(count, 1), 3_default_int, "Last dimer first element should be 3")
      if (allocated(error)) then
         deallocate (monomers, polymers)
         return
      end if

      call check(error, polymers(count, 2), 4_default_int, "Last dimer second element should be 4")

      deallocate (monomers, polymers)
   end subroutine test_generate_dimers

   subroutine test_generate_trimers(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), allocatable :: monomers(:)
      integer(default_int), allocatable :: polymers(:, :)
      integer(default_int) :: count, expected_count
      integer(default_int) :: max_level, n_monomers

      n_monomers = 5_default_int
      max_level = 3_default_int

      ! Create monomer list
      allocate (monomers(n_monomers))
      call create_monomer_list(monomers)

      ! Calculate expected fragments (dimers + trimers)
      expected_count = binomial(n_monomers, 2_default_int) + binomial(n_monomers, 3_default_int)

      ! Allocate output array
      allocate (polymers(expected_count, max_level))
      polymers = 0
      count = 0

      ! Generate fragments
      call generate_fragment_list(monomers, max_level, polymers, count)

      ! Check that we got the right number of fragments
      ! C(5,2) + C(5,3) = 10 + 10 = 20
      call check(error, count, expected_count, "Should generate 20 fragments (10 dimers + 10 trimers)")
      if (allocated(error)) then
         deallocate (monomers, polymers)
         return
      end if

      ! The first 10 should be dimers (with third column = 0)
      call check(error, polymers(1, 3), 0_default_int, "First fragment should be a dimer (3rd col = 0)")
      if (allocated(error)) then
         deallocate (monomers, polymers)
         return
      end if

      ! Fragment 11 should be the first trimer [1, 2, 3]
      call check(error, polymers(11, 1), 1_default_int, "First trimer element 1")
      if (allocated(error)) then
         deallocate (monomers, polymers)
         return
      end if

      call check(error, polymers(11, 2), 2_default_int, "First trimer element 2")
      if (allocated(error)) then
         deallocate (monomers, polymers)
         return
      end if

      call check(error, polymers(11, 3), 3_default_int, "First trimer element 3")

      deallocate (monomers, polymers)
   end subroutine test_generate_trimers

   subroutine test_generate_mixed(error)
      type(error_type), allocatable, intent(out) :: error
      integer(default_int), allocatable :: monomers(:)
      integer(default_int), allocatable :: polymers(:, :)
      integer(default_int) :: count, expected_count
      integer(default_int) :: max_level, n_monomers

      n_monomers = 3_default_int
      max_level = 3_default_int

      ! Create monomer list
      allocate (monomers(n_monomers))
      call create_monomer_list(monomers)

      ! Calculate expected fragments (dimers + trimers)
      ! C(3,2) + C(3,3) = 3 + 1 = 4
      expected_count = binomial(n_monomers, 2_default_int) + binomial(n_monomers, 3_default_int)

      ! Allocate output array
      allocate (polymers(expected_count, max_level))
      polymers = 0
      count = 0

      ! Generate fragments
      call generate_fragment_list(monomers, max_level, polymers, count)

      call check(error, count, 4_default_int, "Should generate 4 fragments from 3 monomers")
      if (allocated(error)) then
         deallocate (monomers, polymers)
         return
      end if

      ! Last fragment should be [1, 2, 3]
      call check(error, polymers(4, 1), 1_default_int, "Full trimer element 1")
      if (allocated(error)) then
         deallocate (monomers, polymers)
         return
      end if

      call check(error, polymers(4, 2), 2_default_int, "Full trimer element 2")
      if (allocated(error)) then
         deallocate (monomers, polymers)
         return
      end if

      call check(error, polymers(4, 3), 3_default_int, "Full trimer element 3")

      deallocate (monomers, polymers)
   end subroutine test_generate_mixed

end module test_mqc_frag_utils

program tester
   use iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_frag_utils, only: collect_mqc_frag_utils_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_frag_utils", collect_mqc_frag_utils_tests) &
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
