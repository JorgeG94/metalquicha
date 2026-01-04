module test_mqc_mbe
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_mbe, only: compute_mbe_energy
   use mqc_result_types, only: calculation_result_t
   use pic_types, only: dp, int64
   implicit none
   private
   public :: collect_mqc_mbe_tests

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_mbe_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("mbe_monomers_only", test_mbe_monomers_only), &
                  new_unittest("mbe_simple_dimer", test_mbe_simple_dimer), &
                  new_unittest("mbe_sorted_order", test_mbe_sorted_order), &
                  new_unittest("mbe_reverse_order", test_mbe_reverse_order), &
                  new_unittest("mbe_random_order", test_mbe_random_order) &
                  ]
   end subroutine collect_mqc_mbe_tests

   subroutine test_mbe_monomers_only(error)
      !! Test MBE energy with monomers only (nlevel=1)
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), allocatable :: results(:)
      integer, allocatable :: polymers(:, :)
      real(dp) :: total_energy
      integer(int64) :: fragment_count
      integer :: max_level

      ! Three monomers with energies -10.0, -15.0, -20.0
      fragment_count = 3
      max_level = 1
      allocate (polymers(fragment_count, max_level))
      allocate (results(fragment_count))

      ! Monomers in order: [1], [2], [3]
      polymers(1, 1) = 1
      polymers(2, 1) = 2
      polymers(3, 1) = 3

      ! Set monomer energies
      results(1)%energy%scf = -10.0_dp
      results(1)%has_energy = .true.
      results(2)%energy%scf = -15.0_dp
      results(2)%has_energy = .true.
      results(3)%energy%scf = -20.0_dp
      results(3)%has_energy = .true.

      ! Compute MBE energy
      call compute_mbe_energy(polymers, fragment_count, max_level, results, total_energy)

      ! Total should be -10 + (-15) + (-20) = -45
      call check(error, total_energy, -45.0_dp, thr=1.0e-10_dp, &
                 message="MBE monomers only should equal sum of monomer energies")
      if (allocated(error)) goto 100

100   continue
      if (allocated(results)) deallocate (results)
      if (allocated(polymers)) deallocate (polymers)

   end subroutine test_mbe_monomers_only

   subroutine test_mbe_simple_dimer(error)
      !! Test MBE with 2 monomers and 1 dimer (nlevel=2)
      !! E_total = E(1) + E(2) + deltaE(1,2)
      !! where deltaE(1,2) = E(1,2) - E(1) - E(2)
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), allocatable :: results(:)
      integer, allocatable :: polymers(:, :)
      real(dp) :: total_energy
      integer(int64) :: fragment_count
      integer :: max_level

      ! 2 monomers + 1 dimer
      fragment_count = 3
      max_level = 2
      allocate (polymers(fragment_count, max_level))
      allocate (results(fragment_count))
      polymers = 0

      ! Fragment 1: monomer [1]
      polymers(1, 1) = 1
      results(1)%energy%scf = -10.0_dp
      results(1)%has_energy = .true.

      ! Fragment 2: monomer [2]
      polymers(2, 1) = 2
      results(2)%energy%scf = -15.0_dp
      results(2)%has_energy = .true.

      ! Fragment 3: dimer [1, 2]
      polymers(3, 1) = 1
      polymers(3, 2) = 2
      results(3)%energy%scf = -26.0_dp  ! Slightly more stable than sum
      results(3)%has_energy = .true.

      ! Compute MBE energy
      call compute_mbe_energy(polymers, fragment_count, max_level, results, total_energy)

      ! deltaE(1,2) = E(1,2) - E(1) - E(2) = -26 - (-10) - (-15) = -26 + 10 + 15 = -1
      ! Total = E(1) + E(2) + deltaE(1,2) = -10 + (-15) + (-1) = -26
      call check(error, total_energy, -26.0_dp, thr=1.0e-10_dp, &
                 message="MBE with simple dimer")
      if (allocated(error)) goto 200

200   continue
      if (allocated(results)) deallocate (results)
      if (allocated(polymers)) deallocate (polymers)

   end subroutine test_mbe_simple_dimer

   subroutine test_mbe_sorted_order(error)
      !! Test MBE with fragments in size order (monomers first)
      !! 3 monomers + 3 dimers + 1 trimer
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), allocatable :: results(:)
      integer, allocatable :: polymers(:, :)
      real(dp) :: total_energy
      integer(int64) :: fragment_count
      integer :: max_level

      fragment_count = 7  ! 3 monomers + 3 dimers + 1 trimer
      max_level = 3
      allocate (polymers(fragment_count, max_level))
      allocate (results(fragment_count))
      polymers = 0

      ! Monomers: [1], [2], [3]
      polymers(1, 1) = 1
      polymers(2, 1) = 2
      polymers(3, 1) = 3
      results(1)%energy%scf = -10.0_dp
      results(2)%energy%scf = -12.0_dp
      results(3)%energy%scf = -11.0_dp
      results(1:3)%has_energy = .true.

      ! Dimers: [1,2], [1,3], [2,3]
      polymers(4, 1:2) = [1, 2]
      polymers(5, 1:2) = [1, 3]
      polymers(6, 1:2) = [2, 3]
      results(4)%energy%scf = -22.5_dp  ! E(1) + E(2) = -22, so delta = -0.5
      results(5)%energy%scf = -21.3_dp  ! E(1) + E(3) = -21, so delta = -0.3
      results(6)%energy%scf = -23.4_dp  ! E(2) + E(3) = -23, so delta = -0.4
      results(4:6)%has_energy = .true.

      ! Trimer: [1,2,3]
      polymers(7, 1:3) = [1, 2, 3]
      results(7)%energy%scf = -33.5_dp
      results(7)%has_energy = .true.

      call compute_mbe_energy(polymers, fragment_count, max_level, results, total_energy)

      ! 1-body: -10 + -12 + -11 = -33
      ! 2-body deltas: -0.5 + -0.3 + -0.4 = -1.2
      ! 3-body delta: E(1,2,3) - E(1) - E(2) - E(3) - delta(1,2) - delta(1,3) - delta(2,3)
      !             = -33.5 - (-10) - (-12) - (-11) - (-0.5) - (-0.3) - (-0.4)
      !             = -33.5 + 10 + 12 + 11 + 0.5 + 0.3 + 0.4 = 0.7
      ! Total = -33 + (-1.2) + 0.7 = -33.5
      call check(error, total_energy, -33.5_dp, thr=1.0e-10_dp, &
                 message="MBE with fragments in sorted order")
      if (allocated(error)) goto 300

300   continue
      if (allocated(results)) deallocate (results)
      if (allocated(polymers)) deallocate (polymers)

   end subroutine test_mbe_sorted_order

   subroutine test_mbe_reverse_order(error)
      !! Test MBE with fragments in REVERSE size order (trimer first, monomers last)
      !! This tests that internal sorting works correctly
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), allocatable :: results(:)
      integer, allocatable :: polymers(:, :)
      real(dp) :: total_energy
      integer(int64) :: fragment_count
      integer :: max_level

      fragment_count = 7  ! 3 monomers + 3 dimers + 1 trimer
      max_level = 3
      allocate (polymers(fragment_count, max_level))
      allocate (results(fragment_count))
      polymers = 0

      ! REVERSE ORDER: Trimer first
      polymers(1, 1:3) = [1, 2, 3]
      results(1)%energy%scf = -33.5_dp
      results(1)%has_energy = .true.

      ! Dimers: [1,2], [1,3], [2,3]
      polymers(2, 1:2) = [1, 2]
      polymers(3, 1:2) = [1, 3]
      polymers(4, 1:2) = [2, 3]
      results(2)%energy%scf = -22.5_dp
      results(3)%energy%scf = -21.3_dp
      results(4)%energy%scf = -23.4_dp
      results(2:4)%has_energy = .true.

      ! Monomers last: [1], [2], [3]
      polymers(5, 1) = 1
      polymers(6, 1) = 2
      polymers(7, 1) = 3
      results(5)%energy%scf = -10.0_dp
      results(6)%energy%scf = -12.0_dp
      results(7)%energy%scf = -11.0_dp
      results(5:7)%has_energy = .true.

      call compute_mbe_energy(polymers, fragment_count, max_level, results, total_energy)

      ! Should get same answer as sorted order test
      call check(error, total_energy, -33.5_dp, thr=1.0e-10_dp, &
                 message="MBE with fragments in reverse order (tests internal sorting)")
      if (allocated(error)) goto 400

400   continue
      if (allocated(results)) deallocate (results)
      if (allocated(polymers)) deallocate (polymers)

   end subroutine test_mbe_reverse_order

   subroutine test_mbe_random_order(error)
      !! Test MBE with fragments in mixed/random order
      !! Order: [dimer], [monomer], [trimer], [monomer], [dimer], [dimer], [monomer]
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), allocatable :: results(:)
      integer, allocatable :: polymers(:, :)
      real(dp) :: total_energy
      integer(int64) :: fragment_count
      integer :: max_level

      fragment_count = 7
      max_level = 3
      allocate (polymers(fragment_count, max_level))
      allocate (results(fragment_count))
      polymers = 0

      ! Mixed order
      ! Fragment 1: dimer [1,2]
      polymers(1, 1:2) = [1, 2]
      results(1)%energy%scf = -22.5_dp
      results(1)%has_energy = .true.

      ! Fragment 2: monomer [1]
      polymers(2, 1) = 1
      results(2)%energy%scf = -10.0_dp
      results(2)%has_energy = .true.

      ! Fragment 3: trimer [1,2,3]
      polymers(3, 1:3) = [1, 2, 3]
      results(3)%energy%scf = -33.5_dp
      results(3)%has_energy = .true.

      ! Fragment 4: monomer [2]
      polymers(4, 1) = 2
      results(4)%energy%scf = -12.0_dp
      results(4)%has_energy = .true.

      ! Fragment 5: dimer [1,3]
      polymers(5, 1:2) = [1, 3]
      results(5)%energy%scf = -21.3_dp
      results(5)%has_energy = .true.

      ! Fragment 6: dimer [2,3]
      polymers(6, 1:2) = [2, 3]
      results(6)%energy%scf = -23.4_dp
      results(6)%has_energy = .true.

      ! Fragment 7: monomer [3]
      polymers(7, 1) = 3
      results(7)%energy%scf = -11.0_dp
      results(7)%has_energy = .true.

      call compute_mbe_energy(polymers, fragment_count, max_level, results, total_energy)

      ! Should get same answer regardless of input order
      call check(error, total_energy, -33.5_dp, thr=1.0e-10_dp, &
                 message="MBE with fragments in random/mixed order")
      if (allocated(error)) goto 500

500   continue
      if (allocated(results)) deallocate (results)
      if (allocated(polymers)) deallocate (polymers)

   end subroutine test_mbe_random_order

end module test_mqc_mbe

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mbe, only: collect_mqc_mbe_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_mbe", collect_mqc_mbe_tests) &
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
