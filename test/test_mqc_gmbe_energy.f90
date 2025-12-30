module test_mqc_gmbe_energy
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_mbe, only: compute_gmbe_energy
   use mqc_result_types, only: calculation_result_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_gmbe_energy_tests

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_gmbe_energy_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("gmbe_no_intersections", test_gmbe_no_intersections), &
                  new_unittest("gmbe_single_intersection", test_gmbe_single_intersection), &
                  new_unittest("gmbe_multiple_intersections", test_gmbe_multiple_intersections) &
                  ]
   end subroutine collect_mqc_gmbe_energy_tests

   subroutine test_gmbe_no_intersections(error)
      !! Test GMBE energy with no overlapping fragments (should equal sum of monomers)
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), allocatable :: monomer_results(:)
      type(calculation_result_t), allocatable :: intersection_results(:)
      integer, allocatable :: monomers(:), intersection_pairs(:, :)
      real(dp) :: total_energy
      integer :: n_monomers, n_intersections

      ! Two non-overlapping monomers with energies -10.0 and -15.0
      n_monomers = 2
      n_intersections = 0
      allocate (monomers(n_monomers))
      allocate (monomer_results(n_monomers))
      monomers = [1, 2]

      ! Set monomer energies
      monomer_results(1)%energy%scf = -10.0_dp
      monomer_results(1)%has_energy = .true.
      monomer_results(2)%energy%scf = -15.0_dp
      monomer_results(2)%has_energy = .true.

      ! No intersections
      allocate (intersection_results(0))
      allocate (intersection_pairs(2, 0))

      ! Compute GMBE energy
      call compute_gmbe_energy(monomers, n_monomers, monomer_results, &
                               n_intersections, intersection_results, &
                               intersection_pairs, total_energy)

      ! Total should be -10 + (-15) = -25
      call check(error, total_energy, -25.0_dp, thr=1.0e-10_dp, &
                 message="GMBE with no intersections should equal sum of monomers")
      if (allocated(error)) goto 100

100   continue
      if (allocated(monomer_results)) deallocate (monomer_results)
      if (allocated(intersection_results)) deallocate (intersection_results)
      if (allocated(monomers)) deallocate (monomers)
      if (allocated(intersection_pairs)) deallocate (intersection_pairs)

   end subroutine test_gmbe_no_intersections

   subroutine test_gmbe_single_intersection(error)
      !! Test GMBE energy with single intersection
      !! E_total = E(monomer1) + E(monomer2) - E(intersection)
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), allocatable :: monomer_results(:)
      type(calculation_result_t), allocatable :: intersection_results(:)
      integer, allocatable :: monomers(:), intersection_pairs(:, :)
      real(dp) :: total_energy
      integer :: n_monomers, n_intersections

      ! Two overlapping monomers
      n_monomers = 2
      n_intersections = 1
      allocate (monomers(n_monomers))
      allocate (monomer_results(n_monomers))
      allocate (intersection_results(n_intersections))
      allocate (intersection_pairs(2, n_intersections))
      monomers = [1, 2]

      ! Set monomer energies: gly0 = -100.0, gly1 = -95.0
      monomer_results(1)%energy%scf = -100.0_dp
      monomer_results(1)%has_energy = .true.
      monomer_results(2)%energy%scf = -95.0_dp
      monomer_results(2)%has_energy = .true.

      ! Set intersection energy: backbone overlap = -8.0
      intersection_results(1)%energy%scf = -8.0_dp
      intersection_results(1)%has_energy = .true.

      ! Intersection is between monomers 1 and 2
      intersection_pairs(:, 1) = [1, 2]

      ! Compute GMBE energy
      call compute_gmbe_energy(monomers, n_monomers, monomer_results, &
                               n_intersections, intersection_results, &
                               intersection_pairs, total_energy)

      ! Total = -100 + (-95) - (-8) = -100 - 95 + 8 = -187
      call check(error, total_energy, -187.0_dp, thr=1.0e-10_dp, &
                 message="GMBE with one intersection should be sum(monomers) - sum(intersections)")
      if (allocated(error)) goto 200

200   continue
      if (allocated(monomer_results)) deallocate (monomer_results)
      if (allocated(intersection_results)) deallocate (intersection_results)
      if (allocated(monomers)) deallocate (monomers)
      if (allocated(intersection_pairs)) deallocate (intersection_pairs)

   end subroutine test_gmbe_single_intersection

   subroutine test_gmbe_multiple_intersections(error)
      !! Test GMBE energy with multiple intersections (tripeptide-like chain)
      !! E_total = E(gly0) + E(gly1) + E(gly2) - E(int01) - E(int12)
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), allocatable :: monomer_results(:)
      type(calculation_result_t), allocatable :: intersection_results(:)
      integer, allocatable :: monomers(:), intersection_pairs(:, :)
      real(dp) :: total_energy
      integer :: n_monomers, n_intersections

      ! Three overlapping monomers in a chain
      n_monomers = 3
      n_intersections = 2  ! int(0,1) and int(1,2)
      allocate (monomers(n_monomers))
      allocate (monomer_results(n_monomers))
      allocate (intersection_results(n_intersections))
      allocate (intersection_pairs(2, n_intersections))
      monomers = [1, 2, 3]

      ! Set monomer energies
      monomer_results(1)%energy%scf = -100.0_dp  ! gly0
      monomer_results(1)%has_energy = .true.
      monomer_results(2)%energy%scf = -95.0_dp   ! gly1
      monomer_results(2)%has_energy = .true.
      monomer_results(3)%energy%scf = -98.0_dp   ! gly2
      monomer_results(3)%has_energy = .true.

      ! Set intersection energies
      intersection_results(1)%energy%scf = -8.0_dp   ! int(0,1)
      intersection_results(1)%has_energy = .true.
      intersection_results(2)%energy%scf = -7.5_dp   ! int(1,2)
      intersection_results(2)%has_energy = .true.

      ! Intersection pairs
      intersection_pairs(:, 1) = [1, 2]  ! gly0-gly1
      intersection_pairs(:, 2) = [2, 3]  ! gly1-gly2

      ! Compute GMBE energy
      call compute_gmbe_energy(monomers, n_monomers, monomer_results, &
                               n_intersections, intersection_results, &
                               intersection_pairs, total_energy)

      ! Total = -100 + (-95) + (-98) - (-8) - (-7.5)
      !       = -100 - 95 - 98 + 8 + 7.5
      !       = -277.5
      call check(error, total_energy, -277.5_dp, thr=1.0e-10_dp, &
                 message="GMBE with multiple intersections in chain")
      if (allocated(error)) goto 300

300   continue
      if (allocated(monomer_results)) deallocate (monomer_results)
      if (allocated(intersection_results)) deallocate (intersection_results)
      if (allocated(monomers)) deallocate (monomers)
      if (allocated(intersection_pairs)) deallocate (intersection_pairs)

   end subroutine test_gmbe_multiple_intersections

end module test_mqc_gmbe_energy

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_gmbe_energy, only: collect_mqc_gmbe_energy_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_gmbe_energy", collect_mqc_gmbe_energy_tests) &
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
