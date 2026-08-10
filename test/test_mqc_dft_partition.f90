!! Unit tests for the Becke partition of space into atomic cells
module test_mqc_dft_partition
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_dft_partition, only: becke_partition_weights, partition_scheme_name, &
                                PARTITION_BECKE, PARTITION_STRATMANN, &
                                ADJUST_NONE, ADJUST_BECKE, ADJUST_TREUTLER
   implicit none
   private

   public :: collect_mqc_dft_partition

   integer, parameter :: N_DIM = 3

contains

   subroutine collect_mqc_dft_partition(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("partition_weights_sum_to_one", test_unity), &
                  new_unittest("partition_of_one_atom_is_everything", test_single_atom), &
                  new_unittest("partition_splits_identical_atoms_evenly", test_symmetry), &
                  new_unittest("partition_favours_the_nearer_atom", test_nearer), &
                  new_unittest("partition_size_adjustment_favours_the_larger", test_adjust), &
                  new_unittest("partition_stratmann_saturates_far_away", test_stratmann_cutoff), &
                  new_unittest("partition_scheme_names", test_names), &
                  new_unittest("partition_rejects_bad_input", test_bad_input) &
                  ]
   end subroutine collect_mqc_dft_partition

   subroutine water(coords, z)
      real(dp), allocatable, intent(out) :: coords(:, :)
      integer, allocatable, intent(out) :: z(:)
      allocate (coords(N_DIM, 3), z(3))
      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                        0.0_dp, -1.4308_dp, 1.1078_dp, &
                        0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3])
      z = [8, 1, 1]
   end subroutine water

   subroutine test_unity(error)
      !! The defining property: the cells partition space, so summing the
      !! weight of every atom at any point must give exactly one.
      !!
      !! Strong for the price. It fails if the pair loop is asymmetric, if a
      !! size-adjustment shift has the wrong sign, or if an atom is dropped
      !! from a product -- none of which a single-atom weight would reveal.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: coords(:, :), probes(:, :), w(:), unity(:)
      integer, allocatable :: z(:), owner(:)
      type(error_t) :: err
      integer :: schemes(2), adjusts(3)
      integer :: s, a, ia, i, j, k, n
      real(dp) :: worst

      call water(coords, z)
      schemes = [PARTITION_BECKE, PARTITION_STRATMANN]
      adjusts = [ADJUST_NONE, ADJUST_BECKE, ADJUST_TREUTLER]

      allocate (probes(N_DIM, 125))
      n = 0
      do i = 0, 4
         do j = 0, 4
            do k = 0, 4
               n = n + 1
               probes(:, n) = [-4.0_dp + 2.0_dp*i, -4.0_dp + 2.0_dp*j, -4.0_dp + 2.0_dp*k]
            end do
         end do
      end do
      allocate (w(n), owner(n), unity(n))

      worst = 0.0_dp
      do s = 1, 2
         do a = 1, 3
            unity = 0.0_dp
            do ia = 1, 3
               owner = ia
               call becke_partition_weights(probes, coords, z, owner, &
                                            schemes(s), adjusts(a), w, err)
               call check(error,.not. err%has_error(), "partition must build")
               if (allocated(error)) return
               unity = unity + w
            end do
            worst = max(worst, maxval(abs(unity - 1.0_dp)))
         end do
      end do
      call check(error, worst < 1.0e-13_dp, "weights must sum to one everywhere")
   end subroutine test_unity

   subroutine test_single_atom(error)
      !! One atom owns all of space; the pair loop has nothing to do
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: coords(N_DIM, 1), probes(N_DIM, 3), w(3)
      integer :: z(1), owner(3)
      type(error_t) :: err

      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 1])
      z = [8]
      probes = reshape([0.1_dp, 0.0_dp, 0.0_dp, &
                        5.0_dp, 0.0_dp, 0.0_dp, &
                        0.0_dp, 0.0_dp, -20.0_dp], [N_DIM, 3])
      owner = 1
      call becke_partition_weights(probes, coords, z, owner, PARTITION_BECKE, &
                                   ADJUST_BECKE, w, err)
      call check(error,.not. err%has_error(), "one atom must be allowed")
      if (allocated(error)) return
      call check(error, all(abs(w - 1.0_dp) < 1.0e-15_dp), "one atom owns everything")
   end subroutine test_single_atom

   subroutine test_symmetry(error)
      !! Two identical atoms split their midpoint exactly in half
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: coords(N_DIM, 2), probes(N_DIM, 1), w(1)
      integer :: z(2), owner(1)
      type(error_t) :: err

      coords = reshape([-2.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 2])
      z = [8, 8]
      probes = reshape([0.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 1])

      owner = 1
      call becke_partition_weights(probes, coords, z, owner, PARTITION_BECKE, &
                                   ADJUST_BECKE, w, err)
      call check(error, abs(w(1) - 0.5_dp) < 1.0e-14_dp, &
                 "identical atoms must split the midpoint evenly")
   end subroutine test_symmetry

   subroutine test_nearer(error)
      !! A point beside one nucleus belongs almost entirely to it
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: coords(N_DIM, 2), probes(N_DIM, 1), w(1)
      integer :: z(2), owner(1)
      type(error_t) :: err

      coords = reshape([-2.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 2])
      z = [8, 8]
      probes = reshape([-1.95_dp, 0.0_dp, 0.0_dp], [N_DIM, 1])

      owner = 1
      call becke_partition_weights(probes, coords, z, owner, PARTITION_BECKE, &
                                   ADJUST_BECKE, w, err)
      call check(error, w(1) > 0.99_dp, "a point at a nucleus belongs to that atom")
   end subroutine test_nearer

   subroutine test_adjust(error)
      !! The size adjustment moves the boundary toward the smaller atom
      !!
      !! With no adjustment the midpoint of any pair splits evenly whatever the
      !! elements are. With it, oxygen takes more than half from hydrogen --
      !! which is the whole point of the correction, and its sign is the easiest
      !! thing to get backwards.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: coords(N_DIM, 2), probes(N_DIM, 1), w(1)
      integer :: z(2), owner(1)
      type(error_t) :: err
      real(dp) :: plain, adjusted

      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp, 4.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 2])
      z = [8, 1]
      probes = reshape([2.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 1])
      owner = 1

      call becke_partition_weights(probes, coords, z, owner, PARTITION_BECKE, &
                                   ADJUST_NONE, w, err)
      plain = w(1)
      call becke_partition_weights(probes, coords, z, owner, PARTITION_BECKE, &
                                   ADJUST_BECKE, w, err)
      adjusted = w(1)

      call check(error, abs(plain - 0.5_dp) < 1.0e-14_dp, &
                 "without adjustment the midpoint splits evenly")
      if (allocated(error)) return
      call check(error, adjusted > plain, &
                 "oxygen must take more than half from hydrogen")
   end subroutine test_adjust

   subroutine test_stratmann_cutoff(error)
      !! Stratmann reaches exactly 0 and 1, which is what makes it screenable
      !!
      !! Becke's cutoff is smooth and never quite reaches either, so a point
      !! deep inside one cell has a tiny but non-zero weight for the other.
      !! Stratmann's is exactly 1, and that difference is the reason it exists.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: coords(N_DIM, 2), probes(N_DIM, 1), w(1)
      integer :: z(2), owner(1)
      type(error_t) :: err

      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp, 30.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 2])
      z = [8, 8]
      probes = reshape([0.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 1])
      owner = 1

      call becke_partition_weights(probes, coords, z, owner, PARTITION_STRATMANN, &
                                   ADJUST_NONE, w, err)
      call check(error, w(1) == 1.0_dp, "Stratmann must saturate to exactly one")
   end subroutine test_stratmann_cutoff

   subroutine test_names(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, partition_scheme_name(PARTITION_BECKE), "Becke")
      if (allocated(error)) return
      call check(error, partition_scheme_name(PARTITION_STRATMANN), "Stratmann")
      if (allocated(error)) return
      call check(error, partition_scheme_name(-1), "unknown")
   end subroutine test_names

   subroutine test_bad_input(error)
      !! Mismatched arrays and out-of-range owners are refused
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: coords(N_DIM, 2), probes(N_DIM, 2), w(2)
      integer :: owner(2)
      type(error_t) :: err

      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp], [N_DIM, 2])
      probes = coords
      owner = 1

      call becke_partition_weights(probes, coords, [8], owner, PARTITION_BECKE, &
                                   ADJUST_BECKE, w, err)
      call check(error, err%has_error(), "a short atomic_numbers must be refused")
      if (allocated(error)) return

      call err%clear()
      owner = [1, 99]
      call becke_partition_weights(probes, coords, [8, 8], owner, PARTITION_BECKE, &
                                   ADJUST_BECKE, w, err)
      call check(error, err%has_error(), "an owner past the atom list must be refused")
      if (allocated(error)) return

      call err%clear()
      owner = 1
      call becke_partition_weights(probes, coords, [8, 8], owner, -5, &
                                   ADJUST_BECKE, w, err)
      call check(error, err%has_error(), "an unknown scheme must be refused")
   end subroutine test_bad_input

end module test_mqc_dft_partition

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_dft_partition, only: collect_mqc_dft_partition
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_dft_partition", collect_mqc_dft_partition)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
