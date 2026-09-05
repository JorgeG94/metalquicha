!! The three pure helpers the batched Fock build's decomposition rests on
module test_mqc_czt_fock_blocking
   !! `pair_degeneracy`, `block_density_max` and `pair_work_order` decide,
   !! respectively, what each quartet is weighted by, what the density screen
   !! compares against, and what order the threads take the work in. The first
   !! two change the answer if they are wrong and do it quietly -- a degeneracy
   !! off by two is a Fock matrix wrong on its diagonal blocks only -- and none
   !! of the three had a test.
   !!
   !! They are tested here rather than through an SCF because each is a pure
   !! function of small integer input: a case that pins the contract exactly
   !! beats one that would have to notice a two-fold error in a converged
   !! energy.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_czt_direct, only: pair_degeneracy, block_density_max, pair_work_order
   implicit none
   private

   public :: collect_czt_fock_blocking

contains

   subroutine collect_czt_fock_blocking(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("degeneracy_counts_shells_not_functions", test_degeneracy), &
                  new_unittest("degeneracy_never_exceeds_eightfold", test_degeneracy_range), &
                  new_unittest("block_density_max_is_the_block_maximum", test_block_max), &
                  new_unittest("block_density_max_is_symmetric", test_block_symmetry), &
                  new_unittest("work_order_is_a_permutation", test_order_permutation), &
                  new_unittest("work_order_puts_the_costliest_first", test_order_descending) &
                  ]
   end subroutine collect_czt_fock_blocking

   subroutine test_degeneracy(error)
      !! The four cases the eightfold symmetry actually distinguishes
      type(error_type), allocatable, intent(out) :: error

      ! All four shells distinct: every permutation is a separate quartet.
      call check(error, pair_degeneracy(1, 2, 3, 4) == 8.0_dp, "all distinct is eightfold")
      if (allocated(error)) return
      ! A bra of one shell with itself already holds both orderings.
      call check(error, pair_degeneracy(1, 1, 3, 4) == 4.0_dp, "s1 == s2 halves it")
      if (allocated(error)) return
      call check(error, pair_degeneracy(1, 2, 3, 3) == 4.0_dp, "s3 == s4 halves it")
      if (allocated(error)) return
      ! Bra and ket the same pair: the bra-ket swap is not a new quartet.
      call check(error, pair_degeneracy(1, 2, 1, 2) == 4.0_dp, "bra == ket halves it")
      if (allocated(error)) return
      ! Everything equal, so nothing to count twice.
      call check(error, pair_degeneracy(1, 1, 1, 1) == 1.0_dp, "one shell throughout is unweighted")
   end subroutine test_degeneracy

   subroutine test_degeneracy_range(error)
      !! Never above eight and never below one, over every shape of four shells
      !!
      !! The factors are independent, so a fifth `if` or a wrong comparison
      !! would show up as a weight off the 1/2/4/8 ladder rather than as a
      !! specific case being wrong.
      type(error_type), allocatable, intent(out) :: error

      integer :: a, b, c, d
      real(dp) :: deg

      do a = 1, 3
         do b = 1, 3
            do c = 1, 3
               do d = 1, 3
                  deg = pair_degeneracy(a, b, c, d)
                  call check(error, deg >= 1.0_dp .and. deg <= 8.0_dp, "weight off the ladder")
                  if (allocated(error)) return
                  call check(error, any(abs(deg - [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp]) < 1.0e-14_dp), &
                             "weight is not a power of two")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_degeneracy_range

   subroutine test_block_max(error)
      !! Each entry is the largest |D| in its shell block
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: density(3, 3)
      real(dp), allocatable :: dsh(:, :)
      integer, parameter :: OFFS(2) = [0, 1]   !! shell 1 is AO 1, shell 2 is AOs 2-3
      integer, parameter :: DIMS(2) = [1, 2]

      ! The negative entry is the largest in magnitude, so a missing `abs`
      ! would take 0.4 instead.
      density = reshape([1.0_dp, 0.2_dp, -0.9_dp, &
                         0.2_dp, 0.3_dp, 0.4_dp, &
                         -0.9_dp, 0.4_dp, 0.5_dp], [3, 3])
      call block_density_max(density, 2, OFFS, DIMS, dsh)

      call check(error, abs(dsh(1, 1) - 1.0_dp) < 1.0e-14_dp, "1x1 block")
      if (allocated(error)) return
      call check(error, abs(dsh(2, 2) - 0.5_dp) < 1.0e-14_dp, "2x2 block")
      if (allocated(error)) return
      call check(error, abs(dsh(1, 2) - 0.9_dp) < 1.0e-14_dp, "off-diagonal takes |D|, not D")
   end subroutine test_block_max

   subroutine test_block_symmetry(error)
      !! `dsh` is symmetric even when the density handed in is not
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: density(3, 3)
      real(dp), allocatable :: dsh(:, :)
      integer, parameter :: OFFS(2) = [0, 1]
      integer, parameter :: DIMS(2) = [1, 2]
      integer :: i, j

      ! Deliberately asymmetric: only the lower triangle is read, and the upper
      ! is filled by reflection, so an implementation that read both would
      ! disagree here.
      density = reshape([1.0_dp, 0.7_dp, 0.6_dp, &
                         0.0_dp, 0.3_dp, 0.4_dp, &
                         0.0_dp, 0.0_dp, 0.5_dp], [3, 3])
      call block_density_max(density, 2, OFFS, DIMS, dsh)

      do i = 1, 2
         do j = 1, 2
            call check(error, abs(dsh(i, j) - dsh(j, i)) < 1.0e-14_dp, "dsh must be symmetric")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_block_symmetry

   subroutine test_order_permutation(error)
      !! Every task appears exactly once
      !!
      !! A task dropped or repeated here is a Fock matrix missing a block or
      !! counting one twice, which no other test in the suite would localise.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: PAIR_I(6) = [1, 2, 2, 3, 3, 3]
      integer, parameter :: PAIR_J(6) = [1, 1, 2, 1, 2, 3]
      integer, parameter :: DIMS(3) = [1, 3, 6]   !! s, p, d
      integer, allocatable :: order(:)
      logical :: seen(6)
      integer :: k

      call pair_work_order(PAIR_I, PAIR_J, DIMS, order)

      call check(error, size(order) == 6, "one entry per task")
      if (allocated(error)) return
      seen = .false.
      do k = 1, size(order)
         call check(error, order(k) >= 1 .and. order(k) <= 6, "index out of range")
         if (allocated(error)) return
         call check(error,.not. seen(order(k)), "task listed twice")
         if (allocated(error)) return
         seen(order(k)) = .true.
      end do
      call check(error, all(seen), "task missing from the order")
   end subroutine test_order_permutation

   subroutine test_order_descending(error)
      !! The costliest task is dispatched first, which is the whole point
      !!
      !! `schedule(dynamic)` hands tasks out in the order given, so the most
      !! expensive one going last leaves every other thread idle behind it.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: PAIR_I(6) = [1, 2, 2, 3, 3, 3]
      integer, parameter :: PAIR_J(6) = [1, 1, 2, 1, 2, 3]
      integer, parameter :: DIMS(3) = [1, 3, 6]
      integer, allocatable :: order(:)
      integer :: k
      integer(kind=8) :: cost(6), cumulative, blk

      call pair_work_order(PAIR_I, PAIR_J, DIMS, order)

      cumulative = 0_8
      do k = 1, 6
         blk = int(DIMS(PAIR_I(k)), 8)*int(DIMS(PAIR_J(k)), 8)
         cumulative = cumulative + blk
         cost(k) = blk*cumulative
      end do

      ! The d-d task is last in the natural order and by far the dearest, so a
      ! sort that ran ascending would put it at the end instead of the front.
      call check(error, order(1) == 6, "the costliest task must be dispatched first")
      if (allocated(error)) return
      do k = 2, 6
         call check(error, cost(order(k - 1)) >= cost(order(k)), "order is not descending in cost")
         if (allocated(error)) return
      end do
   end subroutine test_order_descending

end module test_mqc_czt_fock_blocking

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_fock_blocking, only: collect_czt_fock_blocking
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_fock_blocking", collect_czt_fock_blocking)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
