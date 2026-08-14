!! Unit tests for NWChem angular pruning
module test_mqc_dft_prune
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_lebedev, only: lebedev_is_available
   use mqc_dft_radial, only: treutler_ahlrichs_radial
   use mqc_dft_prune, only: prune_angular_orders, prune_scheme_name, &
                            PRUNE_NONE, PRUNE_NWCHEM
   implicit none
   private

   public :: collect_mqc_dft_prune

contains

   subroutine collect_mqc_dft_prune(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("prune_none_is_a_constant_order", test_none), &
                  new_unittest("prune_emits_only_real_lebedev_orders", test_real_orders), &
                  new_unittest("prune_never_exceeds_the_target", test_bounded), &
                  new_unittest("prune_reduces_the_total", test_reduces), &
                  new_unittest("prune_uses_the_target_somewhere", test_uses_target), &
                  new_unittest("prune_leaves_small_grids_alone", test_small_untouched), &
                  new_unittest("prune_scheme_names", test_names), &
                  new_unittest("prune_rejects_bad_input", test_bad_input) &
                  ]
   end subroutine collect_mqc_dft_prune

   subroutine mesh(z, n, r)
      integer, intent(in) :: z, n
      real(dp), allocatable, intent(out) :: r(:)
      real(dp), allocatable :: dr(:)
      type(error_t) :: err
      call treutler_ahlrichs_radial(n, z, r, dr, err)
   end subroutine mesh

   subroutine test_none(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:)
      integer, allocatable :: orders(:)
      type(error_t) :: err

      call mesh(8, 40, r)
      allocate (orders(size(r)))
      call prune_angular_orders(PRUNE_NONE, 8, r, 302, orders, err)
      call check(error,.not. err%has_error(), "no pruning must be allowed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, all(orders == 302), "without pruning every shell keeps the target")
   end subroutine test_none

   subroutine test_real_orders(error)
      !! Every emitted order must be one lebedev_grid can actually build
      !!
      !! The scheme indexes into a sequence, so an off-by-one lands on a
      !! neighbouring order rather than an invalid one -- which would surface
      !! much later, as a grid of the wrong size, not as an error here.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:)
      integer, allocatable :: orders(:)
      type(error_t) :: err
      integer :: targets(4), z_list(3)
      integer :: i, j, k

      targets = [50, 194, 302, 590]
      z_list = [1, 8, 79]

      do i = 1, size(z_list)
         do j = 1, size(targets)
            call mesh(z_list(i), 50, r)
            if (allocated(orders)) deallocate (orders)
            allocate (orders(size(r)))
            call prune_angular_orders(PRUNE_NWCHEM, z_list(i), r, targets(j), orders, err)
            call check(error,.not. err%has_error(), "pruning must succeed: "//err%get_full_trace())
            if (allocated(error)) return
            do k = 1, size(orders)
               call check(error, lebedev_is_available(orders(k)), &
                          "pruning must emit a real Lebedev order")
               if (allocated(error)) return
            end do
            deallocate (r)
         end do
      end do
   end subroutine test_real_orders

   subroutine test_bounded(error)
      !! Pruning only ever removes points, never adds them
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:)
      integer, allocatable :: orders(:)
      type(error_t) :: err

      call mesh(8, 60, r)
      allocate (orders(size(r)))
      call prune_angular_orders(PRUNE_NWCHEM, 8, r, 302, orders, err)
      call check(error, all(orders <= 302), "no shell may exceed the target order")
   end subroutine test_bounded

   subroutine test_reduces(error)
      !! The point of the scheme: fewer points than the unpruned grid
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:)
      integer, allocatable :: orders(:)
      type(error_t) :: err
      integer :: pruned, unpruned

      call mesh(8, 75, r)
      allocate (orders(size(r)))
      call prune_angular_orders(PRUNE_NWCHEM, 8, r, 302, orders, err)
      pruned = sum(orders)
      unpruned = size(r)*302
      call check(error, pruned < unpruned, "pruning must remove points")
      if (allocated(error)) return
      ! Not a token reduction either.
      call check(error, real(pruned, dp) < 0.8_dp*real(unpruned, dp), &
                 "pruning should remove a useful fraction, not a handful")
   end subroutine test_reduces

   subroutine test_uses_target(error)
      !! The valence zone keeps the full order it was asked for
      !!
      !! A scheme that reduced everywhere would integrate the valence region
      !! worse than requested, which is the region that matters.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:)
      integer, allocatable :: orders(:)
      type(error_t) :: err

      call mesh(8, 75, r)
      allocate (orders(size(r)))
      call prune_angular_orders(PRUNE_NWCHEM, 8, r, 302, orders, err)
      call check(error, any(orders == 302), "some shell must keep the target order")
   end subroutine test_uses_target

   subroutine test_small_untouched(error)
      !! Below 50 points there is nothing worth taking away
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:)
      integer, allocatable :: orders(:)
      type(error_t) :: err

      call mesh(1, 30, r)
      allocate (orders(size(r)))
      call prune_angular_orders(PRUNE_NWCHEM, 1, r, 38, orders, err)
      call check(error,.not. err%has_error(), "a small target must be accepted: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, all(orders == 38), "a grid below the floor is left alone")
   end subroutine test_small_untouched

   subroutine test_names(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, prune_scheme_name(PRUNE_NONE), "none")
      if (allocated(error)) return
      call check(error, prune_scheme_name(PRUNE_NWCHEM), "nwchem")
      if (allocated(error)) return
      call check(error, prune_scheme_name(42), "unknown")
   end subroutine test_names

   subroutine test_bad_input(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: r(:)
      integer, allocatable :: orders(:)
      type(error_t) :: err

      call mesh(8, 20, r)

      allocate (orders(5))
      call prune_angular_orders(PRUNE_NWCHEM, 8, r, 302, orders, err)
      call check(error, err%has_error(), "a short orders array must be refused")
      if (allocated(error)) return

      call err%clear()
      deallocate (orders)
      allocate (orders(size(r)))
      call prune_angular_orders(99, 8, r, 302, orders, err)
      call check(error, err%has_error(), "an unknown scheme must be refused")
      if (allocated(error)) return

      call err%clear()
      call prune_angular_orders(PRUNE_NWCHEM, 8, r, 301, orders, err)
      call check(error, err%has_error(), "a target that is not a Lebedev order must be refused")
   end subroutine test_bad_input

end module test_mqc_dft_prune

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_dft_prune, only: collect_mqc_dft_prune
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_dft_prune", collect_mqc_dft_prune)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
