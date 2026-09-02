!! Unit tests for Lebedev-Laikov angular quadrature
module test_mqc_lebedev
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_lebedev, only: lebedev_available_orders, lebedev_is_available, &
                          lebedev_order_at_least, lebedev_has_negative_weights, &
                          lebedev_grid
   use mqc_physical_constants, only: PI
   implicit none
   private

   public :: collect_mqc_lebedev

   !! Weights are exact rational-ish constants summed in a fixed order, so the
   !! sum-to-one check can be held to a few ulp over 5810 terms.
   real(dp), parameter :: SUM_TOL = 1.0e-13_dp

   !! A grid of degree n integrates spherical harmonics up to degree n exactly.
   !! This is the property that makes it a quadrature rather than a point set.
   real(dp), parameter :: EXACT_TOL = 1.0e-12_dp

contains

   subroutine collect_mqc_lebedev(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("lebedev_orders_are_ascending_and_known", test_orders), &
                  new_unittest("lebedev_availability_query", test_availability), &
                  new_unittest("lebedev_order_at_least_rounds_up", test_order_at_least), &
                  new_unittest("lebedev_negative_weight_orders_declared", test_negative_declared), &
                  new_unittest("lebedev_grid_has_the_promised_size", test_sizes), &
                  new_unittest("lebedev_points_lie_on_the_unit_sphere", test_unit_sphere), &
                  new_unittest("lebedev_weights_sum_to_one", test_weight_sum), &
                  new_unittest("lebedev_integrates_harmonics_exactly", test_harmonics), &
                  new_unittest("lebedev_is_invariant_under_inversion", test_inversion), &
                  new_unittest("lebedev_unknown_order_is_refused", test_unknown_order) &
                  ]
   end subroutine collect_mqc_lebedev

   subroutine test_orders(error)
      !! The published set, ascending, 6 through 5810
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: orders(:)
      integer :: i

      orders = lebedev_available_orders()
      call check(error, size(orders), 32)
      if (allocated(error)) return
      call check(error, orders(1), 6)
      if (allocated(error)) return
      call check(error, orders(size(orders)), 5810)
      if (allocated(error)) return

      do i = 2, size(orders)
         call check(error, orders(i) > orders(i - 1), "orders must ascend")
         if (allocated(error)) return
      end do
   end subroutine test_orders

   subroutine test_availability(error)
      type(error_type), allocatable, intent(out) :: error

      call check(error, lebedev_is_available(302), "302 is a Lebedev order")
      if (allocated(error)) return
      call check(error,.not. lebedev_is_available(300), "300 is not")
      if (allocated(error)) return
      call check(error,.not. lebedev_is_available(0), "0 is not")
      if (allocated(error)) return
      call check(error,.not. lebedev_is_available(-7), "a negative count is not")
   end subroutine test_availability

   subroutine test_order_at_least(error)
      !! Rounds up to an available order, and saturates rather than failing
      type(error_type), allocatable, intent(out) :: error

      call check(error, lebedev_order_at_least(300), 302)
      if (allocated(error)) return
      call check(error, lebedev_order_at_least(302), 302)
      if (allocated(error)) return
      call check(error, lebedev_order_at_least(1), 6)
      if (allocated(error)) return
      ! Asking for more than exists gives the largest rather than an error:
      ! there is nothing finer to offer and refusing would be less useful.
      call check(error, lebedev_order_at_least(99999), 5810)
   end subroutine test_order_at_least

   subroutine test_negative_declared(error)
      !! Orders 74, 230 and 266 really do carry negative weights
      !!
      !! Checked both ways. A grid declared to have them must have them, and a
      !! grid not declared must not -- the first version of this asserted only
      !! that weights were positive, which is false for these three and was the
      !! bug rather than the finding.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: points(:, :), weights(:)
      type(error_t) :: err
      integer, allocatable :: orders(:)
      integer :: i
      logical :: has_negative

      orders = lebedev_available_orders()
      do i = 1, size(orders)
         call lebedev_grid(orders(i), points, weights, err)
         call check(error,.not. err%has_error(), "grid must build: "//err%get_full_trace())
         if (allocated(error)) return

         has_negative = any(weights < 0.0_dp)
         call check(error, has_negative .eqv. lebedev_has_negative_weights(orders(i)), &
                    "negative-weight declaration disagrees with the weights")
         deallocate (points, weights)
         if (allocated(error)) return
      end do

      call check(error, lebedev_has_negative_weights(74), "74 has negative weights")
      if (allocated(error)) return
      call check(error,.not. lebedev_has_negative_weights(302), "302 does not")
   end subroutine test_negative_declared

   subroutine test_sizes(error)
      !! Every order returns exactly as many points as its name says
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: points(:, :), weights(:)
      type(error_t) :: err
      integer, allocatable :: orders(:)
      integer :: i

      orders = lebedev_available_orders()
      do i = 1, size(orders)
         call lebedev_grid(orders(i), points, weights, err)
         call check(error,.not. err%has_error(), "grid must build: "//err%get_full_trace())
         if (allocated(error)) return
         call check(error, size(weights), orders(i))
         if (allocated(error)) return
         call check(error, size(points, 2), orders(i))
         deallocate (points, weights)
         if (allocated(error)) return
      end do
   end subroutine test_sizes

   subroutine test_unit_sphere(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: points(:, :), weights(:)
      type(error_t) :: err
      integer, allocatable :: orders(:)
      integer :: i, k
      real(dp) :: worst

      orders = lebedev_available_orders()
      worst = 0.0_dp
      do i = 1, size(orders)
         call lebedev_grid(orders(i), points, weights, err)
         do k = 1, size(weights)
            worst = max(worst, abs(norm2(points(:, k)) - 1.0_dp))
         end do
         deallocate (points, weights)
      end do
      call check(error, worst < 1.0e-14_dp, "points must lie on the unit sphere")
   end subroutine test_unit_sphere

   subroutine test_weight_sum(error)
      !! Weights are normalised to 1, not to 4*pi
      !!
      !! Which convention holds matters more than it looks: the caller supplies
      !! the 4*pi, and a grid normalised the other way would make every
      !! integral wrong by that factor without failing anything here.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: points(:, :), weights(:)
      type(error_t) :: err
      integer, allocatable :: orders(:)
      integer :: i
      real(dp) :: worst

      orders = lebedev_available_orders()
      worst = 0.0_dp
      do i = 1, size(orders)
         call lebedev_grid(orders(i), points, weights, err)
         worst = max(worst, abs(sum(weights) - 1.0_dp))
         deallocate (points, weights)
      end do
      call check(error, worst < SUM_TOL, "weights must sum to one")
   end subroutine test_weight_sum

   subroutine test_harmonics(error)
      !! The quadrature property, which is the only thing that makes it a grid
      !!
      !! A Lebedev grid of degree n integrates every spherical harmonic up to
      !! degree n exactly. Testing a few low-order polynomials in x, y and z is
      !! the same statement in a form that needs no special functions: over the
      !! unit sphere the average of x, of xy, of x^2 - y^2 and of xyz is zero,
      !! while the average of x^2 + y^2 + z^2 is one by construction.
      !!
      !! This is what a point set that merely lies on a sphere would fail.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: p(:, :), w(:)
      type(error_t) :: err
      integer, allocatable :: orders(:)
      integer :: i
      real(dp) :: x_avg, xy_avg, quad_avg, xyz_avg, radius_avg

      orders = lebedev_available_orders()
      do i = 1, size(orders)
         call lebedev_grid(orders(i), p, w, err)

         x_avg = sum(w*p(1, :))
         xy_avg = sum(w*p(1, :)*p(2, :))
         quad_avg = sum(w*(p(1, :)**2 - p(2, :)**2))
         xyz_avg = sum(w*p(1, :)*p(2, :)*p(3, :))
         radius_avg = sum(w*(p(1, :)**2 + p(2, :)**2 + p(3, :)**2))

         call check(error, abs(x_avg) < EXACT_TOL, "average of x must vanish")
         if (allocated(error)) return
         call check(error, abs(xy_avg) < EXACT_TOL, "average of xy must vanish")
         if (allocated(error)) return
         call check(error, abs(quad_avg) < EXACT_TOL, "average of x^2-y^2 must vanish")
         if (allocated(error)) return
         call check(error, abs(xyz_avg) < EXACT_TOL, "average of xyz must vanish")
         if (allocated(error)) return
         call check(error, abs(radius_avg - 1.0_dp) < EXACT_TOL, &
                    "average of x^2+y^2+z^2 must be one")
         deallocate (p, w)
         if (allocated(error)) return
      end do
   end subroutine test_harmonics

   subroutine test_inversion(error)
      !! Every point has its antipode, with the same weight
      !!
      !! The orbits are octahedral, so the set is closed under r -> -r. This
      !! catches a sign dropped in one branch of the orbit expansion, which a
      !! count and a weight sum would both survive.
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: p(:, :), w(:)
      type(error_t) :: err
      integer :: k, j, n
      logical :: found

      call lebedev_grid(110, p, w, err)
      call check(error,.not. err%has_error(), "grid must build: "//err%get_full_trace())
      if (allocated(error)) return

      n = size(w)
      do k = 1, n
         found = .false.
         do j = 1, n
            if (maxval(abs(p(:, j) + p(:, k))) < 1.0e-12_dp) then
               found = abs(w(j) - w(k)) < 1.0e-14_dp
               exit
            end if
         end do
         call check(error, found, "every point needs its antipode at equal weight")
         if (allocated(error)) return
      end do
   end subroutine test_inversion

   subroutine test_unknown_order(error)
      !! An order that does not exist is refused, not silently rounded
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: points(:, :), weights(:)
      type(error_t) :: err

      call lebedev_grid(300, points, weights, err)
      call check(error, err%has_error(), "300 is not a Lebedev order and must be refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "300") > 0, &
                 "the message should name the order asked for")
   end subroutine test_unknown_order

end module test_mqc_lebedev

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_lebedev, only: collect_mqc_lebedev
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_lebedev", collect_mqc_lebedev)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
