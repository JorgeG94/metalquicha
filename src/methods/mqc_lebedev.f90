module mqc_lebedev
   !! Lebedev-Laikov angular quadrature on the unit sphere
   !!
   !! A Lebedev grid is a union of orbits of the octahedral group: pick a seed
   !! point, apply every coordinate permutation and sign flip, and every point
   !! produced carries the same weight. So the whole grid -- 5810 points at the
   !! largest order -- is stored as 1287 (code, a, b, weight) records in
   !! `mqc_lebedev_data`, and expanded here.
   !!
   !! The six orbit types differ only in their seed:
   !!
   !!    code 0   (1, 0, 0)                        6 points
   !!    code 1   (0, s, s)      s = sqrt(1/2)    12 points
   !!    code 2   (t, t, t)      t = sqrt(1/3)     8 points
   !!    code 3   (a, a, b)      b = sqrt(1-2a^2) 24 points
   !!    code 4   (a, b, 0)      b = sqrt(1-a^2)  24 points
   !!    code 5   (a, b, c)      c = sqrt(1-a^2-b^2)  48 points
   !!
   !! Rather than tabulate each orbit's points, `expand_orbit` generates them
   !! from the seed and discards duplicates -- which is where the smaller counts
   !! come from, a seed with a repeated or zero coordinate having a smaller
   !! orbit than the 48 of a general point. The expected size is asserted, so a
   !! seed that collapsed further than it should cannot pass silently.
   !!
   !! Weights are normalised to sum to 1: the 4*pi of the sphere's surface is
   !! not included and must be supplied by the caller if it is wanted.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_lebedev_data, only: LEBEDEV_ORDERS, ORDER_OFFSET, ORBIT_CODE, &
                               ORBIT_A, ORBIT_B, ORBIT_W
   implicit none
   private

   public :: lebedev_available_orders  !! The orders that exist
   public :: lebedev_is_available    !! Whether a given point count is one of them
   public :: lebedev_order_at_least  !! Smallest available order >= a request
   public :: lebedev_has_negative_weights  !! Whether an order carries negative weights
   public :: lebedev_grid            !! Points and weights for an order

   !! Number of orbit types, codes 0 to N_ORBIT_CODES-1
   integer, parameter :: N_ORBIT_CODES = 6

   !! Largest orbit: a seed with three distinct non-zero coordinates reaches
   !! all 6 permutations in all 8 sign combinations.
   integer, parameter :: MAX_ORBIT_POINTS = 48

   !! Permutations of three coordinates
   integer, parameter :: N_PERMUTATIONS = 6

   !! Spatial dimensions
   integer, parameter :: N_DIM = 3

   !! Points in the orbit of each code, before duplicates are removed
   integer, parameter :: ORBIT_SIZE(0:N_ORBIT_CODES - 1) = [6, 12, 8, 24, 24, MAX_ORBIT_POINTS]

   !! Orders whose weights are not all positive.
   !!
   !! These are valid quadratures -- they integrate to the stated accuracy, and
   !! their weights still sum to 1 -- but a negative weight can drive a
   !! numerically integrated density negative in a low-density region, so DFT
   !! grids conventionally avoid them. Checked against PySCF: the same three
   !! orders, with 8, 6 and 18 negative weights respectively.
   integer, parameter :: NEGATIVE_WEIGHT_ORDERS(3) = [74, 230, 266]

   !! Two generated points closer than this are the same point. Seeds differ in
   !! their parameters by far more than this, so no distinct pair can merge.
   real(dp), parameter :: DUPLICATE_TOL = 1.0e-12_dp

contains

   pure function lebedev_available_orders() result(orders)
      !! Available Lebedev orders, ascending
      integer :: orders(size(LEBEDEV_ORDERS))
      orders = LEBEDEV_ORDERS
   end function lebedev_available_orders

   pure function lebedev_is_available(n_points) result(available)
      !! Whether `n_points` is one of the Lebedev orders
      integer, intent(in) :: n_points
      logical :: available
      available = any(LEBEDEV_ORDERS == n_points)
   end function lebedev_is_available

   pure function lebedev_order_at_least(n_points) result(order)
      !! Smallest available order with at least `n_points` points
      !!
      !! Returns the largest order, 5810, if more than that is asked for; there
      !! is nothing finer to offer and refusing would be less useful.
      integer, intent(in) :: n_points
      integer :: order
      integer :: i

      order = LEBEDEV_ORDERS(size(LEBEDEV_ORDERS))
      do i = 1, size(LEBEDEV_ORDERS)
         if (LEBEDEV_ORDERS(i) >= n_points) then
            order = LEBEDEV_ORDERS(i)
            exit
         end if
      end do
   end function lebedev_order_at_least

   pure function lebedev_has_negative_weights(n_points) result(negative)
      !! Whether this order carries negative weights
      integer, intent(in) :: n_points
      logical :: negative
      negative = any(NEGATIVE_WEIGHT_ORDERS == n_points)
   end function lebedev_has_negative_weights

   subroutine lebedev_grid(n_points, points, weights, error)
      !! Unit-sphere quadrature points and weights for one Lebedev order
      integer, intent(in) :: n_points                     !! Must be an available order
      real(dp), allocatable, intent(out) :: points(:, :)  !! (3, n_points), on the unit sphere
      real(dp), allocatable, intent(out) :: weights(:)    !! Sums to 1
      type(error_t), intent(inout) :: error

      real(dp) :: orbit(N_DIM, MAX_ORBIT_POINTS)
      integer :: order_index, record, first, last, n_written, n_orbit, i
      character(len=64) :: message

      order_index = 0
      do i = 1, size(LEBEDEV_ORDERS)
         if (LEBEDEV_ORDERS(i) == n_points) then
            order_index = i
            exit
         end if
      end do

      if (order_index == 0) then
         write (message, "(a,i0,a)") "no Lebedev grid of ", n_points, " points"
         call error%set(ERROR_VALIDATION, trim(message))
         return
      end if

      allocate (points(3, n_points))
      allocate (weights(n_points))

      first = ORDER_OFFSET(order_index)
      last = ORDER_OFFSET(order_index + 1) - 1
      n_written = 0

      do record = first, last
         call expand_orbit(ORBIT_CODE(record), ORBIT_A(record), ORBIT_B(record), &
                           orbit, n_orbit)

         ! An orbit that came out the wrong size means the seed collapsed in a
         ! way the code does not describe, and every later point would be
         ! misplaced. Stop rather than return a grid that is quietly short.
         if (n_orbit /= ORBIT_SIZE(ORBIT_CODE(record))) then
            write (message, "(a,i0,a,i0,a,i0)") "Lebedev orbit ", record, &
               " gave ", n_orbit, " points, expected ", ORBIT_SIZE(ORBIT_CODE(record))
            call error%set(ERROR_VALIDATION, trim(message))
            return
         end if

         points(:, n_written + 1:n_written + n_orbit) = orbit(:, 1:n_orbit)
         weights(n_written + 1:n_written + n_orbit) = ORBIT_W(record)
         n_written = n_written + n_orbit
      end do

      if (n_written /= n_points) then
         write (message, "(a,i0,a,i0)") "Lebedev order filled ", n_written, &
            " of ", n_points
         call error%set(ERROR_VALIDATION, trim(message))
      end if
   end subroutine lebedev_grid

   pure subroutine expand_orbit(code, a, b, points, n_points)
      !! Every distinct image of one orbit's seed under the octahedral group
      integer, intent(in) :: code            !! Orbit type, 0 to 5
      real(dp), intent(in) :: a, b           !! Seed parameters; unused for codes 0 to 2
      real(dp), intent(out) :: points(N_DIM, MAX_ORBIT_POINTS)  !! Distinct images, first `n_points` valid
      integer, intent(out) :: n_points

      !! The six permutations of three coordinates
      integer, parameter :: PERM(N_DIM, N_PERMUTATIONS) = reshape([1, 2, 3, 1, 3, 2, 2, 1, 3, &
                                                                   2, 3, 1, 3, 1, 2, 3, 2, 1], [N_DIM, N_PERMUTATIONS])
      real(dp) :: seed(3), candidate(3)
      integer :: p, sx, sy, sz, i
      logical :: seen

      seed = orbit_seed(code, a, b)
      n_points = 0

      do p = 1, N_PERMUTATIONS
         do sx = -1, 1, 2
            do sy = -1, 1, 2
               do sz = -1, 1, 2
                  candidate(1) = real(sx, dp)*seed(PERM(1, p))
                  candidate(2) = real(sy, dp)*seed(PERM(2, p))
                  candidate(3) = real(sz, dp)*seed(PERM(3, p))

                  ! A seed with repeated or zero coordinates reaches the same
                  ! point by several routes; keep it once.
                  seen = .false.
                  do i = 1, n_points
                     if (maxval(abs(points(:, i) - candidate)) < DUPLICATE_TOL) then
                        seen = .true.
                        exit
                     end if
                  end do

                  if (.not. seen) then
                     n_points = n_points + 1
                     points(:, n_points) = candidate
                  end if
               end do
            end do
         end do
      end do
   end subroutine expand_orbit

   pure function orbit_seed(code, a, b) result(seed)
      !! The seed point of an orbit, on the unit sphere by construction
      integer, intent(in) :: code
      real(dp), intent(in) :: a, b
      real(dp) :: seed(3)

      select case (code)
      case (0)
         seed = [1.0_dp, 0.0_dp, 0.0_dp]
      case (1)
         seed = [0.0_dp, sqrt(0.5_dp), sqrt(0.5_dp)]
      case (2)
         seed = [sqrt(1.0_dp/3.0_dp), sqrt(1.0_dp/3.0_dp), sqrt(1.0_dp/3.0_dp)]
      case (3)
         seed = [a, a, sqrt(max(0.0_dp, 1.0_dp - 2.0_dp*a*a))]
      case (4)
         seed = [a, sqrt(max(0.0_dp, 1.0_dp - a*a)), 0.0_dp]
      case (5)
         seed = [a, b, sqrt(max(0.0_dp, 1.0_dp - a*a - b*b))]
      case default
         seed = 0.0_dp
      end select
   end function orbit_seed

end module mqc_lebedev
