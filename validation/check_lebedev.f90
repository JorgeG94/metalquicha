program check_lebedev
   !! Dump every Lebedev order so it can be checked against an outside reference
   !!
   !! Writes one line per point -- order, x, y, z, weight -- at full precision.
   !! The self-consistent properties (weights summing to 1, points on the unit
   !! sphere) are checked here; agreement with PySCF's grids is checked by
   !! compare_lebedev.py, which reads this output.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_lebedev, only: lebedev_available_orders, lebedev_grid, &
                          lebedev_has_negative_weights
   implicit none

   real(dp), allocatable :: points(:, :), weights(:)
   type(error_t) :: error
   integer, allocatable :: orders(:)
   integer :: unit, i, k, n_bad
   real(dp) :: weight_sum, radius_dev, worst_weight, worst_radius

   orders = lebedev_available_orders()
   worst_weight = 0.0_dp
   worst_radius = 0.0_dp
   n_bad = 0

   open (newunit=unit, file="lebedev_points.txt", status="replace", action="write")

   do i = 1, size(orders)
      call lebedev_grid(orders(i), points, weights, error)
      if (error%has_error()) then
         write (*, "(a,i0,a,a)") "FAIL order ", orders(i), ": ", error%get_message()
         n_bad = n_bad + 1
         cycle
      end if

      weight_sum = sum(weights)
      radius_dev = 0.0_dp
      do k = 1, size(weights)
         radius_dev = max(radius_dev, abs(norm2(points(:, k)) - 1.0_dp))
         write (unit, "(i0,4(1x,es24.16e3))") orders(i), points(:, k), weights(k)
      end do

      worst_weight = max(worst_weight, abs(weight_sum - 1.0_dp))
      worst_radius = max(worst_radius, radius_dev)

      ! Orders 74, 230 and 266 genuinely have negative weights -- this is a
      ! property of those Lebedev grids, not a transcription error, and PySCF
      ! carries the same ones. Report them, since a DFT grid usually wants to
      ! avoid them, but do not treat them as a failure.
      if (any(weights < 0.0_dp)) then
         if (.not. lebedev_has_negative_weights(orders(i))) then
            write (*, "(a,i0,a)") "FAIL order ", orders(i), &
               " has negative weights but is not declared as such"
            n_bad = n_bad + 1
         else
            write (*, "(a,i0,a,i0,a)") "note: order ", orders(i), " has ", &
               count(weights < 0.0_dp), " negative weights, as expected"
         end if
      else if (lebedev_has_negative_weights(orders(i))) then
         write (*, "(a,i0,a)") "FAIL order ", orders(i), &
            " is declared to have negative weights but does not"
         n_bad = n_bad + 1
      end if

      deallocate (points, weights)
   end do

   close (unit)

   write (*, "(a,i0,a)") "orders written: ", size(orders), ""
   write (*, "(a,es10.3)") "worst |sum(w) - 1|      : ", worst_weight
   write (*, "(a,es10.3)") "worst ||r|| - 1|        : ", worst_radius

   if (n_bad > 0 .or. worst_weight > 1.0e-13_dp .or. worst_radius > 1.0e-13_dp) then
      write (*, "(a)") "FAILED"
      stop 1
   end if
   write (*, "(a)") "self-consistency OK -- now run compare_lebedev.py"
end program check_lebedev
