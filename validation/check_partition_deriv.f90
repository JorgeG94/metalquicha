!! The Becke partition weight derivatives, on their own
program check_partition_deriv
   !! Checks `becke_partition_derivatives` against central differences of
   !! `becke_partition_weights`, at fixed grid points.
   !!
   !! Checked here rather than only through the assembled DFT gradient, for the
   !! reason the density-fitting derivative integrals were: a mistake in this
   !! term produces a gradient of the right magnitude that disagrees with finite
   !! differences for reasons that could be anywhere in the exchange-correlation
   !! contribution. Differencing the partition function itself puts the error
   !! where it happened.
   !!
   !! The partition is a smooth, purely geometric function of the nuclear
   !! positions, so a central difference of it is good to about 1e-10 -- far
   !! tighter than anything the assembled gradient can be checked to, which is
   !! what makes this the sharper test of the two.
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_dft_partition, only: becke_partition_weights, becke_partition_derivatives, &
                                PARTITION_BECKE, PARTITION_STRATMANN, &
                                ADJUST_NONE, ADJUST_TREUTLER
   implicit none

   integer :: n_bad

   n_bad = 0

   call one_case("water, Becke + Treutler", PARTITION_BECKE, ADJUST_TREUTLER, n_bad)
   call one_case("water, Becke, no adjustment", PARTITION_BECKE, ADJUST_NONE, n_bad)
   call one_case("water, Stratmann + Treutler", PARTITION_STRATMANN, ADJUST_TREUTLER, n_bad)

   write (*, "(a)") ""
   if (n_bad == 0) then
      write (*, "(a)") "all partition derivative checks passed"
   else
      write (*, "(a,i0,a)") "FAILED: ", n_bad, " check(s)"
      error stop 1
   end if

contains

   subroutine one_case(label, scheme, adjust, n_bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: scheme, adjust
      integer, intent(inout) :: n_bad

      integer, parameter :: N_ATOMS = 3
      integer, parameter :: N_POINTS = 7
      real(dp), parameter :: STEP = 1.0e-5_dp

      real(dp) :: coords(3, N_ATOMS), points(3, N_POINTS)
      real(dp) :: shifted(3, N_ATOMS)
      integer :: numbers(N_ATOMS), owner(N_POINTS)
      real(dp) :: analytic(3, N_ATOMS, N_POINTS)
      real(dp) :: w_plus(N_POINTS), w_minus(N_POINTS)
      real(dp) :: numeric, worst, diff
      integer :: ia, ic, k
      type(error_t) :: error

      numbers = [8, 1, 1]
      coords = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                        0.0_dp, -1.4308_dp, 1.1078_dp, &
                        0.0_dp, 1.4308_dp, 1.1078_dp], [3, N_ATOMS])

      ! Points scattered through the region where the partition actually
      ! varies. A point far outside every cell has a weight of one or zero and
      ! a derivative of zero, which would pass without testing anything.
      points = reshape([0.3_dp, 0.2_dp, 0.4_dp, &
                        0.0_dp, -0.7_dp, 0.6_dp, &
                        0.1_dp, 0.9_dp, 0.5_dp, &
                        -0.4_dp, 0.0_dp, 1.3_dp, &
                        0.6_dp, -1.1_dp, 0.9_dp, &
                        0.2_dp, 1.2_dp, 1.4_dp, &
                        0.0_dp, 0.0_dp, 2.1_dp], [3, N_POINTS])
      owner = [1, 2, 3, 1, 2, 3, 1]

      call becke_partition_derivatives(points, coords, numbers, owner, scheme, adjust, &
                                       analytic, error)
      if (error%has_error()) then
         write (*, "(a,a)") "FAIL: ", error%get_message()
         n_bad = n_bad + 1
         return
      end if

      worst = 0.0_dp
      do ia = 1, N_ATOMS
         do ic = 1, 3
            shifted = coords
            shifted(ic, ia) = coords(ic, ia) + STEP
            call becke_partition_weights(points, shifted, numbers, owner, scheme, adjust, &
                                         w_plus, error)
            if (error%has_error()) then
               write (*, "(a,a)") "FAIL: ", error%get_message()
               n_bad = n_bad + 1
               return
            end if

            shifted = coords
            shifted(ic, ia) = coords(ic, ia) - STEP
            call becke_partition_weights(points, shifted, numbers, owner, scheme, adjust, &
                                         w_minus, error)
            if (error%has_error()) then
               write (*, "(a,a)") "FAIL: ", error%get_message()
               n_bad = n_bad + 1
               return
            end if

            do k = 1, N_POINTS
               numeric = (w_plus(k) - w_minus(k))/(2.0_dp*STEP)
               diff = abs(numeric - analytic(ic, ia, k))
               worst = max(worst, diff)
            end do
         end do
      end do

      write (*, "(a,a,a,es12.4)") "== ", label, "   largest deviation: ", worst
      if (worst > 1.0e-8_dp) then
         write (*, "(a)") "  FAIL: the analytic partition derivative is wrong"
         n_bad = n_bad + 1
      end if
   end subroutine one_case

end program check_partition_deriv
