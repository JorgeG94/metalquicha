module test_mqc_scf_gdm
   !! The numerical core of geometric direct minimisation (mqc_scf_gdm): the
   !! L-BFGS two-loop recursion against a dense BFGS inverse Hessian built by
   !! the textbook update, the curvature guard, convergence on a convex
   !! quadratic with exact line search, and the preconditioner, step-control
   !! and descent helpers. No SCF, no reference program: every expected value
   !! is computed here by a different (dense) route.
   use pic_types, only: dp
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_error, only: error_t
   use mqc_scf_gdm, only: scf_gdm_history_t, scf_gdm_preconditioner_build, scf_gdm_step_judge, &
                          scf_gdm_step_limit, scf_gdm_is_descent, scf_gdm_scale_min, &
                          scf_gdm_curvature_tol, scf_gdm_descent_margin
   implicit none
   private
   public :: collect_mqc_scf_gdm_tests

   real(dp) :: nan_value, inf_value

contains

   subroutine collect_mqc_scf_gdm_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
      inf_value = ieee_value(1.0_dp, ieee_positive_inf)
      testsuite = [ &
                  new_unittest("empty history is preconditioned steepest descent", test_empty), &
                  new_unittest("two-loop equals dense BFGS", test_dense), &
                  new_unittest("wrap-around keeps the last m pairs", test_wrap), &
                  new_unittest("curvature guard", test_guard), &
                  new_unittest("quadratic converges, diagonal preconditioner helps", test_quadratic), &
                  new_unittest("reset and release", test_reset_release), &
                  new_unittest("preconditioner values and index order", test_precond), &
                  new_unittest("step control", test_step_control), &
                  new_unittest("breakdown returns scale 1", test_breakdown_scale), &
                  new_unittest("curvature guard bounds the direction", test_guard_blowup), &
                  new_unittest("step limit", test_step_limit), &
                  new_unittest("misuse is reported, not rejected", test_misuse), &
                  new_unittest("storage is not initialised", test_no_touch), &
                  new_unittest("descent check", test_descent) &
                  ]
   end subroutine collect_mqc_scf_gdm_tests

   subroutine make_spd(n, lowest, highest, a)
      !! A = Q D Q^T with Q the Householder reflector of a fixed vector and D
      !! geometrically spaced from `lowest` to `highest`: SPD, dense, with a
      !! diagonal that differs from the eigenvalues.
      integer, intent(in) :: n
      real(dp), intent(in) :: lowest, highest
      real(dp), intent(out) :: a(n, n)
      real(dp) :: q(n, n), u(n), dmat(n, n)
      integer :: i, j
      do i = 1, n
         u(i) = 1.0_dp + 0.37_dp*real(i, dp) + 0.11_dp*real(mod(7*i, 5), dp)
      end do
      dmat = 0.0_dp
      do i = 1, n
         dmat(i, i) = lowest*(highest/lowest)**(real(mod(5*i, n), dp)/real(max(n - 1, 1), dp))
      end do
      q = 0.0_dp
      do j = 1, n
         do i = 1, n
            q(i, j) = -2.0_dp*u(i)*u(j)/dot_product(u, u)
         end do
         q(j, j) = q(j, j) + 1.0_dp
      end do
      a = matmul(q, matmul(dmat, transpose(q)))
      a = 0.5_dp*(a + transpose(a))
   end subroutine make_spd

   subroutine dense_bfgs_apply(h0_inv, s, y, g, d)
      !! d = -H g with H from the textbook BFGS inverse-Hessian update,
      !! H <- (I - rho s y^T) H (I - rho y s^T) + rho s s^T, applied for the
      !! columns of s, y in order, starting from diag(h0_inv).
      real(dp), intent(in) :: h0_inv(:), s(:, :), y(:, :), g(:)
      real(dp), intent(out) :: d(size(g))
      real(dp) :: h(size(g), size(g)), v(size(g), size(g)), rho
      integer :: k, i, n
      n = size(g)
      h = 0.0_dp
      do i = 1, n
         h(i, i) = h0_inv(i)
      end do
      do k = 1, size(s, 2)
         rho = 1.0_dp/dot_product(s(:, k), y(:, k))
         ! v = I - rho s y^T : element (i,j) = delta - rho s_i y_j
         v = -rho*spread(s(:, k), 2, n)*spread(y(:, k), 1, n)
         do i = 1, n
            v(i, i) = v(i, i) + 1.0_dp
         end do
         h = matmul(v, matmul(h, transpose(v))) &
             + rho*spread(s(:, k), 2, n)*spread(s(:, k), 1, n)
      end do
      d = -matmul(h, g)
   end subroutine dense_bfgs_apply

   subroutine make_pairs(a, count, s, y)
      !! Deterministic distinct steps s_k and y_k = A s_k.
      real(dp), intent(in) :: a(:, :)
      integer, intent(in) :: count
      real(dp), intent(out) :: s(size(a, 1), count), y(size(a, 1), count)
      integer :: k, i
      do k = 1, count
         do i = 1, size(a, 1)
            s(i, k) = sin(1.3_dp*real(i, dp) + 2.1_dp*real(k, dp)) + 0.2_dp*real(k, dp)/real(i, dp)
         end do
      end do
      y = matmul(a, s)
   end subroutine make_pairs

   real(dp) function rel_diff(x, ref) result(r)
      real(dp), intent(in) :: x(:), ref(:)
      r = sqrt(sum((x - ref)**2))/sqrt(sum(ref**2))
   end function rel_diff

   subroutine compare_to_dense(n, m, npush, error)
      !! Push `npush` pairs into a history of capacity m; the direction must
      !! equal the dense BFGS built from only the last min(npush, m) pairs.
      integer, intent(in) :: n, m, npush
      type(error_type), allocatable, intent(out) :: error
      type(scf_gdm_history_t) :: hist
      type(error_t) :: err
      real(dp) :: a(n, n), s(n, npush), y(n, npush), h0_inv(n), g(n), d(n), ref(n)
      integer :: k, first, i
      logical :: ok

      call make_spd(n, 1.0_dp, 50.0_dp, a)
      call make_pairs(a, npush, s, y)
      do i = 1, n
         h0_inv(i) = 1.0_dp/(1.0_dp + 0.5_dp*real(i, dp))
         g(i) = cos(0.7_dp*real(i, dp)) + 0.3_dp
      end do
      call hist%init(n, m, err)
      do k = 1, npush
         call hist%push(s(:, k), y(:, k), ok, err)
         call check(error, ok, "SPD pair accepted")
         if (allocated(error)) then
            call hist%release()
            return
         end if
      end do
      first = max(1, npush - m + 1)
      call check(error, hist%npairs(), min(npush, m), "stored pair count")
      if (.not. allocated(error)) then
         call hist%direction(g, h0_inv, d, err)
         call dense_bfgs_apply(h0_inv, s(:, first:npush), y(:, first:npush), g, ref)
         call check(error, rel_diff(d, ref) < 1.0e-12_dp, "two-loop direction equals dense BFGS")
      end if
      call hist%release()
   end subroutine compare_to_dense

   subroutine test_empty(error)
      type(error_type), allocatable, intent(out) :: error
      type(scf_gdm_history_t) :: hist
      type(error_t) :: err
      real(dp) :: g(5), h0_inv(5), d(5)
      integer :: i
      do i = 1, 5
         g(i) = 0.4_dp*real(i, dp) - 1.1_dp
         h0_inv(i) = 0.25_dp*real(i, dp)
      end do
      call hist%init(5, 4, err)
      call hist%direction(g, h0_inv, d, err)
      call check(error, all(d == -h0_inv*g), "d == -h0_inv*g exactly")
      call hist%release()
   end subroutine test_empty

   subroutine test_dense(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_to_dense(6, 5, 1, error)
      if (allocated(error)) return
      call compare_to_dense(6, 5, 3, error)
      if (allocated(error)) return
      call compare_to_dense(6, 5, 5, error)
   end subroutine test_dense

   subroutine test_wrap(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_to_dense(6, 3, 4, error)
      if (allocated(error)) return
      call compare_to_dense(6, 3, 8, error)
      if (allocated(error)) return
      call compare_to_dense(6, 4, 9, error)
   end subroutine test_wrap

   subroutine test_guard(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: n = 6
      type(scf_gdm_history_t) :: hist
      type(error_t) :: err
      real(dp) :: a(n, n), s(n, 2), y(n, 2), h0_inv(n), g(n), d_before(n), d_after(n), bad(n)
      logical :: ok, rejected(6), kept_above
      integer :: i, count_after
      real(dp) :: weak(n)

      call make_spd(n, 1.0_dp, 50.0_dp, a)
      call make_pairs(a, 2, s, y)
      do i = 1, n
         h0_inv(i) = 1.0_dp/(1.0_dp + 0.5_dp*real(i, dp))
         g(i) = cos(0.7_dp*real(i, dp)) + 0.3_dp
      end do
      call hist%init(n, 4, err)
      call hist%push(s(:, 1), y(:, 1), ok, err)
      call hist%push(s(:, 2), y(:, 2), ok, err)
      call hist%direction(g, h0_inv, d_before, err)

      ! Every pair below is an ordinary rejection; none touches the history.
      ! y = -s: s.y < 0
      call hist%push(s(:, 1), -s(:, 1), ok, err)
      rejected(1) = .not. ok
      ! y orthogonal to s: s.y = 0
      bad = 0.0_dp
      bad(1) = 1.0_dp
      call hist%push(bad, [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], ok, err)
      rejected(2) = .not. ok
      ! NaN and Inf pairs
      bad = s(:, 1)
      bad(3) = nan_value
      call hist%push(bad, y(:, 1), ok, err)
      rejected(3) = .not. ok
      call hist%push(s(:, 1), bad, ok, err)
      rejected(4) = .not. ok
      bad = y(:, 1)
      bad(2) = inf_value
      call hist%push(s(:, 1), bad, ok, err)
      rejected(5) = .not. ok
      ! Positive but noisy curvature: cosine 0.5*curvature_tol
      bad = 0.0_dp
      bad(1) = 1.0_dp
      weak = 0.0_dp
      weak(1) = 0.5_dp*scf_gdm_curvature_tol
      weak(2) = sqrt(1.0_dp - weak(1)**2)
      call hist%push(bad, weak, ok, err)
      rejected(6) = .not. ok
      count_after = hist%npairs()
      call hist%direction(g, h0_inv, d_after, err)
      ! ... and just above it is a legitimate pair
      weak(1) = 2.0_dp*scf_gdm_curvature_tol
      weak(2) = sqrt(1.0_dp - weak(1)**2)
      call hist%push(bad, weak, ok, err)
      kept_above = ok .and. hist%npairs() == count_after + 1
      call hist%release()

      call check(error, rejected(1), "negative curvature rejected")
      if (allocated(error)) return
      call check(error, rejected(2), "zero curvature rejected")
      if (allocated(error)) return
      call check(error, rejected(3), "NaN step rejected")
      if (allocated(error)) return
      call check(error, rejected(4), "NaN gradient change rejected")
      if (allocated(error)) return
      call check(error, rejected(5), "Inf gradient change rejected")
      if (allocated(error)) return
      call check(error, rejected(6), "cosine below the curvature tolerance rejected")
      if (allocated(error)) return
      call check(error, kept_above, "cosine above the curvature tolerance stored")
      if (allocated(error)) return
      call check(error,.not. err%has_error(), "rejections are not errors")
      if (allocated(error)) return
      call check(error, count_after == 2, "count unchanged after rejections")
      if (allocated(error)) return
      call check(error, all(d_after == d_before), "direction identical after rejected pushes")
      if (allocated(error)) return
      call check(error, scf_gdm_is_descent(g, d_after), "still a descent direction")
   end subroutine test_guard

   subroutine test_quadratic(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: n = 20
      real(dp) :: a(n, n), diagonal(n), ones(n), scaling(n)
      integer :: iterations_plain, iterations_diag, i

      ! A = S A0 S: a mildly conditioned A0 under a diagonal scaling S that
      ! spans a factor 30, the situation an orbital-energy preconditioner is for.
      call make_spd(n, 1.0_dp, 10.0_dp, a)
      do i = 1, n
         scaling(i) = 30.0_dp**(real(mod(3*i, n), dp)/real(n - 1, dp))
      end do
      a = spread(scaling, 2, n)*a*spread(scaling, 1, n)
      ones = 1.0_dp
      diagonal = [(1.0_dp/a(i, i), i=1, n)]
      call minimise(a, ones, iterations_plain)
      call minimise(a, diagonal, iterations_diag)
      call check(error, iterations_plain >= 0 .and. iterations_plain <= n + 4, &
                 "identity preconditioner converges within n + 4 iterations")
      if (allocated(error)) return
      call check(error, iterations_diag >= 0 .and. iterations_diag <= n + 4, &
                 "diagonal preconditioner converges within n + 4 iterations")
      if (allocated(error)) return
      ! Measured: 20 iterations with the identity, 15 with the diagonal.
      call check(error, iterations_diag < iterations_plain, &
                 "exact diagonal needs fewer iterations than the identity")
   end subroutine test_quadratic

   subroutine minimise(a, h0_inv, iterations)
      !! Minimise 1/2 x^T A x - b^T x by L-BFGS (m = n) with exact line
      !! search; iterations = -1 if |g| < 1e-10 is not reached in 3n steps.
      real(dp), intent(in) :: a(:, :), h0_inv(:)
      integer, intent(out) :: iterations
      type(scf_gdm_history_t) :: hist
      type(error_t) :: err
      real(dp) :: x(size(a, 1)), b(size(a, 1)), g(size(a, 1)), gnew(size(a, 1))
      real(dp) :: d(size(a, 1)), s(size(a, 1)), ad(size(a, 1)), step
      integer :: n, i, it
      logical :: ok

      n = size(a, 1)
      do i = 1, n
         b(i) = 1.0_dp + 0.1_dp*real(i, dp)
      end do
      x = 0.0_dp
      g = matmul(a, x) - b
      call hist%init(n, n, err)
      iterations = -1
      do it = 1, 3*n
         if (sqrt(dot_product(g, g)) < 1.0e-10_dp) then
            iterations = it - 1
            exit
         end if
         call hist%direction(g, h0_inv, d, err)
         if (.not. scf_gdm_is_descent(g, d)) then
            call hist%reset()
            d = -h0_inv*g
         end if
         ad = matmul(a, d)
         step = -dot_product(g, d)/dot_product(d, ad)
         s = step*d
         x = x + s
         gnew = matmul(a, x) - b
         call hist%push(s, gnew - g, ok, err)
         g = gnew
      end do
      call hist%release()
   end subroutine minimise

   subroutine test_reset_release(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: n = 6
      type(scf_gdm_history_t) :: hist
      type(error_t) :: err
      real(dp) :: a(n, n), s(n, 3), y(n, 3), h0_inv(n), g(n), d(n)
      logical :: ok
      integer :: k

      call make_spd(n, 1.0_dp, 50.0_dp, a)
      call make_pairs(a, 3, s, y)
      h0_inv = 0.5_dp
      g = 1.0_dp
      call hist%init(n, 3, err)
      do k = 1, 3
         call hist%push(s(:, k), y(:, k), ok, err)
      end do
      call hist%reset()
      call check(error, hist%npairs(), 0, "reset empties the history")
      if (allocated(error)) return
      call hist%direction(g, h0_inv, d, err)
      call check(error, all(d == -h0_inv*g), "after reset: steepest descent")
      if (allocated(error)) return
      call hist%release()
      call check(error, hist%npairs(), 0, "released history is empty")
      if (allocated(error)) return
      call hist%init(n, 2, err)
      do k = 1, 3
         call hist%push(s(:, k), y(:, k), ok, err)
      end do
      call check(error, hist%npairs(), 2, "re-init with a new capacity works")
      if (.not. allocated(error)) then
         call compare_to_dense(n, 2, 3, error)
      end if
      call hist%release()
   end subroutine test_reset_release

   subroutine test_precond(error)
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      real(dp) :: e_occ(2), e_vir(3), h0_inv(6), expected(6)
      integer :: i, a

      e_occ = [-1.0_dp, -0.5_dp]
      e_vir = [0.25_dp, 0.75_dp, 2.0_dp]
      call scf_gdm_preconditioner_build(e_occ, e_vir, 4.0_dp, 0.1_dp, h0_inv, err)
      do i = 1, 2
         do a = 1, 3
            expected((i - 1)*3 + a) = 1.0_dp/(4.0_dp*(e_vir(a) - e_occ(i)))
         end do
      end do
      call check(error, all(abs(h0_inv - expected) < 1.0e-15_dp), "scale 4 values, a fast and i slow")
      if (allocated(error)) return
      ! element (i=2, a=1) is the fourth entry: 1/(4*(0.25+0.5))
      call check(error, abs(h0_inv(4) - 1.0_dp/3.0_dp) < 1.0e-15_dp, "(i=2, a=1) sits at position 4")
      if (allocated(error)) return
      call scf_gdm_preconditioner_build(e_occ, e_vir, 2.0_dp, 0.1_dp, h0_inv, err)
      call check(error, abs(h0_inv(3) - 1.0_dp/(2.0_dp*3.0_dp)) < 1.0e-15_dp, "scale 2 value")
      if (allocated(error)) return
      ! floor: a gap of 0.02 and a negative gap both give 1/floor
      call scf_gdm_preconditioner_build([0.0_dp, 1.0_dp], [0.02_dp, 3.0_dp], 4.0_dp, 0.5_dp, h0_inv(1:4), err)
      call check(error, abs(h0_inv(1) - 2.0_dp) < 1.0e-15_dp, "small gap floored")
      if (allocated(error)) return
      call check(error, abs(h0_inv(2) - 1.0_dp/12.0_dp) < 1.0e-15_dp, "large gap untouched")
      if (allocated(error)) return
      call check(error, abs(h0_inv(3) - 2.0_dp) < 1.0e-15_dp, "negative gap floored")
      if (allocated(error)) return
      call check(error,.not. err%has_error(), "valid arguments are not an error")
      if (allocated(error)) return
      ! The floor acts on the Hessian estimate factor*gap, not on the gap:
      ! gap 0.1 with factor 4 is 0.4, below a floor of 0.5.
      call scf_gdm_preconditioner_build([0.0_dp], [0.1_dp], 4.0_dp, 0.5_dp, h0_inv(1:1), err)
      call check(error, abs(h0_inv(1) - 2.0_dp) < 1.0e-15_dp, "floor is in units of the Hessian estimate")
      if (allocated(error)) return
      ! Misuse: non-positive floor, non-positive factor, wrong output length.
      h0_inv = -7.0_dp
      call scf_gdm_preconditioner_build(e_occ, e_vir, 4.0_dp, 0.0_dp, h0_inv, err)
      call check(error, err%has_error() .and. all(h0_inv == -7.0_dp), "floor 0 reported, output untouched")
      if (allocated(error)) return
      call err%clear()
      call scf_gdm_preconditioner_build(e_occ, e_vir, 4.0_dp, nan_value, h0_inv, err)
      call check(error, err%has_error(), "NaN floor reported")
      if (allocated(error)) return
      call err%clear()
      call scf_gdm_preconditioner_build(e_occ, e_vir, -4.0_dp, 0.1_dp, h0_inv, err)
      call check(error, err%has_error(), "negative factor reported")
      if (allocated(error)) return
      call err%clear()
      call scf_gdm_preconditioner_build(e_occ, e_vir, 4.0_dp, 0.1_dp, h0_inv(1:5), err)
      call check(error, err%has_error(), "wrong output length reported")
   end subroutine test_precond

   subroutine test_step_control(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: scale
      logical :: accepted, broken
      integer :: count

      scale = 1.0_dp
      call scf_gdm_step_judge(-10.0_dp, -10.5_dp, 1.0e-10_dp, scale, accepted, broken)
      call check(error, accepted .and. .not. broken .and. scale == 1.0_dp, "decrease accepted")
      if (allocated(error)) return
      call scf_gdm_step_judge(-10.0_dp, -10.0_dp + 5.0e-11_dp, 1.0e-10_dp, scale, accepted, broken)
      call check(error, accepted, "increase within slack accepted")
      if (allocated(error)) return
      call scf_gdm_step_judge(-10.0_dp, -9.0_dp, 1.0e-10_dp, scale, accepted, broken)
      call check(error,.not. accepted .and. .not. broken .and. scale == 0.5_dp, "increase rejected, scale halved")
      if (allocated(error)) return
      call scf_gdm_step_judge(-10.0_dp, -11.0_dp, 1.0e-10_dp, scale, accepted, broken)
      call check(error, accepted .and. scale == 1.0_dp, "acceptance restores the scale")
      if (allocated(error)) return
      call scf_gdm_step_judge(-10.0_dp, nan_value, 1.0e-10_dp, scale, accepted, broken)
      call check(error,.not. accepted .and. scale == 0.5_dp, "NaN energy rejected")
      if (allocated(error)) return
      call scf_gdm_step_judge(-10.0_dp, -inf_value, 1.0e-10_dp, scale, accepted, broken)
      call check(error,.not. accepted, "-Inf energy rejected")
      if (allocated(error)) return
      ! repeated rejection must report breakdown after a finite number of halvings
      scale = 1.0_dp
      count = 0
      broken = .false.
      do while (.not. broken .and. count < 100)
         call scf_gdm_step_judge(-10.0_dp, -9.0_dp, 0.0_dp, scale, accepted, broken)
         count = count + 1
      end do
      call check(error, broken .and. count == 11 .and. scale == 1.0_dp, &
                 "breakdown after 11 rejections")
   end subroutine test_step_control

   subroutine test_breakdown_scale(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: step_scale
      logical :: accepted, broken
      integer :: count

      step_scale = 1.0_dp
      broken = .false.
      count = 0
      do while (.not. broken .and. count < 100)
         call scf_gdm_step_judge(-10.0_dp, -9.0_dp, 0.0_dp, step_scale, accepted, broken)
         count = count + 1
      end do
      call check(error, broken .and. count == 11 .and. step_scale == 1.0_dp, &
                 "breakdown on the 11th rejection returns scale 1")
      if (allocated(error)) return
      ! A restart with that scale gets a full set of retries, not one.
      call scf_gdm_step_judge(-10.0_dp, -9.0_dp, 0.0_dp, step_scale, accepted, broken)
      call check(error,.not. accepted .and. .not. broken .and. step_scale == 0.5_dp, &
                 "first rejection after a restart does not break down")
   end subroutine test_breakdown_scale

   subroutine test_guard_blowup(error)
      !! With H0 = I and |s| = |y| a pair of cosine c amplifies the gradient by
      !! about 1/c^2 (measured here), so the guard's own cosine bounds it at
      !! about 1e8 (the caller's step cap handles the rest) and an unguarded
      !! 2e-10 would give 1e19.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: n = 4
      type(scf_gdm_history_t) :: hist
      type(error_t) :: err
      real(dp) :: s(n), y(n), g(n), h0_inv(n), d(n)
      real(dp), parameter :: cosine = 2.0_dp*scf_gdm_curvature_tol
      logical :: ok, ok_noisy

      s = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      g = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      h0_inv = 1.0_dp
      call hist%init(n, 3, err)
      y = [cosine, sqrt(1.0_dp - cosine**2), 0.0_dp, 0.0_dp]
      call hist%push(s, y, ok, err)
      call hist%direction(g, h0_inv, d, err)
      call check(error, ok .and. sqrt(sum(d**2)) < 2.0_dp/cosine**2, "legitimate pair: amplification below 2/c^2")
      if (allocated(error)) then
         call hist%release()
         return
      end if
      call hist%reset()
      y = [2.0e-10_dp, sqrt(1.0_dp - 4.0e-20_dp), 0.0_dp, 0.0_dp]
      call hist%push(s, y, ok_noisy, err)
      call hist%direction(g, h0_inv, d, err)
      call hist%release()
      call check(error,.not. ok_noisy .and. sqrt(sum(d**2)) == sqrt(sum(g**2)), &
                 "cosine 2e-10 is not stored: steepest descent")
   end subroutine test_guard_blowup

   subroutine test_step_limit(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: factor, d(3)
      logical :: valid

      d = [0.5_dp, -0.25_dp, 0.1_dp]
      call scf_gdm_step_limit(d, 1.0_dp, factor, valid)
      call check(error, valid .and. factor == 1.0_dp, "fits: factor 1")
      if (allocated(error)) return
      call scf_gdm_step_limit(d, 0.5_dp, factor, valid)
      call check(error, valid .and. factor == 1.0_dp, "exact boundary: factor 1")
      if (allocated(error)) return
      d = [4.0_dp, -8.0_dp, 1.0_dp]
      call scf_gdm_step_limit(d, 1.0_dp, factor, valid)
      call check(error, valid .and. abs(factor - 0.125_dp) < 1.0e-15_dp .and. maxval(abs(factor*d)) <= 1.0_dp, &
                 "needs scaling: factor 1/8, largest element capped")
      if (allocated(error)) return
      ! m/L rounds up here: (m/L)*L = m + 1 ulp without the nearest() correction
      d = [2.9493945741755203_dp, -1.0_dp, 0.1_dp]
      call scf_gdm_step_limit(d, 0.47990491802444557_dp, factor, valid)
      call check(error, valid .and. factor > 0.0_dp .and. factor < 1.0_dp .and. &
                 maxval(abs(factor*d)) <= 0.47990491802444557_dp, "hard cap holds after rounding")
      if (allocated(error)) return
      call scf_gdm_step_limit([0.0_dp, 0.0_dp], 1.0_dp, factor, valid)
      call check(error, valid .and. factor == 1.0_dp, "zero vector: factor 1")
      if (allocated(error)) return
      call scf_gdm_step_limit(d(1:0), 1.0_dp, factor, valid)
      call check(error, valid .and. factor == 1.0_dp, "empty vector: factor 1")
      if (allocated(error)) return
      call scf_gdm_step_limit([1.0_dp, nan_value], 1.0_dp, factor, valid)
      call check(error,.not. valid .and. factor == 1.0_dp, "NaN reported, not scaled")
      if (allocated(error)) return
      call scf_gdm_step_limit([inf_value, 1.0_dp], 1.0_dp, factor, valid)
      call check(error,.not. valid .and. factor == 1.0_dp, "Inf reported, not scaled")
      if (allocated(error)) return
      call scf_gdm_step_limit([1.0_dp], 0.0_dp, factor, valid)
      call check(error,.not. valid, "non-positive cap reported")
   end subroutine test_step_limit

   subroutine test_misuse(error)
      type(error_type), allocatable, intent(out) :: error
      type(scf_gdm_history_t) :: hist
      type(error_t) :: err
      real(dp) :: d(3), keep(3)
      logical :: ok

      ! Uninitialised history.
      ok = .true.
      call hist%push([1.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 0.0_dp, 0.0_dp], ok, err)
      call check(error, err%has_error() .and. .not. ok, "push on an uninitialised history is an error")
      if (allocated(error)) return
      call err%clear()
      call hist%init(0, 3, err)
      call check(error, err%has_error(), "init with n = 0 is an error")
      if (allocated(error)) return
      call err%clear()
      call hist%init(3, 2, err)
      call check(error,.not. err%has_error(), "valid init is not")
      if (allocated(error)) return
      ! Wrong lengths: an error, not an ordinary rejection.
      call hist%push([1.0_dp, 0.0_dp], [1.0_dp, 0.0_dp], ok, err)
      call check(error, err%has_error() .and. .not. ok, "push with a wrong length is an error")
      if (allocated(error)) return
      call err%clear()
      call hist%push([1.0_dp, 0.0_dp, 0.0_dp], [-1.0_dp, 0.0_dp, 0.0_dp], ok, err)
      call check(error,.not. err%has_error() .and. .not. ok, "a curvature rejection is not an error")
      if (allocated(error)) return
      keep = 7.0_dp
      d = keep
      call hist%direction([1.0_dp, 1.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], d, err)
      call check(error, err%has_error() .and. all(d == keep), "direction with a short gradient is an error")
      if (allocated(error)) return
      call err%clear()
      call hist%direction([1.0_dp, 1.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], d, err)
      call check(error, err%has_error(), "direction with a short h0_inv is an error")
      if (allocated(error)) return
      call err%clear()
      call hist%direction([1.0_dp, 1.0_dp, 1.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], d(1:2), err)
      call check(error, err%has_error(), "direction with a short output is an error")
      call hist%release()
   end subroutine test_misuse

   subroutine test_no_touch(error)
      !! `init` allocates and nothing more: a history that is only partly
      !! filled works, and the direction reads only stored slots.
      type(error_type), allocatable, intent(out) :: error
      type(scf_gdm_history_t) :: hist
      type(error_t) :: err
      real(dp) :: d(3)
      logical :: ok
      call hist%init(3, 50, err)
      call hist%push([1.0_dp, 0.0_dp, 0.0_dp], [2.0_dp, 0.0_dp, 0.0_dp], ok, err)
      call hist%direction([1.0_dp, 1.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 1.0_dp], d, err)
      call hist%release()
      ! H = diag(0.5, 1, 1) after one pair along x with curvature 2.
      call check(error, ok .and. abs(d(1) + 0.5_dp) < 1.0e-14_dp .and. abs(d(2) + 1.0_dp) < 1.0e-14_dp, &
                 "partly filled history gives the BFGS direction")
   end subroutine test_no_touch

   subroutine test_descent(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: g(3)
      g = [1.0_dp, -2.0_dp, 0.5_dp]
      call check(error, scf_gdm_is_descent(g, -g), "-g is descent")
      if (allocated(error)) return
      call check(error,.not. scf_gdm_is_descent(g, g), "+g is not")
      if (allocated(error)) return
      call check(error,.not. scf_gdm_is_descent(g, [2.0_dp, 1.0_dp, 0.0_dp]), "orthogonal is not descent")
      if (allocated(error)) return
      call check(error,.not. scf_gdm_is_descent(g, [nan_value, 0.0_dp, 0.0_dp]), "NaN is not descent")
      if (allocated(error)) return
      ! The margin: g.d = -1e-300 and an orthogonal-to-rounding direction are not descent
      call check(error,.not. scf_gdm_is_descent([1.0e-150_dp, 0.0_dp], [-1.0e-150_dp, 1.0_dp]), &
                 "g.d = -1e-300 is not descent")
      if (allocated(error)) return
      call check(error,.not. scf_gdm_is_descent([1.0_dp, 0.0_dp], [-0.5_dp*scf_gdm_descent_margin, 1.0_dp]), &
                 "inside the relative margin is not descent")
      if (allocated(error)) return
      call check(error, scf_gdm_is_descent([1.0_dp, 0.0_dp], [-2.0_dp*scf_gdm_descent_margin, 1.0_dp]), &
                 "outside the relative margin is descent")
   end subroutine test_descent

end module test_mqc_scf_gdm

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_scf_gdm, only: collect_mqc_scf_gdm_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_scf_gdm", collect_mqc_scf_gdm_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
