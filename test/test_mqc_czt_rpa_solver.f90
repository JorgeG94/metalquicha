!! The paired response solver, against a dense reduction of the same problem
module test_mqc_czt_rpa_solver
   !! `mqc_czt_rpa_solver` knows nothing about excitations, orbitals or
   !! integrals: it takes two symmetric operators as a pair of matrix-vector
   !! products and returns the lowest `w` of `[A B; B A]`. So it is tested the
   !! way it is written -- on a matrix small enough to reduce by hand.
   !!
   !! `(A-B)^{1/2}(A+B)(A-B)^{1/2}` is built explicitly with LAPACK and
   !! diagonalised, which carries no solver tolerance of its own, and the
   !! iterative answer is compared with that. The same fixture then feeds the
   !! three ways the solve can fail to produce an answer -- an indefinite
   !! difference, a tolerance nothing can reach, and a subspace small enough
   !! to have to collapse -- so the error paths are exercised on a problem
   !! whose right answer is known rather than on a molecule.
   !!
   !! **The matrices are deterministic.** A congruential generator, not
   !! `random_number`, because a solver test that fails once in fifty runs on
   !! an unlucky draw is worse than no test: the numbers below are the same on
   !! every machine and every compiler, so a failure is a change in the code.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp, int64
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t
   use mqc_czt_rpa_solver, only: paired_operator_t, rpa_solve, RPA_REASON_NONE, &
                                 RPA_REASON_OTHER, RPA_REASON_UNSTABLE_MINUS, &
                                 RPA_REASON_UNSTABLE_PLUS
   implicit none
   private

   public :: collect_mqc_czt_rpa_solver_tests

   integer, parameter :: N_DIM = 40
      !! Big enough that the subspace has to grow and collapse, small enough
      !! that the explicit reduction is instant.
   integer, parameter :: N_ROOTS = 5

   real(dp), parameter :: TOL_SOLVER = 1.0e-10_dp
      !! The iterative roots against the dense ones. Roots are converged on a
      !! residual of 1e-11 and the eigenvalue error near an eigenvector is
      !! second order in the vector error, so this is a bound the solver
      !! clears by construction rather than one fitted to it.

   real(dp), parameter :: TOL_NORM = 1.0e-12_dp
      !! `|X|^2 - |Y|^2` against one. This is an algebraic identity the
      !! solver imposes by division, not a converged quantity, so it holds to
      !! round-off however far the roots themselves have got.

   type, extends(paired_operator_t) :: dense_paired_t
      !! Two explicit matrices, as the pair of products the solver wants
      real(dp), allocatable :: aplus(:, :), aminus(:, :)
      integer :: n_products = 0
   contains
      procedure :: apply_plus => dense_apply_plus
      procedure :: apply_minus => dense_apply_minus
   end type dense_paired_t

contains

   subroutine collect_mqc_czt_rpa_solver_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("the_paired_solver_matches_the_dense_reduction", &
                               test_dense_agreement), &
                  new_unittest("a_collapsing_subspace_finds_the_same_roots", &
                               test_collapse), &
                  new_unittest("an_indefinite_difference_is_named_an_instability", &
                               test_instability), &
                  new_unittest("an_indefinite_sum_is_named_an_instability", &
                               test_plus_instability), &
                  new_unittest("the_failure_reason_tells_the_two_halves_apart", &
                               test_failure_reason), &
                  new_unittest("an_unreachable_tolerance_stops_and_says_so", &
                               test_stagnation) &
                  ]
   end subroutine collect_mqc_czt_rpa_solver_tests

   subroutine dense_apply_plus(this, vectors, images, error)
      !! `(A+B)v`, as a gemm
      class(dense_paired_t), intent(inout) :: this
      real(dp), intent(in) :: vectors(:, :)
      real(dp), intent(out) :: images(:, :)
      type(error_t), intent(inout) :: error

      if (error%has_error()) return
      call pic_gemm(this%aplus, vectors, images)
      this%n_products = this%n_products + size(vectors, 2)
   end subroutine dense_apply_plus

   subroutine dense_apply_minus(this, vectors, images, error)
      !! `(A-B)v`, as a gemm
      class(dense_paired_t), intent(inout) :: this
      real(dp), intent(in) :: vectors(:, :)
      real(dp), intent(out) :: images(:, :)
      type(error_t), intent(inout) :: error

      if (error%has_error()) return
      call pic_gemm(this%aminus, vectors, images)
      this%n_products = this%n_products + size(vectors, 2)
   end subroutine dense_apply_minus

   function noise(n, seed, scale) result(s)
      !! A deterministic symmetric matrix of the given scale
      !!
      !! A multiplicative congruential generator written out, so the entries
      !! are the same integers on every compiler -- `random_number` is neither
      !! reproducible across implementations nor the same twice without a
      !! seed, and a solver test that samples a different matrix each run
      !! cannot distinguish a regression from an unlucky draw.
      integer, intent(in) :: n, seed
      real(dp), intent(in) :: scale
      real(dp), allocatable :: s(:, :)

      integer, parameter :: MODULUS = 2147483647
      integer, parameter :: MULTIPLIER = 48271
      integer(int64) :: state
      integer :: i, j

      allocate (s(n, n))
      state = int(seed, int64)
      do j = 1, n
         do i = 1, n
            state = mod(state*int(MULTIPLIER, int64), int(MODULUS, int64))
            s(i, j) = scale*(2.0_dp*real(state, dp)/real(MODULUS, dp) - 1.0_dp)
         end do
      end do
      s = 0.5_dp*(s + transpose(s))
   end function noise

   subroutine make_problem(operator, unstable, plus_unstable)
      !! A paired problem whose difference is positive definite, or is not
      !!
      !! `A` is a rising diagonal with a little symmetric noise on it and `B`
      !! is noise alone, which is the shape of a real linear-response problem:
      !! the orbital-energy gaps dominate and the two-electron coupling is a
      !! correction. The diagonal dominance is what makes `A-B` positive
      !! definite, and `unstable` breaks it by pulling one diagonal element
      !! below zero -- which is exactly what an unstable reference does.
      type(dense_paired_t), intent(out) :: operator
      logical, intent(in) :: unstable
      logical, intent(in), optional :: plus_unstable
         !! Break `A+B` instead, leaving `A-B` positive definite. That is the
         !! ordinary real-singlet instability, and it is the half `square_root`
         !! cannot see -- it only ever tests `A-B`. Absent is `.false.`.

      real(dp), allocatable :: a(:, :), b(:, :)
      integer :: i
      logical :: break_plus

      break_plus = .false.
      if (present(plus_unstable)) break_plus = plus_unstable

      a = noise(N_DIM, 20260919, 0.02_dp)
      b = noise(N_DIM, 771131, 0.01_dp)
      do i = 1, N_DIM
         a(i, i) = a(i, i) + 1.0_dp + 0.1_dp*real(i, dp)
      end do
      if (unstable) a(1, 1) = a(1, 1) - 1.6_dp
      ! Pulling `b(1,1)` down pushes `A+B` negative there and `A-B` further
      ! positive, so the two halves of the instability are separated.
      if (break_plus) b(1, 1) = b(1, 1) - 2.6_dp

      operator%aplus = a + b
      operator%aminus = a - b
      deallocate (a, b)
   end subroutine make_problem

   function reduced_spectrum(operator, ok) result(values)
      !! Every `w` of the problem, from the explicit reduction
      !!
      !! `(A-B)^{1/2}(A+B)(A-B)^{1/2}`, diagonalised whole. No iteration and
      !! no tolerance, so what the solver is compared against is the answer
      !! and not a second approximation to it.
      type(dense_paired_t), intent(in) :: operator
      logical, intent(out) :: ok
      real(dp), allocatable :: values(:)

      real(dp), allocatable :: vectors(:, :), w(:), half(:, :), scaled(:, :)
      real(dp), allocatable :: work(:, :), reduced(:, :)
      integer :: info, k

      ok = .false.
      allocate (values(N_DIM))
      vectors = operator%aminus
      allocate (w(N_DIM))
      call pic_syev(vectors, w, jobz="V", uplo="U", info=info)
      if (info /= 0 .or. minval(w) <= 0.0_dp) return

      allocate (scaled(N_DIM, N_DIM), half(N_DIM, N_DIM))
      do k = 1, N_DIM
         scaled(:, k) = vectors(:, k)*sqrt(w(k))
      end do
      call pic_gemm(scaled, vectors, half, transb="T")

      allocate (work(N_DIM, N_DIM), reduced(N_DIM, N_DIM))
      call pic_gemm(operator%aplus, half, work)
      call pic_gemm(half, work, reduced)
      reduced = 0.5_dp*(reduced + transpose(reduced))
      call pic_syev(reduced, values, jobz="N", uplo="U", info=info)
      if (info /= 0) return
      values = sqrt(values)
      ok = .true.
   end function reduced_spectrum

   subroutine test_dense_agreement(error)
      !! The five lowest roots, and the norm the paired problem conserves
      type(error_type), allocatable, intent(out) :: error

      type(dense_paired_t) :: operator
      type(error_t) :: err
      real(dp), allocatable :: reference(:), omega(:), xpy(:, :), xmy(:, :)
      real(dp), allocatable :: residuals(:)
      real(dp) :: worst
      integer :: iterations, products, k
      logical :: ok, converged

      call make_problem(operator, .false.)
      reference = reduced_spectrum(operator, ok)
      call check(error, ok, "the dense reduction of the test problem failed, so the "// &
                 "fixture is wrong rather than the solver")
      if (allocated(error)) return

      call rpa_solve(operator, diagonal_of(operator), N_ROOTS, omega, xpy, xmy, &
                     residuals, iterations, products, converged, err, &
                     tolerance=1.0e-11_dp, max_iterations=100)
      call check(error,.not. err%has_error(), "the paired solve failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, converged, "the paired solve did not converge a problem of "// &
                 "forty dimensions in a hundred iterations")
      if (allocated(error)) return

      worst = maxval(abs(omega - reference(1:N_ROOTS)))
      call check(error, worst < TOL_SOLVER, "a paired root disagrees with the dense "// &
                 "reduction of the same two matrices")
      if (allocated(error)) return

      ! `R . L` is `|X|^2 - |Y|^2`, which is the only norm the paired problem
      ! conserves and the one every amplitude downstream is scaled by.
      do k = 1, N_ROOTS
         call check(error, abs(dot_product(xpy(:, k), xmy(:, k)) - 1.0_dp) < TOL_NORM, &
                    "a paired root came back with |X|^2 - |Y|^2 away from one")
         if (allocated(error)) return
      end do

      ! Every trial vector costs one of each half, and the solver is asked
      ! for both of every vector it adds and no more.
      call check(error, products == operator%n_products, "the solver counted a "// &
                 "different number of products than the operator was asked for")
   end subroutine test_dense_agreement

   subroutine test_collapse(error)
      !! A subspace too small to hold the solve finds the same roots
      !!
      !! The collapse is where the stored products have to survive: every
      !! restart vector is a linear combination of vectors whose images are
      !! already in hand, and a collapse that rebuilt them from scratch -- or
      !! carried the wrong ones -- would show up as a root that moved.
      type(error_type), allocatable, intent(out) :: error

      type(dense_paired_t) :: operator
      type(error_t) :: err
      real(dp), allocatable :: reference(:), omega(:), xpy(:, :), xmy(:, :)
      real(dp), allocatable :: residuals(:)
      integer :: iterations, products
      logical :: ok, converged

      call make_problem(operator, .false.)
      reference = reduced_spectrum(operator, ok)
      call check(error, ok, "the dense reduction of the test problem failed")
      if (allocated(error)) return

      call rpa_solve(operator, diagonal_of(operator), N_ROOTS, omega, xpy, xmy, &
                     residuals, iterations, products, converged, err, &
                     tolerance=1.0e-11_dp, max_iterations=100, max_subspace=30)
      call check(error,.not. err%has_error() .and. converged, &
                 "the paired solve with a capped subspace failed: "//err%get_message())
      if (allocated(error)) return
      call check(error, maxval(abs(omega - reference(1:N_ROOTS))) < TOL_SOLVER, &
                 "a root moved when the subspace was made small enough to collapse")
   end subroutine test_collapse

   subroutine test_instability(error)
      !! An indefinite `(A-B)` is refused by name, not by a square root of a
      !! negative number
      type(error_type), allocatable, intent(out) :: error

      type(dense_paired_t) :: operator
      type(error_t) :: err
      real(dp), allocatable :: omega(:), xpy(:, :), xmy(:, :), residuals(:)
      integer :: iterations, products
      logical :: converged

      call make_problem(operator, .true.)
      call rpa_solve(operator, diagonal_of(operator), N_ROOTS, omega, xpy, xmy, &
                     residuals, iterations, products, converged, err, &
                     tolerance=1.0e-8_dp, max_iterations=50)
      call check(error, err%has_error(), "a paired problem whose difference is "// &
                 "indefinite was solved rather than refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "unstable") > 0, &
                 "an indefinite (A-B) was reported without naming the instability: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error,.not. converged, "an indefinite problem reported convergence")
   end subroutine test_instability

   subroutine test_plus_instability(error)
      !! An indefinite `(A+B)` is refused, not stepped over
      !!
      !! This is the half `square_root` cannot see. `(A-B)` stays positive
      !! definite, so the reduction runs and produces a negative `w^2`; the
      !! solver used to take the lowest `n_roots` of the *positive* squared
      !! frequencies and return them, which is a converged, plausible spectrum
      !! with the lowest state missing and every index below it shifted up.
      !! No error, no warning, on the default method.
      type(error_type), allocatable, intent(out) :: error

      type(dense_paired_t) :: operator
      type(error_t) :: err
      real(dp), allocatable :: omega(:), xpy(:, :), xmy(:, :), residuals(:)
      real(dp), allocatable :: values(:)
      integer :: iterations, products, info
      logical :: converged, ok

      call make_problem(operator, .false., plus_unstable=.true.)

      ! The fixture has to actually have the shape the test is about: `A-B`
      ! positive definite, `A+B` not. Asserted rather than assumed, because a
      ! fixture that breaks both would pass through `square_root` instead and
      ! the test would prove nothing.
      allocate (values(N_DIM))
      call spectrum_of(operator%aminus, values, info)
      call check(error, info == 0 .and. minval(values) > 0.0_dp, &
                 "the fixture was meant to leave (A-B) positive definite")
      if (allocated(error)) return
      call spectrum_of(operator%aplus, values, info)
      call check(error, info == 0 .and. minval(values) < 0.0_dp, &
                 "the fixture was meant to make (A+B) indefinite")
      if (allocated(error)) return

      call rpa_solve(operator, diagonal_of(operator), N_ROOTS, omega, xpy, xmy, &
                     residuals, iterations, products, converged, err, &
                     tolerance=1.0e-8_dp, max_iterations=50)
      call check(error, err%has_error(), "a paired problem whose sum is indefinite "// &
                 "was answered from the roots above the imaginary ones rather than "// &
                 "refused")
      if (allocated(error)) return
      call check(error, index(err%get_message(), "unstable") > 0, &
                 "an indefinite (A+B) was reported without naming the instability: "// &
                 err%get_message())
   end subroutine test_plus_instability

   subroutine test_failure_reason(error)
      !! Why the solve failed comes back as a code, not as prose
      !!
      !! The two instabilities both put the word "unstable" in the message
      !! and they do not mean the same thing. `(A-B)` carries only exchange
      !! and is the same operator in every spin manifold; `(A+B)` carries the
      !! Coulomb term and the kernel, and is where a singlet and a triplet
      !! differ. A caller that matched on the text would call both of them a
      !! triplet instability and tell the user to converge an unrestricted
      !! reference on the strength of the manifold they happened to ask for.
      !!
      !! All four outcomes are checked in one case, because what is being
      !! asserted is that the code *separates* them -- a test of one value
      !! would pass on a routine that returned it unconditionally.
      type(error_type), allocatable, intent(out) :: error

      type(dense_paired_t) :: operator
      type(error_t) :: err
      real(dp), allocatable :: omega(:), xpy(:, :), xmy(:, :), residuals(:)
      integer :: iterations, products, why
      logical :: converged

      call make_problem(operator, .false.)
      call rpa_solve(operator, diagonal_of(operator), N_ROOTS, omega, xpy, xmy, &
                     residuals, iterations, products, converged, err, &
                     tolerance=1.0e-10_dp, max_iterations=100, reason=why)
      call check(error,.not. err%has_error() .and. why == RPA_REASON_NONE, &
                 "a solve that succeeded did not report RPA_REASON_NONE")
      if (allocated(error)) return

      call make_problem(operator, .true.)
      call rpa_solve(operator, diagonal_of(operator), N_ROOTS, omega, xpy, xmy, &
                     residuals, iterations, products, converged, err, &
                     tolerance=1.0e-8_dp, max_iterations=50, reason=why)
      call check(error, why == RPA_REASON_UNSTABLE_MINUS, &
                 "an indefinite (A-B) was not reported as RPA_REASON_UNSTABLE_MINUS")
      if (allocated(error)) return

      call err%clear()
      call make_problem(operator, .false., plus_unstable=.true.)
      call rpa_solve(operator, diagonal_of(operator), N_ROOTS, omega, xpy, xmy, &
                     residuals, iterations, products, converged, err, &
                     tolerance=1.0e-8_dp, max_iterations=50, reason=why)
      call check(error, why == RPA_REASON_UNSTABLE_PLUS, &
                 "an indefinite (A+B) was not reported as RPA_REASON_UNSTABLE_PLUS")
      if (allocated(error)) return

      ! A stable problem asked for more roots than it has dimensions: a
      ! failure that says nothing at all about the reference, and the one a
      ! text match on "unstable" would also get right. It is here for the
      ! other direction -- that an ordinary failure is not promoted.
      call err%clear()
      call make_problem(operator, .false.)
      call rpa_solve(operator, diagonal_of(operator), N_DIM + 1, omega, xpy, xmy, &
                     residuals, iterations, products, converged, err, reason=why)
      call check(error, err%has_error() .and. why == RPA_REASON_OTHER, &
                 "an ordinary failure was not reported as RPA_REASON_OTHER")
   end subroutine test_failure_reason

   subroutine spectrum_of(matrix, values, info)
      !! Eigenvalues of a symmetric matrix, for a fixture assertion
      real(dp), intent(in) :: matrix(:, :)
      real(dp), intent(out) :: values(:)
      integer, intent(out) :: info

      real(dp), allocatable :: vectors(:, :)

      vectors = matrix
      call pic_syev(vectors, values, jobz="N", uplo="U", info=info)
   end subroutine spectrum_of

   subroutine test_stagnation(error)
      !! A tolerance nothing can reach returns, with the residual it got to
      !!
      !! Below round-off there is no correction left to add: every residual
      !! direction is already in the subspace, which by then is the whole
      !! space. The solver has to say so rather than spend its iteration
      !! budget adding vectors that project to nothing.
      type(error_type), allocatable, intent(out) :: error

      type(dense_paired_t) :: operator
      type(error_t) :: err
      real(dp), allocatable :: omega(:), xpy(:, :), xmy(:, :), residuals(:)
      integer :: iterations, products
      logical :: converged

      call make_problem(operator, .false.)
      call rpa_solve(operator, diagonal_of(operator), N_ROOTS, omega, xpy, xmy, &
                     residuals, iterations, products, converged, err, &
                     tolerance=1.0e-30_dp, max_iterations=200)
      call check(error, err%has_error(), "a tolerance below round-off was reported "// &
                 "as met")
      if (allocated(error)) return
      call check(error,.not. converged, "a solve that could not reach its tolerance "// &
                 "reported convergence")
      if (allocated(error)) return
      call check(error, iterations < 200, "the solver spent its whole iteration "// &
                 "budget instead of noticing it had stopped improving")
      if (allocated(error)) return
      ! The point of the named error is that it carries the residual reached,
      ! so a caller can decide whether the answer is good enough anyway.
      call check(error, index(err%get_message(), "residual") > 0, &
                 "a stalled solve did not report the residual it reached: "// &
                 err%get_message())
   end subroutine test_stagnation

   function diagonal_of(operator) result(diagonal)
      !! What the solver preconditions on: the diagonal of `(A-B)`
      type(dense_paired_t), intent(in) :: operator
      real(dp), allocatable :: diagonal(:)

      integer :: i

      allocate (diagonal(N_DIM))
      do i = 1, N_DIM
         diagonal(i) = operator%aminus(i, i)
      end do
   end function diagonal_of

end module test_mqc_czt_rpa_solver

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_rpa_solver, only: collect_mqc_czt_rpa_solver_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_rpa_solver", collect_mqc_czt_rpa_solver_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
