!! Solving the coupled-perturbed equations, whatever the perturbation is
module mqc_czt_response
   !! One iterative solver, and a type per kind of perturbation.
   !!
   !! **The two perturbations do not have the same shape.** An electric field
   !! leaves the basis functions where they are, so the overlap is unchanged,
   !! the orbitals stay orthonormal for free, and the response is a
   !! virtual-by-occupied rectangle. Displacing a nucleus drags its functions
   !! with it: orthonormality has to be maintained while that happens, which
   !! fixes the occupied-occupied block of the response at minus half the
   !! overlap derivative, and the vector stops being a rectangle.
   !!
   !! So the solver is written against a flat vector and an abstract operator,
   !! and the two cases differ only in what they put behind it -- the same
   !! arrangement `mqc_davidson` uses. Dispatch costs nothing: `apply` is called
   !! once per iteration and each call does a Fock build over the whole basis.
   use pic_types, only: dp, default_int, int64
   use mqc_error, only: error_t, ERROR_VALIDATION
   use pic_lapack_interfaces, only: pic_getrf, pic_getrs
   use pic_blas_interfaces, only: pic_gemm, pic_gemv
   use pic_logger, only: logger => global_logger
   use mqc_program_limits, only: MAX_LINE_LENGTH
   use mqc_calculation_defaults, only: DEFAULT_RESPONSE_TOL, DEFAULT_RESPONSE_MAX_ITER
   implicit none
   private

   public :: response_operator_t
   public :: solve_response

   type, abstract :: response_operator_t
      !! Whatever can multiply a trial response by the electronic Hessian
      !!
      !! The solver asks for one thing: given a trial vector, the two-electron
      !! response it induces, in the same flat layout. What that layout means --
      !! which rows are virtual, whether there are occupied rows at all -- is
      !! the operator's business and not the solver's.
   contains
      procedure(apply_response), deferred :: apply
      procedure(response_length), deferred :: length
      procedure :: note_header => response_note_header
      procedure :: note => response_note
      procedure :: denominators => response_denominators
   end type response_operator_t

   abstract interface
      subroutine apply_response(this, vector, image, error)
         !! The induced response to a trial vector, in the same layout
         import :: response_operator_t, dp, error_t
         implicit none
         class(response_operator_t), intent(inout) :: this
         real(dp), intent(in) :: vector(:)
         real(dp), intent(out) :: image(:)
         type(error_t), intent(inout) :: error
      end subroutine apply_response

      pure function response_length(this) result(n)
         !! How long a trial vector is, which only the operator knows
         import :: response_operator_t
         implicit none
         class(response_operator_t), intent(in) :: this
         integer :: n
      end function response_length
   end interface

contains

   function response_note_header(this) result(text)
      !! Column headings for what `note` reports, empty unless overridden
      class(response_operator_t), intent(in) :: this
      character(len=:), allocatable :: text

      text = ""
   end function response_note_header

   function response_note(this) result(text)
      !! What the last `apply` did, for the solver's per-cycle line
      !!
      !! Empty by default. An operator that counts something worth seeing --
      !! the nuclear one counts integrals -- overrides both this and
      !! `note_header`, and the solver prints whatever comes back.
      class(response_operator_t), intent(in) :: this
      character(len=:), allocatable :: text

      text = ""
   end function response_note

   function response_denominators(this) result(den)
      !! The diagonal `Delta` that `apply` divides by, for an operator whose
      !! `apply` is `Delta^-1 K` with `K` symmetric on the elements where
      !! `Delta` is positive. Zero marks an element `apply` reads but does
      !! not produce, which the fixed point leaves at its right-hand side.
      !! Zero-size when the operator has none, which rules the
      !! conjugate-gradient solve out.
      class(response_operator_t), intent(in) :: this
      real(dp), allocatable :: den(:)

      allocate (den(0))
   end function response_denominators

   subroutine solve_response(operator, rhs, seed, solution, error, max_iter, tol, &
                             iterations, label, columns, conjugate_gradient)
      !! Iterate `x = seed + K(x)` until it stops moving
      !!
      !! The coupled-perturbed equations put the two-electron response of the
      !! solution back into its own right-hand side, which is what makes them
      !! coupled and what makes this a fixed point rather than a linear solve
      !! against a matrix anyone builds. `seed` is the uncoupled answer -- the
      !! right-hand side divided by the orbital energy differences -- and each
      !! pass adds the response the current estimate induces.
      !!
      !! **A Krylov subspace method, because a pass is what costs.** Every pass
      !! is a Fock build over the whole basis, and a subspace method needs far
      !! fewer of them than straight geometric iteration.
      !!
      !! The equation is `(1 - K) x = rhs`. Each cycle orthogonalises the current
      !! residual against the basis built so far, spends its one operator
      !! application on that new direction, and then solves the projected problem
      !! in the whole subspace exactly. `K x` is the same combination of the
      !! images the basis vectors already produced, so the residual is assembled
      !! from vectors in hand rather than from a second pass.
      !!
      !! **Block form when the vector stacks several perturbations.** With
      !! `columns` greater than one, each column is its own right-hand side of
      !! the same operator, each contributes a direction per cycle to one
      !! shared subspace, and each picks its own coefficients in it. A cycle is
      !! still one operator application -- the new directions travel together
      !! -- so the pass count can only fall. With one column it is the plain
      !! method above, direction for direction.
      !!
      !! `basis` and `images` are each `n/columns` by `columns*max_iter`, so
      !! the memory is two vectors per direction *allowed* rather than per
      !! direction used. Bound it by restarting the subspace, not by shortening
      !! it: a shorter one converges more slowly and each cycle is a Fock build.
      !!
      !! Near a triplet instability `1 - K` stops being positive definite and the
      !! projected system goes singular, which is reported rather than iterated
      !! through. A subspace method does not rescue a reference that is not a
      !! minimum.
      class(response_operator_t), intent(inout) :: operator
      real(dp), intent(in) :: rhs(:)
         !! The uncoupled solution, already scaled by the energy denominators.
         !! The two operators build it differently.
      real(dp), intent(in) :: seed(:)
         !! Where to start. Usually `rhs` itself; a previous solution when a
         !! caller is walking a sequence of related perturbations.
      real(dp), allocatable, intent(out) :: solution(:)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: iterations
      character(len=*), intent(in), optional :: label
         !! Present, one line per cycle is logged at `performance` level under
         !! this heading: cycle, subspace size, columns still improving,
         !! largest residual, and whatever the operator's `note` adds. Absent,
         !! the solve is silent.
      integer, intent(in), optional :: columns
         !! How many independent perturbations the flat vector stacks, each
         !! `n/columns` long, column fastest. One, the default, treats the
         !! vector as a single system.

      logical, intent(in), optional :: conjugate_gradient
         !! Solve by preconditioned conjugate gradient instead of the shared
         !! Krylov subspace. Needs an operator with `denominators`; the seed
         !! is then ignored. The search direction is applied at its own size,
         !! so an absolute density screen inside the operator skips more with
         !! every iteration as the residual shrinks, and nothing is divided
         !! by that size afterwards to amplify what the screen dropped.

      real(dp), allocatable :: basis(:, :), images(:, :)
      real(dp), allocatable :: rhs_c(:, :), sol(:, :), resid(:, :), direction(:, :)
      real(dp), allocatable :: probe(:), image(:), probe_c(:, :)
      real(dp), allocatable :: smat(:, :), factored(:, :), coeff(:, :)
      real(dp), allocatable :: col_change(:), overlaps(:)
      integer(default_int), allocatable :: pivots(:)
      integer(default_int) :: info
      integer :: n, m, ncol, cycles, it, k, kmax, i, c, pass, added, active, kk
      real(dp) :: threshold, change, norm
      integer(int64) :: clock_last, clock_now, clock_rate
      character(len=MAX_LINE_LENGTH) :: line

      if (error%has_error()) return

      if (present(conjugate_gradient)) then
         if (conjugate_gradient) then
            call solve_response_cg(operator, rhs, solution, error, max_iter, tol, &
                                   iterations, label, columns)
            return
         end if
      end if

      n = operator%length()
      if (size(rhs) /= n .or. size(seed) /= n) then
         call error%set(ERROR_VALIDATION, "the coupled-perturbed right-hand side is "// &
                        "not the length this operator works on, so one of them "// &
                        "describes a different orbital space.")
         return
      end if
      ncol = 1
      if (present(columns)) ncol = columns
      if (ncol < 1 .or. mod(n, ncol) /= 0) then
         call error%set(ERROR_VALIDATION, "the coupled-perturbed vector does not divide "// &
                        "into the number of columns the caller says it stacks.")
         return
      end if
      m = n/ncol

      cycles = DEFAULT_RESPONSE_MAX_ITER
      if (present(max_iter)) cycles = max_iter
      threshold = DEFAULT_RESPONSE_TOL
      if (present(tol)) threshold = tol

      kmax = cycles*ncol
      allocate (solution(n), probe(n), image(n))
      allocate (basis(m, kmax), images(m, kmax), smat(kmax, kmax))
      allocate (rhs_c(m, ncol), sol(m, ncol), resid(m, ncol), direction(m, ncol))
      allocate (probe_c(m, ncol), col_change(ncol), overlaps(kmax))

      rhs_c = reshape(rhs, [m, ncol])
      ! The first directions are the seed -- the uncoupled answer, which is
      ! already most of the coupled one, so the subspace starts pointed the
      ! right way rather than at the right-hand side alone.
      direction = reshape(seed, [m, ncol])
      sol = direction
      col_change = huge(1.0_dp)
      k = 0

      if (present(label)) then
         write (line, "(a,es9.1)") "    "//trim(label)//", residual threshold ", threshold
         call logger%performance(trim(line))
         call logger%performance("      cycle  subspace  improving  max residual   seconds" &
                                 //operator%note_header())
      end if
      call system_clock(clock_last, clock_rate)

      do it = 1, cycles
         if (present(iterations)) iterations = it

         ! A converged column contributes no direction: its residual is below
         ! threshold and would only add noise to the subspace.
         !
         ! **Orthogonalised twice.** Once is not enough in floating point when
         ! the new direction is nearly dependent on the subspace, which is
         ! precisely what happens as this converges: the residual shrinks
         ! toward the span and a single pass leaves a component of rounding
         ! error that then gets normalised up to unit length. The directions
         ! accepted earlier in this cycle are part of the span too.
         added = 0
         active = 0
         do c = 1, ncol
            if (col_change(c) < threshold) cycle
            active = active + 1
            ! Classical Gram-Schmidt as two matrix-vector products rather
            ! than one dot product per basis vector: the projection onto the
            ! whole span, then its removal. Twice, which is what makes the
            ! classical form as good as the modified one.
            kk = k + added
            do pass = 1, 2
               if (kk == 0) exit
               call pic_gemv(basis(:, 1:kk), direction(:, c), overlaps(1:kk), trans_a="T")
               call pic_gemv(basis(:, 1:kk), overlaps(1:kk), direction(:, c), &
                             alpha=-1.0_dp, beta=1.0_dp)
            end do
            norm = sqrt(dot_product(direction(:, c), direction(:, c)))
            ! Nothing left that the subspace does not already span.
            if (norm < 1.0e-14_dp) cycle
            added = added + 1
            basis(:, k + added) = direction(:, c)/norm
         end do
         ! No column has anything new: the current solution is the answer,
         ! not a failure to find one.
         if (added == 0) then
            solution = reshape(sol, [n])
            return
         end if

         ! One pass for every new direction. Columns with nothing new carry a
         ! zero, whose image is zero and costs the contraction and not the
         ! integrals.
         probe_c = 0.0_dp
         probe_c(:, 1:added) = basis(:, k + 1:k + added)
         probe = reshape(probe_c, [n])
         call operator%apply(probe, image, error)
         if (error%has_error()) return
         probe_c = reshape(image, [m, ncol])
         images(:, k + 1:k + added) = probe_c(:, 1:added)

         ! `smat` holds `-<v_i|K|v_j>` and the identity is added when the
         ! system is factored, so this stays a plain accumulation of what each
         ! new vector contributes -- a row and a column per new direction.
         ! The new columns and rows of the projected operator, `-B^T K B`
         ! restricted to what this cycle added, as two matrix products.
         kk = k + added
         call pic_gemm(basis(:, 1:kk), images(:, k + 1:kk), smat(1:kk, k + 1:kk), &
                       transa="T", alpha=-1.0_dp, beta=0.0_dp)
         call pic_gemm(basis(:, k + 1:kk), images(:, 1:kk), smat(k + 1:kk, 1:kk), &
                       transa="T", alpha=-1.0_dp, beta=0.0_dp)
         k = kk

         allocate (factored(k, k), pivots(k), coeff(k, ncol))
         factored = smat(1:k, 1:k)
         do i = 1, k
            factored(i, i) = factored(i, i) + 1.0_dp
         end do
         call pic_gemm(basis(:, 1:k), rhs_c, coeff, transa="T")
         call pic_getrf(factored, pivots, info)
         if (info /= 0) then
            deallocate (factored, pivots, coeff)
            call error%set(ERROR_VALIDATION, "the projected coupled-perturbed system "// &
                           "is singular. That is what an instability of the reference "// &
                           "looks like from inside the solver: 1 - K has stopped being "// &
                           "positive definite, so a reference that is not a minimum is "// &
                           "the first thing to check.")
            return
         end if
         call pic_getrs(factored, pivots, coeff, info=info)

         ! The subspace solution of every column, and the residual it leaves.
         ! `K x` is the same combination of the images already computed, so
         ! this costs no pass.
         call pic_gemm(basis(:, 1:k), coeff, sol)
         resid = rhs_c
         call pic_gemm(images(:, 1:k), coeff, resid, beta=1.0_dp)
         resid = resid - sol
         deallocate (factored, pivots, coeff)

         do c = 1, ncol
            col_change(c) = maxval(abs(resid(:, c)))
         end do
         change = maxval(col_change)
         if (present(label)) then
            call system_clock(clock_now)
            write (line, "(a,i7,i10,i11,es14.3,f10.2,a)") "    ", it, k, active, change, &
               real(clock_now - clock_last, dp)/real(clock_rate, dp), operator%note()
            clock_last = clock_now
            call logger%performance(trim(line))
         end if
         if (change < threshold) then
            solution = reshape(sol, [n])
            return
         end if
         direction = resid
      end do

      call error%set(ERROR_VALIDATION, "the coupled-perturbed equations did not "// &
                     "converge. The response is solved in a Krylov subspace, which "// &
                     "stops being well conditioned near an instability of the "// &
                     "reference, so a reference that is not a minimum is the first "// &
                     "thing to check.")
   end subroutine solve_response

   subroutine solve_response_cg(operator, rhs, solution, error, max_iter, tol, &
                                iterations, label, columns)
      !! Preconditioned conjugate gradient on `(Delta - K) x = Delta rhs`
      !!
      !! The fixed point `x = rhs + apply(x)` with `apply = -Delta^-1 G`
      !! is `(Delta + G) x = Delta rhs`, symmetric positive definite while
      !! the reference is a minimum, so conjugate gradient with `Delta` as
      !! the preconditioner applies. `Delta p + G p` is
      !! `Delta * (p - apply(p))`, one operator application per iteration
      !! for every column at once.
      !!
      !! **The direction is applied at its own size.** It is a combination of
      !! preconditioned residuals and shrinks as the solve converges, so the
      !! absolute density screen inside the operator skips more each
      !! iteration; and since the update uses `A p` directly rather than
      !! through a normalised basis, what the screen drops is an absolute
      !! error of the screen's own size, never divided by a small norm.
      !! Three vectors per column, no growing subspace, nothing to
      !! orthogonalise.
      !!
      !! Elements with a zero denominator are ones `apply` reads but never
      !! produces -- the occupied-occupied block a nuclear perturbation fixes
      !! through the overlap constraint. The fixed point leaves them at the
      !! right-hand side, so they are set there once, their mean field is
      !! applied once and moved to the right-hand side, and the iteration
      !! runs on the rest, where the operator is symmetric.
      !!
      !! Convergence is on the same quantity as the subspace solver, the
      !! largest element of `Delta^-1 r`, which is the residual of the fixed
      !! point.
      class(response_operator_t), intent(inout) :: operator
      real(dp), intent(in) :: rhs(:)
      real(dp), allocatable, intent(out) :: solution(:)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: iterations
      character(len=*), intent(in), optional :: label
      integer, intent(in), optional :: columns

      real(dp), allocatable :: den(:, :), x(:, :), r(:, :), z(:, :), p(:, :), q(:, :)
      real(dp), allocatable :: probe(:), image(:), rz(:), col_change(:), b(:, :)
      logical, allocatable :: free(:, :)
      integer :: n, m, ncol, cycles, it, c, active
      real(dp) :: threshold, change, pq, alpha, beta, rz_new
      integer(int64) :: clock_last, clock_now, clock_rate
      character(len=MAX_LINE_LENGTH) :: line

      if (error%has_error()) return

      n = operator%length()
      if (size(rhs) /= n) then
         call error%set(ERROR_VALIDATION, "the coupled-perturbed right-hand side is "// &
                        "not the length this operator works on, so one of them "// &
                        "describes a different orbital space.")
         return
      end if
      ncol = 1
      if (present(columns)) ncol = columns
      if (ncol < 1 .or. mod(n, ncol) /= 0) then
         call error%set(ERROR_VALIDATION, "the coupled-perturbed vector does not divide "// &
                        "into the number of columns the caller says it stacks.")
         return
      end if
      m = n/ncol

      allocate (probe(n), image(n))
      probe = operator%denominators()
      if (size(probe) /= n) then
         call error%set(ERROR_VALIDATION, "the conjugate-gradient solve was asked of a "// &
                        "response operator that carries no preconditioner diagonal, "// &
                        "so its symmetric form is not available.")
         return
      end if
      allocate (den(m, ncol))
      den = reshape(probe, [m, ncol])

      cycles = DEFAULT_RESPONSE_MAX_ITER
      if (present(max_iter)) cycles = max_iter
      threshold = DEFAULT_RESPONSE_TOL
      if (present(tol)) threshold = tol

      allocate (solution(n), x(m, ncol), r(m, ncol), z(m, ncol), p(m, ncol), q(m, ncol))
      allocate (rz(ncol), col_change(ncol), b(m, ncol), free(m, ncol))
      free = den > 0.0_dp

      if (present(label)) then
         write (line, "(a,es9.1)") "    "//trim(label)//", residual threshold ", threshold
         call logger%performance(trim(line))
         call logger%performance("      cycle  subspace  improving  max residual   seconds" &
                                 //operator%note_header())
      end if
      call system_clock(clock_last, clock_rate)

      ! The fixed elements take their right-hand side; their mean field,
      ! applied once, joins the right-hand side of the free ones. A full
      ! pass, and reported as one: it is the difference between the passes
      ! the operator counts and the cycles below.
      b = reshape(rhs, [m, ncol])
      x = 0.0_dp
      where (.not. free) x = b
      probe = reshape(x, [n])
      call operator%apply(probe, image, error)
      if (error%has_error()) return
      q = reshape(image, [m, ncol])
      where (free) b = b + q
      if (present(label)) then
         call system_clock(clock_now)
         write (line, "(a,f10.2,a)") "      fixed occupied block, one pass         ", &
            real(clock_now - clock_last, dp)/real(clock_rate, dp), operator%note()
         call logger%performance(trim(line))
         clock_last = clock_now
      end if

      ! From zero on the free elements: the first direction is then the
      ! uncoupled response, which is what the subspace solver starts from.
      z = 0.0_dp
      where (free) z = b
      r = den*z
      p = z
      do c = 1, ncol
         rz(c) = dot_product(r(:, c), z(:, c))
         col_change(c) = maxval(abs(z(:, c)))
         if (col_change(c) < threshold) p(:, c) = 0.0_dp
      end do

      do it = 1, cycles
         if (present(iterations)) iterations = it
         if (maxval(col_change) < threshold) then
            solution = reshape(x, [n])
            return
         end if

         probe = reshape(p, [n])
         call operator%apply(probe, image, error)
         if (error%has_error()) return
         q = reshape(image, [m, ncol])
         q = den*(p - q)
         where (.not. free) q = 0.0_dp

         active = 0
         do c = 1, ncol
            if (col_change(c) < threshold) cycle
            active = active + 1
            pq = dot_product(p(:, c), q(:, c))
            if (pq <= 0.0_dp) then
               call error%set(ERROR_VALIDATION, "the coupled-perturbed operator is not "// &
                              "positive definite along a conjugate-gradient direction. "// &
                              "That is what an instability of the reference looks like "// &
                              "from inside the solver, so a reference that is not a "// &
                              "minimum is the first thing to check.")
               return
            end if
            alpha = rz(c)/pq
            x(:, c) = x(:, c) + alpha*p(:, c)
            r(:, c) = r(:, c) - alpha*q(:, c)
            where (free(:, c)) z(:, c) = r(:, c)/den(:, c)
            col_change(c) = maxval(abs(z(:, c)))
            if (col_change(c) < threshold) then
               p(:, c) = 0.0_dp
               cycle
            end if
            rz_new = dot_product(r(:, c), z(:, c))
            beta = rz_new/rz(c)
            rz(c) = rz_new
            p(:, c) = z(:, c) + beta*p(:, c)
         end do
         change = maxval(col_change)
         if (present(label)) then
            call system_clock(clock_now)
            write (line, "(a,i7,i10,i11,es14.3,f10.2,a)") "    ", it, it, active, change, &
               real(clock_now - clock_last, dp)/real(clock_rate, dp), operator%note()
            clock_last = clock_now
            call logger%performance(trim(line))
         end if
      end do
      if (maxval(col_change) < threshold) then
         solution = reshape(x, [n])
         return
      end if

      call error%set(ERROR_VALIDATION, "the coupled-perturbed equations did not "// &
                     "converge by conjugate gradient. That operator is positive definite "// &
                     "only while the reference is a minimum, so a reference that is "// &
                     "not one is the first thing to check.")
   end subroutine solve_response_cg

end module mqc_czt_response
