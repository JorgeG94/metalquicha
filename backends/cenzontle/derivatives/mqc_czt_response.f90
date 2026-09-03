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
   use pic_types, only: dp, default_int
   use mqc_error, only: error_t, ERROR_VALIDATION
   use pic_lapack_interfaces, only: pic_getrf, pic_getrs
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

   subroutine solve_response(operator, rhs, seed, solution, error, max_iter, tol, &
                             iterations)
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
      !! `basis` and `images` are each `n` by `max_iter`, so the memory is two
      !! vectors per cycle *allowed* rather than per cycle used. Bound it by
      !! restarting the subspace, not by shortening it: a shorter one converges
      !! more slowly and each cycle is a Fock build.
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

      real(dp), allocatable :: basis(:, :), images(:, :)
      real(dp), allocatable :: direction(:), image(:), residual(:)
      real(dp), allocatable :: smat(:, :), factored(:, :), coeff(:, :)
      integer(default_int), allocatable :: pivots(:)
      integer(default_int) :: info
      integer :: n, cycles, it, k, m, pass
      real(dp) :: threshold, change, norm, overlap

      if (error%has_error()) return

      n = operator%length()
      if (size(rhs) /= n .or. size(seed) /= n) then
         call error%set(ERROR_VALIDATION, "the coupled-perturbed right-hand side is "// &
                        "not the length this operator works on, so one of them "// &
                        "describes a different orbital space.")
         return
      end if

      cycles = 50
      if (present(max_iter)) cycles = max_iter
      threshold = 1.0e-9_dp
      if (present(tol)) threshold = tol

      allocate (solution(n), direction(n), image(n), residual(n))
      allocate (basis(n, cycles), images(n, cycles))
      allocate (smat(cycles, cycles), coeff(cycles, 1))

      ! The first direction is the seed -- the uncoupled answer, which is
      ! already most of the coupled one, so the subspace starts pointed the
      ! right way rather than at the right-hand side alone.
      solution = seed
      direction = seed
      k = 0

      do it = 1, cycles
         if (present(iterations)) iterations = it

         ! **Orthogonalised twice.** Once is not enough in floating point when
         ! the new direction is nearly dependent on the subspace, which is
         ! precisely what happens as this converges: the residual shrinks
         ! toward the span and a single pass leaves a component of rounding
         ! error that then gets normalised up to unit length.
         do pass = 1, 2
            do m = 1, k
               overlap = dot_product(basis(:, m), direction)
               direction = direction - overlap*basis(:, m)
            end do
         end do

         norm = sqrt(dot_product(direction, direction))
         ! Nothing left that the subspace does not already span: the current
         ! solution is the answer, not a failure to find one.
         if (norm < 1.0e-14_dp) return

         k = k + 1
         basis(:, k) = direction/norm
         call operator%apply(basis(:, k), image, error)
         if (error%has_error()) return
         images(:, k) = image

         ! `smat` holds `-<v_i|K|v_j>` and the identity is added when the
         ! system is factored, so this stays a plain accumulation of what each
         ! new vector contributes -- one row and one column per cycle.
         do m = 1, k
            smat(m, k) = -dot_product(basis(:, m), images(:, k))
            smat(k, m) = -dot_product(basis(:, k), images(:, m))
         end do
         do m = 1, k
            coeff(m, 1) = dot_product(basis(:, m), rhs)
         end do

         allocate (factored(k, k), pivots(k))
         factored = smat(1:k, 1:k)
         do m = 1, k
            factored(m, m) = factored(m, m) + 1.0_dp
         end do
         call pic_getrf(factored, pivots, info)
         if (info /= 0) then
            deallocate (factored, pivots)
            call error%set(ERROR_VALIDATION, "the projected coupled-perturbed system "// &
                           "is singular. That is what an instability of the reference "// &
                           "looks like from inside the solver: 1 - K has stopped being "// &
                           "positive definite, so a reference that is not a minimum is "// &
                           "the first thing to check.")
            return
         end if
         call pic_getrs(factored, pivots, coeff(1:k, :), info=info)
         deallocate (factored, pivots)

         ! The subspace solution, and the residual it leaves. `K x` is the same
         ! combination of the images already computed, so this costs no pass.
         solution = 0.0_dp
         residual = rhs
         do m = 1, k
            solution = solution + coeff(m, 1)*basis(:, m)
            residual = residual + coeff(m, 1)*images(:, m)
         end do
         residual = residual - solution

         change = maxval(abs(residual))
         if (change < threshold) return
         direction = residual
      end do

      call error%set(ERROR_VALIDATION, "the coupled-perturbed equations did not "// &
                     "converge. The response is solved in a Krylov subspace, which "// &
                     "stops being well conditioned near an instability of the "// &
                     "reference, so a reference that is not a minimum is the first "// &
                     "thing to check.")
   end subroutine solve_response

end module mqc_czt_response
