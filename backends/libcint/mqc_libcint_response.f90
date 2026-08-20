!! Solving the coupled-perturbed equations, whatever the perturbation is
module mqc_libcint_response
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
   !! and the two cases differ only in what they put behind it. That is the
   !! arrangement `mqc_davidson` already uses for the same reason -- a
   !! restricted active space cannot shape its vector either -- and following it
   !! makes this read as the house pattern rather than as an invention.
   !!
   !! Dispatch costs nothing here. `apply` is called once per iteration, a few
   !! dozen times in total, and each call does a Fock build over the whole
   !! basis. One indirect call against that is not measurable.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
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
      !! **Straight iteration, not a Krylov method.** The operator is a
      !! contraction for the systems this is used on, so this converges
      !! geometrically at a rate set by how far the reference is from an
      !! instability; the cost is one Fock build per pass either way, and a
      !! subspace method would save passes rather than making them cheaper.
      !! Where it will not do is near a triplet instability, and the failure is
      !! loud rather than silent: the residual stops falling and the iteration
      !! count runs out.
      class(response_operator_t), intent(inout) :: operator
      real(dp), intent(in) :: rhs(:)
         !! The uncoupled solution, already scaled by the energy denominators.
         !! Named for what it is to this routine rather than for where it came
         !! from, since the two operators build it differently.
      real(dp), intent(in) :: seed(:)
         !! Where to start. Usually `rhs` itself; a previous solution when a
         !! caller is walking a sequence of related perturbations.
      real(dp), allocatable, intent(out) :: solution(:)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: iterations

      real(dp), allocatable :: image(:), previous(:)
      integer :: n, cycles, it
      real(dp) :: threshold, change

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

      allocate (solution(n), image(n), previous(n))
      solution = seed

      do it = 1, cycles
         previous = solution
         call operator%apply(solution, image, error)
         if (error%has_error()) return
         solution = rhs + image
         change = maxval(abs(solution - previous))
         if (present(iterations)) iterations = it
         if (change < threshold) return
      end do

      call error%set(ERROR_VALIDATION, "the coupled-perturbed equations did not "// &
                     "converge. The response is iterated as a fixed point, which "// &
                     "stops contracting near an instability of the reference, so a "// &
                     "reference that is not a minimum is the first thing to check.")
   end subroutine solve_response

end module mqc_libcint_response
