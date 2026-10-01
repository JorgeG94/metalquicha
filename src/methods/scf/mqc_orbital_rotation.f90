!! The orbital-rotation step machinery, shared by every second-order optimizer
module mqc_orbital_rotation
   !! `C -> C exp(kappa)`, and the level-shifted trust-region Newton step that
   !! chooses `kappa`.
   !!
   !! Lifted out of `mqc_czt_mcscf`, which had both to itself. Nothing here
   !! knows what the energy is a function of -- not a CI vector, not a density
   !! -- so the same step control serves CASSCF's orbital optimisation and the
   !! second-order SCF, and the constants below are decided once rather than
   !! twice. The CASSCF numbers are pinned by validated tests, so this module
   !! is a *move* and not a rewrite: `rotation_matrix` and `level_shifted_step`
   !! are the bodies that were there, with the flat-parameter mapping left
   !! behind in the caller that owns the pair list.
   !!
   !! ## Why the parametrisation is an exponential
   !!
   !! `exp(kappa)` with `kappa` antisymmetric is orthogonal for any `kappa`, so
   !! a step can be as long as it likes without the orbitals drifting out of
   !! orthonormality. A bad step is then only a bad step, never an invalid
   !! state, which is what lets the trust region be enforced by backtracking on
   !! the energy rather than by projecting the step.
   !!
   !! ## The sign convention, which is not a matter of taste
   !!
   !! `kappa(p, q)` for `p > q` is the parameter, `kappa(q, p) = -kappa(p, q)`,
   !! and the new orbitals are `C exp(kappa)`. Under that convention the
   !! gradient whose components this expects is
   !!
   !!     g_pq = dE / d kappa_pq
   !!
   !! and a descent step is `-H^-1 g`, which is what `level_shifted_step`
   !! returns. Both `C exp(kappa)` and `C exp(-kappa)` appear in the
   !! literature, and taking the wrong one gives an optimiser that climbs;
   !! `test_mqc_mcscf.f90` fixes it by finite differences rather than by
   !! assertion, and `test_mqc_czt_soscf.f90` fixes it again for the SCF.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm, pic_gemv
   use pic_io, only: to_char
   use pic_lapack_interfaces, only: pic_syevd
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: rotation_matrix
   public :: level_shifted_step
   public :: rotation_hessian_t
   public :: subspace_newton_step
   public :: MIN_CURVATURE, SADDLE_CURVATURE, MAX_ROTATION, MIN_ROTATION
   public :: ENERGY_RESOLUTION, TRUST_GROWTH

   real(dp), parameter :: MIN_CURVATURE = 1.0e-3_dp
      !! Smallest curvature the Newton equations are allowed to divide by, in
      !! hartree per radian squared.
      !!
      !! Every eigenvalue is raised until the smallest reaches this, so a mode
      !! that is genuinely stiff keeps its own curvature and only the soft and
      !! the inverted ones are regularised. It has to sit below the softest real
      !! mode and above zero; larger damps directions that did not need it and
      !! costs the quadratic convergence the Hessian was built for.
   real(dp), parameter :: SADDLE_CURVATURE = -1.0e-5_dp
      !! How negative an eigenvalue has to be before the point is called a
      !! saddle rather than a minimum.
      !!
      !! Deliberately loose. A nearly redundant rotation -- an active orbital
      !! whose occupation has gone to two, say -- has a nearly zero eigenvalue
      !! whose sign is decided by how well the CI converged, and chasing those
      !! buys nothing. A real symmetry-breaking mode is orders of magnitude
      !! below this.
   real(dp), parameter :: MAX_ROTATION = 0.2_dp
      !! Largest rotation angle the trust radius is allowed to grow back to, in
      !! radians. Roughly 11 degrees.
   real(dp), parameter :: MIN_ROTATION = 1.0e-6_dp
      !! If backtracking has shrunk the trust radius this far and the energy
      !! still rises, the step direction is not a descent direction and halving
      !! it again will not help. Better to stop and say so.
   real(dp), parameter :: ENERGY_RESOLUTION = 1.0e-12_dp
      !! An energy change this small is not a change, in hartree.
      !!
      !! About fifteen units in the last place of a molecular energy, so it is
      !! the resolution of the arithmetic rather than a convergence criterion.
      !! Near the solution a Newton step is worth `g^2/H`, which drops below
      !! what the energy can report while the gradient is still shrinking, so a
      !! step whose *predicted* gain is below this is taken without being
      !! tested rather than rejected on the noise of the objective.
   real(dp), parameter :: TRUST_GROWTH = 1.3_dp
      !! How fast the trust radius recovers after a successful step. Slower than
      !! it shrinks: an over-long step costs a wasted objective evaluation, an
      !! over-short one only an iteration.

   integer, parameter :: DEFAULT_SUBSPACE = 20
      !! Largest Krylov subspace `subspace_newton_step` builds unless told
   real(dp), parameter :: DEFAULT_SUBSPACE_TOLERANCE = 1.0e-3_dp
      !! Residual of the projected Newton equations, relative to the gradient,
      !! at which `subspace_newton_step` stops expanding unless told
   real(dp), parameter :: PRECONDITION_FLOOR = 1.0e-6_dp
      !! Smallest diagonal element the preconditioner divides by
   real(dp), parameter :: LINEAR_DEPENDENCE = 1.0e-8_dp
      !! Norm below which an orthogonalised direction is taken as already spanned

   type, abstract :: rotation_hessian_t
      !! An orbital Hessian that can be applied to a vector of rotation
      !! parameters, in energy per radian squared, without being written down
   contains
      procedure(apply_rotation_hessian), deferred :: apply
   end type rotation_hessian_t

   abstract interface
      subroutine apply_rotation_hessian(this, x, hx, error)
         !! `H x`
         import :: rotation_hessian_t, dp, error_t
         implicit none
         class(rotation_hessian_t), intent(inout) :: this
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: hx(:)
         type(error_t), intent(inout) :: error
      end subroutine apply_rotation_hessian
   end interface

contains

   subroutine rotation_matrix(kappa, rotation)
      !! `exp(kappa)` for antisymmetric `kappa`, by scaling and squaring
      !!
      !! Orthogonal to machine precision for any step size, so a large step is a
      !! bad step but never an invalid one and no reorthogonalisation is needed
      !! after it. Scaled by halving until the norm is below 1/2, because a
      !! plain Taylor series loses accuracy once the step is not small.
      real(dp), intent(in) :: kappa(:, :)
      real(dp), allocatable, intent(out) :: rotation(:, :)

      real(dp), allocatable :: scaled(:, :), term(:, :), next_term(:, :)
      real(dp) :: norm
      integer :: n, squarings, k, i
      integer, parameter :: TAYLOR_TERMS = 18

      n = size(kappa, 1)
      allocate (scaled(n, n), term(n, n), next_term(n, n), rotation(n, n))

      norm = maxval(abs(kappa))
      squarings = 0
      do while (norm > 0.5_dp)
         norm = 0.5_dp*norm
         squarings = squarings + 1
      end do
      scaled = kappa/real(2**squarings, dp)

      rotation = 0.0_dp
      term = 0.0_dp
      do i = 1, n
         rotation(i, i) = 1.0_dp
         term(i, i) = 1.0_dp
      end do
      do k = 1, TAYLOR_TERMS
         call pic_gemm(term, scaled, next_term)
         term = next_term/real(k, dp)
         rotation = rotation + term
         if (maxval(abs(term)) < 1.0e-18_dp) exit
      end do

      do k = 1, squarings
         call pic_gemm(rotation, rotation, next_term)
         rotation = next_term
      end do

      deallocate (scaled, term, next_term)
   end subroutine rotation_matrix

   subroutine level_shifted_step(hessian, gradient, escape, step, lowest, &
                                 predicted, error)
      !! The Newton step, level shifted, and pushed off a saddle when it is on one
      !!
      !! In the eigenbasis of the Hessian the step is one division per mode, and
      !! two things go wrong there.
      !!
      !! A mode with small or negative curvature would divide by nearly nothing,
      !! so every eigenvalue is raised by a single shift until the smallest
      !! reaches `MIN_CURVATURE`. One shift for all modes rather than a per-mode
      !! floor, because that is the exact solution of the trust-region
      !! subproblem and leaves the well-conditioned modes untouched.
      !!
      !! **A mode with negative curvature and no gradient on it is a saddle**,
      !! and it is where an optimiser built from the gradient alone stops and
      !! reports success. It happens whenever the starting orbitals carry a
      !! symmetry the solution does not. The division gives nothing to work with
      !! there, so such a mode is displaced by `escape` instead; either sign
      !! descends, and the caller backtracks if the step was too long.
      !!
      !! Flat in and flat out: the caller owns whatever mapping there is between
      !! its parameters and orbital pairs, so nothing here needs to know one.
      real(dp), intent(in) :: hessian(:, :)
         !! (n_param, n_param), symmetric, in the same units as `gradient`
      real(dp), intent(in) :: gradient(:)
         !! (n_param), `dE/d kappa` on each parameter
      real(dp), intent(in) :: escape
         !! How far to displace a mode that has negative curvature and no
         !! gradient, in radians
      real(dp), allocatable, intent(out) :: step(:)
         !! (n_param), the displacement to take
      real(dp), intent(out) :: lowest
         !! Smallest Hessian eigenvalue. Negative means the point the step was
         !! taken from is not a minimum, whatever the gradient says.
      real(dp), intent(out) :: predicted
         !! What the quadratic model says the step is worth, as a positive
         !! energy decrease. The caller uses it to tell a step that is too small
         !! to measure from one that is too long to take.
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: vectors(:, :), values(:), projected(:), amplitude(:)
      real(dp) :: shift
      integer :: n_param, k, info

      lowest = 0.0_dp
      predicted = 0.0_dp
      n_param = size(gradient)
      allocate (step(n_param))
      step = 0.0_dp
      if (error%has_error()) return
      if (n_param == 0) return
      if (size(hessian, 1) /= n_param .or. size(hessian, 2) /= n_param) then
         call error%set(ERROR_VALIDATION, "the orbital Hessian is not square in the "// &
                        "number of rotation parameters its gradient carries")
         return
      end if

      allocate (vectors(n_param, n_param), values(n_param))
      allocate (projected(n_param), amplitude(n_param))
      vectors = hessian
      call pic_syevd(vectors, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "the orbital Hessian could not be "// &
                        "diagonalized (info = "//to_char(info)//")")
         return
      end if
      lowest = values(1)

      shift = 0.0_dp
      if (values(1) < MIN_CURVATURE) shift = MIN_CURVATURE - values(1)

      call pic_gemv(vectors, gradient, projected, trans_a="T")

      do k = 1, n_param
         amplitude(k) = -projected(k)/(values(k) + shift)
         if (values(k) < SADDLE_CURVATURE .and. abs(amplitude(k)) < escape) then
            amplitude(k) = sign(escape, amplitude(k))
         end if
      end do

      do k = 1, n_param
         predicted = predicted - projected(k)*amplitude(k) &
                     - 0.5_dp*values(k)*amplitude(k)**2
      end do

      call pic_gemv(vectors, amplitude, step)

      deallocate (vectors, values, projected, amplitude)
   end subroutine level_shifted_step

   subroutine subspace_newton_step(hessian, diagonal, gradient, escape, step, lowest, &
                                   predicted, products, error, max_subspace, tolerance)
      !! The trust-region Newton step, from Hessian-vector products alone
      !!
      !! Builds a Krylov subspace off the preconditioned gradient, projects the
      !! Hessian and the gradient into it, and hands the small dense problem to
      !! `level_shifted_step`, so the level shift, the saddle escape and the
      !! predicted gain are the ones an explicit Hessian gets, applied to a
      !! projected matrix instead of the full one.
      !!
      !! The softest diagonal direction is seeded alongside the gradient. Two
      !! reasons, and the second is the important one: it is where an
      !! instability lives, so the bound on `lowest` is tight; and at a saddle
      !! the gradient is zero and the preconditioned-gradient seed is the zero
      !! vector, leaving nothing to build a subspace from.
      class(rotation_hessian_t), intent(inout) :: hessian
      real(dp), intent(in) :: diagonal(:)
         !! (n_param) the Hessian's diagonal or an approximation to it, for
         !! preconditioning, in the same units
      real(dp), intent(in) :: gradient(:)
         !! (n_param), `dE/d kappa` on each parameter
      real(dp), intent(in) :: escape
         !! How far to displace a mode with negative curvature and no gradient,
         !! in radians
      real(dp), allocatable, intent(out) :: step(:)
         !! (n_param), the displacement to take
      real(dp), intent(out) :: lowest
         !! An **upper bound** on the smallest curvature: the smallest
         !! eigenvalue of the projected Hessian
      real(dp), intent(out) :: predicted
         !! What the quadratic model says the step is worth, as a positive
         !! energy decrease
      integer, intent(out) :: products
         !! Hessian-vector products spent
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_subspace
         !! Default `DEFAULT_SUBSPACE`
      real(dp), intent(in), optional :: tolerance
         !! Default `DEFAULT_SUBSPACE_TOLERANCE`

      real(dp), allocatable :: basis(:, :), image(:, :)
      real(dp), allocatable :: small(:, :), small_gradient(:), amplitude(:)
      real(dp), allocatable :: residual(:), direction(:)
      real(dp) :: norm, overlap, gradient_norm, stop_at
      integer :: n_param, nmax, nsub, i, j, pass

      lowest = 0.0_dp
      predicted = 0.0_dp
      products = 0
      n_param = size(gradient)
      allocate (step(n_param))
      step = 0.0_dp
      if (error%has_error()) return
      if (n_param < 1) return

      nmax = DEFAULT_SUBSPACE
      if (present(max_subspace)) nmax = max_subspace
      nmax = max(1, min(nmax, n_param))
      stop_at = DEFAULT_SUBSPACE_TOLERANCE
      if (present(tolerance)) stop_at = tolerance

      allocate (basis(n_param, nmax), image(n_param, nmax))
      allocate (residual(n_param), direction(n_param))
      gradient_norm = sqrt(dot_product(gradient, gradient))

      ! ---- the seeds ------------------------------------------------------
      ! A negative diagonal element is floored to `PRECONDITION_FLOOR`, which
      ! weights that mode up in the normalised direction rather than flipping
      ! it. Only the subspace depends on this: the step comes from the
      ! projected Hessian, so a soft or negative mode drawn in is what the
      ! saddle escape needs.
      nsub = 0
      direction = -gradient/max(diagonal, PRECONDITION_FLOOR)
      call add_direction(basis, nsub, direction)
      direction = 0.0_dp
      direction(minloc(diagonal, 1)) = 1.0_dp
      call add_direction(basis, nsub, direction)
      if (nsub == 0) then
         ! Unreachable for `n_param >= 1`: the second seed is a unit vector and
         ! there is nothing yet to project it against. Refused rather than
         ! returned, because `lowest = 0` would read as a marginal curvature.
         call error%set(ERROR_VALIDATION, "subspace_newton_step: no direction to "// &
                        "build a subspace from.")
         return
      end if
      do i = 1, nsub
         call hessian%apply(basis(:, i), image(:, i), error)
         products = products + 1
         if (error%has_error()) return
      end do

      ! ---- solve, expand, repeat ------------------------------------------
      do
         if (allocated(small)) deallocate (small, small_gradient)
         allocate (small(nsub, nsub), small_gradient(nsub))
         do j = 1, nsub
            do i = 1, nsub
               small(i, j) = dot_product(basis(:, i), image(:, j))
            end do
            small_gradient(j) = dot_product(basis(:, j), gradient)
         end do
         ! Symmetric to rounding only; the average is what a symmetric solver
         ! would read anyway.
         small = 0.5_dp*(small + transpose(small))

         if (allocated(amplitude)) deallocate (amplitude)
         call level_shifted_step(small, small_gradient, escape, amplitude, lowest, &
                                 predicted, error)
         if (error%has_error()) return

         step = 0.0_dp
         residual = gradient
         do i = 1, nsub
            step = step + amplitude(i)*basis(:, i)
            residual = residual + amplitude(i)*image(:, i)
         end do

         if (nsub >= nmax) exit
         norm = sqrt(dot_product(residual, residual))
         if (norm <= stop_at*max(gradient_norm, tiny(1.0_dp))) exit

         direction = -residual/max(diagonal, PRECONDITION_FLOOR)
         norm = sqrt(dot_product(direction, direction))
         if (norm < tiny(1.0_dp)) exit
         direction = direction/norm
         do pass = 1, 2
            do j = 1, nsub
               overlap = dot_product(basis(:, j), direction)
               direction = direction - overlap*basis(:, j)
            end do
         end do
         norm = sqrt(dot_product(direction, direction))
         if (norm < LINEAR_DEPENDENCE) exit

         nsub = nsub + 1
         basis(:, nsub) = direction/norm
         call hessian%apply(basis(:, nsub), image(:, nsub), error)
         products = products + 1
         if (error%has_error()) return
      end do

      deallocate (basis, image, residual, direction)
      if (allocated(small)) deallocate (small, small_gradient)
      if (allocated(amplitude)) deallocate (amplitude)
   end subroutine subspace_newton_step

   subroutine add_direction(basis, nsub, direction)
      !! Orthonormalise a candidate against the subspace and keep it if anything survives
      real(dp), intent(inout) :: basis(:, :)
      integer, intent(inout) :: nsub
      real(dp), intent(in) :: direction(:)

      real(dp), allocatable :: work(:)
      real(dp) :: norm, overlap
      integer :: j, pass

      work = direction
      norm = sqrt(dot_product(work, work))
      if (norm < tiny(1.0_dp)) return
      work = work/norm
      do pass = 1, 2
         do j = 1, nsub
            overlap = dot_product(basis(:, j), work)
            work = work - overlap*basis(:, j)
         end do
      end do
      norm = sqrt(dot_product(work, work))
      if (norm < LINEAR_DEPENDENCE) return
      if (nsub >= size(basis, 2)) return
      nsub = nsub + 1
      basis(:, nsub) = work/norm
   end subroutine add_direction

end module mqc_orbital_rotation
