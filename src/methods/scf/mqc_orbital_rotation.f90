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
   use pic_blas_interfaces, only: pic_gemm
   use pic_io, only: to_char
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: rotation_matrix
   public :: level_shifted_step
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
      integer :: n_param, k, l, info

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
      call pic_syev(vectors, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "the orbital Hessian could not be "// &
                        "diagonalized (info = "//to_char(info)//")")
         return
      end if
      lowest = values(1)

      shift = 0.0_dp
      if (values(1) < MIN_CURVATURE) shift = MIN_CURVATURE - values(1)

      do k = 1, n_param
         projected(k) = 0.0_dp
         do l = 1, n_param
            projected(k) = projected(k) + vectors(l, k)*gradient(l)
         end do
      end do

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

      do l = 1, n_param
         step(l) = dot_product(vectors(l, :), amplitude)
      end do

      deallocate (vectors, values, projected, amplitude)
   end subroutine level_shifted_step

end module mqc_orbital_rotation
