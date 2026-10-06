!! The second-order step a closed-shell SCF takes once DIIS has got it close
module mqc_czt_soscf
   !! Trust-region Newton in the orbital-rotation space, matrix free.
   !!
   !! The pieces of a second-order SCF that are *not* here: the Fock build, the
   !! energy and the convergence rule all belong to `mqc_czt_rhf`, which owns
   !! the operators and the iteration table, and the level-shifted step control
   !! and the matrix exponential belong to `mqc_orbital_rotation`, which shares
   !! them with the CASSCF optimiser. What is left, and what this module is, is
   !! the three things specific to taking a Newton step on an SCF:
   !!
   !! * `soscf_gradient` -- the orbital gradient out of the Fock matrix.
   !! * `soscf_semicanonicalize` -- making the orbital energies mean something
   !!   again after a rotation, without moving the density.
   !! * `soscf_newton_step` -- the Newton step when the Hessian can only be
   !!   applied to a vector, never written down.
   !!
   !! ## The units, which are a factor of four and must not be guessed
   !!
   !! Parametrise the orbitals as `C -> C exp(kappa)` with `kappa`
   !! antisymmetric, as `mqc_orbital_rotation` and `mqc_czt_mcscf` do. For a
   !! closed shell the energy through second order in the non-redundant
   !! rotations `kappa_ai` is
   !!
   !!     E(kappa) = E0 + sum_ai g_ai kappa_ai
   !!                   + (1/2) sum_aibj kappa_ai H_ai,bj kappa_bj
   !!
   !!     g_ai   = 4 F_ai
   !!     H_ai,bj = 4 [ Delta_ai delta_ab delta_ij
   !!                   + 4(ai|bj) - (ab|ij) - (aj|ib) ]
   !!
   !! with `F` the closed-shell Fock matrix in the molecular basis and
   !! `Delta_ai = e_a - e_i`. The bracket is exactly `mqc_czt_ov_hessian`'s
   !! `(A+B)`, so
   !!
   !!     H = HESSIAN_SCALE * (A+B),   g = HESSIAN_SCALE * F_ai
   !!
   !! and the same factor multiplies both. It therefore **cancels out of the
   !! Newton step** `-H^-1 g` -- which is why a sign or factor error here is so
   !! easy to miss -- but it does not cancel out of the curvature, and
   !! `MIN_CURVATURE` and `SADDLE_CURVATURE` are thresholds on a curvature in
   !! hartree per radian squared. Working in the true units is what makes the
   !! constants shared with CASSCF mean the same thing on both paths. The
   !! factor is derived independently of this comment by
   !! `test_mqc_czt_soscf.f90`, which finite-differences the energy.
   !!
   !! The two-electron part of `H` is only the Hessian **when the Fock matrix
   !! is diagonal within the occupied block and within the virtual block** --
   !! otherwise the first term is `delta_ij F_ab - delta_ab F_ij` and not
   !! `Delta_ai`. A rotation destroys that, so every iteration
   !! semicanonicalises before it differentiates. That costs two small
   !! diagonalisations and no Fock build, and it leaves the density untouched,
   !! so it is not a step.
   !!
   !! ## Why the Newton equations are solved in a Krylov subspace
   !!
   !! `(A+B)` can be applied to a vector for one Fock build and cannot be
   !! written down for less than `n_occ*n_vir` of them, which is the whole
   !! matrix. So the step is taken in a small subspace built from the
   !! preconditioned gradient: project `H` and `g` into it, solve the
   !! level-shifted trust-region subproblem there with the *same*
   !! `level_shifted_step` the CASSCF optimiser uses on its dense Hessian, and
   !! expand along the preconditioned residual until the residual is small or
   !! the subspace is full.
   !!
   !! **`lowest` is an upper bound on the smallest curvature, not the smallest
   !! curvature.** A Rayleigh-Ritz value from a subspace always sits above the
   !! true eigenvalue. A negative one is therefore proof of negative curvature
   !! and a positive one is not proof of a minimum; the subspace is seeded with
   !! the softest diagonal direction so that the bound is usually tight, but
   !! `keywords.scf.stability` is what actually answers the question and this
   !! does not pretend to.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_ov_hessian, only: ov_hessian_t
   use mqc_orbital_rotation, only: rotation_hessian_t, subspace_newton_step
   implicit none
   private

   public :: soscf_gradient
   public :: soscf_semicanonicalize
   public :: soscf_newton_step
   public :: HESSIAN_SCALE
   public :: DEFAULT_SOSCF_SUBSPACE

   real(dp), parameter :: HESSIAN_SCALE = 4.0_dp
      !! `(A+B)` to the second derivative of the energy. The module header
      !! derives it; it multiplies the gradient by the same amount.
   integer, parameter :: DEFAULT_SOSCF_SUBSPACE = 20
      !! Krylov vectors the Newton equations are solved in, at most.
      !!
      !! Each one is a Fock build, so this is the per-iteration cost of the
      !! method and the reason a second-order SCF is not free. Twenty is enough
      !! for the preconditioned residual to fall two or three orders on the
      !! systems this was measured on, and the residual test below usually
      !! stops well short of it.
   real(dp), parameter :: SUBSPACE_TOLERANCE = 1.0e-3_dp
      !! How far the Newton equations are solved, relative to the gradient.
      !!
      !! Deliberately loose. The step is a direction to try, checked against
      !! the energy by the caller's backtracking, so solving the quadratic
      !! model to ten digits buys nothing -- the model itself is only good to
      !! the trust radius. Tightening it spends Fock builds and saves
      !! iterations, which is the wrong trade for the one measure that matters.

   type, extends(rotation_hessian_t) :: scaled_ov_hessian_t
      !! `ov_hessian_t` in true energy units: `HESSIAN_SCALE` times `(A+B)`
      type(ov_hessian_t), pointer :: ov => null()
         !! The caller's operator, viewed rather than owned: nothing here
         !! allocates or destroys it, so it must outlive this wrapper
   contains
      procedure :: apply => scaled_ov_apply
   end type scaled_ov_hessian_t

contains

   subroutine soscf_gradient(fock_mo, n_occ, gradient)
      !! `dE/d kappa_ai = 4 F_ai`, flattened virtual fastest
      !!
      !! The same layout `ov_hessian_t` acts on, so the gradient and the
      !! Hessian-vector products are the same vector space with no index
      !! translation in between.
      real(dp), intent(in) :: fock_mo(:, :)
         !! (n_mo, n_mo), the Fock matrix in the molecular basis
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: gradient(:)

      integer :: n_mo, n_vir, i, a

      n_mo = size(fock_mo, 1)
      n_vir = n_mo - n_occ
      allocate (gradient(max(n_vir*n_occ, 0)))
      do i = 1, n_occ
         do a = 1, n_vir
            gradient(a + (i - 1)*n_vir) = HESSIAN_SCALE*fock_mo(n_occ + a, i)
         end do
      end do
   end subroutine soscf_gradient

   subroutine soscf_semicanonicalize(fock_ao, coeff, n_occ, fock_mo, energies, error)
      !! Diagonalise the Fock matrix within the occupied and virtual blocks
      !!
      !! A rotation leaves the Fock matrix with off-diagonal elements inside
      !! each block, and the Hessian expression this module works in needs
      !! them gone -- see the header. Rotating *within* the occupied orbitals
      !! does not change a closed-shell density and rotating within the
      !! virtuals does not change anything at all, so this is a change of
      !! representation and not a step: the energy before and after is the same
      !! number, and no Fock matrix has to be rebuilt.
      !!
      !! What comes out is the Fock matrix in the new molecular basis, whose
      !! occupied-virtual block is the gradient and whose diagonal is the
      !! orbital energies. At convergence that block is zero and the orbitals
      !! are canonical, which is why `run_czt_rhf` can report them as such.
      real(dp), intent(in) :: fock_ao(:, :)
      real(dp), intent(inout) :: coeff(:, :)
         !! (n_ao, n_mo). Rotated in place, within each block.
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: fock_mo(:, :)
         !! (n_mo, n_mo), after the rotation
      real(dp), allocatable, intent(out) :: energies(:)
         !! (n_mo), its diagonal: occupied ascending, then virtual ascending
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: work(:, :), block_matrix(:, :), block_values(:)
      real(dp), allocatable :: rotated(:, :)
      integer :: n_ao, n_mo, n_vir, info, k

      if (error%has_error()) return
      n_ao = size(coeff, 1)
      n_mo = size(coeff, 2)
      n_vir = n_mo - n_occ
      if (n_occ < 1 .or. n_vir < 1) then
         call error%set(ERROR_VALIDATION, "a second-order SCF needs at least one "// &
                        "occupied and one virtual orbital; there are no rotations to "// &
                        "take")
         return
      end if

      allocate (work(n_ao, n_mo), fock_mo(n_mo, n_mo), energies(n_mo))
      allocate (rotated(n_ao, n_mo))

      ! The occupied block first, then the virtual one, each diagonalised on
      ! its own and used to rotate its own columns of C.
      call pic_gemm(fock_ao, coeff, work, beta=0.0_dp)
      call pic_gemm(coeff, work, fock_mo, transa="T", beta=0.0_dp)

      allocate (block_matrix(n_occ, n_occ), block_values(n_occ))
      block_matrix = fock_mo(1:n_occ, 1:n_occ)
      call pic_syev(block_matrix, block_values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "the occupied block of the Fock matrix "// &
                        "could not be diagonalized (info = "//to_char(info)//")")
         return
      end if
      call pic_gemm(coeff(:, 1:n_occ), block_matrix, rotated(:, 1:n_occ), beta=0.0_dp)
      deallocate (block_matrix, block_values)

      allocate (block_matrix(n_vir, n_vir), block_values(n_vir))
      block_matrix = fock_mo(n_occ + 1:n_mo, n_occ + 1:n_mo)
      call pic_syev(block_matrix, block_values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "the virtual block of the Fock matrix "// &
                        "could not be diagonalized (info = "//to_char(info)//")")
         return
      end if
      call pic_gemm(coeff(:, n_occ + 1:n_mo), block_matrix, rotated(:, n_occ + 1:n_mo), &
                    beta=0.0_dp)
      deallocate (block_matrix, block_values)

      coeff = rotated
      ! Rebuilt in the new basis rather than assembled from the two block
      ! rotations: the occupied-virtual block is the gradient and has to be the
      ! transformed one, not the one from before the rotation.
      call pic_gemm(fock_ao, coeff, work, beta=0.0_dp)
      call pic_gemm(coeff, work, fock_mo, transa="T", beta=0.0_dp)
      do k = 1, n_mo
         energies(k) = fock_mo(k, k)
      end do

      deallocate (work, rotated)
   end subroutine soscf_semicanonicalize

   subroutine soscf_newton_step(hessian, gradient, escape, kappa, lowest, predicted, &
                                products, error, max_subspace)
      !! The trust-region Newton step, from Hessian-vector products alone
      !!
      !! `mqc_orbital_rotation`'s `subspace_newton_step` on the scaled
      !! occupied-virtual Hessian, with the step laid out as `kappa`.
      type(ov_hessian_t), intent(inout), target :: hessian
      real(dp), intent(in) :: gradient(:)
         !! (n_vir*n_occ), in true energy units -- `soscf_gradient`'s
      real(dp), intent(in) :: escape
         !! How far to displace a mode with negative curvature and no gradient,
         !! in radians
      real(dp), allocatable, intent(out) :: kappa(:, :)
         !! (n_mo, n_mo), antisymmetric, ready for `rotation_matrix`
      real(dp), intent(out) :: lowest
         !! An **upper bound** on the smallest curvature, in hartree per radian
         !! squared. See the module header.
      real(dp), intent(out) :: predicted
         !! What the quadratic model says the step is worth, as a positive
         !! energy decrease
      integer, intent(out) :: products
         !! Hessian-vector products spent, each one a Fock build
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_subspace

      type(scaled_ov_hessian_t) :: scaled
      real(dp), allocatable :: step(:)
      integer :: n_ov, n_mo, n_occ, n_vir, nmax, i, a

      lowest = 0.0_dp
      predicted = 0.0_dp
      products = 0
      n_occ = hessian%n_occ
      n_vir = hessian%n_vir
      n_mo = hessian%n_mo
      allocate (kappa(n_mo, n_mo))
      kappa = 0.0_dp
      if (error%has_error()) return

      n_ov = hessian%length()
      if (size(gradient) /= n_ov) then
         call error%set(ERROR_VALIDATION, "the orbital gradient handed to the "// &
                        "second-order step is not the length of the non-redundant "// &
                        "occupied-virtual space")
         return
      end if
      if (n_ov < 1) return

      nmax = DEFAULT_SOSCF_SUBSPACE
      if (present(max_subspace)) nmax = max_subspace

      scaled%ov => hessian
      call subspace_newton_step(scaled, HESSIAN_SCALE*hessian%diagonal(), gradient, escape, &
                                step, lowest, predicted, products, error, &
                                max_subspace=nmax, tolerance=SUBSPACE_TOLERANCE)
      nullify (scaled%ov)
      if (error%has_error()) return

      do i = 1, n_occ
         do a = 1, n_vir
            kappa(n_occ + a, i) = step(a + (i - 1)*n_vir)
            kappa(i, n_occ + a) = -step(a + (i - 1)*n_vir)
         end do
      end do
   end subroutine soscf_newton_step

   subroutine scaled_ov_apply(this, x, hx, error)
      !! `H x` in true energy units: `HESSIAN_SCALE` times `(A+B) x`
      class(scaled_ov_hessian_t), intent(inout) :: this
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: hx(:)
      type(error_t), intent(inout) :: error

      if (error%has_error()) return
      call this%ov%apply(x, hx)
      if (this%ov%error%has_error()) then
         error = this%ov%error
         call error%add_context("applying the electronic Hessian for the "// &
                                "second-order SCF step")
         hx = 0.0_dp
         return
      end if
      hx = HESSIAN_SCALE*hx
   end subroutine scaled_ov_apply

end module mqc_czt_soscf
