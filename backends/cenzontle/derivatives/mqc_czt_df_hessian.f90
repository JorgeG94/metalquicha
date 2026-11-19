!! The two-electron derivatives a density-fitted reference's Hessian needs
module mqc_czt_df_hessian
   !! A fitted reference replaces `(mu nu|lambda sigma)` by
   !! `V_mn^T M^-1 V_ls`, with `V_mn,P = (mu nu|P)` and `M_PQ = (P|Q)`, and its
   !! Hessian differentiates that expression and not the integral it stands
   !! for. Two things come out of it, and this module builds both:
   !!
   !! * the Fock matrix's first derivative at fixed density, which is the
   !!   two-electron part of the coupled-perturbed right-hand side;
   !! * the energy's second derivative at fixed density.
   !!
   !! Everything is written with `c = M^-1 V`, the fitting coefficients, and
   !! with one object per perturbation `x`,
   !!
   !!     U^x = V^x - M^x c,
   !!
   !! the derivative of the three-index integrals less what the metric's own
   !! derivative takes back. In those terms, over the pair density
   !! `Gamma = D D - k/2 (D D)_exchange`,
   !!
   !!     (mn|ls)^x  = V^x_mn . c_ls + c_mn . U^x_ls
   !!     E^xy       = sum Gamma [ V^xy . c - 1/2 c . M^xy . c ]
   !!                + sum Gamma   U^x . M^-1 . U^y
   !!
   !! The first line of `E^xy` is the second-derivative integrals against two
   !! fixed densities, which `hess_3c2c_contract` does a shell triple at a
   !! time. The second is a product of first derivatives, and is where the
   !! fit's response to a displaced auxiliary basis lives.
   !!
   !! **Dense, and meant to be.** Each perturbation holds a full
   !! `(nao^2, naux)` array and the exchange term is `naux nao^3` per
   !! perturbation. That is the readable form the leaner one will be held to.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, atom_ao_blocks, three_centre, &
                                two_centre, metric_inverse_sqrt
   use mqc_czt_gradient, only: three_centre_deriv, two_centre_deriv
   use mqc_czt_hess_ints, only: hess_3c2c_contract
   implicit none
   private

   public :: df_hessian_terms_t
   public :: df_two_electron_derivatives
   public :: df_metric_dropped

   type :: df_hessian_terms_t
      !! What a fitted reference contributes to a Hessian, as the assembly reads it
      real(dp), allocatable :: b(:, :)
         !! The fitted tensor `V M^-1/2`, (nao^2, naux); the response operator's
         !! Coulomb and exchange are built from it.
      real(dp), allocatable :: h1(:, :, :, :)
         !! `d(J - k K/2)/dR_{A,x}` at fixed density, (nao, nao, 3, natm)
      real(dp), allocatable :: hess(:, :, :, :)
         !! Second derivative of `1/2 Tr[D (J - k K/2)]` at fixed density,
         !! (3, 3, natm, natm)
   contains
      procedure :: destroy => df_hessian_terms_destroy
   end type df_hessian_terms_t

   integer, parameter :: DERIV_ORBITAL = 1
      !! `three_centre_deriv`'s selector for `(nabla mu nu|P)`
   integer, parameter :: DERIV_AUXILIARY = 2
      !! and for `(mu nu|nabla P)`

contains

   subroutine df_hessian_terms_destroy(this)
      !! Release the three arrays
      class(df_hessian_terms_t), intent(inout) :: this

      if (allocated(this%b)) deallocate (this%b)
      if (allocated(this%h1)) deallocate (this%h1)
      if (allocated(this%hess)) deallocate (this%hess)
   end subroutine df_hessian_terms_destroy

   subroutine df_metric_dropped(aux, n_dropped, error)
      !! Modes the inverse of `aux`'s metric drops
      !!
      !! `-M^-1 M^x M^-1` is the derivative of the inverse only when none is
      !! dropped, so a nonzero count means the analytic Hessian is not that of
      !! the energy the SCF reported.
      type(czt_molecule_t), intent(in) :: aux
      integer, intent(out) :: n_dropped
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: metric(:, :), half(:, :)

      n_dropped = 0
      if (error%has_error()) return
      call two_centre(aux, metric)
      call metric_inverse_sqrt(metric, half, error, n_dropped=n_dropped)
   end subroutine df_metric_dropped

   subroutine df_two_electron_derivatives(orb, aux, density, c_occ, terms, error, k_scale)
      !! The fitted two-electron terms of a closed-shell Hessian at fixed density
      !!
      !! Fills all three components of `terms`. Refuses an auxiliary basis whose metric drops a mode; see
      !! `df_metric_dropped`.
      type(czt_molecule_t), intent(in) :: orb, aux
      real(dp), intent(in) :: density(:, :)   !! Closed shell, carrying its factor of two
      real(dp), intent(in) :: c_occ(:, :)     !! (nao, n_occ)
      type(df_hessian_terms_t), intent(out) :: terms
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: k_scale
         !! Exact-exchange fraction; one (Hartree-Fock) when absent.

      real(dp), allocatable :: v(:, :), metric(:, :), mih(:, :), minv(:, :), c(:, :)
      real(dp), allocatable :: d(:), dvec(:), dc(:, :, :), dcd(:, :), t(:, :)
      real(dp), allocatable :: gamma3(:, :, :), gamma2(:, :)
      real(dp), allocatable :: ip1(:, :, :, :), ip2(:, :, :, :), ip2c(:, :, :)
      real(dp), allocatable :: vx(:, :, :), mx(:, :), u(:, :), ubar(:), jx(:)
      real(dp), allocatable :: kx(:, :), work(:, :), half(:, :), y(:, :)
      real(dp), allocatable :: ut(:, :), yt(:, :, :)
      integer, allocatable :: orb_first(:), orb_count(:), aux_first(:), aux_count(:)
      real(dp) :: kf
      integer :: nao, naux, natm, n_occ, n_dropped, ia, ja, x, yy, p, q, mu, nu, ip, jp
      integer :: m0, m1, a0, a1

      if (error%has_error()) return
      kf = 1.0_dp
      if (present(k_scale)) kf = k_scale
      nao = orb%nao
      naux = aux%nao
      natm = orb%natm
      n_occ = size(c_occ, 2)
      if (aux%natm /= natm) then
         call error%set(ERROR_VALIDATION, "the fitted Hessian needs the auxiliary "// &
                        "basis on the same atoms as the orbital one.")
         return
      end if

      ! ---- the fit ---------------------------------------------------------
      call three_centre(orb, aux, v)
      call two_centre(aux, metric)
      call metric_inverse_sqrt(metric, mih, error, n_dropped=n_dropped)
      if (error%has_error()) return
      if (n_dropped > 0) then
         call error%set(ERROR_VALIDATION, "the auxiliary metric drops a mode, and "// &
                        "-M^-1 dM M^-1 is the derivative of its inverse only when it "// &
                        "drops none.")
         return
      end if
      allocate (minv(naux, naux), c(nao*nao, naux), terms%b(nao*nao, naux))
      call pic_gemm(mih, mih, minv)
      call pic_gemm(v, minv, c)
      call pic_gemm(v, mih, terms%b)
      deallocate (metric, minv)

      ! ---- the two fixed densities of the second-derivative integrals ------
      ! gamma3 = D d_P - k/2 D c_P D, gamma2 = -1/2 [d d^T - k/2 Tr(D c_P D c_Q)].
      allocate (dvec(nao*nao), d(naux))
      dvec = reshape(density, [nao*nao])
      d = matmul(dvec, c)

      allocate (dc(nao, nao, naux), dcd(nao*nao, naux), work(nao, nao))
      do p = 1, naux
         call pic_gemm(density, reshape(c(:, p), [nao, nao]), dc(:, :, p))
         call pic_gemm(dc(:, :, p), density, work)
         dcd(:, p) = reshape(work, [nao*nao])
      end do
      allocate (t(naux, naux))
      call pic_gemm(dcd, c, t, transa="T")

      allocate (gamma3(nao, nao, naux), gamma2(naux, naux))
      do p = 1, naux
         gamma3(:, :, p) = d(p)*density - 0.5_dp*kf*reshape(dcd(:, p), [nao, nao])
      end do
      do q = 1, naux
         do p = 1, naux
            gamma2(p, q) = -0.5_dp*(d(p)*d(q) - 0.5_dp*kf*t(p, q))
         end do
      end do
      deallocate (dcd, t)

      allocate (terms%hess(3, 3, natm, natm))
      terms%hess = 0.0_dp
      call hess_3c2c_contract(orb, aux, terms%hess, error, gamma3=gamma3, gamma2=gamma2)
      if (error%has_error()) return
      deallocate (gamma3, gamma2)

      ! ---- first derivatives, one perturbation at a time -------------------
      call three_centre_deriv(orb, aux, DERIV_ORBITAL, ip1)
      call three_centre_deriv(orb, aux, DERIV_AUXILIARY, ip2)
      call two_centre_deriv(aux, ip2c)

      allocate (orb_first(natm), orb_count(natm), aux_first(natm), aux_count(natm))
      call atom_ao_blocks(orb, orb_first, orb_count)
      call atom_ao_blocks(aux, aux_first, aux_count)

      allocate (terms%h1(nao, nao, 3, natm))
      allocate (vx(nao, nao, naux), mx(naux, naux), u(nao*nao, naux), ubar(naux))
      allocate (jx(nao*nao), kx(nao, nao), half(n_occ, nao), y(n_occ*n_occ, naux))
      allocate (ut(naux, 3*natm), yt(n_occ*n_occ, naux, 3*natm))

      do ia = 1, natm
         m0 = orb_first(ia) + 1
         m1 = orb_first(ia) + orb_count(ia)
         a0 = aux_first(ia) + 1
         a1 = aux_first(ia) + aux_count(ia)
         do x = 1, 3
            ip = 3*(ia - 1) + x

            ! `V^x`. The library's nabla acts on the electron, so the derivative
            ! with respect to a centre is minus it, for whichever of the three
            ! functions sits on this atom.
            vx = 0.0_dp
            do p = 1, naux
               do nu = 1, nao
                  do mu = m0, m1
                     vx(mu, nu, p) = vx(mu, nu, p) - ip1(mu, nu, p, x)
                     vx(nu, mu, p) = vx(nu, mu, p) - ip1(mu, nu, p, x)
                  end do
               end do
            end do
            if (a1 >= a0) vx(:, :, a0:a1) = vx(:, :, a0:a1) - ip2(:, :, a0:a1, x)

            ! `M^x`, the same rule for the metric's two functions.
            mx = 0.0_dp
            do q = 1, naux
               do p = a0, a1
                  mx(p, q) = mx(p, q) - ip2c(p, q, x)
                  mx(q, p) = mx(q, p) - ip2c(p, q, x)
               end do
            end do

            ! `U^x = V^x - M^x c`.
            u = reshape(vx, [nao*nao, naux])
            call pic_gemm(c, mx, u, alpha=-1.0_dp, beta=1.0_dp)
            ubar = matmul(dvec, u)

            ! Coulomb: `J^x = V^x d + c (U^x : D)`.
            jx = matmul(reshape(vx, [nao*nao, naux]), d) + matmul(c, ubar)

            ! Exchange: `K^x = sum_P [ V^x_P (D c_P) + (D c_P)^T U^x_P ]`.
            kx = 0.0_dp
            do p = 1, naux
               call pic_gemm(vx(:, :, p), dc(:, :, p), kx, alpha=1.0_dp, beta=1.0_dp)
               call pic_gemm(dc(:, :, p), reshape(u(:, p), [nao, nao]), kx, transa="T", &
                             alpha=1.0_dp, beta=1.0_dp)
            end do
            terms%h1(:, :, x, ia) = reshape(jx, [nao, nao]) - 0.5_dp*kf*kx

            ! What the product term of `E^xy` needs of this perturbation, with
            ! `M^-1/2` already on it: `U^x : D`, and `U^x` over occupied pairs.
            ut(:, ip) = matmul(mih, ubar)
            do p = 1, naux
               call pic_gemm(c_occ, reshape(u(:, p), [nao, nao]), half, transa="T")
               y(:, p) = reshape(matmul(half, c_occ), [n_occ*n_occ])
            end do
            call pic_gemm(y, mih, yt(:, :, ip))
         end do
      end do

      ! ---- sum Gamma U^x M^-1 U^y ------------------------------------------
      ! `D = 2 C C^T`, so the exchange trace is four times the sum over occupied
      ! pairs, and with its -k/2 that is -2k.
      do ja = 1, natm
         do ia = 1, natm
            do yy = 1, 3
               do x = 1, 3
                  ip = 3*(ia - 1) + x
                  jp = 3*(ja - 1) + yy
                  terms%hess(x, yy, ia, ja) = terms%hess(x, yy, ia, ja) &
                                              + dot_product(ut(:, ip), ut(:, jp)) &
                                              - 2.0_dp*kf*sum(yt(:, :, ip)*yt(:, :, jp))
               end do
            end do
         end do
      end do
   end subroutine df_two_electron_derivatives

end module mqc_czt_df_hessian
