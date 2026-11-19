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
   !! Written twice. `df_two_electron_derivatives_dense` is those formulas as
   !! they stand, a full `(nao^2, naux)` array per perturbation; it is the
   !! reference. `df_two_electron_derivatives` is what the Hessian calls: the
   !! same numbers with every exchange-type contraction taken through the
   !! occupied orbitals, and each atom touching only the integrals it owns.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_memory, only: memory_budget
   use mqc_czt_gemm_threads, only: gemm_over_columns, gemm_over_inner
   use mqc_czt_integrals, only: czt_molecule_t, atom_ao_blocks, three_centre, &
                                two_centre, metric_inverse_sqrt
   use mqc_czt_gradient, only: three_centre_deriv, two_centre_deriv
   use mqc_czt_hess_ints, only: hess_3c2c_contract
   implicit none
   private

   public :: df_hessian_terms_t
   public :: df_two_electron_derivatives
   public :: df_two_electron_derivatives_dense
   public :: df_fit_to_mo
   public :: df_response_mo
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
      real(dp), allocatable :: b_mo(:, :, :)
         !! `C^T B_P C`, (n_mo, n_mo, naux). Allocated by `df_fit_to_mo`, not by
         !! the derivative routines.
   contains
      procedure :: destroy => df_hessian_terms_destroy
   end type df_hessian_terms_t

   integer, parameter :: DERIV_ORBITAL = 1
      !! `three_centre_deriv`'s selector for `(nabla mu nu|P)`
   integer, parameter :: DERIV_AUXILIARY = 2
      !! and for `(mu nu|nabla P)`

   real(dp), parameter :: PRODUCT_BLIND_BYTES = 2.0_dp*1024.0_dp**3
      !! What the product term's batch plans on where the machine reports nothing
   real(dp), parameter :: PRODUCT_SHARE = 0.25_dp
      !! and the share of available memory it plans on where it does

contains

   subroutine df_hessian_terms_destroy(this)
      !! Release the three arrays
      class(df_hessian_terms_t), intent(inout) :: this

      if (allocated(this%b)) deallocate (this%b)
      if (allocated(this%h1)) deallocate (this%h1)
      if (allocated(this%hess)) deallocate (this%hess)
      if (allocated(this%b_mo)) deallocate (this%b_mo)
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

   subroutine df_two_electron_derivatives_dense(orb, aux, density, c_occ, terms, error, k_scale)
      !! `df_two_electron_derivatives`, one dense `(nao^2, naux)` array per perturbation
      !!
      !! Same arguments and the same three components. `naux nao^3` per
      !! perturbation and unthreaded: for checking the other against.
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
   end subroutine df_two_electron_derivatives_dense

   subroutine df_two_electron_derivatives(orb, aux, density, c_occ, terms, error, k_scale)
      !! The fitted two-electron terms of a closed-shell Hessian at fixed density
      !!
      !! Fills `b`, `h1` and `hess` of `terms`. Refuses an auxiliary basis whose
      !! metric drops a mode; see `df_metric_dropped`. `density` must be
      !! `2 c_occ c_occ^T`: the exchange terms are built from `c_occ` alone.
      type(czt_molecule_t), intent(in) :: orb, aux
      real(dp), intent(in) :: density(:, :)   !! Closed shell, carrying its factor of two
      real(dp), intent(in) :: c_occ(:, :)     !! (nao, n_occ)
      type(df_hessian_terms_t), intent(out) :: terms
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: k_scale
         !! Exact-exchange fraction; one (Hartree-Fock) when absent.

      ! Index order throughout: an AO first, then an occupied orbital, then the
      ! auxiliary function, so that `(n; o naux)` and `(n o; naux)` are both
      ! contiguous matrix views and an atom's auxiliary functions are a column
      ! range of either.
      real(dp), allocatable :: v(:, :), metric(:, :), mih(:, :), d(:)
      real(dp), allocatable, target :: b(:, :), half(:, :, :), cc(:, :, :), coo(:, :, :)
      real(dp), allocatable, target :: gamma3(:, :, :), gamma2(:, :), h1(:, :, :, :)
      real(dp), allocatable, target :: ip1(:, :, :, :), ip2(:, :, :, :), ip2c(:, :, :)
      real(dp), allocatable, target :: ip1c(:, :, :, :), f(:, :, :, :), g1(:, :, :)
      real(dp), allocatable, target :: it(:, :, :, :), pm(:, :, :), f_oo(:, :, :, :)
      real(dp), allocatable :: m1(:, :), w2(:, :), i1d(:, :, :), ubar(:, :), ut(:, :)
      real(dp), allocatable :: bx(:, :, :), kx(:, :), ca(:, :), gram(:, :)
      real(dp), pointer :: p2(:, :), q2(:, :), p3(:, :, :)
      integer, allocatable :: orb_first(:), orb_count(:), aux_first(:), aux_count(:)
      real(dp) :: kf
      integer :: nao, naux, natm, n_occ, n_dropped, ia, ja, x, yy, p, q, nu, i, ip, jp
      integer :: m0, m1a, a0, a1, na, xa

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
      allocate (b(nao*nao, naux))
      call gemm_over_columns(v, mih, b)
      deallocate (v, metric)

      ! `cc(mu, i, P) = sum_nu c_P(mu, nu) C(nu, i)`, the fitting coefficients
      ! `c = B M^-1/2` with one index on the occupied orbitals, and `coo` with
      ! both.
      allocate (half(nao, n_occ, naux), cc(nao, n_occ, naux), coo(n_occ, n_occ, naux))
      p3(1:nao, 1:nao, 1:naux) => b
      call contract_second(p3, c_occ, half)
      p2(1:nao*n_occ, 1:naux) => half
      q2(1:nao*n_occ, 1:naux) => cc
      call gemm_over_columns(p2, mih, q2)
      deallocate (half)
      p2(1:nao, 1:n_occ*naux) => cc
      q2(1:n_occ, 1:n_occ*naux) => coo
      call gemm_over_columns(c_occ, p2, q2, transa="T")

      ! `d_P = c_P : D`.
      allocate (d(naux))
      do p = 1, naux
         d(p) = 0.0_dp
         do i = 1, n_occ
            d(p) = d(p) + 2.0_dp*coo(i, i, p)
         end do
      end do

      ! ---- the two fixed densities of the second-derivative integrals ------
      ! gamma3 = D d_P - k/2 D c_P D and gamma2 = -1/2 [d d^T - k/2 Tr(D c_P D
      ! c_Q)], with `D c_P D = 4 C coo_P C^T`.
      allocate (gamma3(nao, nao, naux), gamma2(naux, naux))
      call occupied_to_ao(coo, c_occ, gamma3)
      !$omp parallel do default(none) shared(gamma3, density, d, kf, naux) private(p) &
      !$omp    schedule(static)
      do p = 1, naux
         gamma3(:, :, p) = d(p)*density - 2.0_dp*kf*gamma3(:, :, p)
      end do
      !$omp end parallel do
      p2(1:n_occ*n_occ, 1:naux) => coo
      call gemm_over_inner(p2, p2, gamma2)
      do q = 1, naux
         do p = 1, naux
            gamma2(p, q) = kf*gamma2(p, q) - 0.5_dp*d(p)*d(q)
         end do
      end do

      allocate (terms%hess(3, 3, natm, natm))
      terms%hess = 0.0_dp
      call hess_3c2c_contract(orb, aux, terms%hess, error, gamma3=gamma3, gamma2=gamma2)
      if (error%has_error()) return
      deallocate (gamma3, gamma2)

      allocate (orb_first(natm), orb_count(natm), aux_first(natm), aux_count(natm))
      call atom_ao_blocks(orb, orb_first, orb_count)
      call atom_ao_blocks(aux, aux_first, aux_count)
      allocate (h1(nao, nao, 3, natm), ubar(naux, 3*natm))
      h1 = 0.0_dp
      ubar = 0.0_dp

      ! ---- the metric's derivative -----------------------------------------
      ! `m1 = (nabla P|Q) d`, `pm = (nabla P|Q) M^-1/2`, and
      ! `f(mu, i, P) = sum_Q cc(mu, i, Q) (nabla P|Q)`, from which the
      ! auxiliary-centre integrals are subtracted below.
      call two_centre_deriv(aux, ip2c)
      allocate (m1(naux, 3), pm(naux, naux, 3), f(nao, n_occ, naux, 3))
      do x = 1, 3
         m1(:, x) = matmul(ip2c(:, :, x), d)
         call gemm_over_columns(ip2c(:, :, x), mih, pm(:, :, x))
      end do

      ! ---- (mu nu|nabla P), which only its own atom's perturbation sees -----
      call three_centre_deriv(orb, aux, DERIV_AUXILIARY, ip2)
      allocate (w2(naux, 3))
      do x = 1, 3
         do ia = 1, natm
            do p = aux_first(ia) + 1, aux_first(ia) + aux_count(ia)
               h1(:, :, x, ia) = h1(:, :, x, ia) - d(p)*ip2(:, :, p, x)
            end do
         end do
         call contract_second(ip2(:, :, :, x), c_occ, f(:, :, :, x))
         do p = 1, naux
            w2(p, x) = 2.0_dp*sum(f(:, :, p, x)*c_occ)
         end do
         p2(1:nao*n_occ, 1:naux) => cc
         q2(1:nao*n_occ, 1:naux) => f(:, :, :, x)
         call gemm_over_columns(p2, ip2c(:, :, x), q2, transb="T", alpha=1.0_dp, beta=-1.0_dp)
      end do
      deallocate (ip2)
      allocate (f_oo(n_occ, n_occ, naux, 3))
      do x = 1, 3
         p2(1:nao, 1:n_occ*naux) => f(:, :, :, x)
         q2(1:n_occ, 1:n_occ*naux) => f_oo(:, :, :, x)
         call gemm_over_columns(c_occ, p2, q2, transa="T")
      end do

      ! ---- (nabla mu nu|P) --------------------------------------------------
      ! `i1d = (nabla mu nu|P) d_P`, `ip1c` with `nu` on the occupied orbitals,
      ! `g1 = sum_P (nabla mu lambda|P) (D c_P)(lambda, nu)`, and `it`, which is
      ! `ip1c` with `M^-1/2` on its auxiliary index.
      call three_centre_deriv(orb, aux, DERIV_ORBITAL, ip1)
      allocate (i1d(nao, nao, 3), ip1c(nao, n_occ, naux, 3), g1(nao, nao, 3), &
                it(nao, n_occ, naux, 3))
      do x = 1, 3
         !$omp parallel do default(none) shared(i1d, ip1, d, x, nao, naux) private(nu, p) &
         !$omp    schedule(static)
         do nu = 1, nao
            i1d(:, nu, x) = 0.0_dp
            do p = 1, naux
               i1d(:, nu, x) = i1d(:, nu, x) + d(p)*ip1(:, nu, p, x)
            end do
         end do
         !$omp end parallel do
         call contract_second(ip1(:, :, :, x), c_occ, ip1c(:, :, :, x))
         p2(1:nao, 1:n_occ*naux) => ip1c(:, :, :, x)
         q2(1:nao, 1:n_occ*naux) => cc
         call gemm_over_columns(p2, q2, g1(:, :, x), transb="T", alpha=2.0_dp)
         p2(1:nao*n_occ, 1:naux) => ip1c(:, :, :, x)
         q2(1:nao*n_occ, 1:naux) => it(:, :, :, x)
         call gemm_over_columns(p2, mih, q2)
      end do

      ! ---- U^x : D for every perturbation, and the Coulomb part of h1 -------
      do ia = 1, natm
         m0 = orb_first(ia) + 1
         m1a = orb_first(ia) + orb_count(ia)
         a0 = aux_first(ia) + 1
         a1 = aux_first(ia) + aux_count(ia)
         do x = 1, 3
            ip = 3*(ia - 1) + x
            do p = 1, naux
               ubar(p, ip) = -4.0_dp*sum(ip1c(m0:m1a, :, p, x)*c_occ(m0:m1a, :)) &
                             + dot_product(ip2c(a0:a1, p, x), d(a0:a1))
            end do
            ubar(a0:a1, ip) = ubar(a0:a1, ip) - (w2(a0:a1, x) - m1(a0:a1, x))

            ! `V^x d`: the differentiated function is on this atom.
            h1(m0:m1a, :, x, ia) = h1(m0:m1a, :, x, ia) - i1d(m0:m1a, :, x)
            h1(:, m0:m1a, x, ia) = h1(:, m0:m1a, x, ia) - transpose(i1d(m0:m1a, :, x))
         end do
      end do
      deallocate (ip1c, i1d, w2, m1)
      allocate (ut(naux, 3*natm))
      call pic_gemm(mih, ubar, ut)
      p2(1:nao*nao, 1:3*natm) => h1
      call gemm_over_columns(b, ut, p2, alpha=1.0_dp, beta=1.0_dp)

      ! ---- exchange part of h1, an atom at a time ---------------------------
      ! K^x = kx + kx^T, with kx = -[mu in A] g1 - bx + 2 f_A cc_A^T and
      ! bx = sum_{lambda in A, P} (nabla lambda mu|P) (D c_P)(lambda, nu).
      allocate (bx(nao, nao, 3), kx(nao, nao))
      do ia = 1, natm
         na = orb_count(ia)
         xa = aux_count(ia)
         m0 = orb_first(ia) + 1
         m1a = orb_first(ia) + na
         a0 = aux_first(ia) + 1
         a1 = aux_first(ia) + xa
         bx = 0.0_dp
         if (na > 0) then
            allocate (ca(na, n_occ))
            ca = c_occ(m0:m1a, :)
            call owned_exchange(ip1, cc, ca, m0, bx)
            deallocate (ca)
         end if
         do x = 1, 3
            kx = -bx(:, :, x)
            kx(m0:m1a, :) = kx(m0:m1a, :) - g1(m0:m1a, :, x)
            if (xa > 0) then
               p2(1:nao, 1:n_occ*xa) => f(:, :, a0:a1, x)
               q2(1:nao, 1:n_occ*xa) => cc(:, :, a0:a1)
               call gemm_over_columns(p2, q2, kx, transb="T", alpha=2.0_dp, beta=1.0_dp)
            end if
            h1(:, :, x, ia) = h1(:, :, x, ia) - 0.5_dp*kf*(kx + transpose(kx))
         end do
      end do
      deallocate (bx, kx, ip1, g1, f, cc)

      ! ---- sum Gamma U^x M^-1 U^y ------------------------------------------
      ! `D = 2 C C^T`, so the exchange trace is four times the sum over occupied
      ! pairs, and with its -k/2 that is -2k.
      allocate (gram(3*natm, 3*natm))
      call occupied_product(c_occ, it, f_oo, coo, mih, pm, orb_first, orb_count, &
                            aux_first, aux_count, gram)
      do ja = 1, natm
         do ia = 1, natm
            do yy = 1, 3
               do x = 1, 3
                  ip = 3*(ia - 1) + x
                  jp = 3*(ja - 1) + yy
                  terms%hess(x, yy, ia, ja) = terms%hess(x, yy, ia, ja) &
                                              + dot_product(ut(:, ip), ut(:, jp)) &
                                              - 2.0_dp*kf*gram(ip, jp)
               end do
            end do
         end do
      end do

      call move_alloc(b, terms%b)
      call move_alloc(h1, terms%h1)
   end subroutine df_two_electron_derivatives

   subroutine contract_second(t, c_occ, out)
      !! `out(mu, i, P) = sum_nu t(mu, nu, P) c_occ(nu, i)`
      real(dp), intent(in) :: t(:, :, :)        !! (nao, nao, naux)
      real(dp), intent(in) :: c_occ(:, :)       !! (nao, n_occ)
      real(dp), intent(out) :: out(:, :, :)     !! (nao, n_occ, naux)

      integer :: p

      !$omp parallel do default(none) shared(t, c_occ, out) private(p) schedule(static)
      do p = 1, size(t, 3)
         call pic_gemm(t(:, :, p), c_occ, out(:, :, p))
      end do
      !$omp end parallel do
   end subroutine contract_second

   subroutine occupied_to_ao(coo, c_occ, out)
      !! `out(:, :, P) = c_occ coo(:, :, P) c_occ^T`
      real(dp), intent(in) :: coo(:, :, :)      !! (n_occ, n_occ, naux)
      real(dp), intent(in) :: c_occ(:, :)       !! (nao, n_occ)
      real(dp), intent(out) :: out(:, :, :)     !! (nao, nao, naux)

      real(dp), allocatable :: work(:, :)
      integer :: p

      !$omp parallel default(none) shared(coo, c_occ, out) private(p, work)
      allocate (work(size(c_occ, 1), size(c_occ, 2)))
      !$omp do schedule(static)
      do p = 1, size(coo, 3)
         call pic_gemm(c_occ, coo(:, :, p), work)
         call pic_gemm(work, c_occ, out(:, :, p), transb="T")
      end do
      !$omp end do
      deallocate (work)
      !$omp end parallel
   end subroutine occupied_to_ao

   subroutine owned_exchange(ip1, cc, ca, m0, bx)
      !! `bx(mu, nu, x) = sum_{lambda in A, P} (nabla_x lambda mu|P) (D c_P)(lambda, nu)`
      !!
      !! `A` is the atom whose functions are rows `m0` onward of `ip1`, as many
      !! as `ca` has rows. `nao^2 n_A` per auxiliary function, threaded over them.
      real(dp), intent(in) :: ip1(:, :, :, :)   !! (nao, nao, naux, 3)
      real(dp), intent(in) :: cc(:, :, :)       !! (nao, n_occ, naux)
      real(dp), intent(in) :: ca(:, :)          !! `c_occ`'s rows on the atom, (n_A, n_occ)
      integer, intent(in) :: m0
      real(dp), intent(out) :: bx(:, :, :)      !! (nao, nao, 3)

      real(dp), allocatable :: dca(:, :), slab(:, :), part(:, :, :)
      integer :: nao, na, naux, p, x

      nao = size(ip1, 1)
      naux = size(ip1, 3)
      na = size(ca, 1)
      bx = 0.0_dp
      !$omp parallel default(none) shared(ip1, cc, ca, bx, nao, na, naux, m0) &
      !$omp    private(p, x, dca, slab, part)
      allocate (dca(na, nao), slab(na, nao), part(nao, nao, 3))
      part = 0.0_dp
      !$omp do schedule(static)
      do p = 1, naux
         call pic_gemm(ca, cc(:, :, p), dca, transb="T", alpha=2.0_dp)
         do x = 1, 3
            slab = ip1(m0:m0 + na - 1, :, p, x)
            call pic_gemm(slab, dca, part(:, :, x), transa="T", alpha=1.0_dp, beta=1.0_dp)
         end do
      end do
      !$omp end do
      !$omp critical (df_hessian_owned_exchange)
      bx = bx + part
      !$omp end critical (df_hessian_owned_exchange)
      deallocate (dca, slab, part)
      !$omp end parallel
   end subroutine owned_exchange

   subroutine occupied_product(c_occ, it, fmo, coo, mih, pm, orb_first, orb_count, &
                               aux_first, aux_count, gram)
      !! `gram(x, y) = sum_{ij,R} Yt^x(i, j, R) Yt^y(i, j, R)` over all perturbations
      !!
      !! `Yt^x = (C^T U^x_P C) M^-1/2`, assembled for a batch of its auxiliary
      !! index at a time so that no more than a share of memory holds it:
      !!
      !!     Yt^x(i,j,R) = - sum_{mu in A} [C(mu,i) it(mu,j,R) + C(mu,j) it(mu,i,R)]
      !!                 + sum_{P in A} [fmo(i,j,P) mih(P,R) + coo(i,j,P) pm(P,R)]
      !!
      !! with `it`, `fmo` and `pm` those of the perturbation's direction.
      real(dp), intent(in) :: c_occ(:, :)             !! (nao, n_occ)
      real(dp), intent(in) :: it(:, :, :, :)          !! (nao, n_occ, naux, 3)
      real(dp), intent(in), contiguous, target :: fmo(:, :, :, :)
         !! (n_occ, n_occ, naux, 3)
      real(dp), intent(in), contiguous, target :: coo(:, :, :)
         !! (n_occ, n_occ, naux)
      real(dp), intent(in) :: mih(:, :)               !! (naux, naux)
      real(dp), intent(in) :: pm(:, :, :)             !! (naux, naux, 3)
      integer, intent(in) :: orb_first(:), orb_count(:), aux_first(:), aux_count(:)
      real(dp), intent(out) :: gram(:, :)             !! (3 natm, 3 natm)

      real(dp), allocatable, target :: yt(:, :, :, :)
      real(dp), allocatable, target :: slab(:, :, :), qq(:, :, :)
      real(dp), allocatable :: ca(:, :), right(:, :), part(:, :)
      real(dp), pointer :: p2(:, :), q2(:, :), y2(:, :)
      integer :: nao, n_occ, naux, natm, npert, width, r0, r1, w, ip, ia, x, i, j, r
      integer :: m0, na, a0, xa

      nao = size(c_occ, 1)
      n_occ = size(c_occ, 2)
      naux = size(mih, 1)
      natm = size(orb_first)
      npert = 3*natm
      ! Clamped as a real: a small molecule on a large machine gives a figure
      ! no default integer holds.
      width = int(min(real(naux, dp), &
                      max(1.0_dp, memory_budget(PRODUCT_BLIND_BYTES, PRODUCT_SHARE) &
                          /(8.0_dp*real(n_occ, dp)**2*real(npert, dp)))))

      gram = 0.0_dp
      allocate (part(npert, npert))
      do r0 = 1, naux, width
         r1 = min(r0 + width - 1, naux)
         w = r1 - r0 + 1
         allocate (yt(n_occ, n_occ, w, npert))
         !$omp parallel default(none) &
         !$omp    shared(c_occ, it, fmo, coo, mih, pm, orb_first, orb_count, aux_first) &
         !$omp    shared(aux_count, yt, nao, n_occ, natm, npert, r0, r1, w) &
         !$omp    private(ip, ia, x, i, j, r, m0, na, a0, xa, slab, qq, ca, right, p2, q2, y2)
         !$omp do schedule(dynamic)
         do ip = 1, npert
            ia = (ip - 1)/3 + 1
            x = ip - 3*(ia - 1)
            na = orb_count(ia)
            m0 = orb_first(ia) + 1
            xa = aux_count(ia)
            a0 = aux_first(ia) + 1
            yt(:, :, :, ip) = 0.0_dp
            if (na > 0) then
               allocate (slab(na, n_occ, w), qq(n_occ, n_occ, w), ca(na, n_occ))
               ca = c_occ(m0:m0 + na - 1, :)
               slab = it(m0:m0 + na - 1, :, r0:r1, x)
               p2(1:na, 1:n_occ*w) => slab
               q2(1:n_occ, 1:n_occ*w) => qq
               call pic_gemm(ca, p2, q2, transa="T")
               do r = 1, w
                  do j = 1, n_occ
                     do i = 1, n_occ
                        yt(i, j, r, ip) = -(qq(i, j, r) + qq(j, i, r))
                     end do
                  end do
               end do
               deallocate (slab, qq, ca)
            end if
            if (xa > 0) then
               allocate (right(xa, w))
               y2(1:n_occ*n_occ, 1:w) => yt(:, :, :, ip)
               right = mih(a0:a0 + xa - 1, r0:r1)
               p2(1:n_occ*n_occ, 1:xa) => fmo(:, :, a0:a0 + xa - 1, x)
               call pic_gemm(p2, right, y2, alpha=1.0_dp, beta=1.0_dp)
               right = pm(a0:a0 + xa - 1, r0:r1, x)
               p2(1:n_occ*n_occ, 1:xa) => coo(:, :, a0:a0 + xa - 1)
               call pic_gemm(p2, right, y2, alpha=1.0_dp, beta=1.0_dp)
               deallocate (right)
            end if
         end do
         !$omp end do
         !$omp end parallel
         y2(1:n_occ*n_occ*w, 1:npert) => yt
         call gemm_over_inner(y2, y2, part)
         gram = gram + part
         deallocate (yt)
      end do
   end subroutine occupied_product

   subroutine df_fit_to_mo(terms, orbitals)
      !! Fill `terms%b_mo` with `C^T B_P C`, from `terms%b`
      type(df_hessian_terms_t), intent(inout), target :: terms
      real(dp), intent(in) :: orbitals(:, :)   !! (nao, n_mo)

      real(dp), allocatable :: work(:, :)
      real(dp), pointer :: b3(:, :, :)
      integer :: nao, n_mo, naux, p

      nao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      naux = size(terms%b, 2)
      if (allocated(terms%b_mo)) deallocate (terms%b_mo)
      allocate (terms%b_mo(n_mo, n_mo, naux))
      b3(1:nao, 1:nao, 1:naux) => terms%b
      !$omp parallel default(none) shared(terms, orbitals, b3, nao, n_mo, naux) &
      !$omp    private(p, work)
      allocate (work(nao, n_mo))
      !$omp do schedule(static)
      do p = 1, naux
         call pic_gemm(b3(:, :, p), orbitals, work)
         call pic_gemm(orbitals, work, terms%b_mo(:, :, p), transa="T")
      end do
      !$omp end do
      deallocate (work)
      !$omp end parallel
   end subroutine df_fit_to_mo

   subroutine df_response_mo(b_mo, mo1, gmo, k_scale)
      !! `C^T (J - k K/2)[D1] C_occ` for the densities a batch of trial vectors implies
      !!
      !! `D1 = 2 C (mo1 E^T + E mo1^T) C^T`, the symmetrised closed-shell
      !! response density with `E` the occupied columns of the identity. Nothing
      !! is taken to the atomic basis: `n_mo^2 n_occ` per auxiliary function and
      !! vector, against the `n^3` of a general density.
      real(dp), intent(in) :: b_mo(:, :, :)   !! `C^T B_P C`, (n_mo, n_mo, naux)
      real(dp), intent(in) :: mo1(:, :, :)    !! (n_mo, n_occ, n_vec), occupied rows included
      real(dp), intent(out) :: gmo(:, :, :)   !! (n_mo, n_occ, n_vec)
      real(dp), intent(in), optional :: k_scale
         !! Exact-exchange fraction; one when absent.

      real(dp), allocatable :: bo(:, :), boo(:, :), z(:, :), w(:, :), part(:, :)
      real(dp) :: kf, c_p
      integer :: n_mo, n_occ, n_vec, naux, p, m, i0, i1

      kf = 1.0_dp
      if (present(k_scale)) kf = k_scale
      n_mo = size(mo1, 1)
      n_occ = size(mo1, 2)
      n_vec = size(mo1, 3)
      naux = size(b_mo, 3)
      gmo = 0.0_dp

      ! With `bo` the occupied columns of `B_P` and `boo` its occupied block,
      !     J = 4 (bo : mo1) bo,     K = 2 [B_P (mo1 boo) + bo (mo1^T bo)].
      ! `z` and `part` hold every vector side by side, `n_occ` columns each, so
      ! the one `n_mo^2` product per auxiliary function is a single wide one.
      !$omp parallel default(none) shared(b_mo, mo1, gmo, kf, n_mo, n_occ, n_vec, naux) &
      !$omp    private(p, m, i0, i1, c_p, bo, boo, z, w, part)
      allocate (bo(n_mo, n_occ), boo(n_occ, n_occ), z(n_mo, n_occ*n_vec), &
                w(n_occ, n_occ), part(n_mo, n_occ*n_vec))
      part = 0.0_dp
      !$omp do schedule(static)
      do p = 1, naux
         bo = b_mo(:, 1:n_occ, p)
         boo = bo(1:n_occ, :)
         do m = 1, n_vec
            i0 = (m - 1)*n_occ + 1
            i1 = m*n_occ
            c_p = 4.0_dp*sum(bo*mo1(:, :, m))
            part(:, i0:i1) = part(:, i0:i1) + c_p*bo
            call pic_gemm(mo1(:, :, m), boo, z(:, i0:i1))
            call pic_gemm(mo1(:, :, m), bo, w, transa="T")
            call pic_gemm(bo, w, part(:, i0:i1), alpha=-kf, beta=1.0_dp)
         end do
         call pic_gemm(b_mo(:, :, p), z, part, alpha=-kf, beta=1.0_dp)
      end do
      !$omp end do
      !$omp critical (df_hessian_response_mo)
      do m = 1, n_vec
         gmo(:, :, m) = gmo(:, :, m) + part(:, (m - 1)*n_occ + 1:m*n_occ)
      end do
      !$omp end critical (df_hessian_response_mo)
      deallocate (bo, boo, z, w, part)
      !$omp end parallel
   end subroutine df_response_mo

end module mqc_czt_df_hessian
