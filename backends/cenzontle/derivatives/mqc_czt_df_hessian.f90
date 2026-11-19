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
   !! None of `V^x`, `M^x` or `U^x` is ever formed. Every exchange-type
   !! contraction is taken through the occupied orbitals, and each atom touches
   !! only the integrals whose differentiated function it owns. The formulas as
   !! they stand, one `(nao^2, naux)` array per perturbation, are written out in
   !! `test_mqc_czt_df_hessian`, which holds this module to them.
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
   public :: df_fit_to_mo
   public :: df_operator_t
   public :: df_operator_setup
   public :: df_potential_derivative
   public :: df_cross_second_derivative
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

   type :: df_operator_t
      !! The fit and its first-derivative integrals, held for a caller that
      !! differentiates `J - k K/2` of several densities in turn
      integer :: nao = 0
      integer :: naux = 0
      integer :: natm = 0
      real(dp), allocatable :: c(:, :, :)
         !! Fitting coefficients `M^-1 (P|mu nu)`, (nao, nao, naux)
      real(dp), allocatable :: mih(:, :)
         !! `M^-1/2`, (naux, naux)
      real(dp), allocatable :: ip1(:, :, :, :)
         !! `(nabla mu nu|P)`, (nao, nao, naux, 3)
      real(dp), allocatable :: ip2(:, :, :, :)
         !! `(mu nu|nabla P)`, (nao, nao, naux, 3)
      real(dp), allocatable :: ip2c(:, :, :)
         !! `(nabla P|Q)`, (naux, naux, 3)
      integer, allocatable :: orb_first(:), orb_count(:)
         !! Each atom's orbital functions: zero-based offset and count
      integer, allocatable :: aux_first(:), aux_count(:)
         !! and its auxiliary ones
   contains
      procedure :: destroy => df_operator_destroy
   end type df_operator_t

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

   subroutine df_operator_destroy(this)
      !! Release every array
      class(df_operator_t), intent(inout) :: this

      if (allocated(this%c)) deallocate (this%c)
      if (allocated(this%mih)) deallocate (this%mih)
      if (allocated(this%ip1)) deallocate (this%ip1)
      if (allocated(this%ip2)) deallocate (this%ip2)
      if (allocated(this%ip2c)) deallocate (this%ip2c)
      if (allocated(this%orb_first)) deallocate (this%orb_first, this%orb_count)
      if (allocated(this%aux_first)) deallocate (this%aux_first, this%aux_count)
   end subroutine df_operator_destroy

   subroutine df_operator_setup(orb, aux, ops, error)
      !! Fill `ops` for the pair of basis sets
      !!
      !! Refuses an auxiliary basis whose metric drops a mode; see
      !! `df_metric_dropped`. Seven arrays of `nao^2 naux` stay allocated.
      type(czt_molecule_t), intent(in) :: orb, aux
      type(df_operator_t), intent(out) :: ops
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: v(:, :), metric(:, :), minv(:, :)
      real(dp), allocatable, target :: c2(:, :)
      real(dp), pointer :: c3(:, :, :)
      integer :: n_dropped

      if (error%has_error()) return
      ops%nao = orb%nao
      ops%naux = aux%nao
      ops%natm = orb%natm
      if (aux%natm /= orb%natm) then
         call error%set(ERROR_VALIDATION, "the fitted operator's derivatives need the "// &
                        "auxiliary basis on the same atoms as the orbital one.")
         return
      end if
      call three_centre(orb, aux, v)
      call two_centre(aux, metric)
      call metric_inverse_sqrt(metric, ops%mih, error, n_dropped=n_dropped)
      if (error%has_error()) return
      if (n_dropped > 0) then
         call error%set(ERROR_VALIDATION, "the auxiliary metric drops a mode, and "// &
                        "-M^-1 dM M^-1 is the derivative of its inverse only when it "// &
                        "drops none.")
         return
      end if
      allocate (minv(ops%naux, ops%naux), c2(ops%nao*ops%nao, ops%naux))
      call pic_gemm(ops%mih, ops%mih, minv)
      call gemm_over_columns(v, minv, c2)
      deallocate (v, metric, minv)
      c3(1:ops%nao, 1:ops%nao, 1:ops%naux) => c2
      ops%c = c3
      deallocate (c2)

      call three_centre_deriv(orb, aux, DERIV_ORBITAL, ops%ip1)
      call three_centre_deriv(orb, aux, DERIV_AUXILIARY, ops%ip2)
      call two_centre_deriv(aux, ops%ip2c)
      allocate (ops%orb_first(ops%natm), ops%orb_count(ops%natm), &
                ops%aux_first(ops%natm), ops%aux_count(ops%natm))
      call atom_ao_blocks(orb, ops%orb_first, ops%orb_count)
      call atom_ao_blocks(aux, ops%aux_first, ops%aux_count)
   end subroutine df_operator_setup

   subroutine fit_response_vectors(ops, dens, d, ubar)
      !! `d_P = c_P : D` and `ubar(P, x) = U^x_P : D` for every perturbation
      type(df_operator_t), intent(in) :: ops
      real(dp), intent(in) :: dens(:, :)      !! Symmetric, (nao, nao)
      real(dp), intent(out) :: d(:)           !! (naux)
      real(dp), intent(out) :: ubar(:, :)     !! (naux, 3 natm)

      real(dp), allocatable :: w1(:, :), w2(:), m1(:)
      integer :: ia, x, p, ip, m0, m1a, a0, a1

      !$omp parallel do default(none) shared(ops, dens, d) private(p) schedule(static)
      do p = 1, ops%naux
         d(p) = sum(ops%c(:, :, p)*dens)
      end do
      !$omp end parallel do

      allocate (w1(ops%nao, ops%naux), w2(ops%naux), m1(ops%naux))
      do x = 1, 3
         ! `w1(mu, P) = sum_nu (nabla mu nu|P) D(mu, nu)`, `w2 = (mu nu|nabla P) : D`
         ! and `m1 = (nabla P|Q) d`.
         !$omp parallel do default(none) shared(ops, dens, w1, w2, x) private(p) &
         !$omp    schedule(static)
         do p = 1, ops%naux
            w1(:, p) = sum(ops%ip1(:, :, p, x)*dens, dim=2)
            w2(p) = sum(ops%ip2(:, :, p, x)*dens)
         end do
         !$omp end parallel do
         m1 = matmul(ops%ip2c(:, :, x), d)
         do ia = 1, ops%natm
            ip = 3*(ia - 1) + x
            m0 = ops%orb_first(ia) + 1
            m1a = ops%orb_first(ia) + ops%orb_count(ia)
            a0 = ops%aux_first(ia) + 1
            a1 = ops%aux_first(ia) + ops%aux_count(ia)
            do p = 1, ops%naux
               ubar(p, ip) = -2.0_dp*sum(w1(m0:m1a, p)) &
                             + dot_product(ops%ip2c(a0:a1, p, x), d(a0:a1))
            end do
            ubar(a0:a1, ip) = ubar(a0:a1, ip) - (w2(a0:a1) - m1(a0:a1))
         end do
      end do
   end subroutine fit_response_vectors

   subroutine df_potential_derivative(ops, dens, h1, k_scale)
      !! `d(J - k K/2)[D]/dR` for every nuclear coordinate, at fixed `D`
      !!
      !! `D` is any symmetric matrix in the atomic basis: a relaxed or a
      !! correlation density as well as the reference's own. About two dozen
      !! `nao^3 naux` products, threaded.
      type(df_operator_t), intent(in), target :: ops
      real(dp), intent(in) :: dens(:, :)
      real(dp), intent(out), contiguous, target :: h1(:, :, :, :)
         !! (nao, nao, 3, natm)
      real(dp), intent(in), optional :: k_scale
         !! Exact-exchange fraction; one when absent.

      real(dp), allocatable, target :: dc(:, :, :), de(:, :, :)
      real(dp), allocatable :: d(:), ubar(:, :), i1d(:, :), g1(:, :), kx(:, :), bx(:, :)
      real(dp), pointer :: p2(:, :), q2(:, :)
      real(dp) :: kf
      integer :: nao, naux, natm, ia, x, p, nu, m0, m1a, a0, a1, na

      kf = 1.0_dp
      if (present(k_scale)) kf = k_scale
      nao = ops%nao
      naux = ops%naux
      natm = ops%natm
      allocate (d(naux), ubar(naux, 3*natm))
      call fit_response_vectors(ops, dens, d, ubar)

      ! ---- Coulomb: J^x = V^x d + c (U^x : D) --------------------------------
      h1 = 0.0_dp
      allocate (i1d(nao, nao))
      do x = 1, 3
         !$omp parallel do default(none) shared(ops, d, i1d, x, nao, naux) private(nu, p) &
         !$omp    schedule(static)
         do nu = 1, nao
            i1d(:, nu) = 0.0_dp
            do p = 1, naux
               i1d(:, nu) = i1d(:, nu) + d(p)*ops%ip1(:, nu, p, x)
            end do
         end do
         !$omp end parallel do
         do ia = 1, natm
            m0 = ops%orb_first(ia) + 1
            m1a = ops%orb_first(ia) + ops%orb_count(ia)
            h1(m0:m1a, :, x, ia) = h1(m0:m1a, :, x, ia) - i1d(m0:m1a, :)
            h1(:, m0:m1a, x, ia) = h1(:, m0:m1a, x, ia) - transpose(i1d(m0:m1a, :))
            do p = ops%aux_first(ia) + 1, ops%aux_first(ia) + ops%aux_count(ia)
               h1(:, :, x, ia) = h1(:, :, x, ia) - d(p)*ops%ip2(:, :, p, x)
            end do
         end do
      end do
      p2(1:nao*nao, 1:naux) => ops%c
      q2(1:nao*nao, 1:3*natm) => h1
      call gemm_over_columns(p2, ubar, q2, alpha=1.0_dp, beta=1.0_dp)
      deallocate (i1d, ubar, d)
      if (kf == 0.0_dp) return

      ! ---- exchange: K^x = kx + kx^T ----------------------------------------
      ! With dc_P = D c_P and de_P = sum_Q (nabla P|Q) dc_Q,
      !   kx = - [mu in A] sum_P (nabla mu .|P) dc_P
      !        - sum_{lambda in A, P} (nabla lambda mu|P) dc_P(lambda, nu)
      !        - sum_{P in A} (. .|nabla P) dc_P + sum_{P in A} c_P de_P.
      allocate (dc(nao, nao, naux), de(nao, nao, naux), g1(nao, nao), kx(nao, nao), &
                bx(nao, nao))
      !$omp parallel do default(none) shared(ops, dens, dc, naux) private(p) schedule(static)
      do p = 1, naux
         call pic_gemm(dens, ops%c(:, :, p), dc(:, :, p))
      end do
      !$omp end parallel do
      do x = 1, 3
         call slab_products(ops%ip1(:, :, :, x), dc, 1, naux, g1)
         p2(1:nao*nao, 1:naux) => dc
         q2(1:nao*nao, 1:naux) => de
         call gemm_over_columns(p2, ops%ip2c(:, :, x), q2, transb="T")
         do ia = 1, natm
            na = ops%orb_count(ia)
            m0 = ops%orb_first(ia) + 1
            m1a = ops%orb_first(ia) + na
            a0 = ops%aux_first(ia) + 1
            a1 = ops%aux_first(ia) + ops%aux_count(ia)
            kx = 0.0_dp
            if (na > 0) then
               call owned_rows_product(ops%ip1(:, :, :, x), dc, m0, na, bx)
               kx = -bx
               kx(m0:m1a, :) = kx(m0:m1a, :) - g1(m0:m1a, :)
            end if
            if (a1 >= a0) then
               call slab_products(ops%ip2(:, :, :, x), dc, a0, a1, bx)
               kx = kx - bx
               call slab_products(ops%c, de, a0, a1, bx)
               kx = kx + bx
            end if
            h1(:, :, x, ia) = h1(:, :, x, ia) - 0.5_dp*kf*(kx + transpose(kx))
         end do
      end do
   end subroutine df_potential_derivative

   subroutine slab_products(a, b, p0, p1, out)
      !! `out = sum_{P = p0..p1} a(:, :, P) b(:, :, P)`, threaded over `P`
      real(dp), intent(in) :: a(:, :, :), b(:, :, :)
      integer, intent(in) :: p0, p1
      real(dp), intent(out) :: out(:, :)

      real(dp), allocatable :: part(:, :)
      integer :: p

      out = 0.0_dp
      !$omp parallel default(none) shared(a, b, out, p0, p1) private(p, part)
      allocate (part(size(out, 1), size(out, 2)))
      part = 0.0_dp
      !$omp do schedule(static)
      do p = p0, p1
         call pic_gemm(a(:, :, p), b(:, :, p), part, alpha=1.0_dp, beta=1.0_dp)
      end do
      !$omp end do
      !$omp critical (df_hessian_slab_products)
      out = out + part
      !$omp end critical (df_hessian_slab_products)
      deallocate (part)
      !$omp end parallel
   end subroutine slab_products

   subroutine owned_rows_product(a, b, m0, na, out)
      !! `out(mu, nu) = sum_P sum_{lambda = m0..m0+na-1} a(lambda, mu, P) b(lambda, nu, P)`
      real(dp), intent(in) :: a(:, :, :), b(:, :, :)
      integer, intent(in) :: m0, na
      real(dp), intent(out) :: out(:, :)

      real(dp), allocatable :: part(:, :), sa(:, :), sb(:, :)
      integer :: p

      out = 0.0_dp
      !$omp parallel default(none) shared(a, b, out, m0, na) private(p, part, sa, sb)
      allocate (part(size(out, 1), size(out, 2)), sa(na, size(a, 2)), sb(na, size(b, 2)))
      part = 0.0_dp
      !$omp do schedule(static)
      do p = 1, size(a, 3)
         sa = a(m0:m0 + na - 1, :, p)
         sb = b(m0:m0 + na - 1, :, p)
         call pic_gemm(sa, sb, part, transa="T", alpha=1.0_dp, beta=1.0_dp)
      end do
      !$omp end do
      !$omp critical (df_hessian_owned_rows)
      out = out + part
      !$omp end critical (df_hessian_owned_rows)
      deallocate (part, sa, sb)
      !$omp end parallel
   end subroutine owned_rows_product

   subroutine df_cross_second_derivative(ops, orb, aux, dens, c_occ, hess, error, k_scale)
      !! Add the second derivative of `A . (J - k K/2)[B]` at fixed `A` and `B`
      !!
      !! `A` is `dens`, any symmetric matrix; `B = 2 c_occ c_occ^T`, a
      !! closed-shell density given by its orbitals. Symmetric in the two
      !! nuclear coordinates. With `A = B` it is twice the `hess` that
      !! `df_two_electron_derivatives` returns, which is the derivative of half
      !! this.
      type(df_operator_t), intent(in), target :: ops
      type(czt_molecule_t), intent(in) :: orb, aux
      real(dp), intent(in) :: dens(:, :)        !! `A`, (nao, nao)
      real(dp), intent(in) :: c_occ(:, :)       !! (nao, n_occ)
      real(dp), intent(inout) :: hess(:, :, :, :)   !! (3, 3, natm, natm), accumulated
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: k_scale
         !! Exact-exchange fraction; one when absent.

      real(dp), allocatable, target :: cc(:, :, :), acc(:, :, :), gamma3(:, :, :)
      real(dp), allocatable, target :: f(:, :, :, :), pm(:, :, :)
      real(dp), allocatable :: dens_b(:, :), d_a(:), d_b(:), ubar_a(:, :), ubar_b(:, :)
      real(dp), allocatable :: ut_a(:, :), ut_b(:, :), gamma2(:, :), work(:, :), gram(:, :)
      real(dp), pointer :: p2(:, :), q2(:, :)
      real(dp) :: kf
      integer :: nao, naux, natm, n_occ, npert, p, q, x, ia, ja, yy, ip, jp

      if (error%has_error()) return
      kf = 1.0_dp
      if (present(k_scale)) kf = k_scale
      nao = ops%nao
      naux = ops%naux
      natm = ops%natm
      npert = 3*natm
      n_occ = size(c_occ, 2)

      allocate (dens_b(nao, nao), d_a(naux), d_b(naux), ubar_a(naux, npert), &
                ubar_b(naux, npert))
      call pic_gemm(c_occ, c_occ, dens_b, transb="T", alpha=2.0_dp)
      call fit_response_vectors(ops, dens, d_a, ubar_a)
      call fit_response_vectors(ops, dens_b, d_b, ubar_b)

      ! `cc_P = c_P C` and `acc_P = A c_P C`.
      allocate (cc(nao, n_occ, naux), acc(nao, n_occ, naux))
      call contract_second(ops%c, c_occ, cc)
      p2(1:nao, 1:n_occ*naux) => cc
      q2(1:nao, 1:n_occ*naux) => acc
      call gemm_over_columns(dens, p2, q2)

      ! ---- second-derivative integrals against two fixed densities ----------
      ! gamma3 = A d_B + B d_A - k/2 (A c_P B + B c_P A),
      ! gamma2 = -(d_A d_B^T + d_B d_A^T)/2 + k/2 Tr(c_P A c_Q B),
      ! and `A c_P B = 2 acc_P C^T`.
      allocate (gamma3(nao, nao, naux), gamma2(naux, naux), work(nao, nao))
      !$omp parallel do default(none) &
      !$omp    shared(gamma3, dens, dens_b, d_a, d_b, acc, c_occ, kf, naux) &
      !$omp    private(p, work) schedule(static)
      do p = 1, naux
         call pic_gemm(acc(:, :, p), c_occ, work, transb="T")
         gamma3(:, :, p) = d_b(p)*dens + d_a(p)*dens_b - kf*(work + transpose(work))
      end do
      !$omp end parallel do
      p2(1:nao*n_occ, 1:naux) => cc
      q2(1:nao*n_occ, 1:naux) => acc
      call gemm_over_inner(p2, q2, gamma2)
      do q = 1, naux
         do p = 1, q
            gamma2(p, q) = 0.5_dp*kf*(gamma2(p, q) + gamma2(q, p)) &
                           - 0.5_dp*(d_a(p)*d_b(q) + d_b(p)*d_a(q))
            gamma2(q, p) = gamma2(p, q)
         end do
      end do
      call hess_3c2c_contract(orb, aux, hess, error, gamma3=gamma3, gamma2=gamma2)
      if (error%has_error()) return
      deallocate (gamma3, gamma2, work, acc)

      ! ---- the product of first derivatives ---------------------------------
      ! Coulomb: ut_A^x . ut_B^y + ut_A^y . ut_B^x. Exchange:
      ! sum_PQ M^-1_PQ [Tr(U^x_P A U^y_Q B) + x <-> y] = 4 gram(x, y).
      allocate (ut_a(naux, npert), ut_b(naux, npert), gram(npert, npert))
      call pic_gemm(ops%mih, ubar_a, ut_a)
      call pic_gemm(ops%mih, ubar_b, ut_b)
      gram = 0.0_dp
      if (kf /= 0.0_dp) then
         ! `f(mu, i, P) = sum_Q cc(mu, i, Q) (nabla P|Q) - sum_nu (mu nu|nabla P) C(nu, i)`
         ! and `pm = (nabla P|Q) M^-1/2`.
         allocate (f(nao, n_occ, naux, 3), pm(naux, naux, 3))
         do x = 1, 3
            call contract_second(ops%ip2(:, :, :, x), c_occ, f(:, :, :, x))
            p2(1:nao*n_occ, 1:naux) => cc
            q2(1:nao*n_occ, 1:naux) => f(:, :, :, x)
            call gemm_over_columns(p2, ops%ip2c(:, :, x), q2, transb="T", alpha=1.0_dp, &
                                   beta=-1.0_dp)
            call gemm_over_columns(ops%ip2c(:, :, x), ops%mih, pm(:, :, x))
         end do
         call cross_exchange_product(ops, dens, c_occ, cc, f, pm, gram)
         deallocate (f, pm)
      end if
      do ja = 1, natm
         do ia = 1, natm
            do yy = 1, 3
               do x = 1, 3
                  ip = 3*(ia - 1) + x
                  jp = 3*(ja - 1) + yy
                  hess(x, yy, ia, ja) = hess(x, yy, ia, ja) &
                                        + dot_product(ut_a(:, ip), ut_b(:, jp)) &
                                        + dot_product(ut_a(:, jp), ut_b(:, ip)) &
                                        - 2.0_dp*kf*gram(ip, jp)
               end do
            end do
         end do
      end do
   end subroutine df_cross_second_derivative

   subroutine cross_exchange_product(ops, dens, c_occ, cc, f, pm, gram)
      !! `gram(x, y) = sum_{i, R} Xt^x(., i, R)^T A Xt^y(., i, R)`
      !!
      !! `Xt^x(mu, i, R) = sum_P (U^x_P C)(mu, i) M^-1/2(P, R)`, built a batch
      !! of `R` at a time for all perturbations:
      !!
      !!     Xt^x = - sum_{lambda in A} t(lambda, mu, R) C(lambda, i)
      !!            - [mu in A] sum_lambda t(mu, lambda, R) C(lambda, i)
      !!            + sum_{P in A} [f(mu, i, P) mih(P, R) + cc(mu, i, P) pm(P, R)]
      !!
      !! with `t = (nabla . .|P) M^-1/2`.
      type(df_operator_t), intent(in), target :: ops
      real(dp), intent(in) :: dens(:, :)                !! `A`
      real(dp), intent(in) :: c_occ(:, :)
      real(dp), intent(in), contiguous, target :: cc(:, :, :)     !! (nao, n_occ, naux)
      real(dp), intent(in), contiguous, target :: f(:, :, :, :)   !! (nao, n_occ, naux, 3)
      real(dp), intent(in) :: pm(:, :, :)               !! (naux, naux, 3)
      real(dp), intent(out) :: gram(:, :)               !! (3 natm, 3 natm)

      real(dp), allocatable, target :: xt(:, :, :, :), axt(:, :, :, :), t(:, :, :)
      real(dp), allocatable :: ca(:, :), sa(:, :), right(:, :), part(:, :), mcol(:, :)
      real(dp), pointer :: p2(:, :), q2(:, :), y2(:, :)
      integer :: nao, n_occ, naux, natm, npert, width, r0, r1, w, ip, ia, x, r
      integer :: m0, na, a0, xa

      nao = ops%nao
      naux = ops%naux
      natm = ops%natm
      npert = 3*natm
      n_occ = size(c_occ, 2)
      ! Two stacks of (nao, n_occ, w, npert) and the `nao^2 w` of `t`.
      width = int(min(real(naux, dp), &
                      max(1.0_dp, memory_budget(PRODUCT_BLIND_BYTES, PRODUCT_SHARE) &
                          /(8.0_dp*real(nao, dp)*(2.0_dp*real(n_occ, dp)*real(npert, dp) &
                                                  + real(nao, dp))))))
      gram = 0.0_dp
      allocate (part(npert, npert))
      do r0 = 1, naux, width
         r1 = min(r0 + width - 1, naux)
         w = r1 - r0 + 1
         allocate (xt(nao, n_occ, w, npert), axt(nao, n_occ, w, npert), t(nao, nao, w), &
                   mcol(naux, w))
         mcol = ops%mih(:, r0:r1)
         do x = 1, 3
            p2(1:nao*nao, 1:naux) => ops%ip1(:, :, :, x)
            q2(1:nao*nao, 1:w) => t
            call gemm_over_columns(p2, mcol, q2)
            !$omp parallel default(none) &
            !$omp    shared(ops, c_occ, cc, f, pm, xt, t, nao, n_occ, natm, r0, r1, w, x) &
            !$omp    private(ia, ip, r, m0, na, a0, xa, ca, sa, right, p2, y2)
            !$omp do schedule(dynamic)
            do ia = 1, natm
               ip = 3*(ia - 1) + x
               na = ops%orb_count(ia)
               m0 = ops%orb_first(ia) + 1
               xa = ops%aux_count(ia)
               a0 = ops%aux_first(ia) + 1
               xt(:, :, :, ip) = 0.0_dp
               if (na > 0) then
                  allocate (ca(na, n_occ), sa(na, nao))
                  ca = c_occ(m0:m0 + na - 1, :)
                  do r = 1, w
                     sa = t(m0:m0 + na - 1, :, r)
                     call pic_gemm(sa, ca, xt(:, :, r, ip), transa="T", alpha=-1.0_dp)
                     call pic_gemm(sa, c_occ, ca, alpha=-1.0_dp)
                     xt(m0:m0 + na - 1, :, r, ip) = xt(m0:m0 + na - 1, :, r, ip) + ca
                     ca = c_occ(m0:m0 + na - 1, :)
                  end do
                  deallocate (ca, sa)
               end if
               if (xa > 0) then
                  allocate (right(xa, w))
                  y2(1:nao*n_occ, 1:w) => xt(:, :, :, ip)
                  right = ops%mih(a0:a0 + xa - 1, r0:r1)
                  p2(1:nao*n_occ, 1:xa) => f(:, :, a0:a0 + xa - 1, x)
                  call pic_gemm(p2, right, y2, alpha=1.0_dp, beta=1.0_dp)
                  right = pm(a0:a0 + xa - 1, r0:r1, x)
                  p2(1:nao*n_occ, 1:xa) => cc(:, :, a0:a0 + xa - 1)
                  call pic_gemm(p2, right, y2, alpha=1.0_dp, beta=1.0_dp)
                  deallocate (right)
               end if
            end do
            !$omp end do
            !$omp end parallel
         end do
         p2(1:nao, 1:n_occ*w*npert) => xt
         q2(1:nao, 1:n_occ*w*npert) => axt
         call gemm_over_columns(dens, p2, q2)
         p2(1:nao*n_occ*w, 1:npert) => axt
         q2(1:nao*n_occ*w, 1:npert) => xt
         call gemm_over_inner(p2, q2, part)
         gram = gram + part
         deallocate (xt, axt, t, mcol)
      end do
   end subroutine cross_exchange_product

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
