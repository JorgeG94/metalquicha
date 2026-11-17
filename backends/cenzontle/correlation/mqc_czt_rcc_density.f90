!! The reduced densities of closed-shell CCSD in spatial orbitals
module mqc_czt_rcc_density
   !! What the amplitudes and the Lambda amplitudes together say about the
   !! electrons: the density a property is the trace of, and the density the
   !! gradient contracts with derivative integrals.
   !!
   !! These are the *correlation* parts. The reference determinant's own
   !! density -- two electrons in each occupied orbital -- is not included, so
   !! the one-particle density here has no trace and is what multiplies the
   !! Fock matrix in the correlation energy.
   !!
   !! **Equations.** PySCF's `cc/ccsd_rdm.py`, term for term, on arrays laid out
   !! as that module lays them out. The amplitude and Lambda equations next
   !! door follow the same package, so one set of conventions runs through all
   !! three.
   !!
   !! Unrelaxed: the orbitals are those of the reference and do not respond.
   !! The response is added by whoever needs it, because it depends on what is
   !! being differentiated.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_tensor_contract, only: contract, transposed
   use mqc_czt_rcc, only: rcc_state_t
   use mqc_czt_rcc_lambda, only: rcc_lambda_t
   use mqc_czt_rcc_triples_terms, only: rcc_triples_terms_t
   implicit none
   private

   public :: rccsd_one_particle_density
   public :: rccsd_fitted_two_particle_density
   public :: rcc_active_fit

contains

   subroutine rccsd_one_particle_density(state, lambda, gamma, triples)
      !! The correlation part of the spin-summed one-particle density
      !!
      !! Over the active orbitals, occupied first and then virtual, and
      !! symmetric. The occupied block traces to minus the virtual block's
      !! trace: correlation moves electrons and makes none.
      type(rcc_state_t), intent(in) :: state
      type(rcc_lambda_t), intent(in) :: lambda
      real(dp), allocatable, intent(out) :: gamma(:, :)   !! (no + nv, no + nv)
      type(rcc_triples_terms_t), intent(in), optional :: triples
         !! Present, the density of CCSD(T). The correction adds to the diagonal
         !! only: it depends on the orbital energies and on nothing else of the
         !! Fock matrix.

      real(dp), allocatable :: doo(:, :), dvv(:, :), dvo(:, :)
      real(dp), allocatable :: xtv(:, :), xt1(:, :), xt2(:, :), theta(:, :, :, :)
      integer :: no, nv, n, p

      no = state%no
      nv = state%nv
      n = no + nv
      associate (t1 => state%t1, t2 => state%t2, l1 => lambda%l1, l2 => lambda%l2)
         allocate (doo(no, no), dvv(nv, nv), xtv(no, no), xt1(no, no), xt2(nv, nv), &
                   source=0.0_dp)
         call contract("ja,ia->ij", t1, l1, doo, -1.0_dp)
         call contract("ia,ib->ab", t1, l1, dvv, 1.0_dp)
         call contract("ie,me->im", t1, l1, xtv, 1.0_dp)
         dvo = transpose(t1)
         call contract("im,ma->ai", xtv, t1, dvo, -1.0_dp)

         theta = 2.0_dp*t2 - transposed(t2, [1, 2, 4, 3])
         call contract("jkab,ikab->ij", theta, l2, doo, -1.0_dp)
         call contract("jica,jicb->ab", theta, l2, dvv, 1.0_dp)
         call contract("mnef,inef->mi", l2, theta, xt1, 1.0_dp)
         call contract("mnaf,mnef->ea", l2, theta, xt2, 1.0_dp)
         call contract("imae,me->ai", theta, l1, dvo, 1.0_dp)
         call contract("mi,ma->ai", xt1, t1, dvo, -1.0_dp)
         call contract("ie,ae->ai", t1, xt2, dvo, -1.0_dp)

         allocate (gamma(n, n))
         gamma(1:no, 1:no) = doo + transpose(doo)
         gamma(1:no, no + 1:n) = l1 + transpose(dvo)
         gamma(no + 1:n, 1:no) = transpose(gamma(1:no, no + 1:n))
         gamma(no + 1:n, no + 1:n) = dvv + transpose(dvv)
      end associate
      if (present(triples)) then
         do p = 1, no
            gamma(p, p) = gamma(p, p) + triples%gamma_occ(p)
         end do
         do p = 1, nv
            gamma(no + p, no + p) = gamma(no + p, no + p) + triples%gamma_vir(p)
         end do
      end if
   end subroutine rccsd_one_particle_density

   subroutine rcc_active_fit(state, b, error)
      !! The fitted three-index integrals over every pair of active orbitals
      !!
      !! `(pq|rs)` is `sum_P b(p,q,P) b(r,s,P)`. Symmetric in `p` and `q`,
      !! occupied orbitals first.
      type(rcc_state_t), intent(in) :: state
      real(dp), allocatable, intent(out) :: b(:, :, :)   !! (no + nv, no + nv, n_aux)
      type(error_t), intent(inout) :: error

      integer :: no, nv, n, n_aux, k

      if (.not. (allocated(state%eris%b_oo) .and. allocated(state%eris%b_ov) &
                 .and. allocated(state%eris%b_vv))) then
         call error%set(ERROR_VALIDATION, "the coupled cluster integrals were not fitted, "// &
                        "so there is no three-index tensor to hand on")
         return
      end if
      no = state%no
      nv = state%nv
      n = no + nv
      n_aux = size(state%eris%b_vv, 2)
      allocate (b(n, n, n_aux))
      do k = 1, n_aux
         b(1:no, 1:no, k) = reshape(state%eris%b_oo(:, k), [no, no])
         b(1:no, no + 1:n, k) = reshape(state%eris%b_ov(:, k), [no, nv])
         b(no + 1:n, 1:no, k) = transpose(b(1:no, no + 1:n, k))
         b(no + 1:n, no + 1:n, k) = reshape(state%eris%b_vv(:, k), [nv, nv])
      end do
   end subroutine rcc_active_fit

   subroutine rccsd_fitted_two_particle_density(state, lambda, b, g3, error, triples)
      !! The two-particle density contracted once with the fit
      !!
      !!     g3(p,q,P) = sum_rs Gamma(p,q,r,s) b(r,s,P)
      !!
      !! for the correlation part `Gamma` of the spin-summed two-particle
      !! density, symmetrised over `p` and `q`, over `r` and `s`, and over the
      !! exchange of the two pairs. The two-electron part of the correlation
      !! energy is half the sum of `g3` times `b`.
      !!
      !! The block of `Gamma` with four virtual indices is never held: its
      !! contraction is taken straight from the amplitudes, one auxiliary
      !! function at a time.
      type(rcc_state_t), intent(in) :: state
      type(rcc_lambda_t), intent(in) :: lambda
      real(dp), intent(in) :: b(:, :, :)
         !! (no + nv, no + nv, n_aux) the fit over the active orbitals, in
         !! whichever factorisation of the metric the caller will go on to use
      real(dp), allocatable, intent(out) :: g3(:, :, :)   !! (no + nv, no + nv, n_aux)
      type(error_t), intent(inout) :: error
      type(rcc_triples_terms_t), intent(in), optional :: triples
         !! Present, the density of CCSD(T)

      real(dp), allocatable :: t(:, :, :), tp(:, :, :)
      real(dp), allocatable :: dovov(:, :, :, :), doooo(:, :, :, :), doovv(:, :, :, :)
      real(dp), allocatable :: dovvo(:, :, :, :), dovvv(:, :, :, :), dooov(:, :, :, :)
      real(dp), allocatable :: tau(:, :, :, :), work(:, :, :, :)
      real(dp), allocatable :: bm(:, :), y(:, :), yp(:, :)
      integer :: no, nv, n, n_aux, k, i, j, a, c

      no = state%no
      nv = state%nv
      n = no + nv
      if (.not. allocated(state%t2) .or. .not. allocated(lambda%l2)) then
         call error%set(ERROR_VALIDATION, "the two-particle density needs converged "// &
                        "amplitudes and Lambda amplitudes")
         return
      end if
      if (size(b, 1) /= n .or. size(b, 2) /= n) then
         call error%set(ERROR_VALIDATION, "the two-particle density was handed a fit "// &
                        "that does not span the active orbitals")
         return
      end if
      n_aux = size(b, 3)
      allocate (t(n, n, n_aux), tp(n, n, n_aux), source=0.0_dp)

      call gamma2_blocks(state, lambda, dovov, doooo, doovv, dovvo, dovvv, dooov)
      if (present(triples)) then
         dovov = dovov + triples%dovov
         dovvv = dovvv + triples%dovvv
         dooov = dooov + triples%dooov
      end if

      ! Each block of the density, and each of the blocks its symmetry makes of
      ! it, contracted with the fit on its second pair and on its first.
      work = dovov + transposed(dovov, [3, 4, 1, 2])
      call accumulate(work, 0, no, 0, no)
      call accumulate(transposed(work, [2, 1, 4, 3]), no, 0, no, 0)
      work = doovv + transposed(doovv, [2, 1, 4, 3])
      call accumulate(work, 0, 0, no, no)
      call accumulate(transposed(work, [3, 4, 1, 2]), no, no, 0, 0)
      work = dovvo + transposed(dovvo, [4, 3, 2, 1])
      call accumulate(work, 0, no, no, 0)
      call accumulate(transposed(work, [2, 1, 4, 3]), no, 0, 0, no)
      work = 2.0_dp*(doooo + transposed(doooo, [2, 1, 4, 3]))
      call accumulate(work, 0, 0, 0, 0)
      call accumulate(dovvv, 0, no, no, no)
      call accumulate(transposed(dovvv, [3, 4, 1, 2]), no, no, 0, no)
      call accumulate(transposed(dovvv, [4, 3, 2, 1]), no, no, no, 0)
      call accumulate(transposed(dovvv, [2, 1, 4, 3]), no, 0, no, no)
      call accumulate(dooov, 0, 0, 0, no)
      call accumulate(transposed(dooov, [3, 4, 1, 2]), 0, no, 0, 0)
      call accumulate(transposed(dooov, [2, 1, 4, 3]), 0, 0, no, 0)
      call accumulate(transposed(dooov, [4, 3, 2, 1]), no, 0, 0, 0)
      deallocate (work, dovov, doooo, doovv, dovvo, dovvv, dooov)

      ! The all-virtual block. With g(a,b,c,d) = sum_ij l2(i,j,a,b) tau(i,j,c,d)
      ! the block is 2 (d + d transposed in both pairs) for
      ! d(p,q,r,s) = g(p,r,q,s) - g(p,r,s,q)/2.
      allocate (tau(no, no, nv, nv))
      do c = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  tau(i, j, a, c) = state%t2(i, j, a, c) + state%t1(i, a)*state%t1(j, c)
               end do
            end do
         end do
      end do
      allocate (work(no, no, nv, nv), y(nv, nv), yp(nv, nv))
      do k = 1, n_aux
         bm = b(no + 1:n, no + 1:n, k)
         work = 0.0_dp
         y = 0.0_dp
         call contract("ijpr,rs->ijps", lambda%l2, bm, work, 1.0_dp)
         call contract("ijps,ijqs->pq", work, tau, y, 1.0_dp)
         call contract("ijps,ijsq->pq", work, tau, y, -0.5_dp)
         work = 0.0_dp
         yp = 0.0_dp
         call contract("ijpr,pq->ijrq", lambda%l2, bm, work, 1.0_dp)
         call contract("ijrq,ijqs->rs", work, tau, yp, 1.0_dp)
         call contract("ijrq,ijsq->rs", work, tau, yp, -0.5_dp)
         t(no + 1:n, no + 1:n, k) = t(no + 1:n, no + 1:n, k) + 2.0_dp*(y + transpose(y))
         tp(no + 1:n, no + 1:n, k) = tp(no + 1:n, no + 1:n, k) + 2.0_dp*(yp + transpose(yp))
      end do

      allocate (g3(n, n, n_aux))
      do k = 1, n_aux
         g3(:, :, k) = 0.25_dp*(t(:, :, k) + tp(:, :, k) &
                                + transpose(t(:, :, k)) + transpose(tp(:, :, k)))
      end do

   contains

      subroutine accumulate(d, w0, x0, y0, z0)
         !! Add one block `d(w,x,y,z)` of the density to both contractions
         real(dp), intent(in) :: d(:, :, :, :)
         integer, intent(in) :: w0, x0, y0, z0
            !! Where each index starts among the active orbitals, less one:
            !! zero for an occupied index and `no` for a virtual one

         real(dp), allocatable :: dm(:, :), bs(:, :), out(:, :)
         integer :: nw, nx, ny, nz

         nw = size(d, 1)
         nx = size(d, 2)
         ny = size(d, 3)
         nz = size(d, 4)
         dm = reshape(d, [nw*nx, ny*nz])

         bs = reshape(b(y0 + 1:y0 + ny, z0 + 1:z0 + nz, :), [ny*nz, n_aux])
         allocate (out(nw*nx, n_aux))
         call pic_gemm(dm, bs, out)
         t(w0 + 1:w0 + nw, x0 + 1:x0 + nx, :) = t(w0 + 1:w0 + nw, x0 + 1:x0 + nx, :) &
                                                + reshape(out, [nw, nx, n_aux])
         deallocate (out)

         bs = reshape(b(w0 + 1:w0 + nw, x0 + 1:x0 + nx, :), [nw*nx, n_aux])
         allocate (out(ny*nz, n_aux))
         call pic_gemm(dm, bs, out, transa="T")
         tp(y0 + 1:y0 + ny, z0 + 1:z0 + nz, :) = tp(y0 + 1:y0 + ny, z0 + 1:z0 + nz, :) &
                                                 + reshape(out, [ny, nz, n_aux])
      end subroutine accumulate

   end subroutine rccsd_fitted_two_particle_density

   subroutine gamma2_blocks(state, lambda, dovov, doooo, doovv, dovvo, dovvv, dooov)
      !! `_gamma2_outcore` of the reference, without its all-virtual block
      !!
      !! Each block is named and laid out as the reference names and lays it
      !! out, and is the part of the density the reference stores, before the
      !! symmetry that fills in the rest.
      type(rcc_state_t), intent(in) :: state
      type(rcc_lambda_t), intent(in) :: lambda
      real(dp), allocatable, intent(out) :: dovov(:, :, :, :), doooo(:, :, :, :)
      real(dp), allocatable, intent(out) :: doovv(:, :, :, :), dovvo(:, :, :, :)
      real(dp), allocatable, intent(out) :: dovvv(:, :, :, :), dooov(:, :, :, :)

      real(dp), allocatable :: pv_x(:, :, :, :), pv_d(:, :, :, :)
      ! The reference's pvOOv and pvoOV, both (nv, no, no, nv)
      real(dp), allocatable :: theta(:, :, :, :), tau(:, :, :, :), tau_half(:, :, :, :)
      real(dp), allocatable :: goooo(:, :, :, :), gooov(:, :, :, :), goovv(:, :, :, :)
      real(dp), allocatable :: g_x(:, :, :, :), g_d(:, :, :, :), gvovv(:, :, :, :)
      ! `g_x` and `g_d` are the reference's gOvvO and govVO
      real(dp), allocatable :: work(:, :, :, :)
      real(dp), allocatable :: moo(:, :), mvv(:, :), mia(:, :), mab(:, :), mij(:, :)
      real(dp), allocatable :: tmp_vv(:, :)
      integer :: no, nv, i, j, l, a, b, d

      no = state%no
      nv = state%nv
      associate (t1 => state%t1, t2 => state%t2, l1 => lambda%l1, l2 => lambda%l2)

         allocate (pv_x(nv, no, no, nv), pv_d(nv, no, no, nv), source=0.0_dp)
         call contract("ikca,jkcb->aijb", l2, t2, pv_x, 1.0_dp)
         theta = 2.0_dp*t2 - transposed(t2, [1, 2, 4, 3])
         call contract("ikca,jkbc->aijb", l2, t2, pv_d, -1.0_dp)
         call contract("ikac,jkbc->aijb", l2, theta, pv_d, 1.0_dp)

         ! 'dljd->jl' and 'blld->db' of twice the first plus the second.
         allocate (moo(no, no), mvv(nv, nv), source=0.0_dp)
         do d = 1, nv
            do j = 1, no
               do l = 1, no
                  moo(j, l) = moo(j, l) + 2.0_dp*pv_x(d, l, j, d) + pv_d(d, l, j, d)
               end do
            end do
         end do
         do d = 1, nv
            do l = 1, no
               do b = 1, nv
                  mvv(d, b) = mvv(d, b) + 2.0_dp*pv_x(b, l, l, d) + pv_d(b, l, l, d)
               end do
            end do
         end do

         allocate (gooov(no, no, no, nv), source=0.0_dp)
         call contract("kc,cija->jkia", t1, pv_x, gooov, 1.0_dp)
         call contract("jc,cika->jkia", t1, pv_d, gooov, -1.0_dp)

         allocate (mia(no, nv), mab(nv, nv), source=0.0_dp)
         call contract("kc,ikac->ia", l1, t2, mia, 2.0_dp)
         call contract("kc,ikca->ia", l1, t2, mia, -1.0_dp)
         call contract("kc,kb->cb", l1, t1, mab, 1.0_dp)
         mij = 0.5_dp*moo
         call contract("kc,jc->jk", l1, t1, mij, 1.0_dp)

         allocate (tau(no, no, nv, nv))
         do b = 1, nv
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     tau(i, j, a, b) = t2(i, j, a, b) + t1(i, a)*t1(j, b)
                  end do
               end do
            end do
         end do
         allocate (goooo(no, no, no, no), source=0.0_dp)
         call contract("ijab,klab->ijkl", tau, l2, goooo, 0.5_dp)
         doooo = 2.0_dp*transposed(goooo, [1, 3, 2, 4]) - transposed(goooo, [1, 4, 2, 3])

         call contract("ji,ka->jkia", moo, t1, gooov, -0.5_dp)
         call contract("la,jkil->jkia", t1, goooo, gooov, 2.0_dp)
         call contract("ib,jkba->jkia", l1, tau, gooov, -1.0_dp)
         call contract("jkba,ib->jkia", l2, t1, gooov, -1.0_dp)
         dooov = 2.0_dp*transposed(gooov, [1, 3, 2, 4]) - transposed(gooov, [2, 3, 1, 4])
         deallocate (gooov)

         goovv = 0.5_dp*tau + 0.5_dp*l2
         call contract("ia,jb->ijab", mia, t1, goovv, 1.0_dp)
         call contract("ijkl,klab->ijab", goooo, tau, goovv, 1.0_dp)
         call contract("jk,ikab->ijab", mij, tau, goovv, -1.0_dp)
         call contract("cb,ijac->ijab", mab, t2, goovv, -1.0_dp)
         call contract("bd,ijad->ijab", mvv, tau, goovv, -0.5_dp)

         ! 'kiac,jc,kb->iabj' and 'ikac,jc,kb->iabj', each in two steps.
         g_x = transposed(pv_x, [2, 1, 4, 3])
         g_d = transposed(pv_d, [2, 1, 4, 3])
         allocate (work(no, no, nv, no), source=0.0_dp)
         call contract("kiac,jc->kiaj", l2, t1, work, 1.0_dp)
         call contract("kiaj,kb->iabj", work, t1, g_x, 1.0_dp)
         work = 0.0_dp
         call contract("ikac,jc->ikaj", l2, t1, work, 1.0_dp)
         call contract("ikaj,kb->iabj", work, t1, g_d, -1.0_dp)
         deallocate (work)
         call contract("ia,jb->iabj", l1, t1, g_d, 1.0_dp)
         dovvo = 2.0_dp*g_d + g_x
         doovv = transposed(-2.0_dp*g_x - g_d, [4, 1, 2, 3])
         deallocate (g_x, g_d)

         tau_half = tau - 0.5_dp*t2
         call contract("dlib,jlda->ijab", pv_x, tau_half, goovv, 1.0_dp)
         call contract("dlia,jldb->ijab", pv_d, tau_half, goovv, -1.0_dp)
         call contract("dlia,jlbd->ijab", pv_d + 0.5_dp*pv_x, t2, goovv, 1.0_dp)
         dovov = 2.0_dp*transposed(goovv, [1, 3, 2, 4]) - transposed(goovv, [2, 3, 1, 4])
         deallocate (goovv, goooo, tau_half)

         ! The reference contracts its all-virtual g with t1 here. Summing
         ! over the second virtual index first leaves nothing of fourth order
         ! in the virtuals to hold.
         allocate (gvovv(nv, no, nv, nv), source=0.0_dp)
         allocate (work(no, no, nv, no), source=0.0_dp)
         call contract("jkad,id->jkai", l2, t1, work, 1.0_dp)
         call contract("jkai,jkbc->aibc", work, tau, gvovv, -1.0_dp)
         deallocate (work)
         call contract("akic,kb->aibc", pv_d, t1, gvovv, 1.0_dp)
         call contract("akib,kc->aibc", pv_x, t1, gvovv, -1.0_dp)
         call contract("ja,jibc->aibc", l1, t2, gvovv, 1.0_dp)
         allocate (tmp_vv(nv, nv), source=0.0_dp)
         call contract("ja,jb->ab", l1, t1, tmp_vv, 1.0_dp)
         call contract("ab,ic->aibc", tmp_vv, t1, gvovv, 1.0_dp)
         call contract("ba,ic->aibc", mvv, t1, gvovv, 0.5_dp)
         call contract("ja,jibc->aibc", t1, l2, gvovv, 1.0_dp)
         dovvv = 2.0_dp*transposed(gvovv, [2, 4, 1, 3]) - transposed(gvovv, [2, 3, 1, 4])
      end associate
   end subroutine gamma2_blocks

end module mqc_czt_rcc_density
