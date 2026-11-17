!! What the perturbative triples add to the CCSD Lambda equations and densities
module mqc_czt_rcc_triples_terms
   !! The (T) correction is a function of the converged amplitudes, of the
   !! integrals and, through its denominators, of the orbital energies. Its
   !! gradient therefore reaches the rest of the machinery in three places:
   !!
   !! * the Lambda equations gain a source in each of their two residuals, the
   !!   correction's derivative with respect to `t1` and `t2`;
   !! * the one-particle density gains a diagonal, the derivative with respect
   !!   to the orbital energies; and
   !! * three blocks of the two-particle density gain the derivative with
   !!   respect to the integrals.
   !!
   !! **Equations.** PySCF's `cc/ccsd_t_lambda_slow.py` and
   !! `cc/ccsd_t_rdm_slow.py`, term for term, with `f_ov = 0`.
   !!
   !! **Cost.** One set of three virtuals at a time, as the energy is
   !! evaluated: the six-index intermediates the reference holds whole are
   !! built and consumed an `n_occ^3` block at a time, so nothing here is
   !! larger than the four-index results. The work is `n_occ^3 n_vir^4`,
   !! threaded.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use omp_lib, only: omp_lock_kind, omp_init_lock, omp_set_lock, omp_unset_lock, &
                      omp_destroy_lock
   use mqc_tensor_contract, only: transposed
   use mqc_czt_rcc, only: rcc_state_t
   implicit none
   private

   public :: rcc_triples_terms_t
   public :: rccsd_triples_terms

   integer, parameter :: N_ORDERINGS = 6
      !! Orderings of three virtuals
   integer, parameter :: ORDERINGS(3, N_ORDERINGS) = reshape( &
                         [1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1], &
                         [3, N_ORDERINGS])
      !! Each column is one: which member of the set is first, second and third

   type :: rcc_triples_terms_t
      !! The (T) contributions, each in the reference's layout
      real(dp), allocatable :: l1_source(:, :)          !! (no, nv), added to every `l1` update
      real(dp), allocatable :: l2_source(:, :, :, :)    !! (no, no, nv, nv)
      real(dp), allocatable :: gamma_occ(:)
         !! (no) added to the diagonal of the one-particle density
      real(dp), allocatable :: gamma_vir(:)             !! (nv)
      real(dp), allocatable :: dovov(:, :, :, :)        !! Added to the block of that name
      real(dp), allocatable :: dooov(:, :, :, :)
      real(dp), allocatable :: dovvv(:, :, :, :)
   end type rcc_triples_terms_t

contains

   subroutine rccsd_triples_terms(state, terms)
      !! Everything (T) adds, from converged CCSD amplitudes
      type(rcc_state_t), intent(in) :: state
      type(rcc_triples_terms_t), intent(out) :: terms

      real(dp), allocatable :: tvoo(:, :, :, :), tjk(:, :, :, :), ovvv_p(:, :, :, :)
      real(dp), allocatable :: ovoo_p(:, :, :, :), ovoo_q(:, :, :, :), joovv(:, :, :, :)
      real(dp), allocatable :: w(:, :, :), v(:, :, :), wp(:, :, :), vp(:, :, :)
      real(dp), allocatable :: rw(:, :, :), x(:, :, :), m(:, :, :), blk(:, :, :)
      real(dp), allocatable :: out_ov(:, :), out_oo(:, :), out_ooo(:, :, :)
      real(dp), allocatable :: j_add(:, :, :), dov_add(:, :), l1_add(:), gamma_occ(:)
      integer(omp_lock_kind), allocatable :: locks(:)
      integer :: no, nv, i, j, k, a, b, c, p, q, virt(3), now(3), inverse(3)
      real(dp) :: acc, d3, gv_add
      logical :: seen

      no = state%no
      nv = state%nv

      ! The integrals and amplitudes laid out so that what one virtual triple
      ! reads is contiguous.
      tvoo = transposed(state%t2, [4, 1, 2, 3])            ! (f, k, j, z) = t2(k, j, z, f)
      tjk = transposed(state%t2, [2, 1, 4, 3])             ! (j, k, f, c) = t2(k, j, c, f)
      ovvv_p = transposed(state%eris%ovvv, [1, 3, 2, 4])   ! (i, f, x, y) = (ix|fy)
      ovoo_p = transposed(state%eris%ooov, [3, 1, 2, 4])   ! (i, j, m, x) = (ix|jm)
      ovoo_q = transposed(state%eris%ooov, [1, 3, 2, 4])   ! (m, n, j, c) = (nc|mj)

      allocate (terms%l1_source(no, nv), terms%gamma_vir(nv), source=0.0_dp)
      allocate (terms%dovov(no, nv, no, nv), terms%dooov(no, no, no, nv), source=0.0_dp)
      allocate (terms%dovvv(no, nv, nv, nv), joovv(no, no, nv, nv), source=0.0_dp)
      allocate (gamma_occ(no), source=0.0_dp)
      allocate (locks(nv))
      do a = 1, nv
         call omp_init_lock(locks(a))
      end do

      ! One set of three virtuals at a time. The symmetrised intermediate is
      ! built once for the set, at six products, and each ordering of the set
      ! is that block with its occupied indices permuted to match. Every
      ! result but one is indexed by an ordering's first virtual, which any
      ! thread may be holding, so each ordering computes into its own arrays
      ! and adds them under that virtual's lock. The occupied diagonal is the
      ! exception and is reduced.
      !$omp parallel default(none) &
      !$omp    shared(state, terms, joovv, tvoo, tjk, ovvv_p, ovoo_p, ovoo_q, no, nv, locks) &
      !$omp    private(a, b, c, i, j, k, p, q, virt, now, inverse, seen, d3, acc, gv_add) &
      !$omp    private(w, v, wp, vp, rw, x, m, blk, out_ov, out_oo, out_ooo) &
      !$omp    private(j_add, dov_add, l1_add) &
      !$omp    reduction(+:gamma_occ)
      allocate (w(no, no, no), v(no, no, no), wp(no, no, no), vp(no, no, no))
      allocate (rw(no, no, no), x(no, no, no), m(no, no, no), blk(no, no, no))
      allocate (out_ov(no, nv), out_oo(no, no), out_ooo(no, no, no))
      allocate (j_add(no, no, nv), dov_add(no, no), l1_add(no))
      ! Descending, because the first index of a set carries work in
      ! proportion to its square and the large ones should be handed out
      ! first.
      !$omp do schedule(dynamic)
      do a = nv, 1, -1
         do b = 1, a
            do c = 1, b
               ! The reference's p6: the six simultaneous permutations of the
               ! three (occupied, virtual) pairs.
               w = 0.0_dp
               call bare_w(a, b, c, blk)
               call add_permuted(blk, [1, 2, 3], w)
               call bare_w(a, c, b, blk)
               call add_permuted(blk, [1, 3, 2], w)
               call bare_w(b, a, c, blk)
               call add_permuted(blk, [2, 1, 3], w)
               call bare_w(c, a, b, blk)
               call add_permuted(blk, [3, 1, 2], w)
               call bare_w(b, c, a, blk)
               call add_permuted(blk, [2, 3, 1], w)
               call bare_w(c, b, a, blk)
               call add_permuted(blk, [3, 2, 1], w)
               do k = 1, no
                  do j = 1, no
                     do i = 1, no
                        d3 = state%eps_o(i) + state%eps_o(j) + state%eps_o(k) &
                             - state%eps_v(a) - state%eps_v(b) - state%eps_v(c)
                        w(i, j, k) = w(i, j, k)/d3
                        ! 'iajb,kc->ijkabc' through the same six permutations.
                        v(i, j, k) = (state%eris%ovov(i, a, j, b)*state%t1(k, c) &
                                      + state%eris%ovov(i, a, k, c)*state%t1(j, b) &
                                      + state%eris%ovov(j, b, i, a)*state%t1(k, c) &
                                      + state%eris%ovov(k, c, i, a)*state%t1(j, b) &
                                      + state%eris%ovov(j, b, k, c)*state%t1(i, a) &
                                      + state%eris%ovov(k, c, j, b)*state%t1(i, a))/d3
                     end do
                  end do
               end do

               virt = [a, b, c]
               do p = 1, N_ORDERINGS
                  ! This ordering, unless an earlier one of a set with a
                  ! repeated virtual was the same triple.
                  now = virt(ORDERINGS(:, p))
                  seen = .false.
                  do q = 1, p - 1
                     if (all(virt(ORDERINGS(:, q)) == now)) seen = .true.
                  end do
                  if (seen) cycle
                  ! The occupied index that goes with each virtual goes with it.
                  do q = 1, 3
                     inverse(ORDERINGS(q, p)) = q
                  end do
                  wp = 0.0_dp
                  vp = 0.0_dp
                  call add_permuted(w, inverse, wp)
                  call add_permuted(v, inverse, vp)

                  ! ---- with r6(w) ------------------------------------------
                  call spin_sum(wp, rw)
                  ! l1 source: 'jbkc,ijkabc->ia', halved; and
                  ! 'kc,ijkabc->iajb', halved, into dovov.
                  do i = 1, no
                     acc = 0.0_dp
                     do k = 1, no
                        do j = 1, no
                           acc = acc + state%eris%ovov(j, now(2), k, now(3))*rw(i, j, k)
                        end do
                     end do
                     l1_add(i) = 0.5_dp*acc
                  end do
                  do j = 1, no
                     do i = 1, no
                        acc = 0.0_dp
                        do k = 1, no
                           acc = acc + state%t1(k, now(3))*rw(i, j, k)
                        end do
                        dov_add(i, j) = 0.5_dp*acc
                     end do
                  end do
                  ! The one-particle diagonals, of (w + v/2) with r6(w).
                  gv_add = 0.0_dp
                  do k = 1, no
                     do j = 1, no
                        do i = 1, no
                           acc = (wp(i, j, k) + 0.5_dp*vp(i, j, k))*rw(i, j, k)
                           gamma_occ(i) = gamma_occ(i) - acc
                           gv_add = gv_add + acc
                        end do
                     end do
                  end do

                  ! ---- with 2 w + v/2 --------------------------------------
                  x = 2.0_dp*wp + 0.5_dp*vp
                  ! l2 source. 'kfbe,ijkaef->ijab': this ordering is (a, e, f),
                  ! so it reaches every b through (kf|be). Then
                  ! 'ncmj,imnabc->ijab'.
                  call spin_pair(x, m)
                  call product_of(m, no*no, no, ovvv_p(:, :, now(3), now(2)), nv, j_add)
                  call pic_gemm(reshape(m, [no, no*no]), &
                                reshape(ovoo_q(:, :, :, now(3)), [no*no, no]), out_oo)

                  ! dooov -= 'mkbc,ijkabc->jmia' and dovvv += 'kjcf,ijkabc->iafb'.
                  call spin_sum(x, rw)
                  call fold_ooo(rw, state%t2(:, :, now(2), now(3)), out_ooo)
                  call pic_gemm(reshape(rw, [no, no*no]), &
                                reshape(tjk(:, :, :, now(3)), [no*no, nv]), out_ov)

                  call omp_set_lock(locks(now(1)))
                  terms%l1_source(:, now(1)) = terms%l1_source(:, now(1)) + l1_add
                  terms%dovov(:, now(1), :, now(2)) = terms%dovov(:, now(1), :, now(2)) + dov_add
                  terms%gamma_vir(now(1)) = terms%gamma_vir(now(1)) + gv_add
                  joovv(:, :, now(1), :) = joovv(:, :, now(1), :) + j_add
                  joovv(:, :, now(1), now(2)) = joovv(:, :, now(1), now(2)) - out_oo
                  do k = 1, no
                     do j = 1, no
                        do i = 1, no
                           terms%dooov(j, k, i, now(1)) = terms%dooov(j, k, i, now(1)) &
                                                          - out_ooo(i, j, k)
                        end do
                     end do
                  end do
                  terms%dovvv(:, now(1), :, now(2)) = terms%dovvv(:, now(1), :, now(2)) + out_ov
                  call omp_unset_lock(locks(now(1)))
               end do
            end do
         end do
      end do
      !$omp end do
      deallocate (w, v, wp, vp, rw, x, m, blk, out_ov, out_oo, out_ooo, j_add, dov_add, l1_add)
      !$omp end parallel
      call move_alloc(gamma_occ, terms%gamma_occ)
      do a = 1, nv
         call omp_destroy_lock(locks(a))
      end do

      do a = 1, nv
         do i = 1, no
            terms%l1_source(i, a) = terms%l1_source(i, a)/(state%eps_o(i) - state%eps_v(a))
         end do
      end do
      joovv = joovv + transposed(joovv, [2, 1, 4, 3])
      allocate (terms%l2_source(no, no, nv, nv))
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  terms%l2_source(i, j, a, b) = joovv(i, j, a, b) &
                                                /(state%eps_o(i) - state%eps_v(a) + state%eps_o(j) - state%eps_v(b))
               end do
            end do
         end do
      end do

   contains

      subroutine bare_w(p, q, r, out)
         !! `out(i,j,k)` = 'iafb,kjcf->ijkabc' - 'iajm,mkbc->ijkabc' at the
         !! virtual triple `(p, q, r)`, before the denominator
         integer, intent(in) :: p, q, r
         real(dp), intent(out) :: out(:, :, :)

         real(dp) :: first(no, no, no)
         integer :: ii, jj, kk

         ! sum_f (ip|fq) t2(k,j,r,f), which comes out as (i, k, j).
         call first_term(ovvv_p(:, :, p, q), tvoo(:, :, :, r), first)
         ! sum_m (ip|jm) t2(m,k,q,r), already (i, j, k).
         call second_term(ovoo_p(:, :, :, p), state%t2(:, :, q, r), out)
         do kk = 1, no
            do jj = 1, no
               do ii = 1, no
                  out(ii, jj, kk) = first(ii, kk, jj) - out(ii, jj, kk)
               end do
            end do
         end do
      end subroutine bare_w

   end subroutine rccsd_triples_terms

   subroutine first_term(ovv, tvo, out)
      !! `out(i, k, j) = sum_f ovv(i, f) tvo(f, k, j)`
      real(dp), intent(in) :: ovv(:, :), tvo(:, :, :)
      real(dp), intent(out) :: out(:, :, :)

      real(dp), allocatable :: flat(:, :)

      allocate (flat(size(out, 1), size(out, 2)*size(out, 3)))
      call pic_gemm(ovv, reshape(tvo, [size(tvo, 1), size(tvo, 2)*size(tvo, 3)]), flat)
      out = reshape(flat, shape(out))
   end subroutine first_term

   subroutine second_term(ooo, t, out)
      !! `out(i, j, k) = sum_m ooo(i, j, m) t(m, k)`
      real(dp), intent(in) :: ooo(:, :, :), t(:, :)
      real(dp), intent(out) :: out(:, :, :)

      real(dp), allocatable :: flat(:, :)

      allocate (flat(size(out, 1)*size(out, 2), size(out, 3)))
      call pic_gemm(reshape(ooo, [size(ooo, 1)*size(ooo, 2), size(ooo, 3)]), t, flat)
      out = reshape(flat, shape(out))
   end subroutine second_term

   subroutine fold_ooo(r, t, out)
      !! `out(i, j, m) = sum_k r(i, j, k) t(m, k)`
      real(dp), intent(in) :: r(:, :, :), t(:, :)
      real(dp), intent(out) :: out(:, :, :)

      real(dp), allocatable :: flat(:, :)

      allocate (flat(size(out, 1)*size(out, 2), size(out, 3)))
      call pic_gemm(reshape(r, [size(r, 1)*size(r, 2), size(r, 3)]), t, flat, transb="T")
      out = reshape(flat, shape(out))
   end subroutine fold_ooo

   subroutine product_of(m, n_rows, n_inner, right, n_cols, out)
      !! `out(i, j, b) = sum_k m(i, j, k) right(k, b)`
      integer, intent(in) :: n_rows, n_inner, n_cols
      real(dp), intent(in) :: m(n_rows, n_inner), right(:, :)
      real(dp), intent(out) :: out(n_rows, n_cols)

      call pic_gemm(m, right, out)
   end subroutine product_of

   pure subroutine add_permuted(blk, order, acc)
      !! `acc(i, j, k) += blk` at the occupied indices taken in `order`: with
      !! `order = [2, 3, 1]`, `blk(j, k, i)`
      real(dp), intent(in) :: blk(:, :, :)
      integer, intent(in) :: order(3)
      real(dp), intent(inout) :: acc(:, :, :)

      integer :: idx(3), i, j, k, no

      no = size(blk, 1)
      do k = 1, no
         do j = 1, no
            do i = 1, no
               idx = [i, j, k]
               acc(i, j, k) = acc(i, j, k) + blk(idx(order(1)), idx(order(2)), idx(order(3)))
            end do
         end do
      end do
   end subroutine add_permuted

   pure subroutine spin_sum(x, r)
      !! The reference's `r6` within one virtual triple, where each
      !! permutation of the virtuals is a permutation of the occupied indices
      real(dp), intent(in) :: x(:, :, :)
      real(dp), intent(out) :: r(:, :, :)

      integer :: i, j, k, no

      no = size(x, 1)
      do k = 1, no
         do j = 1, no
            do i = 1, no
               r(i, j, k) = 4.0_dp*x(i, j, k) + x(k, i, j) + x(j, k, i) &
                            - 2.0_dp*x(k, j, i) - 2.0_dp*x(i, k, j) - 2.0_dp*x(j, i, k)
            end do
         end do
      end do
   end subroutine spin_sum

   pure subroutine spin_pair(x, r)
      !! The reference's `as_r6` within one virtual triple
      real(dp), intent(in) :: x(:, :, :)
      real(dp), intent(out) :: r(:, :, :)

      integer :: i, j, k, no

      no = size(x, 1)
      do k = 1, no
         do j = 1, no
            do i = 1, no
               r(i, j, k) = 2.0_dp*x(i, j, k) - x(k, j, i) - x(i, k, j)
            end do
         end do
      end do
   end subroutine spin_pair

end module mqc_czt_rcc_triples_terms
