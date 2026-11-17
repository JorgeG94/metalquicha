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
   !! **Cost.** As the reference does it: the two six-index intermediates are
   !! held whole, `n_occ^3 n_vir^3` each, and every contraction is a loop. That
   !! is the form that can be checked against the reference line by line, and
   !! it is not the form a large molecule can afford.
   use pic_types, only: dp
   use mqc_tensor_contract, only: transposed
   use mqc_czt_rcc, only: rcc_state_t
   implicit none
   private

   public :: rcc_triples_terms_t
   public :: rccsd_triples_terms

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

      real(dp), allocatable :: w(:, :, :, :, :, :), v(:, :, :, :, :, :)
      real(dp), allocatable :: rw(:, :, :, :, :, :), m(:, :, :, :, :, :)
      real(dp), allocatable :: ovoo(:, :, :, :), joovv(:, :, :, :)
      integer :: no, nv, i, j, k, l, a, b, c, e, f
      real(dp) :: acc, d3

      no = state%no
      nv = state%nv
      ovoo = transposed(state%eris%ooov, [3, 4, 1, 2])   ! (ia|jk) from (jk|ia)

      associate (t1 => state%t1, t2 => state%t2, ovov => state%eris%ovov, &
                 ovvv => state%eris%ovvv, eo => state%eps_o, ev => state%eps_v)

         ! w = 'iafb,kjcf->ijkabc' - 'iajm,mkbc->ijkabc', and v = 'iajb,kc->ijkabc',
         ! both over the three-particle denominator.
         allocate (w(no, no, no, nv, nv, nv), v(no, no, no, nv, nv, nv))
         ! `default(shared)` in this routine's regions: the amplitudes and the
         ! integrals are associate names here, and a data-sharing clause may
         ! not list one.
         !$omp parallel do collapse(3) default(shared) &
         !$omp    private(a, b, c, i, j, k, l, f, acc, d3)
         do c = 1, nv
            do b = 1, nv
               do a = 1, nv
                  do k = 1, no
                     do j = 1, no
                        do i = 1, no
                           d3 = eo(i) - ev(a) + eo(j) - ev(b) + eo(k) - ev(c)
                           acc = 0.0_dp
                           do f = 1, nv
                              acc = acc + ovvv(i, a, f, b)*t2(k, j, c, f)
                           end do
                           do l = 1, no
                              acc = acc - ovoo(i, a, j, l)*t2(l, k, b, c)
                           end do
                           w(i, j, k, a, b, c) = acc/d3
                           v(i, j, k, a, b, c) = ovov(i, a, j, b)*t1(k, c)/d3
                        end do
                     end do
                  end do
               end do
            end do
         end do
         !$omp end parallel do
         call symmetrise_p6(w)
         call symmetrise_p6(v)

         ! ---- the Lambda sources ---------------------------------------------
         call spin_sum_r6(w, rw)
         allocate (terms%l1_source(no, nv))
         do a = 1, nv
            do i = 1, no
               acc = 0.0_dp
               do c = 1, nv
                  do b = 1, nv
                     do k = 1, no
                        do j = 1, no
                           acc = acc + ovov(j, b, k, c)*rw(i, j, k, a, b, c)
                        end do
                     end do
                  end do
               end do
               terms%l1_source(i, a) = 0.5_dp*acc/(eo(i) - ev(a))
            end do
         end do

         ! 'kfbe,ijkaef->ijab' - 'ncmj,imnabc->ijab' with m = as_r6(2 w + v/2).
         call spin_pair_r6(2.0_dp*w + 0.5_dp*v, m)
         allocate (joovv(no, no, nv, nv))
         !$omp parallel do collapse(2) default(shared) &
         !$omp    private(a, b, i, j, k, l, c, e, f, acc)
         do b = 1, nv
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     acc = 0.0_dp
                     do f = 1, nv
                        do e = 1, nv
                           do k = 1, no
                              acc = acc + ovvv(k, f, b, e)*m(i, j, k, a, e, f)
                           end do
                        end do
                     end do
                     do c = 1, nv
                        do l = 1, no
                           do k = 1, no
                              ! n -> l, m -> k of the reference's string
                              acc = acc - ovoo(l, c, k, j)*m(i, k, l, a, b, c)
                           end do
                        end do
                     end do
                     joovv(i, j, a, b) = acc
                  end do
               end do
            end do
         end do
         !$omp end parallel do
         deallocate (m)
         joovv = joovv + transposed(joovv, [2, 1, 4, 3])
         allocate (terms%l2_source(no, no, nv, nv))
         do b = 1, nv
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     terms%l2_source(i, j, a, b) = joovv(i, j, a, b) &
                                                   /(eo(i) - ev(a) + eo(j) - ev(b))
                  end do
               end do
            end do
         end do
         deallocate (joovv)

         ! ---- the one-particle density ---------------------------------------
         ! goo = 'iklabc,jklabc->ij' and gvv = 'ijkacd,ijkbcd->ab' of w + v/2
         ! with r6(w). Only their diagonals are a derivative of the energy; the
         ! gradient settles the rest of both blocks from the Lagrangian.
         allocate (terms%gamma_occ(no), terms%gamma_vir(nv), source=0.0_dp)
         do c = 1, nv
            do b = 1, nv
               do a = 1, nv
                  do k = 1, no
                     do j = 1, no
                        do i = 1, no
                           acc = (w(i, j, k, a, b, c) + 0.5_dp*v(i, j, k, a, b, c)) &
                                 *rw(i, j, k, a, b, c)
                           terms%gamma_occ(i) = terms%gamma_occ(i) - acc
                           terms%gamma_vir(a) = terms%gamma_vir(a) + acc
                        end do
                     end do
                  end do
               end do
            end do
         end do

         ! ---- the two-particle density ---------------------------------------
         ! dovov += 'kc,ijkabc->iajb' of t1 with r6(w), halved.
         allocate (terms%dovov(no, nv, no, nv))
         do b = 1, nv
            do j = 1, no
               do a = 1, nv
                  do i = 1, no
                     acc = 0.0_dp
                     do c = 1, nv
                        do k = 1, no
                           acc = acc + t1(k, c)*rw(i, j, k, a, b, c)
                        end do
                     end do
                     terms%dovov(i, a, j, b) = 0.5_dp*acc
                  end do
               end do
            end do
         end do

         ! dooov -= 'mkbc,ijkabc->jmia' and dovvv += 'kjcf,ijkabc->iafb', both
         ! of t2 with r6(2 w + v/2).
         call spin_sum_r6(2.0_dp*w + 0.5_dp*v, rw)
         deallocate (w, v)
         allocate (terms%dooov(no, no, no, nv), terms%dovvv(no, nv, nv, nv))
         !$omp parallel do collapse(2) default(shared) &
         !$omp    private(a, i, j, l, k, b, c, acc)
         do a = 1, nv
            do i = 1, no
               do l = 1, no
                  do j = 1, no
                     acc = 0.0_dp
                     do c = 1, nv
                        do b = 1, nv
                           do k = 1, no
                              acc = acc + t2(l, k, b, c)*rw(i, j, k, a, b, c)
                           end do
                        end do
                     end do
                     terms%dooov(j, l, i, a) = -acc
                  end do
               end do
            end do
         end do
         !$omp end parallel do
         !$omp parallel do collapse(2) default(shared) &
         !$omp    private(a, b, f, i, j, k, c, acc)
         do b = 1, nv
            do f = 1, nv
               do a = 1, nv
                  do i = 1, no
                     acc = 0.0_dp
                     do c = 1, nv
                        do k = 1, no
                           do j = 1, no
                              acc = acc + t2(k, j, c, f)*rw(i, j, k, a, b, c)
                           end do
                        end do
                     end do
                     terms%dovvv(i, a, f, b) = acc
                  end do
               end do
            end do
         end do
         !$omp end parallel do
      end associate
   end subroutine rccsd_triples_terms

   subroutine symmetrise_p6(x)
      !! The reference's `p6`: sum over the six simultaneous permutations of
      !! the three (occupied, virtual) index pairs
      real(dp), allocatable, intent(inout) :: x(:, :, :, :, :, :)   !! (i, j, k, a, b, c)

      real(dp), allocatable :: y(:, :, :, :, :, :)
      integer :: no, nv, i, j, k, a, b, c

      no = size(x, 1)
      nv = size(x, 4)
      allocate (y(no, no, no, nv, nv, nv))
      !$omp parallel do collapse(3) default(none) shared(x, y, no, nv) private(a, b, c, i, j, k)
      do c = 1, nv
         do b = 1, nv
            do a = 1, nv
               do k = 1, no
                  do j = 1, no
                     do i = 1, no
                        y(i, j, k, a, b, c) = x(i, j, k, a, b, c) + x(i, k, j, a, c, b) &
                                              + x(j, i, k, b, a, c) + x(k, i, j, c, a, b) &
                                              + x(j, k, i, b, c, a) + x(k, j, i, c, b, a)
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do
      call move_alloc(y, x)
   end subroutine symmetrise_p6

   subroutine spin_sum_r6(x, r)
      !! The reference's `r6`: 4 x plus the two cyclic permutations of the
      !! virtual indices, less twice each of the three transpositions
      real(dp), intent(in) :: x(:, :, :, :, :, :)
      real(dp), allocatable, intent(inout) :: r(:, :, :, :, :, :)

      integer :: no, nv, i, j, k, a, b, c

      no = size(x, 1)
      nv = size(x, 4)
      if (.not. allocated(r)) allocate (r(no, no, no, nv, nv, nv))
      !$omp parallel do collapse(3) default(none) shared(x, r, no, nv) private(a, b, c, i, j, k)
      do c = 1, nv
         do b = 1, nv
            do a = 1, nv
               do k = 1, no
                  do j = 1, no
                     do i = 1, no
                        r(i, j, k, a, b, c) = 4.0_dp*x(i, j, k, a, b, c) &
                                              + x(i, j, k, b, c, a) + x(i, j, k, c, a, b) &
                                              - 2.0_dp*x(i, j, k, c, b, a) &
                                              - 2.0_dp*x(i, j, k, a, c, b) &
                                              - 2.0_dp*x(i, j, k, b, a, c)
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine spin_sum_r6

   subroutine spin_pair_r6(x, r)
      !! The reference's `as_r6`: 2 x less the exchange of the outer virtual
      !! indices and less the exchange of the last two
      real(dp), intent(in) :: x(:, :, :, :, :, :)
      real(dp), allocatable, intent(out) :: r(:, :, :, :, :, :)

      integer :: no, nv, i, j, k, a, b, c

      no = size(x, 1)
      nv = size(x, 4)
      allocate (r(no, no, no, nv, nv, nv))
      !$omp parallel do collapse(3) default(none) shared(x, r, no, nv) private(a, b, c, i, j, k)
      do c = 1, nv
         do b = 1, nv
            do a = 1, nv
               do k = 1, no
                  do j = 1, no
                     do i = 1, no
                        r(i, j, k, a, b, c) = 2.0_dp*x(i, j, k, a, b, c) &
                                              - x(i, j, k, c, b, a) - x(i, j, k, a, c, b)
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine spin_pair_r6

end module mqc_czt_rcc_triples_terms
