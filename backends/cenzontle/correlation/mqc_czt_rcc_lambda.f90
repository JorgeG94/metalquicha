!! The Lambda equations of closed-shell CCSD in spatial orbitals
module mqc_czt_rcc_lambda
   !! The left-hand amplitudes `l1` and `l2` that make the coupled cluster
   !! energy stationary with respect to `t1` and `t2`. Densities, properties
   !! and gradients are built from the pair; the energy needs neither.
   !!
   !! **Equations.** PySCF's `cc/rccsd_lambda.py`, term for term. The amplitude
   !! equations next door in `mqc_czt_rcc` follow the same package's `rccsd.py`,
   !! so the two share one set of conventions and one reference for their
   !! numbers.
   !!
   !! **Notation.** Each contraction is written with `contract` and carries the
   !! reference's index string unchanged, on arrays laid out as the reference
   !! lays them out: `l2(i,j,a,b)`, `ovov(i,a,j,b) = (ia|jb)`,
   !! `ovoo(i,a,j,k) = (ia|jk)`, `ovvo(i,a,b,j) = (ia|bj)`. The integral blocks
   !! this module is handed are permuted into those once.
   !!
   !! **Canonical orbitals are assumed**, as they are next door: the Fock
   !! matrix is diagonal, so every term in `f_ov` is absent and `f_oo` and
   !! `f_vv` are the orbital energies.
   !!
   !! **What is not held.** `(vv|vv)` enters once, contracted with `l2`, through
   !! the same batched ladder the amplitude equations use.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_diis, only: diis_state_t
   use mqc_program_limits, only: MAX_LINE_LENGTH
   use mqc_convergence_report, only: convergence_header, convergence_footer
   use mqc_tensor_contract, only: contract, transposed
   use mqc_czt_rcc, only: rcc_state_t, rcc_ladder
   use mqc_czt_rcc_triples_terms, only: rcc_triples_terms_t
   implicit none
   private

   public :: rcc_lambda_t
   public :: run_czt_rccsd_lambda

   integer, parameter :: DEFAULT_DIIS = 8
      !! Vectors the extrapolation keeps when the caller names no number

   type :: rcc_lambda_t
      !! Converged Lambda amplitudes
      real(dp), allocatable :: l1(:, :)          !! (no, nv)
      real(dp), allocatable :: l2(:, :, :, :)    !! (no, no, nv, nv), as `(i, j, a, b)`
      integer :: iterations = 0
      logical :: converged = .false.
   end type rcc_lambda_t

   type :: lambda_integrals_t
      !! The integral blocks in the reference's layout
      real(dp), allocatable :: oooo(:, :, :, :)   !! (ij|kl)
      real(dp), allocatable :: ovoo(:, :, :, :)   !! (ia|jk)
      real(dp), allocatable :: oovv(:, :, :, :)   !! (ij|ab)
      real(dp), allocatable :: ovov(:, :, :, :)   !! (ia|jb)
      real(dp), allocatable :: ovvo(:, :, :, :)   !! (ia|bj)
      real(dp), allocatable :: ovvv(:, :, :, :)   !! (ia|bc)
   end type lambda_integrals_t

   type :: lambda_intermediates_t
      !! Everything that depends on `t1` and `t2` alone, built once
      real(dp), allocatable :: tau(:, :, :, :)     !! t2 + t1 t1
      real(dp), allocatable :: theta(:, :, :, :)   !! 2 t2 - t2 with a and b exchanged
      real(dp), allocatable :: v1(:, :)            !! (nv, nv)
      real(dp), allocatable :: v2(:, :)            !! (no, no)
      real(dp), allocatable :: v4(:, :)            !! (no, nv)
      real(dp), allocatable :: w3(:, :)            !! (nv, no)
      real(dp), allocatable :: woooo(:, :, :, :)
      real(dp), allocatable :: wovvo(:, :, :, :)   !! (no, nv, nv, no)
      real(dp), allocatable :: wovvo_x(:, :, :, :) !! The reference's `woVVo`
      real(dp), allocatable :: woovo(:, :, :, :)   !! (no, no, nv, no)
      real(dp), allocatable :: wvvvo(:, :, :, :)   !! (nv, nv, nv, no)
   end type lambda_intermediates_t

contains

   subroutine run_czt_rccsd_lambda(state, max_iter, tolerance, verbose, lambda, error, &
                                   diis_vectors, triples)
      !! Solve the Lambda equations for a converged CCSD calculation
      type(rcc_state_t), intent(in) :: state
      integer, intent(in) :: max_iter
      real(dp), intent(in) :: tolerance
         !! On the norm of the change in `l1` and `l2` between iterations
      logical, intent(in) :: verbose
      type(rcc_lambda_t), intent(out) :: lambda
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: diis_vectors
      type(rcc_triples_terms_t), intent(in), optional :: triples
         !! Present, the equations are those of CCSD(T): the perturbative
         !! correction's derivative with respect to the amplitudes is a source
         !! in each residual

      type(lambda_integrals_t) :: ints
      type(lambda_intermediates_t) :: imds
      type(diis_state_t) :: diis
      real(dp), allocatable :: l1(:, :), l2(:, :, :, :), l1n(:, :), l2n(:, :, :, :)
      real(dp), allocatable :: flat(:), step(:)
      character(len=MAX_LINE_LENGTH) :: line
      logical :: extrapolated
      integer :: no, nv, iter, diis_size, n1, n2
      real(dp) :: change

      no = state%no
      nv = state%nv
      if (.not. allocated(state%t2)) then
         call error%set(ERROR_VALIDATION, "RCCSD Lambda: there are no converged "// &
                        "amplitudes to start from")
         return
      end if
      diis_size = DEFAULT_DIIS
      if (present(diis_vectors)) diis_size = diis_vectors

      call reference_layout(state, ints)
      call make_intermediates(state, ints, imds)

      n1 = no*nv
      n2 = no*no*nv*nv
      allocate (l1(no, nv), l2(no, no, nv, nv), l1n(no, nv), l2n(no, no, nv, nv))
      allocate (flat(n1 + n2), step(n1 + n2))
      l1 = state%t1
      l2 = state%t2
      call diis%init(diis_size, size(flat), size(step))

      call convergence_header(verbose, "RCCSD Lambda iterations", &
                              "    iter        |change|   diis", 34)
      do iter = 1, max_iter
         call update_lambda(state, ints, imds, l1, l2, l1n, l2n)
         if (present(triples)) then
            l1n = l1n + triples%l1_source
            l2n = l2n + triples%l2_source
         end if

         flat(1:n1) = reshape(l1n, [n1])
         flat(n1 + 1:) = reshape(l2n, [n2])
         step(1:n1) = reshape(l1n - l1, [n1])
         step(n1 + 1:) = reshape(l2n - l2, [n2])
         change = norm2(step)
         call diis%push(flat, step)
         call diis%extrapolate(flat, extrapolated)
         l1 = reshape(flat(1:n1), [no, nv])
         l2 = reshape(flat(n1 + 1:), [no, no, nv, nv])

         if (verbose) then
            write (line, "(i8,es16.3,i7)") iter, change, diis%count()
            call logger%info(trim(line))
         end if
         lambda%iterations = iter
         if (change < tolerance) then
            lambda%converged = .true.
            exit
         end if
      end do
      call diis%destroy()
      call convergence_footer(verbose, lambda%converged, lambda%iterations, "iterations", 34)

      if (.not. lambda%converged) then
         call error%set(ERROR_VALIDATION, "the RCCSD Lambda equations did not converge")
         return
      end if
      call move_alloc(l1, lambda%l1)
      call move_alloc(l2, lambda%l2)
   end subroutine run_czt_rccsd_lambda

   subroutine reference_layout(state, ints)
      !! The integral blocks permuted into the layout the equations are written in
      type(rcc_state_t), intent(in) :: state
      type(lambda_integrals_t), intent(out) :: ints

      ints%oooo = state%eris%oooo
      ints%oovv = state%eris%oovv
      ints%ovov = state%eris%ovov
      ints%ovvv = state%eris%ovvv
      ! (ia|jk) from (jk|ia), and (ia|bj) from (ia|jb).
      ints%ovoo = transposed(state%eris%ooov, [3, 4, 1, 2])
      ints%ovvo = transposed(state%eris%ovov, [1, 2, 4, 3])
   end subroutine reference_layout

   subroutine make_intermediates(state, ints, imds)
      !! `make_intermediates` of the reference, with `f_ov = 0`
      type(rcc_state_t), intent(in) :: state
      type(lambda_integrals_t), intent(in) :: ints
      type(lambda_intermediates_t), intent(out) :: imds

      real(dp), allocatable :: ovov1(:, :, :, :), ovoo1(:, :, :, :), ovvv1(:, :, :, :)
      real(dp), allocatable :: w_direct(:, :, :, :), w_exchange(:, :, :, :)
      real(dp), allocatable :: w_both(:, :, :, :), w_diff(:, :, :, :)
      real(dp), allocatable :: tmp(:, :, :, :), tmp_vv(:, :), v5(:, :)
      integer :: no, nv, i, j, a, b

      no = state%no
      nv = state%nv
      associate (t1 => state%t1, t2 => state%t2, ovov => ints%ovov, ovoo => ints%ovoo, &
                 ovvv => ints%ovvv)

         allocate (imds%tau(no, no, nv, nv))
         do b = 1, nv
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     imds%tau(i, j, a, b) = t2(i, j, a, b) + t1(i, a)*t1(j, b)
                  end do
               end do
            end do
         end do
         imds%theta = 2.0_dp*t2 - transposed(t2, [1, 2, 4, 3])
         ovov1 = 2.0_dp*ovov - transposed(ovov, [1, 4, 3, 2])
         ovoo1 = 2.0_dp*ovoo - transposed(ovoo, [3, 2, 1, 4])

         allocate (imds%v1(nv, nv), imds%v2(no, no), imds%v4(no, nv), source=0.0_dp)
         do a = 1, nv
            imds%v1(a, a) = state%eps_v(a)
         end do
         do i = 1, no
            imds%v2(i, i) = state%eps_o(i)
         end do
         call contract("jakc,jkbc->ba", ovov1, imds%tau, imds%v1, -1.0_dp)
         call contract("ibkc,jkbc->ij", ovov1, imds%tau, imds%v2, 1.0_dp)
         call contract("kbij,kb->ij", ovoo1, t1, imds%v2, 1.0_dp)
         call contract("jbkc,kc->jb", ovov1, t1, imds%v4, 1.0_dp)

         ! 'kc,kb,jc->bj' in two steps.
         allocate (v5(nv, no), tmp_vv(nv, nv), source=0.0_dp)
         call contract("kc,kb->cb", imds%v4, t1, tmp_vv, 1.0_dp)
         call contract("cb,jc->bj", tmp_vv, t1, v5, 1.0_dp)
         call contract("lckj,klbc->bj", ovoo1, t2, v5, -1.0_dp)

         imds%woooo = ints%oooo
         call contract("icjl,kc->ikjl", ovoo, t1, imds%woooo, 1.0_dp)
         call contract("jcil,kc->iljk", ovoo, t1, imds%woooo, 1.0_dp)
         call contract("icjd,klcd->ikjl", ovov, imds%tau, imds%woooo, 1.0_dp)

         ! The reference's v4OVvo and v4oVVo.
         w_direct = ints%ovvo
         call contract("ldjb,klcd->jbck", ovov1, t2, w_direct, 1.0_dp)
         call contract("ldjb,kldc->jbck", ovov, t2, w_direct, -1.0_dp)
         w_exchange = -transposed(ints%oovv, [1, 4, 3, 2])
         call contract("jdlb,kldc->jbck", ovov, t2, w_exchange, 1.0_dp)
         w_both = 2.0_dp*w_direct + w_exchange
         w_diff = w_direct - w_exchange

         allocate (imds%w3(nv, no), source=0.0_dp)
         call contract("jbck,jb->ck", w_both, t1, imds%w3, 1.0_dp)

         allocate (tmp(no, no, nv, no), source=0.0_dp)
         call contract("ibck,jb->ijck", w_both, t1, tmp, 1.0_dp)
         imds%woovo = tmp - transposed(tmp, [1, 4, 3, 2])
         deallocate (tmp)
         call contract("ibck,jb->ikcj", w_diff, t1, imds%woovo, 1.0_dp)
         imds%woovo = imds%woovo + transposed(ovoo1, [4, 3, 2, 1])
         call contract("lcik,jlbc->ikbj", ovoo1, imds%theta, imds%woovo, 1.0_dp)
         call contract("lcik,jlbc->ijbk", ovoo1, t2, imds%woovo, -1.0_dp)
         call contract("iclk,ljbc->ijbk", ovoo1, t2, imds%woovo, -1.0_dp)

         allocate (tmp(nv, nv, nv, no), source=0.0_dp)
         call contract("jack,jb->back", w_both, t1, tmp, 1.0_dp)
         imds%wvvvo = tmp - transposed(tmp, [3, 2, 1, 4])
         deallocate (tmp)
         call contract("jack,jb->cabk", w_diff, t1, imds%wvvvo, 1.0_dp)
         call contract("lajk,jlbc->cabk", ovoo1, imds%tau, imds%wvvvo, -1.0_dp)
         deallocate (w_both, w_diff)

         ! From here the two are the reference's wOVvo and woVVo.
         ! 'jbld,kd,lc->jbck' in two steps, once for each.
         allocate (tmp(no, nv, no, no), source=0.0_dp)
         call contract("jbld,kd->jblk", ovov, t1, tmp, 1.0_dp)
         call contract("jblk,lc->jbck", tmp, t1, w_direct, -1.0_dp)
         tmp = 0.0_dp
         call contract("jdlb,kd->jblk", ovov, t1, tmp, 1.0_dp)
         call contract("jblk,lc->jbck", tmp, t1, w_exchange, 1.0_dp)
         deallocate (tmp)
         call contract("jblk,lc->jbck", ovoo, t1, w_direct, -1.0_dp)
         call contract("lbjk,lc->jbck", ovoo, t1, w_exchange, 1.0_dp)

         call contract("kacd,kjbd->bacj", ovvv, t2, imds%wvvvo, 1.5_dp)
         call contract("jbcd,kd->jbck", ovvv, t1, w_direct, 1.0_dp)
         call contract("jdcb,kd->jbck", ovvv, t1, w_exchange, -1.0_dp)

         ovvv1 = 2.0_dp*ovvv - transposed(ovvv, [1, 4, 3, 2])
         call contract("jcba,jc->ba", ovvv1, t1, imds%v1, 1.0_dp)
         call contract("kdbc,jkcd->bj", ovvv1, t2, v5, 1.0_dp)
         call contract("idcb,jkdb->ijck", ovvv1, imds%tau, imds%woovo, 1.0_dp)

         allocate (tmp(nv, nv, nv, no), source=0.0_dp)
         call contract("kdca,jkbd->cabj", ovvv1, imds%theta, tmp, 1.0_dp)
         imds%wvvvo = imds%wvvvo - tmp + 0.5_dp*transposed(tmp, [3, 2, 1, 4]) &
                      - transposed(ovvv1, [4, 3, 2, 1])
         deallocate (tmp)

         imds%w3 = imds%w3 + v5
         call contract("cb,jb->cj", imds%v1, t1, imds%w3, 1.0_dp)
         call contract("jk,jb->bk", imds%v2, t1, imds%w3, -1.0_dp)

         imds%wovvo = 2.0_dp*w_direct + w_exchange
         call move_alloc(w_exchange, imds%wovvo_x)
      end associate
   end subroutine make_intermediates

   subroutine update_lambda(state, ints, imds, l1, l2, l1n, l2n)
      !! One pass of `update_lambda` of the reference, with `f_ov = 0`
      type(rcc_state_t), intent(in) :: state
      type(lambda_integrals_t), intent(in) :: ints
      type(lambda_intermediates_t), intent(in) :: imds
      real(dp), intent(in) :: l1(:, :), l2(:, :, :, :)
      real(dp), intent(out) :: l1n(:, :), l2n(:, :, :, :)

      real(dp), allocatable :: mvv(:, :), moo(:, :), mvv1(:, :), moo1(:, :), tmp1(:, :)
      real(dp), allocatable :: m3(:, :, :, :), work(:, :, :, :), ovov1(:, :, :, :)
      real(dp), allocatable :: l2theta(:, :, :, :)
      integer :: no, nv, i, j, a, b
      real(dp) :: gap

      no = state%no
      nv = state%nv
      associate (t1 => state%t1, ovov => ints%ovov, ovoo => ints%ovoo, ovvv => ints%ovvv)

         l1n = 0.0_dp
         allocate (mvv(nv, nv), moo(no, no), source=0.0_dp)
         call contract("klca,klcb->ba", l2, imds%theta, mvv, 1.0_dp)
         call contract("kicd,kjcd->ij", l2, imds%theta, moo, 1.0_dp)
         mvv1 = mvv
         moo1 = moo
         call contract("jc,jb->bc", l1, t1, mvv1, 1.0_dp)
         call contract("ic,kc->ik", l1, t1, moo1, 1.0_dp)

         ! m3(i,j,c,d) = sum_ab l2(i,j,a,b) (ac|bd), then the hole ladder.
         allocate (m3(no, no, nv, nv), source=0.0_dp)
         call rcc_ladder(state%eris, l2, no, nv, m3)
         call contract("klab,ikjl->ijab", l2, imds%woooo, m3, 1.0_dp)
         m3 = 0.5_dp*m3
         allocate (work(no, no, no, no), source=0.0_dp)
         call contract("ijcd,klcd->ijkl", l2, imds%tau, work, 1.0_dp)
         call contract("kalb,ijkl->ijab", ovov, work, m3, 0.5_dp)
         deallocate (work)

         l2n = 0.5_dp*transposed(ovov, [1, 3, 2, 4])
         call contract("ijac,cb->ijab", l2, imds%v1, l2n, 1.0_dp)
         call contract("ikab,jk->ijab", l2, imds%v2, l2n, -1.0_dp)
         call contract("ca,icjb->ijab", mvv1, ovov, l2n, -1.0_dp)
         call contract("ik,kajb->ijab", moo1, ovov, l2n, -1.0_dp)

         call contract("ik,ka->ia", moo, imds%v4, l1n, -1.0_dp)
         call contract("ca,ic->ia", mvv, imds%v4, l1n, -1.0_dp)
         call contract("ia,jb->ijab", l1, imds%v4, l2n, 1.0_dp)

         tmp1 = t1
         call contract("kc,kjcb->jb", l1, imds%theta, tmp1, 1.0_dp)
         call contract("bd,jd->jb", mvv1, t1, tmp1, -1.0_dp)
         call contract("lj,lb->jb", moo, t1, tmp1, -1.0_dp)
         ovov1 = 2.0_dp*ovov - transposed(ovov, [1, 4, 3, 2])
         call contract("jbia,jb->ia", ovov1, tmp1, l1n, 1.0_dp)
         deallocate (ovov1, tmp1)

         call contract("iacb,bc->ia", ovvv, mvv1, l1n, 2.0_dp)
         call contract("ibca,bc->ia", ovvv, mvv1, l1n, -1.0_dp)
         call contract("ic,jbca->jiba", l1, ovvv, l2n, 1.0_dp)
         allocate (work(no, no, nv, no), source=0.0_dp)
         call contract("ijcd,kd->ijck", l2, t1, work, 1.0_dp)
         call contract("kbca,ijck->ijab", ovvv, work, m3, -1.0_dp)
         deallocate (work)

         l2n = l2n + m3
         call contract("ijab,jb->ia", m3, t1, l1n, 2.0_dp)
         call contract("jiba,jb->ia", m3, t1, l1n, 2.0_dp)
         call contract("ijba,jb->ia", m3, t1, l1n, -1.0_dp)
         call contract("jiab,jb->ia", m3, t1, l1n, -1.0_dp)
         deallocate (m3)

         call contract("iajk,kj->ia", ovoo, moo1, l1n, -2.0_dp)
         call contract("jaik,kj->ia", ovoo, moo1, l1n, 1.0_dp)
         call contract("ka,jbik->ijab", l1, ovoo, l2n, -1.0_dp)

         l2theta = 2.0_dp*l2 - transposed(l2, [1, 2, 4, 3])
         call contract("ikac,jbck->ijab", l2theta, imds%wovvo, l2n, 0.5_dp)
         deallocate (l2theta)
         allocate (work(no, no, nv, nv), source=0.0_dp)
         call contract("ikca,jbck->ijab", l2, imds%wovvo_x, work, 1.0_dp)
         l2n = l2n + 0.5_dp*work + transposed(work, [2, 1, 3, 4])
         deallocate (work)

         call contract("ib,ba->ia", l1, imds%v1, l1n, 1.0_dp)
         call contract("ja,ij->ia", l1, imds%v2, l1n, -1.0_dp)
         call contract("jb,iabj->ia", l1, ints%ovvo, l1n, 2.0_dp)
         call contract("jb,ijba->ia", l1, ints%oovv, l1n, -1.0_dp)
         call contract("ijbc,bacj->ia", l2, imds%wvvvo, l1n, -1.0_dp)
         call contract("kjca,ijck->ia", l2, imds%woovo, l1n, -1.0_dp)
         call contract("ijab,bj->ia", l2, imds%w3, l1n, 2.0_dp)
         call contract("ijba,bj->ia", l2, imds%w3, l1n, -1.0_dp)
      end associate

      ! The residual becomes an update: divide by the orbital energy
      ! denominator and add the amplitude it was a residual for, which takes
      ! back the diagonal of `v1` and `v2`.
      do a = 1, nv
         do i = 1, no
            l1n(i, a) = l1n(i, a)/(state%eps_o(i) - state%eps_v(a)) + l1(i, a)
         end do
      end do
      l2n = l2n + transposed(l2n, [2, 1, 4, 3])
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  gap = state%eps_o(i) - state%eps_v(a) + state%eps_o(j) - state%eps_v(b)
                  l2n(i, j, a, b) = l2n(i, j, a, b)/gap + l2(i, j, a, b)
               end do
            end do
         end do
      end do
   end subroutine update_lambda

end module mqc_czt_rcc_lambda
