!! Closed-shell coupled cluster in spatial orbitals
module mqc_libcint_rcc
   !! Spin-adapted (restricted) CCSD over an RHF reference.
   !!
   !! The companion to `mqc_libcint_cc`, which does the same physics in spin
   !! orbitals. That module says of itself that spin orbitals are the reference
   !! formulation and to "revisit if CCSD becomes something people run rather
   !! than something people check against". This is that revisit. Neither
   !! replaces the other: the spin-orbital path stays as the oracle, because for
   !! a closed shell the two must agree to machine precision and that identity
   !! is a far sharper test than any reference energy.
   !!
   !! **Why.** Everything the spin-orbital path indexes is twice as long in each
   !! of four dimensions. `t2` is `nv^2 no^2` where this is `n_v^2 n_o^2` --
   !! sixteen times -- and DIIS keeps sixteen copies of it, which is the single
   !! largest allocation a CCSD run makes and is invisible until measured. The
   !! `n_occ n_vir^3` integral blocks carry the same factor. Benzene/cc-pVDZ
   !! goes from about 17 GB to 2.2 GB, and on the fitted path, where the
   !! `n_act^4` tensor never exists, from 15.9 GB to 0.9 GB.
   !!
   !! **Equations.** Hirata, Podeszwa, Tobita and Bartlett, JCP 120, 2581 (2004),
   !! Eqs. (35)-(45), followed term for term as PySCF's `cc/rccsd.py` and
   !! `cc/rintermediates.py` write them. That source rather than the paper on
   !! purpose: every reference energy in `validation/check_cc.f90` was generated
   !! by PySCF, so matching its equations means matching its answer by
   !! construction, and a disagreement is debuggable against a readable
   !! implementation instead of against a table.
   !!
   !! **Notation.** Chemists', spatial, and each array is indexed in the order
   !! its name reads: `ovov(i,a,j,b)` is `(ia|jb)`. The einsum string from the
   !! reference is quoted above every contraction, because the translation from
   !! row-major numpy to column-major Fortran is the one place this can go
   !! quietly wrong, and a transposed index gives a converged wrong number.
   !!
   !! Real orbitals throughout, which is what lets `(ia|bj)` be read straight out
   !! of `ovov` as `(ia|jb)` and `(ia|jk)` out of `ooov` as `(jk|ia)`. Both
   !! identities need `(pq|rs) = (pq|sr) = (rs|pq)` and would not hold for
   !! complex orbitals; PySCF carries separate `ovvo` and `ovoo` blocks because
   !! it supports those.
   !!
   !! **Canonical orbitals are assumed**, as they are next door. The Fock matrix
   !! is then diagonal, so `f_ov` vanishes and the occupied-occupied and
   !! virtual-virtual blocks are exactly the orbital energies that go into the
   !! denominators. The reference adds `foo`/`fvv` into the intermediates and
   !! subtracts the same diagonal back out a few lines later; here neither
   !! happens, which is the same arithmetic with two fewer chances to get a sign
   !! wrong.
   !!
   !! **What is not held.** `(vv|vv)` is never materialised in the fitted case
   !! and is read once in the conventional one. The particle-particle ladder is
   !! the only consumer, and it works a batch of `(cd)` columns at a time -- the
   !! same argument, and the same batching, as `particle_ladder` next door.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use omp_lib, only: omp_get_max_threads
   use mqc_timing, only: timing_report_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_diis, only: diis_state_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_df_mo_block
   use mqc_libcint_mp2, only: transform_block
   use pic_logger, only: logger => global_logger
   use mqc_convergence_report, only: convergence_header, convergence_footer
   use mqc_program_limits, only: MAX_LINE_LENGTH
   implicit none
   private

   public :: rcc_result_t
   public :: run_libcint_rccsd
   public :: rcc_eris_t              !! Exported so a test can build them directly
   public :: build_rcc_eris_conventional
   public :: rccsd_correlation_energy  !! Exported for the MP2 identity test
   public :: ri_ladder_prefers_direct  !! Exported so a test can assert which ladder ran

   !> Columns of the compound (cd) index assembled per pass of the ladder.
   !> Same role and same reasoning as LADDER_BATCH in the spin-orbital module:
   !> the array held is n_vir^2 by this, so the memory stays O(n_vir^2).
   integer, parameter :: LADDER_BATCH = 256

   character(len=*), parameter :: STAGE_ITER = "RCCSD iterations"

   !> Orderings of a triple of labels. Three distinct virtual orbitals name the
   !> same physical triple six ways, and the (T) energy expression sums over
   !> every pairing of one such ordering with one ordering of the occupied
   !> labels -- so this six is 3! and nothing else, and appears wherever either
   !> set is enumerated.
   integer, parameter :: N_TRIPLE_PERMS = 6

   !> How many of a triple's innermost virtuals one pass of the triples holds.
   !>
   !> The batching that made the particle half of W into a decent matrix product
   !> wants this large; the fact that each thread keeps its own copy of the
   !> result -- n_occ^3 by six by this -- wants it small. At 32 the products are
   !> still several hundred thousand flops apiece, and a forty-electron case
   !> costs 98 MB a thread instead of the 614 MB an unbounded batch would.
   integer, parameter :: TRIPLES_C_BATCH = 32

   type :: rcc_eris_t
      !! Spatial MO integrals in chemists' notation, by block
      !!
      !! Exactly the five blocks with at least one occupied index, plus whichever
      !! of `vvvv` and `b_vv` the ladder is to read. Those two are alternatives
      !! and never both: `vvvv` is the conventional path's, `b_vv` the fitted
      !! one's three-index tensor, and every routine below decides which it has
      !! by which one is allocated.
      real(dp), allocatable :: oooo(:, :, :, :)   !! (ij|kl)
      real(dp), allocatable :: ooov(:, :, :, :)   !! (ij|ka)
      real(dp), allocatable :: oovv(:, :, :, :)   !! (ij|ab)
      real(dp), allocatable :: ovov(:, :, :, :)   !! (ia|jb)
      real(dp), allocatable :: ovvv(:, :, :, :)   !! (ia|bc)
      real(dp), allocatable :: vvvv(:, :, :, :)   !! (ab|cd), conventional only
      real(dp), allocatable :: b_vv(:, :)         !! (ab, P), fitted only
   end type rcc_eris_t

   type :: rcc_result_t
      !! What a converged spin-adapted coupled cluster calculation leaves behind
      real(dp) :: e_singles = 0.0_dp
      real(dp) :: e_doubles = 0.0_dp
      real(dp) :: e_triples = 0.0_dp   !! (T); not implemented on this path yet
      real(dp) :: e_mp2 = 0.0_dp
         !! MP2 from these same spatial integrals -- the first iteration's
         !! energy, so free, and the sharpest check available on the block
         !! layout before any amplitude equation runs.
      real(dp) :: e_correlation = 0.0_dp
      integer :: iterations = 0
      logical :: converged = .false.
   end type rcc_result_t

contains

   pure function rcc_megabytes(no, nv, conventional) result(mb)
      !! What the integral blocks cost, in MB
      integer, intent(in) :: no, nv
      logical, intent(in) :: conventional
      real(dp) :: mb

      real(dp) :: elements

      elements = real(no, dp)**4 &                          ! oooo
                 + real(no, dp)**3*real(nv, dp) &           ! ooov
                 + 2.0_dp*real(no, dp)**2*real(nv, dp)**2 &  ! oovv, ovov
                 + real(no, dp)*real(nv, dp)**3 &           ! ovvv
                 + real(nv, dp)**2*real(LADDER_BATCH, dp)   ! the ladder batch
      if (conventional) elements = elements + real(nv, dp)**4
      mb = elements*8.0_dp/1.0e6_dp
   end function rcc_megabytes

   subroutine build_rcc_eris_conventional(mol, c_occ, c_vir, eris, clk)
      !! The five occupied-bearing blocks and (vv|vv), by direct transform
      !!
      !! Six calls to `transform_block` rather than one whole-active-space
      !! transform and six slices of it. The total arithmetic is the same --
      !! these blocks partition that tensor -- but the `n_act^4` array is never
      !! allocated, and it is the term the conventional path's memory used to be
      !! set by. `(vv|vv)` at `n_vir^4` is what remains, and it is smaller by
      !! roughly `(n_act/n_vir)^4`.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :)
      type(rcc_eris_t), intent(out) :: eris
      type(timing_report_t), intent(inout), optional :: clk
         !! Lapped here rather than around the call, because the AO build and
         !! the transform are separate costs with separate fixes and a single
         !! lap around both cannot tell you which one you are paying for.

      real(dp), allocatable :: eri(:, :)

      if (present(clk)) call clk%begin("AO integrals")
      call mol%eris_packed(eri)
      if (present(clk)) call clk%lap()
      if (present(clk)) call clk%begin("AO->MO transform")
      call transform_block(eri, c_occ, c_occ, c_occ, c_occ, eris%oooo)
      call transform_block(eri, c_occ, c_occ, c_occ, c_vir, eris%ooov)
      call transform_block(eri, c_occ, c_occ, c_vir, c_vir, eris%oovv)
      call transform_block(eri, c_occ, c_vir, c_occ, c_vir, eris%ovov)
      call transform_block(eri, c_occ, c_vir, c_vir, c_vir, eris%ovvv)
      call transform_block(eri, c_vir, c_vir, c_vir, c_vir, eris%vvvv)
      if (present(clk)) call clk%lap()
      deallocate (eri)
   end subroutine build_rcc_eris_conventional

   subroutine build_rcc_eris_fitted(mol, aux, c_occ, c_vir, eris, error)
      !! The same blocks from the three-index fitted tensor
      !!
      !! `b_vv` outlives this routine because the ladder builds `(ac|bd)` from it
      !! a batch at a time; nothing of fourth order in the virtuals is ever
      !! allocated on this path, which is the whole of what RI buys here.
      type(libcint_molecule_t), intent(in) :: mol, aux
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :)
      type(rcc_eris_t), intent(out) :: eris
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: b_oo(:, :), b_ov(:, :), b_vv_pq(:, :)
      integer :: no, nv

      no = size(c_occ, 2)
      nv = size(c_vir, 2)

      ! Every block is consumed as `B B^T`, so the cheap factor is safe.
      call build_df_mo_block(mol, aux, c_occ, c_occ, b_oo, error, fast_factor=.true.)
      if (error%has_error()) return
      call build_df_mo_block(mol, aux, c_occ, c_vir, b_ov, error, fast_factor=.true.)
      if (error%has_error()) return
      call build_df_mo_block(mol, aux, c_vir, c_vir, b_vv_pq, error, fast_factor=.true.)
      if (error%has_error()) return

      allocate (eris%oooo(no, no, no, no), eris%ooov(no, no, no, nv))
      allocate (eris%oovv(no, no, nv, nv), eris%ovov(no, nv, no, nv))
      allocate (eris%ovvv(no, nv, nv, nv))

      call df_block_gemm(b_oo, b_oo, no*no, no*no, eris%oooo)
      call df_block_gemm(b_oo, b_ov, no*no, no*nv, eris%ooov)
      call df_block_gemm(b_oo, b_vv_pq, no*no, nv*nv, eris%oovv)
      call df_block_gemm(b_ov, b_ov, no*nv, no*nv, eris%ovov)
      call df_block_gemm(b_ov, b_vv_pq, no*nv, nv*nv, eris%ovvv)

      deallocate (b_oo, b_ov)

      ! Kept as (ab, P): the ladder slices it by compound virtual pair, which is
      ! the leading dimension here and so contiguous.
      call move_alloc(b_vv_pq, eris%b_vv)
   end subroutine build_rcc_eris_fitted

   subroutine df_block_gemm(bl, br, nl, nr, block)
      !! (left|right) = sum_P B_left(left,P) B_right(right,P)
      !!
      !! Explicit-shape destination so sequence association lets the gemm write
      !! straight into a rank-four array, exactly as the spin-orbital module's
      !! namesake does and for the same reason: at `(ov|vv)` the reshaping
      !! temporary would be the peak allocation.
      real(dp), intent(in) :: bl(:, :), br(:, :)
      integer, intent(in) :: nl, nr
      real(dp), intent(out) :: block(nl, nr)

      call pic_gemm(bl, br, block, transb="T")
   end subroutine df_block_gemm

   !===========================================================================
   ! Intermediates -- Hirata Eqs. (37)-(45), as rintermediates.py writes them
   !
   ! Every routine quotes the einsum it implements. Two index identities are
   ! used throughout to avoid storing blocks that are transposes of ones we
   ! already have, and both need real orbitals:
   !
   !     eris.ovvo(k,c,a,i) = (kc|ai) = (kc|ia) = ovov(k,c,i,a)
   !     eris.ovoo(i,a,j,k) = (ia|jk) = (jk|ia) = ooov(j,k,i,a)
   !===========================================================================

   subroutine build_tau(t1, t2, no, nv, tau)
      !! tau(i,j,a,b) = t2(i,j,a,b) + t1(i,a) t1(j,b)
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: tau(:, :, :, :)

      integer :: i, j, a, b

      !$omp parallel do default(none) shared(tau, t1, t2, no, nv) &
      !$omp    private(i, j, a, b) schedule(static)
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  tau(i, j, a, b) = t2(i, j, a, b) + t1(i, a)*t1(j, b)
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine build_tau

   pure function lvec(ovov, k, c, l, d) result(v)
      !! The recurring 2(kc|ld) - (kd|lc) combination
      !!
      !! It appears in Foo, Fvv, Fov and the energy, and writing it once is what
      !! keeps the "2 minus exchange" spin-adaptation factor from being typed
      !! four times with four chances to swap c and d.
      real(dp), intent(in) :: ovov(:, :, :, :)
      integer, intent(in) :: k, c, l, d
      real(dp) :: v

      v = 2.0_dp*ovov(k, c, l, d) - ovov(k, d, l, c)
   end function lvec

   subroutine cc_foo(eris, t1, tau, no, nv, fki)
      !! Fki = sum_lcd [2(kc|ld) - (kd|lc)] tau(i,l,c,d)
      !!
      !! rintermediates.cc_Foo, with its four terms folded into tau: the two t2
      !! ones and the two t1 t1 ones differ only by which of t2 and t1 t1 they
      !! carry, which is the definition of tau.
      !!
      !! The bare `foo` the reference adds here is absent: for canonical
      !! orbitals it is diagonal and `update_amps` subtracts that same diagonal
      !! back out immediately afterwards.
      type(rcc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: t1(:, :), tau(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: fki(:, :)

      integer :: k, i, l, c, d

      fki = 0.0_dp
      do i = 1, no
         do d = 1, nv
            do c = 1, nv
               do l = 1, no
                  do k = 1, no
                     fki(k, i) = fki(k, i) + lvec(eris%ovov, k, c, l, d)*tau(i, l, c, d)
                  end do
               end do
            end do
         end do
      end do
   end subroutine cc_foo

   subroutine cc_fvv(eris, t1, tau, no, nv, fac)
      !! Fac = -sum_kld [2(kc|ld) - (kd|lc)] tau(k,l,a,d)
      type(rcc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: t1(:, :), tau(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: fac(:, :)

      integer :: a, c, k, l, d

      fac = 0.0_dp
      !$omp parallel do default(none) shared(fac, eris, tau, no, nv) &
      !$omp    private(c, d, l, k, a) schedule(static)
      do c = 1, nv
         do d = 1, nv
            do l = 1, no
               do k = 1, no
                  do a = 1, nv
                     fac(a, c) = fac(a, c) - lvec(eris%ovov, k, c, l, d)*tau(k, l, a, d)
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine cc_fvv

   subroutine cc_fov(eris, t1, no, nv, fkc)
      !! Fkc = sum_ld [2(kc|ld) - (kd|lc)] t1(l,d)
      type(rcc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: t1(:, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: fkc(:, :)

      integer :: k, c, l, d

      fkc = 0.0_dp
      do d = 1, nv
         do l = 1, no
            do c = 1, nv
               do k = 1, no
                  fkc(k, c) = fkc(k, c) + lvec(eris%ovov, k, c, l, d)*t1(l, d)
               end do
            end do
         end do
      end do
   end subroutine cc_fov

   subroutine cc_loo(eris, fki, t1, no, nv, lki)
      !! Lki = Fki + sum_lc [2 ooov(k,i,l,c) - ooov(l,i,k,c)] t1(l,c)
      !!
      !! rintermediates.Loo. Its `2 einsum('lcki,lc->ki', ovoo, t1)` reads
      !! `(lc|ki) = ooov(k,i,l,c)`, and its exchange partner `(kc|li) =
      !! ooov(l,i,k,c)`.
      type(rcc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: fki(:, :), t1(:, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: lki(:, :)

      integer :: k, i, l, c

      lki = fki
      do c = 1, nv
         do l = 1, no
            do i = 1, no
               do k = 1, no
                  lki(k, i) = lki(k, i) &
                              + (2.0_dp*eris%ooov(k, i, l, c) - eris%ooov(l, i, k, c))*t1(l, c)
               end do
            end do
         end do
      end do
   end subroutine cc_loo

   subroutine cc_lvv(eris, fac, t1, no, nv, lac)
      !! Lac = Fac + sum_kd [2 ovvv(k,d,a,c) - ovvv(k,c,a,d)] t1(k,d)
      type(rcc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: fac(:, :), t1(:, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: lac(:, :)

      integer :: a, c, k, d

      lac = fac
      do d = 1, nv
         do k = 1, no
            do c = 1, nv
               do a = 1, nv
                  lac(a, c) = lac(a, c) &
                              + (2.0_dp*eris%ovvv(k, d, a, c) - eris%ovvv(k, c, a, d))*t1(k, d)
               end do
            end do
         end do
      end do
   end subroutine cc_lvv

   subroutine cc_woooo(eris, t1, tau, no, nv, w)
      !! W(k,l,i,j) = (ki|lj) + sum_c ooov(k,i,l,c) t1(j,c) + sum_c ooov(l,j,k,c) t1(i,c)
      !!                      + sum_cd (kc|ld) tau(i,j,c,d)
      type(rcc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: t1(:, :), tau(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: w(:, :, :, :)

      integer :: k, l, i, j, c, d

      do j = 1, no
         do i = 1, no
            do l = 1, no
               do k = 1, no
                  w(k, l, i, j) = eris%oooo(k, i, l, j)
               end do
            end do
         end do
      end do

      !$omp parallel do default(none) shared(w, eris, t1, no, nv) &
      !$omp    private(j, c, i, l, k) schedule(static)
      do j = 1, no
         do c = 1, nv
            do i = 1, no
               do l = 1, no
                  do k = 1, no
                     w(k, l, i, j) = w(k, l, i, j) &
                                     + eris%ooov(k, i, l, c)*t1(j, c) &
                                     + eris%ooov(l, j, k, c)*t1(i, c)
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(w, eris, tau, no, nv) &
      !$omp    private(j, d, c, i, l, k) schedule(static)
      do j = 1, no
         do d = 1, nv
            do c = 1, nv
               do i = 1, no
                  do l = 1, no
                     do k = 1, no
                        w(k, l, i, j) = w(k, l, i, j) + eris%ovov(k, c, l, d)*tau(i, j, c, d)
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine cc_woooo

   subroutine cc_wvoov(eris, t1, t2, no, nv, w)
      !! W(a,k,i,c), rintermediates.cc_Wvoov
      !!
      !!   + sum_d ovvv(k,c,a,d) t1(i,d)
      !!   - sum_l ooov(l,i,k,c) t1(l,a)
      !!   + ovov(k,c,i,a)
      !!   - 1/2 sum_ld ovov(l,d,k,c) t2(i,l,d,a)
      !!   - 1/2 sum_ld ovov(l,c,k,d) t2(i,l,a,d)
      !!   -     sum_ld ovov(l,d,k,c) t1(i,d) t1(l,a)
      !!   +     sum_ld ovov(l,d,k,c) t2(i,l,a,d)
      type(rcc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: w(:, :, :, :)

      integer :: a, k, i, c, d, l

      do c = 1, nv
         do i = 1, no
            do k = 1, no
               do a = 1, nv
                  w(a, k, i, c) = eris%ovov(k, c, i, a)
               end do
            end do
         end do
      end do

      ! `c` is hoisted outermost in both nests so that each thread owns a
      ! disjoint slice of W. The contracted index moves inside, which costs a
      ! little locality and buys the loop being safe to split at all.
      !$omp parallel do default(none) shared(w, eris, t1, no, nv) &
      !$omp    private(c, d, i, k, a) schedule(static)
      do c = 1, nv
         do d = 1, nv
            do i = 1, no
               do k = 1, no
                  do a = 1, nv
                     w(a, k, i, c) = w(a, k, i, c) + eris%ovvv(k, c, a, d)*t1(i, d)
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(w, eris, t1, no, nv) &
      !$omp    private(c, l, i, k, a) schedule(static)
      do c = 1, nv
         do l = 1, no
            do i = 1, no
               do k = 1, no
                  do a = 1, nv
                     w(a, k, i, c) = w(a, k, i, c) - eris%ooov(l, i, k, c)*t1(l, a)
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! The four amplitude terms, as two gemms over the contracted (l,d).
      !
      ! They group by which of the two integral orderings they carry, and the
      ! grouping is the point: written as the reference writes them this is
      ! four passes of O(n_occ^3 n_vir^3) with a strided read of `ovov` in the
      ! innermost position. Grouped, it is two matrix products over compound
      ! indices, and the packing that makes them possible is
      ! O(n_occ^2 n_vir^2).
      block
         real(dp), allocatable :: a1(:, :), a2(:, :), b1(:, :), b2(:, :), r(:, :)
         integer :: nov, kc, ld, ia

         nov = no*nv
         allocate (a1(nov, nov), a2(nov, nov), b1(nov, nov), b2(nov, nov), r(nov, nov))

         do d = 1, nv
            do l = 1, no
               ld = (d - 1)*no + l
               do c = 1, nv
                  do k = 1, no
                     kc = (c - 1)*no + k
                     a1(kc, ld) = eris%ovov(l, d, k, c)
                     a2(kc, ld) = eris%ovov(l, c, k, d)
                  end do
               end do
            end do
         end do

         do a = 1, nv
            do i = 1, no
               ia = (a - 1)*no + i
               do d = 1, nv
                  do l = 1, no
                     ld = (d - 1)*no + l
                     b1(ld, ia) = -0.5_dp*t2(i, l, d, a) - t1(i, d)*t1(l, a) &
                                  + t2(i, l, a, d)
                     b2(ld, ia) = -0.5_dp*t2(i, l, a, d)
                  end do
               end do
            end do
         end do

         call gemm_over_columns(a1, b1, r)
         call accumulate_wvoov(r, no, nv, w)
         call gemm_over_columns(a2, b2, r)
         call accumulate_wvoov(r, no, nv, w)

         deallocate (a1, a2, b1, b2, r)
      end block
   end subroutine cc_wvoov

   subroutine cc_wvovo(eris, t1, t2, no, nv, w)
      !! W(a,k,c,i), rintermediates.cc_Wvovo
      !!
      !!   + sum_d ovvv(k,d,a,c) t1(i,d)
      !!   - sum_l ooov(k,i,l,c) t1(l,a)
      !!   + oovv(k,i,a,c)
      !!   - 1/2 sum_ld ovov(l,c,k,d) t2(i,l,d,a)
      !!   -     sum_ld ovov(l,c,k,d) t1(i,d) t1(l,a)
      type(rcc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: w(:, :, :, :)

      integer :: a, k, c, i, d, l

      do i = 1, no
         do c = 1, nv
            do k = 1, no
               do a = 1, nv
                  w(a, k, c, i) = eris%oovv(k, i, a, c)
               end do
            end do
         end do
      end do

      !$omp parallel do default(none) shared(w, eris, t1, no, nv) &
      !$omp    private(c, d, i, k, a) schedule(static)
      do c = 1, nv
         do d = 1, nv
            do i = 1, no
               do k = 1, no
                  do a = 1, nv
                     w(a, k, c, i) = w(a, k, c, i) + eris%ovvv(k, d, a, c)*t1(i, d)
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(w, eris, t1, no, nv) &
      !$omp    private(c, l, i, k, a) schedule(static)
      do c = 1, nv
         do l = 1, no
            do i = 1, no
               do k = 1, no
                  do a = 1, nv
                     w(a, k, c, i) = w(a, k, c, i) - eris%ooov(k, i, l, c)*t1(l, a)
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! Both amplitude terms carry the same integral ordering, so they are one
      ! gemm over (l,d) rather than two passes.
      block
         real(dp), allocatable :: a2(:, :), b3(:, :), r(:, :)
         integer :: nov, kc, ld, ia

         nov = no*nv
         allocate (a2(nov, nov), b3(nov, nov), r(nov, nov))

         do d = 1, nv
            do l = 1, no
               ld = (d - 1)*no + l
               do c = 1, nv
                  do k = 1, no
                     kc = (c - 1)*no + k
                     a2(kc, ld) = eris%ovov(l, c, k, d)
                  end do
               end do
            end do
         end do

         do a = 1, nv
            do i = 1, no
               ia = (a - 1)*no + i
               do d = 1, nv
                  do l = 1, no
                     ld = (d - 1)*no + l
                     b3(ld, ia) = -0.5_dp*t2(i, l, d, a) - t1(i, d)*t1(l, a)
                  end do
               end do
            end do
         end do

         call gemm_over_columns(a2, b3, r)
         call accumulate_wvovo(r, no, nv, w)

         deallocate (a2, b3, r)
      end block
   end subroutine cc_wvovo

   !===========================================================================
   ! Energy and amplitudes -- Hirata Eqs. (35)-(36), as rccsd.py writes them
   !===========================================================================

   subroutine rccsd_correlation_energy(eris, t1, t2, no, nv, e_singles, e_doubles)
      !! E = sum_ijab [2(ia|jb) - (ib|ja)] (t2(i,j,a,b) + t1(i,a) t1(j,b))
      !!
      !! Split at the tau that produced it rather than by diagram: the doubles
      !! part is what t2 alone gives and the singles part what the t1 t1 product
      !! adds, which is the decomposition the spin-orbital module reports and so
      !! the one the two can be compared term by term with. The bare `2 f_ov t1`
      !! term of the reference is zero for canonical orbitals.
      type(rcc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: e_singles, e_doubles

      integer :: i, j, a, b
      real(dp) :: l

      e_singles = 0.0_dp
      e_doubles = 0.0_dp
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  l = 2.0_dp*eris%ovov(i, a, j, b) - eris%ovov(i, b, j, a)
                  e_doubles = e_doubles + l*t2(i, j, a, b)
                  e_singles = e_singles + l*t1(i, a)*t1(j, b)
               end do
            end do
         end do
      end do
   end subroutine rccsd_correlation_energy

   subroutine mp2_guess(eris, eps_o, eps_v, no, nv, t2, e_mp2)
      !! First-order doubles and the MP2 energy they give
      !!
      !!     t2(i,j,a,b) = (ia|jb) / (e_i + e_j - e_a - e_b)
      !!
      !! The checkpoint before any amplitude equation runs: this must reproduce
      !! `run_libcint_mp2` exactly, and if it does then the block layout, the
      !! denominators and the frozen-core offset are all right.
      type(rcc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: eps_o(:), eps_v(:)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: t2(:, :, :, :)
      real(dp), intent(out) :: e_mp2

      integer :: i, j, a, b
      real(dp) :: d, dummy

      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  d = eps_o(i) + eps_o(j) - eps_v(a) - eps_v(b)
                  t2(i, j, a, b) = eris%ovov(i, a, j, b)/d
               end do
            end do
         end do
      end do

      block
         real(dp), allocatable :: t1zero(:, :)
         allocate (t1zero(no, nv))
         t1zero = 0.0_dp
         call rccsd_correlation_energy(eris, t1zero, t2, no, nv, dummy, e_mp2)
      end block
   end subroutine mp2_guess

   subroutine add_symmetrised(t2new, tmp, no, nv, factor)
      !! t2new(i,j,a,b) += factor * [tmp(i,j,a,b) + tmp(j,i,b,a)]
      !!
      !! The `tmp + tmp.transpose(1,0,3,2)` that appears eight times in the
      !! reference's T2 equation. Written once: the permutation is the same
      !! every time and typing it eight times is eight chances to swap one pair
      !! and not the other, which is a symmetric-looking error that still
      !! converges.
      real(dp), intent(inout) :: t2new(:, :, :, :)
      real(dp), intent(in) :: tmp(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: factor

      integer :: i, j, a, b

      !$omp parallel do default(none) shared(t2new, tmp, no, nv, factor) &
      !$omp    private(i, j, a, b) schedule(static)
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  t2new(i, j, a, b) = t2new(i, j, a, b) &
                                      + factor*(tmp(i, j, a, b) + tmp(j, i, b, a))
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine add_symmetrised

   subroutine rccsd_iteration(eris, eps_o, eps_v, no, nv, t1, t2, t1n, t2n)
      !! One amplitude update: Hirata Eqs. (35) and (36)
      !!
      !! Correctness before speed, deliberately. The contractions below are
      !! written as loops that read as the einsum strings quoted above them,
      !! except the particle-particle ladder, which is a gemm because it is both
      !! the asymptotically dominant term and the one whose memory this module
      !! exists to fix. Turning the O(n_o^3 n_v^3) ring terms into gemms is the
      !! obvious next step and is deliberately a separate change: it is a
      !! rearrangement of arithmetic that must not move a digit, and that is
      !! only checkable against a version already known to be right.
      type(rcc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: eps_o(:), eps_v(:)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      real(dp), intent(out) :: t1n(:, :), t2n(:, :, :, :)

      real(dp), allocatable :: tau(:, :, :, :), tmp(:, :, :, :)
      real(dp), allocatable :: fki(:, :), fac(:, :), fkc(:, :), lki(:, :), lac(:, :)
      real(dp), allocatable :: woooo(:, :, :, :), wvoov(:, :, :, :), wvovo(:, :, :, :)
      real(dp), allocatable :: tmp2(:, :, :, :), tmp2b(:, :, :, :)
      integer :: i, j, k, l, a, b, c, d
      real(dp) :: acc, den

      allocate (tau(no, no, nv, nv))
      call build_tau(t1, t2, no, nv, tau)

      allocate (fki(no, no), fac(nv, nv), fkc(no, nv), lki(no, no), lac(nv, nv))
      call cc_foo(eris, t1, tau, no, nv, fki)
      call cc_fvv(eris, t1, tau, no, nv, fac)
      call cc_fov(eris, t1, no, nv, fkc)
      call cc_loo(eris, fki, t1, no, nv, lki)
      call cc_lvv(eris, fac, t1, no, nv, lac)

      ! ---- T1, Eq. (35) ----------------------------------------------------
      t1n = 0.0_dp

      ! 'ac,ic->ia' and '-ki,ka->ia'
      do c = 1, nv
         do a = 1, nv
            do i = 1, no
               t1n(i, a) = t1n(i, a) + fac(a, c)*t1(i, c)
            end do
         end do
      end do
      do k = 1, no
         do i = 1, no
            do a = 1, nv
               t1n(i, a) = t1n(i, a) - fki(k, i)*t1(k, a)
            end do
         end do
      end do

      ! '2 kc,kica->ia', '-kc,ikca->ia', 'kc,ic,ka->ia'
      do c = 1, nv
         do k = 1, no
            do a = 1, nv
               do i = 1, no
                  t1n(i, a) = t1n(i, a) + fkc(k, c)*(2.0_dp*t2(k, i, c, a) - t2(i, k, c, a)) &
                              + fkc(k, c)*t1(i, c)*t1(k, a)
               end do
            end do
         end do
      end do

      ! '2 kcai,kc->ia' with (kc|ai) = ovov(k,c,i,a); '-kiac,kc->ia'
      do c = 1, nv
         do k = 1, no
            do a = 1, nv
               do i = 1, no
                  t1n(i, a) = t1n(i, a) &
                              + (2.0_dp*eris%ovov(k, c, i, a) - eris%oovv(k, i, a, c))*t1(k, c)
               end do
            end do
         end do
      end do

      ! '2 kdac,ikcd->ia' - 'kcad,ikcd->ia', and the same pair with t2 -> t1 t1
      ! `a` outermost, so each thread owns its own columns of t1n.
      !$omp parallel do default(none) shared(t1n, eris, t1, t2, no, nv) &
      !$omp    private(a, d, c, k, i) schedule(static)
      do a = 1, nv
         do d = 1, nv
            do c = 1, nv
               do k = 1, no
                  do i = 1, no
                     t1n(i, a) = t1n(i, a) &
                                 + 2.0_dp*eris%ovvv(k, d, a, c)*(t2(i, k, c, d) + t1(k, d)*t1(i, c)) &
                                 - eris%ovvv(k, c, a, d)*(t2(i, k, c, d) + t1(k, d)*t1(i, c))
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! '-2 lcki,klac->ia' + 'kcli,klac->ia' with (lc|ki) = ooov(k,i,l,c),
      ! (kc|li) = ooov(l,i,k,c); and the same pair with t2 -> t1 t1
      !$omp parallel do default(none) shared(t1n, eris, t1, t2, no, nv) &
      !$omp    private(a, c, l, k, i) schedule(static)
      do a = 1, nv
         do c = 1, nv
            do l = 1, no
               do k = 1, no
                  do i = 1, no
                     t1n(i, a) = t1n(i, a) &
                                 - 2.0_dp*eris%ooov(k, i, l, c)*(t2(k, l, a, c) + t1(l, c)*t1(k, a)) &
                                 + eris%ooov(l, i, k, c)*(t2(k, l, a, c) + t1(l, c)*t1(k, a))
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do

      do a = 1, nv
         do i = 1, no
            t1n(i, a) = t1n(i, a)/(eps_o(i) - eps_v(a))
         end do
      end do

      ! ---- T2, Eq. (36) ----------------------------------------------------
      allocate (tmp(no, no, nv, nv))
      t2n = 0.0_dp

      ! t2new += ovov(i,a,j,b)   ['eris.ovov.transpose(0,2,1,3)']
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  t2n(i, j, a, b) = eris%ovov(i, a, j, b)
               end do
            end do
         end do
      end do

      ! tmp2(a,b,i,c) = ovvv(i,a,c,b) - sum_k oovv(k,i,b,c) t1(k,a)
      ! tmp(i,j,a,b)  = sum_c tmp2(a,b,i,c) t1(j,c);  t2new += tmp + P(ij,ab) tmp
      allocate (tmp2(nv, nv, no, nv))
      !$omp parallel do default(none) shared(tmp2, eris, t1, no, nv) &
      !$omp    private(c, i, b, a, k, acc) schedule(static)
      do c = 1, nv
         do i = 1, no
            do b = 1, nv
               do a = 1, nv
                  acc = eris%ovvv(i, a, c, b)
                  do k = 1, no
                     acc = acc - eris%oovv(k, i, b, c)*t1(k, a)
                  end do
                  tmp2(a, b, i, c) = acc
               end do
            end do
         end do
      end do
      !$omp end parallel do
      tmp = 0.0_dp
      ! `b` outermost: tmp(:,:,:,b) is one thread's alone.
      !$omp parallel do default(none) shared(tmp, tmp2, t1, no, nv) &
      !$omp    private(b, c, a, j, i) schedule(static)
      do b = 1, nv
         do c = 1, nv
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     tmp(i, j, a, b) = tmp(i, j, a, b) + tmp2(a, b, i, c)*t1(j, c)
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do
      call add_symmetrised(t2n, tmp, no, nv, 1.0_dp)
      deallocate (tmp2)

      ! tmp2b(a,k,i,j) = sum_c ovov(k,c,i,a) t1(j,c) + ooov(j,k,i,a)
      ! tmp(i,j,a,b)   = sum_k tmp2b(a,k,i,j) t1(k,b);  t2new -= tmp + P tmp
      allocate (tmp2b(nv, no, no, no))
      do j = 1, no
         do i = 1, no
            do k = 1, no
               do a = 1, nv
                  acc = eris%ooov(j, k, i, a)
                  do c = 1, nv
                     acc = acc + eris%ovov(k, c, i, a)*t1(j, c)
                  end do
                  tmp2b(a, k, i, j) = acc
               end do
            end do
         end do
      end do
      tmp = 0.0_dp
      do k = 1, no
         do b = 1, nv
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     tmp(i, j, a, b) = tmp(i, j, a, b) + tmp2b(a, k, i, j)*t1(k, b)
                  end do
               end do
            end do
         end do
      end do
      call add_symmetrised(t2n, tmp, no, nv, -1.0_dp)
      deallocate (tmp2b)

      ! 'klij,klab->ijab' with the four-occupied intermediate
      allocate (woooo(no, no, no, no))
      call cc_woooo(eris, t1, tau, no, nv, woooo)
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  acc = 0.0_dp
                  do l = 1, no
                     do k = 1, no
                        acc = acc + woooo(k, l, i, j)*tau(k, l, a, b)
                     end do
                  end do
                  t2n(i, j, a, b) = t2n(i, j, a, b) + acc
               end do
            end do
         end do
      end do
      deallocate (woooo)

      ! 'abcd,ijcd->ijab' -- the particle-particle ladder, never held whole
      call particle_ladder(eris, t1, tau, no, nv, t2n)

      ! 'ac,ijcb->ijab' and '-ki,kjab->ijab', each symmetrised
      tmp = 0.0_dp
      !$omp parallel do default(none) shared(tmp, lac, t2, no, nv) &
      !$omp    private(b, c, a, j, i) schedule(static)
      do b = 1, nv
         do c = 1, nv
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     tmp(i, j, a, b) = tmp(i, j, a, b) + lac(a, c)*t2(i, j, c, b)
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do
      call add_symmetrised(t2n, tmp, no, nv, 1.0_dp)

      tmp = 0.0_dp
      !$omp parallel do default(none) shared(tmp, lki, t2, no, nv) &
      !$omp    private(b, k, a, j, i) schedule(static)
      do b = 1, nv
         do k = 1, no
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     tmp(i, j, a, b) = tmp(i, j, a, b) + lki(k, i)*t2(k, j, a, b)
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do
      call add_symmetrised(t2n, tmp, no, nv, -1.0_dp)

      ! The three ring terms, together
      !
      ! Each contracts over (k,c) and each is O(n_occ^3 n_vir^3), which makes
      ! them jointly the largest thing in this routine after the ladder. Written
      ! as loops they read straight off the reference; written as gemms they
      ! need the operands laid out with the contracted pair as one compound
      ! index, which is what the packing below does. The packs are all
      ! O(n_occ^2 n_vir^2) -- a factor of n_occ n_vir cheaper than the products
      ! they feed -- so the rearrangement is free in the only sense that
      ! matters.
      !
      ! Three of the five operands are shared, which is why they are done as a
      ! group rather than one at a time:
      !
      !     X1 = 2 X2 - X3        so only X2 and X3 are built
      !     term C reuses Y2       with its free virtual named a instead of b
      allocate (wvoov(nv, no, no, nv), wvovo(nv, no, nv, no))
      call cc_wvoov(eris, t1, t2, no, nv, wvoov)
      call cc_wvovo(eris, t1, t2, no, nv, wvovo)

      block
         real(dp), allocatable :: x1(:, :), x2(:, :), x3(:, :)
         real(dp), allocatable :: y1(:, :), y2(:, :), r(:, :)
         integer :: nov, ai, kc, jb

         nov = no*nv
         allocate (x1(nov, nov), x2(nov, nov), x3(nov, nov))
         allocate (y1(nov, nov), y2(nov, nov), r(nov, nov))

         ! X2(ai,kc) = Wvoov(a,k,i,c);  X3(ai,kc) = Wvovo(a,k,c,i)
         do c = 1, nv
            do k = 1, no
               kc = (c - 1)*no + k
               do i = 1, no
                  do a = 1, nv
                     ai = (i - 1)*nv + a
                     x2(ai, kc) = wvoov(a, k, i, c)
                     x3(ai, kc) = wvovo(a, k, c, i)
                  end do
               end do
            end do
         end do
         x1 = 2.0_dp*x2 - x3

         ! Y1(kc,jb) = t2(k,j,c,b);  Y2(kc,jb) = t2(k,j,b,c)
         do b = 1, nv
            do j = 1, no
               jb = (b - 1)*no + j
               do c = 1, nv
                  do k = 1, no
                     kc = (c - 1)*no + k
                     y1(kc, jb) = t2(k, j, c, b)
                     y2(kc, jb) = t2(k, j, b, c)
                  end do
               end do
            end do
         end do

         ! '2 akic,kjcb->ijab' - 'akci,kjcb->ijab'
         call gemm_over_columns(x1, y1, r)
         tmp = 0.0_dp
         call scatter_ring(r, no, nv, .false., tmp)
         call add_symmetrised(t2n, tmp, no, nv, 1.0_dp)

         ! '-akic,kjbc->ijab'
         call gemm_over_columns(x2, y2, r)
         tmp = 0.0_dp
         call scatter_ring(r, no, nv, .false., tmp)
         call add_symmetrised(t2n, tmp, no, nv, -1.0_dp)

         ! '-bkci,kjac->ijab'. The only one whose free indices come out paired
         ! the other way round, (b,i) with (j,a), which `swapped` says.
         call gemm_over_columns(x3, y2, r)
         tmp = 0.0_dp
         call scatter_ring(r, no, nv, .true., tmp)
         call add_symmetrised(t2n, tmp, no, nv, -1.0_dp)

         deallocate (x1, x2, x3, y1, y2, r)
      end block
      deallocate (wvoov, wvovo)

      !$omp parallel do default(none) shared(t2n, eps_o, eps_v, no, nv) &
      !$omp    private(b, a, j, i, den) schedule(static)
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  den = eps_o(i) + eps_o(j) - eps_v(a) - eps_v(b)
                  t2n(i, j, a, b) = t2n(i, j, a, b)/den
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine rccsd_iteration

   subroutine accumulate_wvoov(r, no, nv, w)
      !! w(a,k,i,c) += R((k,c),(i,a))
      real(dp), intent(in) :: r(:, :)
      integer, intent(in) :: no, nv
      real(dp), intent(inout) :: w(:, :, :, :)

      integer :: a, k, i, c

      do c = 1, nv
         do i = 1, no
            do k = 1, no
               do a = 1, nv
                  w(a, k, i, c) = w(a, k, i, c) + r((c - 1)*no + k, (a - 1)*no + i)
               end do
            end do
         end do
      end do
   end subroutine accumulate_wvoov

   subroutine accumulate_wvovo(r, no, nv, w)
      !! w(a,k,c,i) += R((k,c),(i,a))
      !!
      !! Same product, different destination ordering -- the two W intermediates
      !! differ in exactly the placement of their occupied and virtual pair,
      !! which is what their names say and the only thing separating them.
      real(dp), intent(in) :: r(:, :)
      integer, intent(in) :: no, nv
      real(dp), intent(inout) :: w(:, :, :, :)

      integer :: a, k, c, i

      do i = 1, no
         do c = 1, nv
            do k = 1, no
               do a = 1, nv
                  w(a, k, c, i) = w(a, k, c, i) + r((c - 1)*no + k, (a - 1)*no + i)
               end do
            end do
         end do
      end do
   end subroutine accumulate_wvovo

   subroutine scatter_ring(r, no, nv, swapped, tmp)
      !! Add a ring term's gemm result back into (i,j,a,b) order
      !!
      !! The product comes out indexed by the two compound indices the
      !! contraction left free. Two of the three terms leave (a,i) against
      !! (j,b); the third leaves (b,i) against (j,a), which `swapped` selects.
      !! Getting that one wrong gives a result that is still symmetric under the
      !! i<->j, a<->b permutation applied afterwards, and still the right
      !! magnitude.
      real(dp), intent(in) :: r(:, :)
      integer, intent(in) :: no, nv
      logical, intent(in) :: swapped
      real(dp), intent(inout) :: tmp(:, :, :, :)

      integer :: i, j, a, b, row, col

      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  if (swapped) then
                     row = (i - 1)*nv + b
                     col = (a - 1)*no + j
                  else
                     row = (i - 1)*nv + a
                     col = (b - 1)*no + j
                  end if
                  tmp(i, j, a, b) = tmp(i, j, a, b) + r(row, col)
               end do
            end do
         end do
      end do
   end subroutine scatter_ring

   subroutine gemm_over_columns(a, b, c)
      !! C = A B, with C's columns split across threads
      !!
      !! For the products whose result is wide but not tall. Every operand slice
      !! here is a column range, so it is contiguous and nothing is copied:
      !! thread t computes C(:, s0:s1) from all of A and B(:, s0:s1), and the
      !! ranges are disjoint, so there is no reduction and no ordering question.
      !!
      !! This exists because a threaded BLAS does not help at these shapes. The
      !! ring and intermediate products are n_occ n_vir on a side -- a few
      !! hundred -- and measured with a threaded MKL the whole iteration went
      !! from 0.073 s to 0.059 s and stopped improving past four threads. The
      !! parallelism has to come from outside the call.
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), intent(inout) :: c(:, :)

      integer :: n, nchunk, width, chunk, s0, s1

      n = size(c, 2)
      nchunk = omp_get_max_threads()
      if (nchunk < 1) nchunk = 1
      if (nchunk > n) nchunk = n
      width = (n + nchunk - 1)/nchunk

      !$omp parallel do default(none) &
      !$omp    shared(a, b, c, n, nchunk, width) &
      !$omp    private(chunk, s0, s1) schedule(static)
      do chunk = 1, nchunk
         s0 = (chunk - 1)*width + 1
         s1 = min(chunk*width, n)
         if (s0 > s1) cycle
         call pic_gemm(a, b(:, s0:s1), c(:, s0:s1))
      end do
      !$omp end parallel do
   end subroutine gemm_over_columns

   subroutine ladder_accumulate(no2, nv2, nb, tau_cols, wblk, t2n)
      !! t2n(ij, ab) += sum_cd tau(ij, cd) W(ab, cd), over one batch of (cd)
      !!
      !! Explicit-shape dummies so the caller's rank-four arrays arrive as the
      !! matrices their memory already is, by sequence association -- the same
      !! device `ladder_t1_block` uses next door, and for the same reason: an
      !! assumed-shape dummy cannot take a block that way, and reshaping a
      !! n_occ^2 by n_vir^2 array every batch would cost more than the gemm.
      integer, intent(in) :: no2, nv2, nb
      real(dp), intent(in) :: tau_cols(no2, nb)   !! this batch's columns of tau
      real(dp), intent(in) :: wblk(nv2, nb)       !! Wvvvv for the same columns
      real(dp), intent(inout) :: t2n(no2, nv2)

      integer :: chunk, nchunk, width, r0, r1

      ! Split by output column, so each thread owns a disjoint column range of
      ! `t2n` and no reduction is needed. Threaded here rather than inside the
      ! product because the product is n_occ^2 rows tall -- sixty-four of them
      ! at ordinary sizes -- and a BLAS cannot usefully divide that, whatever it
      ! is told about threads.
      !
      ! Slicing `wblk` by row costs a copy, which the compiler makes. It is
      ! n_vir^2 by nb elements across all threads, against the n_occ^2 times
      ! that in flops, so it is paid back sixty-four times over.
      nchunk = omp_get_max_threads()
      if (nchunk < 1) nchunk = 1
      if (nchunk > nv2) nchunk = nv2
      width = (nv2 + nchunk - 1)/nchunk

      !$omp parallel do default(none) &
      !$omp    shared(tau_cols, wblk, t2n, no2, nv2, nb, nchunk, width) &
      !$omp    private(chunk, r0, r1) schedule(static)
      do chunk = 1, nchunk
         r0 = (chunk - 1)*width + 1
         r1 = min(chunk*width, nv2)
         if (r0 > r1) cycle
         call pic_gemm(tau_cols, wblk(r0:r1, 1:nb), t2n(:, r0:r1), &
                       transb="T", beta=1.0_dp)
      end do
      !$omp end parallel do
   end subroutine ladder_accumulate

   pure function ri_ladder_prefers_direct(no, nv, naux) result(direct)
      !! Whether to contract the fitted tensor against tau without forming (ac|bd)
      !!
      !! Two ways to reach the same sum, and neither wins everywhere:
      !!
      !!   assemble   naux n_vir^4 to build (ac|bd), then n_occ^2 n_vir^4 to
      !!              contract it -- what the conventional path does too, with
      !!              the build replacing a read
      !!   direct     2 naux n_vir^3 n_occ^2, never forming (ac|bd) at all
      !!
      !! Direct wins when n_vir exceeds roughly 2 n_occ^2, which is a small
      !! molecule in a large basis; assemble wins with many occupied orbitals,
      !! which is where anyone actually reaches for density fitting. Water in
      !! cc-pVTZ is 1.0 GFlop against 0.6; benzene in cc-pVDZ is 63 against 284
      !! the other way. Choosing on measured counts rather than picking one is
      !! the only honest option, because the ratio spans two orders of magnitude
      !! across ordinary cases.
      !!
      !! The comparison is on operation counts only. Both forms are gemms of
      !! comparable shape, so the counts predict the ordering; what they do not
      !! predict is the constant, which is why the crossover is written as the
      !! full expression rather than the tidy `nv > 2 no^2` it reduces to when
      !! naux dominates n_occ^2.
      integer, intent(in) :: no, nv, naux
      logical :: direct

      real(dp) :: assemble_ops, direct_ops

      assemble_ops = (real(naux, dp) + real(no, dp)**2)*real(nv, dp)**4
      direct_ops = 2.0_dp*real(naux, dp)*real(nv, dp)**3*real(no, dp)**2
      direct = (direct_ops < assemble_ops)
   end function ri_ladder_prefers_direct

   subroutine ri_ladder_pair(nv, gp, tauij, work, accab)
      !! accab(a,b) += sum_cd G(a,c) tau(c,d) G(b,d), one auxiliary function
      !!
      !! `gp` is one column of `b_vv` seen as the n_vir by n_vir matrix its
      !! memory already is: the compound index runs virtual-fastest, so column P
      !! laid out as (a,c) is exactly B^P_ac with nothing moved.
      integer, intent(in) :: nv
      real(dp), intent(in) :: gp(nv, nv)      !! B^P_ac
      real(dp), intent(in) :: tauij(nv, nv)   !! tau(i,j,c,d) at fixed (i,j)
      real(dp), intent(inout) :: work(nv, nv)
      real(dp), intent(inout) :: accab(nv, nv)

      call pic_gemm(gp, tauij, work)
      call pic_gemm(work, gp, accab, transb="T", beta=1.0_dp)
   end subroutine ri_ladder_pair

   subroutine ri_ladder_direct(eris, tau, no, nv, naux, t2n)
      !! The fitted ladder without ever forming (ac|bd)
      !!
      !!     sum_cd (ac|bd) tau(i,j,c,d) = sum_P [ G_P tau_ij G_P^T ](a,b)
      !!
      !! which is two n_vir-cubed products per auxiliary function per occupied
      !! pair, against the n_vir^4 of assembling the integral first. Nothing
      !! held here exceeds n_vir^2, and `G_P` is not held at all -- it is a
      !! column of `b_vv` read in place.
      type(rcc_eris_t), intent(in) :: eris
      integer, intent(in) :: no, nv, naux
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp), intent(inout) :: t2n(no, no, nv, nv)

      real(dp), allocatable :: tauij(:, :), work(:, :), accab(:, :)
      integer :: i, j, a, b, c, d, p

      ! Each occupied pair writes its own (a,b) block of `t2n` and reads
      ! nothing another pair writes, so the pair loop threads with no reduction
      ! and no ordering question. The products inside are n_vir cubed, which is
      ! a shape a threaded BLAS could split -- but there are naux of them per
      ! pair and n_occ^2 pairs, so the outer loop has far more parallelism to
      ! give and costs nothing to take it from.
      !$omp parallel default(none) &
      !$omp    shared(eris, tau, t2n, no, nv, naux) &
      !$omp    private(i, j, a, b, c, d, p, tauij, work, accab)
      allocate (tauij(nv, nv), work(nv, nv), accab(nv, nv))

      !$omp do collapse(2) schedule(static)
      do j = 1, no
         do i = 1, no
            do d = 1, nv
               do c = 1, nv
                  tauij(c, d) = tau(i, j, c, d)
               end do
            end do

            accab = 0.0_dp
            do p = 1, naux
               call ri_ladder_pair(nv, eris%b_vv(1, p), tauij, work, accab)
            end do

            do b = 1, nv
               do a = 1, nv
                  t2n(i, j, a, b) = t2n(i, j, a, b) + accab(a, b)
               end do
            end do
         end do
      end do
      !$omp end do

      deallocate (tauij, work, accab)
      !$omp end parallel
   end subroutine ri_ladder_direct

   subroutine ladder_t1_dressing(eris, t1, tau, no, nv, t2n)
      !! The t1 half of the particle-particle ladder, without forming it
      !!
      !! Wvvvv carries two singles terms beside `(ac|bd)`, and contracting them
      !! with tau in the order they are written costs O(n_occ n_vir^4) -- the
      !! same order as the ladder itself, and measured at 54% of an iteration
      !! once everything around it had been turned into matrix products.
      !!
      !! Reassociating fixes it. Summing over (c,d) first,
      !!
      !!     -sum_cd [sum_k ovvv(k,d,a,c) t1(k,b)] tau(i,j,c,d)
      !!         = -sum_k t1(k,b) Z1(k,a,i,j),  Z1 = sum_cd ovvv(k,d,a,c) tau(i,j,c,d)
      !!
      !! and likewise for the second term with `ovvv(k,c,b,d)`. Z is
      !! O(n_occ^3 n_vir^3) to build and the t1 contraction that follows is
      !! smaller again, so the n_vir^4 disappears entirely.
      !!
      !! Held one virtual at a time, so nothing here is larger than
      !! `n_occ n_vir^2` -- the alternative, packing both permuted copies of
      !! `ovvv` whole, would have tripled the largest block this module holds
      !! and undone the reason it exists.
      type(rcc_eris_t), intent(in) :: eris
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv)
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp), intent(inout) :: t2n(no, no, nv, nv)

      real(dp), allocatable :: taut(:, :), p(:, :), z(:, :), cmat(:, :)
      integer :: i, j, a, b, c, d, k, v, cd, ij

      allocate (taut(nv*nv, no*no))

      ! tau with the virtual pair leading, once per iteration.
      do d = 1, nv
         do c = 1, nv
            cd = (d - 1)*nv + c
            do j = 1, no
               do i = 1, no
                  taut(cd, (j - 1)*no + i) = tau(i, j, c, d)
               end do
            end do
         end do
      end do

      ! Z1 and Z2 are threaded as two loops rather than one.
      !
      ! Within a single v they touch disjoint parts of `t2n` -- Z1 writes
      ! t2n(:,:,v,:) and Z2 writes t2n(:,:,:,v) -- but across two threads
      ! holding v1 and v2 those two writes collide at (a,b) = (v1,v2). Split in
      ! two, each loop's writes are disjoint by construction and the collision
      ! cannot happen. The scratch is per thread; only `taut` is shared, and it
      ! is read-only here.
      !$omp parallel default(none) &
      !$omp    shared(eris, t1, t2n, taut, no, nv) &
      !$omp    private(v, i, j, b, c, d, k, cd, p, z, cmat)
      allocate (p(no, nv*nv), z(no, no*no), cmat(no*no, nv))

      !$omp do schedule(static)
      do v = 1, nv
         ! ---- Z1: the free virtual is the first of the pair ----------------
         do c = 1, nv
            do d = 1, nv
               cd = (d - 1)*nv + c
               do k = 1, no
                  p(k, cd) = eris%ovvv(k, d, v, c)
               end do
            end do
         end do
         call pic_gemm(p, taut, z)
         call pic_gemm(z, t1, cmat, transa="T")
         do b = 1, nv
            do j = 1, no
               do i = 1, no
                  t2n(i, j, v, b) = t2n(i, j, v, b) - cmat((j - 1)*no + i, b)
               end do
            end do
         end do
      end do
      !$omp end do

      !$omp do schedule(static)
      do v = 1, nv
         ! ---- Z2: the free virtual is the second ---------------------------
         do c = 1, nv
            do d = 1, nv
               cd = (d - 1)*nv + c
               do k = 1, no
                  p(k, cd) = eris%ovvv(k, c, v, d)
               end do
            end do
         end do
         call pic_gemm(p, taut, z)
         call pic_gemm(z, t1, cmat, transa="T")
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  t2n(i, j, a, v) = t2n(i, j, a, v) - cmat((j - 1)*no + i, a)
               end do
            end do
         end do
      end do
      !$omp end do

      deallocate (p, z, cmat)
      !$omp end parallel

      deallocate (taut)
   end subroutine ladder_t1_dressing

   subroutine particle_ladder(eris, t1, tau, no, nv, t2n)
      !! t2new(i,j,a,b) += sum_cd Wvvvv(a,b,c,d) tau(i,j,c,d), never holding Wvvvv
      !!
      !!     Wvvvv(a,b,c,d) = (ac|bd) - sum_k ovvv(k,d,a,c) t1(k,b)
      !!                              - sum_k ovvv(k,c,b,d) t1(k,a)
      !!
      !! Only the `(ac|bd)` part is batched here. The two t1 terms are handled
      !! by `ladder_t1_dressing`, which never forms them: contracted against tau
      !! first they cost O(n_occ^3 n_vir^3) instead of the O(n_occ n_vir^4) they
      !! cost inside this loop, and that difference was measured at just over
      !! half of one iteration.
      !!
      !! `n_vir^4` is the largest thing coupled cluster asks for and it is read
      !! exactly once, here, which is what makes batching it possible rather
      !! than merely desirable: a block of `(cd)` columns is built, contracted
      !! and dropped. Same argument, same structure and the same batch size as
      !! `particle_ladder` in the spin-orbital module -- the difference is that
      !! every index here is already spatial, so there is no spin map to hoist.
      !!
      !! **Where (ac|bd) comes from.** Exactly one of `vvvv` and `b_vv` is
      !! allocated. `vvvv` is the conventional path's, read at permuted indices.
      !! `b_vv` is the fitted three-index block laid out `(ab, P)`, and for a
      !! fixed `(c,d)` the matrix `(ac|bd)` over `a` and `b` is one gemm over the
      !! auxiliary index between two contiguous row ranges of it -- which is why
      !! the fitted path never allocates anything of fourth order in the
      !! virtuals at all.
      type(rcc_eris_t), intent(in) :: eris
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(no, nv)
      ! Explicit shape, not assumed: the accumulation below hands a column range
      ! of `tau` and the whole of `t2n` to a rank-two dummy by sequence
      ! association, and only an explicit-shape or allocatable actual may be
      ! sliced that way.
      real(dp), intent(in) :: tau(no, no, nv, nv)
      real(dp), intent(inout) :: t2n(no, no, nv, nv)

      real(dp), allocatable :: wblk(:, :, :), gvv(:, :)
      integer :: nv2, cd0, cd1, nb, col, cd, a, b, c, d, i, j, k, naux
      logical :: fitted
      real(dp) :: acc

      nv2 = nv*nv
      fitted = allocated(eris%b_vv)
      if (fitted) then
         naux = size(eris%b_vv, 2)
         ! Which way round to do the fitted ladder. See `ri_ladder_prefers_direct`.
         if (ri_ladder_prefers_direct(no, nv, naux)) then
            call ri_ladder_direct(eris, tau, no, nv, naux, t2n)
            call ladder_t1_dressing(eris, t1, tau, no, nv, t2n)
            return
         end if
      end if

      allocate (wblk(nv, nv, LADDER_BATCH))

      do cd0 = 1, nv2, LADDER_BATCH
         cd1 = min(cd0 + LADDER_BATCH - 1, nv2)
         nb = cd1 - cd0 + 1

         ! One column of the batch is one (c,d), and the columns are
         ! independent: each writes its own page of `wblk` and reads nothing
         ! another writes.
         !
         ! Threading the loop rather than the products inside it, because of
         ! their shape. On the fitted path this is most of an iteration --
         ! naux n_vir^4 of assembling (ac|bd) -- but it arrives as n_vir^2
         ! separate products with a n_vir by n_vir result, which is far too
         ! small for a threaded BLAS to split and far too many to leave serial.
         ! Measured: a threaded MKL took the fitted iteration from 0.116 s to
         ! 0.101 s and stopped improving past four threads.
         !$omp parallel default(none) &
         !$omp    shared(eris, wblk, nv, nb, cd0, naux, fitted) &
         !$omp    private(col, cd, a, b, c, d, gvv)
         if (fitted) allocate (gvv(nv, nv))
         !$omp do schedule(static)
         do col = 1, nb
            cd = cd0 + col - 1
            c = mod(cd - 1, nv) + 1
            d = (cd - 1)/nv + 1

            if (fitted) then
               ! (ac|bd) over a and b: rows (c-1)nv+1 .. c nv of b_vv against
               ! rows (d-1)nv+1 .. d nv, contracted over the auxiliary index.
               call pic_gemm(eris%b_vv((c - 1)*nv + 1:c*nv, 1:naux), &
                             eris%b_vv((d - 1)*nv + 1:d*nv, 1:naux), gvv, transb="T")
               do b = 1, nv
                  do a = 1, nv
                     wblk(a, b, col) = gvv(a, b)
                  end do
               end do
            else
               do b = 1, nv
                  do a = 1, nv
                     wblk(a, b, col) = eris%vvvv(a, c, b, d)
                  end do
               end do
            end if

         end do
         !$omp end do
         if (allocated(gvv)) deallocate (gvv)
         !$omp end parallel

         ! One gemm for the whole batch. `tau` is (no^2, nv^2) as it lies in
         ! memory and this batch is a contiguous column range of it, starting at
         ! the (c,d) that opens the batch; `wblk` is (nv^2, nb) likewise. So the
         ! accumulation is C += A B^T with nothing copied and nothing reshaped,
         ! which is the whole reason the batch runs over the compound (cd)
         ! rather than over c and d separately.
         c = mod(cd0 - 1, nv) + 1
         d = (cd0 - 1)/nv + 1
         call ladder_accumulate(no*no, nv*nv, nb, tau(1, 1, c, d), wblk, t2n(1, 1, 1, 1))
      end do

      deallocate (wblk)

      call ladder_t1_dressing(eris, t1, tau, no, nv, t2n)
   end subroutine particle_ladder

   !===========================================================================
   ! Driver
   !===========================================================================

   subroutine run_libcint_rccsd(mol, coeff, orbital_energies, n_occ, frozen, &
                                max_iter, energy_tol, want_triples, verbose, &
                                result, error, diis_vectors, aux)
      !! Drive spin-adapted CCSD to convergence
      !!
      !! Signature deliberately identical to `run_libcint_ccsd` so the bridge can
      !! choose between them without knowing which it called, and so the two can
      !! be run back to back on one set of orbitals -- which is what the identity
      !! test does.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ                  !! Spatial occupied count
      integer, intent(in) :: frozen
      integer, intent(in) :: max_iter
      real(dp), intent(in) :: energy_tol
      logical, intent(in) :: want_triples
      logical, intent(in) :: verbose
      type(rcc_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: diis_vectors
      type(libcint_molecule_t), intent(in), optional :: aux

      type(rcc_eris_t) :: eris
      real(dp), allocatable :: c_act(:, :), eps_o(:), eps_v(:)
      real(dp), allocatable :: t1(:, :), t2(:, :, :, :), t1n(:, :), t2n(:, :, :, :)
      real(dp), allocatable :: flat(:), err_flat(:)
      type(diis_state_t) :: diis
      logical :: extrapolated
      integer :: n_ao, n_mo, n_act, no, nv, iter, diis_size, i, a
      real(dp) :: e_corr, e_old, de, t_iter
      character(len=MAX_LINE_LENGTH) :: line
      type(timing_report_t) :: clk

      n_ao = size(coeff, 1)
      n_mo = size(coeff, 2)

      if (frozen < 0 .or. frozen >= n_occ) then
         call error%set(ERROR_VALIDATION, "RCCSD: the frozen core count must leave at "// &
                        "least one occupied orbital")
         return
      end if
      n_act = n_mo - frozen
      no = n_occ - frozen
      nv = n_mo - n_occ
      if (nv < 1) then
         call error%set(ERROR_VALIDATION, "RCCSD: no virtual orbitals -- the basis is "// &
                        "saturated by the occupied space and there is nothing to excite into")
         return
      end if

      diis_size = 8
      if (present(diis_vectors)) diis_size = diis_vectors

      if (verbose) then
         write (line, "(a,i0,a,i0,a,i0,a)") "  spin-adapted coupled cluster: ", n_act, &
            " spatial orbitals, ", no, " occupied, ", nv, " virtual"
         call logger%info(trim(line))
         if (frozen > 0) then
            write (line, "(a,i0,a)") "  frozen core: ", frozen, " spatial orbitals"
            call logger%info(trim(line))
         end if
         ! Against what the spin-orbital path would have taken for the same
         ! system, because that difference is the entire reason this exists.
         write (line, "(a,f0.1,a,f0.1,a)") "  integrals: ", &
            rcc_megabytes(no, nv,.not. present(aux)), " MB in blocks, against ", &
            rcc_megabytes(2*no, 2*nv,.not. present(aux)), " MB in spin orbitals"
         call logger%info(trim(line))
      end if

      call clk%start()
      allocate (c_act(n_ao, n_act))
      c_act = coeff(:, frozen + 1:n_mo)
      if (present(aux)) then
         ! One stage: the fitted path never forms an AO tensor to separate out.
         call clk%begin("B tensor (3c/2c fit)")
         call build_rcc_eris_fitted(mol, aux, c_act(:, 1:no), c_act(:, no + 1:n_act), &
                                    eris, error)
         if (error%has_error()) return
         call clk%lap()
      else
         ! Two stages, lapped inside: the AO build and the transform.
         call build_rcc_eris_conventional(mol, c_act(:, 1:no), c_act(:, no + 1:n_act), &
                                          eris, clk)
      end if
      deallocate (c_act)

      allocate (eps_o(no), eps_v(nv))
      do i = 1, no
         eps_o(i) = orbital_energies(frozen + i)
      end do
      do a = 1, nv
         eps_v(a) = orbital_energies(n_occ + a)
      end do

      ! ---- MP2, the checkpoint before any amplitude equation ----------------
      call clk%begin("MP2 amplitudes")
      allocate (t1(no, nv), t2(no, no, nv, nv))
      t1 = 0.0_dp
      call mp2_guess(eris, eps_o, eps_v, no, nv, t2, result%e_mp2)
      call clk%lap()
      if (verbose) then
         write (line, "(a,f20.12)") "  MP2 (spin adapted) = ", result%e_mp2
         call logger%info(trim(line))
      end if

      ! ---- CCSD -------------------------------------------------------------
      allocate (t1n(no, nv), t2n(no, no, nv, nv))
      allocate (flat(no*nv + no*no*nv*nv), err_flat(no*nv + no*no*nv*nv))
      call diis%init(diis_size, size(flat), size(err_flat))

      e_old = 0.0_dp
      result%converged = .false.
      call convergence_header(verbose, "RCCSD iterations", &
                              "    iter                 E_corr          dE   diis       time", 60)

      do iter = 1, max_iter
         t_iter = clk%seconds_of(STAGE_ITER)
         call rccsd_iteration(eris, eps_o, eps_v, no, nv, t1, t2, t1n, t2n)
         call clk%lap(STAGE_ITER)
         t_iter = clk%seconds_of(STAGE_ITER) - t_iter

         call pack_amplitudes(t1n, t2n, no, nv, flat)
         call pack_step(t1n, t2n, t1, t2, no, nv, err_flat)
         call diis%push(flat, err_flat)
         call diis%extrapolate(flat, extrapolated)
         if (extrapolated) call unpack_amplitudes(flat, no, nv, t1n, t2n)

         t1 = t1n
         t2 = t2n

         call rccsd_correlation_energy(eris, t1, t2, no, nv, result%e_singles, result%e_doubles)
         e_corr = result%e_singles + result%e_doubles
         de = abs(e_corr - e_old)
         if (verbose) then
            write (line, "(i8,f23.12,es12.3,i7,f9.2,a)") &
               iter, e_corr, de, diis%count(), t_iter, " s"
            call logger%info(trim(line))
         end if

         e_old = e_corr
         result%iterations = iter
         if (iter > 1 .and. de < energy_tol) then
            result%converged = .true.
            exit
         end if
      end do
      call diis%destroy()
      call convergence_footer(verbose, result%converged, result%iterations, "iterations", 50)

      if (.not. result%converged) then
         call error%set(ERROR_VALIDATION, "RCCSD did not converge")
         return
      end if

      ! (T) only on a converged reference: it is a correction evaluated once at
      ! the fixed point, and taking it on amplitudes that are still moving would
      ! give a number with no meaning rather than a slightly worse one.
      if (want_triples) then
         call clk%begin("(T) correction")
         call triples_correction(eris, eps_o, eps_v, t1, t2, no, nv, result%e_triples)
         call clk%lap()
         if (verbose) then
            write (line, "(a,f20.12)") "  (T)                = ", result%e_triples
            call logger%info(trim(line))
         end if
      end if

      result%e_correlation = result%e_singles + result%e_doubles + result%e_triples

      ! Same breakdown the spin-orbital driver prints, and under the same stage
      ! names, so the two can be laid side by side without translating.
      call clk%finish()
      call clk%report("RCCSD")
   end subroutine run_libcint_rccsd

   pure function int_text(n) result(text)
      !! An integer as a trimmed string, for message building
      integer, intent(in) :: n
      character(len=16) :: text

      write (text, "(i0)") n
   end function int_text

   subroutine pack_amplitudes(t1, t2, no, nv, flat)
      !! t1 then t2, one contiguous vector for DIIS
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: flat(:)

      integer :: n1

      n1 = no*nv
      flat(1:n1) = reshape(t1, [n1])
      flat(n1 + 1:) = reshape(t2, [no*no*nv*nv])
   end subroutine pack_amplitudes

   subroutine unpack_amplitudes(flat, no, nv, t1, t2)
      !! The inverse of `pack_amplitudes`
      real(dp), intent(in) :: flat(:)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: t1(:, :), t2(:, :, :, :)

      integer :: n1

      n1 = no*nv
      t1 = reshape(flat(1:n1), [no, nv])
      t2 = reshape(flat(n1 + 1:), [no, no, nv, nv])
   end subroutine unpack_amplitudes

   subroutine pack_step(t1n, t2n, t1, t2, no, nv, err_flat)
      !! The step between two iterations, flattened
      !!
      !! DIIS extrapolates the amplitudes against the change in them: the step
      !! vanishes exactly at convergence, so it is the right thing to drive to
      !! zero. Same choice as the spin-orbital path.
      real(dp), intent(in) :: t1n(:, :), t2n(:, :, :, :), t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: err_flat(:)

      integer :: n1

      n1 = no*nv
      err_flat(1:n1) = reshape(t1n - t1, [n1])
      err_flat(n1 + 1:) = reshape(t2n - t2, [no*no*nv*nv])
   end subroutine pack_step

   !===========================================================================
   ! Perturbative triples
   !
   ! JCP 94, 442 (1991), as PySCF's cc/ccsd_t_slow.py writes it. The paper's
   ! Eq. (1) has a known error in its restriction, which that file notes and
   ! corrects to [ia] >= [jb] >= [kc]; the loop below follows the corrected
   ! form rather than the printed one.
   !
   ! Three index identities let the reordered blocks that file builds be read
   ! straight out of the ones already held, so no fourth transposed copy of
   ! `ovvv` -- the largest block -- is ever allocated:
   !
   !     vvov(a,b,i,f) = ovvv(i,a,f,b)     (ia|fb)
   !     vooo(a,i,j,m) = ooov(j,m,i,a)     (ia|jm)
   !     vvoo(a,b,i,j) = ovov(i,a,j,b)     (ia|jb)
   !
   ! Memory is `n_occ^3` per intermediate and twelve of them are live at once,
   ! which is kilobytes -- the triples are a time problem, not a space one, and
   ! deliberately hold no `t3`.
   !===========================================================================

   pure subroutine permuted_triple(p, i, j, k, pi, pj, pk)
      !! The six orderings of (i,j,k), in the same order as those of (a,b,c)
      !!
      !! Pairing the two is the whole trick of the energy expression: the
      !! permutation applied to the virtual labels is mirrored in the occupied
      !! ones, so one index into a table of six covers both.
      integer, intent(in) :: p, i, j, k
      integer, intent(out) :: pi, pj, pk

      !> Which of (i,j,k) lands in each slot, for the six orderings in the
      !> order they are named in the reference: ijk, ikj, jik, jki, kij, kji.
      integer, parameter :: SLOT(3, N_TRIPLE_PERMS) = reshape([ &
                                                              1, 2, 3, &
                                                              1, 3, 2, &
                                                              2, 1, 3, &
                                                              2, 3, 1, &
                                                              3, 1, 2, &
                                                              3, 2, 1], [3, N_TRIPLE_PERMS])

      integer :: t(3)

      t = [i, j, k]
      pi = t(SLOT(1, p))
      pj = t(SLOT(2, p))
      pk = t(SLOT(3, p))
   end subroutine permuted_triple

   subroutine triples_term2(eris, t2, no, u, v, w3, w)
      !! w(i,j,k) -= sum_m (iu|jm) t2(m,k,v,w3), the ring half of W
      !!
      !! Left as a loop nest. It is O(n_occ^4 n_vir^3) against the other half's
      !! O(n_occ^3 n_vir^4), so at any size where the triples are worth running
      !! it is a twenty-sixth of the work, and the gemm that would replace it is
      !! n_occ^2 by n_occ by n_occ -- too small to pay for the call.
      type(rcc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: t2(:, :, :, :)
      integer, intent(in) :: no, u, v, w3
      real(dp), intent(inout) :: w(:, :, :)

      integer :: i, j, k, m
      real(dp) :: acc

      do k = 1, no
         do j = 1, no
            do i = 1, no
               acc = 0.0_dp
               do m = 1, no
                  acc = acc + eris%ooov(j, m, i, u)*t2(m, k, v, w3)
               end do
               w(i, j, k) = w(i, j, k) - acc
            end do
         end do
      end do
   end subroutine triples_term2

   subroutine triples_v(eris, t1, no, a, b, c, v)
      !! v(i,j,k) = (ia|jb) t1(k,c)
      !!
      !! The reference's second term, `t2(i,j,a,b) f_vo(c,k)`, is absent: the
      !! occupied-virtual Fock block vanishes for canonical orbitals, which this
      !! module assumes throughout.
      type(rcc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: t1(:, :)
      integer, intent(in) :: no, a, b, c
      real(dp), intent(out) :: v(:, :, :)

      integer :: i, j, k

      do k = 1, no
         do j = 1, no
            do i = 1, no
               v(i, j, k) = eris%ovov(i, a, j, b)*t1(k, c)
            end do
         end do
      end do
   end subroutine triples_v

   subroutine triples_r3(x, d3, no, z)
      !! z = [4x(i,j,k) + x(k,i,j) + x(j,k,i) - 2x(k,j,i) - 2x(i,k,j) - 2x(j,i,k)] / d3
      real(dp), intent(in) :: x(:, :, :), d3(:, :, :)
      integer, intent(in) :: no
      real(dp), intent(out) :: z(:, :, :)

      integer :: i, j, k

      do k = 1, no
         do j = 1, no
            do i = 1, no
               z(i, j, k) = (4.0_dp*x(i, j, k) + x(k, i, j) + x(j, k, i) &
                             - 2.0_dp*x(k, j, i) - 2.0_dp*x(i, k, j) &
                             - 2.0_dp*x(j, i, k))/d3(i, j, k)
            end do
         end do
      end do
   end subroutine triples_r3

   subroutine triples_pack_t2(t2, no, nv, tt)
      !! tt(f, (c,j,k)) = t2(k,j,c,f), built once for the whole (T)
      !!
      !! The particle half of W contracts over a virtual `f` against t2, and
      !! wants it with `f` leading. n_vir by n_occ^2 n_vir is a few tens of
      !! thousands of doubles, so it is built whole rather than per triple --
      !! which is what lets the batched products below take contiguous slices of
      !! it instead of gathering.
      real(dp), intent(in) :: t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: tt(:, :)

      integer :: c, f, j, k

      do c = 1, nv
         do f = 1, nv
            do j = 1, no
               do k = 1, no
                  tt(f, (c - 1)*no*no + (j - 1)*no + k) = t2(k, j, c, f)
               end do
            end do
         end do
      end do
   end subroutine triples_pack_t2

   subroutine triples_pack_ovvv(eris, no, nv, c0, nc, fixed, fixed_first, mmt)
      !! mmt(f, (c,i)) = (iu|fv) for the four permutations whose pair carries c
      !!
      !! `fixed_first` says which side of the pair the loop's fixed virtual sits
      !! on: true gives (i,fixed|f,c), false gives (i,c|f,fixed). Transposed --
      !! `f` leading rather than `i` -- so that the column range for this batch
      !! is contiguous and the product can take it without a copy.
      type(rcc_eris_t), intent(in) :: eris
      integer, intent(in) :: no, nv, c0, nc, fixed
      logical, intent(in) :: fixed_first
      real(dp), intent(out) :: mmt(:, :)

      integer :: cl, c, i, f

      do cl = 1, nc
         c = c0 + cl - 1
         do f = 1, nv
            do i = 1, no
               if (fixed_first) then
                  mmt(f, (cl - 1)*no + i) = eris%ovvv(i, fixed, f, c)
               else
                  mmt(f, (cl - 1)*no + i) = eris%ovvv(i, c, f, fixed)
               end if
            end do
         end do
      end do
   end subroutine triples_pack_ovvv

   subroutine triples_correction(eris, eps_o, eps_v, t1, t2, no, nv, e_triples)
      !! The (T) correction, batched over the innermost virtual
      !!
      !! The outer loop runs a >= b >= c and weights the denominator to recover
      !! the terms it skipped -- six-fold when all three coincide, two-fold when
      !! exactly two do.
      !!
      !! **Why the c loop is not innermost any more.** The particle half of W is
      !! O(n_occ^3 n_vir^4) and dominates everything else in the triples, but per
      !! (a,b,c) it is only n_occ n_vir by n_occ^2 -- a few thousand flops, far
      !! too small a matrix product to reach any useful fraction of peak, and
      !! called n_vir^3 times so that the call overhead alone is measurable.
      !!
      !! Batching the whole c range of a given (a,b) into one product fixes the
      !! shape without changing the arithmetic. The six permutations split two
      !! ways by where c appears, which is the whole of the bookkeeping below:
      !!
      !!   p = 1,3      c is the third index, so the integral half is fixed and
      !!                the amplitude half stacks   -> (n_occ x n_vir) (n_vir x n_occ^2 nc)
      !!   p = 2,4,5,6  c is in the pair, so the amplitude half is fixed and the
      !!                integral half stacks        -> (nc n_occ x n_vir) (n_vir x n_occ^2)
      !!
      !! Everything downstream -- the denominators, r3, and the thirty-six-term
      !! sum -- is unchanged and still runs one (a,b,c) at a time.
      type(rcc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: eps_o(:), eps_v(:)
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: e_triples

      !> Which W pairs with which occupied permutation, for each Z. Row is the
      !> Z's virtual permutation, column the occupied one, entry the W's. This
      !> is the S3 multiplication table and is written out rather than derived:
      !> deriving it needs a convention for whether a permutation acts on
      !> positions or on labels, and getting that backwards transposes half the
      !> table into a still-symmetric, still-plausible, wrong answer.
      integer, parameter :: WMAP(N_TRIPLE_PERMS, N_TRIPLE_PERMS) = reshape([ &
                                                                           1, 2, 3, 4, 5, 6, &
                                                                           2, 1, 5, 6, 3, 4, &
                                                                           3, 4, 1, 2, 6, 5, &
                                                                           4, 3, 6, 5, 1, 2, &
                                                                           5, 6, 2, 1, 4, 3, &
                                                                           6, 5, 4, 3, 2, 1], &
                                                                           [N_TRIPLE_PERMS, N_TRIPLE_PERMS], order=[2, 1])

      real(dp), allocatable :: w(:, :, :, :, :), z(:, :, :, :), x(:, :, :), d3(:, :, :)
      real(dp), allocatable :: tt(:, :), mfix(:, :), mmt(:, :), r1(:, :), r2(:, :)
      integer, allocatable :: pidx(:, :), alist(:), blist(:)
      real(dp), allocatable :: e_pair(:)
      integer :: a, b, c, cl, c0, c1, i, j, k, p, q, pi, pj, pk, nc, no2
      integer :: idx, npair
      integer :: trip(3), perm(3, N_TRIPLE_PERMS)
      real(dp) :: scale

      no2 = no*no

      ! Only the shared arrays are allocated here. W and the scratch that feeds
      ! it are allocated per thread inside the parallel region below, which is
      ! also where they are sized -- an allocatable named in a `private` clause
      ! arrives unallocated, and allocating it out here would be both wrong and
      ! a per-thread copy of the wrong size.
      allocate (tt(nv, no2*nv))

      ! Where each occupied permutation sends each (i,j,k), once.
      !
      ! The thirty-six-term sum below touches n_occ^3 elements for each of
      ! thirty-six pairings of every triple, so any per-element work there is
      ! multiplied by tens of millions. Resolving the permutation by a table
      ! lookup rather than by a call was worth 35% of the whole correction.
      allocate (pidx(no*no*no, N_TRIPLE_PERMS))
      do q = 1, N_TRIPLE_PERMS
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  call permuted_triple(q, i, j, k, pi, pj, pk)
                  pidx(i + (j - 1)*no + (k - 1)*no2, q) = pi + (pj - 1)*no + (pk - 1)*no2
               end do
            end do
         end do
      end do

      call triples_pack_t2(t2, no, nv, tt)

      ! The (a,b) loop is triangular, so it is flattened before being handed to
      ! OpenMP. Two reasons, and the second is the one that matters: a triangular
      ! nest cannot be collapsed, and the work per outer index runs from nothing
      ! at a = 1 to the whole inner range at a = n_vir, so splitting on `a` alone
      ! would leave the first threads idle. One flat index over the pairs, taken
      ! dynamically, balances itself.
      npair = nv*(nv + 1)/2
      allocate (alist(npair), blist(npair), e_pair(npair))
      idx = 0
      do a = 1, nv
         do b = 1, a
            idx = idx + 1
            alist(idx) = a
            blist(idx) = b
         end do
      end do

      e_pair = 0.0_dp

      ! Everything above this line is shared and read-only for the duration:
      ! the packed amplitudes, the permutation table and the pair list. Each
      ! thread keeps its own W and its own scratch, which is what TRIPLES_C_BATCH
      ! is sized against.
      !$omp parallel default(none) &
      !$omp    shared(eris, t1, t2, eps_o, eps_v, no, nv, no2, tt, pidx, &
      !$omp           alist, blist, npair) &
      !$omp    shared(e_pair) &
      !$omp    private(w, z, x, d3, mfix, mmt, r1, r2, &
      !$omp            a, b, c, cl, c0, c1, nc, i, j, k, p, q, trip, perm, scale)
      allocate (w(no, no, no, N_TRIPLE_PERMS, TRIPLES_C_BATCH))
      allocate (z(no, no, no, N_TRIPLE_PERMS))
      allocate (x(no, no, no), d3(no, no, no))
      allocate (mfix(no, nv), mmt(nv, TRIPLES_C_BATCH*no))
      allocate (r1(no, no2*TRIPLES_C_BATCH), r2(TRIPLES_C_BATCH*no, no2))

      !$omp do schedule(dynamic)
      do idx = 1, npair
         a = alist(idx)
         b = blist(idx)

         do c0 = 1, b, TRIPLES_C_BATCH
            c1 = min(c0 + TRIPLES_C_BATCH - 1, b)
            nc = c1 - c0 + 1

            ! ---- p = 1 and 3: the pair is fixed, the amplitudes stack -------
            do concurrent(i=1:no, k=1:nv)
               mfix(i, k) = eris%ovvv(i, a, k, b)
            end do
            call pic_gemm(mfix, tt(:, (c0 - 1)*no2 + 1:c1*no2), r1(:, 1:no2*nc))
            call scatter_w_third(r1, no, nc, 1, w)

            do concurrent(i=1:no, k=1:nv)
               mfix(i, k) = eris%ovvv(i, b, k, a)
            end do
            call pic_gemm(mfix, tt(:, (c0 - 1)*no2 + 1:c1*no2), r1(:, 1:no2*nc))
            call scatter_w_third(r1, no, nc, 3, w)

            ! ---- p = 2, 4, 5, 6: the amplitudes are fixed, the pair stacks --
            ! The third index is b for p = 2 and 5, a for p = 4 and 6, so each
            ! takes the matching contiguous block of `tt`.
            call triples_pack_ovvv(eris, no, nv, c0, nc, a, .true., mmt)
            call pic_gemm(mmt(:, 1:nc*no), tt(:, (b - 1)*no2 + 1:b*no2), &
                          r2(1:nc*no, :), transa="T")
            call scatter_w_pair(r2, no, nc, 2, w)

            call triples_pack_ovvv(eris, no, nv, c0, nc, b, .true., mmt)
            call pic_gemm(mmt(:, 1:nc*no), tt(:, (a - 1)*no2 + 1:a*no2), &
                          r2(1:nc*no, :), transa="T")
            call scatter_w_pair(r2, no, nc, 4, w)

            call triples_pack_ovvv(eris, no, nv, c0, nc, a, .false., mmt)
            call pic_gemm(mmt(:, 1:nc*no), tt(:, (b - 1)*no2 + 1:b*no2), &
                          r2(1:nc*no, :), transa="T")
            call scatter_w_pair(r2, no, nc, 5, w)

            call triples_pack_ovvv(eris, no, nv, c0, nc, b, .false., mmt)
            call pic_gemm(mmt(:, 1:nc*no), tt(:, (a - 1)*no2 + 1:a*no2), &
                          r2(1:nc*no, :), transa="T")
            call scatter_w_pair(r2, no, nc, 6, w)

            do cl = 1, nc
               c = c0 + cl - 1
               perm(:, 1) = [a, b, c]
               perm(:, 2) = [a, c, b]
               perm(:, 3) = [b, a, c]
               perm(:, 4) = [b, c, a]
               perm(:, 5) = [c, a, b]
               perm(:, 6) = [c, b, a]

               ! The ring half, still one permutation at a time.
               do p = 1, N_TRIPLE_PERMS
                  trip = perm(:, p)
                  call triples_term2(eris, t2, no, trip(1), trip(2), trip(3), &
                                     w(:, :, :, p, cl))
               end do

               scale = 1.0_dp
               if (a == c) then
                  scale = 6.0_dp          ! a == b == c
               else if (a == b .or. b == c) then
                  scale = 2.0_dp
               end if

               do k = 1, no
                  do j = 1, no
                     do i = 1, no
                        d3(i, j, k) = scale*(eps_o(i) + eps_o(j) + eps_o(k) &
                                             - eps_v(a) - eps_v(b) - eps_v(c))
                     end do
                  end do
               end do

               do p = 1, N_TRIPLE_PERMS
                  trip = perm(:, p)
                  call triples_v(eris, t1, no, trip(1), trip(2), trip(3), x)
                  x = w(:, :, :, p, cl) + 0.5_dp*x
                  call triples_r3(x, d3, no, z(:, :, :, p))
               end do

               ! The thirty-six terms: every Z against every occupied
               ! permutation of the W the table pairs it with.
               do p = 1, N_TRIPLE_PERMS
                  do q = 1, N_TRIPLE_PERMS
                     e_pair(idx) = e_pair(idx) &
                                   + triples_dot(no*no*no, pidx(1, q), &
                                                 w(1, 1, 1, WMAP(p, q), cl), z(1, 1, 1, p))
                  end do
               end do
            end do
         end do
      end do
      !$omp end do

      deallocate (w, z, x, d3, mfix, mmt, r1, r2)
      !$omp end parallel

      ! Summed here, in pair order, rather than by an OpenMP reduction. Each
      ! pair's contribution is written by whichever thread took it and read by
      ! nobody else, so there is no race either way -- but a reduction combines
      ! the per-thread partials in whatever order the threads happen to finish,
      ! and dynamic scheduling means that order changes between runs. Summing a
      ! fixed array in a fixed order makes the correction independent of the
      ! thread count and of the scheduling, which is worth a n_vir^2/2 array:
      ! this energy is differenced against others in a many-body expansion, and
      ! there a wandering last digit is not noise but a term.
      e_triples = 2.0_dp*sum(e_pair)

      deallocate (tt, pidx, alist, blist, e_pair)
   end subroutine triples_correction

   pure function triples_dot(no3, permuted, wblk, zblk) result(acc)
      !! sum over one occupied permutation of W against Z
      !!
      !! Both blocks are n_occ^3 contiguous doubles -- one permutation's W and
      !! one Z, taken by sequence association out of the arrays that hold all
      !! six of each. `permuted` is where the permutation sends each flat index,
      !! so the permutation costs a load rather than a call.
      integer, intent(in) :: no3
      integer, intent(in) :: permuted(no3)
      real(dp), intent(in) :: wblk(no3), zblk(no3)
      real(dp) :: acc

      integer :: n

      acc = 0.0_dp
      do n = 1, no3
         acc = acc + wblk(permuted(n))*zblk(n)
      end do
   end function triples_dot

   subroutine scatter_w_third(r, no, nc, p, w)
      !! R(i, (c,j,k)) -> w(i,j,k,p,c), for the permutations where c is third
      real(dp), intent(in) :: r(:, :)
      integer, intent(in) :: no, nc, p
      real(dp), intent(inout) :: w(:, :, :, :, :)

      integer :: c, i, j, k

      do c = 1, nc
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  w(i, j, k, p, c) = r(i, (c - 1)*no*no + (j - 1)*no + k)
               end do
            end do
         end do
      end do
   end subroutine scatter_w_third

   subroutine scatter_w_pair(r, no, nc, p, w)
      !! R((c,i), (j,k)) -> w(i,j,k,p,c), for the permutations where c is paired
      real(dp), intent(in) :: r(:, :)
      integer, intent(in) :: no, nc, p
      real(dp), intent(inout) :: w(:, :, :, :, :)

      integer :: c, i, j, k

      do c = 1, nc
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  w(i, j, k, p, c) = r((c - 1)*no + i, (j - 1)*no + k)
               end do
            end do
         end do
      end do
   end subroutine scatter_w_pair

end module mqc_libcint_rcc
