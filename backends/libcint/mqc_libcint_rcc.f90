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

   !> Columns of the compound (cd) index assembled per pass of the ladder.
   !> Same role and same reasoning as LADDER_BATCH in the spin-orbital module:
   !> the array held is n_vir^2 by this, so the memory stays O(n_vir^2).
   integer, parameter :: LADDER_BATCH = 256

   character(len=*), parameter :: STAGE_ITER = "RCCSD iterations"

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

   subroutine build_rcc_eris_conventional(mol, c_occ, c_vir, eris)
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

      real(dp), allocatable :: eri(:, :)

      call mol%eris_packed(eri)
      call transform_block(eri, c_occ, c_occ, c_occ, c_occ, eris%oooo)
      call transform_block(eri, c_occ, c_occ, c_occ, c_vir, eris%ooov)
      call transform_block(eri, c_occ, c_occ, c_vir, c_vir, eris%oovv)
      call transform_block(eri, c_occ, c_vir, c_occ, c_vir, eris%ovov)
      call transform_block(eri, c_occ, c_vir, c_vir, c_vir, eris%ovvv)
      call transform_block(eri, c_vir, c_vir, c_vir, c_vir, eris%vvvv)
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

      call build_df_mo_block(mol, aux, c_occ, c_occ, b_oo, error)
      if (error%has_error()) return
      call build_df_mo_block(mol, aux, c_occ, c_vir, b_ov, error)
      if (error%has_error()) return
      call build_df_mo_block(mol, aux, c_vir, c_vir, b_vv_pq, error)
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

      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  tau(i, j, a, b) = t2(i, j, a, b) + t1(i, a)*t1(j, b)
               end do
            end do
         end do
      end do
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
      do d = 1, nv
         do l = 1, no
            do k = 1, no
               do c = 1, nv
                  do a = 1, nv
                     fac(a, c) = fac(a, c) - lvec(eris%ovov, k, c, l, d)*tau(k, l, a, d)
                  end do
               end do
            end do
         end do
      end do
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

      do c = 1, nv
         do j = 1, no
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

      do d = 1, nv
         do c = 1, nv
            do j = 1, no
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

      do d = 1, nv
         do c = 1, nv
            do i = 1, no
               do k = 1, no
                  do a = 1, nv
                     w(a, k, i, c) = w(a, k, i, c) + eris%ovvv(k, c, a, d)*t1(i, d)
                  end do
               end do
            end do
         end do
      end do

      do l = 1, no
         do c = 1, nv
            do i = 1, no
               do k = 1, no
                  do a = 1, nv
                     w(a, k, i, c) = w(a, k, i, c) - eris%ooov(l, i, k, c)*t1(l, a)
                  end do
               end do
            end do
         end do
      end do

      do d = 1, nv
         do l = 1, no
            do c = 1, nv
               do i = 1, no
                  do k = 1, no
                     do a = 1, nv
                        w(a, k, i, c) = w(a, k, i, c) &
                                        - 0.5_dp*eris%ovov(l, d, k, c)*t2(i, l, d, a) &
                                        - 0.5_dp*eris%ovov(l, c, k, d)*t2(i, l, a, d) &
                                        - eris%ovov(l, d, k, c)*t1(i, d)*t1(l, a) &
                                        + eris%ovov(l, d, k, c)*t2(i, l, a, d)
                     end do
                  end do
               end do
            end do
         end do
      end do
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

      do d = 1, nv
         do i = 1, no
            do c = 1, nv
               do k = 1, no
                  do a = 1, nv
                     w(a, k, c, i) = w(a, k, c, i) + eris%ovvv(k, d, a, c)*t1(i, d)
                  end do
               end do
            end do
         end do
      end do

      do l = 1, no
         do i = 1, no
            do c = 1, nv
               do k = 1, no
                  do a = 1, nv
                     w(a, k, c, i) = w(a, k, c, i) - eris%ooov(k, i, l, c)*t1(l, a)
                  end do
               end do
            end do
         end do
      end do

      do d = 1, nv
         do l = 1, no
            do i = 1, no
               do c = 1, nv
                  do k = 1, no
                     do a = 1, nv
                        w(a, k, c, i) = w(a, k, c, i) &
                                        - 0.5_dp*eris%ovov(l, c, k, d)*t2(i, l, d, a) &
                                        - eris%ovov(l, c, k, d)*t1(i, d)*t1(l, a)
                     end do
                  end do
               end do
            end do
         end do
      end do
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
      do d = 1, nv
         do c = 1, nv
            do k = 1, no
               do a = 1, nv
                  do i = 1, no
                     t1n(i, a) = t1n(i, a) &
                                 + 2.0_dp*eris%ovvv(k, d, a, c)*(t2(i, k, c, d) + t1(k, d)*t1(i, c)) &
                                 - eris%ovvv(k, c, a, d)*(t2(i, k, c, d) + t1(k, d)*t1(i, c))
                  end do
               end do
            end do
         end do
      end do

      ! '-2 lcki,klac->ia' + 'kcli,klac->ia' with (lc|ki) = ooov(k,i,l,c),
      ! (kc|li) = ooov(l,i,k,c); and the same pair with t2 -> t1 t1
      do c = 1, nv
         do l = 1, no
            do k = 1, no
               do a = 1, nv
                  do i = 1, no
                     t1n(i, a) = t1n(i, a) &
                                 - 2.0_dp*eris%ooov(k, i, l, c)*(t2(k, l, a, c) + t1(l, c)*t1(k, a)) &
                                 + eris%ooov(l, i, k, c)*(t2(k, l, a, c) + t1(l, c)*t1(k, a))
                  end do
               end do
            end do
         end do
      end do

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
      tmp = 0.0_dp
      do c = 1, nv
         do b = 1, nv
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     tmp(i, j, a, b) = tmp(i, j, a, b) + tmp2(a, b, i, c)*t1(j, c)
                  end do
               end do
            end do
         end do
      end do
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
      do c = 1, nv
         do b = 1, nv
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     tmp(i, j, a, b) = tmp(i, j, a, b) + lac(a, c)*t2(i, j, c, b)
                  end do
               end do
            end do
         end do
      end do
      call add_symmetrised(t2n, tmp, no, nv, 1.0_dp)

      tmp = 0.0_dp
      do k = 1, no
         do b = 1, nv
            do a = 1, nv
               do j = 1, no
                  do i = 1, no
                     tmp(i, j, a, b) = tmp(i, j, a, b) + lki(k, i)*t2(k, j, a, b)
                  end do
               end do
            end do
         end do
      end do
      call add_symmetrised(t2n, tmp, no, nv, -1.0_dp)

      ! The three ring terms
      allocate (wvoov(nv, no, no, nv), wvovo(nv, no, nv, no))
      call cc_wvoov(eris, t1, t2, no, nv, wvoov)
      call cc_wvovo(eris, t1, t2, no, nv, wvovo)

      ! '2 akic,kjcb->ijab' - 'akci,kjcb->ijab'
      tmp = 0.0_dp
      do c = 1, nv
         do k = 1, no
            do b = 1, nv
               do a = 1, nv
                  do j = 1, no
                     do i = 1, no
                        tmp(i, j, a, b) = tmp(i, j, a, b) &
                                          + (2.0_dp*wvoov(a, k, i, c) - wvovo(a, k, c, i)) &
                                          *t2(k, j, c, b)
                     end do
                  end do
               end do
            end do
         end do
      end do
      call add_symmetrised(t2n, tmp, no, nv, 1.0_dp)

      ! '-akic,kjbc->ijab'
      tmp = 0.0_dp
      do c = 1, nv
         do k = 1, no
            do b = 1, nv
               do a = 1, nv
                  do j = 1, no
                     do i = 1, no
                        tmp(i, j, a, b) = tmp(i, j, a, b) + wvoov(a, k, i, c)*t2(k, j, b, c)
                     end do
                  end do
               end do
            end do
         end do
      end do
      call add_symmetrised(t2n, tmp, no, nv, -1.0_dp)

      ! '-bkci,kjac->ijab'
      tmp = 0.0_dp
      do c = 1, nv
         do k = 1, no
            do b = 1, nv
               do a = 1, nv
                  do j = 1, no
                     do i = 1, no
                        tmp(i, j, a, b) = tmp(i, j, a, b) + wvovo(b, k, c, i)*t2(k, j, a, c)
                     end do
                  end do
               end do
            end do
         end do
      end do
      call add_symmetrised(t2n, tmp, no, nv, -1.0_dp)
      deallocate (wvoov, wvovo)

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
   end subroutine rccsd_iteration

   subroutine particle_ladder(eris, t1, tau, no, nv, t2n)
      !! t2new(i,j,a,b) += sum_cd Wvvvv(a,b,c,d) tau(i,j,c,d), never holding Wvvvv
      !!
      !!     Wvvvv(a,b,c,d) = (ac|bd) - sum_k ovvv(k,d,a,c) t1(k,b)
      !!                              - sum_k ovvv(k,c,b,d) t1(k,a)
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
      real(dp), intent(in) :: t1(:, :), tau(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(inout) :: t2n(:, :, :, :)

      real(dp), allocatable :: wblk(:, :, :), gvv(:, :)
      integer :: nv2, cd0, cd1, nb, col, cd, a, b, c, d, i, j, k, naux
      logical :: fitted
      real(dp) :: acc

      nv2 = nv*nv
      fitted = allocated(eris%b_vv)
      if (fitted) naux = size(eris%b_vv, 2)

      allocate (wblk(nv, nv, LADDER_BATCH))
      if (fitted) allocate (gvv(nv, nv))

      do cd0 = 1, nv2, LADDER_BATCH
         cd1 = min(cd0 + LADDER_BATCH - 1, nv2)
         nb = cd1 - cd0 + 1

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

            ! The t1 dressing. n_occ n_vir^4 an iteration, the same order as the
            ! contraction it feeds.
            do b = 1, nv
               do a = 1, nv
                  acc = 0.0_dp
                  do k = 1, no
                     acc = acc - eris%ovvv(k, d, a, c)*t1(k, b) &
                           - eris%ovvv(k, c, b, d)*t1(k, a)
                  end do
                  wblk(a, b, col) = wblk(a, b, col) + acc
               end do
            end do
         end do

         do col = 1, nb
            cd = cd0 + col - 1
            c = mod(cd - 1, nv) + 1
            d = (cd - 1)/nv + 1
            do b = 1, nv
               do a = 1, nv
                  acc = wblk(a, b, col)
                  if (acc == 0.0_dp) cycle
                  do j = 1, no
                     do i = 1, no
                        t2n(i, j, a, b) = t2n(i, j, a, b) + acc*tau(i, j, c, d)
                     end do
                  end do
               end do
            end do
         end do
      end do

      deallocate (wblk)
      if (allocated(gvv)) deallocate (gvv)
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

      ! (T) is not on this path yet. Refused rather than silently returning zero:
      ! a CCSD(T) deck that got CCSD and a zero triples term would report a
      ! number that is wrong by exactly the thing that was asked for.
      if (want_triples) then
         call error%set(ERROR_VALIDATION, "RCCSD: the perturbative triples are not "// &
                        "implemented in the spatial-orbital path yet. Run CCSD, or take "// &
                        "CCSD(T) through the spin-orbital path.")
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
         call build_rcc_eris_fitted(mol, aux, c_act(:, 1:no), c_act(:, no + 1:n_act), &
                                    eris, error)
         if (error%has_error()) return
      else
         call build_rcc_eris_conventional(mol, c_act(:, 1:no), c_act(:, no + 1:n_act), eris)
      end if
      deallocate (c_act)
      call clk%lap("AO->MO integrals")

      allocate (eps_o(no), eps_v(nv))
      do i = 1, no
         eps_o(i) = orbital_energies(frozen + i)
      end do
      do a = 1, nv
         eps_v(a) = orbital_energies(n_occ + a)
      end do

      ! ---- MP2, the checkpoint before any amplitude equation ----------------
      allocate (t1(no, nv), t2(no, no, nv, nv))
      t1 = 0.0_dp
      call mp2_guess(eris, eps_o, eps_v, no, nv, t2, result%e_mp2)
      call clk%lap("MP2 amplitudes")
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

      result%e_correlation = result%e_singles + result%e_doubles + result%e_triples

      if (.not. result%converged) then
         call error%set(ERROR_VALIDATION, "RCCSD did not converge")
         return
      end if
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

end module mqc_libcint_rcc
