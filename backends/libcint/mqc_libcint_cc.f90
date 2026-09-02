!! Coupled cluster on the CPU, in the spin-orbital basis
module mqc_libcint_cc
   !! CCSD and CCSD(T) over a restricted or an unrestricted reference,
   !! written in spin orbitals.
   !!
   !! Every index is a spin orbital, and which spatial orbital each one belongs
   !! to is the table `so_map_t`. Equal occupied and virtual counts interleave --
   !! alpha at `2p-1`, beta at `2p` -- and unequal counts block as
   !! `[occ_a, occ_b, vir_a, vir_b]`, because interleaving keeps the occupied
   !! space contiguous only while the counts agree. The branch is on the counts
   !! and not on the reference, so a closed-shell unrestricted deck interleaves
   !! too. See `build_so_map`.
   !!
   !! The antisymmetrised integrals are
   !!
   !!     <pq||rs> = (pr|qs) - (ps|qr)
   !!
   !! with each spatial integral surviving only if the spins on the two indices
   !! it pairs agree. The transform produces chemists' notation and the equations
   !! are written in physicists'; that conversion happens in `asym` alone.
   !!
   !! The amplitude equations follow Stanton, Gauss, Watts and Bartlett (JCP 94,
   !! 4334 (1991)) with their one- and two-body intermediates.
   !!
   !! **Frozen core** drops orbitals from the active space before anything else
   !! happens, so nothing below this line knows they existed.
   !!
   !! Every antisymmetrised block is materialised in spin orbitals, sixteen times
   !! its spatial size, so `<ma||ef>` and `<ab||ej>` at n_occ n_vir^3 are the
   !! wall. `<ab||cd>` is not materialised -- `particle_ladder` assembles it in
   !! batches -- and on the fitted path neither is the n_act^4 spatial tensor:
   !! RI-CCSD holds five spatial blocks with an occupied index and the
   !! three-index `b_vv`.
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

   public :: cc_result_t
   public :: run_libcint_ccsd
   public :: spin_orbital_integrals   !! Exported so a test can check antisymmetry
   public :: so_map_t                 !! ... which needs the ordering it was built with
   public :: so_source_t              !! ... and the tensors it reads
   public :: build_so_map

   integer, parameter :: LADDER_BATCH = 256
      !! Columns of the compound (ef) index assembled per pass of the
      !! particle-particle ladder. The block held is n_vir^2 by this, so it has
      !! to stay a fixed count rather than a fraction of n_vir^2.

   character(len=*), parameter :: STAGE_ITER = "CCSD iterations"
      !! Timing stage the amplitude update accumulates under; the per-iteration
      !! row reads it back with `seconds_of`.

   type :: so_map_t
      !! Which spatial orbital, and which spin, each spin orbital is
      !!
      !! Occupied before virtual. Equal alpha and beta counts interleave, spin
      !! orbitals `2p-1` and `2p` sharing spatial orbital `p`; unequal counts
      !! block, because interleaving would put beta virtuals among the alpha
      !! occupieds and every `do i = 1, no` below would walk into the virtual
      !! space -- converging, to a wrong number.
      integer, allocatable :: spat(:)   !! spin orbital -> spatial index, within its spin
      logical, allocatable :: is_a(:)   !! spin orbital -> alpha?
      integer :: n_so = 0
   end type so_map_t

   type :: so_source_t
      !! The spatial MO tensors the spin-orbital integrals are read out of
      !!
      !! One tensor for a restricted reference, three for an unrestricted one:
      !! `(aa|aa)`, `(bb|bb)` and the mixed `(aa|bb)`. There is no fourth --
      !! `(bb|aa)` is `(aa|bb)` with the electron pairs swapped, which is what
      !! `chem` does rather than storing it twice. When `restricted`, only `aa`
      !! is allocated and every spin case reads it.
      logical :: restricted = .true.
      real(dp), allocatable :: aa(:, :, :, :)
      real(dp), allocatable :: bb(:, :, :, :)
      real(dp), allocatable :: ab(:, :, :, :)
   end type so_source_t

   type :: cc_eris_t
      !! Antisymmetrised integrals, by block
      !!
      !! No `vvvv`: at n_vir^4 it is the largest block, and it is read in exactly
      !! one place -- `particle_ladder`, which assembles it in batches.
      real(dp), allocatable :: oooo(:, :, :, :)   !! <ij||kl>
      real(dp), allocatable :: ooov(:, :, :, :)   !! <mn||ie>
      real(dp), allocatable :: oovv(:, :, :, :)   !! <ij||ab>
      real(dp), allocatable :: ovov(:, :, :, :)   !! <ia||jb>
      real(dp), allocatable :: ovvo(:, :, :, :)   !! <mb||ej>
      real(dp), allocatable :: ovoo(:, :, :, :)   !! <mb||ij>
      real(dp), allocatable :: ovvv(:, :, :, :)   !! <ma||ef>
      real(dp), allocatable :: vvvo(:, :, :, :)   !! <ab||ej>
   end type cc_eris_t

   type :: cc_chem_t
      !! Fitted spatial integrals, in chemists' notation, by block
      !!
      !! Only the blocks with at least one occupied index. `(vv|vv)` is absent
      !! for the same reason `cc_eris_t` has no `vvvv`: `particle_ladder`
      !! assembles it, here from `b_vv` rather than from a stored tensor.
      real(dp), allocatable :: oooo(:, :, :, :)   !! (ij|kl)
      real(dp), allocatable :: ooov(:, :, :, :)   !! (ij|ka)
      real(dp), allocatable :: oovv(:, :, :, :)   !! (ij|ab)
      real(dp), allocatable :: ovov(:, :, :, :)   !! (ia|jb)
      real(dp), allocatable :: ovvv(:, :, :, :)   !! (ia|bc)
   end type cc_chem_t

   type :: cc_result_t
      !! What a converged coupled cluster calculation leaves behind
      real(dp) :: e_singles = 0.0_dp    !! The t1 t1 part of the correlation energy
      real(dp) :: e_doubles = 0.0_dp    !! The t2 part
      real(dp) :: e_triples = 0.0_dp    !! (T), zero unless it was asked for
      real(dp) :: e_mp2 = 0.0_dp
         !! MP2 from these same spin-orbital integrals, which must reproduce
         !! `run_libcint_mp2` exactly. It is the first iteration's energy, so it
         !! costs nothing, and it checks the antisymmetrisation and the
         !! denominators.
      real(dp) :: e_correlation = 0.0_dp  !! singles + doubles + triples
      integer :: iterations = 0
      logical :: converged = .false.
   end type cc_result_t

contains

   pure function integral_megabytes(no, nv) result(mb)
      !! What the integral blocks and the ladder batch cost, in MB
      !!
      !! `<ab||ef>` is never held, so its contribution is one batch, n_vir^2 by
      !! `LADDER_BATCH`; that is why this is n_occ n_vir^3 and not n_vir^4.
      integer, intent(in) :: no, nv
      real(dp) :: mb

      real(dp) :: elements

      elements = real(no, dp)**4 &                        ! oooo
                 + 2.0_dp*real(no, dp)**3*real(nv, dp) &  ! ooov, ovoo
                 + 3.0_dp*real(no, dp)**2*real(nv, dp)**2 &  ! oovv, ovov, ovvo
                 + 2.0_dp*real(no, dp)*real(nv, dp)**3 &  ! ovvv, vvvo
                 + real(nv, dp)**2*real(LADDER_BATCH, dp)   ! the ladder batch
      mb = elements*8.0_dp/1.0e6_dp
   end function integral_megabytes

   subroutine build_so_map(n_occ_a, n_occ_b, n_vir_a, n_vir_b, map)
      !! The spin orbital ordering, occupied first
      !!
      !! Equal counts (`n_occ_a == n_occ_b` and `n_vir_a == n_vir_b`) interleave:
      !! spin orbital `2p-1` is alpha spatial `p`, `2p` is beta spatial `p`.
      !! Unequal counts block as `[occ_a, occ_b, vir_a, vir_b]`, the only ordering
      !! that keeps the occupied space contiguous.
      integer, intent(in) :: n_occ_a, n_occ_b, n_vir_a, n_vir_b
      type(so_map_t), intent(out) :: map

      integer :: s, p

      map%n_so = n_occ_a + n_occ_b + n_vir_a + n_vir_b
      allocate (map%spat(map%n_so), map%is_a(map%n_so))

      if (n_occ_a == n_occ_b .and. n_vir_a == n_vir_b) then
         do s = 1, map%n_so
            map%spat(s) = (s + 1)/2
            map%is_a(s) = mod(s, 2) == 1
         end do
         return
      end if

      s = 0
      do p = 1, n_occ_a
         s = s + 1
         map%spat(s) = p
         map%is_a(s) = .true.
      end do
      do p = 1, n_occ_b
         s = s + 1
         map%spat(s) = p
         map%is_a(s) = .false.
      end do
      do p = 1, n_vir_a
         s = s + 1
         map%spat(s) = n_occ_a + p
         map%is_a(s) = .true.
      end do
      do p = 1, n_vir_b
         s = s + 1
         map%spat(s) = n_occ_b + p
         map%is_a(s) = .false.
      end do
   end subroutine build_so_map

   pure function spatial(map, s) result(p)
      !! Which spatial orbital a spin orbital belongs to
      type(so_map_t), intent(in) :: map
      integer, intent(in) :: s
      integer :: p
      p = map%spat(s)
   end function spatial

   pure function same_spin(map, s, t) result(same)
      !! Whether two spin orbitals carry the same spin
      type(so_map_t), intent(in) :: map
      integer, intent(in) :: s, t
      logical :: same
      same = map%is_a(s) .eqv. map%is_a(t)
   end function same_spin

   pure function tensor(src, s1, s2, a, b, c, d) result(v)
      !! `(ab|cd)` from spatial indices whose spins the caller has already resolved
      !!
      !! `s1` is the spin of the first electron's pair, `s2` the second's; 1 is
      !! alpha and 0 beta. `chem` makes the same choice from spin orbitals, for
      !! callers that still hold them.
      type(so_source_t), intent(in) :: src
      integer, intent(in) :: s1, s2, a, b, c, d
      real(dp) :: v

      if (src%restricted) then
         v = src%aa(a, b, c, d)
      else if (s1 == 1 .and. s2 == 1) then
         v = src%aa(a, b, c, d)
      else if (s1 == 0 .and. s2 == 0) then
         v = src%bb(a, b, c, d)
      else if (s1 == 1) then
         v = src%ab(a, b, c, d)
      else
         v = src%ab(c, d, a, b)
      end if
   end function tensor

   pure function chem(src, map, a, b, c, d) result(v)
      !! The spatial integral `(ab|cd)`, for spin orbitals whose spins already agree
      !!
      !! Callers guarantee `spin(a) == spin(b)` and `spin(c) == spin(d)`, so this
      !! only chooses which of the three tensors holds that combination. A
      !! beta-alpha pairing is stored as alpha-beta with the electron pairs
      !! swapped, which `(ab|cd) = (cd|ab)` makes legitimate.
      type(so_source_t), intent(in) :: src
      type(so_map_t), intent(in) :: map
      integer, intent(in) :: a, b, c, d
      real(dp) :: v

      if (src%restricted) then
         v = src%aa(map%spat(a), map%spat(b), map%spat(c), map%spat(d))
      else if (map%is_a(a) .and. map%is_a(c)) then
         v = src%aa(map%spat(a), map%spat(b), map%spat(c), map%spat(d))
      else if (.not. map%is_a(a) .and. .not. map%is_a(c)) then
         v = src%bb(map%spat(a), map%spat(b), map%spat(c), map%spat(d))
      else if (map%is_a(a)) then
         v = src%ab(map%spat(a), map%spat(b), map%spat(c), map%spat(d))
      else
         v = src%ab(map%spat(c), map%spat(d), map%spat(a), map%spat(b))
      end if
   end function chem

   pure function asym(src, map, p, q, r, s) result(v)
      !! <pq||rs> from the spatial MO tensors, in spin orbitals
      !!
      !! The one place chemists' notation becomes physicists':
      !!
      !!     <pq||rs> = (pr|qs) - (ps|qr)
      !!
      !! A spatial integral contributes only when both index pairs it couples
      !! share a spin -- (pr|qs) needs spin(p) == spin(r) and spin(q) == spin(s),
      !! the Kronecker delta the spin integration leaves behind. Those guards are
      !! what lets `chem` assume the spins within each pair already agree.
      type(so_source_t), intent(in) :: src
      type(so_map_t), intent(in) :: map
      integer, intent(in) :: p, q, r, s
      real(dp) :: v

      v = 0.0_dp
      if (same_spin(map, p, r) .and. same_spin(map, q, s)) then
         v = v + chem(src, map, p, r, q, s)
      end if
      if (same_spin(map, p, s) .and. same_spin(map, q, r)) then
         v = v - chem(src, map, p, s, q, r)
      end if
   end function asym

   subroutine build_cc_eris(src, map, no, nv, eris)
      !! The antisymmetrised blocks the amplitude equations actually read
      !!
      !! Every block is a slice of the single (2 n_act)^4 tensor. The slice that
      !! is missing -- `vvvv` -- is the one that dominates; `particle_ladder`
      !! assembles it in batches instead, which takes the scaling from n_vir^4 to
      !! n_occ n_vir^3.
      !!
      !! `spin_orbital_integrals` is kept for the tests, which check antisymmetry
      !! over the whole tensor. Nothing in the calculation builds it.
      type(so_source_t), intent(in) :: src
      type(so_map_t), intent(in) :: map
      integer, intent(in) :: no, nv
      type(cc_eris_t), intent(out) :: eris

      integer :: i, j, k, l, a, b, e, f, m, n

      allocate (eris%oooo(no, no, no, no), eris%ooov(no, no, no, nv))
      allocate (eris%oovv(no, no, nv, nv), eris%ovov(no, nv, no, nv))
      allocate (eris%ovvo(no, nv, nv, no), eris%ovoo(no, nv, no, no))
      allocate (eris%ovvv(no, nv, nv, nv), eris%vvvo(nv, nv, nv, no))

      !$omp parallel do default(none) shared(src, map, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do l = 1, no
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  eris%oooo(i, j, k, l) = asym(src, map, i, j, k, l)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(src, map, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do e = 1, nv
         do i = 1, no
            do n = 1, no
               do m = 1, no
                  eris%ooov(m, n, i, e) = asym(src, map, m, n, i, no + e)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(src, map, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  eris%oovv(i, j, a, b) = asym(src, map, i, j, no + a, no + b)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(src, map, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do b = 1, nv
         do j = 1, no
            do a = 1, nv
               do i = 1, no
                  eris%ovov(i, a, j, b) = asym(src, map, i, no + a, j, no + b)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(src, map, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do j = 1, no
         do e = 1, nv
            do b = 1, nv
               do m = 1, no
                  eris%ovvo(m, b, e, j) = asym(src, map, m, no + b, no + e, j)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(src, map, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do j = 1, no
         do i = 1, no
            do b = 1, nv
               do m = 1, no
                  eris%ovoo(m, b, i, j) = asym(src, map, m, no + b, i, j)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(src, map, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do f = 1, nv
         do e = 1, nv
            do a = 1, nv
               do m = 1, no
                  eris%ovvv(m, a, e, f) = asym(src, map, m, no + a, no + e, no + f)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(src, map, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do j = 1, no
         do e = 1, nv
            do b = 1, nv
               do a = 1, nv
                  eris%vvvo(a, b, e, j) = asym(src, map, no + a, no + b, no + e, j)
               end do
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine build_cc_eris

   pure function asym2(direct, keep_direct, exchange, keep_exchange) result(v)
      !! <pq||rs> = (pr|qs) - (ps|qr), given the two chemists' integrals
      !!
      !! `asym` for the fitted path, where which block an integral lives in
      !! depends on which of the four indices are occupied and is known only at
      !! the call site. The caller fetches both terms; this applies the signs.
      real(dp), intent(in) :: direct, exchange
      logical, intent(in) :: keep_direct, keep_exchange
         !! Whether the spins of the pair each term couples agree
      real(dp) :: v

      v = 0.0_dp
      if (keep_direct) v = v + direct
      if (keep_exchange) v = v - exchange
   end function asym2

   subroutine build_cc_eris_fitted(chem, map, no, nv, eris)
      !! The antisymmetrised blocks, from the fitted spatial ones
      !!
      !! Block for block the same result as `build_cc_eris`, reached without a
      !! full spatial tensor to slice. Each of the eight is written against
      !! whichever chemists' block its two terms fall in, using the pair symmetry
      !! (pq|rs) = (qp|rs) = (pq|sr) = (rs|pq); the mapping is stated above each
      !! one.
      !!
      !! **This routine indexes `map` with within-space indices, and only the
      !! interleaved ordering makes that legal.** `same_spin(map, n, e)` below
      !! passes a virtual's index `e` in `1..nv` where `map` is indexed by spin
      !! orbital in `1..n_so`, so the lookup lands in the occupied block. Under
      !! the interleaving that is harmless: spins alternate and `no` is even, so
      !! `is_a(e)` and `is_a(no+e)` agree. A blocked map reaches here only
      !! because `run_libcint_ccsd` refuses an unrestricted reference on the
      !! fitted path. **Enabling fitted UCCSD means offsetting every index in
      !! this routine first**; forgetting converges to a wrong number.
      type(so_map_t), intent(in) :: map
      type(cc_chem_t), intent(in) :: chem
      integer, intent(in) :: no, nv
      type(cc_eris_t), intent(out) :: eris

      integer :: i, j, k, l, a, b, e, f, m, n

      allocate (eris%oooo(no, no, no, no), eris%ooov(no, no, no, nv))
      allocate (eris%oovv(no, no, nv, nv), eris%ovov(no, nv, no, nv))
      allocate (eris%ovvo(no, nv, nv, no), eris%ovoo(no, nv, no, no))
      allocate (eris%ovvv(no, nv, nv, nv), eris%vvvo(nv, nv, nv, no))

      ! <ij||kl> = (ik|jl) - (il|jk)
      !$omp parallel do default(none) shared(chem, map, eris, no, nv) &
      !$omp    private(i, j, k, l) schedule(static) collapse(2)
      do l = 1, no
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  eris%oooo(i, j, k, l) = asym2( &
                                          chem%oooo(spatial(map, i), spatial(map, k), spatial(map, j), spatial(map, l)), &
                                          same_spin(map, i, k) .and. same_spin(map, j, l), &
                                          chem%oooo(spatial(map, i), spatial(map, l), spatial(map, j), spatial(map, k)), &
                                          same_spin(map, i, l) .and. same_spin(map, j, k))
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! <mn||ie> = (mi|ne) - (me|ni), and (me|ni) = (ni|me)
      !$omp parallel do default(none) shared(chem, map, eris, no, nv) &
      !$omp    private(m, n, i, e) schedule(static) collapse(2)
      do e = 1, nv
         do i = 1, no
            do n = 1, no
               do m = 1, no
                  eris%ooov(m, n, i, e) = asym2( &
                                          chem%ooov(spatial(map, m), spatial(map, i), spatial(map, n), spatial(map, e)), &
                                          same_spin(map, m, i) .and. same_spin(map, n, e), &
                                          chem%ooov(spatial(map, n), spatial(map, i), spatial(map, m), spatial(map, e)), &
                                          same_spin(map, m, e) .and. same_spin(map, n, i))
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! <ij||ab> = (ia|jb) - (ib|ja)
      !$omp parallel do default(none) shared(chem, map, eris, no, nv) &
      !$omp    private(i, j, a, b) schedule(static) collapse(2)
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  eris%oovv(i, j, a, b) = asym2( &
                                          chem%ovov(spatial(map, i), spatial(map, a), spatial(map, j), spatial(map, b)), &
                                          same_spin(map, i, a) .and. same_spin(map, j, b), &
                                          chem%ovov(spatial(map, i), spatial(map, b), spatial(map, j), spatial(map, a)), &
                                          same_spin(map, i, b) .and. same_spin(map, j, a))
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! <ia||jb> = (ij|ab) - (ib|aj), and (ib|aj) = (ib|ja)
      !$omp parallel do default(none) shared(chem, map, eris, no, nv) &
      !$omp    private(i, j, a, b) schedule(static) collapse(2)
      do b = 1, nv
         do j = 1, no
            do a = 1, nv
               do i = 1, no
                  eris%ovov(i, a, j, b) = asym2( &
                                          chem%oovv(spatial(map, i), spatial(map, j), spatial(map, a), spatial(map, b)), &
                                          same_spin(map, i, j) .and. same_spin(map, a, b), &
                                          chem%ovov(spatial(map, i), spatial(map, b), spatial(map, j), spatial(map, a)), &
                                          same_spin(map, i, b) .and. same_spin(map, a, j))
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! <mb||ej> = (me|bj) - (mj|be), and (me|bj) = (me|jb)
      !$omp parallel do default(none) shared(chem, map, eris, no, nv) &
      !$omp    private(m, b, e, j) schedule(static) collapse(2)
      do j = 1, no
         do e = 1, nv
            do b = 1, nv
               do m = 1, no
                  eris%ovvo(m, b, e, j) = asym2( &
                                          chem%ovov(spatial(map, m), spatial(map, e), spatial(map, j), spatial(map, b)), &
                                          same_spin(map, m, e) .and. same_spin(map, b, j), &
                                          chem%oovv(spatial(map, m), spatial(map, j), spatial(map, b), spatial(map, e)), &
                                          same_spin(map, m, j) .and. same_spin(map, b, e))
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! <mb||ij> = (mi|bj) - (mj|bi), and (mi|bj) = (mi|jb)
      !$omp parallel do default(none) shared(chem, map, eris, no, nv) &
      !$omp    private(m, b, i, j) schedule(static) collapse(2)
      do j = 1, no
         do i = 1, no
            do b = 1, nv
               do m = 1, no
                  eris%ovoo(m, b, i, j) = asym2( &
                                          chem%ooov(spatial(map, m), spatial(map, i), spatial(map, j), spatial(map, b)), &
                                          same_spin(map, m, i) .and. same_spin(map, b, j), &
                                          chem%ooov(spatial(map, m), spatial(map, j), spatial(map, i), spatial(map, b)), &
                                          same_spin(map, m, j) .and. same_spin(map, b, i))
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! <ma||ef> = (me|af) - (mf|ae)
      !$omp parallel do default(none) shared(chem, map, eris, no, nv) &
      !$omp    private(m, a, e, f) schedule(static) collapse(2)
      do f = 1, nv
         do e = 1, nv
            do a = 1, nv
               do m = 1, no
                  eris%ovvv(m, a, e, f) = asym2( &
                                          chem%ovvv(spatial(map, m), spatial(map, e), spatial(map, a), spatial(map, f)), &
                                          same_spin(map, m, e) .and. same_spin(map, a, f), &
                                          chem%ovvv(spatial(map, m), spatial(map, f), spatial(map, a), spatial(map, e)), &
                                          same_spin(map, m, f) .and. same_spin(map, a, e))
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! <ab||ej> = (ae|bj) - (aj|be), and (ae|bj) = (jb|ae), (aj|be) = (ja|be)
      !$omp parallel do default(none) shared(chem, map, eris, no, nv) &
      !$omp    private(a, b, e, j) schedule(static) collapse(2)
      do j = 1, no
         do e = 1, nv
            do b = 1, nv
               do a = 1, nv
                  eris%vvvo(a, b, e, j) = asym2( &
                                          chem%ovvv(spatial(map, j), spatial(map, b), spatial(map, a), spatial(map, e)), &
                                          same_spin(map, a, e) .and. same_spin(map, b, j), &
                                          chem%ovvv(spatial(map, j), spatial(map, a), spatial(map, b), spatial(map, e)), &
                                          same_spin(map, a, j) .and. same_spin(map, b, e))
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine build_cc_eris_fitted

   subroutine spin_orbital_integrals(src, map, n_so, w)
      !! The full antisymmetrised spin-orbital tensor <pq||rs>
      !!
      !! Kept for the antisymmetry test. The calculation reads the blocks
      !! `build_cc_eris` carves instead.
      type(so_source_t), intent(in) :: src
      type(so_map_t), intent(in) :: map
      integer, intent(in) :: n_so
      real(dp), allocatable, intent(out) :: w(:, :, :, :)

      integer :: p, q, r, s

      allocate (w(n_so, n_so, n_so, n_so))
      !$omp parallel do default(none) shared(w, src, map, n_so) private(p, q, r, s) &
      !$omp    schedule(static) collapse(2)
      do s = 1, n_so
         do r = 1, n_so
            do q = 1, n_so
               do p = 1, n_so
                  w(p, q, r, s) = asym(src, map, p, q, r, s)
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine spin_orbital_integrals

   subroutine run_libcint_ccsd(mol, coeff, orbital_energies, n_occ, frozen, &
                               max_iter, energy_tol, want_triples, verbose, &
                               result, error, diis_vectors, aux, &
                               coeff_b, orbital_energies_b, n_occ_b)
      !! Drive CCSD, and optionally (T), to convergence
      !!
      !! **Memory.** The spin-orbital blocks are sixteen times their spatial
      !! size, and the two n_occ n_vir^3 ones set the ceiling. The conventional
      !! path adds the n_act^4 spatial tensor, which the ladder reads (ae|bf) out
      !! of; the fitted path has no such tensor -- the ladder builds what it needs
      !! from `b_vv`, so RI-CCSD reaches further than CCSD rather than merely
      !! differing by the fitting error. (T) holds no triples amplitude either:
      !! one occupied triple at a time, n_vir^3 at once.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)              !! MO coefficients, (n_ao, n_mo)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ                     !! Spatial occupied count
      integer, intent(in) :: frozen                    !! Spatial orbitals to freeze
      integer, intent(in) :: max_iter
      real(dp), intent(in) :: energy_tol
      logical, intent(in) :: want_triples
      logical, intent(in) :: verbose
      type(cc_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: diis_vectors
      type(libcint_molecule_t), intent(in), optional :: aux
         !! Auxiliary basis. Present means every MO integral is density-fitted
         !! rather than transformed exactly -- RI-CCSD. *Every* class is fitted,
         !! not only the particle-particle ladder, which is what PySCF's `dfccsd`
         !! does and what the comparison against it requires.
      real(dp), intent(in), optional :: coeff_b(:, :)
         !! Beta MO coefficients. Present means an unrestricted reference, and
         !! the other two beta arguments must come with it.
      real(dp), intent(in), optional :: orbital_energies_b(:)
      integer, intent(in), optional :: n_occ_b

      real(dp), allocatable :: eri(:, :), mo(:, :, :, :), b_vv(:, :)
      type(cc_eris_t) :: eris
      type(cc_chem_t) :: chem
      real(dp), allocatable :: c_act(:, :), eps(:)
      real(dp), allocatable :: t1(:, :), t2(:, :, :, :), t1n(:, :), t2n(:, :, :, :)
      real(dp), allocatable :: flat(:), err_flat(:)
      type(diis_state_t) :: diis
      logical :: extrapolated
      integer :: n_ao, n_mo, n_act, n_o, n_v, no, nv, n_so, iter, diis_size
      ! TODO(mqc): `i` and `a` are declared and never used here.
      integer :: i, a, s
      real(dp) :: e_corr, e_old, de, t_iter

      character(len=MAX_LINE_LENGTH) :: line

      type(timing_report_t) :: clk
      type(so_map_t) :: map
      type(so_source_t) :: src
      logical :: unrestricted
      integer :: n_ob, n_vb
      real(dp), allocatable :: cb_act(:, :)

      n_ao = size(coeff, 1)
      n_mo = size(coeff, 2)

      if (frozen < 0 .or. frozen >= n_occ) then
         call error%set(ERROR_VALIDATION, "CCSD: the frozen core count must leave at "// &
                        "least one occupied orbital")
         return
      end if
      n_act = n_mo - frozen
      unrestricted = present(coeff_b)
      if (unrestricted .and. present(aux)) then
         ! Refused rather than run. The fitted path builds its blocks and its
         ! ladder straight off `b_vv`, which has no spin blocks, and both still
         ! assume alpha and beta share a spatial orbital. Run that way it would
         ! not fail -- it would converge, to the restricted answer for the alpha
         ! orbitals, which is a wrong number that looks like a right one.
         call error%set(ERROR_VALIDATION, "UCCSD: density fitting is not available "// &
                        "over an unrestricted reference. Run without an auxiliary "// &
                        "basis, or use a restricted reference.")
         return
      end if
      if (unrestricted) then
         if (.not. present(orbital_energies_b) .or. .not. present(n_occ_b)) then
            call error%set(ERROR_VALIDATION, "UCCSD: beta coefficients given without "// &
                           "beta orbital energies or a beta occupation count")
            return
         end if

         ! The beta arrays are sized from the alpha ones below -- `n_mo` comes
         ! from `coeff` and slices `coeff_b` -- so a beta block of a different
         ! width is read past its end rather than reported. The two spins span
         ! the same basis, so this is a mismatch and not a case to support.
         if (size(coeff_b, 1) /= n_ao .or. size(coeff_b, 2) /= n_mo) then
            call error%set(ERROR_VALIDATION, "UCCSD: the beta coefficients are not the "// &
                           "same shape as the alpha ones. Both spins span the same "// &
                           "basis, so this is a bookkeeping error upstream.")
            return
         end if
         if (size(orbital_energies_b) /= n_mo) then
            call error%set(ERROR_VALIDATION, "UCCSD: there is one beta orbital energy "// &
                           "per beta orbital, and the count does not match the "// &
                           "coefficients.")
            return
         end if
         if (n_occ_b < 0 .or. n_occ_b > n_mo) then
            call error%set(ERROR_VALIDATION, "UCCSD: the beta occupation count is "// &
                           "outside the orbital space.")
            return
         end if

         ! **Zero active beta orbitals is legitimate and negative is not.** A
         ! frozen core removes the same *spatial* orbitals from both spins, and
         ! a high-spin system can have no beta electrons left outside it --
         ! frozen-core lithium is one valence electron and it is alpha. So the
         ! bound is `>`, not `>=`.
         !
         ! Past it the failure is not a wrong number. `build_so_map` sizes
         ! itself `n_o + n_ob + n_v + n_vb` and then fills `n_o + n_v + n_vb`
         ! entries, because the beta occupied loop does not run at a negative
         ! count -- so it writes past the end of its own arrays. The alpha
         ! guard above does not catch it: Li+ at multiplicity 3 has two alpha
         ! and no beta electrons, and one frozen core orbital.
         if (frozen > n_occ_b) then
            call error%set(ERROR_VALIDATION, "UCCSD: the frozen core takes more "// &
                           "orbitals than the beta space has occupied. A frozen core "// &
                           "removes the same spatial orbitals from both spins, so it "// &
                           "cannot exceed the smaller occupation. Reduce "// &
                           "keywords.correlation.n_frozen_core, or turn freeze_core off.")
            return
         end if

         n_ob = n_occ_b - frozen
         n_vb = n_mo - n_occ_b
      else
         n_ob = n_occ - frozen
         n_vb = n_mo - n_occ
      end if
      n_o = n_occ - frozen
      n_v = n_mo - n_occ
      if (n_v < 1 .or. n_vb < 1) then
         call error%set(ERROR_VALIDATION, "CCSD: no virtual orbitals -- the basis is "// &
                        "saturated by the occupied space and there is nothing to excite into")
         return
      end if

      no = n_o + n_ob
      nv = n_v + n_vb
      n_so = no + nv

      ! Restricted: the two spins have the same orbitals and the same counts, so
      ! the map interleaves. Unrestricted: blocked, because interleaving stops
      ! keeping occupied contiguous the moment the counts differ.
      call build_so_map(n_o, n_ob, n_v, n_vb, map)

      diis_size = 8
      if (present(diis_vectors)) diis_size = diis_vectors

      if (verbose) then
         write (line, "(a,i0,a,i0,a,i0,a)") "  coupled cluster: ", n_so, " spin orbitals, ", no, &
            " occupied, ", nv, " virtual"
         call logger%info(trim(line))
         if (frozen > 0) then
            write (line, "(a,i0,a)") "  frozen core: ", frozen, " spatial orbitals"
            call logger%info(trim(line))
         end if
         ! Memory is the interesting constraint here, and the comparison worth
         ! reporting is against the whole spin-orbital tensor, n_so^4.
         write (line, "(a,f0.1,a,f0.1,a)") "  integrals: ", integral_megabytes(no, nv), &
            " MB in blocks, against ", real(n_so, dp)**4*8.0_dp/1.0e6_dp, " MB for the full tensor"
         call logger%info(trim(line))
         ! The spatial tensor is the conventional path's alone, and it is the
         ! difference between the two paths' ceilings -- both report the same
         ! block total above.
         if (.not. present(aux)) then
            write (line, "(a,f0.1,a)") "  plus the spatial tensor: ", &
               real(n_act, dp)**4*8.0_dp/1.0e6_dp, " MB, which density fitting does not build"
            call logger%info(trim(line))
         end if
      end if

      ! ---- integrals --------------------------------------------------------
      !
      ! The two paths differ in what they can avoid materialising. Conventional
      ! has to form the n_act^4 spatial tensor -- the ladder reads (ae|bf) out of
      ! it and there is nowhere else to get it. Fitted does not: the five spatial
      ! blocks with an occupied index cover every antisymmetrised block, and the
      ! ladder builds its own from `b_vv`. So `mo` and `b_vv` are alternatives,
      ! never both, and everything downstream is told which by which one is
      ! allocated.
      call clk%start()
      call clk%begin("AO->MO integrals")
      allocate (c_act(n_ao, n_act))
      c_act = coeff(:, frozen + 1:n_mo)
      if (present(aux)) then
         call fitted_chem_blocks(mol, aux, c_act(:, 1:n_o), c_act(:, n_o + 1:n_act), &
                                 chem, b_vv, error)
         if (error%has_error()) return
      else
         call mol%eris_packed(eri)
         if (unrestricted) then
            ! Three tensors, not one: `(aa|aa)`, `(bb|bb)` and the mixed
            ! `(aa|bb)`. `(bb|aa)` is not built -- `tensor` and `chem` reach it
            ! by swapping the electron pairs. This is three times the memory of
            ! the restricted path, and that is inherent: the numbers genuinely
            ! differ once the two spins have different orbitals.
            allocate (cb_act(n_ao, n_act))
            cb_act = coeff_b(:, frozen + 1:n_mo)
            src%restricted = .false.
            call transform_block(eri, c_act, c_act, c_act, c_act, src%aa)
            call transform_block(eri, cb_act, cb_act, cb_act, cb_act, src%bb)
            call transform_block(eri, c_act, c_act, cb_act, cb_act, src%ab)
            deallocate (cb_act)
         else
            call transform_block(eri, c_act, c_act, c_act, c_act, mo)
            ! Moved, not copied: this tensor is the conventional path's memory
            ! ceiling and duplicating it to fill a container would double it.
            ! `mo` is unallocated afterwards and the ladder reads `src%aa`.
            src%restricted = .true.
            call move_alloc(mo, src%aa)
         end if
         deallocate (eri)
      end if
      deallocate (c_act)

      ! Blocks rather than the whole (2 n_act)^4 tensor. On the conventional path
      ! `src` is kept afterwards because the ladder assembles its own from it in
      ! batches; on the fitted one the spatial blocks are finished with here and
      ! only `b_vv` outlives them.
      call clk%lap()
      call clk%begin("CC integral blocks")
      if (present(aux)) then
         call build_cc_eris_fitted(chem, map, no, nv, eris)
         deallocate (chem%oooo, chem%ooov, chem%oovv, chem%ovov, chem%ovvv)
      else
         call build_cc_eris(src, map, no, nv, eris)
      end if
      call clk%lap()

      ! Spin-orbital energies. Canonical, so the Fock matrix is this diagonal and
      ! nothing else -- which is what lets every denominator below be a sum of
      ! four of these rather than a matrix element.
      allocate (eps(n_so))
      do s = 1, n_so
         if (map%is_a(s) .or. .not. unrestricted) then
            eps(s) = orbital_energies(frozen + spatial(map, s))
         else
            eps(s) = orbital_energies_b(frozen + spatial(map, s))
         end if
      end do

      ! ---- MP2, as the checkpoint before any amplitude equation --------------
      call clk%begin("MP2 amplitudes")
      allocate (t1(nv, no), t2(nv, nv, no, no))
      t1 = 0.0_dp
      call mp2_amplitudes(eris, eps, no, nv, t2, result%e_mp2)
      call clk%lap()
      if (verbose) then
         write (line, "(a,f20.12)") "  MP2 (spin orbital) = ", result%e_mp2
         call logger%info(trim(line))
      end if

      ! ---- CCSD -------------------------------------------------------------
      allocate (t1n(nv, no), t2n(nv, nv, no, no))
      allocate (flat(nv*no + nv*nv*no*no), err_flat(nv*no + nv*nv*no*no))
      call diis%init(diis_size, size(flat), size(err_flat))

      e_old = 0.0_dp
      result%converged = .false.
      ! The amplitude update is timed per iteration and printed with the row.
      ! Every iteration costs the same here -- there is no screening that
      ! tightens as things converge -- so the first row is already the answer to
      ! "how long is this going to take", which on a system large enough to be
      ! worth asking about is the only thing anyone wants from the table.
      call convergence_header(verbose, "CCSD iterations", &
                              "    iter                 E_corr          dE   diis       time", 60)
      do iter = 1, max_iter
         t_iter = clk%seconds_of(STAGE_ITER)
         call ccsd_iteration(src, map, b_vv, eris, eps, no, nv, t1, t2, t1n, t2n)
         call clk%lap(STAGE_ITER)
         t_iter = clk%seconds_of(STAGE_ITER) - t_iter

         ! Extrapolate on the amplitudes with the step as the error vector, which
         ! is the same trick `mqc_diis` already plays on the Fock matrix: the
         ! step vanishes exactly at convergence, so it is the right thing to
         ! drive to zero.
         call pack_amplitudes(t1n, t2n, no, nv, flat)
         call pack_step(t1n, t2n, t1, t2, no, nv, err_flat)
         call diis%push(flat, err_flat)
         call diis%extrapolate(flat, extrapolated)
         if (extrapolated) call unpack_amplitudes(flat, no, nv, t1n, t2n)

         t1 = t1n
         t2 = t2n

         call ccsd_energy(eris, t1, t2, no, nv, result%e_singles, result%e_doubles)
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
         call error%set(ERROR_VALIDATION, "CCSD did not converge")
         return
      end if

      ! ---- (T) --------------------------------------------------------------
      if (want_triples) then
         call clk%begin("(T) correction")
         call triples_correction(eris, eps, t1, t2, no, nv, result%e_triples)
         call clk%lap()
         if (verbose) then
            write (line, "(a,f20.12)") "  (T) = ", result%e_triples
            call logger%info(trim(line))
         end if
      end if

      call clk%finish()
      call clk%report("CCSD")

      result%e_correlation = result%e_singles + result%e_doubles + result%e_triples
   end subroutine run_libcint_ccsd

   subroutine df_block_gemm(bl, br, nl, nr, block)
      !! (pq|rs) = sum_P B^P_pq B^P_rs for one pair of fitted blocks
      !!
      !! The product is a matrix over the two compound indices and the caller
      !! wants it with four; the explicit-shape `block` takes the destination by
      !! sequence association, so the gemm writes straight into the rank-four
      !! array rather than into a temporary that is then reshaped.
      real(dp), intent(in) :: bl(:, :), br(:, :)
      integer, intent(in) :: nl, nr
      real(dp), intent(out) :: block(nl, nr)

      call pic_gemm(bl, br, block, transb="T")
   end subroutine df_block_gemm

   subroutine fitted_chem_blocks(mol, aux, c_occ, c_vir, chem, b_vv, error)
      !! The fitted spatial integrals RI-CCSD actually reads
      !!
      !! Every antisymmetrised block the amplitude equations want carries at
      !! least one occupied index -- the sole exception is `<ab||ef>`, which
      !! `particle_ladder` assembles for itself from `b_vv`. So the five blocks
      !! below are the whole of what has to be materialised, and the largest of
      !! them is `(ov|vv)` at n_occ n_vir^3, smaller than the conventional path's
      !! n_act^4 tensor by roughly n_vir/n_occ.
      type(libcint_molecule_t), intent(in) :: mol, aux
      real(dp), intent(in) :: c_occ(:, :), c_vir(:, :)
      type(cc_chem_t), intent(out) :: chem
      real(dp), allocatable, intent(out) :: b_vv(:, :)  !! (P, ab), kept for the ladder
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: b_oo(:, :), b_ov(:, :), b_vv_pq(:, :)
      integer :: n_o, n_v

      n_o = size(c_occ, 2)
      n_v = size(c_vir, 2)

      ! `build_df_mo_block` lays the compound index out left-fastest, so
      ! `b_ov` is (i,a) with i fastest -- which is the ordering every block
      ! below is indexed in.
      ! Every block is consumed as `B B^T`, so the cheap factor is safe.
      call build_df_mo_block(mol, aux, c_occ, c_occ, b_oo, error, fast_factor=.true.)
      if (error%has_error()) return
      call build_df_mo_block(mol, aux, c_occ, c_vir, b_ov, error, fast_factor=.true.)
      if (error%has_error()) return
      call build_df_mo_block(mol, aux, c_vir, c_vir, b_vv_pq, error, fast_factor=.true.)
      if (error%has_error()) return

      allocate (chem%oooo(n_o, n_o, n_o, n_o), chem%ooov(n_o, n_o, n_o, n_v))
      allocate (chem%oovv(n_o, n_o, n_v, n_v), chem%ovov(n_o, n_v, n_o, n_v))
      allocate (chem%ovvv(n_o, n_v, n_v, n_v))

      call df_block_gemm(b_oo, b_oo, n_o*n_o, n_o*n_o, chem%oooo)
      call df_block_gemm(b_oo, b_ov, n_o*n_o, n_o*n_v, chem%ooov)
      call df_block_gemm(b_oo, b_vv_pq, n_o*n_o, n_v*n_v, chem%oovv)
      call df_block_gemm(b_ov, b_ov, n_o*n_v, n_o*n_v, chem%ovov)
      call df_block_gemm(b_ov, b_vv_pq, n_o*n_v, n_v*n_v, chem%ovvv)

      deallocate (b_oo, b_ov)

      ! Transposed because the ladder is the only thing that survives this
      ! routine, and it wants the auxiliary index leading: it slices `b_vv` by
      ! spatial orbital, which is contiguous in (P, ab) and strided the other
      ! way round.
      b_vv = transpose(b_vv_pq)
      deallocate (b_vv_pq)
   end subroutine fitted_chem_blocks

   subroutine mp2_amplitudes(eris, eps, no, nv, t2, e_mp2)
      !! First-order doubles, and the MP2 energy they give
      !!
      !!     t_ij^ab = <ij||ab> / D_ij^ab,   E = 1/4 sum <ij||ab> t_ij^ab
      !!
      !! Must reproduce `run_libcint_mp2` exactly. If it does not, the
      !! antisymmetrisation or the spin-orbital index map is wrong and no
      !! amplitude equation after it can be right.
      type(cc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: eps(:)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: t2(:, :, :, :)
      real(dp), intent(out) :: e_mp2

      integer :: i, j, a, b
      real(dp) :: d

      e_mp2 = 0.0_dp
      !$omp parallel do default(none) shared(eris, eps, t2, no, nv) &
      !$omp    private(i, j, a, b, d) reduction(+:e_mp2) schedule(static) collapse(2)
      do j = 1, no
         do i = 1, no
            do b = 1, nv
               do a = 1, nv
                  d = eps(i) + eps(j) - eps(no + a) - eps(no + b)
                  t2(a, b, i, j) = eris%oovv(i, j, a, b)/d
                  e_mp2 = e_mp2 + 0.25_dp*eris%oovv(i, j, a, b)*t2(a, b, i, j)
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine mp2_amplitudes

   subroutine ccsd_energy(eris, t1, t2, no, nv, e_singles, e_doubles)
      !! E_CCSD = 1/4 sum <ij||ab> t_ij^ab + 1/2 sum <ij||ab> t_i^a t_j^b
      !!
      !! Reported as two numbers because `cc_energy_t` carries them separately.
      !! The f_ia t_i^a term of the general expression is absent: the reference is
      !! canonical, so the occupied-virtual Fock block is zero.
      type(cc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: e_singles, e_doubles

      integer :: i, j, a, b
      real(dp) :: v

      e_singles = 0.0_dp
      e_doubles = 0.0_dp
      !$omp parallel do default(none) shared(eris, t1, t2, no, nv) &
      !$omp    private(i, j, a, b, v) reduction(+:e_singles, e_doubles) &
      !$omp    schedule(static) collapse(2)
      do j = 1, no
         do i = 1, no
            do b = 1, nv
               do a = 1, nv
                  v = eris%oovv(i, j, a, b)
                  e_doubles = e_doubles + 0.25_dp*v*t2(a, b, i, j)
                  e_singles = e_singles + 0.5_dp*v*t1(a, i)*t1(b, j)
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine ccsd_energy

   subroutine ccsd_iteration(src, map, b_vv, eris, eps, no, nv, t1, t2, t1n, t2n)
      !! One CCSD amplitude update, through the Stanton et al. intermediates
      type(so_source_t), intent(in) :: src
         !! The spatial tensors, on the conventional path. Read only by
         !! `particle_ladder`, which assembles <ab||ef> in batches rather than
         !! holding it.
      type(so_map_t), intent(in) :: map
      real(dp), allocatable, intent(in) :: b_vv(:, :)
         !! The fitted three-index virtual block, on the RI path. Exactly one of
         !! this and `src%aa` is allocated; the ladder is the only thing that
         !! cares which.
      type(cc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: eps(:)
      integer, intent(in) :: no, nv
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      real(dp), intent(out) :: t1n(:, :), t2n(:, :, :, :)

      ! Every contraction over `tau` runs over either the particle pair or the
      ! hole pair, never a mixture, so `tau` is held as a matrix over those
      ! compound indices from the start -- ab = (b-1) n_v + a, ij = (j-1) n_o + i.
      ! That is what turns the four terms it appears in into four gemms. `taut`
      ! stays four-index: the intermediates it feeds contract over one particle
      ! and both holes, which is not a matrix product in this layout.
      real(dp), allocatable :: tau2(:, :), taut(:, :, :, :)
      real(dp), allocatable :: fae(:, :), fmi(:, :), fme(:, :)
      real(dp), allocatable :: wmnij(:, :)
      real(dp), allocatable :: oovv2(:, :), lad(:, :)
      real(dp), allocatable :: gae(:, :), gmi(:, :)
      ! The ring term, in the layout that makes it a gemm. See the block that
      ! builds them for why the compound indices are (me) and (bj).
      real(dp), allocatable :: w2(:, :)          !! Wmbej as (me, bj)
      real(dp), allocatable :: t2p(:, :)         !! t2(a,e,i,m) as (ai, me)
      real(dp), allocatable :: ring2(:, :)       !! the ring, as (ai, bj)
      real(dp), allocatable :: oovv_menf(:, :)   !! <mn||ef> as (me, nf)
      real(dp), allocatable :: ttnf(:, :)        !! the damped doubles as (nf, bj)
      real(dp), allocatable :: pz(:, :, :, :)    !! sum_e <mb||ej> t1(e,i)
      real(dp), allocatable :: zz(:, :, :, :)    !! the t1 t1 half of the ring
      ! TODO(mqc): `ab` and `mn` are declared and never used -- left behind when
      ! the compound-index loops became gemms.
      integer :: i, j, m, n, a, b, e, f, ab, ij, mn, ef
      integer :: ai, aj, bi, bj
      real(dp) :: acc, d

      allocate (tau2(nv*nv, no*no), taut(nv, nv, no, no))
      call build_tau(t1, t2, no, nv, tau2, taut)

      ! <mn||ef> as a matrix over the two pairs, which is the shape three of the
      ! gemms below want it in.
      allocate (oovv2(no*no, nv*nv))
      !$omp parallel do default(none) shared(eris, oovv2, no, nv) &
      !$omp    private(m, n, e, f, mn, ef) schedule(static) collapse(2)
      do f = 1, nv
         do e = 1, nv
            ef = (f - 1)*nv + e
            do n = 1, no
               do m = 1, no
                  oovv2((n - 1)*no + m, ef) = eris%oovv(m, n, e, f)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! ---- one-body intermediates (Stanton 3-5) -----------------------------
      allocate (fae(nv, nv), fmi(no, no), fme(no, nv))

      !$omp parallel do default(none) shared(eris, t1, taut, fae, no, nv) &
      !$omp    private(a, e, m, f, n, acc) schedule(static)
      do e = 1, nv
         do a = 1, nv
            acc = 0.0_dp
            do m = 1, no
               do f = 1, nv
                  acc = acc + t1(f, m)*eris%ovvv(m, a, f, e)
               end do
            end do
            do n = 1, no
               do m = 1, no
                  do f = 1, nv
                     acc = acc - 0.5_dp*taut(a, f, m, n)*eris%oovv(m, n, e, f)
                  end do
               end do
            end do
            fae(a, e) = acc
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(eris, t1, taut, fmi, no, nv) &
      !$omp    private(m, i, e, n, f, acc) schedule(static)
      do i = 1, no
         do m = 1, no
            acc = 0.0_dp
            do n = 1, no
               do e = 1, nv
                  acc = acc + t1(e, n)*eris%ooov(m, n, i, e)
               end do
            end do
            do n = 1, no
               do f = 1, nv
                  do e = 1, nv
                     acc = acc + 0.5_dp*taut(e, f, i, n)*eris%oovv(m, n, e, f)
                  end do
               end do
            end do
            fmi(m, i) = acc
         end do
      end do
      !$omp end parallel do

      do e = 1, nv
         do m = 1, no
            acc = 0.0_dp
            do n = 1, no
               do f = 1, nv
                  acc = acc + t1(f, n)*eris%oovv(m, n, e, f)
               end do
            end do
            fme(m, e) = acc
         end do
      end do

      ! ---- two-body intermediates (Stanton 6-8) -----------------------------
      allocate (wmnij(no*no, no*no))

      ! The quadratic term of each ladder intermediate is a gemm, so it goes in
      ! first and the linear terms are added on top of it.
      call pic_gemm(oovv2, tau2, wmnij, alpha=0.25_dp, beta=0.0_dp)
      !$omp parallel do default(none) shared(eris, t1, wmnij, no, nv) &
      !$omp    private(m, n, i, j, e, acc, ij) schedule(static) collapse(2)
      do j = 1, no
         do i = 1, no
            ij = (j - 1)*no + i
            do n = 1, no
               do m = 1, no
                  acc = eris%oooo(m, n, i, j)
                  do e = 1, nv
                     acc = acc + t1(e, j)*eris%ooov(m, n, i, e) - t1(e, i)*eris%ooov(m, n, j, e)
                  end do
                  wmnij((n - 1)*no + m, ij) = wmnij((n - 1)*no + m, ij) + acc
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! ---- both ladders ------------------------------------------------------
      !
      ! The two terms that dominate CCSD: n_o^4 n_v^2 for the holes and
      ! n_v^4 n_o^2 for the particles. Both are gemms; the difference is that the
      ! hole one contracts an intermediate small enough to hold, and the particle
      ! one does not.
      allocate (lad(nv*nv, no*no))
      call pic_gemm(tau2, wmnij, lad, alpha=0.5_dp, beta=0.0_dp)
      call particle_ladder(src, map, b_vv, eris, t1, tau2, oovv2, no, nv, lad)

      ! ---- the ring intermediate (Stanton 8) ---------------------------------
      !
      ! Held as W(me, bj) rather than W(m,b,e,j). Both the term that dominates
      ! building it and the term that reads it contract over the particle-hole
      ! pair (m,e), so in this layout both are gemms; written with the four
      ! indices apart they are n_occ^3 n_vir^3 of scalar strided loads each, and
      ! together they were most of a CCSD iteration.
      allocate (w2(no*nv, nv*no))
      allocate (oovv_menf(no*nv, no*nv), ttnf(no*nv, nv*no))

      !$omp parallel do default(none) shared(eris, oovv_menf, no, nv) &
      !$omp    private(m, n, e, f) schedule(static) collapse(2)
      do f = 1, nv
         do e = 1, nv
            do n = 1, no
               do m = 1, no
                  oovv_menf((e - 1)*no + m, (f - 1)*no + n) = eris%oovv(m, n, e, f)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(t1, t2, ttnf, no, nv) &
      !$omp    private(b, j, n, f) schedule(static) collapse(2)
      do j = 1, no
         do b = 1, nv
            do f = 1, nv
               do n = 1, no
                  ttnf((f - 1)*no + n, (j - 1)*nv + b) = &
                     0.5_dp*t2(f, b, j, n) + t1(f, j)*t1(b, n)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! The quadratic term, which is the expensive one, as a single gemm.
      call pic_gemm(oovv_menf, ttnf, w2, alpha=-1.0_dp, beta=0.0_dp)
      deallocate (oovv_menf, ttnf)

      ! The three linear terms on top of it. Both inner sums are one order of
      ! n_vir cheaper than the one above, so they stay loops.
      !$omp parallel do default(none) shared(eris, t1, w2, no, nv) &
      !$omp    private(m, b, e, j, f, n, acc) schedule(static) collapse(2)
      do j = 1, no
         do b = 1, nv
            do e = 1, nv
               do m = 1, no
                  acc = eris%ovvo(m, b, e, j)
                  do f = 1, nv
                     acc = acc + t1(f, j)*eris%ovvv(m, b, e, f)
                  end do
                  do n = 1, no
                     ! <mn||ej> = <mn||je> with the last pair swapped, hence the
                     ! sign: -(-ooov) = +.
                     acc = acc + t1(b, n)*eris%ooov(m, n, j, e)
                  end do
                  w2((e - 1)*no + m, (j - 1)*nv + b) = &
                     w2((e - 1)*no + m, (j - 1)*nv + b) + acc
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! ---- the ring contraction itself ---------------------------------------
      !
      !     R(a,i,b,j) = sum_me t2(a,e,i,m) W(m,b,e,j)
      !                - sum_me t1(e,i) t1(a,m) <mb||ej>
      !
      ! The T2 equation wants this at four permutations of (ij) and (ab), and all
      ! four are the same object read at permuted indices -- the trick
      ! `triples_block` already plays over the virtuals. Contracting once and
      ! indexing four times is a factor four off the arithmetic before the gemm
      ! is worth anything at all.
      allocate (t2p(nv*no, no*nv), ring2(nv*no, nv*no))

      !$omp parallel do default(none) shared(t2, t2p, no, nv) &
      !$omp    private(m, e, i, a) schedule(static) collapse(2)
      do m = 1, no
         do e = 1, nv
            do i = 1, no
               do a = 1, nv
                  t2p((i - 1)*nv + a, (e - 1)*no + m) = t2(a, e, i, m)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      call pic_gemm(t2p, w2, ring2)
      deallocate (t2p, w2)

      ! The t1 t1 half factorises into two contractions of one index each, which
      ! is n_occ^3 n_vir^2 rather than the n_occ^3 n_vir^3 it costs carried
      ! through the quadruple loop.
      allocate (pz(no, nv, no, no), zz(nv, nv, no, no))

      !$omp parallel do default(none) shared(eris, t1, pz, no, nv) &
      !$omp    private(i, j, b, m, e, acc) schedule(static) collapse(2)
      do i = 1, no
         do j = 1, no
            do b = 1, nv
               do m = 1, no
                  acc = 0.0_dp
                  do e = 1, nv
                     acc = acc + eris%ovvo(m, b, e, j)*t1(e, i)
                  end do
                  pz(m, b, j, i) = acc
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(t1, pz, zz, no, nv) &
      !$omp    private(i, j, b, a, m, acc) schedule(static) collapse(2)
      do i = 1, no
         do j = 1, no
            do b = 1, nv
               do a = 1, nv
                  acc = 0.0_dp
                  do m = 1, no
                     acc = acc + t1(a, m)*pz(m, b, j, i)
                  end do
                  zz(a, b, j, i) = acc
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(ring2, zz, no, nv) &
      !$omp    private(i, j, a, b) schedule(static) collapse(2)
      do j = 1, no
         do b = 1, nv
            do i = 1, no
               do a = 1, nv
                  ring2((i - 1)*nv + a, (j - 1)*nv + b) = &
                     ring2((i - 1)*nv + a, (j - 1)*nv + b) - zz(a, b, j, i)
               end do
            end do
         end do
      end do
      !$omp end parallel do
      deallocate (pz, zz)

      ! ---- T1 (Stanton 1) ---------------------------------------------------
      !$omp parallel do default(none) shared(eris, eps, t1, t2, fae, fmi, fme, t1n, no, nv) &
      !$omp    private(a, i, e, m, n, f, acc, d) schedule(static)
      do i = 1, no
         do a = 1, nv
            acc = 0.0_dp
            do e = 1, nv
               acc = acc + t1(e, i)*fae(a, e)
            end do
            do m = 1, no
               acc = acc - t1(a, m)*fmi(m, i)
            end do
            do m = 1, no
               do e = 1, nv
                  acc = acc + t2(a, e, i, m)*fme(m, e)
               end do
            end do
            do n = 1, no
               do f = 1, nv
                  acc = acc - t1(f, n)*eris%ovov(n, a, i, f)
               end do
            end do
            do m = 1, no
               do f = 1, nv
                  do e = 1, nv
                     acc = acc - 0.5_dp*t2(e, f, i, m)*eris%ovvv(m, a, e, f)
                  end do
               end do
            end do
            do n = 1, no
               do m = 1, no
                  do e = 1, nv
                     ! <nm||ei> = <mn||ie>, two pair swaps.
                     acc = acc - 0.5_dp*t2(a, e, m, n)*eris%ooov(m, n, i, e)
                  end do
               end do
            end do
            d = eps(i) - eps(no + a)
            t1n(a, i) = acc/d
         end do
      end do
      !$omp end parallel do

      ! ---- T2 (Stanton 2) ---------------------------------------------------
      !
      ! The two damped one-body intermediates are hoisted: they depend on the
      ! index being permuted, not on the amplitude being solved for, so building
      ! them once outside is worth n_o n_v of arithmetic each.
      allocate (gae(nv, nv), gmi(no, no))
      do e = 1, nv
         do a = 1, nv
            acc = fae(a, e)
            do m = 1, no
               acc = acc - 0.5_dp*t1(a, m)*fme(m, e)
            end do
            gae(a, e) = acc
         end do
      end do
      do i = 1, no
         do m = 1, no
            acc = fmi(m, i)
            do e = 1, nv
               acc = acc + 0.5_dp*t1(e, i)*fme(m, e)
            end do
            gmi(m, i) = acc
         end do
      end do

      !$omp parallel do default(none) &
      !$omp    shared(eris, eps, t1, t2, lad, ring2, gae, gmi, t2n, no, nv) &
      !$omp    private(a, b, i, j, e, m, acc, d, ai, aj, bi, bj) &
      !$omp    schedule(static) collapse(2)
      do j = 1, no
         do i = 1, no
            do b = 1, nv
               do a = 1, nv
                  acc = eris%oovv(i, j, a, b)

                  ! P(ab) sum_e t2(a,e,i,j) Gae(b,e)
                  do e = 1, nv
                     acc = acc + t2(a, e, i, j)*gae(b, e) - t2(b, e, i, j)*gae(a, e)
                  end do

                  ! -P(ij) sum_m t2(a,b,i,m) Gmi(m,j)
                  do m = 1, no
                     acc = acc - t2(a, b, i, m)*gmi(m, j) + t2(a, b, j, m)*gmi(m, i)
                  end do

                  ! Both ladders, already contracted above.
                  acc = acc + lad((b - 1)*nv + a, (j - 1)*no + i)

                  ! Rings: P(ij)P(ab) over the four index pairings, already
                  ! contracted into `ring2` and only read here.
                  ai = (i - 1)*nv + a
                  aj = (j - 1)*nv + a
                  bi = (i - 1)*nv + b
                  bj = (j - 1)*nv + b
                  acc = acc + ring2(ai, bj) - ring2(aj, bi) &
                        - ring2(bi, aj) + ring2(bj, ai)

                  ! P(ij) sum_e t1(e,i) <ab||ej>
                  do e = 1, nv
                     acc = acc + t1(e, i)*eris%vvvo(a, b, e, j) &
                           - t1(e, j)*eris%vvvo(a, b, e, i)
                  end do

                  ! -P(ab) sum_m t1(a,m) <mb||ij>
                  do m = 1, no
                     acc = acc - t1(a, m)*eris%ovoo(m, b, i, j) + t1(b, m)*eris%ovoo(m, a, i, j)
                  end do

                  d = eps(i) + eps(j) - eps(no + a) - eps(no + b)
                  t2n(a, b, i, j) = acc/d
               end do
            end do
         end do
      end do
      !$omp end parallel do

      deallocate (tau2, taut, fae, fmi, fme, wmnij, gae, gmi)
      deallocate (oovv2, lad, ring2)
   end subroutine ccsd_iteration

   subroutine particle_ladder(src, map, b_vv, eris, t1, tau2, oovv2, no, nv, lad)
      !! The particle-particle ladder, never holding <ab||ef>
      !!
      !!     lad(ab,ij) += 1/2 sum_ef Wabef(ab,ef) tau(ef,ij)
      !!
      !! Wabef is n_vir^4, the largest array coupled cluster asks for, and it is
      !! read exactly once -- by this contraction -- so a block of `ef` columns is
      !! built, contracted and dropped.
      !!
      !! The batching is over `ef` and every operand stays contiguous, which is
      !! the only reason this is still gemms:
      !!
      !!   * `oovv2(:, ef0:ef1)` is a column range of an (mn, ef) matrix;
      !!   * `tau2` is passed whole for the quadratic term;
      !!   * the accumulation needs tau(ef, ij) over the batch, which is a *row*
      !!     range of tau2 and so not contiguous -- hence `tau2t`, the transpose,
      !!     built once per iteration and sliced by column instead.
      !!
      !! The arithmetic is unchanged: the energies are bit-identical to the
      !! version that held the whole array.
      !!
      !! **Where <ab||ef> comes from.** Exactly one of `src%aa` and `b_vv` is
      !! allocated. `src` holds the conventional path's spatial tensors, read at
      !! permuted indices. `b_vv` is the fitted three-index block, and the
      !! spatial (a e | b f) for one pair of spatial (e, f) is then a gemm over
      !! the auxiliary index. Those gemms cost naux n_vir^4 an iteration against
      !! the n_occ^2 n_vir^4 of the contraction they feed.
      type(so_source_t), intent(in) :: src
      type(so_map_t), intent(in) :: map
      real(dp), allocatable, intent(in) :: b_vv(:, :)   !! (P, ab), spatial virtuals
      type(cc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: t1(:, :), tau2(:, :), oovv2(:, :)
      integer, intent(in) :: no, nv
      real(dp), intent(inout) :: lad(:, :)

      real(dp), allocatable :: wblk(:, :), tau2t(:, :), tblk(:, :)
      real(dp), allocatable :: gvv(:, :), gvt(:, :)
      integer :: nv2, no2, ef0, ef1, nb, col, ef, a, b, e, f
      integer :: nocc2, nvs, pa, pb, pe, pf, sa, se, sf, sb, e0, f0
      integer :: pe_have, pf_local, x, y
      logical :: fitted
      real(dp) :: acc

      nv2 = nv*nv
      no2 = no*no
      ! The occupied spin orbitals are the first `no`, and `no` is even, so a
      ! virtual spin orbital `a` carries spin `mod(a,2)` and sits, counting
      ! virtuals only as `b_vv` does, in spatial orbital `(a+1)/2`.
      ! TODO(mqc): `nocc2` is assigned here and read nowhere -- a dead local left
      ! behind when the fitted branch stopped offsetting past the occupied
      ! spatial orbitals.
      nocc2 = no/2
      nvs = nv/2
      fitted = allocated(b_vv)

      allocate (tau2t(no2, nv2))
      tau2t = transpose(tau2)

      allocate (wblk(nv2, LADDER_BATCH))
      allocate (tblk(nv, nv*LADDER_BATCH))

      do ef0 = 1, nv2, LADDER_BATCH
         ef1 = min(ef0 + LADDER_BATCH - 1, nv2)
         nb = ef1 - ef0 + 1

         ! Quadratic term first, so the linear ones add onto it.
         call pic_gemm(tau2, oovv2(:, ef0:ef1), wblk(:, 1:nb), &
                       alpha=0.25_dp, beta=0.0_dp)

         ! The t1 dressing, sum_m t1(b,m) <ma||ef>, is n_occ n_vir^4 per
         ! iteration -- the same order as the contraction this feeds -- and it is
         ! a gemm. `ovvv` is laid out (m, a, e, f) and the batch is a contiguous
         ! range of the compound (ef), so the block starting at (1, 1, e0, f0) is
         ! exactly the (no, nv*nb) matrix the product wants; `ladder_t1_block`
         ! exists only to give that block its rank.
         e0 = mod(ef0 - 1, nv) + 1
         f0 = (ef0 - 1)/nv + 1
         call ladder_t1_block(no, nv, nb, t1, eris%ovvv(1, 1, e0, f0), tblk)

         ! Two loops rather than one with a test in it. They differ only in
         ! where <ab||ef> comes from, but that test sits under n_vir^4
         ! iterations and the two paths want different things hoisted, so the
         ! duplication is the cheaper of the two mistakes.
         if (fitted) then
            !$omp parallel default(none) &
            !$omp    shared(b_vv, tblk, wblk, nv, nb, ef0, nvs) &
            !$omp    private(col, ef, a, b, e, f, pa, pb, pe, pf, sa, se, sf, sb, acc) &
            !$omp    private(gvv, gvt, pe_have, pf_local, x, y)
            ! The fitted (x pe | y pf) for this column, built by the thread
            ! that needs it rather than by a serial pass before the loop:
            ! it is O(naux n_vir^3) an iteration and everything around it is
            ! threaded, so a single gemm for it runs one core against
            ! thirteen idle ones. Each column reads only its own (pe, pf)
            ! block, so there is nothing to share.
            !
            ! `gvt` is its transpose. The direct integral is contiguous in
            ! `a`; the exchange one is the same block with `a` and `b` the
            ! other way round, which read in place is a stride of n_vir on
            ! every one of n_vir^4 elements.
            !
            ! Consecutive columns are consecutive `e`, and two spin orbitals
            ! share a spatial one, so the pair is rebuilt every second
            ! column and not every column.
            allocate (gvv(nvs, nvs), gvt(nvs, nvs))
            pe_have = 0
            pf_local = 0
            !$omp do schedule(static)
            do col = 1, nb
               ef = ef0 + col - 1
               e = mod(ef - 1, nv) + 1
               f = (ef - 1)/nv + 1
               se = mod(e, 2)
               sf = mod(f, 2)
               pe = (e + 1)/2
               pf = (f + 1)/2
               if (pe /= pe_have .or. pf /= pf_local) then
                  call pic_gemm(b_vv(:, (pe - 1)*nvs + 1:pe*nvs), &
                                b_vv(:, (pf - 1)*nvs + 1:pf*nvs), gvv, transa="T")
                  ! Written out rather than `transpose`, which allocates a
                  ! temporary and this runs n_vir^2 times an iteration.
                  do y = 1, nvs
                     do x = 1, nvs
                        gvt(y, x) = gvv(x, y)
                     end do
                  end do
                  pe_have = pe
                  pf_local = pf
               end if
               do b = 1, nv
                  sb = mod(b, 2)
                  pb = (b + 1)/2
                  do a = 1, nv
                     ! <ab||ef> = (ae|bf) - (af|be), the second being
                     ! (be|af): the same block with the two virtuals swapped.
                     pa = (a + 1)/2
                     acc = 0.0_dp
                     if (mod(a, 2) == se .and. sb == sf) acc = acc + gvv(pa, pb)
                     if (mod(a, 2) == sf .and. sb == se) acc = acc - gvt(pa, pb)
                     ! <am||ef> = -<ma||ef>, and likewise for b.
                     acc = acc + tblk(b, (col - 1)*nv + a) - tblk(a, (col - 1)*nv + b)
                     wblk((b - 1)*nv + a, col) = wblk((b - 1)*nv + a, col) + acc
                  end do
               end do
            end do
            !$omp end do
            deallocate (gvv, gvt)
            !$omp end parallel
         else
            !$omp parallel do default(none) &
            !$omp    shared(src, map, tblk, wblk, nv, nb, ef0, no) &
            !$omp    private(col, ef, a, b, e, f, pa, pb, pe, pf, sa, se, sf, sb, acc) &
            !$omp    schedule(static)
            do col = 1, nb
               ef = ef0 + col - 1
               e = mod(ef - 1, nv) + 1
               f = (ef - 1)/nv + 1
               ! Virtual `e` here is spin orbital `no + e`; the map carries both
               ! its spin and its spatial index, the latter already offset past
               ! that spin's occupied orbitals.
               se = merge(1, 0, map%is_a(no + e))
               sf = merge(1, 0, map%is_a(no + f))
               pe = map%spat(no + e)
               pf = map%spat(no + f)
               do b = 1, nv
                  sb = merge(1, 0, map%is_a(no + b))
                  pb = map%spat(no + b)
                  do a = 1, nv
                     ! The one place <ab||ef> is needed, straight from the
                     ! spatial tensor. This is `asym` with the spin
                     ! bookkeeping lifted out of the innermost loop -- the
                     ! same two terms, the same signs.
                     pa = map%spat(no + a)
                     sa = merge(1, 0, map%is_a(no + a))
                     acc = 0.0_dp
                     if (sa == se .and. sb == sf) then
                        acc = acc + tensor(src, sa, sb, pa, pe, pb, pf)
                     end if
                     if (sa == sf .and. sb == se) then
                        acc = acc - tensor(src, sa, sb, pa, pf, pb, pe)
                     end if
                     ! <am||ef> = -<ma||ef>, and likewise for b.
                     acc = acc + tblk(b, (col - 1)*nv + a) - tblk(a, (col - 1)*nv + b)
                     wblk((b - 1)*nv + a, col) = wblk((b - 1)*nv + a, col) + acc
                  end do
               end do
            end do
            !$omp end parallel do
         end if

         call pic_gemm(wblk(:, 1:nb), tau2t(:, ef0:ef1), lad, &
                       transb="T", alpha=0.5_dp, beta=1.0_dp)
      end do

      deallocate (wblk, tau2t, tblk)
   end subroutine particle_ladder

   subroutine ladder_t1_block(no, nv, nb, t1, ovvv_blk, tmat)
      !! sum_m t1(b,m) <ma||ef> over one batch of the compound (ef)
      !!
      !! One gemm. `eris%ovvv` is rank four and the product wants it as the
      !! (no, nv*nb) matrix its memory already is; the explicit-shape dummy takes
      !! the block by sequence association, so nothing is copied.
      integer, intent(in) :: no, nv, nb
      real(dp), intent(in) :: t1(nv, no)
      real(dp), intent(in) :: ovvv_blk(no, nv*nb)   !! <ma||ef>, (m, (a,ef))
      real(dp), intent(out) :: tmat(nv, nv*nb)      !! (b, (a,ef))

      call pic_gemm(t1, ovvv_blk, tmat)
   end subroutine ladder_t1_block

   subroutine build_tau(t1, t2, no, nv, tau2, taut)
      !! The two effective doubles the intermediates are written in
      !!
      !!     tau      = t2 + t1 t1 - t1 t1   (transposed)
      !!     tau~     = t2 + 1/2 (t1 t1 - t1 t1)
      !!
      !! They differ only in that half, and mixing them up still converges -- to
      !! the wrong energy. Built together so the difference is visible in one
      !! place.
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: tau2(:, :)          !! (ab, ij)
      real(dp), intent(out) :: taut(:, :, :, :)

      integer :: i, j, a, b
      real(dp) :: cross

      !$omp parallel do default(none) shared(t1, t2, tau2, taut, no, nv) &
      !$omp    private(i, j, a, b, cross) schedule(static) collapse(2)
      do j = 1, no
         do i = 1, no
            do b = 1, nv
               do a = 1, nv
                  cross = t1(a, i)*t1(b, j) - t1(b, i)*t1(a, j)
                  tau2((b - 1)*nv + a, (j - 1)*no + i) = t2(a, b, i, j) + cross
                  taut(a, b, i, j) = t2(a, b, i, j) + 0.5_dp*cross
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine build_tau

   subroutine triples_correction(eris, eps, t1, t2, no, nv, e_triples)
      !! (T), the perturbative triples correction
      !!
      !! Non-iterative, so one pass after CCSD converges, and n_occ^3 n_vir^4.
      !! No triples amplitude is stored: the two n_vir^3 blocks are rebuilt for
      !! each occupied triple, where a full t3 would be n_o^3 n_v^3.
      !!
      !!     E(T) = 1/36 sum t3c D (t3c + t3d)
      !!
      !! with the connected and disconnected triples
      !!
      !!     t3d D = P(i/jk) P(a/bc) t1(a,i) <jk||bc>
      !!     t3c D = P(i/jk) P(a/bc) [ sum_e t2(a,e,j,k) <ei||bc>
      !!                               - sum_m t2(b,c,i,m) <ma||jk> ]
      !!
      !! Both numerators are fully antisymmetric in (i,j,k), so the summand is
      !! symmetric under permuting them and vanishes when two coincide. Only
      !! i > j > k is visited, weighted by six. The virtual permutations are
      !! *not* restricted: P(a/bc) is applied inside `triples_block`, which needs
      !! the whole (a,b,c) cube to read at permuted indices.
      type(cc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: eps(:), t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: e_triples

      real(dp), allocatable :: t3c(:, :, :), t3d(:, :, :)
      real(dp), allocatable :: ovvv_p(:, :, :), t2bc(:, :, :)
      real(dp), allocatable :: ovoo_p(:, :, :, :), oovv_p(:, :, :)
      real(dp), allocatable :: acc(:, :), bcc(:, :), gg(:, :, :), dd(:, :, :)
      integer :: i, j, k, a, b, c
      real(dp) :: d, d_occ, e_local

      e_triples = 0.0_dp

      ! Four integral blocks and one amplitude block, repacked so the
      ! contractions below are gemms over contiguous memory rather than strided
      ! reads out of the full spin-orbital tensor. The largest is n_o n_v^3 and
      ! each is touched n_occ^3 times, so packing once is the whole point.
      call pack_triples_blocks(eris, t2, no, nv, ovvv_p, t2bc, ovoo_p, oovv_p)

      !$omp parallel default(none) &
      !$omp    shared(eps, t1, t2, no, nv, ovvv_p, t2bc, ovoo_p, oovv_p) &
      !$omp    private(i, j, k, a, b, c, d, d_occ, t3c, t3d, e_local, acc, bcc, gg, dd) &
      !$omp    reduction(+:e_triples)
      allocate (t3c(nv, nv, nv), t3d(nv, nv, nv))
      allocate (acc(nv, nv*nv), bcc(nv*nv, nv), gg(nv, nv, nv), dd(nv, nv, nv))
      !$omp do schedule(dynamic) collapse(2)
      do i = 1, no
         do j = 1, no
            if (j >= i) cycle
            do k = 1, j - 1
               call triples_block(t1, t2, no, nv, i, j, k, ovvv_p, t2bc, ovoo_p, &
                                  oovv_p, acc, bcc, gg, dd, t3c, t3d)
               e_local = 0.0_dp
               d_occ = eps(i) + eps(j) + eps(k)
               do c = 1, nv
                  do b = 1, nv
                     do a = 1, nv
                        d = d_occ - eps(no + a) - eps(no + b) - eps(no + c)
                        ! t3c and t3d arrive as numerators; dividing here rather
                        ! than twice inside the block keeps the block additive.
                        e_local = e_local + t3c(a, b, c)*(t3c(a, b, c) + t3d(a, b, c))/d
                     end do
                  end do
               end do
               ! Six, for the six orderings of {i,j,k} this one stands in for.
               e_triples = e_triples + 6.0_dp*e_local/36.0_dp
            end do
         end do
      end do
      !$omp end do
      deallocate (t3c, t3d, acc, bcc, gg, dd)
      !$omp end parallel

      deallocate (ovvv_p, t2bc, ovoo_p, oovv_p)
   end subroutine triples_correction

   subroutine pack_triples_blocks(eris, t2, no, nv, ovvv_p, t2bc, ovoo_p, oovv_p)
      !! Repack what (T) contracts, so each contraction is a gemm
      !!
      !! The compound index is always `bc = (c-1) n_v + b`, which puts it in the
      !! position a gemm wants: leading for the amplitude block it contracts over,
      !! trailing for the integral block it indexes.
      type(cc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), allocatable, intent(out) :: ovvv_p(:, :, :)    !! (e, bc, i) = <ie||bc>
      real(dp), allocatable, intent(out) :: t2bc(:, :, :)      !! (bc, m, i) = t2(b,c,i,m)
      real(dp), allocatable, intent(out) :: ovoo_p(:, :, :, :)  !! (m, a, j, k) = <ma||jk>
      real(dp), allocatable, intent(out) :: oovv_p(:, :, :)    !! (bc, j, k) = <jk||bc>

      integer :: i, j, k, m, a, b, c, bc

      allocate (ovvv_p(nv, nv*nv, no), t2bc(nv*nv, no, no))
      allocate (ovoo_p(no, nv, no, no), oovv_p(nv*nv, no, no))

      !$omp parallel do default(none) shared(eris, ovvv_p, no, nv) &
      !$omp    private(i, b, c, bc) schedule(static)
      do i = 1, no
         do c = 1, nv
            do b = 1, nv
               bc = (c - 1)*nv + b
               ovvv_p(:, bc, i) = eris%ovvv(i, :, b, c)
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(t2, t2bc, no, nv) &
      !$omp    private(i, m, b, c, bc) schedule(static) collapse(2)
      do i = 1, no
         do m = 1, no
            do c = 1, nv
               do b = 1, nv
                  bc = (c - 1)*nv + b
                  t2bc(bc, m, i) = t2(b, c, i, m)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(eris, ovoo_p, oovv_p, no, nv) &
      !$omp    private(j, k, m, a, b, c, bc) schedule(static) collapse(2)
      do k = 1, no
         do j = 1, no
            do a = 1, nv
               do m = 1, no
                  ovoo_p(m, a, j, k) = eris%ovoo(m, a, j, k)
               end do
            end do
            do c = 1, nv
               do b = 1, nv
                  bc = (c - 1)*nv + b
                  oovv_p(bc, j, k) = eris%oovv(j, k, b, c)
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine pack_triples_blocks

   subroutine triples_block(t1, t2, no, nv, i, j, k, ovvv_p, t2bc, ovoo_p, oovv_p, &
                            acc, bcc, gg, dd, t3c, t3d)
      !! The connected and disconnected triples numerators for one (i,j,k)
      !!
      !! P(i/jk) is the three-term permutation f(ijk) - f(jik) - f(kji), and
      !! P(a/bc) likewise, so each numerator is nine signed contributions. For a
      !! fixed occupied permutation the quantity permuted over the virtuals,
      !!
      !!     G(a,b,c) = -sum_e t2(a,e,j,k) <ie||bc> - sum_m t2(b,c,i,m) <ma||jk>
      !!
      !! is one object over all (a,b,c), and P(a/bc) then only reads it at
      !! permuted indices. So three occupied permutations give three pairs of
      !! gemms, and the virtual permutations cost nothing but the indexing.
      !!
      !! The scratch arrays come in from the caller rather than being allocated
      !! here, because this runs n_occ^3 times.
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv, i, j, k
      real(dp), intent(in) :: ovvv_p(:, :, :), t2bc(:, :, :)
      real(dp), intent(in) :: ovoo_p(:, :, :, :), oovv_p(:, :, :)
      real(dp), intent(inout) :: acc(:, :), bcc(:, :)    !! Scratch, (nv, nv^2) and (nv^2, nv)
      real(dp), intent(inout) :: gg(:, :, :), dd(:, :, :)  !! Scratch, (nv, nv, nv)
      real(dp), intent(out) :: t3c(:, :, :), t3d(:, :, :)

      integer, parameter :: N_PERM = 3
      integer :: oi(N_PERM)
      integer :: po, a, b, c, bc, pi, pj, pk
      real(dp) :: sign_o

      t3c = 0.0_dp
      t3d = 0.0_dp

      do po = 1, N_PERM
         ! P(i/jk): identity, then i<->j, then i<->k, the last two with a minus.
         select case (po)
         case (1)
            oi = [i, j, k]
            sign_o = 1.0_dp
         case (2)
            oi = [j, i, k]
            sign_o = -1.0_dp
         case default
            oi = [k, j, i]
            sign_o = -1.0_dp
         end select
         pi = oi(1)
         pj = oi(2)
         pk = oi(3)

         ! -sum_e t2(a,e,j,k) <ie||bc>, as (nv,nv) x (nv,nv^2). The amplitude
         ! slice t2(:,:,pj,pk) is already contiguous, so it goes in as it lies.
         call pic_gemm(t2(:, :, pj, pk), ovvv_p(:, :, pi), acc, alpha=-1.0_dp, beta=0.0_dp)

         ! -sum_m t2(b,c,i,m) <ma||jk>, as (nv^2,no) x (no,nv).
         call pic_gemm(t2bc(:, :, pi), ovoo_p(:, :, pj, pk), bcc, &
                       alpha=-1.0_dp, beta=0.0_dp)

         do c = 1, nv
            do b = 1, nv
               bc = (c - 1)*nv + b
               do a = 1, nv
                  gg(a, b, c) = acc(a, bc) + bcc(bc, a)
                  dd(a, b, c) = t1(a, pi)*oovv_p(bc, pj, pk)
               end do
            end do
         end do

         ! P(a/bc) reads the same two arrays at permuted indices.
         do c = 1, nv
            do b = 1, nv
               do a = 1, nv
                  t3c(a, b, c) = t3c(a, b, c) + sign_o &
                                 *(gg(a, b, c) - gg(b, a, c) - gg(c, b, a))
                  t3d(a, b, c) = t3d(a, b, c) + sign_o &
                                 *(dd(a, b, c) - dd(b, a, c) - dd(c, b, a))
               end do
            end do
         end do
      end do
   end subroutine triples_block

   subroutine pack_amplitudes(t1, t2, no, nv, flat)
      !! t1 and t2 laid end to end, for DIIS
      real(dp), intent(in) :: t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: flat(:)

      integer :: n1

      n1 = nv*no
      flat(1:n1) = reshape(t1, [n1])
      flat(n1 + 1:n1 + nv*nv*no*no) = reshape(t2, [nv*nv*no*no])
   end subroutine pack_amplitudes

   subroutine unpack_amplitudes(flat, no, nv, t1, t2)
      !! The inverse of `pack_amplitudes`
      real(dp), intent(in) :: flat(:)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: t1(:, :), t2(:, :, :, :)

      integer :: n1

      n1 = nv*no
      t1 = reshape(flat(1:n1), [nv, no])
      t2 = reshape(flat(n1 + 1:n1 + nv*nv*no*no), [nv, nv, no, no])
   end subroutine unpack_amplitudes

   subroutine pack_step(t1n, t2n, t1, t2, no, nv, flat)
      !! The change in the amplitudes, which is the DIIS error vector
      !!
      !! It vanishes exactly at a fixed point of the amplitude equations, so it
      !! is the right quantity to extrapolate against -- the same argument as the
      !! FDS - SDF commutator in the SCF.
      real(dp), intent(in) :: t1n(:, :), t2n(:, :, :, :), t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: flat(:)

      integer :: n1

      n1 = nv*no
      flat(1:n1) = reshape(t1n - t1, [n1])
      flat(n1 + 1:n1 + nv*nv*no*no) = reshape(t2n - t2, [nv*nv*no*no])
   end subroutine pack_step

end module mqc_libcint_cc
