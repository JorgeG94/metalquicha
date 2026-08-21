!! Coupled cluster on the CPU, in the spin-orbital basis
module mqc_libcint_cc
   !! CCSD and CCSD(T) over an RHF reference, written in spin orbitals.
   !!
   !! Spin orbitals rather than a spin-adapted closed-shell formulation, and the
   !! arithmetic argument loses on purpose. This is the reference path -- the same
   !! one where conventional MP2 exists so RI-MP2 has something exact to be
   !! checked against -- so one set of equations out of any textbook, correct and
   !! checkable term by term, is worth more than four times the speed on the half
   !! of cases that are closed shell. Revisit if CCSD becomes something people
   !! run rather than something people check against.
   !!
   !! The consequence is that everything here is indexed by spin orbital. Active
   !! spatial orbital `p` becomes spin orbitals `2p-1` (alpha) and `2p` (beta), so
   !! `spatial(s) = (s+1)/2` and `is_alpha(s) = mod(s,2) == 1`. Interleaving
   !! rather than blocking the spins is what keeps occupied and virtual
   !! contiguous: the first `2*n_occ` spin orbitals are occupied, whatever the
   !! spin pattern.
   !!
   !! The antisymmetrised integrals are
   !!
   !!     <pq||rs> = (pr|qs) - (ps|qr)
   !!
   !! with each spatial integral surviving only if the spins on the two indices
   !! it pairs agree. Note the index reordering: the transform produces chemists'
   !! notation and the equations are written in physicists'. That conversion is
   !! the single most likely place for a sign or index error in the whole module,
   !! so it happens in exactly one function -- `asym` -- and that function is unit
   !! tested on its own antisymmetry.
   !!
   !! The amplitude equations follow Stanton, Gauss, Watts and Bartlett (JCP 94,
   !! 4334 (1991)) with their one- and two-body intermediates. Written without
   !! them the equations are both unreadable and slow; written with them they are
   !! neither, and the table in that paper is worth following literally.
   !!
   !! **Frozen core** drops orbitals from the active space before anything else
   !! happens, so nothing below this line knows they existed.
   !!
   !! **What this is not.** Every antisymmetrised block is materialised in spin
   !! orbitals, so `<ma||ef>` and `<ab||ej>` at n_occ n_vir^3 are the wall --
   !! sixteen times their spatial size, which is the price of the formulation
   !! above and is paid on both paths. What is *not* materialised is `<ab||cd>`,
   !! assembled in batches by `particle_ladder`, and on the fitted path the
   !! n_act^4 spatial tensor with it: RI-CCSD holds five spatial blocks with an
   !! occupied index and the three-index `b_vv`, and nothing of fourth order in
   !! the virtuals. See the memory note in `run_libcint_ccsd`.
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

   !> Columns of the compound (ef) index assembled per pass of the
   !> particle-particle ladder. The array held is n_vir^2 by this, so the memory
   !> is O(n_vir^2) rather than the O(n_vir^4) of holding <ab||ef> whole -- which
   !> is the entire point, and only true if this stays a fixed count rather than a
   !> fraction of n_vir^2.
   !>
   !> A larger value opens fewer parallel regions and at small n_vir that is
   !> visible, because each pass then has too little arithmetic to amortise the
   !> barrier. It stays small anyway: the sizes where the memory matters are also
   !> the sizes where a pass has plenty of work.
   integer, parameter :: LADDER_BATCH = 256

   !> The stage name the amplitude update accumulates under. Named because the
   !> per-iteration row reads it back with `seconds_of` to print one iteration's
   !> time, and a typo between the two would silently print zero.
   character(len=*), parameter :: STAGE_ITER = "CCSD iterations"

   type :: cc_eris_t
      !! Antisymmetrised integrals, by block
      !!
      !! `vvvv` is deliberately absent. It is n_vir^4, it is read in exactly one
      !! place -- the particle-particle ladder -- and it is assembled there in
      !! batches straight from the spatial tensor. Keeping it out of this type is
      !! what makes the memory argument in `build_cc_eris` true.
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
      !! for the same reason `cc_eris_t` has no `vvvv`: it is read in one place
      !! and assembled there, here from `b_vv` rather than from a stored tensor.
      !! Together those two absences are what keeps RI-CCSD off the n_act^4
      !! ceiling the conventional path sits on.
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
         !! MP2 from these same spin-orbital integrals. Free -- it is the first
         !! iteration's energy -- and the sharpest check available on the
         !! antisymmetrisation and the denominators, because `run_libcint_mp2`
         !! already agrees with PySCF and this must reproduce it exactly.
      real(dp) :: e_correlation = 0.0_dp  !! singles + doubles + triples
      integer :: iterations = 0
      logical :: converged = .false.
   end type cc_result_t

contains

   pure function integral_megabytes(no, nv) result(mb)
      !! What the integral blocks and the ladder batch cost, in MB
      !!
      !! Everything except the particle-particle ladder, which is never held --
      !! its contribution is the batch, n_vir^2 by EF_BATCH, and that is the whole
      !! reason this number is n_occ n_vir^3 rather than n_vir^4.
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

   pure function spatial(s) result(p)
      !! Which spatial orbital a spin orbital belongs to
      integer, intent(in) :: s
      integer :: p
      p = (s + 1)/2
   end function spatial

   pure function same_spin(s, t) result(same)
      !! Whether two spin orbitals carry the same spin
      integer, intent(in) :: s, t
      logical :: same
      same = mod(s, 2) == mod(t, 2)
   end function same_spin

   pure function asym(mo, p, q, r, s) result(v)
      !! <pq||rs> from the spatial MO tensor, in spin orbitals
      !!
      !! The one place chemists' notation becomes physicists'. `mo` is
      !! (pq|rs) as the transform produced it, so
      !!
      !!     <pq||rs> = (pr|qs) - (ps|qr)
      !!
      !! and a spatial integral contributes only when both index pairs it
      !! couples share a spin -- (pr|qs) needs spin(p) == spin(r) and
      !! spin(q) == spin(s), which is the Kronecker delta the spin integration
      !! leaves behind.
      real(dp), intent(in) :: mo(:, :, :, :)
      integer, intent(in) :: p, q, r, s
      real(dp) :: v

      v = 0.0_dp
      if (same_spin(p, r) .and. same_spin(q, s)) then
         v = v + mo(spatial(p), spatial(r), spatial(q), spatial(s))
      end if
      if (same_spin(p, s) .and. same_spin(q, r)) then
         v = v - mo(spatial(p), spatial(s), spatial(q), spatial(r))
      end if
   end function asym

   subroutine build_cc_eris(mo, no, nv, eris)
      !! The antisymmetrised blocks the amplitude equations actually read
      !!
      !! Replaces the single (2 n_act)^4 tensor. Every block below is a slice of
      !! it, and the slice that is missing -- `vvvv` -- is the one that dominates:
      !! at ten occupied and thirty-eight virtual spin orbitals the whole tensor
      !! is 42 MB, of which the blocks here are 13 MB and `vvvv` alone is 17 MB
      !! again as `wabef`. Assembling `vvvv` in batches instead (see
      !! `ladder_from_batches`) takes the total from 59 MB to 14, and changes the
      !! scaling from n_vir^4 to n_occ n_vir^3 -- a factor n_vir/n_occ, which is
      !! the ratio that grows with basis size at fixed electron count.
      !!
      !! `spin_orbital_integrals` is kept for the tests, which check antisymmetry
      !! over the whole tensor. Nothing in the calculation builds it any more.
      real(dp), intent(in) :: mo(:, :, :, :)
      integer, intent(in) :: no, nv
      type(cc_eris_t), intent(out) :: eris

      integer :: i, j, k, l, a, b, e, f, m, n

      allocate (eris%oooo(no, no, no, no), eris%ooov(no, no, no, nv))
      allocate (eris%oovv(no, no, nv, nv), eris%ovov(no, nv, no, nv))
      allocate (eris%ovvo(no, nv, nv, no), eris%ovoo(no, nv, no, no))
      allocate (eris%ovvv(no, nv, nv, nv), eris%vvvo(nv, nv, nv, no))

      !$omp parallel do default(none) shared(mo, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do l = 1, no
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  eris%oooo(i, j, k, l) = asym(mo, i, j, k, l)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(mo, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do e = 1, nv
         do i = 1, no
            do n = 1, no
               do m = 1, no
                  eris%ooov(m, n, i, e) = asym(mo, m, n, i, no + e)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(mo, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  eris%oovv(i, j, a, b) = asym(mo, i, j, no + a, no + b)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(mo, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do b = 1, nv
         do j = 1, no
            do a = 1, nv
               do i = 1, no
                  eris%ovov(i, a, j, b) = asym(mo, i, no + a, j, no + b)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(mo, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do j = 1, no
         do e = 1, nv
            do b = 1, nv
               do m = 1, no
                  eris%ovvo(m, b, e, j) = asym(mo, m, no + b, no + e, j)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(mo, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do j = 1, no
         do i = 1, no
            do b = 1, nv
               do m = 1, no
                  eris%ovoo(m, b, i, j) = asym(mo, m, no + b, i, j)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(mo, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do f = 1, nv
         do e = 1, nv
            do a = 1, nv
               do m = 1, no
                  eris%ovvv(m, a, e, f) = asym(mo, m, no + a, no + e, no + f)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(mo, eris, no, nv) &
      !$omp    private(i, j, k, l, a, b, e, f, m, n) schedule(static) collapse(2)
      do j = 1, no
         do e = 1, nv
            do b = 1, nv
               do a = 1, nv
                  eris%vvvo(a, b, e, j) = asym(mo, no + a, no + b, no + e, j)
               end do
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine build_cc_eris

   pure function asym2(direct, keep_direct, exchange, keep_exchange) result(v)
      !! <pq||rs> = (pr|qs) - (ps|qr), given the two chemists' integrals
      !!
      !! `asym` for the fitted path. It cannot look the integrals up itself --
      !! which block they live in depends on which of the four indices are
      !! occupied, and that is known at the call site and nowhere else -- so the
      !! caller fetches and this decides. What stays here is the part worth
      !! having in one place: which term is the direct one, which the exchange,
      !! and that the exchange carries the minus.
      real(dp), intent(in) :: direct, exchange
      logical, intent(in) :: keep_direct, keep_exchange
         !! Whether the spins of the pair each term couples agree
      real(dp) :: v

      v = 0.0_dp
      if (keep_direct) v = v + direct
      if (keep_exchange) v = v - exchange
   end function asym2

   subroutine build_cc_eris_fitted(chem, no, nv, eris)
      !! The antisymmetrised blocks, from the fitted spatial ones
      !!
      !! Block for block the same result as `build_cc_eris`, reached without a
      !! full spatial tensor to slice. Each of the eight is written against
      !! whichever chemists' block its two terms fall in, using the pair
      !! symmetry (pq|rs) = (qp|rs) = (pq|sr) = (rs|pq) to put the occupied
      !! indices where the stored block has them; the mapping is stated above
      !! each one because it is the only place it can be checked by reading.
      !!
      !! Occupied spin orbital `i` sits in spatial occupied orbital `(i+1)/2`
      !! and virtual spin orbital `a` in spatial virtual orbital `(a+1)/2`, both
      !! with spin `mod(index,2)` -- `no` is even, so a virtual's spin is the
      !! same whether it is counted from one or from `no+1`, which is what lets
      !! `same_spin` be used on the within-space indices directly.
      type(cc_chem_t), intent(in) :: chem
      integer, intent(in) :: no, nv
      type(cc_eris_t), intent(out) :: eris

      integer :: i, j, k, l, a, b, e, f, m, n

      allocate (eris%oooo(no, no, no, no), eris%ooov(no, no, no, nv))
      allocate (eris%oovv(no, no, nv, nv), eris%ovov(no, nv, no, nv))
      allocate (eris%ovvo(no, nv, nv, no), eris%ovoo(no, nv, no, no))
      allocate (eris%ovvv(no, nv, nv, nv), eris%vvvo(nv, nv, nv, no))

      ! <ij||kl> = (ik|jl) - (il|jk)
      !$omp parallel do default(none) shared(chem, eris, no, nv) &
      !$omp    private(i, j, k, l) schedule(static) collapse(2)
      do l = 1, no
         do k = 1, no
            do j = 1, no
               do i = 1, no
                  eris%oooo(i, j, k, l) = asym2( &
                                          chem%oooo(spatial(i), spatial(k), spatial(j), spatial(l)), &
                                          same_spin(i, k) .and. same_spin(j, l), &
                                          chem%oooo(spatial(i), spatial(l), spatial(j), spatial(k)), &
                                          same_spin(i, l) .and. same_spin(j, k))
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! <mn||ie> = (mi|ne) - (me|ni), and (me|ni) = (ni|me)
      !$omp parallel do default(none) shared(chem, eris, no, nv) &
      !$omp    private(m, n, i, e) schedule(static) collapse(2)
      do e = 1, nv
         do i = 1, no
            do n = 1, no
               do m = 1, no
                  eris%ooov(m, n, i, e) = asym2( &
                                          chem%ooov(spatial(m), spatial(i), spatial(n), spatial(e)), &
                                          same_spin(m, i) .and. same_spin(n, e), &
                                          chem%ooov(spatial(n), spatial(i), spatial(m), spatial(e)), &
                                          same_spin(m, e) .and. same_spin(n, i))
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! <ij||ab> = (ia|jb) - (ib|ja)
      !$omp parallel do default(none) shared(chem, eris, no, nv) &
      !$omp    private(i, j, a, b) schedule(static) collapse(2)
      do b = 1, nv
         do a = 1, nv
            do j = 1, no
               do i = 1, no
                  eris%oovv(i, j, a, b) = asym2( &
                                          chem%ovov(spatial(i), spatial(a), spatial(j), spatial(b)), &
                                          same_spin(i, a) .and. same_spin(j, b), &
                                          chem%ovov(spatial(i), spatial(b), spatial(j), spatial(a)), &
                                          same_spin(i, b) .and. same_spin(j, a))
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! <ia||jb> = (ij|ab) - (ib|aj), and (ib|aj) = (ib|ja)
      !$omp parallel do default(none) shared(chem, eris, no, nv) &
      !$omp    private(i, j, a, b) schedule(static) collapse(2)
      do b = 1, nv
         do j = 1, no
            do a = 1, nv
               do i = 1, no
                  eris%ovov(i, a, j, b) = asym2( &
                                          chem%oovv(spatial(i), spatial(j), spatial(a), spatial(b)), &
                                          same_spin(i, j) .and. same_spin(a, b), &
                                          chem%ovov(spatial(i), spatial(b), spatial(j), spatial(a)), &
                                          same_spin(i, b) .and. same_spin(a, j))
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! <mb||ej> = (me|bj) - (mj|be), and (me|bj) = (me|jb)
      !$omp parallel do default(none) shared(chem, eris, no, nv) &
      !$omp    private(m, b, e, j) schedule(static) collapse(2)
      do j = 1, no
         do e = 1, nv
            do b = 1, nv
               do m = 1, no
                  eris%ovvo(m, b, e, j) = asym2( &
                                          chem%ovov(spatial(m), spatial(e), spatial(j), spatial(b)), &
                                          same_spin(m, e) .and. same_spin(b, j), &
                                          chem%oovv(spatial(m), spatial(j), spatial(b), spatial(e)), &
                                          same_spin(m, j) .and. same_spin(b, e))
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! <mb||ij> = (mi|bj) - (mj|bi), and (mi|bj) = (mi|jb)
      !$omp parallel do default(none) shared(chem, eris, no, nv) &
      !$omp    private(m, b, i, j) schedule(static) collapse(2)
      do j = 1, no
         do i = 1, no
            do b = 1, nv
               do m = 1, no
                  eris%ovoo(m, b, i, j) = asym2( &
                                          chem%ooov(spatial(m), spatial(i), spatial(j), spatial(b)), &
                                          same_spin(m, i) .and. same_spin(b, j), &
                                          chem%ooov(spatial(m), spatial(j), spatial(i), spatial(b)), &
                                          same_spin(m, j) .and. same_spin(b, i))
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! <ma||ef> = (me|af) - (mf|ae)
      !$omp parallel do default(none) shared(chem, eris, no, nv) &
      !$omp    private(m, a, e, f) schedule(static) collapse(2)
      do f = 1, nv
         do e = 1, nv
            do a = 1, nv
               do m = 1, no
                  eris%ovvv(m, a, e, f) = asym2( &
                                          chem%ovvv(spatial(m), spatial(e), spatial(a), spatial(f)), &
                                          same_spin(m, e) .and. same_spin(a, f), &
                                          chem%ovvv(spatial(m), spatial(f), spatial(a), spatial(e)), &
                                          same_spin(m, f) .and. same_spin(a, e))
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! <ab||ej> = (ae|bj) - (aj|be), and (ae|bj) = (jb|ae), (aj|be) = (ja|be)
      !$omp parallel do default(none) shared(chem, eris, no, nv) &
      !$omp    private(a, b, e, j) schedule(static) collapse(2)
      do j = 1, no
         do e = 1, nv
            do b = 1, nv
               do a = 1, nv
                  eris%vvvo(a, b, e, j) = asym2( &
                                          chem%ovvv(spatial(j), spatial(b), spatial(a), spatial(e)), &
                                          same_spin(a, e) .and. same_spin(b, j), &
                                          chem%ovvv(spatial(j), spatial(a), spatial(b), spatial(e)), &
                                          same_spin(a, j) .and. same_spin(b, e))
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine build_cc_eris_fitted

   subroutine spin_orbital_integrals(mo, n_so, w)
      !! The full antisymmetrised spin-orbital tensor <pq||rs>
      !!
      !! Built whole rather than as blocks. The blocks below are slices of it, and
      !! carving them out of one tensor cannot disagree about a sign the way six
      !! separate constructions could.
      real(dp), intent(in) :: mo(:, :, :, :)
      integer, intent(in) :: n_so
      real(dp), allocatable, intent(out) :: w(:, :, :, :)

      integer :: p, q, r, s

      allocate (w(n_so, n_so, n_so, n_so))
      !$omp parallel do default(none) shared(w, mo, n_so) private(p, q, r, s) &
      !$omp    schedule(static) collapse(2)
      do s = 1, n_so
         do r = 1, n_so
            do q = 1, n_so
               do p = 1, n_so
                  w(p, q, r, s) = asym(mo, p, q, r, s)
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine spin_orbital_integrals

   subroutine run_libcint_ccsd(mol, coeff, orbital_energies, n_occ, frozen, &
                               max_iter, energy_tol, want_triples, verbose, &
                               result, error, diis_vectors, aux)
      !! Drive CCSD, and optionally (T), to convergence
      !!
      !! **Memory.** The spin-orbital blocks are sixteen times their spatial
      !! size, and the two n_occ n_vir^3 ones set the ceiling. On the
      !! conventional path the n_act^4 spatial tensor sits on top of them,
      !! because the ladder reads (ae|bf) out of it and there is nowhere else to
      !! get it; on the fitted path there is no such tensor at all -- the ladder
      !! builds what it needs from `b_vv` for a few per cent of its own cost, so
      !! RI-CCSD reaches further than CCSD rather than merely differing from it
      !! by the fitting error. The (T) step holds no triples amplitude either --
      !! it works one occupied triple at a time, n_vir^3 at once.
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
         !! rather than transformed exactly -- RI-CCSD.
         !!
         !! *Every* class, not just the particle-particle ladder, because that is
         !! what PySCF's `dfccsd` does and the comparison is only possible against
         !! the same approximation. Its `ao2mo` returns `_make_df_eris`, which
         !! builds oooo, ovov, ovvo, oovv and ovvv from the fitted tensor and
         !! keeps `vvL` for the ladder -- so "it overrides `_add_vvvv`" describes
         !! where the ladder is contracted, not which integrals are fitted. Fitting
         !! only the ladder would leave us permanently ~1e-5 from the reference,
         !! which is close enough to a convergence threshold to be argued about
         !! rather than diagnosed.

      real(dp), allocatable :: eri(:, :), mo(:, :, :, :), b_vv(:, :)
      type(cc_eris_t) :: eris
      type(cc_chem_t) :: chem
      real(dp), allocatable :: c_act(:, :), eps(:)
      real(dp), allocatable :: t1(:, :), t2(:, :, :, :), t1n(:, :), t2n(:, :, :, :)
      real(dp), allocatable :: flat(:), err_flat(:)
      type(diis_state_t) :: diis
      logical :: extrapolated
      integer :: n_ao, n_mo, n_act, n_o, n_v, no, nv, n_so, iter, diis_size
      integer :: i, a, s
      real(dp) :: e_corr, e_old, de, t_iter

      character(len=MAX_LINE_LENGTH) :: line

      type(timing_report_t) :: clk

      n_ao = size(coeff, 1)
      n_mo = size(coeff, 2)

      if (frozen < 0 .or. frozen >= n_occ) then
         call error%set(ERROR_VALIDATION, "CCSD: the frozen core count must leave at "// &
                        "least one occupied orbital")
         return
      end if
      n_act = n_mo - frozen
      n_o = n_occ - frozen
      n_v = n_mo - n_occ
      if (n_v < 1) then
         call error%set(ERROR_VALIDATION, "CCSD: no virtual orbitals -- the basis is "// &
                        "saturated by the occupied space and there is nothing to excite into")
         return
      end if

      no = 2*n_o
      nv = 2*n_v
      n_so = no + nv

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
         ! Reported because the memory is the interesting constraint here, and
         ! because the number makes the claim checkable rather than asserted. The
         ! comparison worth knowing is the whole spin-orbital tensor, n_so^4, which
         ! is what this used to build and no longer does.
         write (line, "(a,f0.1,a,f0.1,a)") "  integrals: ", integral_megabytes(no, nv), &
            " MB in blocks, against ", real(n_so, dp)**4*8.0_dp/1.0e6_dp, " MB for the full tensor"
         call logger%info(trim(line))
         ! The spatial tensor is the conventional path's alone. Said out loud
         ! because it is the difference between the two paths' ceilings and it
         ! is otherwise invisible -- both report the same block total above.
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
      allocate (c_act(n_ao, n_act))
      c_act = coeff(:, frozen + 1:n_mo)
      if (present(aux)) then
         call fitted_chem_blocks(mol, aux, c_act(:, 1:n_o), c_act(:, n_o + 1:n_act), &
                                 chem, b_vv, error)
         if (error%has_error()) return
      else
         call mol%eris_packed(eri)
         call transform_block(eri, c_act, c_act, c_act, c_act, mo)
         deallocate (eri)
      end if
      deallocate (c_act)

      ! Blocks rather than the whole (2 n_act)^4 tensor. On the conventional path
      ! `mo` is kept afterwards because the ladder assembles its own from it in
      ! batches; on the fitted one the spatial blocks are finished with here and
      ! only `b_vv` outlives them.
      call clk%lap("AO->MO integrals")
      if (present(aux)) then
         call build_cc_eris_fitted(chem, no, nv, eris)
         deallocate (chem%oooo, chem%ooov, chem%oovv, chem%ovov, chem%ovvv)
      else
         call build_cc_eris(mo, no, nv, eris)
      end if
      call clk%lap("CC integral blocks")

      ! Spin-orbital energies. Canonical, so the Fock matrix is this diagonal and
      ! nothing else -- which is what lets every denominator below be a sum of
      ! four of these rather than a matrix element.
      allocate (eps(n_so))
      do s = 1, n_so
         eps(s) = orbital_energies(frozen + spatial(s))
      end do

      ! ---- MP2, as the checkpoint before any amplitude equation --------------
      allocate (t1(nv, no), t2(nv, nv, no, no))
      t1 = 0.0_dp
      call mp2_amplitudes(eris, eps, no, nv, t2, result%e_mp2)
      call clk%lap("MP2 amplitudes")
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
         call ccsd_iteration(mo, b_vv, eris, eps, no, nv, t1, t2, t1n, t2n)
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
         if (verbose) call logger%info("  CCSD iterations done! beginning the perturbative triples (T)")
         call triples_correction(eris, eps, t1, t2, no, nv, result%e_triples)
         call clk%lap("(T) correction")
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
      !! wants it with four; an explicit-shape dummy takes the destination by
      !! sequence association, so the gemm writes straight into the rank-four
      !! array rather than into a temporary that is then reshaped. At `(ov|vv)`,
      !! n_occ n_vir^3, that temporary was the peak.
      real(dp), intent(in) :: bl(:, :), br(:, :)
      integer, intent(in) :: nl, nr
      real(dp), intent(out) :: block(nl, nr)

      call pic_gemm(bl, br, block, transb="T")
   end subroutine df_block_gemm

   subroutine fitted_chem_blocks(mol, aux, c_occ, c_vir, chem, b_vv, error)
      !! The fitted spatial integrals RI-CCSD actually reads
      !!
      !! Every antisymmetrised block the amplitude equations want carries at
      !! least one occupied index -- the sole exception is `<ab||ef>`, which the
      !! particle-particle ladder assembles for itself from `b_vv`. So the five
      !! blocks below are the whole of what has to be materialised, and the
      !! largest of them is `(ov|vv)` at n_occ n_vir^3.
      !!
      !! This is what RI buys, and until now it bought nothing: the previous
      !! version formed the same n_act^4 tensor the conventional path forms, so
      !! the ceiling was identical and only the numbers moved, by the fitting
      !! error. Against that tensor these blocks are smaller by roughly
      !! n_vir/n_occ -- and the one-time n_act^4 gemm that built it, naux n_act^4,
      !! is gone with it.
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
      call build_df_mo_block(mol, aux, c_occ, c_occ, b_oo, error)
      if (error%has_error()) return
      call build_df_mo_block(mol, aux, c_occ, c_vir, b_ov, error)
      if (error%has_error()) return
      call build_df_mo_block(mol, aux, c_vir, c_vir, b_vv_pq, error)
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
      !! This is the free checkpoint the plan asks for: it must reproduce
      !! `run_libcint_mp2` exactly, and that already agrees with PySCF. If it
      !! does not, the antisymmetrisation or the spin-orbital index map is wrong
      !! and no amplitude equation after it can be right.
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
      !! Reported as two numbers because `cc_energy_t` carries them separately and
      !! a lumped correlation energy cannot be taken apart afterwards. The f_ia
      !! t_i^a term of the general expression is absent: the reference is
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

   subroutine ccsd_iteration(mo, b_vv, eris, eps, no, nv, t1, t2, t1n, t2n)
      !! One CCSD amplitude update, through the Stanton et al. intermediates
      real(dp), allocatable, intent(in) :: mo(:, :, :, :)
         !! The spatial tensor, on the conventional path. Read only by the
         !! particle-particle ladder, which assembles <ab||ef> in batches rather
         !! than holding it -- see `particle_ladder`.
      real(dp), allocatable, intent(in) :: b_vv(:, :)
         !! The fitted three-index virtual block, on the RI path. Exactly one of
         !! this and `mo` is allocated; the ladder is the only thing that cares
         !! which.
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
      call particle_ladder(mo, b_vv, eris, t1, tau2, oovv2, no, nv, lad)

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

   subroutine particle_ladder(mo, b_vv, eris, t1, tau2, oovv2, no, nv, lad)
      !! The particle-particle ladder, never holding <ab||ef>
      !!
      !!     lad(ab,ij) += 1/2 sum_ef Wabef(ab,ef) tau(ef,ij)
      !!
      !! Wabef is n_vir^4 and is the largest array coupled cluster asks for -- 17
      !! MB at thirty-eight virtual spin orbitals, and the term that decided
      !! whether CCSD was usable here. It is also read exactly once, by this
      !! contraction, which is what makes batching it possible rather than merely
      !! desirable: a block of `ef` columns is built, contracted, and dropped.
      !!
      !! The batching is over `ef` and every operand stays contiguous, which is
      !! the only reason this is still gemms:
      !!
      !!   * `oovv2(:, ef0:ef1)` is a column range of an (mn, ef) matrix;
      !!   * `tau2` is passed whole for the quadratic term;
      !!   * the accumulation needs tau(ef, ij) over the batch, which is a *row*
      !!     range of tau2 and so not contiguous -- hence `tau2t`, the transpose,
      !!     built once per iteration for 1.4 MB and sliced by column instead.
      !!
      !! Nothing about the arithmetic changes. The energies are bit-identical to
      !! the version that held the whole array.
      !!
      !! **Where <ab||ef> comes from.** Exactly one of `mo` and `b_vv` is
      !! allocated. `mo` is the conventional path's spatial tensor, read at
      !! permuted indices. `b_vv` is the fitted three-index block, and the
      !! spatial (a e | b f) for one pair of spatial (e, f) is then a gemm over
      !! the auxiliary index -- which is why the fitted path never needs a
      !! spatial tensor at all. Those gemms cost naux n_vir^4 an iteration
      !! against the n_occ^2 n_vir^4 of the contraction they feed, so a few per
      !! cent of the ladder buys the whole n_act^4 array.
      real(dp), allocatable, intent(in) :: mo(:, :, :, :)
      real(dp), allocatable, intent(in) :: b_vv(:, :)   !! (P, ab), spatial virtuals
      type(cc_eris_t), intent(in) :: eris
      real(dp), intent(in) :: t1(:, :), tau2(:, :), oovv2(:, :)
      integer, intent(in) :: no, nv
      real(dp), intent(inout) :: lad(:, :)

      !> Columns of (ef) per pass. n_vir^2 by this, so 12 MB at thirty-eight
      !> virtuals against the 17 MB the whole array would take -- and the ratio
      !> improves as n_vir grows, which is the point. Chosen large because the
      !> assembly is cheap and it is the *number of passes* that costs: each one
      !> opens a parallel region, and at 256 the barrier time was larger than all
      !> the arithmetic in the routine.
      real(dp), allocatable :: wblk(:, :), tau2t(:, :), tblk(:, :)
      real(dp), allocatable :: gvv(:, :), gvt(:, :)
      integer :: nv2, no2, ef0, ef1, nb, col, ef, a, b, e, f
      integer :: nocc2, nvs, pa, pb, pe, pf, se, sf, sb, e0, f0
      integer :: pe_have, pf_local, x, y
      logical :: fitted
      real(dp) :: acc

      nv2 = nv*nv
      no2 = no*no
      ! The occupied spin orbitals are the first `no`, and `no` is even, so a
      ! virtual spin orbital `a` carries spin `mod(a,2)` and sits in spatial
      ! orbital `nocc2 + (a+1)/2` -- or, counting virtuals only as `b_vv` does,
      ! in `(a+1)/2`. Hoisting that out of `asym` is what lets the loop below
      ! read plain 2-D slices instead of recomputing four index maps and two spin
      ! tests for every one of the n_vir^4 elements.
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
            !$omp    private(col, ef, a, b, e, f, pa, pb, pe, pf, se, sf, sb, acc) &
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
            !$omp    shared(mo, tblk, wblk, nv, nb, ef0, nocc2) &
            !$omp    private(col, ef, a, b, e, f, pa, pb, pe, pf, se, sf, sb, acc) &
            !$omp    schedule(static)
            do col = 1, nb
               ef = ef0 + col - 1
               e = mod(ef - 1, nv) + 1
               f = (ef - 1)/nv + 1
               se = mod(e, 2)
               sf = mod(f, 2)
               pe = nocc2 + (e + 1)/2
               pf = nocc2 + (f + 1)/2
               do b = 1, nv
                  sb = mod(b, 2)
                  pb = nocc2 + (b + 1)/2
                  do a = 1, nv
                     ! The one place <ab||ef> is needed, straight from the
                     ! spatial tensor. This is `asym` with the spin
                     ! bookkeeping lifted out of the innermost loop -- the
                     ! same two terms, the same signs.
                     pa = nocc2 + (a + 1)/2
                     acc = 0.0_dp
                     if (mod(a, 2) == se .and. sb == sf) acc = acc + mo(pa, pe, pb, pf)
                     if (mod(a, 2) == sf .and. sb == se) acc = acc - mo(pa, pf, pb, pe)
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
      !! The whole of this routine is one gemm. It exists because `eris%ovvv` is
      !! rank four and the product wants it as the (no, nv*nb) matrix its memory
      !! already is: an explicit-shape dummy takes the block by sequence
      !! association, which an assumed-shape one cannot do. Nothing is copied.
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
      !! Two of them, differing only in that half, and mixing them up is a
      !! mistake that still converges -- to the wrong energy by a few
      !! millihartree. Built together so the difference is visible in one place.
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
      !! Non-iterative, so one pass after CCSD converges, and n_occ^3 n_vir^4 --
      !! the wall. No triples amplitude is stored: the two n_vir^3 blocks are
      !! rebuilt for each occupied triple, which is what keeps this inside memory
      !! when a full t3 would be n_o^3 n_v^3.
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
      !! symmetric under permuting them and vanishes when two of them coincide.
      !! Only i > j > k is therefore visited, weighted by six: exactly the same
      !! number, a sixth of the work. The virtual permutations are *not*
      !! restricted -- P(a/bc) is applied inside `triples_block`, which needs the
      !! whole (a,b,c) cube to read at permuted indices, so there is nothing to
      !! skip there.
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
      ! reads out of the full spin-orbital tensor. All of them are small -- the
      ! largest is n_o n_v^3, 4.4 MB at 38 virtuals -- and each is touched
      ! n_occ^3 times, so packing once is the whole optimisation.
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
      !! P(a/bc) likewise, so each numerator is nine signed contributions. The
      !! nine are not nine pieces of work, though, and that is the whole point of
      !! this arrangement: for a fixed occupied permutation the quantity being
      !! permuted over the virtuals,
      !!
      !!     G(a,b,c) = -sum_e t2(a,e,j,k) <ie||bc> - sum_m t2(b,c,i,m) <ma||jk>
      !!
      !! is one object over all (a,b,c), and P(a/bc) then only reads it at
      !! permuted indices. So three occupied permutations give three pairs of
      !! gemms, and the virtual permutations cost nothing but the indexing.
      !!
      !! Written the direct way -- nine calls to a function computing one element
      !! -- this was 118 of 260 seconds of the whole program, because every one of
      !! the nine redid an O(n_vir) sum for every (a,b,c). The equations are
      !! unchanged; only the order of the loops is.
      !!
      !! The scratch arrays come in from the caller rather than being allocated
      !! here: this runs n_occ^3 times, and four allocations per call was itself
      !! measurable.
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
      !! It vanishes exactly at a fixed point of the amplitude equations, which
      !! is what convergence means here, so it is the right quantity to
      !! extrapolate against -- the same argument as the FDS - SDF commutator in
      !! the SCF.
      real(dp), intent(in) :: t1n(:, :), t2n(:, :, :, :), t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: flat(:)

      integer :: n1

      n1 = nv*no
      flat(1:n1) = reshape(t1n - t1, [n1])
      flat(n1 + 1:n1 + nv*nv*no*no) = reshape(t2n - t2, [nv*nv*no*no])
   end subroutine pack_step

end module mqc_libcint_cc
