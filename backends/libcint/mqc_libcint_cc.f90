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
   !! **What this is not.** Every integral block is materialised, including
   !! `<ab||cd>` at n_vir^4 in spin orbitals. That is the wall, and it is
   !! deliberate: correctness first, on systems small enough to check against
   !! PySCF, with the particle-particle ladder as the obvious first thing to fit
   !! or batch when this needs to be fast. See the memory note in
   !! `run_libcint_ccsd`.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_diis, only: diis_state_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_df_mo_block
   use mqc_libcint_mp2, only: transform_block
   implicit none
   private

   public :: cc_result_t
   public :: run_libcint_ccsd
   public :: spin_orbital_integrals   !! Exported so a test can check antisymmetry

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
      !! **Memory.** The spin-orbital tensor is (2 n_act)^4, which is sixteen
      !! times the spatial one: 42 MB at 24 active orbitals, 1.1 GB at 48, and
      !! out of reach not long after. That ceiling is why this is a reference
      !! implementation and is checked on molecules that fit inside it. The (T)
      !! step is the exception and holds no triples amplitude -- it works one
      !! occupied triple at a time, n_vir^3 at once.
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

      real(dp), allocatable :: eri(:, :), mo(:, :, :, :), w(:, :, :, :)
      real(dp), allocatable :: c_act(:, :), eps(:)
      real(dp), allocatable :: t1(:, :), t2(:, :, :, :), t1n(:, :), t2n(:, :, :, :)
      real(dp), allocatable :: flat(:), err_flat(:)
      type(diis_state_t) :: diis
      logical :: extrapolated
      integer :: n_ao, n_mo, n_act, n_o, n_v, no, nv, n_so, iter, diis_size
      integer :: i, a, s
      real(dp) :: e_corr, e_old, de

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
         write (*, "(a,i0,a,i0,a,i0)") "  coupled cluster: ", n_so, " spin orbitals, ", &
            no, " occupied, ", nv, " virtual"
         if (frozen > 0) write (*, "(a,i0,a)") "  frozen core: ", frozen, " spatial orbitals"
      end if

      ! ---- integrals --------------------------------------------------------
      allocate (c_act(n_ao, n_act))
      c_act = coeff(:, frozen + 1:n_mo)
      if (present(aux)) then
         call fitted_mo_tensor(mol, aux, c_act, mo, error)
         if (error%has_error()) return
      else
         call mol%eris_packed(eri)
         call transform_block(eri, c_act, c_act, c_act, c_act, mo)
         deallocate (eri)
      end if
      deallocate (c_act)

      call spin_orbital_integrals(mo, n_so, w)
      deallocate (mo)

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
      call mp2_amplitudes(w, eps, no, nv, t2, result%e_mp2)
      if (verbose) write (*, "(a,f20.12)") "  MP2 (spin orbital) = ", result%e_mp2

      ! ---- CCSD -------------------------------------------------------------
      allocate (t1n(nv, no), t2n(nv, nv, no, no))
      allocate (flat(nv*no + nv*nv*no*no), err_flat(nv*no + nv*nv*no*no))
      call diis%init(diis_size, size(flat), size(err_flat))

      e_old = 0.0_dp
      result%converged = .false.
      do iter = 1, max_iter
         call ccsd_iteration(w, eps, no, nv, t1, t2, t1n, t2n)

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

         call ccsd_energy(w, t1, t2, no, nv, result%e_singles, result%e_doubles)
         e_corr = result%e_singles + result%e_doubles
         de = abs(e_corr - e_old)
         if (verbose) write (*, "(a,i4,a,f20.12,a,es10.3,a,i0)") &
            "  cc iter ", iter, "  E_corr = ", e_corr, "  dE = ", de, &
            "  diis ", diis%count()

         e_old = e_corr
         result%iterations = iter
         if (iter > 1 .and. de < energy_tol) then
            result%converged = .true.
            exit
         end if
      end do
      call diis%destroy()

      if (.not. result%converged) then
         call error%set(ERROR_VALIDATION, "CCSD did not converge")
         return
      end if

      ! ---- (T) --------------------------------------------------------------
      if (want_triples) then
         call triples_correction(w, eps, t1, t2, no, nv, result%e_triples)
         if (verbose) write (*, "(a,f20.12)") "  (T) = ", result%e_triples
      end if

      result%e_correlation = result%e_singles + result%e_doubles + result%e_triples
   end subroutine run_libcint_ccsd

   subroutine fitted_mo_tensor(mol, aux, c_act, mo, error)
      !! The active MO tensor from the fitted three-index one
      !!
      !!     (pq|rs) = sum_P B^P_pq B^P_rs
      !!
      !! One gemm, because `build_df_mo_block` lays B out with the compound index
      !! leading. That is the whole of RI-CCSD's integral step.
      !!
      !! **This does not save memory yet, and is not meant to.** The n_act^4
      !! tensor is still formed, so the ceiling is where it was; what changes is
      !! the number, by the fitting error. Never forming it is the optimisation
      !! that follows, and doing it in this order means the fitting can be
      !! verified against PySCF before the loop structure changes underneath it --
      !! the same order the conventional path was built in.
      type(libcint_molecule_t), intent(in) :: mol, aux
      real(dp), intent(in) :: c_act(:, :)
      real(dp), allocatable, intent(out) :: mo(:, :, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: b(:, :), mo2(:, :)
      integer :: n_act

      n_act = size(c_act, 2)
      call build_df_mo_block(mol, aux, c_act, c_act, b, error)
      if (error%has_error()) return

      allocate (mo2(n_act*n_act, n_act*n_act))
      call pic_gemm(b, b, mo2, transb="T")
      deallocate (b)

      ! Reshaped rather than aliased: `asym` reads the tensor with four indices,
      ! and it is the one function in this module worth not touching.
      allocate (mo(n_act, n_act, n_act, n_act))
      mo = reshape(mo2, [n_act, n_act, n_act, n_act])
      deallocate (mo2)
   end subroutine fitted_mo_tensor

   subroutine mp2_amplitudes(w, eps, no, nv, t2, e_mp2)
      !! First-order doubles, and the MP2 energy they give
      !!
      !!     t_ij^ab = <ij||ab> / D_ij^ab,   E = 1/4 sum <ij||ab> t_ij^ab
      !!
      !! This is the free checkpoint the plan asks for: it must reproduce
      !! `run_libcint_mp2` exactly, and that already agrees with PySCF. If it
      !! does not, the antisymmetrisation or the spin-orbital index map is wrong
      !! and no amplitude equation after it can be right.
      real(dp), intent(in) :: w(:, :, :, :), eps(:)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: t2(:, :, :, :)
      real(dp), intent(out) :: e_mp2

      integer :: i, j, a, b
      real(dp) :: d

      e_mp2 = 0.0_dp
      !$omp parallel do default(none) shared(w, eps, t2, no, nv) &
      !$omp    private(i, j, a, b, d) reduction(+:e_mp2) schedule(static) collapse(2)
      do j = 1, no
         do i = 1, no
            do b = 1, nv
               do a = 1, nv
                  d = eps(i) + eps(j) - eps(no + a) - eps(no + b)
                  t2(a, b, i, j) = w(i, j, no + a, no + b)/d
                  e_mp2 = e_mp2 + 0.25_dp*w(i, j, no + a, no + b)*t2(a, b, i, j)
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine mp2_amplitudes

   subroutine ccsd_energy(w, t1, t2, no, nv, e_singles, e_doubles)
      !! E_CCSD = 1/4 sum <ij||ab> t_ij^ab + 1/2 sum <ij||ab> t_i^a t_j^b
      !!
      !! Reported as two numbers because `cc_energy_t` carries them separately and
      !! a lumped correlation energy cannot be taken apart afterwards. The f_ia
      !! t_i^a term of the general expression is absent: the reference is
      !! canonical, so the occupied-virtual Fock block is zero.
      real(dp), intent(in) :: w(:, :, :, :), t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: e_singles, e_doubles

      integer :: i, j, a, b
      real(dp) :: v

      e_singles = 0.0_dp
      e_doubles = 0.0_dp
      !$omp parallel do default(none) shared(w, t1, t2, no, nv) &
      !$omp    private(i, j, a, b, v) reduction(+:e_singles, e_doubles) &
      !$omp    schedule(static) collapse(2)
      do j = 1, no
         do i = 1, no
            do b = 1, nv
               do a = 1, nv
                  v = w(i, j, no + a, no + b)
                  e_doubles = e_doubles + 0.25_dp*v*t2(a, b, i, j)
                  e_singles = e_singles + 0.5_dp*v*t1(a, i)*t1(b, j)
               end do
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine ccsd_energy

   subroutine ccsd_iteration(w, eps, no, nv, t1, t2, t1n, t2n)
      !! One CCSD amplitude update, through the Stanton et al. intermediates
      real(dp), intent(in) :: w(:, :, :, :), eps(:)
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
      real(dp), allocatable :: wmnij(:, :), wabef(:, :), wmbej(:, :, :, :)
      real(dp), allocatable :: oovv2(:, :), lad(:, :)
      real(dp), allocatable :: gae(:, :), gmi(:, :)
      integer :: i, j, m, n, a, b, e, f, ab, ij, mn, ef
      real(dp) :: acc, d, term

      allocate (tau2(nv*nv, no*no), taut(nv, nv, no, no))
      call build_tau(t1, t2, no, nv, tau2, taut)

      ! <mn||ef> as a matrix over the two pairs, which is the shape three of the
      ! gemms below want it in.
      allocate (oovv2(no*no, nv*nv))
      !$omp parallel do default(none) shared(w, oovv2, no, nv) &
      !$omp    private(m, n, e, f, mn, ef) schedule(static) collapse(2)
      do f = 1, nv
         do e = 1, nv
            ef = (f - 1)*nv + e
            do n = 1, no
               do m = 1, no
                  oovv2((n - 1)*no + m, ef) = w(m, n, no + e, no + f)
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! ---- one-body intermediates (Stanton 3-5) -----------------------------
      allocate (fae(nv, nv), fmi(no, no), fme(no, nv))

      !$omp parallel do default(none) shared(w, t1, taut, fae, no, nv) &
      !$omp    private(a, e, m, f, n, acc) schedule(static)
      do e = 1, nv
         do a = 1, nv
            acc = 0.0_dp
            do m = 1, no
               do f = 1, nv
                  acc = acc + t1(f, m)*w(m, no + a, no + f, no + e)
               end do
            end do
            do n = 1, no
               do m = 1, no
                  do f = 1, nv
                     acc = acc - 0.5_dp*taut(a, f, m, n)*w(m, n, no + e, no + f)
                  end do
               end do
            end do
            fae(a, e) = acc
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(none) shared(w, t1, taut, fmi, no, nv) &
      !$omp    private(m, i, e, n, f, acc) schedule(static)
      do i = 1, no
         do m = 1, no
            acc = 0.0_dp
            do n = 1, no
               do e = 1, nv
                  acc = acc + t1(e, n)*w(m, n, i, no + e)
               end do
            end do
            do n = 1, no
               do f = 1, nv
                  do e = 1, nv
                     acc = acc + 0.5_dp*taut(e, f, i, n)*w(m, n, no + e, no + f)
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
                  acc = acc + t1(f, n)*w(m, n, no + e, no + f)
               end do
            end do
            fme(m, e) = acc
         end do
      end do

      ! ---- two-body intermediates (Stanton 6-8) -----------------------------
      allocate (wmnij(no*no, no*no), wabef(nv*nv, nv*nv), wmbej(no, nv, nv, no))

      ! The quadratic term of each ladder intermediate is a gemm, so it goes in
      ! first and the linear terms are added on top of it.
      call pic_gemm(oovv2, tau2, wmnij, alpha=0.25_dp, beta=0.0_dp)
      !$omp parallel do default(none) shared(w, t1, wmnij, no, nv) &
      !$omp    private(m, n, i, j, e, acc, ij) schedule(static) collapse(2)
      do j = 1, no
         do i = 1, no
            ij = (j - 1)*no + i
            do n = 1, no
               do m = 1, no
                  acc = w(m, n, i, j)
                  do e = 1, nv
                     acc = acc + t1(e, j)*w(m, n, i, no + e) - t1(e, i)*w(m, n, j, no + e)
                  end do
                  wmnij((n - 1)*no + m, ij) = wmnij((n - 1)*no + m, ij) + acc
               end do
            end do
         end do
      end do
      !$omp end parallel do

      call pic_gemm(tau2, oovv2, wabef, alpha=0.25_dp, beta=0.0_dp)
      !$omp parallel do default(none) shared(w, t1, wabef, no, nv) &
      !$omp    private(a, b, e, f, m, acc, ef) schedule(static) collapse(2)
      do f = 1, nv
         do e = 1, nv
            ef = (f - 1)*nv + e
            do b = 1, nv
               do a = 1, nv
                  acc = w(no + a, no + b, no + e, no + f)
                  do m = 1, no
                     ! <am||ef> = -<ma||ef>, and likewise for b. Writing it this
                     ! way keeps every integral read out of the same block.
                     acc = acc + t1(b, m)*w(m, no + a, no + e, no + f) &
                           - t1(a, m)*w(m, no + b, no + e, no + f)
                  end do
                  wabef((b - 1)*nv + a, ef) = wabef((b - 1)*nv + a, ef) + acc
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! ---- both ladders, one gemm each --------------------------------------
      !
      ! The two terms that dominate CCSD: n_o^4 n_v^2 for the holes and
      ! n_v^4 n_o^2 for the particles. As scalar loops they also had the worst
      ! access pattern in the module -- for fixed (a,b) the particle ladder walks
      ! `wabef` down strides of n_v^2 and n_v^3 -- and moving them here took the
      ! T2 equation from 90 seconds of CPU time to a gemm.
      allocate (lad(nv*nv, no*no))
      call pic_gemm(tau2, wmnij, lad, alpha=0.5_dp, beta=0.0_dp)
      call pic_gemm(wabef, tau2, lad, alpha=0.5_dp, beta=1.0_dp)

      !$omp parallel do default(none) shared(w, t1, t2, wmbej, no, nv) &
      !$omp    private(m, b, e, j, f, n, acc) schedule(static) collapse(2)
      do j = 1, no
         do e = 1, nv
            do b = 1, nv
               do m = 1, no
                  acc = w(m, no + b, no + e, j)
                  do f = 1, nv
                     acc = acc + t1(f, j)*w(m, no + b, no + e, no + f)
                  end do
                  do n = 1, no
                     ! <mn||ej> = <mn||je> with the last pair swapped, hence the
                     ! sign: -(-ooov) = +.
                     acc = acc + t1(b, n)*w(m, n, j, no + e)
                  end do
                  do n = 1, no
                     do f = 1, nv
                        acc = acc - (0.5_dp*t2(f, b, j, n) + t1(f, j)*t1(b, n)) &
                              *w(m, n, no + e, no + f)
                     end do
                  end do
                  wmbej(m, b, e, j) = acc
               end do
            end do
         end do
      end do
      !$omp end parallel do

      ! ---- T1 (Stanton 1) ---------------------------------------------------
      !$omp parallel do default(none) shared(w, eps, t1, t2, fae, fmi, fme, t1n, no, nv) &
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
                  acc = acc - t1(f, n)*w(n, no + a, i, no + f)
               end do
            end do
            do m = 1, no
               do f = 1, nv
                  do e = 1, nv
                     acc = acc - 0.5_dp*t2(e, f, i, m)*w(m, no + a, no + e, no + f)
                  end do
               end do
            end do
            do n = 1, no
               do m = 1, no
                  do e = 1, nv
                     ! <nm||ei> = <mn||ie>, two pair swaps.
                     acc = acc - 0.5_dp*t2(a, e, m, n)*w(m, n, i, no + e)
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
      !$omp    shared(w, eps, t1, t2, lad, wmbej, gae, gmi, t2n, no, nv) &
      !$omp    private(a, b, i, j, e, f, m, n, acc, d, term) schedule(static) collapse(2)
      do j = 1, no
         do i = 1, no
            do b = 1, nv
               do a = 1, nv
                  acc = w(i, j, no + a, no + b)

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

                  ! Rings: P(ij)P(ab) over the four index pairings.
                  do m = 1, no
                     do e = 1, nv
                        term = t2(a, e, i, m)*wmbej(m, b, e, j) &
                               - t1(e, i)*t1(a, m)*w(m, no + b, no + e, j)
                        acc = acc + term
                        term = t2(a, e, j, m)*wmbej(m, b, e, i) &
                               - t1(e, j)*t1(a, m)*w(m, no + b, no + e, i)
                        acc = acc - term
                        term = t2(b, e, i, m)*wmbej(m, a, e, j) &
                               - t1(e, i)*t1(b, m)*w(m, no + a, no + e, j)
                        acc = acc - term
                        term = t2(b, e, j, m)*wmbej(m, a, e, i) &
                               - t1(e, j)*t1(b, m)*w(m, no + a, no + e, i)
                        acc = acc + term
                     end do
                  end do

                  ! P(ij) sum_e t1(e,i) <ab||ej>
                  do e = 1, nv
                     acc = acc + t1(e, i)*w(no + a, no + b, no + e, j) &
                           - t1(e, j)*w(no + a, no + b, no + e, i)
                  end do

                  ! -P(ab) sum_m t1(a,m) <mb||ij>
                  do m = 1, no
                     acc = acc - t1(a, m)*w(m, no + b, i, j) + t1(b, m)*w(m, no + a, i, j)
                  end do

                  d = eps(i) + eps(j) - eps(no + a) - eps(no + b)
                  t2n(a, b, i, j) = acc/d
               end do
            end do
         end do
      end do
      !$omp end parallel do

      deallocate (tau2, taut, fae, fmi, fme, wmnij, wabef, wmbej, gae, gmi)
      deallocate (oovv2, lad)
   end subroutine ccsd_iteration

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

   subroutine triples_correction(w, eps, t1, t2, no, nv, e_triples)
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
      !! The 1/36 is the price of summing every permutation of both triples
      !! rather than restricting to i<j<k and a<b<c. Restricting is the obvious
      !! optimisation and is deliberately not done here -- the unrestricted form
      !! is what the equations say, and this is the pass that has to be right
      !! before anything is made fast.
      real(dp), intent(in) :: w(:, :, :, :), eps(:), t1(:, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), intent(out) :: e_triples

      real(dp), allocatable :: t3c(:, :, :), t3d(:, :, :)
      real(dp), allocatable :: ovvv_p(:, :, :), t2bc(:, :, :)
      real(dp), allocatable :: ovoo_p(:, :, :, :), oovv_p(:, :, :)
      real(dp), allocatable :: acc(:, :), bcc(:, :), gg(:, :, :), dd(:, :, :)
      integer :: i, j, k, a, b, c
      real(dp) :: d, e_local

      e_triples = 0.0_dp

      ! Four integral blocks and one amplitude block, repacked so the
      ! contractions below are gemms over contiguous memory rather than strided
      ! reads out of the full spin-orbital tensor. All of them are small -- the
      ! largest is n_o n_v^3, 4.4 MB at 38 virtuals -- and each is touched
      ! n_occ^3 times, so packing once is the whole optimisation.
      call pack_triples_blocks(w, t2, no, nv, ovvv_p, t2bc, ovoo_p, oovv_p)

      !$omp parallel default(none) &
      !$omp    shared(eps, t1, t2, no, nv, ovvv_p, t2bc, ovoo_p, oovv_p) &
      !$omp    private(i, j, k, a, b, c, d, t3c, t3d, e_local, acc, bcc, gg, dd) &
      !$omp    reduction(+:e_triples)
      allocate (t3c(nv, nv, nv), t3d(nv, nv, nv))
      allocate (acc(nv, nv*nv), bcc(nv*nv, nv), gg(nv, nv, nv), dd(nv, nv, nv))
      !$omp do schedule(dynamic) collapse(2)
      do i = 1, no
         do j = 1, no
            do k = 1, no
               call triples_block(t1, t2, no, nv, i, j, k, ovvv_p, t2bc, ovoo_p, &
                                  oovv_p, acc, bcc, gg, dd, t3c, t3d)
               e_local = 0.0_dp
               do c = 1, nv
                  do b = 1, nv
                     do a = 1, nv
                        d = eps(i) + eps(j) + eps(k) &
                            - eps(no + a) - eps(no + b) - eps(no + c)
                        ! t3c and t3d arrive as numerators; dividing here rather
                        ! than twice inside the block keeps the block additive.
                        e_local = e_local + t3c(a, b, c)*(t3c(a, b, c) + t3d(a, b, c))/d
                     end do
                  end do
               end do
               e_triples = e_triples + e_local/36.0_dp
            end do
         end do
      end do
      !$omp end do
      deallocate (t3c, t3d, acc, bcc, gg, dd)
      !$omp end parallel

      deallocate (ovvv_p, t2bc, ovoo_p, oovv_p)
   end subroutine triples_correction

   subroutine pack_triples_blocks(w, t2, no, nv, ovvv_p, t2bc, ovoo_p, oovv_p)
      !! Repack what (T) contracts, so each contraction is a gemm
      !!
      !! The compound index is always `bc = (c-1) n_v + b`, which puts it in the
      !! position a gemm wants: leading for the amplitude block it contracts over,
      !! trailing for the integral block it indexes.
      real(dp), intent(in) :: w(:, :, :, :), t2(:, :, :, :)
      integer, intent(in) :: no, nv
      real(dp), allocatable, intent(out) :: ovvv_p(:, :, :)    !! (e, bc, i) = <ie||bc>
      real(dp), allocatable, intent(out) :: t2bc(:, :, :)      !! (bc, m, i) = t2(b,c,i,m)
      real(dp), allocatable, intent(out) :: ovoo_p(:, :, :, :)  !! (m, a, j, k) = <ma||jk>
      real(dp), allocatable, intent(out) :: oovv_p(:, :, :)    !! (bc, j, k) = <jk||bc>

      integer :: i, j, k, m, a, b, c, bc

      allocate (ovvv_p(nv, nv*nv, no), t2bc(nv*nv, no, no))
      allocate (ovoo_p(no, nv, no, no), oovv_p(nv*nv, no, no))

      !$omp parallel do default(none) shared(w, ovvv_p, no, nv) &
      !$omp    private(i, b, c, bc) schedule(static)
      do i = 1, no
         do c = 1, nv
            do b = 1, nv
               bc = (c - 1)*nv + b
               ovvv_p(:, bc, i) = w(i, no + 1:no + nv, no + b, no + c)
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

      !$omp parallel do default(none) shared(w, ovoo_p, oovv_p, no, nv) &
      !$omp    private(j, k, m, a, b, c, bc) schedule(static) collapse(2)
      do k = 1, no
         do j = 1, no
            do a = 1, nv
               do m = 1, no
                  ovoo_p(m, a, j, k) = w(m, no + a, j, k)
               end do
            end do
            do c = 1, nv
               do b = 1, nv
                  bc = (c - 1)*nv + b
                  oovv_p(bc, j, k) = w(j, k, no + b, no + c)
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
