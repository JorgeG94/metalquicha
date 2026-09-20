!! The response operator's mean field, in one place for every route that applies it
module mqc_czt_response_product
   !! `(A+B)u` and `(A-B)u`, batched, over whichever integral source is in hand.
   !!
   !! Static coupled-perturbed solves, frequency-dependent ones, and the
   !! analytic Hessian's fixed point all apply the same operator to a trial
   !! rotation; they differ in what the trial vector is indexed by and in what
   !! is done with the image, not in the physics between. Written out three
   !! times, they drifted: the static route grew the exchange-correlation
   !! kernel and no long-range exchange, the frequency-dependent one grew
   !! neither, and the Hessian grew both. This module is the one copy, and the
   !! three callers are the wrappers that remain.
   !!
   !! Two levels, because the Hessian's trial vector is `(n_mo, n_occ)` and
   !! carries occupied rows while everything else is `(n_vir, n_occ)`:
   !!
   !! * `response_mean_field` takes AO response densities and returns `G(D')`
   !!   -- Coulomb, exact exchange at whatever fraction the reference kept, the
   !!   attenuated second pass of a range-separated functional, the
   !!   semilocal kernel and the non-local one.
   !! * `response_product` is that plus the transforms either side of it and
   !!   the orbital-energy diagonal, which is the operator itself.
   !!
   !! **What `minus` changes.** `A - B` acts on the *antisymmetrised* response
   !! density. The Coulomb term vanishes there identically -- the integral is
   !! symmetric in its ket pair and the density is not -- and so does the
   !! exchange-correlation kernel, whose response is to a density change that
   !! an antisymmetric matrix does not make. So `minus` is exchange only, both
   !! passes of it, and the grid is never touched.
   !!
   !! **What `triplet` changes.** The Coulomb term goes, and the semilocal
   !! kernel becomes the spin difference `(f_aa - f_ab)/2`. Exchange does not
   !! move -- the same integrals at the same coefficients, both ranges of them
   !! -- because a triplet transition density is alpha minus beta and exchange
   !! is same-spin, while Coulomb sees only the sum, which a triplet leaves
   !! unchanged. VV10 goes for that reason too, and is dropped rather than
   !! computed. `triplet` with `minus` is `(A - B)`, one operator for both
   !! spins.
   use pic_types, only: dp, int64
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_xc, only: xc_context_t, xc_kernel_apply_many, vv10_kernel_apply, &
                         xc_kernel_cache_t, xc_kernel_apply_uks_many, &
                         xc_kernel_cache_uks_t
   use mqc_czt_direct, only: build_fock, build_fock_direct_many, &
                             build_fock_direct_uhf_many, direct_stats_t
   implicit none
   private

   public :: response_mean_field
   public :: response_product
   public :: response_mean_field_df
   public :: response_mean_field_uhf
   public :: response_product_uhf

contains

   subroutine response_mean_field(mol, dens, zero_h, g, error, minus, direct, eri, &
                                  bounds, k_scale, xc, reference, rs_k_lr, rs_omega, &
                                  cache, triplet, density_screen, screen_floor, stats)
      !! `G(D')` for a batch of response densities, over one pass of the integrals
      !!
      !! The Coulomb and exchange terms come from whichever source the caller
      !! has -- the direct build or a stored tensor -- and the kernel terms go
      !! on top of either, since the choice of integral source is orthogonal to
      !! the functional.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(inout) :: dens(:, :, :)
         !! `(n_ao, n_ao, n_set)`, already symmetrised or antisymmetrised.
         !! **Consumed**: with `screen_floor` given it is rescaled in place, and
         !! no caller may read it back afterwards.
      real(dp), intent(in) :: zero_h(:, :)
         !! Added to every set, so a zero matrix returns `G` alone -- which is
         !! what every response caller wants and passes.
      real(dp), allocatable, intent(out) :: g(:, :, :)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: minus
         !! The densities are antisymmetric and this is the `A - B` half:
         !! exchange only, no Coulomb pass and no grid pass. Off by default.
      logical, intent(in), optional :: direct
         !! Recompute the integrals rather than read them from `eri`. On by
         !! default; `eri` is read when it is off.
      real(dp), intent(in), optional :: eri(:, :, :, :)
      real(dp), intent(in), optional :: bounds(:, :)   !! From `schwarz_bounds`
      real(dp), intent(in), optional :: k_scale
         !! The exact-exchange fraction the *reference* kept: one for
         !! Hartree-Fock, zero for a pure functional, the mixing fraction for a
         !! hybrid. Building full exchange for a pure functional makes the
         !! operator indefinite and the solver then reports a saddle point that
         !! is not there. Absent, it is `xc%exx_fraction` when an `xc` context
         !! was given and one otherwise, so a caller holding the context need
         !! not restate what the context already says.
      type(xc_context_t), intent(inout), optional :: xc
         !! Present, the exchange-correlation kernel is added. Needs
         !! `reference`. Ignored when `minus`, where it does not contribute.
      real(dp), intent(in), optional :: reference(:, :)
         !! The converged reference density the kernel is evaluated at, which
         !! is not the trial density being contracted against it.
      real(dp), intent(in), optional :: rs_k_lr
      real(dp), intent(in), optional :: rs_omega
         !! A range-separated functional's second exchange term: `rs_k_lr` of
         !! the exchange built against `erf(rs_omega r)/r`, with no Coulomb,
         !! since the full-range pass already supplied it. Zero or absent
         !! `rs_omega` is the ordinary case and costs nothing. Absent, both
         !! come from `xc` when it is given and says it is range separated --
         !! which is what makes a CAM-B3LYP response carry its long-range
         !! exchange without every caller in the chain knowing about it.
      type(xc_kernel_cache_t), intent(in), optional :: cache
         !! The reference's kernel coefficients over the whole grid, from
         !! `xc_kernel_cache_fill`. Forwarded to the kernel contraction, which
         !! then skips the libxc pass it would otherwise make on every call --
         !! and a response solve makes hundreds. A cache that was never filled
         !! is ignored, so a caller holding one as a plain component may pass
         !! it unconditionally.
      logical, intent(in), optional :: triplet
         !! The **triplet** mean field: no Coulomb term at all, and the
         !! exchange-correlation kernel taken as `(f_aa - f_ab)/2`. Exchange,
         !! of both ranges, is untouched -- the two-electron part of `A_T` is
         !! `-c_x (ab|ij)`, the same integral at the same coefficient as the
         !! singlet's, with only `J` gone. Off by default, and immaterial when
         !! `minus`, whose product is exchange only and the same for the two
         !! spins.
      logical, intent(in), optional :: density_screen
         !! Weight the Schwarz bound by the largest trial-density element a
         !! quartet touches. Off by default, which keeps a batch bit-for-bit
         !! what the same densities give one at a time.
      real(dp), intent(in), optional :: screen_floor
         !! **A tiny trial density is applied at unit scale and its image
         !! scaled back.** The density screen is absolute and a conjugate
         !! gradient direction shrinks as the solve converges; once its largest
         !! element is near the screening tolerance the surviving quartets are
         !! an arbitrary subset of the operator and the image is noise, which
         !! the solver reads as an operator that is not positive definite.
         !! Everything applied here is linear in the density, so scaling a set
         !! up to this floor, screening at the same tolerance and scaling the
         !! image back is the same operator with a screen relative to the
         !! direction. Absent, or with the screen off, nothing is rescaled.
      type(direct_stats_t), intent(out), optional :: stats
         !! Quartets computed and skipped, summed over every integral pass made.

      real(dp), allocatable :: g_lr(:, :, :), vnl(:, :, :), scale(:)
      type(direct_stats_t) :: pass
      real(dp) :: kf, k_lr, omega, dmax
      integer :: n_ao, n_set, p
      logical :: anti, is_direct, screen, use_cache, spin_flip

      if (error%has_error()) return

      n_ao = size(dens, 1)
      n_set = size(dens, 3)
      anti = .false.
      if (present(minus)) anti = minus
      is_direct = .true.
      if (present(direct)) is_direct = direct
      screen = .false.
      if (present(density_screen)) screen = density_screen
      use_cache = .false.
      if (present(cache)) use_cache = cache%filled
      spin_flip = .false.
      if (present(triplet)) spin_flip = triplet
      ! The context first, an explicit coefficient over it. The analytic
      ! Hessian carries these as scalars on its operator and passes them; the
      ! coupled-perturbed routes pass only the context.
      kf = 1.0_dp
      k_lr = 0.0_dp
      omega = 0.0_dp
      if (present(xc)) then
         kf = xc%exx_fraction
         if (xc%range_separated) then
            k_lr = xc%rs_k_lr
            omega = xc%rs_omega
         end if
      end if
      if (present(k_scale)) kf = k_scale
      if (present(rs_k_lr)) k_lr = rs_k_lr
      if (present(rs_omega)) omega = rs_omega

      if (present(xc) .and. .not. present(reference)) then
         call error%set(ERROR_VALIDATION, "the response operator was given an "// &
                        "exchange-correlation context but no reference density to "// &
                        "evaluate its kernel at")
         return
      end if
      if (is_direct .and. .not. present(bounds)) then
         call error%set(ERROR_VALIDATION, "the response mean field was asked for an "// &
                        "integral-direct build without the Schwarz bounds it screens on")
         return
      end if
      if (.not. is_direct .and. .not. present(eri)) then
         call error%set(ERROR_VALIDATION, "the response mean field was asked to read a "// &
                        "stored two-electron tensor that was not given")
         return
      end if
      if (omega > 0.0_dp .and. .not. is_direct) then
         ! The stored tensor is full-range and there is no attenuated twin of
         ! it, so the long-range term could only be dropped -- which is the bug
         ! this module exists to remove, not one to keep quietly.
         call error%set(ERROR_VALIDATION, "a range-separated functional's response "// &
                        "needs the integral-direct build: the stored two-electron "// &
                        "tensor is full-range and carries no attenuated companion")
         return
      end if

      allocate (scale(n_set))
      scale = 1.0_dp
      if (present(screen_floor) .and. screen) then
         if (screen_floor > 0.0_dp) then
            do p = 1, n_set
               dmax = maxval(abs(dens(:, :, p)))
               if (dmax > 0.0_dp .and. dmax < screen_floor) scale(p) = screen_floor/dmax
               if (scale(p) /= 1.0_dp) dens(:, :, p) = scale(p)*dens(:, :, p)
            end do
         end if
      end if

      if (present(stats)) then
         stats%quartets_total = 0_int64
         stats%quartets_computed = 0_int64
         stats%quartets_screened = 0_int64
      end if

      if (is_direct) then
         ! Announced antisymmetric, the fast build's folded accumulation is
         ! exact: it drops the Coulomb term, which vanishes, and antisymmetrises
         ! instead of symmetrising. `build_fock_direct_nosym` writes the same
         ! permutations out at several times the cost and is not needed here.
         !
         ! `j_scale` is the whole two-electron difference a triplet makes: the
         ! Coulomb term is the response to a change in the *total* density, and
         ! the alpha and beta halves of a triplet transition density cancel
         ! there exactly.
         call build_fock_direct_many(mol, zero_h, dens, bounds, g, pass, error, &
                                     k_scale=kf, j_scale=j_fraction(spin_flip), &
                                     antisymmetric=anti, &
                                     density_screen=screen)
         if (error%has_error()) return
         call add_stats(stats, pass)
         if (omega > 0.0_dp) then
            call build_fock_direct_many(mol, zero_h, dens, bounds, g_lr, pass, error, &
                                        k_scale=k_lr, j_scale=0.0_dp, omega=omega, &
                                        antisymmetric=anti, density_screen=screen)
            if (error%has_error()) return
            g = g + g_lr
            call add_stats(stats, pass)
            deallocate (g_lr)
         end if
      else
         ! The plain four-index contraction needs no announcement: nothing in
         ! it assumes a symmetry, and the Coulomb term vanishes of its own
         ! accord on an antisymmetric density.
         allocate (g(n_ao, n_ao, n_set))
         do p = 1, n_set
            call build_fock(zero_h, eri, dens(:, :, p), g(:, :, p), k_scale=kf, &
                            j_scale=j_fraction(spin_flip))
         end do
      end if

      ! The kernel, for a Kohn-Sham reference, on top of whatever built `g`.
      ! Leaving it out does not fail; it converges to the wrong orbital
      ! response. One grid pass serves the whole batch.
      if (present(xc) .and. .not. anti) then
         if (use_cache) then
            call xc_kernel_apply_many(xc, mol, reference, dens, g, error, cache=cache, &
                                      triplet=spin_flip)
         else
            call xc_kernel_apply_many(xc, mol, reference, dens, g, error, &
                                      triplet=spin_flip)
         end if
         if (error%has_error()) return
         ! The non-local kernel, once for the batch rather than per set:
         ! `vv10_kernel_apply`'s pair sweep is O(npts^2) whether it carries one
         ! trial density or a dozen. It accumulates, hence the zeroed buffer.
         !
         ! Not for a triplet. VV10 is a functional of the total density alone,
         ! and a triplet transition density changes that by nothing, so its
         ! response is identically zero rather than merely small; applying it
         ! would contract a spin density against a kernel with no spin channel.
         if ((xc%nlc_b /= 0.0_dp .or. xc%nlc_c /= 0.0_dp) .and. .not. spin_flip) then
            allocate (vnl(n_ao, n_ao, n_set))
            vnl = 0.0_dp
            call vv10_kernel_apply(xc, mol, reference, dens, vnl, error)
            if (error%has_error()) return
            g = g + vnl
            deallocate (vnl)
         end if
      end if

      do p = 1, n_set
         if (scale(p) /= 1.0_dp) g(:, :, p) = g(:, :, p)/scale(p)
      end do

      deallocate (scale)
   end subroutine response_mean_field

   pure function j_fraction(triplet) result(jf)
      !! One for the ordinary response, zero for a triplet
      !!
      !! Written out because it is the single arithmetic difference between the
      !! two two-electron builds, and a bare `0.0` at the call site would say
      !! nothing about which spin it belonged to.
      logical, intent(in) :: triplet
      real(dp) :: jf

      jf = 1.0_dp
      if (triplet) jf = 0.0_dp
   end function j_fraction

   subroutine add_stats(total, pass)
      !! Sum one integral pass's quartet counts into the running total
      type(direct_stats_t), intent(inout), optional :: total
      type(direct_stats_t), intent(in) :: pass

      if (.not. present(total)) return
      total%quartets_total = total%quartets_total + pass%quartets_total
      total%quartets_computed = total%quartets_computed + pass%quartets_computed
      total%quartets_screened = total%quartets_screened + pass%quartets_screened
   end subroutine add_stats

   subroutine response_product(mol, c_occ, c_vir, gaps, zero_h, u, idx, nact, minus, &
                               au, error, direct, eri, bounds, k_scale, xc, reference, &
                               rs_k_lr, rs_omega, cache, triplet, bmat, density_screen, &
                               t_dens, t_fock, t_back)
      !! `(A+B)u` or `(A-B)u` for many trial rotations in one integral pass
      !!
      !! Writing the trial rotation as a density,
      !!
      !!     Dt = C_vir u C_occ^T  +/-  transpose
      !!
      !! and contracting it as an SCF contracts a density gives half the
      !! bracket, so the two-electron part of the image is `2 C_vir^T G(Dt)
      !! C_occ` and the whole operator is
      !!
      !!     (A +/- B) u = (eps_a - eps_i) u + 2 C_vir^T G(Dt) C_occ
      !!
      !! `idx(1:nact)` selects which of the trial vectors are still wanted, so
      !! a converged system stops costing anything rather than riding along,
      !! and the images of the others are left untouched.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: c_occ(:, :)   !! (n_ao, n_occ)
      real(dp), intent(in) :: c_vir(:, :)   !! (n_ao, n_vir)
      real(dp), intent(in) :: gaps(:, :)
         !! `eps_a - eps_i`, (n_vir, n_occ), the operator's diagonal
      real(dp), intent(in) :: zero_h(:, :)
         !! Zero, so the two-electron build returns `G` alone
      real(dp), intent(in) :: u(:, :, :)    !! (n_vir, n_occ, n_sys)
      integer, intent(in) :: idx(:)         !! Which systems are active
      integer, intent(in) :: nact           !! How many of `idx` to read
      logical, intent(in) :: minus          !! `A - B` rather than `A + B`
      real(dp), intent(inout) :: au(:, :, :)
         !! (n_vir, n_occ, n_sys); only the active systems are written
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: direct
      real(dp), intent(in), optional :: eri(:, :, :, :)
      real(dp), intent(in), optional :: bounds(:, :)
      real(dp), intent(in), optional :: k_scale
      type(xc_context_t), intent(inout), optional :: xc
      real(dp), intent(in), optional :: reference(:, :)
      real(dp), intent(in), optional :: rs_k_lr, rs_omega
      type(xc_kernel_cache_t), intent(in), optional :: cache
         !! The reference's kernel coefficients, filled once and reused. See
         !! `response_mean_field`, which is where it ends up on the exact
         !! route; the fitted one below reads it too.
      logical, intent(in), optional :: triplet
         !! The triplet operator: no Coulomb term and the spin-difference
         !! kernel. Refused on the fitted route, which has no Coulomb-free
         !! build. See `response_mean_field`.
      real(dp), intent(in), optional :: bmat(:, :)
         !! The fitted tensor `B(mu nu, P)`, in place of any four-index
         !! integrals. Not a storage choice: it makes the operator the fitted
         !! reference's own, so give it when and only when the reference was
         !! itself fitted.
      logical, intent(in), optional :: density_screen
      real(dp), intent(inout), optional :: t_dens, t_fock, t_back
         !! Seconds accumulated in the three stages, for a caller that reports
         !! where a matrix-free solve went.

      real(dp), allocatable :: dens(:, :, :), g(:, :, :), half(:, :, :), work(:, :)
      real(dp) :: t0, t1, kf
      integer :: n_ao, n_occ, m, j
      logical :: use_cache, spin_flip

      if (error%has_error()) return
      if (nact <= 0) return

      use_cache = .false.
      if (present(cache)) use_cache = cache%filled
      spin_flip = .false.
      if (present(triplet)) spin_flip = triplet

      ! Only the fitted branch needs this here; the others let
      ! `response_mean_field` resolve it the same way.
      kf = 1.0_dp
      if (present(xc)) kf = xc%exx_fraction
      if (present(k_scale)) kf = k_scale

      n_ao = size(c_occ, 1)
      n_occ = size(c_occ, 2)
      allocate (dens(n_ao, n_ao, nact), half(n_ao, n_occ, nact), work(n_ao, n_occ))
      call cpu_time(t0)

      do m = 1, nact
         j = idx(m)
         call pic_gemm(c_vir, u(:, :, j), half(:, :, m))
         call pic_gemm(half(:, :, m), c_occ, dens(:, :, m), transb="T")
         if (minus) then
            dens(:, :, m) = dens(:, :, m) - transpose(dens(:, :, m))
         else
            dens(:, :, m) = dens(:, :, m) + transpose(dens(:, :, m))
         end if
      end do

      call cpu_time(t1)
      if (present(t_dens)) t_dens = t_dens + (t1 - t0)
      t0 = t1

      if (present(bmat)) then
         if (spin_flip) then
            call error%set(ERROR_VALIDATION, "a triplet response cannot be applied "// &
                           "through the fitted tensor: that build assembles its "// &
                           "Coulomb term unconditionally, and a triplet has none")
            return
         end if
         if (present(xc)) then
            if (xc%range_separated) then
               call error%set(ERROR_VALIDATION, "a range-separated functional's "// &
                              "response cannot be applied through the fitted tensor: "// &
                              "there is no attenuated fitting to build its long-range "// &
                              "exchange from")
               return
            end if
         end if
         ! `half` is C_vir u, which is exactly the factor the fitted build
         ! wants: it never assembles Dt at all, and only reads it for Coulomb.
         allocate (g(n_ao, n_ao, nact))
         do m = 1, nact
            call response_mean_field_df(bmat, half(:, :, m), c_occ, dens(:, :, m), &
                                        g(:, :, m), minus, k_scale=kf)
         end do
         if (present(xc) .and. .not. minus) then
            if (.not. present(reference)) then
               call error%set(ERROR_VALIDATION, "the response operator was given an "// &
                              "exchange-correlation context but no reference density "// &
                              "to evaluate its kernel at")
               return
            end if
            ! TODO(mqc): the non-local VV10 kernel is not applied on the
            ! fitted route, only on the two exact ones, so a fitted reference
            ! with VV10 gets a response missing that term. It was missing
            ! before this routine existed too.
            if (use_cache) then
               call xc_kernel_apply_many(xc, mol, reference, dens, g, error, cache=cache)
            else
               call xc_kernel_apply_many(xc, mol, reference, dens, g, error)
            end if
            if (error%has_error()) return
         end if
      else
         call response_mean_field(mol, dens, zero_h, g, error, minus=minus, &
                                  direct=direct, eri=eri, bounds=bounds, &
                                  k_scale=k_scale, xc=xc, reference=reference, &
                                  rs_k_lr=rs_k_lr, rs_omega=rs_omega, cache=cache, &
                                  triplet=spin_flip, density_screen=density_screen)
         if (error%has_error()) return
      end if

      call cpu_time(t1)
      if (present(t_fock)) t_fock = t_fock + (t1 - t0)
      t0 = t1

      do m = 1, nact
         j = idx(m)
         call pic_gemm(g(:, :, m), c_occ, work)
         call pic_gemm(c_vir, work, au(:, :, j), transa="T")
         au(:, :, j) = gaps*u(:, :, j) + 2.0_dp*au(:, :, j)
      end do

      call cpu_time(t1)
      if (present(t_back)) t_back = t_back + (t1 - t0)

      deallocate (dens, g, half, work)
   end subroutine response_product

   subroutine response_mean_field_df(b, x, c_occ, dtilde, g, minus, k_scale)
      !! `J - k K/2` for a response density, from the fitted tensor
      !!
      !! Not `build_fock_df`: that one assumes an idempotent SCF density to get
      !! its occupied orbitals from, and a response density is symmetric but
      !! *indefinite*, so it has none. It arrives already factored instead,
      !!
      !!     Dt = X C_occ^T +/- C_occ X^T,     X = C_vir u
      !!
      !! which turns exchange into
      !!
      !!     K = sum_P [ (B_P X)(B_P C_occ)^T +/- (B_P C_occ)(B_P X)^T ]
      !!
      !! two n^2 n_occ products per auxiliary function. `Dt` is still passed,
      !! but only for the Coulomb term, where any density will do -- and the
      !! antisymmetric combination does not need it at all.
      real(dp), intent(in) :: b(:, :)        !! `B(mu nu, P)`, (n_ao^2, naux)
      real(dp), intent(in) :: x(:, :)        !! `C_vir u`, (n_ao, n_occ)
      real(dp), intent(in) :: c_occ(:, :)    !! (n_ao, n_occ)
      real(dp), intent(in) :: dtilde(:, :)   !! The assembled response density, for J
      real(dp), intent(out) :: g(:, :)
      logical, intent(in) :: minus
         !! The antisymmetric combination: the second exchange term changes
         !! sign and there is no Coulomb term.
      real(dp), intent(in), optional :: k_scale
         !! The exchange fraction the reference kept. Absent is all of it.

      real(dp), allocatable :: coul(:, :), exch(:, :), bx(:, :), bc(:, :), b_p(:, :)
      real(dp), allocatable :: coul_t(:, :), exch_t(:, :)
      real(dp) :: c_p, kf, sign_two
      integer :: n, n_occ, naux, p

      kf = 1.0_dp
      if (present(k_scale)) kf = k_scale
      sign_two = 1.0_dp
      if (minus) sign_two = -1.0_dp

      n = size(c_occ, 1)
      n_occ = size(c_occ, 2)
      naux = size(b, 2)
      allocate (coul(n, n), exch(n, n))

      coul = 0.0_dp
      exch = 0.0_dp
      ! Threaded over the auxiliary functions, each thread with its own pair of
      ! `n^2` accumulators, because the BLAS is sequential: the Z-vector
      ! iterations of a fitted-reference gradient ran on one core for 36 s here.
      !$omp parallel default(none) &
      !$omp    shared(b, x, c_occ, dtilde, coul, exch, n, n_occ, naux, minus, sign_two) &
      !$omp    private(p, c_p, b_p, bx, bc, coul_t, exch_t)
      allocate (b_p(n, n), bx(n, n_occ), bc(n, n_occ), coul_t(n, n), exch_t(n, n))
      coul_t = 0.0_dp
      exch_t = 0.0_dp
      !$omp do schedule(static)
      do p = 1, naux
         b_p = reshape(b(:, p), [n, n])
         if (.not. minus) then
            c_p = sum(b_p*dtilde)
            coul_t = coul_t + c_p*b_p
         end if
         call pic_gemm(b_p, x, bx)
         call pic_gemm(b_p, c_occ, bc)
         call pic_gemm(bx, bc, exch_t, transb="T", alpha=1.0_dp, beta=1.0_dp)
         call pic_gemm(bc, bx, exch_t, transb="T", alpha=sign_two, beta=1.0_dp)
      end do
      !$omp end do
      !$omp critical
      coul = coul + coul_t
      exch = exch + exch_t
      !$omp end critical
      deallocate (b_p, bx, bc, coul_t, exch_t)
      !$omp end parallel

      g = coul - 0.5_dp*kf*exch
      deallocate (coul, exch)
   end subroutine response_mean_field_df

   subroutine response_mean_field_uhf(mol, dens_a, dens_b, zero_h, bounds, g_a, g_b, &
                                      error, minus, k_scale, xc, ref_a, ref_b, &
                                      rs_k_lr, rs_omega, cache, stats)
      !! `G_sigma(D')` for a batch of response density pairs, over one integral pass
      !!
      !! The unrestricted counterpart of `response_mean_field`:
      !!
      !!     G_a = J[D'_a + D'_b] - c_x K[D'_a] + sum_t f_xc^(a,t) . drho_t
      !!
      !! and the same with the spins exchanged. **Full** same-spin exchange,
      !! not the half a closed-shell build carries, and the Coulomb term
      !! reading the sum of the two -- which is what makes the two spin blocks
      !! couple through `J` and through the kernel and through nothing else.
      !!
      !! **Integral-direct only.** There is no stored-tensor branch and no
      !! fitted one: a fitted reference is refused before the SCF runs, and
      !! the in-core path has no unrestricted batched build to read. A caller
      !! that wants either has to add it rather than have this silently route
      !! around the missing term.
      !!
      !! **What `minus` changes** is what it changes on the restricted side:
      !! `A - B` acts on the antisymmetrised pair, where the Coulomb term
      !! vanishes identically and so does the kernel, whose response is to a
      !! density change an antisymmetric matrix does not make. So `minus` is
      !! exchange only, both ranges of it, and the grid is never touched --
      !! and, cross-spin exchange being nothing, `(A - B)` is block diagonal
      !! in the spin.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: dens_a(:, :, :), dens_b(:, :, :)
         !! `(n_ao, n_ao, n_set)` each, already symmetrised or antisymmetrised
      real(dp), intent(in) :: zero_h(:, :)
         !! Added to every set, so a zero matrix returns `G` alone
      real(dp), intent(in) :: bounds(:, :)
         !! From `schwarz_bounds`. Required, unlike the restricted twin's:
         !! there is no stored-tensor branch here for an absent one to mean.
      real(dp), allocatable, intent(out) :: g_a(:, :, :), g_b(:, :, :)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: minus
         !! The densities are antisymmetric and this is the `A - B` half.
         !! Off by default.
      real(dp), intent(in), optional :: k_scale
         !! The exact-exchange fraction the reference kept. Absent, it is
         !! `xc%exx_fraction` when an `xc` context was given and one otherwise.
      type(xc_context_t), intent(inout), optional :: xc
         !! Present, the spin-resolved kernel is added. Needs `ref_a` and
         !! `ref_b`, and a **polarised** context. Ignored when `minus`.
      real(dp), intent(in), optional :: ref_a(:, :), ref_b(:, :)
         !! The converged reference spin densities the kernel is evaluated at
      real(dp), intent(in), optional :: rs_k_lr, rs_omega
         !! A range-separated functional's attenuated second exchange pass.
         !! Absent, both come from `xc` where it says it is range separated.
      type(xc_kernel_cache_uks_t), intent(in), optional :: cache
         !! The reference's polarised kernel coefficients, from
         !! `xc_kernel_cache_uks_fill`. A cache that was never filled is
         !! ignored, so a caller holding one as a plain component may pass it
         !! unconditionally.
      type(direct_stats_t), intent(out), optional :: stats
         !! Quartets computed and skipped, summed over every integral pass made.

      real(dp), allocatable :: ga_lr(:, :, :), gb_lr(:, :, :)
      type(direct_stats_t) :: pass
      real(dp) :: kf, k_lr, omega
      logical :: anti, use_cache

      if (error%has_error()) return

      anti = .false.
      if (present(minus)) anti = minus
      use_cache = .false.
      if (present(cache)) use_cache = cache%filled

      kf = 1.0_dp
      k_lr = 0.0_dp
      omega = 0.0_dp
      if (present(xc)) then
         kf = xc%exx_fraction
         if (xc%range_separated) then
            k_lr = xc%rs_k_lr
            omega = xc%rs_omega
         end if
      end if
      if (present(k_scale)) kf = k_scale
      if (present(rs_k_lr)) k_lr = rs_k_lr
      if (present(rs_omega)) omega = rs_omega

      if (present(xc) .and. .not. (present(ref_a) .and. present(ref_b))) then
         call error%set(ERROR_VALIDATION, "the unrestricted response operator was "// &
                        "given an exchange-correlation context but not both spin "// &
                        "densities to evaluate its kernel at")
         return
      end if

      if (present(stats)) then
         stats%quartets_total = 0_int64
         stats%quartets_computed = 0_int64
         stats%quartets_screened = 0_int64
      end if

      ! Announced antisymmetric, the folded accumulation is exact: it drops
      ! the Coulomb term, which vanishes, and antisymmetrises the result.
      call build_fock_direct_uhf_many(mol, zero_h, dens_a, dens_b, bounds, g_a, g_b, &
                                      pass, error, k_scale=kf, antisymmetric=anti)
      if (error%has_error()) return
      call add_stats(stats, pass)
      if (omega > 0.0_dp) then
         call build_fock_direct_uhf_many(mol, zero_h, dens_a, dens_b, bounds, ga_lr, &
                                         gb_lr, pass, error, k_scale=k_lr, &
                                         j_scale=0.0_dp, omega=omega, antisymmetric=anti)
         if (error%has_error()) return
         g_a = g_a + ga_lr
         g_b = g_b + gb_lr
         call add_stats(stats, pass)
         deallocate (ga_lr, gb_lr)
      end if

      ! The kernel, for a Kohn-Sham reference, on top of whatever built `G`.
      ! One grid pass serves the whole batch and both spins.
      !
      ! No VV10: a non-local reference is refused before the SCF, because its
      ! kernel here would be a term silently left out rather than one this
      ! routine could supply.
      if (present(xc) .and. .not. anti) then
         if (use_cache) then
            call xc_kernel_apply_uks_many(xc, mol, ref_a, ref_b, dens_a, dens_b, &
                                          g_a, g_b, error, cache=cache)
         else
            call xc_kernel_apply_uks_many(xc, mol, ref_a, ref_b, dens_a, dens_b, &
                                          g_a, g_b, error)
         end if
         if (error%has_error()) return
      end if
   end subroutine response_mean_field_uhf

   subroutine response_product_uhf(mol, c_occ_a, c_vir_a, c_occ_b, c_vir_b, gaps_a, &
                                   gaps_b, zero_h, bounds, u_a, u_b, minus, au_a, &
                                   au_b, error, k_scale, xc, ref_a, ref_b, rs_k_lr, &
                                   rs_omega, cache)
      !! `(A+B)u` or `(A-B)u` for many spin-blocked trial rotations, one pass
      !!
      !! The unrestricted counterpart of `response_product`. Writing each
      !! spin's trial rotation as a density,
      !!
      !!     D'_s = C_vir,s u_s C_occ,s^T  +/-  transpose
      !!
      !! the image is
      !!
      !!     (A +/- B) u|_s = (eps_a - eps_i)_s u_s + C_vir,s^T G_s(D') C_occ,s
      !!
      !! -- **no factor of two** in front of the mean field, where the
      !! restricted product carries one. That two is the two electrons a
      !! closed-shell spatial orbital holds; here each spin is its own, and
      !! the same two appears instead inside `G` as full rather than half
      !! same-spin exchange. Feed this `u_a = u_b` on a closed shell and the
      !! restricted `(A+B)` comes back exactly, which is what the cross-check
      !! test asserts.
      !!
      !! **A spin with no rotations is allowed**, which is what a reference
      !! with no beta electrons has: that spin's trial and image rectangles
      !! are then empty, its response density is zero, and the other spin
      !! still sees it through `J`. The transforms are skipped rather than
      !! run at zero extent, so no BLAS call is made with a vanishing inner
      !! dimension.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: c_occ_a(:, :), c_vir_a(:, :)   !! (n_ao, n_occ_a), (n_ao, n_vir_a)
      real(dp), intent(in) :: c_occ_b(:, :), c_vir_b(:, :)
      real(dp), intent(in) :: gaps_a(:, :), gaps_b(:, :)
         !! `eps_a - eps_i` per spin, `(n_vir_s, n_occ_s)`
      real(dp), intent(in) :: zero_h(:, :)
         !! Zero, so the two-electron build returns `G` alone
      real(dp), intent(in) :: bounds(:, :)   !! From `schwarz_bounds`
      real(dp), intent(in) :: u_a(:, :, :), u_b(:, :, :)
         !! (n_vir_s, n_occ_s, n_set) the trial rotations
      logical, intent(in) :: minus          !! `A - B` rather than `A + B`
      real(dp), intent(out) :: au_a(:, :, :), au_b(:, :, :)
         !! (n_vir_s, n_occ_s, n_set) the images
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: k_scale
      type(xc_context_t), intent(inout), optional :: xc
      real(dp), intent(in), optional :: ref_a(:, :), ref_b(:, :)
      real(dp), intent(in), optional :: rs_k_lr, rs_omega
      type(xc_kernel_cache_uks_t), intent(in), optional :: cache

      real(dp), allocatable :: da(:, :, :), db(:, :, :), ga(:, :, :), gb(:, :, :)
      real(dp), allocatable :: half_a(:, :), half_b(:, :), work_a(:, :), work_b(:, :)
      integer :: n_ao, n_occ_a, n_occ_b, n_set, m
      logical :: has_alpha, has_beta

      if (error%has_error()) return

      n_ao = size(c_occ_a, 1)
      n_occ_a = size(c_occ_a, 2)
      n_occ_b = size(c_occ_b, 2)
      n_set = size(u_a, 3)
      if (n_set <= 0) return
      has_alpha = n_occ_a > 0 .and. size(c_vir_a, 2) > 0
      has_beta = n_occ_b > 0 .and. size(c_vir_b, 2) > 0

      allocate (da(n_ao, n_ao, n_set), db(n_ao, n_ao, n_set))
      allocate (half_a(n_ao, n_occ_a), half_b(n_ao, n_occ_b))
      allocate (work_a(n_ao, n_occ_a), work_b(n_ao, n_occ_b))

      do m = 1, n_set
         if (has_alpha) then
            call pic_gemm(c_vir_a, u_a(:, :, m), half_a)
            call pic_gemm(half_a, c_occ_a, da(:, :, m), transb="T")
         else
            da(:, :, m) = 0.0_dp
         end if
         if (has_beta) then
            call pic_gemm(c_vir_b, u_b(:, :, m), half_b)
            call pic_gemm(half_b, c_occ_b, db(:, :, m), transb="T")
         else
            db(:, :, m) = 0.0_dp
         end if
         if (minus) then
            da(:, :, m) = da(:, :, m) - transpose(da(:, :, m))
            db(:, :, m) = db(:, :, m) - transpose(db(:, :, m))
         else
            da(:, :, m) = da(:, :, m) + transpose(da(:, :, m))
            db(:, :, m) = db(:, :, m) + transpose(db(:, :, m))
         end if
      end do

      call response_mean_field_uhf(mol, da, db, zero_h, bounds, ga, gb, error, &
                                   minus=minus, k_scale=k_scale, xc=xc, &
                                   ref_a=ref_a, ref_b=ref_b, rs_k_lr=rs_k_lr, &
                                   rs_omega=rs_omega, cache=cache)
      if (error%has_error()) then
         deallocate (da, db, half_a, half_b, work_a, work_b)
         return
      end if

      do m = 1, n_set
         if (has_alpha) then
            call pic_gemm(ga(:, :, m), c_occ_a, work_a)
            call pic_gemm(c_vir_a, work_a, au_a(:, :, m), transa="T")
            au_a(:, :, m) = gaps_a*u_a(:, :, m) + au_a(:, :, m)
         end if
         if (has_beta) then
            call pic_gemm(gb(:, :, m), c_occ_b, work_b)
            call pic_gemm(c_vir_b, work_b, au_b(:, :, m), transa="T")
            au_b(:, :, m) = gaps_b*u_b(:, :, m) + au_b(:, :, m)
         end if
      end do

      deallocate (da, db, ga, gb, half_a, half_b, work_a, work_b)
   end subroutine response_product_uhf

end module mqc_czt_response_product
