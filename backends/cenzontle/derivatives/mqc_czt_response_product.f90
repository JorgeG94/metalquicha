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
   use pic_types, only: dp, int64
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_xc, only: xc_context_t, xc_kernel_apply_many, vv10_kernel_apply, &
                         xc_kernel_cache_t
   use mqc_czt_direct, only: build_fock, build_fock_direct, build_fock_direct_many, &
                             direct_stats_t
   implicit none
   private

   public :: response_mean_field
   public :: response_product
   ! `response_mean_field_df` stays private, as `response_operator_df` was
   ! before it: the fitted build is reached through `response_product`, which
   ! is where the factored form it needs is assembled.

contains

   subroutine response_mean_field(mol, dens, zero_h, g, error, minus, direct, eri, &
                                  bounds, k_scale, xc, reference, rs_k_lr, rs_omega, &
                                  cache, density_screen, screen_floor, stats)
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
      logical :: anti, is_direct, screen, use_cache

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
         ! TODO(mqc): unreachable today, because the only route that reaches
         ! here with a stored tensor is a double-hybrid Z-vector solve and all
         ! three supported double hybrids are global hybrids. Add a
         ! range-separated one and this turns a working gradient into an abort:
         ! give the stored route an attenuated companion tensor, or send a
         ! range-separated reference down the direct path, before that lands.
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
         call direct_pass(mol, zero_h, dens, bounds, g, pass, error, k_scale=kf, &
                          j_scale=1.0_dp, omega=0.0_dp, antisymmetric=anti, &
                          density_screen=screen)
         if (error%has_error()) return
         call add_stats(stats, pass)
         if (omega > 0.0_dp) then
            call direct_pass(mol, zero_h, dens, bounds, g_lr, pass, error, &
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
            call build_fock(zero_h, eri, dens(:, :, p), g(:, :, p), k_scale=kf)
         end do
      end if

      ! The kernel, for a Kohn-Sham reference, on top of whatever built `g`.
      ! Leaving it out does not fail; it converges to the wrong orbital
      ! response. One grid pass serves the whole batch.
      if (present(xc) .and. .not. anti) then
         if (use_cache) then
            call xc_kernel_apply_many(xc, mol, reference, dens, g, error, cache=cache)
         else
            call xc_kernel_apply_many(xc, mol, reference, dens, g, error)
         end if
         if (error%has_error()) return
         ! The non-local kernel, once for the batch rather than per set:
         ! `vv10_kernel_apply`'s pair sweep is O(npts^2) whether it carries one
         ! trial density or a dozen. It accumulates, hence the zeroed buffer.
         !
         ! TODO(mqc): this one is not cached, so on a VV10 reference the kernel
         ! cache removes the smaller of the two grid costs an application pays
         ! and leaves the larger. What is reusable is the reference half --
         ! `vv10_nlc`'s U..C pair sums and the omega and kappa derivatives,
         ! plus `rho`, `sigma` and `rho_grad` -- thirteen more grid-sized
         ! arrays over the NLC grid, which is not the grid `xc_kernel_cache_t`
         ! is filled over, so it would be a second cache with its own budget
         ! rather than four more components on that one. The trial half,
         ! `vv10_hessian_kernel`, cannot be cached at all: it is
         ! O(npts^2 n_set) and moves with the densities. So the saving is one
         ! pair sweep out of two at `n_set = 1` and one out of `n_set + 1` as
         ! the batch widens -- most to a coupled-perturbed solve, least to the
         ! Hessian's wide batches.
         if (xc%nlc_b /= 0.0_dp .or. xc%nlc_c /= 0.0_dp) then
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

   subroutine direct_pass(mol, zero_h, dens, bounds, g, stats, error, k_scale, &
                          j_scale, omega, antisymmetric, density_screen)
      !! One integral-direct build over the batch, by the cheaper of the two routes
      !!
      !! `build_fock_direct_many` is bit-for-bit equal to `build_fock_direct` on
      !! a single symmetric density, and not free: its accumulator is indexed by
      !! the set, so every update is a length-one vector operation, and it
      !! carries its tiling, its lock array and its window buffers whatever the
      !! width. Measured on the double-hybrid Hessian's Z-vector solves, one
      !! density at a time through the batched build costs 16 per cent more at
      !! 19 basis functions and 12 per cent at 29, against the same answer to
      !! the last bit. The coupled-perturbed solvers apply the operator to one
      !! trial vector at a time, so that is their whole matvec.
      !!
      !! An antisymmetric single density still goes through the batched build:
      !! folding that symmetry in is what `antisymmetric` does and the
      !! single-density routine has no such argument.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: zero_h(:, :)
      real(dp), intent(in) :: dens(:, :, :)
      real(dp), intent(in) :: bounds(:, :)
      real(dp), allocatable, intent(out) :: g(:, :, :)
      type(direct_stats_t), intent(out) :: stats
      type(error_t), intent(inout) :: error
      real(dp), intent(in) :: k_scale, j_scale, omega
      logical, intent(in) :: antisymmetric, density_screen

      if (size(dens, 3) == 1 .and. .not. antisymmetric) then
         allocate (g(size(dens, 1), size(dens, 2), 1))
         call build_fock_direct(mol, zero_h, dens(:, :, 1), bounds, g(:, :, 1), &
                                stats, error, k_scale=k_scale, j_scale=j_scale, &
                                omega=omega, density_screen=density_screen)
      else
         call build_fock_direct_many(mol, zero_h, dens, bounds, g, stats, error, &
                                     k_scale=k_scale, j_scale=j_scale, omega=omega, &
                                     antisymmetric=antisymmetric, &
                                     density_screen=density_screen)
      end if
   end subroutine direct_pass

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
                               rs_k_lr, rs_omega, cache, bmat, density_screen, &
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
      logical :: use_cache

      if (error%has_error()) return
      if (nact <= 0) return

      use_cache = .false.
      if (present(cache)) use_cache = cache%filled

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
         if (present(xc)) then
            ! TODO(mqc): unreachable today for the same reason the stored-tensor
            ! refusal in `response_mean_field` is -- a fitted reference reaches
            ! here only from a double-hybrid Z-vector solve, and all three
            ! supported double hybrids are global hybrids. A range-separated
            ! one would turn a working gradient into an abort rather than
            ! quietly dropping a term, which is the right order, but it still
            ! needs an attenuated fitting before it can be offered.
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
                                  density_screen=density_screen)
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

end module mqc_czt_response_product
