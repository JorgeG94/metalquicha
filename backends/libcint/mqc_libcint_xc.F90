!! Exchange-correlation, as one thing the SCF can be handed
module mqc_libcint_xc
   !! Everything the Kohn-Sham SCF needs that Hartree-Fock does not, behind one
   !! derived type and one call.
   !!
   !! The SCF gains one optional argument of type `xc_context_t` and one call
   !! inside its Fock build; Hartree-Fock is the case where that argument is
   !! absent, rather than a separate path.
   !!
   !! The context owns the grid, the resolved libxc functionals and their
   !! weights, and the exact-exchange fraction. It does **not** own basis
   !! function values: those are rebuilt per block of grid points on every
   !! iteration rather than cached, because cached they are `n_points` by
   !! `n_ao`.
   !!
   !! LDA, GGA and meta-GGA are all evaluated, restricted and spin-polarised.
   !! There are two evaluators, `xc_add_potential` and `xc_add_potential_uks`,
   !! because libxc fixes the spin channel when a functional is initialised and
   !! the polarised arrays are spin-interleaved, so the two cases take different
   !! array shapes. What they share -- the three terms that turn a pointwise
   !! potential into a matrix -- is `accumulate_xc_matrix`.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_vv10, only: vv10_nlc, vv10_hessian_kernel
   use mqc_dft_grid, only: dft_grid_t, build_dft_grid, DEFAULT_GRID_LEVEL
   use mqc_xc_spec, only: xc_spec_t, xc_spec_from_name, MAX_XC_COMPONENTS
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_libcint_ao, only: eval_ao_block, eval_rho, AO_POINT_BLOCK, &
                             shell_extents, block_significant_aos
#ifdef MQC_WITH_LIBXC
   use xc_f03_lib_m, only: xc_f03_func_t, xc_f03_func_init, xc_f03_func_end, &
                           xc_f03_lda_exc_vxc, xc_f03_func_get_info, &
                           xc_f03_func_info_get_family, xc_f03_hyb_exx_coef, &
                           xc_f03_hyb_cam_coef, xc_f03_nlc_coef, &
                           xc_f03_functional_get_number, xc_f03_func_info_t, &
                           xc_f03_gga_exc_vxc, xc_f03_mgga_exc_vxc, &
                           xc_f03_lda_fxc, xc_f03_gga_fxc, xc_f03_mgga_fxc, &
                           xc_f03_lda_kxc, xc_f03_gga_kxc, &
                           xc_f03_func_info_get_flags, XC_FLAGS_NEEDS_LAPLACIAN, &
                           XC_UNPOLARIZED, XC_POLARIZED, &
                           XC_FAMILY_LDA, XC_FAMILY_HYB_LDA, &
                           XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, &
                           XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA
#endif
   implicit none
   private

   real(dp), parameter :: AO_SCREEN_TOL = 1.0e-12_dp
      !! The AO value below which a shell is dropped from a block.
      !!
      !! Tighter than the 1e-10 often quoted: the quantity being converged is a
      !! total energy at the microhartree level over ~10^6 points, and the extra
      !! digit costs a couple of bohr of radius on the most diffuse shells.

   integer, parameter, public :: NLC_GRID_LEVEL = 1
      !! Default level for the non-local inner grid.
      !!
      !! Not the exchange grid's level: PySCF defaults its `nlcgrids` to the
      !! main grid's, which costs the full product, and this is one step down.
      !! A deck can raise it back.

   public :: xc_context_t
   public :: xc_context_create
   public :: xc_add_potential
   public :: xc_add_potential_uks
   public :: xc_available
   public :: xc_grid_lda_quantities
   public :: xc_grid_gga_quantities
   public :: xc_kernel_apply
   public :: xc_kernel2_apply
   public :: xc_grid_kernel_quantities
   public :: KERNEL_RHO_FLOOR   !! Where the kernel's divergence is cut off
   public :: ensure_nlc_grid
   public :: vv10_add_potential
   public :: vv10_kernel_apply

   real(dp), parameter :: KERNEL_RHO_FLOOR = 1.0e-10_dp
      !! Grid points below this density contribute no *second* derivative.
      !!
      !! Not a tolerance -- a necessity. The LDA kernel is `d2e/drho2`, which for
      !! the exchange part goes as `rho^(-2/3)` and diverges at the tail of every
      !! atomic grid. Left in, those points swamp the real ones and the response
      !! operator stops being positive definite, so the conjugate-gradient solver
      !! reports a saddle point that is not there. `vrho` and `exc` stay finite
      !! where `v2rho2` does not, which is why nothing before the kernel needed a
      !! floor.

   type :: xc_context_t
      !! A functional, a grid, and the fraction of exact exchange to keep
      logical :: active = .false.
      real(dp) :: exx_fraction = 0.0_dp
         !! How much Fock exchange the functional wants. Zero for pure DFT, one
         !! for Hartree-Fock, in between for a hybrid. Read from libxc where libxc
         !! owns the functional, and from `mqc_xc_spec` only where it does not.
      logical :: range_separated = .false.
         !! Whether the exchange splits into short and long range. The Fock matrix
         !! then needs two exchange passes rather than one, so the SCF has to know.
      real(dp) :: rs_omega = 0.0_dp
         !! The range-separation parameter, for libcint's `env(PTR_RANGE_OMEGA)`.
      real(dp) :: rs_k_lr = 0.0_dp
         !! Coefficient of the long-range exchange matrix. `exx_fraction` carries
         !! the coefficient of the full one, so that
         !!
         !!     K_eff = exx_fraction * K_full + rs_k_lr * K_lr(rs_omega)
         !!
         !! libxc reports (omega, alpha, beta) with alpha the long-range
         !! coefficient and alpha+beta the short-range one, so `exx_fraction` is
         !! alpha+beta and this is -beta.
      real(dp) :: pt2_fraction = 0.0_dp
         !! MP2 correlation fraction, for a double hybrid. Nothing in this module
         !! acts on it -- perturbative correlation is not a grid quantity.
      type(dft_grid_t) :: nlc_grid
         !! A second, coarser quadrature for VV10's inner sum only.
         !!
         !! The non-local term is a double integral, so its cost goes as the
         !! product of the two grids' sizes while everything else here is linear
         !! in one.
      integer :: nlc_grid_level = NLC_GRID_LEVEL
      real(dp) :: nlc_b = 0.0_dp
      real(dp) :: nlc_c = 0.0_dp
         !! VV10's two parameters, as libxc reports them. Both zero means the
         !! functional carries no non-local correlation, which is every
         !! functional here except the `-V` family.
      real(dp) :: screen_tol = AO_SCREEN_TOL
         !! The AO value below which a shell is dropped from a grid block.
         !! From `keywords.dft.screening_tolerance`; the default is the constant.
      integer :: point_block = AO_POINT_BLOCK
         !! Grid points per block. From `keywords.dft.block_size`.
         !!
         !! Two things at once. A smaller block is spatially tighter, so fewer
         !! shells reach it and the screen keeps less; and the loop over blocks
         !! *is* the OpenMP loop, so the block count is the thread granularity.
      logical :: polarized = .false.
         !! Whether the functionals were initialised spin-polarised.
         !!
         !! libxc fixes this at initialisation, and the two cases take different
         !! array shapes and return a different number of potential components.
         !! So a context belongs to a restricted calculation or to an
         !! unrestricted one, and `xc_add_potential` and `xc_add_potential_uks`
         !! each refuse the other kind rather than reinterpreting the arrays --
         !! which would be a silent read of the wrong element, not a crash.
      type(dft_grid_t) :: grid
      integer :: n_func = 0
      integer :: family(MAX_XC_COMPONENTS) = 0
         !! libxc's family per component. A composition may mix them -- a hybrid
         !! GGA exchange beside an LDA correlation -- so the block loop asks per
         !! component rather than once for the whole functional.
      logical :: any_gga = .false.
         !! Whether anything here needs density gradients, and so whether the AO
         !! gradients are worth evaluating at all.
      logical :: any_mgga = .false.
         !! Whether anything here needs the kinetic energy density too. Separate
         !! from `any_gga` because a meta-GGA needs both, and the gradients are the
         !! expensive half.
      real(dp) :: weight(MAX_XC_COMPONENTS) = 0.0_dp
#ifdef MQC_WITH_LIBXC
      type(xc_f03_func_t) :: func(MAX_XC_COMPONENTS)
#endif
   contains
      procedure :: destroy => xc_context_destroy
   end type xc_context_t

contains

   pure function xc_available() result(available)
      !! Whether this build can evaluate a functional at all
      logical :: available
#ifdef MQC_WITH_LIBXC
      available = .true.
#else
      available = .false.
#endif
   end function xc_available

   subroutine xc_context_create(mol, functional, ctx, error, level, polarized, &
                                screen_tol, point_block, nlc_level, allow_half)
      !! Resolve a functional name and build the grid it will be integrated on
      type(libcint_molecule_t), intent(in) :: mol
      character(len=*), intent(in) :: functional
      type(xc_context_t), intent(out) :: ctx
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: level
      logical, intent(in), optional :: polarized
         !! Initialise the functionals spin-polarised, for an unrestricted
         !! calculation. Default restricted. Fixed here because libxc fixes it at
         !! initialisation, so it cannot be decided later by whoever evaluates.
      integer, intent(in), optional :: nlc_level
         !! Level for VV10's own quadrature. Absent or negative keeps
         !! `NLC_GRID_LEVEL`; a deck reaches this through
         !! `keywords.dft.nlc_grid_level`.
      real(dp), intent(in), optional :: screen_tol
         !! AO screening threshold; zero or negative disables the screen and
         !! evaluates the whole basis, which is the way to check what it costs.
      integer, intent(in), optional :: point_block
         !! Grid points per block. Non-positive keeps the default.
      logical, intent(in), optional :: allow_half
         !! Accept a libxc name carrying only exchange or only correlation.
         !! Passed to `xc_spec_from_name`, which explains why it exists and why
         !! no deck may set it. The derivative tests are the only caller.

      type(xc_spec_t) :: spec
      integer :: grid_level, i, id, family
      integer, allocatable :: numbers(:)
#ifdef MQC_WITH_LIBXC
      type(xc_f03_func_info_t) :: info
      real(dp) :: libxc_exx
      real(dp) :: cam_omega, cam_alpha, cam_beta, nlc_b, nlc_c
#endif

      if (.not. xc_available()) then
         call error%set(ERROR_VALIDATION, "this build has no exchange-correlation "// &
                        "functionals: configure with -DMQC_ENABLE_LIBXC=ON")
         return
      end if

      call xc_spec_from_name(functional, spec, error, allow_half=allow_half)
      if (error%has_error()) return

      grid_level = DEFAULT_GRID_LEVEL
      if (present(level)) grid_level = level

      ctx%nlc_grid_level = NLC_GRID_LEVEL
      if (present(nlc_level)) then
         if (nlc_level >= 0) ctx%nlc_grid_level = nlc_level
      end if

      allocate (numbers(mol%natm))
      ! The *element*, not the charge it presents. `mol%charges` has the ECP's
      ! core subtracted, and every table the grid reaches for is an element
      ! property: the radial count and Lebedev order come from the period, the
      ! Treutler-Ahlrichs xi and the Becke radii are per element. The reduced
      ! charge builds iodine's grid as though it were manganese, 53 - 28 = 25.
      ! `core_electrons` is zero for an all-electron atom and for a ghost, so
      ! this is the identity everywhere an ECP is not involved.
      numbers = nint(mol%charges) + mol%core_electrons
      call build_dft_grid(mol%coords, numbers, ctx%grid, error, level=grid_level)
      if (error%has_error()) return
      deallocate (numbers)
      if (error%has_error()) return

      ctx%n_func = spec%n_components
      ctx%pt2_fraction = spec%pt2_fraction
      ctx%exx_fraction = spec%exx_fraction
      if (present(polarized)) ctx%polarized = polarized

#ifdef MQC_WITH_LIBXC
      do i = 1, spec%n_components
         ctx%weight(i) = spec%component(i)%weight
         id = xc_f03_functional_get_number(trim(spec%component(i)%name))
         if (id <= 0) then
            call error%set(ERROR_VALIDATION, "libxc does not know the functional '"// &
                           trim(spec%component(i)%name)//"'")
            return
         end if
         if (ctx%polarized) then
            call xc_f03_func_init(ctx%func(i), id, XC_POLARIZED)
         else
            call xc_f03_func_init(ctx%func(i), id, XC_UNPOLARIZED)
         end if

         info = xc_f03_func_get_info(ctx%func(i))
         family = xc_f03_func_info_get_family(info)
         ctx%family(i) = family
         select case (family)
         case (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
            ! nothing extra
         case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
            ctx%any_gga = .true.
         case (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
            ! A meta-GGA needs the gradients as well, since sigma appears alongside
            ! tau in every one of them.
            ctx%any_gga = .true.
            ctx%any_mgga = .true.
            ! Some meta-GGAs want the density Laplacian, which is a second
            ! derivative of every basis function and not implemented. libxc says
            ! which, so this is a refusal on the functional's own account rather
            ! than a guess about which ones are safe.
            if (iand(xc_f03_func_info_get_flags(info), XC_FLAGS_NEEDS_LAPLACIAN) /= 0) then
               call error%set(ERROR_VALIDATION, "'"//trim(spec%component(i)%name)// &
                              "' needs the density Laplacian, which the CPU path does "// &
                              "not compute. Refused rather than passed zeros, which "// &
                              "would return a different functional under its name.")
               return
            end if
         case default
            call error%set(ERROR_VALIDATION, "'"//trim(spec%component(i)%name)// &
                           "' is beyond meta-GGA, which the CPU path does not "// &
                           "implement.")
            return
         end select

         ! **Range separation, detected rather than asked about.** libxc 7.1.2's
         ! Fortran bindings expose no `xc_hyb_type`, so the test is on the
         ! coefficients: a global hybrid reports omega = 0. Detecting it matters
         ! because `hyb_exx_coef` does not mean the same thing for a
         ! range-separated functional as for a global one, and taking it at face
         ! value converges to a badly wrong energy without a warning.
         !
         ! The attenuated integrals need no new integral code: the same entry
         ! points with `env(PTR_RANGE_OMEGA)` set.
         call xc_f03_hyb_cam_coef(ctx%func(i), cam_omega, cam_alpha, cam_beta)
         if (cam_omega /= 0.0_dp .and. cam_beta /= 0.0_dp) then
            if (ctx%range_separated) then
               call error%set(ERROR_VALIDATION, "two range-separated components in "// &
                              "one functional is not something this assembles: their "// &
                              "omegas would have to agree and nothing checks that.")
               return
            end if
            ctx%range_separated = .true.
            ctx%rs_omega = cam_omega
            ! alpha is the long-range coefficient and alpha+beta the short-range
            ! one; K_full carries the short-range value and the long-range matrix
            ! carries the difference.
            ctx%exx_fraction = ctx%exx_fraction &
                               + spec%component(i)%weight*(cam_alpha + cam_beta)
            ctx%rs_k_lr = ctx%rs_k_lr - spec%component(i)%weight*cam_beta
         end if

         ! Non-local correlation: VV10 is a double integral over the density,
         ! not a functional of it at a point, so libxc supplies only the
         ! semilocal half and hands back the two parameters. Kept here and
         ! evaluated by `mqc_libcint_vv10` where the potential is assembled.
         !
         ! A composition mixing two VV10 functionals with different parameters
         ! is refused rather than averaged: `b` and `c` are not linear in the
         ! energy, so there is no meaningful blend of them.
         call xc_f03_nlc_coef(ctx%func(i), nlc_b, nlc_c)
         if (nlc_b /= 0.0_dp .or. nlc_c /= 0.0_dp) then
            if (ctx%nlc_b /= 0.0_dp .and. &
                (ctx%nlc_b /= nlc_b .or. ctx%nlc_c /= nlc_c)) then
               call error%set(ERROR_VALIDATION, "this composition mixes two non-local "// &
                              "correlation functionals with different VV10 parameters, "// &
                              "which cannot be combined into one kernel.")
               return
            end if
            ctx%nlc_b = nlc_b
            ctx%nlc_c = nlc_c
         end if

         ! libxc owns a hybrid's fraction, so ask rather than assume -- of every
         ! component, not only of a functional libxc carries whole. SCAN0 is the
         ! case that separates the two: it is a composition here, because libxc
         ! pairs no correlation with it, but its exchange component
         ! `hyb_mgga_x_scan0` carries the quarter of exact exchange itself.
         !
         ! Safe to ask always: every semilocal component of every composition in
         ! `mqc_xc_spec` reports zero here. A spec's own `exx_fraction` remains
         ! for exchange libxc knows nothing about, which is the double hybrids
         ! and only them.
         if (.not. (cam_omega /= 0.0_dp .and. cam_beta /= 0.0_dp)) then
            libxc_exx = xc_f03_hyb_exx_coef(ctx%func(i))
            ctx%exx_fraction = ctx%exx_fraction + spec%component(i)%weight*libxc_exx
         end if
      end do
#endif

      if (present(screen_tol)) ctx%screen_tol = screen_tol
      if (present(point_block)) then
         if (point_block > 0) ctx%point_block = point_block
      end if
      ctx%active = .true.
   end subroutine xc_context_create

   subroutine xc_context_destroy(this)
      class(xc_context_t), intent(inout) :: this
      integer :: i

#ifdef MQC_WITH_LIBXC
      do i = 1, this%n_func
         call xc_f03_func_end(this%func(i))
      end do
#endif
      call this%grid%destroy()
      call this%nlc_grid%destroy()
      this%n_func = 0
      this%active = .false.
      this%exx_fraction = 0.0_dp
      this%pt2_fraction = 0.0_dp
      this%range_separated = .false.
      this%rs_omega = 0.0_dp
      this%rs_k_lr = 0.0_dp
   end subroutine xc_context_destroy

   subroutine xc_grid_lda_quantities(ctx, mol, density, rho, exc, vrho, error, &
                                     density_beta, rho_beta, vrho_beta)
      !! rho, eps_xc and v_xc at every grid point, for an LDA functional
      !!
      !! The gradient needs these three on the grid rather than contracted into a
      !! matrix: one term contracts `vrho` against the moving basis functions and
      !! another contracts `rho*exc` against the moving quadrature weights, and
      !! neither can be recovered from the Fock-matrix contribution
      !! `xc_add_potential` returns.
      !!
      !! LDA only, and anything else is refused rather than returned incomplete.
      ! TODO(mqc): this loop, and the ones in `xc_grid_gga_quantities` and
      ! `xc_grid_kernel_quantities`, block on the `AO_POINT_BLOCK` constant and
      ! apply no AO screen, where every other loop in this module blocks on
      ! `ctx%point_block`. So `keywords.dft.block_size` and
      ! `keywords.dft.screening_tolerance` are silently ignored on the gradient
      ! and kernel-quantity paths.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)         !! Total density, or alpha
      real(dp), allocatable, intent(out) :: rho(:)  !! Total density on the grid
      real(dp), allocatable, intent(out) :: exc(:)  !! Energy per electron
      real(dp), allocatable, intent(out) :: vrho(:)  !! dE_xc/drho, or the alpha part
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: density_beta(:, :)
      real(dp), allocatable, intent(out), optional :: rho_beta(:)
      real(dp), allocatable, intent(out), optional :: vrho_beta(:)

      real(dp), allocatable :: ao(:, :), rho_blk(:), exc_i(:), vrho_i(:)
      real(dp), allocatable :: rho_a_blk(:), rho_b_blk(:)
      integer :: g0, g1, nb, i, ig, npts
      logical :: unrestricted

      unrestricted = present(density_beta)
      npts = ctx%grid%n_points

      allocate (rho(npts), exc(npts), vrho(npts))
      rho = 0.0_dp
      exc = 0.0_dp
      vrho = 0.0_dp
      if (unrestricted) then
         allocate (rho_beta(npts), vrho_beta(npts))
         rho_beta = 0.0_dp
         vrho_beta = 0.0_dp
      end if

      if (.not. ctx%active) return
      if (.not. xc_available()) then
         call error%set(ERROR_VALIDATION, "no libxc in this build")
         return
      end if
      if (ctx%any_gga .or. ctx%any_mgga) then
         call error%set(ERROR_VALIDATION, "the exchange-correlation gradient is "// &
                        "implemented for LDA functionals only; this one needs the "// &
                        "density gradient")
         return
      end if
      if (unrestricted .neqv. ctx%polarized) then
         call error%set(ERROR_VALIDATION, "xc_grid_lda_quantities: the spin case does "// &
                        "not match how this context was built")
         return
      end if

#ifdef MQC_WITH_LIBXC
      ! `ctx%point_block`, not the module constant: `keywords.dft.block_size`
      ! sets the former and every other loop in this module reads it, so
      ! blocking on the constant here made the keyword a no-op on this path.
      ! These three routines still apply no AO screen, so
      ! `keywords.dft.screening_tolerance` remains inert for them -- adding it
      ! means resizing every downstream array to the significant-AO subset, the
      ! way the potential loops do, which is a change to numerics and belongs on
      ! its own.
      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error)
         if (error%has_error()) return

         if (allocated(exc_i)) deallocate (exc_i, vrho_i)

         if (unrestricted) then
            call eval_rho(ao, density, rho_a_blk)
            call eval_rho(ao, density_beta, rho_b_blk)
            ! libxc's polarised arrays are spin-interleaved.
            if (allocated(rho_blk)) deallocate (rho_blk)
            allocate (rho_blk(2*nb), exc_i(nb), vrho_i(2*nb))
            do ig = 1, nb
               rho_blk(2*ig - 1) = rho_a_blk(ig)
               rho_blk(2*ig) = rho_b_blk(ig)
            end do
            do i = 1, ctx%n_func
               exc_i = 0.0_dp
               vrho_i = 0.0_dp
               call xc_f03_lda_exc_vxc(ctx%func(i), int(nb, 8), rho_blk, exc_i, vrho_i)
               do ig = 1, nb
                  exc(g0 + ig - 1) = exc(g0 + ig - 1) + ctx%weight(i)*exc_i(ig)
                  vrho(g0 + ig - 1) = vrho(g0 + ig - 1) + ctx%weight(i)*vrho_i(2*ig - 1)
                  vrho_beta(g0 + ig - 1) = vrho_beta(g0 + ig - 1) &
                                           + ctx%weight(i)*vrho_i(2*ig)
               end do
            end do
            do ig = 1, nb
               rho(g0 + ig - 1) = rho_a_blk(ig)
               rho_beta(g0 + ig - 1) = rho_b_blk(ig)
            end do
         else
            call eval_rho(ao, density, rho_blk)
            allocate (exc_i(nb), vrho_i(nb))
            do i = 1, ctx%n_func
               exc_i = 0.0_dp
               vrho_i = 0.0_dp
               call xc_f03_lda_exc_vxc(ctx%func(i), int(nb, 8), rho_blk, exc_i, vrho_i)
               do ig = 1, nb
                  exc(g0 + ig - 1) = exc(g0 + ig - 1) + ctx%weight(i)*exc_i(ig)
                  vrho(g0 + ig - 1) = vrho(g0 + ig - 1) + ctx%weight(i)*vrho_i(ig)
               end do
            end do
            do ig = 1, nb
               rho(g0 + ig - 1) = rho_blk(ig)
            end do
         end if
      end do
#else
      call error%set(ERROR_VALIDATION, "no libxc in this build")
#endif
   end subroutine xc_grid_lda_quantities

   subroutine xc_grid_gga_quantities(ctx, mol, density, rho, exc, vrho, grad_coeff, error, &
                                     density_beta, rho_beta, vrho_beta, grad_coeff_beta, &
                                     vtau, sigma_out)
      !! rho, eps_xc, dE/drho and dE/d(grad rho) at every grid point, for a GGA
      !!
      !! The counterpart of `xc_grid_lda_quantities` for a functional that reads
      !! the density gradient, returning one thing that routine does not:
      !! `grad_coeff`, the quantity conjugate to `grad rho`.
      !!
      !! **`grad_coeff` rather than `vsigma` and `rho_grad` separately.** libxc
      !! parametrises a GGA by sigma = |grad rho|^2, so what multiplies
      !! d(grad rho)/dR is 2 vsigma grad rho, and in the polarised case the spins
      !! couple through sigma_ab:
      !!
      !!     dE/d(grad rho_a) = 2 vsigma_aa grad rho_a + vsigma_ab grad rho_b
      !!
      !! Resolving it here means the caller contracts a single vector per spin
      !! and never has to know how many sigma channels there were. It is also the
      !! combination `xc_add_potential` builds for the Fock matrix, so the two
      !! paths cannot disagree about the chain rule.
      !!
      !! LDA functionals are accepted and return `grad_coeff` zero. Meta-GGA
      !! needs `vtau`, and is refused without it.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)         !! Total density, or alpha
      real(dp), allocatable, intent(out) :: rho(:)
      real(dp), allocatable, intent(out) :: exc(:)
      real(dp), allocatable, intent(out) :: vrho(:)
      real(dp), allocatable, intent(out) :: grad_coeff(:, :)  !! (npts, 3)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: density_beta(:, :)
      real(dp), allocatable, intent(out), optional :: rho_beta(:)
      real(dp), allocatable, intent(out), optional :: vrho_beta(:)
      real(dp), allocatable, intent(out), optional :: grad_coeff_beta(:, :)
      real(dp), allocatable, intent(out), optional :: vtau(:)
         !! `dE/dtau` per point, for a meta-GGA. Absent, a meta-GGA is refused
         !! rather than evaluated without it: the dispatch further down would
         !! otherwise treat one as an LDA.
      real(dp), allocatable, intent(out), optional :: sigma_out(:)
         !! |grad rho|^2 on the whole grid. Formed per block here already and
         !! otherwise discarded; VV10 needs it over every point at once, since
         !! its inner sum is not blockable.

      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :)
      real(dp), allocatable :: rho_blk(:), sigma(:), exc_i(:), vrho_i(:), vsigma_i(:)
      real(dp), allocatable :: vsigma(:)
      real(dp), allocatable :: tau_blk(:), lapl(:), vlapl(:), vtau_i(:)
      real(dp), allocatable :: rho_a_blk(:), rho_b_blk(:), grad_a(:, :), grad_b(:, :)
      real(dp), allocatable :: rho_grad(:, :)
      integer :: g0, g1, nb, i, ig, id, npts
      logical :: unrestricted

      unrestricted = present(density_beta)
      npts = ctx%grid%n_points

      allocate (rho(npts), exc(npts), vrho(npts), grad_coeff(npts, 3))
      if (present(sigma_out)) then
         allocate (sigma_out(npts))
         sigma_out = 0.0_dp
      end if
      rho = 0.0_dp
      exc = 0.0_dp
      vrho = 0.0_dp
      grad_coeff = 0.0_dp
      if (unrestricted) then
         allocate (rho_beta(npts), vrho_beta(npts), grad_coeff_beta(npts, 3))
         rho_beta = 0.0_dp
         vrho_beta = 0.0_dp
         grad_coeff_beta = 0.0_dp
      end if

      if (.not. ctx%active) return
      if (.not. xc_available()) then
         call error%set(ERROR_VALIDATION, "no libxc in this build")
         return
      end if
      if (ctx%any_mgga .and. .not. present(vtau)) then
         call error%set(ERROR_VALIDATION, "this functional needs the kinetic energy "// &
                        "density, and the caller asked for no `vtau`. Refused rather "// &
                        "than evaluated without it: the dispatch below would fall "// &
                        "through to the LDA branch and return a converged, wrong "// &
                        "answer.")
         return
      end if
      if (ctx%any_mgga .and. unrestricted) then
         call error%set(ERROR_VALIDATION, "the meta-GGA gradient is implemented for a "// &
                        "restricted reference only; the unrestricted kinetic energy "// &
                        "density carries a term per spin that is not built here.")
         return
      end if
      if (present(vtau)) then
         allocate (vtau(npts))
         vtau = 0.0_dp
      end if
      if (unrestricted .neqv. ctx%polarized) then
         call error%set(ERROR_VALIDATION, "xc_grid_gga_quantities: the spin case does "// &
                        "not match how this context was built")
         return
      end if

#ifdef MQC_WITH_LIBXC
      ! `ctx%point_block`, not the module constant: `keywords.dft.block_size`
      ! sets the former and every other loop in this module reads it, so
      ! blocking on the constant here made the keyword a no-op on this path.
      ! These three routines still apply no AO screen, so
      ! `keywords.dft.screening_tolerance` remains inert for them -- adding it
      ! means resizing every downstream array to the significant-AO subset, the
      ! way the potential loops do, which is a change to numerics and belongs on
      ! its own.
      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, grad=ao_grad)
         if (error%has_error()) return

         if (allocated(exc_i)) deallocate (exc_i)
         if (allocated(vrho_i)) deallocate (vrho_i)
         if (allocated(vsigma_i)) deallocate (vsigma_i)
         if (allocated(vsigma)) deallocate (vsigma)
         if (allocated(sigma)) deallocate (sigma)

         if (unrestricted) then
            call eval_rho(ao, density, rho_a_blk, ao_grad=ao_grad, rho_grad=grad_a)
            call eval_rho(ao, density_beta, rho_b_blk, ao_grad=ao_grad, rho_grad=grad_b)

            ! libxc's polarised arrays are spin-interleaved, and sigma runs
            ! (aa, ab, bb) per point.
            if (allocated(rho_blk)) deallocate (rho_blk)
            allocate (rho_blk(2*nb), sigma(3*nb), exc_i(nb))
            allocate (vrho_i(2*nb), vsigma_i(3*nb), vsigma(3*nb))
            vsigma = 0.0_dp
            do ig = 1, nb
               rho_blk(2*ig - 1) = rho_a_blk(ig)
               rho_blk(2*ig) = rho_b_blk(ig)
               sigma(3*ig - 2) = dot_product(grad_a(ig, :), grad_a(ig, :))
               sigma(3*ig - 1) = dot_product(grad_a(ig, :), grad_b(ig, :))
               sigma(3*ig) = dot_product(grad_b(ig, :), grad_b(ig, :))
            end do

            do i = 1, ctx%n_func
               exc_i = 0.0_dp
               vrho_i = 0.0_dp
               vsigma_i = 0.0_dp
               ! Same dispatch the Fock build uses. A composition may mix an
               ! LDA component with a GGA one -- a hybrid's correlation part
               ! commonly is -- so the family is per functional rather than per
               ! context, and `any_gga` only says at least one needs sigma.
               select case (ctx%family(i))
               case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                  call xc_f03_gga_exc_vxc(ctx%func(i), int(nb, 8), rho_blk, sigma, &
                                          exc_i, vrho_i, vsigma_i)
                  vsigma = vsigma + ctx%weight(i)*vsigma_i
               case default
                  call xc_f03_lda_exc_vxc(ctx%func(i), int(nb, 8), rho_blk, exc_i, vrho_i)
               end select
               do ig = 1, nb
                  exc(g0 + ig - 1) = exc(g0 + ig - 1) + ctx%weight(i)*exc_i(ig)
                  vrho(g0 + ig - 1) = vrho(g0 + ig - 1) + ctx%weight(i)*vrho_i(2*ig - 1)
                  vrho_beta(g0 + ig - 1) = vrho_beta(g0 + ig - 1) &
                                           + ctx%weight(i)*vrho_i(2*ig)
               end do
            end do

            do ig = 1, nb
               rho(g0 + ig - 1) = rho_a_blk(ig)
               rho_beta(g0 + ig - 1) = rho_b_blk(ig)
               do id = 1, 3
                  grad_coeff(g0 + ig - 1, id) = 2.0_dp*vsigma(3*ig - 2)*grad_a(ig, id) &
                                                + vsigma(3*ig - 1)*grad_b(ig, id)
                  grad_coeff_beta(g0 + ig - 1, id) = 2.0_dp*vsigma(3*ig)*grad_b(ig, id) &
                                                     + vsigma(3*ig - 1)*grad_a(ig, id)
               end do
            end do
         else
            if (ctx%any_mgga) then
               call eval_rho(ao, density, rho_blk, ao_grad=ao_grad, rho_grad=rho_grad, &
                             tau=tau_blk)
            else
               call eval_rho(ao, density, rho_blk, ao_grad=ao_grad, rho_grad=rho_grad)
            end if

            allocate (sigma(nb), exc_i(nb), vrho_i(nb), vsigma_i(nb), vsigma(nb))
            vsigma = 0.0_dp
            if (ctx%any_mgga) then
               allocate (lapl(nb), vlapl(nb), vtau_i(nb))
               ! Zero, and read only by functionals that do not need it: the
               ! Laplacian-dependent ones were refused at construction.
               lapl = 0.0_dp
            end if
            do ig = 1, nb
               sigma(ig) = rho_grad(ig, 1)**2 + rho_grad(ig, 2)**2 + rho_grad(ig, 3)**2
               if (present(sigma_out)) sigma_out(g0 + ig - 1) = sigma(ig)
            end do

            do i = 1, ctx%n_func
               exc_i = 0.0_dp
               vrho_i = 0.0_dp
               vsigma_i = 0.0_dp
               ! Same dispatch the Fock build uses. A composition may mix an
               ! LDA component with a GGA one -- a hybrid's correlation part
               ! commonly is -- so the family is per functional rather than per
               ! context, and `any_gga` only says at least one needs sigma.
               select case (ctx%family(i))
               case (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
                  vtau_i = 0.0_dp
                  call xc_f03_mgga_exc_vxc(ctx%func(i), int(nb, 8), rho_blk, sigma, &
                                           lapl, tau_blk, exc_i, vrho_i, vsigma_i, &
                                           vlapl, vtau_i)
                  vsigma = vsigma + ctx%weight(i)*vsigma_i
                  do ig = 1, nb
                     vtau(g0 + ig - 1) = vtau(g0 + ig - 1) + ctx%weight(i)*vtau_i(ig)
                  end do
               case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
                  call xc_f03_gga_exc_vxc(ctx%func(i), int(nb, 8), rho_blk, sigma, &
                                          exc_i, vrho_i, vsigma_i)
                  vsigma = vsigma + ctx%weight(i)*vsigma_i
               case default
                  call xc_f03_lda_exc_vxc(ctx%func(i), int(nb, 8), rho_blk, exc_i, vrho_i)
               end select
               do ig = 1, nb
                  exc(g0 + ig - 1) = exc(g0 + ig - 1) + ctx%weight(i)*exc_i(ig)
                  vrho(g0 + ig - 1) = vrho(g0 + ig - 1) + ctx%weight(i)*vrho_i(ig)
               end do
            end do
            if (ctx%any_mgga) deallocate (lapl, vlapl, vtau_i)

            do ig = 1, nb
               rho(g0 + ig - 1) = rho_blk(ig)
               do id = 1, 3
                  grad_coeff(g0 + ig - 1, id) = 2.0_dp*vsigma(ig)*rho_grad(ig, id)
               end do
            end do
         end if
      end do
#else
      call error%set(ERROR_VALIDATION, "no libxc in this build")
#endif
   end subroutine xc_grid_gga_quantities

   subroutine xc_grid_kernel_quantities(ctx, mol, density, rho, rho_grad, vrho, &
                                        vsigma, frr, frs, fss, error, &
                                        tau, vtau, frt, fst, ftt, &
                                        grrr, grrs, grss, gsss)
      !! First *and* second functional derivatives at every grid point
      !!
      !! What `xc_grid_gga_quantities` is to the exchange-correlation gradient,
      !! this is to the derivative of the exchange-correlation *potential*. The
      !! double hybrid gradient needs
      !!
      !!     d/dR Tr(P V_xc[D])
      !!
      !! for a density `P` that is not the reference -- the Z-vector's -- and
      !! since `V_xc` is itself a first derivative, differentiating it brings
      !! second ones.
      !!
      !! **Returned raw rather than pre-combined**, unlike next door: there are
      !! two densities and four independent combinations of them, and which is
      !! wanted depends on whether the moving object is the reference's density
      !! or the other one.
      !!
      !! Restricted only, and a meta-GGA is served through the optional tau
      !! outputs.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), allocatable, intent(out) :: rho(:)
      real(dp), allocatable, intent(out) :: rho_grad(:, :)   !! (npts, 3)
      real(dp), allocatable, intent(out) :: vrho(:), vsigma(:)
      real(dp), allocatable, intent(out) :: frr(:), frs(:), fss(:)
         !! d2E/drho2, d2E/drho dsigma and d2E/dsigma2
      real(dp), allocatable, intent(out), optional :: tau(:)
      real(dp), allocatable, intent(out), optional :: vtau(:)
      real(dp), allocatable, intent(out), optional :: frt(:), fst(:), ftt(:)
      real(dp), allocatable, intent(out), optional :: grrr(:), grrs(:), grss(:), gsss(:)
         !! Third functional derivatives: `v3rho3`, `v3rho2sigma`, `v3rhosigma2`
         !! and `v3sigma3`, weighted and summed over components exactly as the
         !! second derivatives beside them.
         !!
         !! Asked for only by the double-hybrid Hessian, which is the first
         !! thing here to differentiate the *kernel*. A double-hybrid gradient
         !! stops at `f_xc`; its Hessian does not, because the perturbed
         !! Z-vector's operator is the differentiated kernel and the pass-one
         !! term `Tr(D_rel V_xc^(XY))` carries `g_xc rho^X rho^Y` at fixed
         !! density. Omitting it is a 400 per cent error on that term, not a
         !! correction to it.
         !!
         !! LDA and GGA. A meta-GGA would need the tau channels of the third
         !! derivative too, and no double hybrid shipped here is one.
         !! The kinetic-energy-density channel: `tau` itself, `dE/dtau`, and
         !! the three second derivatives that touch it. Absent for anything but
         !! a meta-GGA, where they are zero and the caller has no use for them.
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :)
      real(dp), allocatable :: rho_blk(:), grad_blk(:, :), sigma(:)
      real(dp), allocatable :: exc_i(:), vrho_i(:), vsigma_i(:)
      real(dp), allocatable :: frr_i(:), frs_i(:), fss_i(:)
      real(dp), allocatable :: vtau_i(:), frt_i(:), fst_i(:), ftt_i(:)
      real(dp), allocatable :: lapl_k(:), lscr(:), tau_k(:)
      ! One buffer per unwanted output. libxc's derivatives are `intent(out)`
      ! dummies through the bind(C) interface, and handing the same array to
      ! four of them at once is aliasing the standard does not allow -- it
      ! happens to work today because nothing wanted shares the buffer, which
      ! is not a property to depend on.
      real(dp), allocatable :: lscr_rl(:), lscr_sl(:), lscr_ll(:), lscr_lt(:)
      logical :: want_tau
      logical :: want_kxc
      real(dp), allocatable :: grrr_i(:), grrs_i(:), grss_i(:), gsss_i(:)
      integer :: g0, g1, nb, i, ig, id, npts, n_kxc

      npts = ctx%grid%n_points
      allocate (rho(npts), rho_grad(npts, 3), vrho(npts), vsigma(npts), &
                frr(npts), frs(npts), fss(npts))
      ! All five tau outputs travel together: a caller that wants one wants the
      ! set, and a meta-GGA is the only thing that makes any of them non-zero.
      want_tau = present(tau) .and. present(vtau) .and. present(frt) &
                 .and. present(fst) .and. present(ftt)
      ! All four together or none: a caller with three of them has a bug, and
      ! silently handing back an unallocated fourth would surface as a segfault
      ! somewhere else entirely.
      !
      ! **Enforced, not merely stated.** Written as an `and` alone, three of
      ! four made `want_kxc` false and every requested output came back
      ! unallocated with no error -- so the failure detached from the call that
      ! caused it and arrived later as a segfault in whatever touched the array.
      ! That is precisely the shape this comment was warning about, left
      ! undefended.
      n_kxc = 0
      if (present(grrr)) n_kxc = n_kxc + 1
      if (present(grrs)) n_kxc = n_kxc + 1
      if (present(grss)) n_kxc = n_kxc + 1
      if (present(gsss)) n_kxc = n_kxc + 1
      if (n_kxc /= 0 .and. n_kxc /= 4) then
         call error%set(ERROR_VALIDATION, "the third functional derivatives are "// &
                        "all four or none: grrr, grrs, grss and gsss are one "// &
                        "quantity in four pieces, and a caller asking for a "// &
                        "subset has a bug rather than a preference.")
         return
      end if
      want_kxc = n_kxc == 4
      ! Refused here rather than left to the caller. The meta-GGA branch below
      ! calls no third-derivative routine, so asking for one and getting a
      ! silent array of zeros is the shape of the failure this file refuses
      ! everywhere else: a converged, plausible, wrong answer with nothing in
      ! the output to say a term was missing. The one consumer today checks
      ! this itself, which is exactly why the guard belongs on this side -- the
      ! second consumer will not know to.
      !
      ! **`any_mgga` rather than the family constants**, which is not a style
      ! choice: `XC_FAMILY_MGGA` and `XC_FAMILY_HYB_MGGA` come from libxc's own
      ! module and are only in scope under `MQC_WITH_LIBXC`, while this routine
      ! is compiled either way. Naming them here broke every build without
      ! libxc. The flag is set from those two families and nothing else, so it
      ! says the same thing and says it in both configurations --
      ! `xc_kernel2_apply` already asks the question this way.
      if (want_kxc .and. ctx%any_mgga) then
         call error%set(ERROR_VALIDATION, "third functional derivatives "// &
                        "are not implemented for meta-GGAs: the tau channels "// &
                        "of the third derivative would be needed and are not "// &
                        "evaluated here.")
         return
      end if
      if (want_kxc) then
         allocate (grrr(npts), grrs(npts), grss(npts), gsss(npts))
         grrr = 0.0_dp
         grrs = 0.0_dp
         grss = 0.0_dp
         gsss = 0.0_dp
      end if
      if (want_tau) then
         allocate (tau(npts), vtau(npts), frt(npts), fst(npts), ftt(npts))
         tau = 0.0_dp
         vtau = 0.0_dp
         frt = 0.0_dp
         fst = 0.0_dp
         ftt = 0.0_dp
      end if
      rho = 0.0_dp
      rho_grad = 0.0_dp
      vrho = 0.0_dp
      vsigma = 0.0_dp
      frr = 0.0_dp
      frs = 0.0_dp
      fss = 0.0_dp

      if (.not. ctx%active) return
      if (.not. xc_available()) then
         call error%set(ERROR_VALIDATION, "no libxc in this build")
         return
      end if
      ! Meta-GGA is served through the optional tau outputs. A caller that does
      ! not ask for them gets the LDA and GGA channels and would be silently
      ! missing tau, so asking is the thing that is checked rather than the
      ! functional family.
      if (ctx%any_mgga .and. .not. want_tau) then
         call error%set(ERROR_VALIDATION, "this is a meta-GGA, so the kinetic-energy "// &
                        "density channel is needed: ask for tau, vtau, frt, fst and ftt. "// &
                        "Refused rather than returned with those terms missing.")
         return
      end if
      if (ctx%polarized) then
         call error%set(ERROR_VALIDATION, "the derivative of the exchange-correlation "// &
                        "potential is implemented for a restricted reference only")
         return
      end if

#ifdef MQC_WITH_LIBXC
      ! `ctx%point_block`, not the module constant: `keywords.dft.block_size`
      ! sets the former and every other loop in this module reads it, so
      ! blocking on the constant here made the keyword a no-op on this path.
      ! These three routines still apply no AO screen, so
      ! `keywords.dft.screening_tolerance` remains inert for them -- adding it
      ! means resizing every downstream array to the significant-AO subset, the
      ! way the potential loops do, which is a change to numerics and belongs on
      ! its own.
      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         ! The AO gradients are needed whenever anything here is a GGA, and
         ! harmless otherwise -- but they are the expensive half of the
         ! evaluation, so a pure LDA does not pay for them.
         if (ctx%any_gga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, grad=ao_grad)
            if (error%has_error()) return
            if (want_tau) then
               call eval_rho(ao, density, rho_blk, ao_grad=ao_grad, &
                             rho_grad=grad_blk, tau=tau_k)
            else
               call eval_rho(ao, density, rho_blk, ao_grad=ao_grad, rho_grad=grad_blk)
            end if
         else
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error)
            if (error%has_error()) return
            call eval_rho(ao, density, rho_blk)
            if (allocated(grad_blk)) deallocate (grad_blk)
            allocate (grad_blk(nb, 3))
            grad_blk = 0.0_dp
         end if

         if (allocated(sigma)) deallocate (sigma, exc_i, vrho_i, vsigma_i, &
                                           frr_i, frs_i, fss_i)
         allocate (sigma(nb), exc_i(nb), vrho_i(nb), vsigma_i(nb), &
                   frr_i(nb), frs_i(nb), fss_i(nb))
         if (want_kxc) then
            if (allocated(grrr_i)) deallocate (grrr_i, grrs_i, grss_i, gsss_i)
            allocate (grrr_i(nb), grrs_i(nb), grss_i(nb), gsss_i(nb))
            grrr_i = 0.0_dp
            grrs_i = 0.0_dp
            grss_i = 0.0_dp
            gsss_i = 0.0_dp
         end if
         if (want_tau) then
            if (allocated(vtau_i)) deallocate (vtau_i, frt_i, fst_i, ftt_i, lapl_k, lscr, &
                                               lscr_rl, lscr_sl, lscr_ll, lscr_lt)
            allocate (vtau_i(nb), frt_i(nb), fst_i(nb), ftt_i(nb), lapl_k(nb), lscr(nb), &
                      lscr_rl(nb), lscr_sl(nb), lscr_ll(nb), lscr_lt(nb))
            lapl_k = 0.0_dp
         end if
         do ig = 1, nb
            sigma(ig) = grad_blk(ig, 1)**2 + grad_blk(ig, 2)**2 + grad_blk(ig, 3)**2
         end do

         do i = 1, ctx%n_func
            select case (ctx%family(i))
            case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
               call xc_f03_gga_exc_vxc(ctx%func(i), int(nb, 8), rho_blk, sigma, &
                                       exc_i, vrho_i, vsigma_i)
               call xc_f03_gga_fxc(ctx%func(i), int(nb, 8), rho_blk, sigma, &
                                   frr_i, frs_i, fss_i)
               if (want_kxc) then
                  call xc_f03_gga_kxc(ctx%func(i), int(nb, 8), rho_blk, sigma, &
                                      grrr_i, grrs_i, grss_i, gsss_i)
               end if
               do ig = 1, nb
                  vsigma(g0 + ig - 1) = vsigma(g0 + ig - 1) + ctx%weight(i)*vsigma_i(ig)
                  frs(g0 + ig - 1) = frs(g0 + ig - 1) + ctx%weight(i)*frs_i(ig)
                  fss(g0 + ig - 1) = fss(g0 + ig - 1) + ctx%weight(i)*fss_i(ig)
               end do
               if (want_kxc) then
                  do ig = 1, nb
                     grrs(g0 + ig - 1) = grrs(g0 + ig - 1) + ctx%weight(i)*grrs_i(ig)
                     grss(g0 + ig - 1) = grss(g0 + ig - 1) + ctx%weight(i)*grss_i(ig)
                     gsss(g0 + ig - 1) = gsss(g0 + ig - 1) + ctx%weight(i)*gsss_i(ig)
                  end do
               end if
            case (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
               call xc_f03_mgga_exc_vxc(ctx%func(i), int(nb, 8), rho_blk, sigma, lapl_k, &
                                        tau_k, exc_i, vrho_i, vsigma_i, lscr, vtau_i)
               call xc_f03_mgga_fxc(ctx%func(i), int(nb, 8), rho_blk, sigma, lapl_k, tau_k, &
                                    frr_i, frs_i, lscr_rl, frt_i, &
                                    fss_i, lscr_sl, fst_i, &
                                    lscr_ll, lscr_lt, ftt_i)
               do ig = 1, nb
                  vsigma(g0 + ig - 1) = vsigma(g0 + ig - 1) + ctx%weight(i)*vsigma_i(ig)
                  frs(g0 + ig - 1) = frs(g0 + ig - 1) + ctx%weight(i)*frs_i(ig)
                  fss(g0 + ig - 1) = fss(g0 + ig - 1) + ctx%weight(i)*fss_i(ig)
                  vtau(g0 + ig - 1) = vtau(g0 + ig - 1) + ctx%weight(i)*vtau_i(ig)
                  frt(g0 + ig - 1) = frt(g0 + ig - 1) + ctx%weight(i)*frt_i(ig)
                  fst(g0 + ig - 1) = fst(g0 + ig - 1) + ctx%weight(i)*fst_i(ig)
                  ftt(g0 + ig - 1) = ftt(g0 + ig - 1) + ctx%weight(i)*ftt_i(ig)
               end do
            case default
               call xc_f03_lda_exc_vxc(ctx%func(i), int(nb, 8), rho_blk, exc_i, vrho_i)
               call xc_f03_lda_fxc(ctx%func(i), int(nb, 8), rho_blk, frr_i)
               if (want_kxc) then
                  call xc_f03_lda_kxc(ctx%func(i), int(nb, 8), rho_blk, grrr_i)
               end if
            end select
            do ig = 1, nb
               vrho(g0 + ig - 1) = vrho(g0 + ig - 1) + ctx%weight(i)*vrho_i(ig)
               frr(g0 + ig - 1) = frr(g0 + ig - 1) + ctx%weight(i)*frr_i(ig)
            end do
            if (want_kxc) then
               do ig = 1, nb
                  grrr(g0 + ig - 1) = grrr(g0 + ig - 1) + ctx%weight(i)*grrr_i(ig)
               end do
            end if
         end do

         do ig = 1, nb
            rho(g0 + ig - 1) = rho_blk(ig)
            do id = 1, 3
               rho_grad(g0 + ig - 1, id) = grad_blk(ig, id)
            end do
            ! The same divergence the kernel floors, and for the same reason:
            ! `v2rho2` goes as rho^(-2/3) at the tail of every atomic grid. Only
            ! the *second* derivatives are floored -- `vrho` and `vsigma` stay
            ! finite there, and they are what the reference's own potential is
            ! built from, so zeroing them would perturb a converged quantity.
            if (want_tau) tau(g0 + ig - 1) = tau_k(ig)
            if (rho_blk(ig) < KERNEL_RHO_FLOOR) then
               frr(g0 + ig - 1) = 0.0_dp
               frs(g0 + ig - 1) = 0.0_dp
               fss(g0 + ig - 1) = 0.0_dp
               if (want_tau) then
                  frt(g0 + ig - 1) = 0.0_dp
                  fst(g0 + ig - 1) = 0.0_dp
                  ftt(g0 + ig - 1) = 0.0_dp
               end if
               ! The third derivatives go the same way, and they have to: they
               ! diverge faster than the second ones do, so a consumer contracting
               ! them at the tail would be handed the largest numbers on the grid
               ! from the region contributing least to anything. Flooring both
               ! orders together is also what lets one be differenced against the
               ! other: across this boundary an unfloored order meets a step.
               if (want_kxc) then
                  grrr(g0 + ig - 1) = 0.0_dp
                  grrs(g0 + ig - 1) = 0.0_dp
                  grss(g0 + ig - 1) = 0.0_dp
                  gsss(g0 + ig - 1) = 0.0_dp
               end if
            end if
         end do
      end do
#endif
   end subroutine xc_grid_kernel_quantities

   subroutine xc_add_potential(ctx, mol, density, v_xc, e_xc, n_elec, error)
      !! The exchange-correlation potential and energy for one density
      !!
      !!     E_xc = sum_g w_g rho_g eps_xc(rho_g)
      !!     V_uv = sum_g w_g v_xc(rho_g) chi_u(r_g) chi_v(r_g)
      !!
      !! `v_xc` comes back on its own rather than added into a Fock matrix,
      !! because the Kohn-Sham energy takes E_xc directly and not half the trace
      !! of D V_xc, so the caller needs the Fock matrix *without* it.
      !!
      !! `n_elec` is the integrated density, returned because it is the cheapest
      !! check that the grid and the density agree.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(out) :: v_xc(:, :)
      real(dp), intent(out) :: e_xc
      real(dp), intent(out) :: n_elec
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: ao(:, :), rho(:), exc(:), vrho(:)
      real(dp), allocatable :: exc_i(:), vrho_i(:), grad_coeff(:, :)
      real(dp), allocatable :: ao_grad(:, :, :), rho_grad(:, :), sigma(:)
      real(dp), allocatable :: vsigma(:), vsigma_i(:)
      real(dp) :: e_nl_total
      real(dp), allocatable :: tau(:), vtau(:), vtau_i(:), lapl(:), vlapl(:)
      real(dp), allocatable :: v_local(:, :), v_sig(:, :), d_sig(:, :)
      real(dp), allocatable :: extents(:)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_list(:), ao_offset(:)
      integer :: n_sig, ia, ja
      real(dp) :: e_local, n_local
      type(error_t) :: local_error
      logical :: failed
      integer :: g0, g1, nb, i, ig, id

      v_xc = 0.0_dp
      e_xc = 0.0_dp
      n_elec = 0.0_dp

      if (.not. ctx%active) return
      if (.not. xc_available()) then
         call error%set(ERROR_VALIDATION, "no libxc in this build")
         return
      end if
      if (ctx%polarized) then
         call error%set(ERROR_VALIDATION, "this exchange-correlation context was built "// &
                        "spin-polarised and cannot be evaluated on a single density; "// &
                        "the unrestricted path wants xc_add_potential_uks")
         return
      end if

#ifdef MQC_WITH_LIBXC
      ! One thread per block of grid points. This loop is where a density
      ! functional run spends nearly all of its time, and the blocks are
      ! independent: each evaluates the basis on its own points, calls libxc on
      ! its own densities, and contributes to `v_xc`, `e_xc` and `n_elec` only by
      ! accumulation. `schedule(dynamic)` because a block's cost depends on how
      ! many basis functions reach its points.
      !
      ! Threaded here rather than by the BLAS underneath, so that the level which
      ! respects `omp_set_num_threads` is this one. The regions inside
      ! `eval_ao_block` nest within this one and so run serially, nesting being
      ! off by default.
      !
      ! `local_error` is **firstprivate and not private**. A private copy of a
      ! derived type is not required to pick up its default initialisation, and
      ! here it did not: every thread started holding an error. Clearing it
      ! inside the region is not the fix either -- `clear` deallocates the
      ! message, and doing that to an uninitialised copy frees a pointer that was
      ! never allocated.
      failed = .false.

      ! One bound per shell for the whole molecule, computed once and read by
      ! every block of every iteration. The tolerance is on the AO *value*, and a
      ! shell dropped beyond its radius contributes less than that at every point
      ! of the block, so this is a truncation with a stated bound.
      call shell_extents(mol, ctx%screen_tol, extents)

      !$omp parallel default(none) &
      !$omp    shared(ctx, mol, density, v_xc, e_xc, n_elec, error, failed, extents) &
      !$omp    private(g0, g1, nb, i, ig, id, ao, rho, exc, vrho, exc_i, vrho_i, &
      !$omp            grad_coeff, ao_grad, rho_grad, sigma, vsigma, vsigma_i, &
      !$omp            tau, vtau, vtau_i, lapl, vlapl) &
      !$omp    private(v_local, e_local, n_local, v_sig, d_sig, shell_mask, &
      !$omp            ao_list, ao_offset, n_sig, ia, ja) &
      !$omp    firstprivate(local_error)
      allocate (v_local(size(v_xc, 1), size(v_xc, 2)))
      v_local = 0.0_dp
      ! Sized to the whole basis once rather than to `n_sig` per block: the
      ! leading sub-block is what gets used, and re-allocating inside the loop
      ! is the allocator traffic this change exists to remove.
      allocate (v_sig(mol%nao, mol%nao), d_sig(mol%nao, mol%nao))
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(mol%nao))
      e_local = 0.0_dp
      n_local = 0.0_dp

      !$omp do schedule(dynamic)
      do g0 = 1, ctx%grid%n_points, ctx%point_block
         ! A thread that has seen a failure stops doing work, but the loop still
         ! has to be run out: leaving an OpenMP region early is not allowed, and
         ! the barrier at its end has to be reached by everybody.
         if (failed) cycle

         g1 = min(g0 + ctx%point_block - 1, ctx%grid%n_points)
         nb = g1 - g0 + 1

         ! Which shells reach this block at all. Everything downstream -- the
         ! basis evaluation, the density contraction, the potential gemm -- runs
         ! on that subset: both gemms are n_points * n_ao^2, so halving n_ao
         ! quarters them.
         call block_significant_aos(mol, ctx%grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig)
         if (n_sig == 0) cycle          ! empty space; no basis function reaches it

         ! The density restricted to the same subset. Gathered per block, which
         ! is n_sig^2 of copying against n_points * n_sig^2 of arithmetic.
         do ja = 1, n_sig
            do ia = 1, n_sig
               d_sig(ia, ja) = density(ao_list(ia), ao_list(ja))
            end do
         end do

         if (ctx%any_mgga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, local_error, &
                               grad=ao_grad, shell_mask=shell_mask, &
                               ao_offset=ao_offset, n_ao_out=n_sig)
            call eval_rho(ao, d_sig(1:n_sig, 1:n_sig), rho, ao_grad=ao_grad, &
                          rho_grad=rho_grad, tau=tau)
         else if (ctx%any_gga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, local_error, &
                               grad=ao_grad, shell_mask=shell_mask, &
                               ao_offset=ao_offset, n_ao_out=n_sig)
            call eval_rho(ao, d_sig(1:n_sig, 1:n_sig), rho, ao_grad=ao_grad, &
                          rho_grad=rho_grad)
         else
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, local_error, &
                               shell_mask=shell_mask, ao_offset=ao_offset, &
                               n_ao_out=n_sig)
            call eval_rho(ao, d_sig(1:n_sig, 1:n_sig), rho)
         end if
         if (local_error%has_error()) then
            !$omp critical (xc_failure)
            if (.not. failed) then
               failed = .true.
               error = local_error
            end if
            !$omp end critical (xc_failure)
            cycle
         end if

         if (allocated(exc)) deallocate (exc, vrho, exc_i, vrho_i)
         allocate (exc(nb), vrho(nb), exc_i(nb), vrho_i(nb))
         exc = 0.0_dp
         vrho = 0.0_dp
         if (ctx%any_gga) then
            if (allocated(sigma)) deallocate (sigma, vsigma, vsigma_i)
            allocate (sigma(nb), vsigma(nb), vsigma_i(nb))
            vsigma = 0.0_dp
            do ig = 1, nb
               sigma(ig) = rho_grad(ig, 1)**2 + rho_grad(ig, 2)**2 + rho_grad(ig, 3)**2
            end do
         end if
         if (ctx%any_mgga) then
            if (allocated(vtau)) deallocate (vtau, vtau_i, lapl, vlapl)
            allocate (vtau(nb), vtau_i(nb), lapl(nb), vlapl(nb))
            vtau = 0.0_dp
            ! Zero, and only ever read by functionals that do not need it -- the
            ! ones that do were refused at construction.
            lapl = 0.0_dp
         end if

         ! Each component contributes its weight of the energy density and of the
         ! potential. One functional with weight one is the ordinary case; more
         ! than one is a composition this repository defines.
         do i = 1, ctx%n_func
            select case (ctx%family(i))
            case (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
               call xc_f03_mgga_exc_vxc(ctx%func(i), int(nb, 8), rho, sigma, lapl, tau, &
                                        exc_i, vrho_i, vsigma_i, vlapl, vtau_i)
               vsigma = vsigma + ctx%weight(i)*vsigma_i
               vtau = vtau + ctx%weight(i)*vtau_i
            case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
               call xc_f03_gga_exc_vxc(ctx%func(i), int(nb, 8), rho, sigma, &
                                       exc_i, vrho_i, vsigma_i)
               vsigma = vsigma + ctx%weight(i)*vsigma_i
            case default
               call xc_f03_lda_exc_vxc(ctx%func(i), int(nb, 8), rho, exc_i, vrho_i)
            end select
            exc = exc + ctx%weight(i)*exc_i
            vrho = vrho + ctx%weight(i)*vrho_i
         end do

         n_local = n_local + sum(ctx%grid%weights(g0:g1)*rho)
         e_local = e_local + sum(ctx%grid%weights(g0:g1)*rho*exc)

         ! The gradient coefficient is dE/d(grad rho). Differentiating
         ! sigma = |grad rho|^2 gives 2 vsigma grad rho; the unrestricted case
         ! adds a cross-spin term, which is why this is assembled by the caller
         ! rather than inside `accumulate_xc_matrix`.
         if (ctx%any_gga) then
            if (allocated(grad_coeff)) deallocate (grad_coeff)
            allocate (grad_coeff(nb, 3))
            do id = 1, 3
               do ig = 1, nb
                  grad_coeff(ig, id) = 2.0_dp*vsigma(ig)*rho_grad(ig, id)
               end do
            end do
         end if

         v_sig(1:n_sig, 1:n_sig) = 0.0_dp
         call accumulate_xc_matrix(ctx%grid%weights(g0:g1), ao, vrho, &
                                   v_sig(1:n_sig, 1:n_sig), &
                                   ao_grad=ao_grad, grad_coeff=grad_coeff, &
                                   vtau=vtau, any_gga=ctx%any_gga, &
                                   any_mgga=ctx%any_mgga)
         ! Back into the full matrix. The scatter is n_sig^2 per block against
         ! the gemm's n_points * n_sig^2, so it disappears against what it
         ! bought.
         do ja = 1, n_sig
            do ia = 1, n_sig
               v_local(ao_list(ia), ao_list(ja)) = &
                  v_local(ao_list(ia), ao_list(ja)) + v_sig(ia, ja)
            end do
         end do
      end do
      !$omp end do

      ! One reduction per thread rather than per block. The matrix is n_ao
      ! square -- half a megabyte here -- so the copies cost far less than
      ! contending for the shared one a hundred times each.
      !$omp critical (xc_reduce)
      v_xc = v_xc + v_local
      e_xc = e_xc + e_local
      n_elec = n_elec + n_local
      !$omp end critical (xc_reduce)
      !$omp end parallel

      ! VV10 last, on its own grid, adding into the same matrix and the same
      ! energy. Outside the parallel region because it runs its own two sweeps
      ! over a different quadrature.
      if (ctx%nlc_b /= 0.0_dp .or. ctx%nlc_c /= 0.0_dp) then
         call vv10_add_potential(ctx, mol, density, v_xc, e_nl_total, error)
         if (error%has_error()) return
         e_xc = e_xc + e_nl_total
      end if
#endif
   end subroutine xc_add_potential

   subroutine xc_add_potential_uks(ctx, mol, d_alpha, d_beta, v_alpha, v_beta, &
                                   e_xc, n_elec, error)
      !! The exchange-correlation potential and energy for a pair of spin densities
      !!
      !! The unrestricted counterpart of `xc_add_potential`, with the same
      !! structure. Two things differ, and they are the two places an
      !! unrestricted functional goes wrong.
      !!
      !! **libxc's arrays are spin-interleaved.** A polarised functional takes
      !! `rho` as (rho_a, rho_b) per point, `sigma` as (sigma_aa, sigma_ab,
      !! sigma_bb) and `tau` as (tau_a, tau_b), spin fastest, and returns `vrho`
      !! and `vtau` in the same layout with `vsigma` in three components. The F03
      !! bindings take bare arrays, so a wrong stride here is a silent misread.
      !!
      !! **The gradient term couples the spins.** sigma_ab = grad rho_a . grad rho_b
      !! belongs to both, so
      !!
      !!     dE/dgrad rho_a = 2 vsigma_aa grad rho_a + vsigma_ab grad rho_b
      !!
      !! and dropping that cross term leaves a GGA that converges and is wrong.
      !!
      !! The spin densities are the true ones, not doubled: `d_alpha` holds
      !! C_a C_a^T, so rho_a integrates to the number of alpha electrons and the two
      !! together to the total. `n_elec` is that total.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: d_alpha(:, :), d_beta(:, :)
      real(dp), intent(out) :: v_alpha(:, :), v_beta(:, :)
      real(dp), intent(out) :: e_xc
      real(dp), intent(out) :: n_elec
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :)
      real(dp) :: e_nl_total
      real(dp), allocatable :: v_nl(:, :)
         !! VV10's potential, built once and added to both spins
      real(dp), allocatable :: rho_a(:), rho_b(:), grad_a(:, :), grad_b(:, :)
      real(dp), allocatable :: tau_a(:), tau_b(:)
      real(dp), allocatable :: rho(:), sigma(:), tau(:), lapl(:)
      real(dp), allocatable :: exc(:), vrho(:), vsigma(:), vtau(:)
      real(dp), allocatable :: exc_i(:), vrho_i(:), vsigma_i(:), vtau_i(:), vlapl(:)
      real(dp), allocatable :: vrho_s(:), vtau_s(:), grad_coeff(:, :)
      real(dp), allocatable :: va_local(:, :), vb_local(:, :)
      real(dp), allocatable :: va_sig(:, :), vb_sig(:, :), da_sig(:, :), db_sig(:, :)
      real(dp), allocatable :: extents(:)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_list(:), ao_offset(:)
      integer :: n_sig, ia, ja
      real(dp) :: e_local, n_local
      type(error_t) :: local_error
      logical :: failed
      integer :: g0, g1, nb, i, ig, id

      v_alpha = 0.0_dp
      v_beta = 0.0_dp
      e_xc = 0.0_dp
      n_elec = 0.0_dp

      if (.not. ctx%active) return
      if (.not. xc_available()) then
         call error%set(ERROR_VALIDATION, "no libxc in this build")
         return
      end if
      if (.not. ctx%polarized) then
         call error%set(ERROR_VALIDATION, "this exchange-correlation context was built "// &
                        "spin-restricted and would return a restricted potential for "// &
                        "two spin densities; build it with polarized=.true.")
         return
      end if

#ifdef MQC_WITH_LIBXC
      ! Threaded over blocks exactly as the restricted path is.
      ! `firstprivate(local_error)` rather than `private` -- a private copy of a
      ! derived type need not pick up its default initialisation, and here it
      ! does not.
      failed = .false.

      ! One bound per shell for the whole molecule, computed once and read by
      ! every block of every iteration.
      call shell_extents(mol, ctx%screen_tol, extents)

      !$omp parallel default(none) &
      !$omp    shared(ctx, mol, d_alpha, d_beta, v_alpha, v_beta, e_xc, n_elec, &
      !$omp           error, failed, extents) &
      !$omp    private(g0, g1, nb, i, ig, id, ao, ao_grad, rho_a, rho_b, grad_a, &
      !$omp            grad_b, tau_a, tau_b, rho, sigma, tau, lapl, exc, vrho, &
      !$omp            vsigma, vtau, exc_i, vrho_i, vsigma_i, vtau_i, vlapl, &
      !$omp            vrho_s, vtau_s, grad_coeff) &
      !$omp    private(va_local, vb_local, e_local, n_local) &
      !$omp    private(va_sig, vb_sig, da_sig, db_sig, shell_mask, ao_list, &
      !$omp            ao_offset, n_sig, ia, ja) &
      !$omp    firstprivate(local_error)
      allocate (va_local(size(v_alpha, 1), size(v_alpha, 2)))
      allocate (vb_local(size(v_beta, 1), size(v_beta, 2)))
      va_local = 0.0_dp
      vb_local = 0.0_dp
      ! Sized to the whole basis once rather than to `n_sig` per block: the
      ! leading sub-block is what gets used, and re-allocating inside the loop is
      ! the allocator traffic this exists to remove.
      allocate (va_sig(mol%nao, mol%nao), vb_sig(mol%nao, mol%nao))
      allocate (da_sig(mol%nao, mol%nao), db_sig(mol%nao, mol%nao))
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(mol%nao))
      e_local = 0.0_dp
      n_local = 0.0_dp

      !$omp do schedule(dynamic)
      do g0 = 1, ctx%grid%n_points, ctx%point_block
         ! A thread that has seen a failure stops working, but the loop still has
         ! to run out: leaving an OpenMP region early is not allowed and the
         ! barrier at its end has to be reached by every thread.
         if (failed) cycle

         g1 = min(g0 + ctx%point_block - 1, ctx%grid%n_points)
         nb = g1 - g0 + 1

         ! Which shells reach this block at all. Both spins share the answer --
         ! the test is on the basis and the geometry, and knows nothing about
         ! which density is being contracted -- so one screen serves both, and
         ! the two spin matrices stay on the same index set as each other.
         call block_significant_aos(mol, ctx%grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig)
         if (n_sig == 0) cycle          ! empty space; no basis function reaches it

         do ja = 1, n_sig
            do ia = 1, n_sig
               da_sig(ia, ja) = d_alpha(ao_list(ia), ao_list(ja))
               db_sig(ia, ja) = d_beta(ao_list(ia), ao_list(ja))
            end do
         end do

         ! One AO evaluation for both spins -- the expensive part -- then the
         ! density contraction once per spin.
         !
         ! The evaluation is hoisted out of the spin-family branch so its error is
         ! checked once, before anything reads `ao`. That ordering is load-bearing:
         ! `eval_ao_block` returns without allocating `ao` on its error paths, so
         ! calling `eval_rho` first would hand it an unallocated array.
         if (ctx%any_gga .or. ctx%any_mgga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, local_error, &
                               grad=ao_grad, shell_mask=shell_mask, &
                               ao_offset=ao_offset, n_ao_out=n_sig)
         else
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, local_error, &
                               shell_mask=shell_mask, ao_offset=ao_offset, &
                               n_ao_out=n_sig)
         end if
         if (local_error%has_error()) then
            !$omp critical (xc_uks_failure)
            if (.not. failed) then
               failed = .true.
               error = local_error
            end if
            !$omp end critical (xc_uks_failure)
            cycle
         end if

         if (ctx%any_mgga) then
            call eval_rho(ao, da_sig(1:n_sig, 1:n_sig), rho_a, ao_grad=ao_grad, &
                          rho_grad=grad_a, tau=tau_a)
            call eval_rho(ao, db_sig(1:n_sig, 1:n_sig), rho_b, ao_grad=ao_grad, &
                          rho_grad=grad_b, tau=tau_b)
         else if (ctx%any_gga) then
            call eval_rho(ao, da_sig(1:n_sig, 1:n_sig), rho_a, ao_grad=ao_grad, &
                          rho_grad=grad_a)
            call eval_rho(ao, db_sig(1:n_sig, 1:n_sig), rho_b, ao_grad=ao_grad, &
                          rho_grad=grad_b)
         else
            call eval_rho(ao, da_sig(1:n_sig, 1:n_sig), rho_a)
            call eval_rho(ao, db_sig(1:n_sig, 1:n_sig), rho_b)
         end if

         if (allocated(rho)) deallocate (rho, exc, vrho, exc_i, vrho_i, vrho_s)
         allocate (rho(2*nb), exc(nb), vrho(2*nb), exc_i(nb), vrho_i(2*nb), vrho_s(nb))
         exc = 0.0_dp
         vrho = 0.0_dp
         do ig = 1, nb
            rho(2*ig - 1) = rho_a(ig)
            rho(2*ig) = rho_b(ig)
         end do

         if (ctx%any_gga) then
            if (allocated(sigma)) deallocate (sigma, vsigma, vsigma_i)
            allocate (sigma(3*nb), vsigma(3*nb), vsigma_i(3*nb))
            vsigma = 0.0_dp
            do ig = 1, nb
               sigma(3*ig - 2) = dot_product(grad_a(ig, :), grad_a(ig, :))
               sigma(3*ig - 1) = dot_product(grad_a(ig, :), grad_b(ig, :))
               sigma(3*ig) = dot_product(grad_b(ig, :), grad_b(ig, :))
            end do
         end if
         if (ctx%any_mgga) then
            if (allocated(tau)) deallocate (tau, vtau, vtau_i, vtau_s, lapl, vlapl)
            allocate (tau(2*nb), vtau(2*nb), vtau_i(2*nb), vtau_s(nb), &
                      lapl(2*nb), vlapl(2*nb))
            vtau = 0.0_dp
            ! Zero, and only read by functionals that do not need it -- the ones
            ! that do were refused at construction.
            lapl = 0.0_dp
            do ig = 1, nb
               tau(2*ig - 1) = tau_a(ig)
               tau(2*ig) = tau_b(ig)
            end do
         end if

         do i = 1, ctx%n_func
            select case (ctx%family(i))
            case (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
               call xc_f03_mgga_exc_vxc(ctx%func(i), int(nb, 8), rho, sigma, lapl, tau, &
                                        exc_i, vrho_i, vsigma_i, vlapl, vtau_i)
               vsigma = vsigma + ctx%weight(i)*vsigma_i
               vtau = vtau + ctx%weight(i)*vtau_i
            case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
               call xc_f03_gga_exc_vxc(ctx%func(i), int(nb, 8), rho, sigma, &
                                       exc_i, vrho_i, vsigma_i)
               vsigma = vsigma + ctx%weight(i)*vsigma_i
            case default
               call xc_f03_lda_exc_vxc(ctx%func(i), int(nb, 8), rho, exc_i, vrho_i)
            end select
            exc = exc + ctx%weight(i)*exc_i
            vrho = vrho + ctx%weight(i)*vrho_i
         end do

         ! exc is per particle, so the energy density is the *total* density times
         ! it -- one of the two places a polarised functional is silently halved.
         n_local = n_local + sum(ctx%grid%weights(g0:g1)*(rho_a + rho_b))
         e_local = e_local + sum(ctx%grid%weights(g0:g1)*(rho_a + rho_b)*exc)

         ! Alpha, then beta: the same assembly with that spin's derivatives, and
         ! the cross-spin gradient term pointing at the other spin's gradient.
         if (allocated(grad_coeff)) deallocate (grad_coeff)
         if (ctx%any_gga) allocate (grad_coeff(nb, 3))

         do ig = 1, nb
            vrho_s(ig) = vrho(2*ig - 1)
         end do
         if (ctx%any_gga) then
            do id = 1, 3
               do ig = 1, nb
                  grad_coeff(ig, id) = 2.0_dp*vsigma(3*ig - 2)*grad_a(ig, id) &
                                       + vsigma(3*ig - 1)*grad_b(ig, id)
               end do
            end do
         end if
         if (ctx%any_mgga) then
            do ig = 1, nb
               vtau_s(ig) = vtau(2*ig - 1)
            end do
         end if
         va_sig(1:n_sig, 1:n_sig) = 0.0_dp
         call accumulate_xc_matrix(ctx%grid%weights(g0:g1), ao, vrho_s, &
                                   va_sig(1:n_sig, 1:n_sig), &
                                   ao_grad=ao_grad, grad_coeff=grad_coeff, &
                                   vtau=vtau_s, any_gga=ctx%any_gga, &
                                   any_mgga=ctx%any_mgga)
         do ja = 1, n_sig
            do ia = 1, n_sig
               va_local(ao_list(ia), ao_list(ja)) = &
                  va_local(ao_list(ia), ao_list(ja)) + va_sig(ia, ja)
            end do
         end do

         do ig = 1, nb
            vrho_s(ig) = vrho(2*ig)
         end do
         if (ctx%any_gga) then
            do id = 1, 3
               do ig = 1, nb
                  grad_coeff(ig, id) = 2.0_dp*vsigma(3*ig)*grad_b(ig, id) &
                                       + vsigma(3*ig - 1)*grad_a(ig, id)
               end do
            end do
         end if
         if (ctx%any_mgga) then
            do ig = 1, nb
               vtau_s(ig) = vtau(2*ig)
            end do
         end if
         vb_sig(1:n_sig, 1:n_sig) = 0.0_dp
         call accumulate_xc_matrix(ctx%grid%weights(g0:g1), ao, vrho_s, &
                                   vb_sig(1:n_sig, 1:n_sig), &
                                   ao_grad=ao_grad, grad_coeff=grad_coeff, &
                                   vtau=vtau_s, any_gga=ctx%any_gga, &
                                   any_mgga=ctx%any_mgga)
         do ja = 1, n_sig
            do ia = 1, n_sig
               vb_local(ao_list(ia), ao_list(ja)) = &
                  vb_local(ao_list(ia), ao_list(ja)) + vb_sig(ia, ja)
            end do
         end do
      end do
      !$omp end do

      ! One reduction per thread rather than per block, as on the restricted
      ! side: two n_ao-square matrices copied once each beats contending for the
      ! shared pair on every block.
      !$omp critical (xc_uks_reduce)
      v_alpha = v_alpha + va_local
      v_beta = v_beta + vb_local
      e_xc = e_xc + e_local
      n_elec = n_elec + n_local
      !$omp end critical (xc_uks_reduce)
      !$omp end parallel

      ! VV10 sees only the total density, so one evaluation serves both spins and
      ! its contribution is identical in each: the kernel is built from rho and
      ! |grad rho| of the total density with no spin dependence anywhere in it,
      ! so dE/drho_a = dE/drho_b and likewise for the gradient term.
      if (ctx%nlc_b /= 0.0_dp .or. ctx%nlc_c /= 0.0_dp) then
         ! Evaluated once and added to both, rather than called twice: with no
         ! spin dependence the two calls did identical work, and the non-local
         ! term is the dominant cost of a `-V` SCF iteration.
         if (.not. allocated(v_nl)) allocate (v_nl(size(v_alpha, 1), size(v_alpha, 2)))
         v_nl = 0.0_dp
         call vv10_add_potential(ctx, mol, d_alpha + d_beta, v_nl, e_nl_total, error)
         if (error%has_error()) return
         v_alpha = v_alpha + v_nl
         v_beta = v_beta + v_nl
         ! Added once, not twice: `e_nl` is the energy of the whole density.
         e_xc = e_xc + e_nl_total
      end if
#endif
   end subroutine xc_add_potential_uks

   subroutine xc_kernel_apply(ctx, mol, density, dtilde, v_kernel, error)
      !! The exchange-correlation kernel applied to a response density
      !!
      !! `f_xc` is the second functional derivative, and this returns
      !!
      !!     V_uv = int w(r) chi_u(r) chi_v(r) f_xc(r) drho(r)
      !!
      !! with `drho` the density change the trial rotation makes. It is what
      !! turns a coupled-perturbed *Hartree-Fock* operator into a
      !! coupled-perturbed *Kohn-Sham* one, and without it a response over a
      !! Kohn-Sham reference is missing a term of the same order as the
      !! exchange it does include.
      !!
      !! **An addend, not a fourth branch.** The two-electron part of the
      !! response operator comes three ways -- stored, integral-direct and
      !! density-fitted -- and the kernel is orthogonal to that choice.
      !!
      !! **LDA and GGA; meta-GGA is refused.** For a GGA the response density
      !! has a gradient too, and both of the potential's pieces respond:
      !!
      !!     dv_rho  = f_rr drho + f_rs dsigma
      !!     dv_grad = 2 (f_rs drho + f_ss dsigma) grad rho + 2 v_sigma grad drho
      !!     dsigma  = 2 grad rho . grad drho
      !!
      !! which is the same pair of coefficients `accumulate_xc_matrix` already
      !! turns into a matrix for the potential -- so the kernel is a different
      !! `vrho` and `grad_coeff`, not a different assembly. The last term has no
      !! analogue in the LDA case: `v_sigma` is a *first* derivative, and it
      !! enters because the response density's gradient multiplies it.
      !!
      !! A meta-GGA would add `v2tau2`, `v2rhotau` and `v2sigmatau` and a tau
      !! component of the response density; returning a GGA kernel for one would
      !! be a converged, plausible, wrong response.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)   !! The converged SCF density
      real(dp), intent(in) :: dtilde(:, :)    !! The response density
      real(dp), intent(inout) :: v_kernel(:, :)  !! Accumulated into
      type(error_t), intent(inout) :: error

#ifdef MQC_WITH_LIBXC
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :)
      real(dp), allocatable :: rho(:), rho_grad(:, :), drho(:), drho_grad(:, :)
      real(dp), allocatable :: sigma(:), dsigma(:)
      real(dp), allocatable :: frr(:), frs(:), fss(:), vsig(:)
      real(dp), allocatable :: frr_i(:), frs_i(:), fss_i(:)
      real(dp), allocatable :: exc_i(:), vrho_i(:), vsigma_i(:)
      real(dp), allocatable :: c_rho(:), c_grad(:, :), c_tau(:), no_tau(:)
      real(dp), allocatable :: tau(:), dtau(:), lapl(:)
      real(dp), allocatable :: frt(:), fst(:), ftt(:), vtau(:)
      real(dp), allocatable :: frt_i(:), fst_i(:), ftt_i(:), vtau_i(:)
      real(dp), allocatable :: lapl_scratch(:)
         !! libxc's meta-GGA entry points take the Laplacian and return its
         !! derivatives whether or not the functional uses them. Laplacian
         !! dependent functionals are refused at construction, so these are zeros
         !! in and discarded out.
      real(dp), allocatable :: v_local(:, :), v_sig(:, :), d_sig(:, :), dt_sig(:, :)
      real(dp), allocatable :: extents(:)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_list(:), ao_offset(:)
      integer :: n_sig, ia, ja
      type(error_t) :: local_error
      logical :: failed
      integer :: g0, g1, nb, i, ig, id, npts
      logical :: gga, mgga

      if (.not. ctx%active) return
      if (ctx%polarized) then
         call error%set(ERROR_VALIDATION, "the exchange-correlation kernel is "// &
                        "implemented for a restricted reference only")
         return
      end if

      gga = ctx%any_gga
      mgga = ctx%any_mgga
      ! A meta-GGA needs the AO gradients for the same reason a GGA does, and
      ! then some: tau is built from them too.
      gga = gga .or. mgga
      npts = ctx%grid%n_points
      failed = .false.

      call shell_extents(mol, ctx%screen_tol, extents)

      ! Threaded over blocks like the potential paths. Driven once per CPHF
      ! iteration rather than once per SCF iteration, so a response property or a
      ! Z-vector pays for it repeatedly.
      !$omp parallel default(none) &
      !$omp    shared(ctx, mol, density, dtilde, v_kernel, error, failed, &
      !$omp           gga, mgga, npts, extents) &
      !$omp    private(g0, g1, nb, i, ig, id, ao, ao_grad, rho, rho_grad, drho, &
      !$omp            drho_grad, sigma, dsigma, frr, frs, fss, vsig, frr_i, &
      !$omp            frs_i, fss_i, exc_i, vrho_i, vsigma_i, c_rho, c_grad, &
      !$omp            c_tau, no_tau, tau, dtau, lapl, frt, fst, ftt, vtau, &
      !$omp            frt_i, fst_i, ftt_i, vtau_i, lapl_scratch, v_local) &
      !$omp    private(v_sig, d_sig, dt_sig, shell_mask, ao_list, ao_offset, &
      !$omp            n_sig, ia, ja) &
      !$omp    firstprivate(local_error)
      allocate (v_local(size(v_kernel, 1), size(v_kernel, 2)))
      v_local = 0.0_dp
      allocate (v_sig(mol%nao, mol%nao), d_sig(mol%nao, mol%nao), &
                dt_sig(mol%nao, mol%nao))
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(mol%nao))

      !$omp do schedule(dynamic)
      do g0 = 1, npts, ctx%point_block
         ! As in the potential: a thread that has failed stops working, but the
         ! loop still has to run out so every thread reaches the barrier.
         if (failed) cycle

         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         ! One screen for both densities. The reference and the response density
         ! are contracted against the same basis at the same points, so they
         ! share the kept set and the result scatters back through one `ao_list`.
         call block_significant_aos(mol, ctx%grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig)
         if (n_sig == 0) cycle          ! empty space; no basis function reaches it

         do ja = 1, n_sig
            do ia = 1, n_sig
               d_sig(ia, ja) = density(ao_list(ia), ao_list(ja))
               dt_sig(ia, ja) = dtilde(ao_list(ia), ao_list(ja))
            end do
         end do

         ! The reference density is what `f_xc` is evaluated at; the response
         ! density is what it multiplies. Both are ordinary densities on the
         ! grid, so one routine builds them -- and for a GGA both need their
         ! gradients, which is the only reason the AO gradients are asked for.
         ! Hoisted out of the branch so the error is checked once, before
         ! anything reads `ao` -- `eval_ao_block` returns without allocating it
         ! on its error paths.
         if (gga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, local_error, &
                               grad=ao_grad, shell_mask=shell_mask, &
                               ao_offset=ao_offset, n_ao_out=n_sig)
         else
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, local_error, &
                               shell_mask=shell_mask, ao_offset=ao_offset, &
                               n_ao_out=n_sig)
         end if
         if (local_error%has_error()) then
            !$omp critical (xc_kernel_failure)
            if (.not. failed) then
               failed = .true.
               error = local_error
            end if
            !$omp end critical (xc_kernel_failure)
            cycle
         end if

         if (mgga) then
            ! **The tau convention never has to be decided here.** Whatever
            ! `eval_rho` means by tau is what the energy path fed libxc and what
            ! `accumulate_xc_matrix` differentiates, and the response density
            ! goes through the same routine with `dtilde` in place of `density`.
            call eval_rho(ao, d_sig(1:n_sig, 1:n_sig), rho, ao_grad=ao_grad, &
                          rho_grad=rho_grad, tau=tau)
            call eval_rho(ao, dt_sig(1:n_sig, 1:n_sig), drho, ao_grad=ao_grad, &
                          rho_grad=drho_grad, tau=dtau)
         else if (gga) then
            call eval_rho(ao, d_sig(1:n_sig, 1:n_sig), rho, ao_grad=ao_grad, &
                          rho_grad=rho_grad)
            call eval_rho(ao, dt_sig(1:n_sig, 1:n_sig), drho, ao_grad=ao_grad, &
                          rho_grad=drho_grad)
         else
            call eval_rho(ao, d_sig(1:n_sig, 1:n_sig), rho)
            call eval_rho(ao, dt_sig(1:n_sig, 1:n_sig), drho)
         end if

         if (allocated(frr)) deallocate (frr, frs, fss, vsig, frr_i, frs_i, fss_i, &
                                         exc_i, vrho_i, vsigma_i, sigma, dsigma, c_rho)
         allocate (frr(nb), frs(nb), fss(nb), vsig(nb), frr_i(nb), frs_i(nb), fss_i(nb), &
                   exc_i(nb), vrho_i(nb), vsigma_i(nb), sigma(nb), dsigma(nb), c_rho(nb))
         frr = 0.0_dp
         frs = 0.0_dp
         fss = 0.0_dp
         vsig = 0.0_dp
         if (allocated(frt)) deallocate (frt, fst, ftt, vtau, frt_i, fst_i, ftt_i, &
                                         vtau_i, lapl, lapl_scratch, c_tau)
         allocate (frt(nb), fst(nb), ftt(nb), vtau(nb), frt_i(nb), fst_i(nb), ftt_i(nb), &
                   vtau_i(nb), lapl(nb), lapl_scratch(nb), c_tau(nb))
         frt = 0.0_dp
         fst = 0.0_dp
         ftt = 0.0_dp
         vtau = 0.0_dp
         lapl = 0.0_dp
         c_tau = 0.0_dp
         if (.not. mgga) then
            ! Same reason `sigma` is zeroed on the LDA path: the coefficients
            ! that multiply these are zero, and zero times uninitialised is a
            ! NaN rather than nothing.
            if (.not. allocated(tau)) allocate (tau(nb))
            if (.not. allocated(dtau)) allocate (dtau(nb))
            tau = 0.0_dp
            dtau = 0.0_dp
         end if
         ! Zero rather than left alone on the LDA path: their coefficients are
         ! zero there, and `0 * uninitialised` is a NaN rather than nothing.
         sigma = 0.0_dp
         dsigma = 0.0_dp
         if (gga) then
            do ig = 1, nb
               sigma(ig) = rho_grad(ig, 1)**2 + rho_grad(ig, 2)**2 + rho_grad(ig, 3)**2
               dsigma(ig) = 2.0_dp*(rho_grad(ig, 1)*drho_grad(ig, 1) &
                                    + rho_grad(ig, 2)*drho_grad(ig, 2) &
                                    + rho_grad(ig, 3)*drho_grad(ig, 3))
            end do
         end if

         ! Per component, as everywhere else here: a composition may put an LDA
         ! correlation beside a GGA exchange, and `any_gga` only says that at
         ! least one of them needs sigma.
         do i = 1, ctx%n_func
            select case (ctx%family(i))
            case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
               call xc_f03_gga_fxc(ctx%func(i), int(nb, 8), rho, sigma, &
                                   frr_i, frs_i, fss_i)
               ! `v_sigma` is a first derivative and comes from the ordinary
               ! evaluator. It belongs here because the *response* density's
               ! gradient multiplies it -- the one kernel term that is not a
               ! second derivative.
               call xc_f03_gga_exc_vxc(ctx%func(i), int(nb, 8), rho, sigma, &
                                       exc_i, vrho_i, vsigma_i)
               frr = frr + ctx%weight(i)*frr_i
               frs = frs + ctx%weight(i)*frs_i
               fss = fss + ctx%weight(i)*fss_i
               vsig = vsig + ctx%weight(i)*vsigma_i
            case (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
               ! Six of the ten second derivatives are wanted; the four
               ! Laplacian ones are written into scratch and dropped.
               call xc_f03_mgga_fxc(ctx%func(i), int(nb, 8), rho, sigma, lapl, tau, &
                                    frr_i, frs_i, lapl_scratch, frt_i, &
                                    fss_i, lapl_scratch, fst_i, &
                                    lapl_scratch, lapl_scratch, ftt_i)
               ! `v_sigma` and `v_tau` are *first* derivatives and belong here
               ! for the same reason on this rung as on the last: the response
               ! density's gradient multiplies one and its tau the other, so
               ! neither is a second-derivative term but both are kernel terms.
               call xc_f03_mgga_exc_vxc(ctx%func(i), int(nb, 8), rho, sigma, lapl, tau, &
                                        exc_i, vrho_i, vsigma_i, lapl_scratch, vtau_i)
               frr = frr + ctx%weight(i)*frr_i
               frs = frs + ctx%weight(i)*frs_i
               fss = fss + ctx%weight(i)*fss_i
               frt = frt + ctx%weight(i)*frt_i
               fst = fst + ctx%weight(i)*fst_i
               ftt = ftt + ctx%weight(i)*ftt_i
               vsig = vsig + ctx%weight(i)*vsigma_i
               vtau = vtau + ctx%weight(i)*vtau_i
            case default
               call xc_f03_lda_fxc(ctx%func(i), int(nb, 8), rho, frr_i)
               frr = frr + ctx%weight(i)*frr_i
            end select
         end do

         do ig = 1, nb
            if (rho(ig) < KERNEL_RHO_FLOOR) then
               frr(ig) = 0.0_dp
               frs(ig) = 0.0_dp
               fss(ig) = 0.0_dp
               vsig(ig) = 0.0_dp
               frt(ig) = 0.0_dp
               fst(ig) = 0.0_dp
               ftt(ig) = 0.0_dp
               vtau(ig) = 0.0_dp
            end if
         end do

         ! The potential has three pieces on this rung and every one of them
         ! responds to every component of the response density:
         !
         !     dv_rho   = f_rr drho + f_rs dsigma + f_rt dtau
         !     dv_sigma = f_rs drho + f_ss dsigma + f_st dtau
         !     dv_tau   = f_rt drho + f_st dsigma + f_tt dtau
         !
         ! `dtau` is zero on the LDA and GGA paths, so the same three lines
         ! reduce to what they were rather than branching.
         do ig = 1, nb
            c_rho(ig) = frr(ig)*drho(ig) + frs(ig)*dsigma(ig) + frt(ig)*dtau(ig)
         end do
         if (gga) then
            if (allocated(c_grad)) deallocate (c_grad)
            allocate (c_grad(nb, 3))
            do id = 1, 3
               do ig = 1, nb
                  c_grad(ig, id) = 2.0_dp*(frs(ig)*drho(ig) + fss(ig)*dsigma(ig) &
                                           + fst(ig)*dtau(ig))*rho_grad(ig, id) &
                                   + 2.0_dp*vsig(ig)*drho_grad(ig, id)
               end do
            end do
         end if
         if (mgga) then
            do ig = 1, nb
               c_tau(ig) = frt(ig)*drho(ig) + fst(ig)*dsigma(ig) + ftt(ig)*dtau(ig)
            end do
         end if

         ! The same assembly the potential uses, with the kernel's coefficients
         ! in place of the potential's.
         v_sig(1:n_sig, 1:n_sig) = 0.0_dp
         if (mgga) then
            call accumulate_xc_matrix(ctx%grid%weights(g0:g1), ao, c_rho, &
                                      v_sig(1:n_sig, 1:n_sig), &
                                      ao_grad=ao_grad, grad_coeff=c_grad, vtau=c_tau, &
                                      any_gga=gga, any_mgga=.true.)
         else
            call accumulate_xc_matrix(ctx%grid%weights(g0:g1), ao, c_rho, &
                                      v_sig(1:n_sig, 1:n_sig), &
                                      ao_grad=ao_grad, grad_coeff=c_grad, vtau=no_tau, &
                                      any_gga=gga, any_mgga=.false.)
         end if
         do ja = 1, n_sig
            do ia = 1, n_sig
               v_local(ao_list(ia), ao_list(ja)) = &
                  v_local(ao_list(ia), ao_list(ja)) + v_sig(ia, ja)
            end do
         end do
      end do
      !$omp end do

      ! `v_kernel` is accumulated into rather than assigned -- the caller adds
      ! this to a response operator it has already partly built -- so the
      ! reduction adds too.
      !$omp critical (xc_kernel_reduce)
      v_kernel = v_kernel + v_local
      !$omp end critical (xc_kernel_reduce)
      !$omp end parallel
#else
      call error%set(ERROR_VALIDATION, "no libxc in this build")
      if (size(density) < 0 .or. size(dtilde) < 0 .or. size(v_kernel) < 0) return
      if (mol%nao < 0) return
      if (ctx%n_func < 0) return
#endif
   end subroutine xc_kernel_apply

   subroutine xc_kernel2_apply(ctx, mol, density, dtilde_a, dtilde_b, v_kernel, error)
      !! The third functional derivative contracted against two response densities
      !!
      !! `xc_kernel_apply` one rung up: where that routine is the second
      !! derivative against one trial density, this is the third against two,
      !!
      !!     V_uv = int w(r) chi_u(r) chi_v(r) g_xc(r) drho_a(r) drho_b(r)
      !!
      !! with the sigma channels expanded alongside the rho ones exactly as the
      !! kernel expands `f_xc`. Equivalently: the derivative of
      !! `xc_kernel_apply(density, dtilde_b)` when `density` moves by
      !! `dtilde_a` at fixed geometry, which is what the test differences.
      !! The double-hybrid Hessian needs it because the perturbed Z-vector's
      !! operator differentiates the *kernel*, and the kernel's coefficients
      !! -- `f_xc` and `v_sigma` -- are themselves functionals of the density.
      !!
      !! Writing `u = (rho, sigma)`, `u_x = (rho_x, sigma_x)` with
      !! `sigma_x = 2 grad rho . grad rho_x`, and `s_ab = 2 grad rho_a . grad rho_b`,
      !! differentiating the kernel's two coefficient groups gives
      !!
      !!     c_rho  = g_rrr rho_a rho_b + g_rrs (rho_a sigma_b + sigma_a rho_b)
      !!            + g_rss sigma_a sigma_b + f_rs s_ab
      !!     c_grad = 2 [g_rrs rho_a rho_b + g_rss (rho_a sigma_b + sigma_a rho_b)
      !!                 + g_sss sigma_a sigma_b + f_ss s_ab] grad rho
      !!            + 2 (f_rs rho_b + f_ss sigma_b) grad rho_a
      !!            + 2 (f_rs rho_a + f_ss sigma_a) grad rho_b
      !!
      !! symmetric under a <-> b, as a third variation must be. `f_rr` never
      !! appears: its own derivative is the `g_rrr`/`g_rrs` pair, and nothing
      !! multiplies it undifferentiated. Neither does `v_sigma`: its derivative
      !! is the `f_rs rho_a + f_ss sigma_a` factor on `grad rho_b`.
      !!
      !! **LDA and GGA; meta-GGA is refused**, for the same reason the third
      !! derivatives themselves refuse it: the tau channels of `g_xc` are not
      !! evaluated, and returning the GGA channels alone would be a converged,
      !! plausible, wrong operator.
      !!
      !! Both trial densities must be symmetric matrices -- the `grad rho_x`
      !! built here carries the same symmetry factor `eval_rho` documents.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)    !! The converged SCF density
      real(dp), intent(in) :: dtilde_a(:, :)   !! First response density
      real(dp), intent(in) :: dtilde_b(:, :)   !! Second response density
      real(dp), intent(inout) :: v_kernel(:, :)   !! Accumulated into
      type(error_t), intent(inout) :: error

#ifdef MQC_WITH_LIBXC
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :)
      real(dp), allocatable :: rho(:), rho_grad(:, :)
      real(dp), allocatable :: drho_a(:), drho_a_grad(:, :)
      real(dp), allocatable :: drho_b(:), drho_b_grad(:, :)
      real(dp), allocatable :: sigma(:), sig_a(:), sig_b(:), s_ab(:)
      real(dp), allocatable :: frs(:), fss(:)
      real(dp), allocatable :: grrr(:), grrs(:), grss(:), gsss(:)
      real(dp), allocatable :: frr_i(:), frs_i(:), fss_i(:)
      real(dp), allocatable :: grrr_i(:), grrs_i(:), grss_i(:), gsss_i(:)
      real(dp), allocatable :: c_rho(:), c_grad(:, :), no_tau(:)
      real(dp), allocatable :: v_sig(:, :), d_sig(:, :), da_sig(:, :), db_sig(:, :)
      real(dp), allocatable :: extents(:)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_list(:), ao_offset(:)
      integer :: n_sig, ia, ja
      integer :: g0, g1, nb, i, ig, id, npts
      logical :: gga

      if (.not. ctx%active) return
      if (ctx%polarized) then
         call error%set(ERROR_VALIDATION, "the second exchange-correlation kernel is "// &
                        "implemented for a restricted reference only")
         return
      end if
      if (ctx%any_mgga) then
         call error%set(ERROR_VALIDATION, "the second exchange-correlation kernel is "// &
                        "LDA and GGA only: a meta-GGA needs the tau channels of the "// &
                        "third functional derivative, which nothing here provides. "// &
                        "Refused rather than computed with those terms missing.")
         return
      end if

      gga = ctx%any_gga
      npts = ctx%grid%n_points

      call shell_extents(mol, ctx%screen_tol, extents)
      allocate (v_sig(mol%nao, mol%nao), d_sig(mol%nao, mol%nao), &
                da_sig(mol%nao, mol%nao), db_sig(mol%nao, mol%nao))
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(mol%nao))

      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         ! One screen for all three densities, as the kernel does for its two:
         ! every one of them is contracted against the same basis at the same
         ! points, so they share the kept set and one `ao_list` scatters back.
         call block_significant_aos(mol, ctx%grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig)
         if (n_sig == 0) cycle          ! empty space; no basis function reaches it

         do ja = 1, n_sig
            do ia = 1, n_sig
               d_sig(ia, ja) = density(ao_list(ia), ao_list(ja))
               da_sig(ia, ja) = dtilde_a(ao_list(ia), ao_list(ja))
               db_sig(ia, ja) = dtilde_b(ao_list(ia), ao_list(ja))
            end do
         end do

         if (gga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, &
                               grad=ao_grad, shell_mask=shell_mask, &
                               ao_offset=ao_offset, n_ao_out=n_sig)
         else
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, &
                               shell_mask=shell_mask, ao_offset=ao_offset, &
                               n_ao_out=n_sig)
         end if
         if (error%has_error()) return

         ! The reference density is what `g_xc` is evaluated at; the two
         ! response densities are what it multiplies. All three go through the
         ! same builder, so the symmetry factor in every gradient is the one
         ! convention `eval_rho` already owns.
         if (gga) then
            call eval_rho(ao, d_sig(1:n_sig, 1:n_sig), rho, ao_grad=ao_grad, &
                          rho_grad=rho_grad)
            call eval_rho(ao, da_sig(1:n_sig, 1:n_sig), drho_a, ao_grad=ao_grad, &
                          rho_grad=drho_a_grad)
            call eval_rho(ao, db_sig(1:n_sig, 1:n_sig), drho_b, ao_grad=ao_grad, &
                          rho_grad=drho_b_grad)
         else
            call eval_rho(ao, d_sig(1:n_sig, 1:n_sig), rho)
            call eval_rho(ao, da_sig(1:n_sig, 1:n_sig), drho_a)
            call eval_rho(ao, db_sig(1:n_sig, 1:n_sig), drho_b)
         end if

         if (allocated(frs)) deallocate (frs, fss, grrr, grrs, grss, gsss, &
                                         frr_i, frs_i, fss_i, &
                                         grrr_i, grrs_i, grss_i, gsss_i, &
                                         sigma, sig_a, sig_b, s_ab, c_rho)
         allocate (frs(nb), fss(nb), grrr(nb), grrs(nb), grss(nb), gsss(nb), &
                   frr_i(nb), frs_i(nb), fss_i(nb), &
                   grrr_i(nb), grrs_i(nb), grss_i(nb), gsss_i(nb), &
                   sigma(nb), sig_a(nb), sig_b(nb), s_ab(nb), c_rho(nb))
         frs = 0.0_dp
         fss = 0.0_dp
         grrr = 0.0_dp
         grrs = 0.0_dp
         grss = 0.0_dp
         gsss = 0.0_dp
         ! Zero rather than left alone on the LDA path: their coefficients are
         ! zero there, and `0 * uninitialised` is a NaN rather than nothing.
         sigma = 0.0_dp
         sig_a = 0.0_dp
         sig_b = 0.0_dp
         s_ab = 0.0_dp
         if (gga) then
            do ig = 1, nb
               sigma(ig) = rho_grad(ig, 1)**2 + rho_grad(ig, 2)**2 + rho_grad(ig, 3)**2
               sig_a(ig) = 2.0_dp*(rho_grad(ig, 1)*drho_a_grad(ig, 1) &
                                   + rho_grad(ig, 2)*drho_a_grad(ig, 2) &
                                   + rho_grad(ig, 3)*drho_a_grad(ig, 3))
               sig_b(ig) = 2.0_dp*(rho_grad(ig, 1)*drho_b_grad(ig, 1) &
                                   + rho_grad(ig, 2)*drho_b_grad(ig, 2) &
                                   + rho_grad(ig, 3)*drho_b_grad(ig, 3))
               s_ab(ig) = 2.0_dp*(drho_a_grad(ig, 1)*drho_b_grad(ig, 1) &
                                  + drho_a_grad(ig, 2)*drho_b_grad(ig, 2) &
                                  + drho_a_grad(ig, 3)*drho_b_grad(ig, 3))
            end do
         end if

         ! Per component, as everywhere else here: a composition may put an LDA
         ! correlation beside a GGA exchange, and `any_gga` only says that at
         ! least one of them needs sigma.
         do i = 1, ctx%n_func
            select case (ctx%family(i))
            case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
               ! `f_rs` and `f_ss` are second derivatives and come from the
               ! ordinary kernel evaluator. They belong here because the trial
               ! densities' own sigma cross-terms multiply them -- the two
               ! pieces of the third variation that are not third derivatives.
               call xc_f03_gga_fxc(ctx%func(i), int(nb, 8), rho, sigma, &
                                   frr_i, frs_i, fss_i)
               call xc_f03_gga_kxc(ctx%func(i), int(nb, 8), rho, sigma, &
                                   grrr_i, grrs_i, grss_i, gsss_i)
               frs = frs + ctx%weight(i)*frs_i
               fss = fss + ctx%weight(i)*fss_i
               grrr = grrr + ctx%weight(i)*grrr_i
               grrs = grrs + ctx%weight(i)*grrs_i
               grss = grss + ctx%weight(i)*grss_i
               gsss = gsss + ctx%weight(i)*gsss_i
            case default
               call xc_f03_lda_kxc(ctx%func(i), int(nb, 8), rho, grrr_i)
               grrr = grrr + ctx%weight(i)*grrr_i
            end select
         end do

         ! The same floor the kernel and the grid quantities apply, at the same
         ! density: every coefficient here sits at second order or above, and
         ! the third derivatives diverge faster at the tail than the second
         ! ones the floor was introduced for.
         do ig = 1, nb
            if (rho(ig) < KERNEL_RHO_FLOOR) then
               frs(ig) = 0.0_dp
               fss(ig) = 0.0_dp
               grrr(ig) = 0.0_dp
               grrs(ig) = 0.0_dp
               grss(ig) = 0.0_dp
               gsss(ig) = 0.0_dp
            end if
         end do

         do ig = 1, nb
            c_rho(ig) = grrr(ig)*drho_a(ig)*drho_b(ig) &
                        + grrs(ig)*(drho_a(ig)*sig_b(ig) + sig_a(ig)*drho_b(ig)) &
                        + grss(ig)*sig_a(ig)*sig_b(ig) &
                        + frs(ig)*s_ab(ig)
         end do
         if (gga) then
            if (allocated(c_grad)) deallocate (c_grad)
            allocate (c_grad(nb, 3))
            do id = 1, 3
               do ig = 1, nb
                  c_grad(ig, id) = 2.0_dp*(grrs(ig)*drho_a(ig)*drho_b(ig) &
                                           + grss(ig)*(drho_a(ig)*sig_b(ig) &
                                                       + sig_a(ig)*drho_b(ig)) &
                                           + gsss(ig)*sig_a(ig)*sig_b(ig) &
                                           + fss(ig)*s_ab(ig))*rho_grad(ig, id) &
                                   + 2.0_dp*(frs(ig)*drho_b(ig) &
                                             + fss(ig)*sig_b(ig))*drho_a_grad(ig, id) &
                                   + 2.0_dp*(frs(ig)*drho_a(ig) &
                                             + fss(ig)*sig_a(ig))*drho_b_grad(ig, id)
               end do
            end do
         end if

         ! The same assembly the potential and the kernel use, with this rung's
         ! coefficients. Writing a second one would be two copies of the
         ! arithmetic that is hardest to get right here.
         v_sig(1:n_sig, 1:n_sig) = 0.0_dp
         call accumulate_xc_matrix(ctx%grid%weights(g0:g1), ao, c_rho, &
                                   v_sig(1:n_sig, 1:n_sig), &
                                   ao_grad=ao_grad, grad_coeff=c_grad, vtau=no_tau, &
                                   any_gga=gga, any_mgga=.false.)
         do ja = 1, n_sig
            do ia = 1, n_sig
               v_kernel(ao_list(ia), ao_list(ja)) = &
                  v_kernel(ao_list(ia), ao_list(ja)) + v_sig(ia, ja)
            end do
         end do
      end do
#else
      call error%set(ERROR_VALIDATION, "no libxc in this build")
      if (size(density) < 0 .or. size(dtilde_a) < 0 .or. size(dtilde_b) < 0) return
      if (size(v_kernel) < 0) return
      if (mol%nao < 0) return
      if (ctx%n_func < 0) return
#endif
   end subroutine xc_kernel2_apply

   subroutine vv10_add_potential(ctx, mol, density, v_nl, e_nl, error)
      !! VV10's energy and Fock contribution, entirely on its own coarse grid
      !!
      !! **Both grids are the coarse one.** The non-local term is a double
      !! integral, so its cost goes as the product of the two grids' sizes;
      !! coarsening only the inner sum is not enough. PySCF makes the same
      !! choice. It costs a second AO pass, because a potential has to be
      !! contracted on the grid it was evaluated on, and that pass is linear in
      !! points on a grid an order of magnitude smaller than the exchange one.
      !!
      !! Two sweeps over that grid rather than one: the kernel's inner sum runs
      !! over every point, so nothing can be contracted until rho and sigma are
      !! known everywhere. Both are threaded over blocks, as is the double sum
      !! between them.
      !!
      !! No AO screening here, unlike the exchange loop: on a grid this size the
      !! saving is small against the indexing it would need.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(inout) :: v_nl(:, :)
      real(dp), intent(out) :: e_nl
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: rho(:), sigma(:), rho_grad(:, :)
      real(dp), allocatable :: exc(:), vrho(:), vsigma(:)
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :), grad_coeff(:, :)
      real(dp), allocatable :: rho_blk(:), rho_grad_blk(:, :), vtau_none(:)
      real(dp), allocatable :: v_local(:, :)
      integer :: npts, g0, g1, nb, ig, id
      integer :: n_kept
      type(error_t) :: local_error
      logical :: failed

      e_nl = 0.0_dp

      call ensure_nlc_grid(ctx, mol, error)
      if (error%has_error()) return

      npts = ctx%nlc_grid%n_points
      if (npts == 0) return

      allocate (rho(npts), sigma(npts), rho_grad(npts, 3))
      rho = 0.0_dp
      sigma = 0.0_dp
      rho_grad = 0.0_dp

      ! Sweep one: rho and sigma everywhere.
      !
      ! Threaded over blocks, as the semilocal loops are. Every block writes to
      ! its own slice of `rho`, `sigma` and `rho_grad` -- the index is the grid
      ! point, so there is nothing to reduce and no two threads touch the same
      ! element.
      !
      ! `firstprivate(local_error)` rather than `private`, for the reason the
      ! restricted path sets out.
      failed = .false.
      !$omp parallel do default(none) &
      !$omp    shared(ctx, mol, density, npts, rho, sigma, rho_grad, error, failed) &
      !$omp    private(g0, g1, nb, ig, id, ao, ao_grad, rho_blk, rho_grad_blk) &
      !$omp    firstprivate(local_error) &
      !$omp    schedule(dynamic)
      do g0 = 1, npts, ctx%point_block
         if (failed) cycle
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1
         if (allocated(rho_blk)) deallocate (rho_blk, rho_grad_blk)
         allocate (rho_blk(nb), rho_grad_blk(nb, 3))
         call eval_ao_block(mol, ctx%nlc_grid%coords(:, g0:g1), ao, local_error, grad=ao_grad)
         if (local_error%has_error()) then
            !$omp critical (vv10_rho_failure)
            if (.not. failed) then
               failed = .true.
               error = local_error
            end if
            !$omp end critical (vv10_rho_failure)
            cycle
         end if
         call eval_rho(ao, density, rho_blk, ao_grad=ao_grad, rho_grad=rho_grad_blk)
         do ig = 1, nb
            rho(g0 + ig - 1) = rho_blk(ig)
            do id = 1, 3
               rho_grad(g0 + ig - 1, id) = rho_grad_blk(ig, id)
            end do
            sigma(g0 + ig - 1) = rho_grad_blk(ig, 1)**2 + rho_grad_blk(ig, 2)**2 &
                                 + rho_grad_blk(ig, 3)**2
         end do
      end do
      !$omp end parallel do
      if (failed) return

      allocate (exc(npts), vrho(npts), vsigma(npts))
      call vv10_nlc(ctx%nlc_b, ctx%nlc_c, ctx%nlc_grid%coords, rho, sigma, &
                    ctx%nlc_grid%coords, rho, sigma, ctx%nlc_grid%weights, &
                    exc, vrho, vsigma, n_inner_kept=n_kept)

      ! TODO(mqc): this diagnostic is commented out, and it is the only consumer
      ! of `n_kept` and of the `logger` and `to_char` imports at the top of the
      ! module -- all three are now dead. Either restore the line or drop them.
      !call logger%verbose("  VV10: "//to_char(npts)//" grid points, "// &
      !                    to_char(n_kept)//" carry density, "// &
      !                    to_char(real(npts, dp)*real(n_kept, dp)/1.0e6_dp)// &
      !                    " million pairs")

      e_nl = sum(ctx%nlc_grid%weights*rho*exc)

      allocate (vtau_none(0))

      ! Sweep two: contract the potential on the grid it was evaluated on.
      !
      ! This one accumulates into a single matrix, so unlike sweep one it needs a
      ! reduction: each thread fills its own `v_local` and adds it in once at the
      ! end. `v_nl` itself is *not* zeroed here -- the restricted caller passes
      ! the matrix that already holds the semilocal potential.
      failed = .false.
      !$omp parallel default(none) &
      !$omp    shared(ctx, mol, npts, vrho, vsigma, rho_grad, v_nl, vtau_none, &
      !$omp           error, failed) &
      !$omp    private(g0, g1, nb, ig, id, ao, ao_grad, grad_coeff, v_local) &
      !$omp    firstprivate(local_error)
      allocate (v_local(size(v_nl, 1), size(v_nl, 2)))
      v_local = 0.0_dp
      !$omp do schedule(dynamic)
      do g0 = 1, npts, ctx%point_block
         if (failed) cycle
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1
         call eval_ao_block(mol, ctx%nlc_grid%coords(:, g0:g1), ao, local_error, grad=ao_grad)
         if (local_error%has_error()) then
            !$omp critical (vv10_pot_failure)
            if (.not. failed) then
               failed = .true.
               error = local_error
            end if
            !$omp end critical (vv10_pot_failure)
            cycle
         end if
         if (allocated(grad_coeff)) deallocate (grad_coeff)
         allocate (grad_coeff(nb, 3))
         ! dE/d(grad rho) = 2 vsigma grad rho, the same chain rule the
         ! semilocal path applies, so the two cannot disagree about it.
         do id = 1, 3
            do ig = 1, nb
               grad_coeff(ig, id) = 2.0_dp*vsigma(g0 + ig - 1)*rho_grad(g0 + ig - 1, id)
            end do
         end do
         call accumulate_xc_matrix(ctx%nlc_grid%weights(g0:g1), ao, vrho(g0:g1), &
                                   v_local, ao_grad=ao_grad, grad_coeff=grad_coeff, &
                                   vtau=vtau_none, any_gga=.true., any_mgga=.false.)
      end do
      !$omp end do
      !$omp critical (vv10_pot_reduce)
      v_nl = v_nl + v_local
      !$omp end critical (vv10_pot_reduce)
      !$omp end parallel
      if (failed) return
   end subroutine vv10_add_potential

   subroutine vv10_kernel_apply(ctx, mol, density, dtilde, v_kernel, error)
      !! The VV10 kernel applied to a batch of response densities
      !!
      !! `xc_kernel_apply`'s non-local counterpart, on the NLC grid: the second
      !! functional derivative of `E_nl` contracted against each trial density,
      !! which is the term a coupled-perturbed solve over a `-V` functional is
      !! missing without this. The VV10 potential matrix has two pieces,
      !! `f_rho chi_u chi_v` and `2 f_gamma grad rho . grad(chi_u chi_v)`, and
      !! both respond to a trial density:
      !!
      !!     d f_rho, d f_gamma  --  `vv10_hessian_kernel`, the operator the
      !!                             explicit Hessian and the Fock derivative
      !!                             already validated, one pair sweep for the
      !!                             whole batch
      !!     d grad rho          --  `2 f_gamma grad drho . grad(chi_u chi_v)`,
      !!                             the same term `xc_kernel_apply` carries as
      !!                             `2 v_sigma grad drho`
      !!
      !! PySCF's `get_vnlc_resp`, term for term.
      !!
      !! **Batched, unlike `xc_kernel_apply`.** The pair sweep is O(npts^2)
      !! whether it carries one trial or thirty, so applying this per
      !! perturbation would multiply the only expensive part of the routine by
      !! `3*natm` every iteration.
      !!
      !! **Accumulates into `v_kernel`**, as `xc_kernel_apply` does; the caller
      !! zeroes. Each trial must be symmetric, which is what `eval_rho`'s
      !! gradient assumes. Restricted only.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)      !! The converged SCF density
      real(dp), intent(in) :: dtilde(:, :, :)    !! (nao, nao, n_trial), each symmetric
      real(dp), intent(inout) :: v_kernel(:, :, :)  !! (nao, nao, n_trial), accumulated into
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: rho(:), sigma(:), rho_grad(:, :)
      real(dp), allocatable :: rho_blk(:), rho_grad_blk(:, :)
      real(dp), allocatable :: rho_t(:, :), gamma_t(:, :), grad_t(:, :, :)
      real(dp), allocatable :: exc(:), vrho(:), vsigma(:)
      real(dp), allocatable :: pu(:), pw(:), pa(:), pb(:), pc(:)
      real(dp), allocatable :: dodr(:), dodg(:), d2odr2(:), d2odg2(:), d2odrdg(:)
      real(dp), allocatable :: dkdr(:), d2kdr2(:)
      real(dp), allocatable :: f_rho_t(:, :), f_gamma_t(:, :)
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :)
      real(dp), allocatable :: grad_coeff(:, :), vtau_none(:), c_rho(:)
      integer :: npts, n_trial, g0, g1, nb, ig, id, it, g

      if (error%has_error()) return
      if (ctx%nlc_b == 0.0_dp .and. ctx%nlc_c == 0.0_dp) return
      if (ctx%polarized) then
         call error%set(ERROR_VALIDATION, "the VV10 kernel is implemented for "// &
                        "a restricted reference only")
         return
      end if

      call ensure_nlc_grid(ctx, mol, error)
      if (error%has_error()) return
      npts = ctx%nlc_grid%n_points
      if (npts == 0) return
      n_trial = size(dtilde, 3)
      if (n_trial == 0) return

      ! Sweep one: the reference density and every trial's, gradients included,
      ! over the whole NLC grid -- the pair sums read every point before any
      ! output exists. No AO screening, matching `vv10_add_potential`.
      allocate (rho(npts), sigma(npts), rho_grad(npts, 3))
      allocate (rho_t(n_trial, npts), gamma_t(n_trial, npts))
      allocate (grad_t(npts, 3, n_trial))
      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1
         call eval_ao_block(mol, ctx%nlc_grid%coords(:, g0:g1), ao, error, grad=ao_grad)
         if (error%has_error()) return
         call eval_rho(ao, density, rho_blk, ao_grad=ao_grad, rho_grad=rho_grad_blk)
         do ig = 1, nb
            rho(g0 + ig - 1) = rho_blk(ig)
            do id = 1, 3
               rho_grad(g0 + ig - 1, id) = rho_grad_blk(ig, id)
            end do
            sigma(g0 + ig - 1) = rho_grad_blk(ig, 1)**2 + rho_grad_blk(ig, 2)**2 &
                                 + rho_grad_blk(ig, 3)**2
         end do
         do it = 1, n_trial
            call eval_rho(ao, dtilde(:, :, it), rho_blk, ao_grad=ao_grad, &
                          rho_grad=rho_grad_blk)
            do ig = 1, nb
               rho_t(it, g0 + ig - 1) = rho_blk(ig)
               do id = 1, 3
                  grad_t(g0 + ig - 1, id, it) = rho_grad_blk(ig, id)
               end do
            end do
         end do
      end do
      ! gamma_t = 2 grad rho . grad drho -- the derivative of sigma along the
      ! trial, the same contraction `xc_kernel_apply` calls `dsigma`.
      do g = 1, npts
         do it = 1, n_trial
            gamma_t(it, g) = 2.0_dp*(rho_grad(g, 1)*grad_t(g, 1, it) &
                                     + rho_grad(g, 2)*grad_t(g, 2, it) &
                                     + rho_grad(g, 3)*grad_t(g, 3, it))
         end do
      end do

      ! One pair sweep for the potential -- whose `vsigma` *is* PySCF's
      ! `f_gamma` -- and every kernel intermediate the perturbed potential needs.
      allocate (exc(npts), vrho(npts), vsigma(npts))
      allocate (pu(npts), pw(npts), pa(npts), pb(npts), pc(npts))
      allocate (dodr(npts), dodg(npts), d2odr2(npts), d2odg2(npts), d2odrdg(npts))
      allocate (dkdr(npts), d2kdr2(npts))
      call vv10_nlc(ctx%nlc_b, ctx%nlc_c, ctx%nlc_grid%coords, rho, sigma, &
                    ctx%nlc_grid%coords, rho, sigma, ctx%nlc_grid%weights, &
                    exc, vrho, vsigma, &
                    hess_u=pu, hess_w=pw, hess_a=pa, hess_b=pb, hess_c=pc, &
                    domega_drho=dodr, domega_dgamma=dodg, &
                    d2omega_drho2=d2odr2, d2omega_dgamma2=d2odg2, &
                    d2omega_drho_dgamma=d2odrdg, &
                    dkappa_drho=dkdr, d2kappa_drho2=d2kdr2)

      ! The kernel: every trial in, the perturbed potentials `d f_rho` and
      ! `d f_gamma` out, in one pair sweep for the whole batch. The inner
      ! quadrature weight lives inside; the outer one is applied below.
      allocate (f_rho_t(n_trial, npts), f_gamma_t(n_trial, npts))
      call vv10_hessian_kernel(ctx%nlc_b, ctx%nlc_c, ctx%nlc_grid%coords, &
                               rho, sigma, ctx%nlc_grid%weights, &
                               pu, pw, pa, pb, pc, dodr, dodg, dkdr, &
                               d2odr2, d2odg2, d2odrdg, d2kdr2, &
                               rho_t, gamma_t, f_rho_t, f_gamma_t)

      ! Sweep two: contract each trial's perturbed potential on the grid it
      ! was evaluated on, with the potential path's own assembly --
      !
      !     dV_uv = int w [ d f_rho chi_u chi_v
      !                     + 2 (d f_gamma grad rho + f_gamma grad drho)
      !                       . grad(chi_u chi_v) ]
      !
      ! which is `accumulate_xc_matrix` with the kernel's coefficients in
      ! place of the potential's, exactly as the semilocal kernel reuses it.
      allocate (vtau_none(0))
      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1
         call eval_ao_block(mol, ctx%nlc_grid%coords(:, g0:g1), ao, error, grad=ao_grad)
         if (error%has_error()) return
         if (allocated(grad_coeff)) deallocate (grad_coeff, c_rho)
         allocate (grad_coeff(nb, 3), c_rho(nb))
         do it = 1, n_trial
            do ig = 1, nb
               g = g0 + ig - 1
               c_rho(ig) = f_rho_t(it, g)
               do id = 1, 3
                  grad_coeff(ig, id) = 2.0_dp*(f_gamma_t(it, g)*rho_grad(g, id) &
                                               + vsigma(g)*grad_t(g, id, it))
               end do
            end do
            call accumulate_xc_matrix(ctx%nlc_grid%weights(g0:g1), ao, c_rho, &
                                      v_kernel(:, :, it), ao_grad=ao_grad, &
                                      grad_coeff=grad_coeff, vtau=vtau_none, &
                                      any_gga=.true., any_mgga=.false.)
         end do
      end do
   end subroutine vv10_kernel_apply

   subroutine ensure_nlc_grid(ctx, mol, error)
      !! Build the non-local correlation grid, once, on first use
      !!
      !! Built on first use rather than beside the exchange grid, because whether
      !! it is needed is only known once the functional's components have been
      !! read. Shared by the potential and the gradient, so the two integrate and
      !! differentiate the same quadrature.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      type(error_t), intent(inout) :: error

      integer, allocatable :: numbers(:)

      if (ctx%nlc_grid%n_points > 0) return
      allocate (numbers(mol%natm))
      ! The *element*, not the charge it presents -- see the exchange grid above.
      numbers = nint(mol%charges) + mol%core_electrons
      call build_dft_grid(mol%coords, numbers, ctx%nlc_grid, error, &
                          level=ctx%nlc_grid_level)
      deallocate (numbers)
   end subroutine ensure_nlc_grid

   subroutine accumulate_xc_matrix(weights, ao, vrho, v, ao_grad, grad_coeff, vtau, &
                                   any_gga, any_mgga)
      !! One spin's exchange-correlation matrix, over one block of grid points
      !!
      !! Shared by the restricted and unrestricted paths. What differs between the
      !! two callers is not this arithmetic but what goes into `vrho`,
      !! `grad_coeff` and `vtau` -- for a spin-polarised functional those are that
      !! spin's derivatives, including the cross-spin gradient term.
      !!
      !! Accumulates into `v` rather than returning, so the caller's matrix is the
      !! sum over blocks.
      real(dp), intent(in) :: weights(:)      !! Grid weights for this block
      real(dp), intent(in) :: ao(:, :)        !! (n_points, n_ao)
      real(dp), intent(in) :: vrho(:)         !! dE/drho at each point
      real(dp), intent(inout) :: v(:, :)
      real(dp), allocatable, intent(in) :: ao_grad(:, :, :)   !! (n_points, n_ao, 3)
      real(dp), allocatable, intent(in) :: grad_coeff(:, :)   !! dE/d(grad rho), (n_points, 3)
      real(dp), allocatable, intent(in) :: vtau(:)            !! dE/dtau at each point
      logical, intent(in) :: any_gga, any_mgga

      real(dp), allocatable :: scaled(:, :)
      integer :: nb, nao, mu, ig, id

      nb = size(ao, 1)
      nao = size(ao, 2)
      allocate (scaled(nb, nao))

      ! V += (w v_rho chi)^T chi, as a gemm. Scaling the left factor rather than
      ! forming a diagonal matrix keeps this one multiply per element.
      do mu = 1, nao
         do ig = 1, nb
            scaled(ig, mu) = weights(ig)*vrho(ig)*ao(ig, mu)
         end do
      end do
      call pic_gemm(scaled, ao, v, transa="T", alpha=1.0_dp, beta=1.0_dp)

      ! The gradient term,
      !
      !     V_uv += sum_g w_g dE/dgrad rho . (grad chi_u chi_v + chi_u grad chi_v)
      !
      ! whose two halves are transposes of each other, so one gemm plus a
      ! symmetrisation does both.
      if (any_gga) then
         do mu = 1, nao
            do ig = 1, nb
               scaled(ig, mu) = 0.0_dp
               do id = 1, 3
                  scaled(ig, mu) = scaled(ig, mu) &
                                   + weights(ig)*grad_coeff(ig, id)*ao_grad(ig, mu, id)
               end do
            end do
         end do
         ! scaled^T ao gives one half; adding its transpose gives the other.
         call pic_gemm(scaled, ao, v, transa="T", alpha=1.0_dp, beta=1.0_dp)
         call pic_gemm(ao, scaled, v, transa="T", alpha=1.0_dp, beta=1.0_dp)
      end if

      ! The kinetic-energy-density term. d tau / d D_uv is half the sum over
      ! directions of grad chi_u grad chi_v, which is already symmetric in u and
      ! v -- so unlike the sigma term this one needs no transpose added.
      if (any_mgga) then
         do id = 1, 3
            do mu = 1, nao
               do ig = 1, nb
                  scaled(ig, mu) = 0.5_dp*weights(ig)*vtau(ig)*ao_grad(ig, mu, id)
               end do
            end do
            call pic_gemm(scaled, ao_grad(:, :, id), v, transa="T", &
                          alpha=1.0_dp, beta=1.0_dp)
         end do
      end if

      deallocate (scaled)
   end subroutine accumulate_xc_matrix

end module mqc_libcint_xc
