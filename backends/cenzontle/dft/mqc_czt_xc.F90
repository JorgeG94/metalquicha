!! Exchange-correlation, as one thing the SCF can be handed
module mqc_czt_xc
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
   use mqc_memory, only: memory_budget
   use mqc_czt_vv10, only: vv10_nlc, vv10_hessian_kernel
   use mqc_dft_grid, only: dft_grid_t, build_dft_grid, DEFAULT_GRID_LEVEL
   use mqc_xc_spec, only: xc_spec_t, xc_spec_from_name, MAX_XC_COMPONENTS
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_ao, only: eval_ao_block, eval_rho, AO_POINT_BLOCK, &
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

   integer, parameter, public :: NLC_GRID_LEVEL_DEFAULT = 1
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
   public :: xc_kernel_apply_many
   public :: xc_kernel_cache_t
   public :: xc_kernel_cache_fill
   public :: xc_kernel2_apply
   public :: xc_grid_kernel_quantities
   public :: KERNEL_RHO_FLOOR   !! Where the kernel's divergence is cut off
   public :: ensure_nlc_grid
   public :: vv10_add_potential
   public :: vv10_kernel_apply

   real(dp), parameter :: KERNEL_RHO_FLOOR = 1.0e-10_dp
   integer, parameter :: KERNEL_SET_CHUNK = 8
      !! Response densities contracted per gemm in `xc_kernel_apply_many`.
      !! Eight of 250 significant functions makes an inner dimension of two
      !! thousand; the four stacked work arrays are then a few tens of
      !! megabytes per thread.
      !! Grid points below this density contribute no *second* derivative.
      !!
      !! Not a tolerance -- a necessity. The LDA kernel is `d2e/drho2`, which for
      !! the exchange part goes as `rho^(-2/3)` and diverges at the tail of every
      !! atomic grid. Left in, those points swamp the real ones and the response
      !! operator stops being positive definite, so the conjugate-gradient solver
      !! reports a saddle point that is not there. `vrho` and `exc` stay finite
      !! where `v2rho2` does not, which is why nothing before the kernel needed a
      !! floor.

   real(dp), parameter :: KERNEL_CACHE_BUDGET_SHARE = 0.25_dp
      !! The share of the machine's available memory `xc_kernel_cache_fill` may
      !! plan on. Well below the response solver's own share, because the solve
      !! that asked for the cache is holding its trial vectors, its operators
      !! and the reference at the same time, and the cache only has to pay for
      !! itself: declining it costs a libxc pass per application, not the run.

   real(dp), parameter :: KERNEL_CACHE_BLIND_LIMIT = 2.0e9_dp
      !! What the cache may take, in bytes, on a machine that reports nothing
      !! about its memory. Two gigabytes is a meta-GGA grid of twenty-three
      !! million points, which is a few thousand atoms.

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
      integer :: nlc_grid_level = NLC_GRID_LEVEL_DEFAULT
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
      integer :: func_id(MAX_XC_COMPONENTS) = 0
         !! libxc's number for each component, kept because libxc fixes the spin
         !! channel when a functional is initialised and a handle cannot be
         !! asked to change it afterwards. The triplet kernel is a *polarised*
         !! evaluation of the same functional this context holds unpolarised, so
         !! it needs a second handle, and building one needs the number again.
      logical :: polarized_twin = .false.
         !! Whether `func_pol` holds initialised handles.
#ifdef MQC_WITH_LIBXC
      type(xc_f03_func_t) :: func(MAX_XC_COMPONENTS)
      type(xc_f03_func_t) :: func_pol(MAX_XC_COMPONENTS)
         !! The same components again, spin-polarised, for the triplet kernel.
         !! Uninitialised and unread until `polarized_twin`; created once per
         !! context by `ensure_polarized_twins` rather than once per grid block,
         !! which is what Psi4 does and pays for.
#endif
   contains
      procedure :: destroy => xc_context_destroy
   end type xc_context_t

   type :: xc_kernel_cache_t
      !! The kernel's coefficients on the whole grid, evaluated once
      !!
      !! `xc_kernel_apply_many` evaluates the reference density and libxc's
      !! second derivatives per grid block on **every** call, and every one of
      !! those is a property of the converged density alone: a response solve
      !! applying the same kernel a hundred times pays for the same numbers a
      !! hundred times. These are O(n_points) scalars, so holding them costs a
      !! handful of arrays over the grid and saves a libxc pass per call.
      !!
      !! **Not the basis functions**, which stay per block for the reason the
      !! module header gives: those are `n_points` by `n_ao` and a cache of them
      !! is the one thing here that does not fit.
      !!
      !! **What it costs, per grid point.** Only the channels the rung defines
      !! are held, so the figure is 8 bytes on an LDA (`frr` alone), 56 on a GGA
      !! (`rho_grad`, `frs`, `fss`, `vsig` beside it) and 88 on a meta-GGA (the
      !! four tau channels as well). A 200-atom level-5 grid is a few million
      !! points, so a GGA cache over it is a couple of hundred megabytes and a
      !! meta-GGA one half again as much. `xc_kernel_cache_fill` weighs that
      !! against `memory_budget` and declines rather than allocating over it;
      !! a declined cache comes back unfilled, and its consumer takes the
      !! uncached path.
      !!
      !! Filled by `xc_kernel_cache_fill` from the reference density, and handed
      !! back to `xc_kernel_apply` or `xc_kernel_apply_many` as an optional
      !! argument. Absent, those routines evaluate as they always did -- the
      !! same code, since the fill and the uncached path share
      !! `kernel_block_reference` rather than each having their own copy of it.
      !!
      !! A cache belongs to one context, one geometry and one reference
      !! density; nothing here can check the last two, so a consumer that moves
      !! the nuclei or re-converges the density has to fill again.
      logical :: filled = .false.
         !! Whether the arrays below hold anything. A cache that was never
         !! filled is refused rather than read as zeros, which would be a
         !! silently missing kernel term.
      integer :: n_points = 0
         !! Points the cache was filled over, checked against the context's grid.
      logical :: gga = .false.
      logical :: mgga = .false.
         !! What the context was when the fill ran. Checked on use, because a
         !! GGA cache read by a meta-GGA contraction is short three channels.
      integer :: point_block = 0
         !! The block width the fill walked the grid in.
      real(dp) :: screen_tol = 0.0_dp
         !! The AO screening tolerance the fill decided block membership with.
         !!
         !! Both are the context's, and both are checked on use: the fill and
         !! the contraction agree on where the blocks fall and on which
         !! functions reach them only because they read the same two numbers,
         !! and a cache filled against one context and applied against another
         !! that differs in either is misaligned point by point.
      real(dp), allocatable :: rho_grad(:, :)
         !! (n_points, 3) the reference density's gradient, which the
         !! contraction reads: it makes `dsigma` and multiplies the response's
         !! gradient coefficient. Unallocated on an LDA, where it is zero.
      real(dp), allocatable :: frr(:), frs(:), fss(:)
         !! (n_points) `v2rho2`, `v2rhosigma` and `v2sigma2`, summed over the
         !! functional's components with their weights and floored at
         !! `KERNEL_RHO_FLOOR`. Only `frr` is allocated on an LDA.
      real(dp), allocatable :: vsig(:)
         !! (n_points) `vsigma`, a *first* derivative, here because the
         !! response density's gradient multiplies it. Unallocated on an LDA.
      real(dp), allocatable :: frt(:), fst(:), ftt(:)
         !! (n_points) the kinetic-energy-density channels `v2rhotau`,
         !! `v2sigmatau` and `v2tau2`. Unallocated unless the context is a
         !! meta-GGA, where they are zero.
      real(dp), allocatable :: vtau(:)
         !! (n_points) `vtau`. Floored with the rest and carried for symmetry
         !! with `vsig`; the contraction does not read it, because tau is
         !! linear in the density and so has no analogue of the `vsigma` term.
         !! Unallocated off the meta-GGA rung, like the three above.
      logical :: triplet = .false.
         !! Whether the four arrays below were filled as well. A cache filled
         !! for singlets alone and then handed to a triplet contraction is
         !! refused rather than read as zeros, which would be a kernel silently
         !! left out of half the spectrum.
      real(dp), allocatable :: frr_t(:), frs_t(:), fss_t(:), vsig_t(:)
         !! (n_points) the same four coefficients for the **triplet** kernel:
         !! a spin-polarised evaluation of the same functional at
         !! `rho_a = rho_b = rho/2` and `sigma_aa = sigma_ab = sigma_bb =
         !! sigma/4`, combined as `kernel_block_reference` documents. Allocated
         !! where `triplet`, and per rung as the singlet channels are: `frr_t`
         !! always, the other three from the GGA up.
   contains
      procedure :: destroy => xc_kernel_cache_destroy
   end type xc_kernel_cache_t

   type :: xc_kernel_block_t
      !! One grid block's worth of the reference density and the kernel
      !!
      !! What the contraction needs on a block that does not depend on the
      !! response density, in the shapes it indexes: everything from one, not
      !! from the block's first grid point. `kernel_block_reference` evaluates
      !! one and `kernel_block_from_cache` copies one out of a filled cache,
      !! and they are interchangeable -- that they produce the same object is
      !! what makes the cache a cache rather than a second kernel.
      !!
      !! Every component is allocated whichever rung is in play, because the
      !! contraction multiplies all of them; the ones the rung does not define
      !! are zero.
      real(dp), allocatable :: rho(:)
         !! (n_block) the reference density. Evaluated on the uncached path
         !! because libxc is evaluated at it; unallocated on the cached one,
         !! where nothing downstream reads it.
      real(dp), allocatable :: rho_grad(:, :)
         !! (n_block, 3) its gradient. Zero off the GGA and meta-GGA rungs.
      real(dp), allocatable :: frr(:), frs(:), fss(:)
         !! (n_block) `v2rho2`, `v2rhosigma` and `v2sigma2`, weighted over the
         !! functional's components and floored at `KERNEL_RHO_FLOOR`.
      real(dp), allocatable :: vsig(:)
         !! (n_block) `vsigma`, a first derivative, which the response
         !! density's gradient multiplies.
      real(dp), allocatable :: frt(:), fst(:), ftt(:)
         !! (n_block) the kinetic-energy-density channels. Zero off the
         !! meta-GGA rung.
      real(dp), allocatable :: vtau(:)
         !! (n_block) `vtau`, floored with the rest and not read: tau is
         !! linear in the density, so it has no analogue of the `vsigma` term.
   end type xc_kernel_block_t

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
                                screen_tol, point_block, nlc_level, allow_half, &
                                n_radial, n_angular)
      !! Resolve a functional name and build the grid it will be integrated on
      type(czt_molecule_t), intent(in) :: mol
      character(len=*), intent(in) :: functional
      type(xc_context_t), intent(out) :: ctx
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: level
         !! Grid level, `DEFAULT_GRID_LEVEL` when absent. Negative means the
         !! deck gave explicit counts instead, which then have to arrive in
         !! `n_radial` and `n_angular`: a negative level on its own is refused,
         !! since the grid builder would clamp it to level 0 and integrate on
         !! the coarsest grid there is without a word.
      integer, intent(in), optional :: n_radial, n_angular
         !! Explicit per-atom radial shells and Lebedev order, in force when
         !! `level` is negative. `keywords.dft.radial_points` and
         !! `angular_points`.
      logical, intent(in), optional :: polarized
         !! Initialise the functionals spin-polarised, for an unrestricted
         !! calculation. Default restricted. Fixed here because libxc fixes it at
         !! initialisation, so it cannot be decided later by whoever evaluates.
      integer, intent(in), optional :: nlc_level
         !! Level for VV10's own quadrature. Absent or negative keeps
         !! `NLC_GRID_LEVEL_DEFAULT`; a deck reaches this through
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

      ctx%nlc_grid_level = NLC_GRID_LEVEL_DEFAULT
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
      ! The element itself, when the molecule recorded it: `charges` is zero
      ! on a ghost atom -- a counterpoise partner, or a nucleus quantised by
      ! NEO -- and a grid built for "element 0" there integrates a density
      ! that is anything but empty. The charge-plus-core sum is the fallback
      ! for a molecule assembled without the record.
      if (allocated(mol%atomic_numbers)) then
         numbers = mol%atomic_numbers
      else
         numbers = nint(mol%charges) + mol%core_electrons
      end if
      if (grid_level < 0) then
         if (.not. (present(n_radial) .and. present(n_angular))) then
            call error%set(ERROR_VALIDATION, "xc_context_create: a negative grid level "// &
                           "stands for explicit radial and angular counts, and none "// &
                           "were passed")
            return
         end if
         if (n_radial < 1 .or. n_angular < 1) then
            call error%set(ERROR_VALIDATION, "xc_context_create: the explicit grid counts "// &
                           "must both be positive")
            return
         end if
         call build_dft_grid(mol%coords, numbers, ctx%grid, error, &
                             n_radial=n_radial, n_angular=n_angular)
      else
         call build_dft_grid(mol%coords, numbers, ctx%grid, error, level=grid_level)
      end if
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
         ctx%func_id(i) = id
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
         ! evaluated by `mqc_czt_vv10` where the potential is assembled.
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
         if (this%polarized_twin) call xc_f03_func_end(this%func_pol(i))
      end do
#endif
      call this%grid%destroy()
      call this%nlc_grid%destroy()
      this%polarized_twin = .false.
      this%func_id = 0
      this%n_func = 0
      this%active = .false.
      this%exx_fraction = 0.0_dp
      this%pt2_fraction = 0.0_dp
      this%range_separated = .false.
      this%rs_omega = 0.0_dp
      this%rs_k_lr = 0.0_dp
   end subroutine xc_context_destroy

   subroutine xc_kernel_cache_destroy(this)
      !! Release the grid-sized arrays and mark the cache unfilled
      class(xc_kernel_cache_t), intent(inout) :: this

      if (allocated(this%rho_grad)) deallocate (this%rho_grad)
      if (allocated(this%frr)) deallocate (this%frr)
      if (allocated(this%frs)) deallocate (this%frs)
      if (allocated(this%fss)) deallocate (this%fss)
      if (allocated(this%vsig)) deallocate (this%vsig)
      if (allocated(this%frt)) deallocate (this%frt)
      if (allocated(this%fst)) deallocate (this%fst)
      if (allocated(this%ftt)) deallocate (this%ftt)
      if (allocated(this%vtau)) deallocate (this%vtau)
      if (allocated(this%frr_t)) deallocate (this%frr_t)
      if (allocated(this%frs_t)) deallocate (this%frs_t)
      if (allocated(this%fss_t)) deallocate (this%fss_t)
      if (allocated(this%vsig_t)) deallocate (this%vsig_t)
      this%filled = .false.
      this%n_points = 0
      this%gga = .false.
      this%mgga = .false.
      this%triplet = .false.
      this%point_block = 0
      this%screen_tol = 0.0_dp
   end subroutine xc_kernel_cache_destroy

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
      type(czt_molecule_t), intent(in) :: mol
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
      type(error_t) :: local_error
      logical :: failed

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
      !
      ! Threaded over blocks: every output is indexed by the point, so the
      ! blocks write disjoint ranges. `default(shared)` because the spin-beta
      ! arguments are optional dummies, and naming an absent one in a
      ! data-sharing clause is not portable.
      failed = .false.
      !$omp parallel default(shared) &
      !$omp    private(g0, g1, nb, i, ig, ao, rho_blk, exc_i, vrho_i, rho_a_blk, rho_b_blk) &
      !$omp    firstprivate(local_error)
      !$omp do schedule(dynamic)
      do g0 = 1, npts, ctx%point_block
         if (failed) cycle
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, local_error)
         if (local_error%has_error()) then
            !$omp critical (xc_lda_quantities_failure)
            if (.not. failed) then
               failed = .true.
               error = local_error
            end if
            !$omp end critical (xc_lda_quantities_failure)
            cycle
         end if

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
      !$omp end do
      !$omp end parallel
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
      type(czt_molecule_t), intent(in) :: mol
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
      type(error_t) :: local_error
      logical :: failed

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
      !
      ! Threaded over blocks, as `xc_grid_lda_quantities` and for the same
      ! reasons, `default(shared)` included.
      failed = .false.
      !$omp parallel default(shared) &
      !$omp    private(g0, g1, nb, i, ig, id, ao, ao_grad, rho_blk, sigma, exc_i, vrho_i, &
      !$omp            vsigma_i, vsigma, tau_blk, lapl, vlapl, vtau_i, rho_a_blk, rho_b_blk, &
      !$omp            grad_a, grad_b, rho_grad) &
      !$omp    firstprivate(local_error)
      !$omp do schedule(dynamic)
      do g0 = 1, npts, ctx%point_block
         if (failed) cycle
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, local_error, grad=ao_grad)
         if (local_error%has_error()) then
            !$omp critical (xc_gga_quantities_failure)
            if (.not. failed) then
               failed = .true.
               error = local_error
            end if
            !$omp end critical (xc_gga_quantities_failure)
            cycle
         end if

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
      !$omp end do
      !$omp end parallel
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
      type(czt_molecule_t), intent(in) :: mol
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
      type(error_t) :: local_error
      logical :: failed

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
      !
      ! Threaded over blocks: every output is indexed by the point, so the
      ! blocks write disjoint ranges and nothing is reduced. `default(shared)`
      ! rather than the usual `default(none)` because the tau and third-
      ! derivative outputs are optional dummies, and naming an absent one in
      ! a data-sharing clause is not portable.
      failed = .false.
      !$omp parallel default(shared) &
      !$omp    private(g0, g1, nb, i, ig, id, ao, ao_grad, rho_blk, grad_blk, sigma, &
      !$omp            exc_i, vrho_i, vsigma_i, frr_i, frs_i, fss_i, vtau_i, frt_i, fst_i, &
      !$omp            ftt_i, lapl_k, lscr, tau_k, lscr_rl, lscr_sl, lscr_ll, lscr_lt, &
      !$omp            grrr_i, grrs_i, grss_i, gsss_i) &
      !$omp    firstprivate(local_error)
      !$omp do schedule(dynamic)
      do g0 = 1, npts, ctx%point_block
         if (failed) cycle
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         ! The AO gradients are needed whenever anything here is a GGA, and
         ! harmless otherwise -- but they are the expensive half of the
         ! evaluation, so a pure LDA does not pay for them.
         if (ctx%any_gga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, local_error, grad=ao_grad)
         else
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, local_error)
         end if
         if (local_error%has_error()) then
            !$omp critical (xc_kernel_quantities_failure)
            if (.not. failed) then
               failed = .true.
               error = local_error
            end if
            !$omp end critical (xc_kernel_quantities_failure)
            cycle
         end if
         if (ctx%any_gga) then
            if (want_tau) then
               call eval_rho(ao, density, rho_blk, ao_grad=ao_grad, &
                             rho_grad=grad_blk, tau=tau_k)
            else
               call eval_rho(ao, density, rho_blk, ao_grad=ao_grad, rho_grad=grad_blk)
            end if
         else
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
      !$omp end do
      !$omp end parallel
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
      type(czt_molecule_t), intent(in) :: mol
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
      type(czt_molecule_t), intent(in) :: mol
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

   subroutine ensure_polarized_twins(ctx, error)
      !! Initialise a spin-polarised handle beside each of the context's own
      !!
      !! libxc decides at initialisation whether a functional is polarised, so
      !! the triplet kernel -- which is a polarised evaluation at
      !! `rho_a = rho_b = rho/2` -- cannot reuse the handles a restricted
      !! calculation built. It gets its own, once per context: the handles cost
      !! nothing to keep and `xc_context_destroy` ends them with the rest.
      !!
      !! **Not callable from inside a parallel region.** It writes to `ctx`, and
      !! every caller invokes it before the block loop rather than in it.
      type(xc_context_t), intent(inout) :: ctx
      type(error_t), intent(inout) :: error

      integer :: i

      if (ctx%polarized_twin) return
      if (ctx%polarized) then
         ! The twins exist to be polarised beside unpolarised ones. A context
         ! that is already polarised has no restricted kernel to be the triplet
         ! partner of, and the unrestricted response is Layer 6.
         call error%set(ERROR_VALIDATION, "a spin-polarised exchange-correlation "// &
                        "context was asked for the triplet kernel of a closed shell")
         return
      end if
#ifdef MQC_WITH_LIBXC
      ! Every number is checked before any handle is created. `polarized_twin`
      ! is the only record that `func_pol` holds anything, and
      ! `xc_context_destroy` ends the twins all-or-nothing on the strength of
      ! it -- so a run that created some and then refused would leak exactly
      ! the ones it had made.
      do i = 1, ctx%n_func
         if (ctx%func_id(i) <= 0) then
            call error%set(ERROR_VALIDATION, "this exchange-correlation context does "// &
                           "not carry the libxc numbers of its components, so no "// &
                           "polarised twin of them can be built for the triplet kernel")
            return
         end if
      end do
      do i = 1, ctx%n_func
         call xc_f03_func_init(ctx%func_pol(i), ctx%func_id(i), XC_POLARIZED)
      end do
      ctx%polarized_twin = .true.
#else
      call error%set(ERROR_VALIDATION, "no libxc in this build")
      i = 0
#endif
   end subroutine ensure_polarized_twins

#ifdef MQC_WITH_LIBXC
   subroutine kernel_block_reference(ctx, ao, ao_grad, d_sig, blk, error, triplet)
      !! The reference density and the kernel's coefficients on one grid block
      !!
      !! Everything the kernel contraction needs that does not depend on the
      !! response density: the reference `rho` and its gradient, libxc's second
      !! derivatives summed over the functional's components with their
      !! weights, and the two first derivatives a response *gradient*
      !! multiplies.
      !!
      !! **One routine rather than two copies.** `kernel_apply_batch` calls it
      !! per block and `xc_kernel_cache_fill` calls it over the whole grid, and
      !! a cache that rounded differently from the path it replaces would be a
      !! second implementation of the kernel rather than a cache of the first.
      !! `KERNEL_RHO_FLOOR` is applied here for the same reason.
      !!
      !! The rung is read off `ctx` rather than passed, so the two callers
      !! cannot disagree about it either.
      !!
      !! ## The triplet coefficients
      !!
      !! With `triplet`, the four GGA coefficients come out of a *polarised*
      !! evaluation of the same functional at the closed-shell point instead.
      !! They go in the same slots of the same block and are consumed by the
      !! same contraction, so the caller changes nothing but this flag.
      type(xc_context_t), intent(inout) :: ctx
      real(dp), intent(in) :: ao(:, :)   !! (n_block, n_sig), the kept functions
      real(dp), allocatable, intent(in) :: ao_grad(:, :, :)
         !! (n_block, n_sig, 3), and unallocated where the context is an LDA
      real(dp), intent(in) :: d_sig(:, :)
         !! (n_sig, n_sig), the reference density over the kept functions
      type(xc_kernel_block_t), intent(out) :: blk
         !! Every channel allocated over the block, and zero on a rung that
         !! does not define it
      type(error_t), intent(inout) :: error
         !! Set only by the triplet precondition below. Every caller is inside
         !! a parallel region, so this is a thread's own `error_t` that the
         !! region's usual critical section promotes.
      logical, intent(in), optional :: triplet
         !! Return the triplet kernel's coefficients rather than the singlet's.
         !! Needs `ensure_polarized_twins` to have run on `ctx`, which is the
         !! caller's job because this is called from inside a parallel region.
         !! Refused rather than trusted, in `triplet_component`. Off by default.

      logical :: gga, mgga
      real(dp), allocatable :: rho(:), rho_grad(:, :)
      real(dp), allocatable :: frr(:), frs(:), fss(:), vsig(:)
      real(dp), allocatable :: frt(:), fst(:), ftt(:), vtau(:)
      real(dp), allocatable :: sigma(:), tau(:), lapl(:), lapl_scratch(:)
         !! libxc's meta-GGA entry points take the Laplacian and return its
         !! derivatives whether or not the functional uses them. Laplacian
         !! dependent functionals are refused at construction, so these are
         !! zeros in and discarded out.
      real(dp), allocatable :: exc_i(:), vrho_i(:), vsigma_i(:)
      real(dp), allocatable :: frr_i(:), frs_i(:), fss_i(:)
      real(dp), allocatable :: frt_i(:), fst_i(:), ftt_i(:), vtau_i(:)
      integer :: nb, i, ig
      logical :: want_triplet

      nb = size(ao, 1)
      mgga = ctx%any_mgga
      ! A meta-GGA needs the density gradients a GGA needs, and tau on top.
      gga = ctx%any_gga .or. mgga
      want_triplet = .false.
      if (present(triplet)) want_triplet = triplet

      ! The reference density, once for the block: what `f_xc` is evaluated
      ! at. Its tau only where a meta-GGA asks.
      if (mgga) then
         call eval_rho(ao, d_sig, rho, ao_grad=ao_grad, rho_grad=rho_grad, tau=tau)
      else if (gga) then
         call eval_rho(ao, d_sig, rho, ao_grad=ao_grad, rho_grad=rho_grad)
      else
         call eval_rho(ao, d_sig, rho)
      end if
      ! Allocated and zero rather than absent on the rungs that do not define
      ! them: the coefficients multiplying them are zero there, and zero times
      ! uninitialised is a NaN rather than nothing.
      if (.not. allocated(rho_grad)) then
         allocate (rho_grad(nb, 3))
         rho_grad = 0.0_dp
      end if
      if (.not. allocated(tau)) then
         allocate (tau(nb))
         tau = 0.0_dp
      end if

      allocate (frr(nb), frs(nb), fss(nb), vsig(nb), frt(nb), fst(nb), ftt(nb), vtau(nb))
      frr = 0.0_dp
      frs = 0.0_dp
      fss = 0.0_dp
      vsig = 0.0_dp
      frt = 0.0_dp
      fst = 0.0_dp
      ftt = 0.0_dp
      vtau = 0.0_dp
      allocate (sigma(nb), lapl(nb), lapl_scratch(nb), exc_i(nb), vrho_i(nb), &
                vsigma_i(nb), frr_i(nb), frs_i(nb), fss_i(nb), frt_i(nb), fst_i(nb), &
                ftt_i(nb), vtau_i(nb))
      lapl = 0.0_dp
      sigma = 0.0_dp
      if (gga) then
         do ig = 1, nb
            sigma(ig) = rho_grad(ig, 1)**2 + rho_grad(ig, 2)**2 + rho_grad(ig, 3)**2
         end do
      end if

      ! Per component, as everywhere else here: a composition may put an LDA
      ! correlation beside a GGA exchange, and `any_gga` only says that at
      ! least one of them needs sigma.
      do i = 1, ctx%n_func
         if (want_triplet) then
            ! A different functional derivative in the same four slots. The
            ! callers refuse a meta-GGA triplet, so no tau channel reaches here.
            call triplet_component(ctx, i, nb, rho, sigma, frr, frs, fss, vsig, error)
            if (error%has_error()) return
            cycle
         end if
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
            ! density's gradient multiplies one and its tau the other.
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

      call move_alloc(rho, blk%rho)
      call move_alloc(rho_grad, blk%rho_grad)
      call move_alloc(frr, blk%frr)
      call move_alloc(frs, blk%frs)
      call move_alloc(fss, blk%fss)
      call move_alloc(vsig, blk%vsig)
      call move_alloc(frt, blk%frt)
      call move_alloc(fst, blk%fst)
      call move_alloc(ftt, blk%ftt)
      call move_alloc(vtau, blk%vtau)
   end subroutine kernel_block_reference

   subroutine triplet_component(ctx, i, nb, rho, sigma, frr, frs, fss, vsig, error)
      !! One component's contribution to the triplet kernel, accumulated
      !!
      !! The polarised twin of component `i` evaluated at the closed-shell
      !! point -- `rho_a = rho_b = rho/2`, `sigma_aa = sigma_ab = sigma_bb =
      !! sigma/4` -- and its derivatives combined into the four coefficients
      !! the restricted contraction already consumes, times `ctx%weight(i)`.
      !!
      !! `ctx%polarized_twin` has to be set, and the caller runs
      !! `ensure_polarized_twins` outside its parallel region to do it. That
      !! is checked here rather than assumed: libxc leaves `func_pol` as
      !! uninitialised storage until then, and handing one of those to
      !! `xc_f03_gga_fxc` is a segmentation fault, not an error return.
      type(xc_context_t), intent(in) :: ctx
      integer, intent(in) :: i     !! Which component
      integer, intent(in) :: nb    !! Points in the block
      real(dp), intent(in) :: rho(:)      !! (nb) the total reference density
      real(dp), intent(in) :: sigma(:)    !! (nb) `|grad rho|^2`, zero for an LDA
      real(dp), intent(inout) :: frr(:), frs(:), fss(:), vsig(:)
         !! (nb) accumulated into, in the singlet coefficients' own slots
      type(error_t), intent(inout) :: error

      ! ## Where the four combinations come from
      !
      ! The restricted contraction reads a response density `dn` and forms
      !
      !     c_rho  = F_rr dn + F_rs dsigma,          dsigma = 2 grad rho . grad dn
      !     c_grad = 2 (F_rs dn + F_ss dsigma) grad rho + 2 V_sig grad dn
      !
      ! and `c_rho`, `c_grad` are the response of the *alpha* potential. A
      ! singlet perturbation is `drho_a = drho_b = dn/2` and a triplet one is
      ! `drho_a = -drho_b = dn/2`; nothing else about the contraction differs,
      ! so the triplet kernel is the same four slots filled from the polarised
      ! derivatives of the second perturbation.
      !
      ! At the closed-shell point that perturbation gives, with
      ! `s = grad rho . grad dn`,
      !
      !     dsigma_aa = s/2,   dsigma_ab = 0,   dsigma_bb = -s/2
      !
      ! -- the cross term vanishes because `grad rho_a = grad rho_b` while the
      ! two perturbations are opposite, which is the whole reason the triplet
      ! kernel is a *difference* of derivatives. Then
      !
      !     dv_a^rho  = (f_ra,ra - f_ra,rb) dn/2 + (f_ra,saa - f_ra,sbb) s/2
      !     dv_a^grad = de_saa grad rho + (e_saa - e_sab/2) grad dn
      !     de_saa    = (f_saa,ra - f_saa,rb) dn/2 + (f_saa,saa - f_saa,sbb) s/2
      !
      ! using `de_sab = 0`: `sigma_ab` is symmetric under a spin swap, so
      ! `f_sab,ra = f_sab,rb` and `f_sab,saa = f_sab,sbb`, and both differences
      ! are zero. Matching term by term against the two lines at the top, and
      ! remembering that `dsigma = 2s` is what the contraction forms:
      !
      !     F_rr  = (f_ra,ra  - f_ra,rb ) / 2
      !     F_rs  = (f_ra,saa - f_ra,sbb) / 4
      !     F_ss  = (f_saa,saa - f_saa,sbb) / 8
      !     V_sig = e_saa / 2 - e_sab / 4
      !
      ! which is the list Psi4 builds in `libfock/v.cc`, in `RV::compute_Vx`,
      ! derived here in this code's own `(rho, sigma)` representation rather
      ! than ported from PySCF's transformed one. For an LDA only the first survives, and there
      ! it equals the singlet coefficient whenever the functional carries no
      ! opposite-spin correlation -- pure exchange has `f_ra,rb = 0`.
      !
      ! TODO(mqc): no tau channels. A meta-GGA triplet kernel needs
      ! `f_ra,ta - f_ra,tb`, `f_saa,ta - f_saa,tb` and `f_ta,ta - f_ta,tb`
      ! alongside these; every caller refuses a meta-GGA triplet rather than
      ! returning these four and calling it the kernel.

      real(dp), allocatable :: rho_pol(:), sigma_pol(:)
      real(dp), allocatable :: frr_p(:), frs_p(:), fss_p(:)
      real(dp), allocatable :: exc_p(:), vrho_p(:), vsigma_p(:)
      real(dp) :: w
      integer :: ig

      if (.not. ctx%polarized_twin) then
         call error%set(ERROR_VALIDATION, "the triplet exchange-correlation kernel "// &
                        "was asked for against a context whose polarised libxc "// &
                        "handles were never built; ensure_polarized_twins has to "// &
                        "run on it first, outside any parallel region")
         return
      end if

      w = ctx%weight(i)
      allocate (rho_pol(2*nb), sigma_pol(3*nb), frr_p(3*nb))
      do ig = 1, nb
         rho_pol(2*ig - 1) = 0.5_dp*rho(ig)
         rho_pol(2*ig) = 0.5_dp*rho(ig)
         sigma_pol(3*ig - 2) = 0.25_dp*sigma(ig)
         sigma_pol(3*ig - 1) = 0.25_dp*sigma(ig)
         sigma_pol(3*ig) = 0.25_dp*sigma(ig)
      end do

      select case (ctx%family(i))
      case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
         allocate (frs_p(6*nb), fss_p(6*nb), exc_p(nb), vrho_p(2*nb), vsigma_p(3*nb))
         call xc_f03_gga_fxc(ctx%func_pol(i), int(nb, 8), rho_pol, sigma_pol, &
                             frr_p, frs_p, fss_p)
         call xc_f03_gga_exc_vxc(ctx%func_pol(i), int(nb, 8), rho_pol, sigma_pol, &
                                 exc_p, vrho_p, vsigma_p)
         ! libxc's polarised layout, spin fastest: `v2rho2` is (aa, ab, bb),
         ! `v2rhosigma` is (a|aa, a|ab, a|bb, b|aa, b|ab, b|bb) and `v2sigma2`
         ! is (aa|aa, aa|ab, aa|bb, ab|ab, ab|bb, bb|bb).
         do ig = 1, nb
            frr(ig) = frr(ig) + w*0.5_dp*(frr_p(3*ig - 2) - frr_p(3*ig - 1))
            frs(ig) = frs(ig) + w*0.25_dp*(frs_p(6*ig - 5) - frs_p(6*ig - 3))
            fss(ig) = fss(ig) + w*0.125_dp*(fss_p(6*ig - 5) - fss_p(6*ig - 3))
            vsig(ig) = vsig(ig) + w*(0.5_dp*vsigma_p(3*ig - 2) &
                                     - 0.25_dp*vsigma_p(3*ig - 1))
         end do
      case default
         ! An LDA component, and a meta-GGA one only because the caller has
         ! already refused that combination: its sigma and tau channels would
         ! be missing rather than wrong.
         call xc_f03_lda_fxc(ctx%func_pol(i), int(nb, 8), rho_pol, frr_p)
         do ig = 1, nb
            frr(ig) = frr(ig) + w*0.5_dp*(frr_p(3*ig - 2) - frr_p(3*ig - 1))
         end do
      end select
   end subroutine triplet_component
#endif

   subroutine xc_kernel_apply(ctx, mol, density, dtilde, v_kernel, error, cache)
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
      !! **LDA, GGA and meta-GGA.** For a GGA the response density has a
      !! gradient too, and both of the potential's pieces respond:
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
      !! A meta-GGA adds `v2rhotau`, `v2sigmatau` and `v2tau2` and a tau
      !! component of the response density. This routine forwards to
      !! `xc_kernel_apply_many`, which carries all three, so the rung is
      !! whatever the context is rather than a choice made here. The header
      !! said otherwise until the tau channels landed in `_many` and it was
      !! left behind.
      type(xc_context_t), intent(inout) :: ctx
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)   !! The converged SCF density
      real(dp), intent(in) :: dtilde(:, :)    !! The response density
      real(dp), intent(inout) :: v_kernel(:, :)  !! Accumulated into
      type(error_t), intent(inout) :: error
      type(xc_kernel_cache_t), intent(in), optional :: cache
         !! The reference's kernel coefficients, from `xc_kernel_cache_fill`.
         !! Absent, they are evaluated here as they always were.

      real(dp), allocatable :: many_in(:, :, :), many_out(:, :, :)

      if (.not. ctx%active) return
      ! A batch of one. The copies are two matrices, against a grid pass.
      allocate (many_in(size(dtilde, 1), size(dtilde, 2), 1), &
                many_out(size(v_kernel, 1), size(v_kernel, 2), 1))
      many_in(:, :, 1) = dtilde
      many_out = 0.0_dp
      if (present(cache)) then
         call xc_kernel_apply_many(ctx, mol, density, many_in, many_out, error, cache=cache)
      else
         call xc_kernel_apply_many(ctx, mol, density, many_in, many_out, error)
      end if
      if (error%has_error()) return
      v_kernel = v_kernel + many_out(:, :, 1)
   end subroutine xc_kernel_apply

   subroutine xc_kernel_apply_many(ctx, mol, density, dtildes, v_kernels, error, cache, &
                                   triplet)
      !! The exchange-correlation kernel applied to a batch of response densities
      !!
      !! `xc_kernel_apply` for `n_set` densities in one pass over the grid. The
      !! basis functions, the reference density and its kernel are evaluated
      !! once per block and every set is contracted against them; what a set
      !! costs on top of the first is its own density on the block and one
      !! matrix assembly. Accumulates into `v_kernels`, as `xc_kernel_apply`
      !! does into its one.
      !!
      !! **A set is skipped where its density is negligible.** A response to
      !! one atom's displacement is small far from that atom, and the bound
      !! `max|D| (sum_u |chi_u|)^2` on its density over the block, against the
      !! set's largest element scaled by `ctx%screen_tol`, says so before
      !! anything is contracted. Relative to the set itself, so a trial vector
      !! of any size is screened the same way.
      !!
      !! **Sets are contracted in stacks.** The densities of up to
      !! `KERNEL_SET_CHUNK` sets are gathered side by side into one matrix,
      !! contracted against the block's basis in one gemm with a long inner
      !! dimension, and their output matrices assembled by one gemm the same
      !! way. Set by set the same arithmetic streams a half-megabyte density
      !! through a gemm too small to hide it, and on a hundred cores at once
      !! that ran at a third of the rate this reaches. A meta-GGA keeps the
      !! set-by-set path, since its tau needs three more gemms per set.
      !!
      !! **Given a filled `cache`, the reference is not evaluated here at
      !! all.** The basis functions still are, per block, and so is every
      !! response density; what the cache removes is `eval_rho` on the
      !! reference and the libxc calls, which is what a Davidson or a
      !! coupled-perturbed solve repeats identically on every application. The
      !! answer is the same to the bit, because the fill and the uncached path
      !! are the same routine.
      !!
      !! The body is `kernel_apply_batch`; this exists to turn one optional
      !! argument into two ordinary ones. An absent optional dummy cannot
      !! portably be named in an OpenMP data-sharing clause, and the
      !! contraction's parallel region is `default(none)`, which is worth
      !! keeping over a shorter call chain.
      type(xc_context_t), intent(inout) :: ctx
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)       !! The converged SCF density
      real(dp), intent(in) :: dtildes(:, :, :)    !! (n_ao, n_ao, n_set), the response densities
      real(dp), intent(inout) :: v_kernels(:, :, :)  !! (n_ao, n_ao, n_set), accumulated into
      type(error_t), intent(inout) :: error
      type(xc_kernel_cache_t), intent(in), optional :: cache
         !! The reference's kernel coefficients over the whole grid, from
         !! `xc_kernel_cache_fill`. Present, the per-block evaluation of the
         !! reference density and of libxc's derivatives is skipped and these
         !! are read instead; absent, nothing changes.
      logical, intent(in), optional :: triplet
         !! Contract with the **triplet** kernel `(f_aa - f_ab)/2` rather than
         !! the singlet one. Off by default. With a cache, the cache has to have
         !! been filled for triplets too; without one, the polarised evaluation
         !! is made per block.

      type(xc_kernel_cache_t) :: unfilled
         !! Stands in for an absent `cache`, so the worker's argument is never
         !! optional. Holds no arrays, and `have_cache` keeps it unread.
      logical :: want_triplet

      want_triplet = .false.
      if (present(triplet)) want_triplet = triplet
      if (present(cache)) then
         call kernel_apply_batch(ctx, mol, density, dtildes, v_kernels, error, cache, &
                                 .true., want_triplet)
      else
         call kernel_apply_batch(ctx, mol, density, dtildes, v_kernels, error, unfilled, &
                                 .false., want_triplet)
      end if
   end subroutine xc_kernel_apply_many

   subroutine kernel_apply_batch(ctx, mol, density, dtildes, v_kernels, error, &
                                 cache, have_cache, want_triplet)
      !! The batched kernel contraction, with the cache made unconditional
      !!
      !! The body of `xc_kernel_apply_many`, whose docstring carries the
      !! contract: the grid pass, the per-set screen, the stacking of up to
      !! `KERNEL_SET_CHUNK` densities per gemm and what a filled cache removes.
      !! What is only true down here is the pair of arguments -- `cache` is an
      !! ordinary dummy and `have_cache` says whether to read it, because an
      !! absent optional cannot be named in the `default(none)` clause of the
      !! parallel region below.
      !!
      !! **Written under a lock per set.** The block's result for a set is the
      !! significant functions squared, and it goes straight into that set's
      !! output; there is no copy per thread, so the memory is the batch and
      !! nothing times the thread count.
!$    use omp_lib, only: omp_lock_kind, omp_init_lock, omp_set_lock, omp_unset_lock, omp_destroy_lock
      type(xc_context_t), intent(inout) :: ctx
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)       !! The converged SCF density
      real(dp), intent(in) :: dtildes(:, :, :)    !! (n_ao, n_ao, n_set), the response densities
      real(dp), intent(inout) :: v_kernels(:, :, :)  !! (n_ao, n_ao, n_set), accumulated into
      type(error_t), intent(inout) :: error
      type(xc_kernel_cache_t), intent(in) :: cache
         !! Read only where `have_cache`, and unfilled otherwise.
      logical, intent(in) :: have_cache
      logical, intent(in) :: want_triplet
         !! The triplet kernel rather than the singlet one, in the same slots.

#ifdef MQC_WITH_LIBXC
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :)
      real(dp), allocatable :: drho(:), drho_grad(:, :)
      real(dp), allocatable :: dsigma(:)
      real(dp), allocatable :: c_rho(:), c_grad(:, :), c_tau(:), no_tau(:)
      real(dp), allocatable :: dtau(:)
      real(dp), allocatable :: v_sig(:, :), d_sig(:, :), dt_sig(:, :)
      type(xc_kernel_block_t) :: blk
         !! The reference and the kernel over the block one thread is on,
         !! cached or evaluated
      real(dp), allocatable :: extents(:), dmax(:)
      real(dp), allocatable :: dt_all(:, :), x_all(:, :), s_all(:, :), m_all(:, :), wg(:)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_list(:), ao_offset(:), keep(:)
      integer :: n_sig, ia, ja, iset, n_set, nk, ik, nc, c, col0, mu
      type(error_t) :: local_error
      logical :: failed
      integer :: g0, g1, nb, i, ig, id, npts
      real(dp) :: amax, agmax, dmax_blk, s
      logical :: gga, mgga
!$    integer(omp_lock_kind), allocatable :: locks(:)

      if (.not. ctx%active) return
      if (ctx%polarized) then
         call error%set(ERROR_VALIDATION, "the exchange-correlation kernel is "// &
                        "implemented for a restricted reference only")
         return
      end if
      n_set = size(dtildes, 3)
      if (size(v_kernels, 3) /= n_set) then
         call error%set(ERROR_VALIDATION, "xc_kernel_apply_many: as many outputs as "// &
                        "response densities")
         return
      end if
      if (n_set == 0) return

      gga = ctx%any_gga
      mgga = ctx%any_mgga
      ! A meta-GGA needs the AO gradients for the same reason a GGA does, and
      ! then some: tau is built from them too.
      gga = gga .or. mgga
      npts = ctx%grid%n_points
      failed = .false.

      ! A cache belonging to another context is the failure worth catching
      ! here: its arrays are the right shape often enough, and a GGA's
      ! coefficients read by a meta-GGA contraction are three channels short
      ! of the answer with nothing to show for it. The block width and the
      ! screening tolerance are in the test because the fill and this loop
      ! agree on where the blocks fall, and on which functions reach them,
      ! only by reading the same two numbers off the same context.
      if (have_cache) then
         if (.not. cache%filled) then
            call error%set(ERROR_VALIDATION, "the exchange-correlation kernel was "// &
                           "given a cache that was never filled: call "// &
                           "xc_kernel_cache_fill on the reference density first")
            return
         end if
         if (cache%n_points /= npts .or. (cache%gga .neqv. gga) &
             .or. (cache%mgga .neqv. mgga) &
             .or. cache%point_block /= ctx%point_block &
             .or. cache%screen_tol /= ctx%screen_tol) then
            call error%set(ERROR_VALIDATION, "the exchange-correlation kernel cache "// &
                           "was filled for a different grid, a different "// &
                           "functional rung or a different blocking than the "// &
                           "context it was handed with")
            return
         end if
         if (want_triplet .and. .not. cache%triplet) then
            call error%set(ERROR_VALIDATION, "the triplet exchange-correlation kernel "// &
                           "was asked for against a cache filled for singlets only: "// &
                           "call xc_kernel_cache_fill with triplet = .true.")
            return
         end if
      end if
      if (want_triplet) then
         if (mgga) then
            call error%set(ERROR_VALIDATION, "the triplet exchange-correlation "// &
                           "kernel of a meta-GGA is not implemented: its tau "// &
                           "channels have no spin-difference form here")
            return
         end if
         ! Outside the parallel region, and once: the polarised handles are a
         ! property of the functional, not of a grid block.
         if (.not. have_cache) then
            call ensure_polarized_twins(ctx, error)
            if (error%has_error()) return
         end if
      end if

      call shell_extents(mol, ctx%screen_tol, extents)

      ! The largest element of each set, which its screen is relative to.
      allocate (dmax(n_set))
      do iset = 1, n_set
         dmax(iset) = maxval(abs(dtildes(:, :, iset)))
      end do

!$    allocate (locks(n_set))
!$    do iset = 1, n_set
!$       call omp_init_lock(locks(iset))
!$    end do

      !$omp parallel default(none) &
      !$omp    shared(ctx, mol, density, dtildes, v_kernels, error, failed, &
      !$omp           gga, mgga, npts, extents, n_set, dmax, locks, &
      !$omp           cache, have_cache, want_triplet) &
      !$omp    private(g0, g1, nb, i, ig, id, ao, ao_grad, blk, drho, &
      !$omp            drho_grad, dsigma, c_rho, c_grad, &
      !$omp            c_tau, no_tau, dtau) &
      !$omp    private(v_sig, d_sig, dt_sig, shell_mask, ao_list, ao_offset, &
      !$omp            n_sig, ia, ja, iset, amax, agmax, dmax_blk, s, &
      !$omp            dt_all, x_all, s_all, m_all, wg, keep, nk, ik, nc, c, col0, mu) &
      !$omp    firstprivate(local_error)
      allocate (v_sig(mol%nao, mol%nao), d_sig(mol%nao, mol%nao), &
                dt_sig(mol%nao, mol%nao))
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(mol%nao), keep(n_set))

      !$omp do schedule(dynamic)
      do g0 = 1, npts, ctx%point_block
         ! A thread that has failed stops working, but the loop still has to
         ! run out so every thread reaches the barrier.
         if (failed) cycle

         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         ! One screen for every density. The reference and the responses are
         ! contracted against the same basis at the same points, so they share
         ! the kept set and every result scatters back through one `ao_list`.
         call block_significant_aos(mol, ctx%grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig)
         if (n_sig == 0) cycle          ! empty space; no basis function reaches it

         ! The reference's own density over the kept functions, which is what
         ! the kernel is evaluated at. Not gathered when the cache already
         ! holds the answer -- it is read for nothing else.
         if (.not. have_cache) then
            do ja = 1, n_sig
               do ia = 1, n_sig
                  d_sig(ia, ja) = density(ao_list(ia), ao_list(ja))
               end do
            end do
         end if

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
            !$omp critical (xc_kernel_many_failure)
            if (.not. failed) then
               failed = .true.
               error = local_error
            end if
            !$omp end critical (xc_kernel_many_failure)
            cycle
         end if

         ! The reference density's gradient and the kernel's coefficients on
         ! this block: read from the cache where there is one, evaluated where
         ! there is not. The evaluation is one routine shared with the fill,
         ! so the two paths cannot drift -- including where the floor falls.
         if (have_cache) then
            call kernel_block_from_cache(cache, g0, g1, blk, want_triplet)
         else
            call kernel_block_reference(ctx, ao, ao_grad, d_sig(1:n_sig, 1:n_sig), &
                                        blk, local_error, triplet=want_triplet)
            if (local_error%has_error()) then
               !$omp critical (xc_kernel_many_failure)
               if (.not. failed) then
                  failed = .true.
                  error = local_error
               end if
               !$omp end critical (xc_kernel_many_failure)
               cycle
            end if
         end if

         ! The response side's own scratch, which no cache can hold: it is the
         ! trial density's, not the reference's.
         if (allocated(c_rho)) deallocate (c_rho, dsigma, c_tau, dtau)
         allocate (c_rho(nb), dsigma(nb), c_tau(nb), dtau(nb))
         c_tau = 0.0_dp
         ! Zero on the LDA and GGA paths, where the coefficient multiplying it
         ! is zero too and zero times uninitialised is a NaN rather than nothing.
         dtau = 0.0_dp

         ! How large the basis is on this block, for the per-set screen: the
         ! largest over the points of `sum_u |chi_u|`, and of the same sum over
         ! the gradient components.
         amax = 0.0_dp
         do ig = 1, nb
            s = 0.0_dp
            do i = 1, n_sig
               s = s + abs(ao(ig, i))
            end do
            amax = max(amax, s)
         end do
         agmax = 0.0_dp
         if (gga) then
            do ig = 1, nb
               s = 0.0_dp
               do id = 1, 3
                  do i = 1, n_sig
                     s = s + abs(ao_grad(ig, i, id))
                  end do
               end do
               agmax = max(agmax, s)
            end do
         end if

         if (mgga) then
            do iset = 1, n_set
               if (dmax(iset) == 0.0_dp) cycle
               dmax_blk = 0.0_dp
               do ja = 1, n_sig
                  do ia = 1, n_sig
                     dt_sig(ia, ja) = dtildes(ao_list(ia), ao_list(ja), iset)
                     dmax_blk = max(dmax_blk, abs(dt_sig(ia, ja)))
                  end do
               end do
               ! `|drho| <= max|D| (sum |chi|)^2` and `|grad drho| <= 2 max|D|
               ! (sum |chi|)(sum |grad chi|)`, both against the set's own scale.
               if (dmax_blk*amax*max(amax, 2.0_dp*agmax) < ctx%screen_tol*dmax(iset)) then
                  cycle
               end if

               if (mgga) then
                  ! **The tau convention never has to be decided here.** Whatever
                  ! `eval_rho` means by tau is what the energy path fed libxc and
                  ! what `accumulate_xc_matrix` differentiates, and the response
                  ! density goes through the same routine.
                  call eval_rho(ao, dt_sig(1:n_sig, 1:n_sig), drho, ao_grad=ao_grad, &
                                rho_grad=drho_grad, tau=dtau)
               else if (gga) then
                  call eval_rho(ao, dt_sig(1:n_sig, 1:n_sig), drho, ao_grad=ao_grad, &
                                rho_grad=drho_grad)
               else
                  call eval_rho(ao, dt_sig(1:n_sig, 1:n_sig), drho)
               end if

               dsigma = 0.0_dp
               if (gga) then
                  do ig = 1, nb
                     dsigma(ig) = 2.0_dp*(blk%rho_grad(ig, 1)*drho_grad(ig, 1) &
                                          + blk%rho_grad(ig, 2)*drho_grad(ig, 2) &
                                          + blk%rho_grad(ig, 3)*drho_grad(ig, 3))
                  end do
               end if

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
                  c_rho(ig) = blk%frr(ig)*drho(ig) + blk%frs(ig)*dsigma(ig) &
                              + blk%frt(ig)*dtau(ig)
               end do
               if (gga) then
                  if (allocated(c_grad)) deallocate (c_grad)
                  allocate (c_grad(nb, 3))
                  do id = 1, 3
                     do ig = 1, nb
                        c_grad(ig, id) = 2.0_dp*(blk%frs(ig)*drho(ig) &
                                                 + blk%fss(ig)*dsigma(ig) &
                                                 + blk%fst(ig)*dtau(ig))*blk%rho_grad(ig, id) &
                                         + 2.0_dp*blk%vsig(ig)*drho_grad(ig, id)
                     end do
                  end do
               end if
               if (mgga) then
                  do ig = 1, nb
                     c_tau(ig) = blk%frt(ig)*drho(ig) + blk%fst(ig)*dsigma(ig) &
                                 + blk%ftt(ig)*dtau(ig)
                  end do
               end if

               ! The same assembly the potential uses, with the kernel's
               ! coefficients in place of the potential's.
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

!$             call omp_set_lock(locks(iset))
               do ja = 1, n_sig
                  do ia = 1, n_sig
                     v_kernels(ao_list(ia), ao_list(ja), iset) = &
                        v_kernels(ao_list(ia), ao_list(ja), iset) + v_sig(ia, ja)
                  end do
               end do
!$             call omp_unset_lock(locks(iset))
            end do
         else
            ! Which sets are worth contracting on this block, their densities
            ! gathered side by side: `|drho| <= max|D| (sum |chi|)^2` and
            ! `|grad drho| <= 2 max|D| (sum |chi|)(sum |grad chi|)`, both
            ! against the set's own scale.
            if (allocated(dt_all)) deallocate (dt_all, x_all, s_all, m_all, wg)
            allocate (dt_all(n_sig, n_sig*min(n_set, KERNEL_SET_CHUNK)), &
                      x_all(nb, n_sig*min(n_set, KERNEL_SET_CHUNK)), &
                      s_all(nb, n_sig*min(n_set, KERNEL_SET_CHUNK)), &
                      m_all(n_sig*min(n_set, KERNEL_SET_CHUNK), n_sig), wg(nb))
            wg = ctx%grid%weights(g0:g1)
            if (allocated(drho)) deallocate (drho)
            allocate (drho(nb))
            if (gga) then
               if (allocated(drho_grad)) deallocate (drho_grad)
               if (allocated(c_grad)) deallocate (c_grad)
               allocate (drho_grad(nb, 3), c_grad(nb, 3))
            end if

            nk = 0
            do iset = 1, n_set
               if (dmax(iset) == 0.0_dp) cycle
               dmax_blk = 0.0_dp
               do ja = 1, n_sig
                  do ia = 1, n_sig
                     dmax_blk = max(dmax_blk, abs(dtildes(ao_list(ia), ao_list(ja), iset)))
                  end do
               end do
               if (dmax_blk*amax*max(amax, 2.0_dp*agmax) < ctx%screen_tol*dmax(iset)) then
                  cycle
               end if
               nk = nk + 1
               keep(nk) = iset
            end do

            do ik = 1, nk, KERNEL_SET_CHUNK
               nc = min(KERNEL_SET_CHUNK, nk - ik + 1)
               do c = 1, nc
                  col0 = (c - 1)*n_sig
                  iset = keep(ik + c - 1)
                  do ja = 1, n_sig
                     do ia = 1, n_sig
                        dt_all(ia, col0 + ja) = dtildes(ao_list(ia), ao_list(ja), iset)
                     end do
                  end do
               end do

               ! `X = chi D` for every set of the stack in one gemm; then per
               ! set the row-wise dots `eval_rho` takes, the kernel
               ! coefficients, and the scaled left factor of the assembly.
               call pic_gemm(ao(1:nb, 1:n_sig), dt_all(1:n_sig, 1:nc*n_sig), &
                             x_all(1:nb, 1:nc*n_sig), beta=0.0_dp)
               do c = 1, nc
                  col0 = (c - 1)*n_sig
                  drho = 0.0_dp
                  do mu = 1, n_sig
                     do ig = 1, nb
                        drho(ig) = drho(ig) + x_all(ig, col0 + mu)*ao(ig, mu)
                     end do
                  end do
                  dsigma = 0.0_dp
                  if (gga) then
                     drho_grad = 0.0_dp
                     do id = 1, 3
                        do mu = 1, n_sig
                           do ig = 1, nb
                              drho_grad(ig, id) = drho_grad(ig, id) &
                                                  + 2.0_dp*x_all(ig, col0 + mu)*ao_grad(ig, mu, id)
                           end do
                        end do
                     end do
                     do ig = 1, nb
                        dsigma(ig) = 2.0_dp*(blk%rho_grad(ig, 1)*drho_grad(ig, 1) &
                                             + blk%rho_grad(ig, 2)*drho_grad(ig, 2) &
                                             + blk%rho_grad(ig, 3)*drho_grad(ig, 3))
                     end do
                  end if

                  !     dv_rho   = f_rr drho + f_rs dsigma
                  !     dv_sigma = f_rs drho + f_ss dsigma
                  !
                  ! and the assembly is `M + M^T` with
                  ! `M = (w dv_rho chi / 2 + w dv_grad . grad chi)^T chi`, as
                  ! `accumulate_xc_matrix` builds it.
                  do ig = 1, nb
                     c_rho(ig) = blk%frr(ig)*drho(ig) + blk%frs(ig)*dsigma(ig)
                  end do
                  if (gga) then
                     do id = 1, 3
                        do ig = 1, nb
                           c_grad(ig, id) = 2.0_dp*(blk%frs(ig)*drho(ig) &
                                                    + blk%fss(ig)*dsigma(ig)) &
                                            *blk%rho_grad(ig, id) &
                                            + 2.0_dp*blk%vsig(ig)*drho_grad(ig, id)
                        end do
                     end do
                     do mu = 1, n_sig
                        do ig = 1, nb
                           s_all(ig, col0 + mu) = wg(ig)*(0.5_dp*c_rho(ig)*ao(ig, mu) &
                                                          + c_grad(ig, 1)*ao_grad(ig, mu, 1) &
                                                          + c_grad(ig, 2)*ao_grad(ig, mu, 2) &
                                                          + c_grad(ig, 3)*ao_grad(ig, mu, 3))
                        end do
                     end do
                  else
                     do mu = 1, n_sig
                        do ig = 1, nb
                           s_all(ig, col0 + mu) = 0.5_dp*wg(ig)*c_rho(ig)*ao(ig, mu)
                        end do
                     end do
                  end if
               end do

               ! `M` for the whole stack in one gemm, then each set's
               ! `M + M^T` into its output under its lock.
               call pic_gemm(s_all(1:nb, 1:nc*n_sig), ao(1:nb, 1:n_sig), &
                             m_all(1:nc*n_sig, 1:n_sig), transa="T", beta=0.0_dp)
               do c = 1, nc
                  col0 = (c - 1)*n_sig
                  iset = keep(ik + c - 1)
!$                call omp_set_lock(locks(iset))
                  do ja = 1, n_sig
                     do ia = 1, n_sig
                        v_kernels(ao_list(ia), ao_list(ja), iset) = &
                           v_kernels(ao_list(ia), ao_list(ja), iset) &
                           + m_all(col0 + ia, ja) + m_all(col0 + ja, ia)
                     end do
                  end do
!$                call omp_unset_lock(locks(iset))
               end do
            end do
         end if
      end do
      !$omp end do
      deallocate (v_sig, d_sig, dt_sig, shell_mask, ao_offset, ao_list, keep)
      if (allocated(dt_all)) deallocate (dt_all, x_all, s_all, m_all, wg)
      !$omp end parallel

!$    do iset = 1, n_set
!$       call omp_destroy_lock(locks(iset))
!$    end do
!$    deallocate (locks)
#else
      call error%set(ERROR_VALIDATION, "no libxc in this build")
      if (size(density) < 0 .or. size(dtildes) < 0 .or. size(v_kernels) < 0) return
      if (mol%nao < 0) return
      if (ctx%n_func < 0) return
      if (have_cache .and. cache%n_points < 0) return
      if (want_triplet) return
#endif
   end subroutine kernel_apply_batch

   subroutine kernel_block_from_cache(cache, g0, g1, blk, triplet)
      !! One block's worth of cached coefficients, in the shapes the block wants
      !!
      !! What `kernel_block_reference` would have evaluated for this block,
      !! read out of a filled cache instead. A copy rather than a slice,
      !! because the contraction indexes everything from one and a pointer
      !! into the middle of a cache array would have to be rebased anyway. It
      !! is ten `n_block` copies against a libxc pass.
      !!
      !! **Every channel comes back allocated, whatever the cache holds.** The
      !! cache only stores the ones its rung defines; the rest are zero
      !! everywhere and are written as zeros here, because the contraction
      !! multiplies them unconditionally and zero times unallocated is not
      !! zero. `blk%rho` is the exception and stays unallocated: the cache
      !! does not hold it and nothing downstream reads it.
      type(xc_kernel_cache_t), intent(in) :: cache
      integer, intent(in) :: g0, g1   !! First and last grid point of the block
      type(xc_kernel_block_t), intent(out) :: blk
      logical, intent(in) :: triplet
         !! Read the triplet coefficients into the same four slots. The caller
         !! has checked that the cache carries them.

      integer :: nb, id

      nb = g1 - g0 + 1
      allocate (blk%rho_grad(nb, 3))
      if (allocated(cache%rho_grad)) then
         do id = 1, 3
            blk%rho_grad(:, id) = cache%rho_grad(g0:g1, id)
         end do
      else
         blk%rho_grad = 0.0_dp
      end if
      allocate (blk%frr(nb), blk%frs(nb), blk%fss(nb), blk%vsig(nb), &
                blk%frt(nb), blk%fst(nb), blk%ftt(nb), blk%vtau(nb))
      ! The triplet coefficients go in the singlet slots: the contraction that
      ! reads them is the same one, and which manifold the four numbers belong
      ! to is the caller's to know.
      if (triplet) then
         call cached_channel(cache%frr_t, g0, g1, blk%frr)
         call cached_channel(cache%frs_t, g0, g1, blk%frs)
         call cached_channel(cache%fss_t, g0, g1, blk%fss)
         call cached_channel(cache%vsig_t, g0, g1, blk%vsig)
      else
         blk%frr = cache%frr(g0:g1)
         call cached_channel(cache%frs, g0, g1, blk%frs)
         call cached_channel(cache%fss, g0, g1, blk%fss)
         call cached_channel(cache%vsig, g0, g1, blk%vsig)
      end if
      call cached_channel(cache%frt, g0, g1, blk%frt)
      call cached_channel(cache%fst, g0, g1, blk%fst)
      call cached_channel(cache%ftt, g0, g1, blk%ftt)
      call cached_channel(cache%vtau, g0, g1, blk%vtau)
   end subroutine kernel_block_from_cache

   subroutine cached_channel(channel, g0, g1, block_values)
      !! One cached channel over one block, or zeros where the rung has none
      real(dp), allocatable, intent(in) :: channel(:)
         !! (n_points) over the whole grid, unallocated off its own rung
      integer, intent(in) :: g0, g1   !! First and last grid point of the block
      real(dp), intent(out) :: block_values(:)   !! (n_block)

      if (allocated(channel)) then
         block_values = channel(g0:g1)
      else
         block_values = 0.0_dp
      end if
   end subroutine cached_channel

   subroutine xc_kernel_cache_fill(ctx, mol, density, cache, error, triplet)
      !! Evaluate the kernel's coefficients over the whole grid, once
      !!
      !! The grid pass `xc_kernel_apply_many` makes on every call, made here
      !! instead and kept. What comes back is a `xc_kernel_cache_t` that any
      !! number of later applications can be handed; each of them then
      !! evaluates the basis functions and the response densities per block
      !! and nothing else.
      !!
      !! **The same blocks, the same screen, the same routine.** The loop below
      !! is the contraction's own prologue -- `ctx%point_block` points at a
      !! time, `block_significant_aos` deciding which functions reach them,
      !! `kernel_block_reference` doing the evaluation -- so a cached
      !! application reproduces an uncached one bit for bit rather than to
      !! within a tolerance. `ctx%point_block` and `ctx%screen_tol` are
      !! recorded in the cache so that the agreement is checked on use rather
      !! than assumed.
      !!
      !! **Only the rung's own channels are held.** An LDA defines `frr` and
      !! nothing else, a GGA adds the density gradient and three more, and the
      !! four tau channels arrive with the meta-GGA -- 8, 56 and 88 bytes per
      !! grid point. The rest are zero everywhere, and
      !! `kernel_block_from_cache` hands the contraction zeros for them.
      !!
      !! **It can decline.** The arrays are O(n_points) and a large grid makes
      !! them hundreds of megabytes, so the total is weighed against
      !! `memory_budget` first. Over budget, nothing is allocated, the cache
      !! comes back with `filled` false and no error is raised: every consumer
      !! already has an uncached path, and taking it costs a libxc pass per
      !! application rather than the run. The decision is logged.
      !!
      !! Restricted only, as the kernel itself is.
      type(xc_context_t), intent(inout) :: ctx
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)   !! The converged SCF density
      type(xc_kernel_cache_t), intent(out) :: cache
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: triplet
         !! Fill the triplet coefficients beside the singlet ones, so a solve
         !! over both manifolds walks the grid once rather than twice. Off by
         !! default. Refused for a meta-GGA, whose triplet kernel is missing
         !! its tau channels.

      logical :: want_triplet
#ifdef MQC_WITH_LIBXC
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :)
      real(dp), allocatable :: d_sig(:, :), extents(:)
      type(xc_kernel_block_t) :: blk
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_list(:), ao_offset(:)
      integer :: g0, g1, ia, ja, id, n_sig, npts, n_channels
      logical :: gga, mgga, failed
      real(dp) :: wanted, budget
      type(error_t) :: local_error
#endif

      want_triplet = .false.
      if (present(triplet)) want_triplet = triplet

      if (.not. ctx%active) return
      if (.not. xc_available()) then
         call error%set(ERROR_VALIDATION, "no libxc in this build")
         return
      end if
      if (ctx%polarized) then
         call error%set(ERROR_VALIDATION, "the exchange-correlation kernel is "// &
                        "implemented for a restricted reference only")
         return
      end if
      if (want_triplet) then
         if (ctx%any_mgga) then
            call error%set(ERROR_VALIDATION, "the triplet exchange-correlation "// &
                           "kernel of a meta-GGA is not implemented: its tau "// &
                           "channels have no spin-difference form here, and "// &
                           "filling the three that do exist would be a kernel "// &
                           "missing a term rather than a meta-GGA one")
            return
         end if
         call ensure_polarized_twins(ctx, error)
         if (error%has_error()) return
      end if

#ifdef MQC_WITH_LIBXC
      gga = ctx%any_gga .or. ctx%any_mgga
      mgga = ctx%any_mgga
      npts = ctx%grid%n_points

      ! What this rung actually holds: `frr` always, the gradient and three
      ! more channels on a GGA, the four tau channels on a meta-GGA. One
      ! `n_points` vector each, and `rho_grad` counts three.
      n_channels = 1
      if (gga) n_channels = n_channels + 6
      if (mgga) n_channels = n_channels + 4
      ! The triplet twins of `frr` and, from the GGA up, of `frs`, `fss` and
      ! `vsig`. Weighed with the rest, so a cache that only just fits for one
      ! manifold is not silently taken for two.
      if (want_triplet) then
         n_channels = n_channels + 1
         if (gga) n_channels = n_channels + 3
      end if
      wanted = real(npts, dp)*real(n_channels, dp)*8.0_dp
      budget = memory_budget(KERNEL_CACHE_BLIND_LIMIT, KERNEL_CACHE_BUDGET_SHARE)
      if (wanted > budget) then
         call logger%info("  xc kernel cache: declined, it wants "// &
                          to_char(nint(wanted/1.048576e6_dp))//" MB against a budget of "// &
                          to_char(nint(budget/1.048576e6_dp))//" MB; the kernel is "// &
                          "re-evaluated on each application instead")
         return
      end if

      allocate (cache%frr(npts))
      ! Zero where no basis function reaches, which is where the contraction
      ! skips the block outright and never reads these.
      cache%frr = 0.0_dp
      if (gga) then
         allocate (cache%rho_grad(npts, 3), cache%frs(npts), cache%fss(npts), &
                   cache%vsig(npts))
         cache%rho_grad = 0.0_dp
         cache%frs = 0.0_dp
         cache%fss = 0.0_dp
         cache%vsig = 0.0_dp
      end if
      if (mgga) then
         allocate (cache%frt(npts), cache%fst(npts), cache%ftt(npts), cache%vtau(npts))
         cache%frt = 0.0_dp
         cache%fst = 0.0_dp
         cache%ftt = 0.0_dp
         cache%vtau = 0.0_dp
      end if
      cache%n_points = npts
      cache%gga = gga
      cache%mgga = mgga
      cache%point_block = ctx%point_block
      cache%screen_tol = ctx%screen_tol
      cache%triplet = want_triplet
      if (want_triplet) then
         ! Per rung, as the singlet channels above are: `frr_t` on every rung,
         ! the other three from the GGA up. A meta-GGA triplet is refused
         ! before it gets here, so there is no tau twin to hold.
         allocate (cache%frr_t(npts))
         cache%frr_t = 0.0_dp
         if (gga) then
            allocate (cache%frs_t(npts), cache%fss_t(npts), cache%vsig_t(npts))
            cache%frs_t = 0.0_dp
            cache%fss_t = 0.0_dp
            cache%vsig_t = 0.0_dp
         end if
      end if

      call shell_extents(mol, ctx%screen_tol, extents)
      failed = .false.

      ! Threaded over blocks and nothing reduced: every point belongs to one
      ! block, so the threads write disjoint slices of each array.
      !$omp parallel default(none) &
      !$omp    shared(ctx, mol, density, cache, error, failed, gga, mgga, npts, &
      !$omp           extents, want_triplet) &
      !$omp    private(g0, g1, ia, ja, id, n_sig, ao, ao_grad, blk, d_sig, &
      !$omp            shell_mask, ao_list, ao_offset) &
      !$omp    firstprivate(local_error)
      allocate (d_sig(mol%nao, mol%nao))
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(mol%nao))

      !$omp do schedule(dynamic)
      do g0 = 1, npts, ctx%point_block
         if (failed) cycle
         g1 = min(g0 + ctx%point_block - 1, npts)

         call block_significant_aos(mol, ctx%grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig)
         if (n_sig == 0) cycle          ! empty space; no basis function reaches it

         do ja = 1, n_sig
            do ia = 1, n_sig
               d_sig(ia, ja) = density(ao_list(ia), ao_list(ja))
            end do
         end do

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
            !$omp critical (xc_kernel_cache_failure)
            if (.not. failed) then
               failed = .true.
               error = local_error
            end if
            !$omp end critical (xc_kernel_cache_failure)
            cycle
         end if

         call kernel_block_reference(ctx, ao, ao_grad, d_sig(1:n_sig, 1:n_sig), &
                                     blk, local_error)
         if (local_error%has_error()) then
            !$omp critical (xc_kernel_cache_failure)
            if (.not. failed) then
               failed = .true.
               error = local_error
            end if
            !$omp end critical (xc_kernel_cache_failure)
            cycle
         end if

         cache%frr(g0:g1) = blk%frr
         if (gga) then
            do id = 1, 3
               cache%rho_grad(g0:g1, id) = blk%rho_grad(:, id)
            end do
            cache%frs(g0:g1) = blk%frs
            cache%fss(g0:g1) = blk%fss
            cache%vsig(g0:g1) = blk%vsig
         end if
         if (mgga) then
            cache%frt(g0:g1) = blk%frt
            cache%fst(g0:g1) = blk%fst
            cache%ftt(g0:g1) = blk%ftt
            cache%vtau(g0:g1) = blk%vtau
         end if

         ! The second, polarised evaluation on the same block. It re-enters
         ! `kernel_block_reference` rather than extending it, so the two
         ! coefficient sets are produced by one routine under one floor; what
         ! is repeated is `eval_rho` and the libxc calls, not the basis
         ! functions, which are the expensive part and are passed back in.
         if (want_triplet) then
            call kernel_block_reference(ctx, ao, ao_grad, d_sig(1:n_sig, 1:n_sig), &
                                        blk, local_error, triplet=.true.)
            if (local_error%has_error()) then
               !$omp critical (xc_kernel_cache_failure)
               if (.not. failed) then
                  failed = .true.
                  error = local_error
               end if
               !$omp end critical (xc_kernel_cache_failure)
               cycle
            end if
            cache%frr_t(g0:g1) = blk%frr
            if (gga) then
               cache%frs_t(g0:g1) = blk%frs
               cache%fss_t(g0:g1) = blk%fss
               cache%vsig_t(g0:g1) = blk%vsig
            end if
         end if
      end do
      !$omp end do
      deallocate (d_sig, shell_mask, ao_offset, ao_list)
      !$omp end parallel

      if (failed) then
         call cache%destroy()
         return
      end if
      cache%filled = .true.
#else
      if (size(density) < 0 .or. mol%nao < 0 .or. want_triplet) return
#endif
   end subroutine xc_kernel_cache_fill

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
      type(czt_molecule_t), intent(in) :: mol
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
      type(czt_molecule_t), intent(in) :: mol
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
      type(czt_molecule_t), intent(in) :: mol
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
      type(czt_molecule_t), intent(in) :: mol
      type(error_t), intent(inout) :: error

      integer, allocatable :: numbers(:)

      if (ctx%nlc_grid%n_points > 0) return
      allocate (numbers(mol%natm))
      ! The *element*, not the charge it presents -- see the exchange grid above.
      ! The element itself, when the molecule recorded it: `charges` is zero
      ! on a ghost atom -- a counterpoise partner, or a nucleus quantised by
      ! NEO -- and a grid built for "element 0" there integrates a density
      ! that is anything but empty. The charge-plus-core sum is the fallback
      ! for a molecule assembled without the record.
      if (allocated(mol%atomic_numbers)) then
         numbers = mol%atomic_numbers
      else
         numbers = nint(mol%charges) + mol%core_electrons
      end if
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

      real(dp), allocatable :: scaled(:, :), half(:, :)
      integer :: nb, nao, mu, nu, ig, id

      nb = size(ao, 1)
      nao = size(ao, 2)
      allocate (scaled(nb, nao))

      ! The density term is `(w v_rho chi)^T chi`, symmetric, and the gradient
      ! term is
      !
      !     V_uv += sum_g w_g dE/dgrad rho . (grad chi_u chi_v + chi_u grad chi_v)
      !
      ! whose two halves are transposes of each other. Both are therefore
      ! `M + M^T` for one `M = (w v_rho chi / 2 + w dE/dgrad rho . grad chi)^T chi`,
      ! so one gemm and a transposed add replace what were three gemms --
      ! and this is the assembly every Kohn-Sham Fock build and every kernel
      ! application ends in. Scaling the left factor rather than forming a
      ! diagonal matrix keeps it one multiply per element.
      if (any_gga) then
         do mu = 1, nao
            do ig = 1, nb
               scaled(ig, mu) = 0.5_dp*weights(ig)*vrho(ig)*ao(ig, mu)
               do id = 1, 3
                  scaled(ig, mu) = scaled(ig, mu) &
                                   + weights(ig)*grad_coeff(ig, id)*ao_grad(ig, mu, id)
               end do
            end do
         end do
         allocate (half(nao, nao))
         call pic_gemm(scaled, ao, half, transa="T", beta=0.0_dp)
         do nu = 1, nao
            do mu = 1, nao
               v(mu, nu) = v(mu, nu) + half(mu, nu) + half(nu, mu)
            end do
         end do
         deallocate (half)
      else
         do mu = 1, nao
            do ig = 1, nb
               scaled(ig, mu) = weights(ig)*vrho(ig)*ao(ig, mu)
            end do
         end do
         call pic_gemm(scaled, ao, v, transa="T", alpha=1.0_dp, beta=1.0_dp)
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

end module mqc_czt_xc
