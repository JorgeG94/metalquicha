!! Exchange-correlation, as one thing the SCF can be handed
module mqc_libcint_xc
   !! Everything the Kohn-Sham SCF needs that Hartree-Fock does not, behind one
   !! derived type and one call.
   !!
   !! The shape is deliberate. `run_libcint_rhf` already takes a dozen arguments,
   !! and exchange-correlation needs a grid, a functional, its components, their
   !! weights, an exact-exchange fraction and a running energy -- six more would
   !! have made the signature unreadable and every future addition worse. So the
   !! SCF gains exactly one optional argument of type `xc_context_t`, and one call
   !! inside its Fock build. Hartree-Fock is then the case where that argument is
   !! absent, rather than a separate path, and there is no second SCF loop to keep
   !! in step with the first.
   !!
   !! What the context owns: the grid, the resolved libxc functionals and their
   !! weights, and the exact-exchange fraction. What it deliberately does not own:
   !! basis function values. Those are rebuilt per block of grid points on every
   !! iteration rather than cached, because cached they are `n_points` by `n_ao` --
   !! 6 MB on water/cc-pVDZ but hundreds at a real size, and the coupled-cluster
   !! work already paid once for retrofitting blocking onto something that assumed
   !! it could hold everything.
   !!
   !! LDA, GGA and meta-GGA are all evaluated, restricted and spin-polarised, and
   !! each rung was added inside the same block loop without the interface moving --
   !! which was the bet this shape was making. There are two evaluators rather than
   !! one, `xc_add_potential` and `xc_add_potential_uks`, because libxc fixes the
   !! spin channel when a functional is initialised and the polarised arrays are
   !! spin-interleaved: the two cases genuinely take different array shapes. What
   !! they do share -- the three terms that turn a pointwise potential into a matrix
   !! -- is factored into `accumulate_xc_matrix`, since that is the arithmetic worth
   !! having exactly one copy of.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_dft_grid, only: dft_grid_t, build_dft_grid, DEFAULT_GRID_LEVEL
   use mqc_xc_spec, only: xc_spec_t, xc_spec_from_name, MAX_XC_COMPONENTS
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_libcint_ao, only: eval_ao_block, eval_rho, AO_POINT_BLOCK
#ifdef MQC_WITH_LIBXC
   use xc_f03_lib_m, only: xc_f03_func_t, xc_f03_func_init, xc_f03_func_end, &
                           xc_f03_lda_exc_vxc, xc_f03_func_get_info, &
                           xc_f03_func_info_get_family, xc_f03_hyb_exx_coef, &
                           xc_f03_hyb_cam_coef, xc_f03_nlc_coef, &
                           xc_f03_functional_get_number, xc_f03_func_info_t, &
                           xc_f03_gga_exc_vxc, xc_f03_mgga_exc_vxc, &
                           xc_f03_func_info_get_flags, XC_FLAGS_NEEDS_LAPLACIAN, &
                           XC_UNPOLARIZED, XC_POLARIZED, &
                           XC_FAMILY_LDA, XC_FAMILY_HYB_LDA, &
                           XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, &
                           XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA
#endif
   implicit none
   private

   public :: xc_context_t
   public :: xc_context_create
   public :: xc_add_potential
   public :: xc_add_potential_uks
   public :: xc_available
   public :: xc_grid_lda_quantities
   public :: xc_grid_gga_quantities

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
         !! libxc reports (omega, alpha, beta) with alpha the long-range coefficient
         !! and alpha+beta the short-range one, so exx_fraction is alpha+beta and
         !! this is -beta. Checked against PySCF's resolved values rather than read
         !! off the convention: CAM-B3LYP comes out 0.19 short-range and 0.65 long,
         !! and wB97X 0.157706 and 1.0, which are the published numbers.
      real(dp) :: pt2_fraction = 0.0_dp
         !! MP2 correlation fraction, for a double hybrid. Carried here so the
         !! caller can see it without re-parsing the name; nothing in this module
         !! acts on it, because perturbative correlation is not a grid quantity.
      logical :: polarized = .false.
         !! Whether the functionals were initialised spin-polarised.
         !!
         !! Not a detail a caller can be indifferent to: libxc fixes this when a
         !! functional is initialised, and the two cases take different array
         !! shapes and return a different number of potential components. So a
         !! context belongs to a restricted calculation or to an unrestricted one,
         !! and `xc_add_potential` and `xc_add_potential_uks` each refuse the other
         !! kind rather than reinterpreting the arrays -- which would be a silent
         !! read of the wrong element, not a crash.
      type(dft_grid_t) :: grid
      integer :: n_func = 0
      integer :: family(MAX_XC_COMPONENTS) = 0
         !! libxc's family per component. A composition may mix them -- a hybrid
         !! GGA exchange beside an LDA correlation is ordinary -- so the block loop
         !! asks per component rather than once for the whole functional.
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

   subroutine xc_context_create(mol, functional, ctx, error, level, polarized)
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

      call xc_spec_from_name(functional, spec, error)
      if (error%has_error()) return

      grid_level = DEFAULT_GRID_LEVEL
      if (present(level)) grid_level = level

      allocate (numbers(mol%natm))
      numbers = nint(mol%charges)
      call build_dft_grid(mol%coords, numbers, ctx%grid, error, level=grid_level)
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
         ! coefficients: a global hybrid reports omega = 0, and anything else
         ! splits its exchange into short and long range over an erf-attenuated
         ! kernel. Detecting it matters either way, because `hyb_exx_coef` does
         ! not mean for such a functional what it means for a global one --
         ! taking it at face value gave CAM-B3LYP an energy 3.4 Hartree low and
         ! wB97X 6.4 low, both converged and neither flagged.
         !
         ! The attenuated integrals need no new integral code: they are the same
         ! entry points with `env(PTR_RANGE_OMEGA)` set, so a second exchange
         ! pass and these coefficients are the whole of it.
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

         ! Non-local correlation likewise: VV10 is a double integral over the
         ! density, not a functional of it at a point, and a functional carrying it
         ! would silently lose that whole term.
         call xc_f03_nlc_coef(ctx%func(i), nlc_b, nlc_c)
         if (nlc_b /= 0.0_dp .or. nlc_c /= 0.0_dp) then
            call error%set(ERROR_VALIDATION, "'"//trim(spec%component(i)%name)// &
                           "' carries a non-local correlation term (VV10), which the "// &
                           "CPU path does not implement. Refused rather than evaluated "// &
                           "without it.")
            return
         end if

         ! libxc owns a global hybrid's fraction, so ask rather than assume. Only a
         ! composition libxc does not carry may state its own, and `mqc_xc_spec`
         ! leaves that at zero for everything libxc knows -- so the two can never
         ! both be nonzero and disagree.
         if (spec%from_libxc .and. .not. (cam_omega /= 0.0_dp .and. cam_beta /= 0.0_dp)) then
            libxc_exx = xc_f03_hyb_exx_coef(ctx%func(i))
            ctx%exx_fraction = ctx%exx_fraction + spec%component(i)%weight*libxc_exx
         end if
      end do
#endif

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
      !! The gradient needs these three on the grid rather than contracted into
      !! a matrix: the exchange-correlation gradient has one term contracting
      !! `vrho` against the moving basis functions and another contracting
      !! `rho*exc` against the moving quadrature weights, and neither can be
      !! recovered from the Fock-matrix contribution that `xc_add_potential`
      !! returns.
      !!
      !! LDA only, and it refuses anything else rather than returning something
      !! incomplete. A GGA needs `vsigma` and the density gradient as well, and
      !! its gradient needs second derivatives of the basis functions, which
      !! `eval_ao_block` does not produce.
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
      do g0 = 1, npts, AO_POINT_BLOCK
         g1 = min(g0 + AO_POINT_BLOCK - 1, npts)
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
                                     vtau)
      !! rho, eps_xc, dE/drho and dE/d(grad rho) at every grid point, for a GGA
      !!
      !! The counterpart of `xc_grid_lda_quantities` for a functional that reads
      !! the density gradient, and it returns one thing that routine does not:
      !! `grad_coeff`, the quantity conjugate to `grad rho`.
      !!
      !! **It returns that rather than `vsigma` and `rho_grad` separately, on
      !! purpose.** libxc parametrises a GGA by sigma = |grad rho|^2, so what
      !! multiplies d(grad rho)/dR is 2 vsigma grad rho -- and in the polarised
      !! case the spins couple through sigma_ab, giving
      !!
      !!     dE/d(grad rho_a) = 2 vsigma_aa grad rho_a + vsigma_ab grad rho_b
      !!
      !! Resolving that here means the gradient contracts a single vector per
      !! spin and never has to know how many sigma channels there were. It is
      !! also the same combination `xc_add_potential` builds for the Fock
      !! matrix, so the two paths cannot disagree about the chain rule.
      !!
      !! LDA functionals are accepted and return `grad_coeff` zero, so a caller
      !! that has already decided to take the GGA path does not need a second
      !! branch. Meta-GGA is refused: tau brings a term of its own that nothing
      !! here computes.
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
         !! rather than evaluated without it -- see the check below, and note
         !! that the dispatch further down would otherwise treat one as an LDA.

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
      do g0 = 1, npts, AO_POINT_BLOCK
         g1 = min(g0 + AO_POINT_BLOCK - 1, npts)
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

   subroutine xc_add_potential(ctx, mol, density, v_xc, e_xc, n_elec, error)
      !! The exchange-correlation potential and energy for one density
      !!
      !!     E_xc = sum_g w_g rho_g eps_xc(rho_g)
      !!     V_uv = sum_g w_g v_xc(rho_g) chi_u(r_g) chi_v(r_g)
      !!
      !! `v_xc` comes back on its own rather than added into a Fock matrix, because
      !! the energy expression needs the Fock matrix *without* it: the Kohn-Sham
      !! energy takes E_xc directly, not half the trace of D V_xc, so a caller that
      !! could not separate them would have to reconstruct one or the other.
      !!
      !! `n_elec` is the integrated density, returned because it costs nothing and
      !! is the cheapest possible statement that the grid and the density agree --
      !! a caller printing it once per run turns a silent grid problem into an
      !! obvious one.
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
      real(dp), allocatable :: tau(:), vtau(:), vtau_i(:), lapl(:), vlapl(:)
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
      do g0 = 1, ctx%grid%n_points, AO_POINT_BLOCK
         g1 = min(g0 + AO_POINT_BLOCK - 1, ctx%grid%n_points)
         nb = g1 - g0 + 1

         if (ctx%any_mgga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, grad=ao_grad)
            if (error%has_error()) return
            call eval_rho(ao, density, rho, ao_grad=ao_grad, rho_grad=rho_grad, tau=tau)
         else if (ctx%any_gga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, grad=ao_grad)
            if (error%has_error()) return
            call eval_rho(ao, density, rho, ao_grad=ao_grad, rho_grad=rho_grad)
         else
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error)
            if (error%has_error()) return
            call eval_rho(ao, density, rho)
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
         ! than one is a composition this repository defines, and a double hybrid
         ! is the only such case in view.
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

         n_elec = n_elec + sum(ctx%grid%weights(g0:g1)*rho)
         e_xc = e_xc + sum(ctx%grid%weights(g0:g1)*rho*exc)

         ! The gradient coefficient is dE/d(grad rho). Differentiating
         ! sigma = |grad rho|^2 gives 2 vsigma grad rho; the unrestricted case has
         ! a cross-spin term as well, which is the only difference between the two
         ! and the reason this is assembled by the caller rather than inside.
         if (ctx%any_gga) then
            if (allocated(grad_coeff)) deallocate (grad_coeff)
            allocate (grad_coeff(nb, 3))
            do id = 1, 3
               do ig = 1, nb
                  grad_coeff(ig, id) = 2.0_dp*vsigma(ig)*rho_grad(ig, id)
               end do
            end do
         end if

         call accumulate_xc_matrix(ctx%grid%weights(g0:g1), ao, vrho, v_xc, &
                                   ao_grad=ao_grad, grad_coeff=grad_coeff, &
                                   vtau=vtau, any_gga=ctx%any_gga, &
                                   any_mgga=ctx%any_mgga)
      end do
#endif
   end subroutine xc_add_potential

   subroutine xc_add_potential_uks(ctx, mol, d_alpha, d_beta, v_alpha, v_beta, &
                                   e_xc, n_elec, error)
      !! The exchange-correlation potential and energy for a pair of spin densities
      !!
      !! The unrestricted counterpart of `xc_add_potential`, and the same structure:
      !! one pass over blocks of grid points, libxc per component, one matrix
      !! assembled per spin. Two things are genuinely different, and they are the
      !! two places an unrestricted functional goes wrong.
      !!
      !! **libxc's arrays are spin-interleaved.** A polarised functional takes
      !! `rho` as (rho_a, rho_b) per point, `sigma` as (sigma_aa, sigma_ab,
      !! sigma_bb) and `tau` as (tau_a, tau_b), all with spin as the fastest index,
      !! and returns `vrho` and `vtau` in the same layout with `vsigma` in three
      !! components. Nothing in the interface enforces that -- the F03 bindings take
      !! bare arrays -- so a wrong stride here is a silent misread, which is why the
      !! closed-shell-equals-restricted test exists.
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
      real(dp), allocatable :: rho_a(:), rho_b(:), grad_a(:, :), grad_b(:, :)
      real(dp), allocatable :: tau_a(:), tau_b(:)
      real(dp), allocatable :: rho(:), sigma(:), tau(:), lapl(:)
      real(dp), allocatable :: exc(:), vrho(:), vsigma(:), vtau(:)
      real(dp), allocatable :: exc_i(:), vrho_i(:), vsigma_i(:), vtau_i(:), vlapl(:)
      real(dp), allocatable :: vrho_s(:), vtau_s(:), grad_coeff(:, :)
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
      do g0 = 1, ctx%grid%n_points, AO_POINT_BLOCK
         g1 = min(g0 + AO_POINT_BLOCK - 1, ctx%grid%n_points)
         nb = g1 - g0 + 1

         ! One AO evaluation for both spins -- the expensive part -- and then the
         ! density contraction once per spin. `eval_rho` takes a density matrix and
         ! asks nothing about whose it is, so a spin density needs no new code.
         if (ctx%any_mgga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, grad=ao_grad)
            if (error%has_error()) return
            call eval_rho(ao, d_alpha, rho_a, ao_grad=ao_grad, rho_grad=grad_a, tau=tau_a)
            call eval_rho(ao, d_beta, rho_b, ao_grad=ao_grad, rho_grad=grad_b, tau=tau_b)
         else if (ctx%any_gga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, grad=ao_grad)
            if (error%has_error()) return
            call eval_rho(ao, d_alpha, rho_a, ao_grad=ao_grad, rho_grad=grad_a)
            call eval_rho(ao, d_beta, rho_b, ao_grad=ao_grad, rho_grad=grad_b)
         else
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error)
            if (error%has_error()) return
            call eval_rho(ao, d_alpha, rho_a)
            call eval_rho(ao, d_beta, rho_b)
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
         n_elec = n_elec + sum(ctx%grid%weights(g0:g1)*(rho_a + rho_b))
         e_xc = e_xc + sum(ctx%grid%weights(g0:g1)*(rho_a + rho_b)*exc)

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
         call accumulate_xc_matrix(ctx%grid%weights(g0:g1), ao, vrho_s, v_alpha, &
                                   ao_grad=ao_grad, grad_coeff=grad_coeff, &
                                   vtau=vtau_s, any_gga=ctx%any_gga, &
                                   any_mgga=ctx%any_mgga)

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
         call accumulate_xc_matrix(ctx%grid%weights(g0:g1), ao, vrho_s, v_beta, &
                                   ao_grad=ao_grad, grad_coeff=grad_coeff, &
                                   vtau=vtau_s, any_gga=ctx%any_gga, &
                                   any_mgga=ctx%any_mgga)
      end do
#endif
   end subroutine xc_add_potential_uks

   subroutine accumulate_xc_matrix(weights, ao, vrho, v, ao_grad, grad_coeff, vtau, &
                                   any_gga, any_mgga)
      !! One spin's exchange-correlation matrix, over one block of grid points
      !!
      !! Shared by the restricted and unrestricted paths, which is the whole point:
      !! the three terms below are where a functional is wrong by millihartree
      !! while converging perfectly, and there should be one copy of them. What
      !! differs between the two callers is not this arithmetic but what goes into
      !! `vrho`, `grad_coeff` and `vtau` -- for a spin-polarised functional those
      !! are that spin's derivatives, including the cross-spin gradient term, and
      !! this routine neither knows nor needs to.
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
