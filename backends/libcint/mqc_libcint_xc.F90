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
   !! Only LDA is evaluated so far. GGA needs the AO gradients and a second term in
   !! the potential; the block loop below is where both go, and nothing about the
   !! interface changes when they arrive.
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
                           XC_UNPOLARIZED, XC_FAMILY_LDA, XC_FAMILY_HYB_LDA, &
                           XC_FAMILY_GGA, XC_FAMILY_HYB_GGA, &
                           XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA
#endif
   implicit none
   private

   public :: xc_context_t
   public :: xc_context_create
   public :: xc_add_potential
   public :: xc_available

   type :: xc_context_t
      !! A functional, a grid, and the fraction of exact exchange to keep
      logical :: active = .false.
      real(dp) :: exx_fraction = 0.0_dp
         !! How much Fock exchange the functional wants. Zero for pure DFT, one
         !! for Hartree-Fock, in between for a hybrid. Read from libxc where libxc
         !! owns the functional, and from `mqc_xc_spec` only where it does not.
      real(dp) :: pt2_fraction = 0.0_dp
         !! MP2 correlation fraction, for a double hybrid. Carried here so the
         !! caller can see it without re-parsing the name; nothing in this module
         !! acts on it, because perturbative correlation is not a grid quantity.
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

   subroutine xc_context_create(mol, functional, ctx, error, level)
      !! Resolve a functional name and build the grid it will be integrated on
      type(libcint_molecule_t), intent(in) :: mol
      character(len=*), intent(in) :: functional
      type(xc_context_t), intent(out) :: ctx
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: level

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

#ifdef MQC_WITH_LIBXC
      do i = 1, spec%n_components
         ctx%weight(i) = spec%component(i)%weight
         id = xc_f03_functional_get_number(trim(spec%component(i)%name))
         if (id <= 0) then
            call error%set(ERROR_VALIDATION, "libxc does not know the functional '"// &
                           trim(spec%component(i)%name)//"'")
            return
         end if
         call xc_f03_func_init(ctx%func(i), id, XC_UNPOLARIZED)

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

         ! **Range separation is refused, and has to be detected rather than
         ! asked about.** libxc 7.1.2's Fortran bindings expose no `xc_hyb_type`,
         ! so the test is on the coefficients themselves: a global hybrid reports
         ! omega = 0, and anything else splits its exchange into short and long
         ! range with an erf-attenuated kernel. Those integrals are not computed
         ! here, and `hyb_exx_coef` does not mean for such a functional what it
         ! means for a global one -- taking it at face value gave CAM-B3LYP an
         ! energy 3.4 Hartree low and wB97X 6.4 low, both converged and neither
         ! flagged. Exactly what this backend exists to refuse.
         call xc_f03_hyb_cam_coef(ctx%func(i), cam_omega, cam_alpha, cam_beta)
         if (cam_omega /= 0.0_dp .or. cam_beta /= 0.0_dp) then
            call error%set(ERROR_VALIDATION, "'"//trim(spec%component(i)%name)// &
                           "' is a range-separated hybrid, which needs "// &
                           "erf-attenuated two-electron integrals the CPU path does "// &
                           "not compute. Refused rather than treated as a global "// &
                           "hybrid, which would converge to a badly wrong energy.")
            return
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
         if (spec%from_libxc) then
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
   end subroutine xc_context_destroy

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

      real(dp), allocatable :: ao(:, :), rho(:), exc(:), vrho(:), scaled(:, :)
      real(dp), allocatable :: exc_i(:), vrho_i(:)
      real(dp), allocatable :: ao_grad(:, :, :), rho_grad(:, :), sigma(:)
      real(dp), allocatable :: vsigma(:), vsigma_i(:)
      real(dp), allocatable :: tau(:), vtau(:), vtau_i(:), lapl(:), vlapl(:)
      integer :: g0, g1, nb, i, mu, ig, id

      v_xc = 0.0_dp
      e_xc = 0.0_dp
      n_elec = 0.0_dp

      if (.not. ctx%active) return
      if (.not. xc_available()) then
         call error%set(ERROR_VALIDATION, "no libxc in this build")
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

         ! V += (w v_xc chi)^T chi, as a gemm. Scaling the left factor rather than
         ! forming a diagonal matrix keeps this one multiply per element.
         allocate (scaled(nb, mol%nao))
         do mu = 1, mol%nao
            do ig = 1, nb
               scaled(ig, mu) = ctx%grid%weights(g0 + ig - 1)*vrho(ig)*ao(ig, mu)
            end do
         end do
         call pic_gemm(scaled, ao, v_xc, transa="T", alpha=1.0_dp, beta=1.0_dp)

         ! The gradient term. From differentiating sigma = |grad rho|^2 with
         ! grad rho = 2 sum_uv D_uv chi_v grad chi_u,
         !
         !     V_uv += sum_g w_g 2 vsigma_g grad rho_g . (grad chi_u chi_v
         !                                                + chi_u grad chi_v)
         !
         ! and the two halves are transposes of each other, so one gemm plus a
         ! symmetrisation does both. The factor of two and that symmetrisation are
         ! the usual place a GGA is wrong by millihartree while converging
         ! perfectly, which is why the finite-difference test on the gradients
         ! exists separately from any energy.
         if (ctx%any_gga) then
            do mu = 1, mol%nao
               do ig = 1, nb
                  scaled(ig, mu) = 0.0_dp
                  do id = 1, 3
                     scaled(ig, mu) = scaled(ig, mu) &
                                      + 2.0_dp*ctx%grid%weights(g0 + ig - 1) &
                                      *vsigma(ig)*rho_grad(ig, id)*ao_grad(ig, mu, id)
                  end do
               end do
            end do
            ! scaled^T ao gives one half; adding its transpose gives the other.
            call pic_gemm(scaled, ao, v_xc, transa="T", alpha=1.0_dp, beta=1.0_dp)
            call pic_gemm(ao, scaled, v_xc, transa="T", alpha=1.0_dp, beta=1.0_dp)
         end if

         ! The kinetic-energy-density term. d tau / d D_uv is half the sum over
         ! directions of grad chi_u grad chi_v, which is already symmetric in u and
         ! v -- so unlike the sigma term this one needs no transpose added.
         if (ctx%any_mgga) then
            do id = 1, 3
               do mu = 1, mol%nao
                  do ig = 1, nb
                     scaled(ig, mu) = 0.5_dp*ctx%grid%weights(g0 + ig - 1) &
                                      *vtau(ig)*ao_grad(ig, mu, id)
                  end do
               end do
               call pic_gemm(scaled, ao_grad(:, :, id), v_xc, transa="T", &
                             alpha=1.0_dp, beta=1.0_dp)
            end do
         end if
         deallocate (scaled)
      end do
#endif
   end subroutine xc_add_potential

end module mqc_libcint_xc
