!! Density Functional Theory (DFT) method implementation for metalquicha
module mqc_method_dft
   !! Kohn-Sham density functional theory.
   !!
   !! Dispatches to whichever backend `model.backend` names, exactly as
   !! `mqc_method_hf` does; the CPU path is `run_czt_hf`, which takes the
   !! functional as one more setting on the same SCF.
   !!
   !! Density fitting is not optional on cuEST: it has no conventional
   !! four-index ERI path, so J (and K for hybrids) are always fitted and an
   !! auxiliary JKFIT basis is required. `density_fitting` is ignored by that
   !! backend rather than switchable. The amount of exact exchange is queried
   !! from cuEST's XC plan rather than assumed from the functional name, so a
   !! hybrid cannot end up with mismatched Coulomb and XC definitions.
   use pic_types, only: dp
   use mqc_scf_types, only: guess_step_t
   use mqc_method_config, only: scf_options_t, pcm_config_t, properties_config_t
   use mqc_method_base, only: qc_method_t
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_semi_numerical_hessian, only: finite_difference_hessian
   use mqc_cuest_iface, only: apply_properties_settings, apply_scf_settings, cuest_scf_settings_t, parse_backend_name, &
                              BACKEND_CUEST, BACKEND_CZT, BACKEND_TERCO
   use mqc_cuest_bridge, only: run_cuest_scf
   use mqc_terco_bridge, only: run_terco_scf
   use mqc_czt_bridge, only: run_czt_hf
   use mqc_dispersion, only: dispersion_correction
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: dft_method_t, dft_options_t

   type, extends(scf_options_t) :: dft_options_t
      !! DFT calculation options
      character(len=32) :: functional = "b3lyp"
         !! Exchange-correlation functional
      character(len=16) :: grid_type = "medium"
         !! Integration grid quality
      integer :: radial_points = 75
         !! Number of radial grid points per atom
      integer :: angular_points = 302
         !! Number of angular grid points (Lebedev order)
      integer :: grid_level = 3
         !! 0 to 9, the standard tables. Three is the usual default.
      integer :: nlc_grid_level = -1
         !! VV10's quadrature level; negative means the backend default.
      real(dp) :: screening_tolerance = 1.0e-12_dp
         !! AO value below which a shell is dropped from a grid block.
      integer :: block_size = -1
         !! Grid points per block; -1 keeps the backend default.
      logical :: use_dispersion = .false.
         !! Add empirical dispersion correction
      character(len=16) :: dispersion_type = "d3bj"
         !! Which correction, in the spelling `keywords.dft.dispersion` uses.
         !! Only "d3bj" today; see `DISPERSION_KINDS` in `mqc_dispersion_names`.
   end type dft_options_t

   type, extends(qc_method_t) :: dft_method_t
      !! Kohn-Sham DFT with a configurable functional and integration grid
      type(dft_options_t) :: options
   contains
      procedure :: calc_energy => dft_calc_energy
      procedure :: calc_gradient => dft_calc_gradient
      procedure :: calc_hessian => dft_calc_hessian
   end type dft_method_t

contains

   subroutine dft_calc_energy(this, fragment, result)
      !! Electronic energy of a fragment
      class(dft_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call dft_run(this, fragment, result, want_gradient=.false.)
   end subroutine dft_calc_energy

   subroutine dft_run(this, fragment, result, want_gradient, want_hessian)
      !! Run the SCF through whichever backend `options%backend` resolves to
      class(dft_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in) :: want_gradient
      logical, intent(in), optional :: want_hessian
         !! Only the CPU backend has an analytic Hessian, and only for some of
         !! what reaches it; the backend decides, not this routine.

      type(cuest_scf_settings_t) :: settings
      type(error_t) :: backend_error
      type(error_t) :: dispersion_error
      real(dp) :: e_dispersion
      real(dp), allocatable :: g_dispersion(:, :)

      ! The empirical dispersion correction, computed before the SCF rather
      ! than after it.
      !
      ! It costs microseconds and depends on nothing but the nuclei, so running
      ! it first buys the refusals for free: an unknown correction, a functional
      ! with no published damping parameters, or a build with no dispersion
      ! library is reported before an SCF is started rather than after one has
      ! converged. The numbers are held until the backend has returned something
      ! to add them to.
      !
      ! Here rather than in a backend, because the correction is a function of
      ! the nuclei alone: cenzontle, cuEST and terco would each need the same
      ! code and could each get it subtly differently.
      e_dispersion = 0.0_dp
      if (this%options%use_dispersion) then
         if (want_gradient) then
            allocate (g_dispersion(3, fragment%n_atoms), source=0.0_dp)
            call dispersion_correction(this%options%dispersion_type, this%options%functional, &
                                       fragment%element_numbers, fragment%coordinates, &
                                       e_dispersion, g_dispersion, dispersion_error)
         else
            call dispersion_correction(this%options%dispersion_type, this%options%functional, &
                                       fragment%element_numbers, fragment%coordinates, &
                                       e_dispersion, error=dispersion_error)
         end if
         if (dispersion_error%has_error()) then
            call result%error%set(ERROR_VALIDATION, dispersion_error%get_message())
            result%has_error = .true.
            return
         end if
      end if

      call apply_scf_settings(settings, this%options)
      settings%functional = this%options%functional
      call parse_backend_name(this%options%backend, settings%backend, backend_error)
      if (backend_error%has_error()) then
         call result%error%set(ERROR_VALIDATION, backend_error%get_message())
         result%has_error = .true.
         return
      end if
      settings%radial_points = this%options%radial_points
      settings%angular_points = this%options%angular_points
      settings%grid_level = this%options%grid_level
      settings%nlc_grid_level = this%options%nlc_grid_level
      settings%screening_tolerance = this%options%screening_tolerance
      settings%block_size = this%options%block_size
      ! TODO(mqc): `grid_type` is not copied here, and is read nowhere in the
      ! tree, although the factory fills it from `dft.grid_type`. A deck naming
      ! it silently gets whatever `grid_level` and the point counts above chose.
      ! The quasi-atomic bonding analysis is refused, not ignored: it is defined
      ! against a Hartree-Fock or MCSCF wavefunction, and `run_czt_hf`
      ! dispatches on the deck naming an analysis alone, so passing it on would
      ! decompose Kohn-Sham orbitals and report numbers for it. Tested on the
      ! name rather than through `bonding_analysis_kind`, which lives in a
      ! backend module this one must not depend on.
      if (len_trim(this%options%properties%bonding_analysis) > 0 .and. &
          trim(adjustl(this%options%properties%bonding_analysis)) /= "none") then
         call result%error%set(ERROR_VALIDATION, "the quasi-atomic bonding analysis is "// &
                               "not available for a Kohn-Sham reference: it is defined "// &
                               "against a Hartree-Fock or MCSCF wavefunction. Request it "// &
                               "on one of those instead.")
         result%has_error = .true.
         return
      end if
      call apply_properties_settings(settings, this%options%properties)

      ! Which backend, and refuse a request that cannot be honoured. Asking for
      ! cuEST on a CPU-only build reaches the stub `run_cuest_scf`, which
      ! reports the missing build rather than falling through to libcint.
      ! TODO(mqc): the MP2/CC refusal below is copied from `hf_run` and is dead
      ! here -- nothing on this path ever sets `run_mp2` or `run_cc`.
      select case (settings%backend)
      case (BACKEND_CUEST)
         if (settings%cartesian) then
            call result%error%set(ERROR_VALIDATION, "backend 'cuest' was asked for, but "// &
                                  "'model.cartesian' is on and the GPU path builds its "// &
                                  "AO shells spherical whatever the basis says. Running "// &
                                  "it would answer with a different basis than the deck "// &
                                  "asked for and say nothing. Ask for backend 'libcint', "// &
                                  "or drop 'model.cartesian'.")
            result%has_error = .true.
            return
         end if
         if (settings%run_mp2 .or. settings%run_cc) then
            call result%error%set(ERROR_VALIDATION, "backend 'cuest' was asked for, but "// &
                                  "MP2 and coupled cluster have no GPU implementation "// &
                                  "here -- they run through the CPU backend. Ask for "// &
                                  "'auto', or drop the correlated method.")
            result%has_error = .true.
            return
         end if
         if (settings%excited%enabled) then
            call result%error%set(ERROR_VALIDATION, "backend 'cuest' was asked for, but "// &
                                  "keywords.excited_states has no GPU implementation "// &
                                  "here -- the linear-response solver lives on the CPU "// &
                                  "backend. Ask for 'auto' or 'libcint', or drop the "// &
                                  "excited states.")
            result%has_error = .true.
            return
         end if
         call run_cuest_scf(settings, fragment, result, want_gradient)
      case (BACKEND_TERCO)
         ! Every refusal terco needs -- gradients, correlated methods,
         ! spherical d, angular momentum above d, an ECP -- is made inside
         ! the driver, against the basis it actually built. Repeating them
         ! here would be a second copy to keep in step.
         call run_terco_scf(settings, fragment, result, want_gradient)
      case (BACKEND_CZT)
         call run_czt_hf(settings, fragment, result, want_gradient, want_hessian)
      case default
#ifdef MQC_WITH_CUEST
         if (settings%cartesian) then
            call result%error%set(ERROR_VALIDATION, "'model.cartesian' is on and this "// &
                                  "build resolves 'auto' to the GPU backend, which "// &
                                  "builds its AO shells spherical whatever the basis "// &
                                  "says. Ask for backend 'libcint', or drop "// &
                                  "'model.cartesian'.")
            result%has_error = .true.
            return
         end if
         ! The same refusal as the explicit request above. Asking for
         ! excited states on a build that resolves 'auto' to the GPU would
         ! otherwise reach `run_cuest_scf`, which knows nothing about the
         ! block and would answer with a ground state and no complaint.
         if (settings%excited%enabled) then
            call result%error%set(ERROR_VALIDATION, "this build resolves 'auto' to the "// &
                                  "GPU backend, which has no linear-response solver, "// &
                                  "and keywords.excited_states asked for roots. Ask "// &
                                  "for backend 'libcint', or drop the excited states.")
            result%has_error = .true.
            return
         end if
         call run_cuest_scf(settings, fragment, result, want_gradient)
#else
         call run_czt_hf(settings, fragment, result, want_gradient, want_hessian)
#endif
      end select

      call add_dispersion(this%options%use_dispersion, e_dispersion, g_dispersion, result)
   end subroutine dft_run

   subroutine add_dispersion(requested, energy, gradient, result)
      !! Fold the dispersion correction into a result the backend has filled in
      !!
      !! Kept beside the total rather than added into `energy%scf`: the
      !! correction is not a functional of the density and converged nothing,
      !! and a total that hides it cannot be compared with a published DFT-D
      !! number or with the same geometry run without it. `energy%total()` adds
      !! it, so every consumer of the total already has it.
      !!
      !! The gradient sign is the library's, unchanged: s-dftd3 returns dE/dR in
      !! Hartree per Bohr, which is what `result%gradient` holds, so the two
      !! simply add. Over the fragment's own atoms, H-caps included -- the same
      !! atoms the SCF gradient covers, so the cap redistribution downstream
      !! sees one consistent gradient rather than two conventions.
      logical, intent(in) :: requested
      real(dp), intent(in) :: energy
      real(dp), allocatable, intent(in) :: gradient(:, :)
      type(calculation_result_t), intent(inout) :: result

      character(len=80) :: line

      if (.not. requested) return
      if (result%has_error) return

      result%energy%dispersion = energy
      if (result%has_gradient .and. allocated(gradient) .and. allocated(result%gradient)) then
         result%gradient = result%gradient + gradient
      end if

      ! On its own line, and unconditionally. A correction of a few millihartree
      ! that appears only inside a total is the kind of thing a run is later
      ! unable to say whether it had.
      write (line, "(a,f20.12)") "  empirical dispersion ", result%energy%dispersion
      call logger%info(trim(line))
   end subroutine add_dispersion

   subroutine dft_calc_gradient(this, fragment, result)
      !! Energy and nuclear gradient of a fragment
      class(dft_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call dft_run(this, fragment, result, want_gradient=.true.)
   end subroutine dft_calc_gradient

   subroutine dft_calc_hessian(this, fragment, result)
      !! The analytic Hessian where there is one, central differences otherwise
      !!
      !! The same shape as `hf_calc_hessian`: the request goes down and
      !! `has_hessian` comes back, true meaning it was computed and false with
      !! no error meaning it was declined and finite differences take over.
      !! Which functionals qualify is a fact about the backend, not about
      !! anything visible here.
      class(dft_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      ! Dispersion goes straight to finite differences, without asking the
      ! backend whether it has an analytic Hessian. s-dftd3's C API returns an
      ! energy, a gradient and a virial and no second derivative, so an analytic
      ! Kohn-Sham Hessian plus this correction would be a Hessian missing one of
      ! its terms with nothing to say so. Differencing `calc_gradient`, which
      ! does carry the correction, gives the whole thing.
      if (.not. this%options%use_dispersion) then
         call dft_run(this, fragment, result, want_gradient=.true., want_hessian=.true.)
         if (result%has_error) return
         if (result%has_hessian) return
      end if

      ! Declined. The finite-difference path runs its own reference point, so
      ! nothing computed above is reusable.
      call result%destroy()
      call finite_difference_hessian(this, fragment, result, verbose=this%options%verbose, &
                                     displacement_in=this%options%hessian_displacement)
   end subroutine dft_calc_hessian

end module mqc_method_dft
