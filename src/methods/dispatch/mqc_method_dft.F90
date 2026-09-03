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
                              BACKEND_CUEST, BACKEND_CZT
   use mqc_cuest_bridge, only: run_cuest_scf
   use mqc_czt_bridge, only: run_czt_hf
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
      character(len=8) :: dispersion_type = "d3bj"
         !! Dispersion type: "d3", "d3bj", "d4"
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

      if (this%options%use_dispersion) then
         ! Refused rather than dropped: a missing dispersion correction biases
         ! every fragment energy the same way.
         call result%error%set(ERROR_VALIDATION, &
                               "Empirical dispersion ("//trim(this%options%dispersion_type)// &
                               ") is not implemented")
         result%has_error = .true.
         return
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
         call run_cuest_scf(settings, fragment, result, want_gradient)
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
         call run_cuest_scf(settings, fragment, result, want_gradient)
#else
         call run_czt_hf(settings, fragment, result, want_gradient, want_hessian)
#endif
      end select
   end subroutine dft_run

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

      call dft_run(this, fragment, result, want_gradient=.true., want_hessian=.true.)
      if (result%has_error) return
      if (result%has_hessian) return

      ! Declined. The finite-difference path runs its own reference point, so
      ! nothing computed above is reusable.
      call result%destroy()
      call finite_difference_hessian(this, fragment, result, verbose=this%options%verbose, &
                                     displacement_in=this%options%hessian_displacement)
   end subroutine dft_calc_hessian

end module mqc_method_dft
