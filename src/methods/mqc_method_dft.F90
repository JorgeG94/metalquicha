!! Density Functional Theory (DFT) method implementation for metalquicha
module mqc_method_dft
   !! Closed-shell Kohn-Sham DFT.
   !!
   !! With the cuEST backend (CMake, MQC_ENABLE_CUEST=ON) the Coulomb and exact-exchange matrices come from
   !! cuEST's density-fitted engine and the exchange-correlation energy and
   !! potential from its grid engine; the SCF itself runs in `mqc_cuest_scf`.
   !!
   !! Density fitting is not optional here. cuEST has no conventional
   !! four-index ERI path, so J (and K for hybrids) are always fitted and an
   !! auxiliary JKFIT basis is required. `use_density_fitting` is therefore
   !! ignored by this backend rather than switchable.
   !!
   !! The amount of exact exchange is never assumed from the functional name:
   !! it is queried from cuEST's XC plan and handed to the DF plan, so a hybrid
   !! cannot end up with mismatched Coulomb and XC definitions.
   use pic_types, only: dp
   use mqc_config_types, only: guess_step_t
   use mqc_method_config, only: scf_options_t, pcm_config_t, properties_config_t
   use mqc_method_base, only: qc_method_t
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_semi_numerical_hessian, only: finite_difference_hessian
   use mqc_cuest_iface, only: apply_properties_settings, apply_scf_settings, cuest_scf_settings_t, parse_backend_name, &
                              BACKEND_CUEST, BACKEND_LIBCINT
   use mqc_cuest_bridge, only: run_cuest_scf
   use mqc_libcint_bridge, only: run_libcint_hf
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
      integer :: grid_level = 3
      integer :: nlc_grid_level = -1
         !! VV10's quadrature level; negative means the backend default.
         !! 0 to 9, the standard tables. Three is the usual default and what a
         !! production calculation should start from.
      real(dp) :: screening_tolerance = 1.0e-12_dp
         !! AO value below which a shell is dropped from a grid block.
      integer :: block_size = -1
         !! Grid points per block; -1 keeps the backend default.
         !! Number of angular grid points (Lebedev)

      ! Density fitting
      logical :: use_dispersion = .false.
         !! Add empirical dispersion correction
      character(len=8) :: dispersion_type = "d3bj"
         !! Dispersion type: "d3", "d3bj", "d4"

      ! DIIS acceleration
   end type dft_options_t

   type, extends(qc_method_t) :: dft_method_t
      !! DFT method implementation
      !!
      !! Kohn-Sham DFT with configurable exchange-correlation functional,
      !! integration grid, and optional density fitting.
      type(dft_options_t) :: options
   contains
      procedure :: calc_energy => dft_calc_energy
      procedure :: calc_gradient => dft_calc_gradient
      procedure :: calc_hessian => dft_calc_hessian
   end type dft_method_t

contains

   subroutine dft_calc_energy(this, fragment, result)
      !! Calculate electronic energy using Kohn-Sham DFT
      class(dft_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call dft_run(this, fragment, result, want_gradient=.false.)
   end subroutine dft_calc_energy

   subroutine dft_run(this, fragment, result, want_gradient, want_hessian)
      !! Run a density-fitted Kohn-Sham calculation through cuEST
      class(dft_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in) :: want_gradient
      logical, intent(in), optional :: want_hessian
         !! Only the CPU backend has an analytic Hessian, and only for some of
         !! what reaches it. Passed along rather than decided here: which cases
         !! qualify is a fact about that backend.

      type(cuest_scf_settings_t) :: settings
      type(error_t) :: backend_error

      if (this%options%use_dispersion) then
         ! Silently dropping a dispersion correction would bias every fragment
         ! energy the same way, which is exactly the kind of error a many-body
         ! expansion amplifies rather than cancels.
         call result%error%set(ERROR_VALIDATION, &
                               "Empirical dispersion ("//trim(this%options%dispersion_type)// &
                               ") is not implemented")
         result%has_error = .true.
         return
      end if

      call apply_scf_settings(settings, this%options)
      settings%functional = this%options%functional
      ! Resolved here rather than carried as a string, so an unknown name fails
      ! once, before any integrals, instead of at each dispatch.
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
      ! Where the molecule reacts. Absent here until now, so a DFT deck asking
      ! for `properties.fukui` got a normal DFT run and no analysis, with no
      ! error to say why: the bridge gates on this being allocated.
      ! Refused rather than run. The QUAO bonding options are not meant for a
      ! Kohn-Sham reference, and nothing in the backend stops them: the
      ! dispatch is gated only on the deck naming an analysis, so handing these
      ! over would run a decomposition on Kohn-Sham orbitals and report numbers
      ! for it. Before this refactor a Kohn-Sham deck asking for them was
      ! silently ignored, which is the failure this work exists to end -- but
      ! the fix is to say no, not to start obeying a request the method was
      ! never meant to serve. `bonding_analysis` itself stays copied: it selects
      ! *which* analysis and its own dispatch handles "none".
      if (this%options%properties%bonding_energy .or. &
          this%options%properties%bonding_no_sharing .or. &
          this%options%properties%bonding_restrict_localization) then
         call result%error%set(ERROR_VALIDATION, "the QUAO bonding options are not "// &
                               "available for a Kohn-Sham reference: the energy "// &
                               "decomposition and its sharing controls are defined "// &
                               "against a Hartree-Fock or MCSCF wavefunction. Request "// &
                               "them on one of those instead.")
         result%has_error = .true.
         return
      end if
      call apply_properties_settings(settings, this%options%properties)
      ! Set unconditionally, and deliberately not guarded on the backend. cuEST
      ! has no four-index path so it fits regardless and ignores this, per the
      ! note at the top of this module; on the libcint side it is a real choice.
      ! Guarding on `BACKEND_LIBCINT` was tried and was wrong: an unspecified
      ! backend parses to `BACKEND_AUTO` and is only resolved further down, so
      ! the guard was false for every ordinary deck and the flag never arrived.
      !
      ! Until this line existed a Kohn-Sham deck could not turn fitting on
      ! however it asked. That made the fitted path unreachable rather than
      ! wrong, which is the better of the two failures but still a gap.

      ! The same choice Hartree-Fock makes, and made the same way rather than a
      ! second time: cuEST when this build has it, because that is the production
      ! path and the only one with gradients, and the CPU backend otherwise --
      ! which is how a fragmented Kohn-Sham calculation runs on a laptop, and what
      ! gives the GPU path a second implementation to disagree with.
      !
      ! Kohn-Sham reaches the CPU backend through `run_libcint_hf` rather than a
      ! routine of its own. That is deliberate: on that side a functional is one
      ! optional argument to the same SCF, so a separate entry point would be two
      ! names for one code path and two places to keep the settings in step.
      ! Which backend, and refuse a request that cannot be honoured.
      !
      ! `run_cuest_scf` exists either way -- the stub compiled without the
      ! backend reports the missing build and computes nothing -- so asking for
      ! cuEST on a CPU-only build produces that message rather than silently
      ! running on the CPU. That is the point of naming a backend: a deck that
      ! said `cuest` and got libcint would report a provenance that was not true.
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
      case (BACKEND_LIBCINT)
         call run_libcint_hf(settings, fragment, result, want_gradient, want_hessian)
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
         call run_libcint_hf(settings, fragment, result, want_gradient, want_hessian)
#endif
      end select
   end subroutine dft_run

   subroutine dft_calc_gradient(this, fragment, result)
      !! Calculate energy gradient using Kohn-Sham DFT
      class(dft_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call dft_run(this, fragment, result, want_gradient=.true.)
   end subroutine dft_calc_gradient

   subroutine dft_calc_hessian(this, fragment, result)
      !! The analytic Hessian where there is one, central differences otherwise
      !!
      !! The same shape as `hf_calc_hessian`, and for the same reason: whether
      !! an analytic Hessian applies depends on what the backend knows -- which
      !! rung the functional sits on, whether it carries VV10, whether the SCF
      !! converged restricted and unfitted -- and not on anything visible here.
      !! The request goes down and `has_hessian` comes back: true means it was
      !! computed, false with no error means it was declined and finite
      !! differences take over.
      !!
      !! Declining is not failing. A backend that cannot honour a *gradient*
      !! refuses loudly, because there is no other way to get one; a Hessian
      !! has a second way that is correct and merely slower.
      class(dft_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call dft_run(this, fragment, result, want_gradient=.true., want_hessian=.true.)
      if (result%has_error) return
      if (result%has_hessian) return

      ! Declined. Nothing computed above is reusable -- the finite-difference
      ! path runs its own reference point -- so this starts over rather than
      ! trying to keep the SCF that just ran.
      call result%destroy()
      call finite_difference_hessian(this, fragment, result, verbose=this%options%verbose)
   end subroutine dft_calc_hessian

end module mqc_method_dft
