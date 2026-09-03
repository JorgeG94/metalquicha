!! Hartree-Fock method implementation for metalquicha
module mqc_method_hf
   !! Hartree-Fock, and the MP2 and coupled-cluster corrections built on it.
   !!
   !! Dispatches to whichever backend `model.backend` names: the CPU one
   !! through `mqc_libcint_bridge`, or the GPU one through `mqc_cuest_bridge`.
   !! cuEST has no conventional four-index ERI path, so J and K are always
   !! density-fitted there and an auxiliary (JKFIT) basis is required alongside
   !! the orbital basis.
   use pic_types, only: dp
   use mqc_scf_types, only: guess_step_t
   use mqc_method_base, only: qc_method_t
   use mqc_method_config, only: scf_options_t, pcm_config_t, properties_config_t
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

   public :: hf_method_t, hf_options_t

   type, extends(scf_options_t) :: hf_options_t
      !! Hartree-Fock calculation options
      logical :: run_mp2 = .false.
         !! Follow the reference with an MP2 correction. Set by the factory
         !! from the method name, not by a deck keyword.
      logical :: corr_density_fitting = .false.
      real(dp) :: scs_ss = 1.0_dp
      real(dp) :: scs_os = 1.0_dp
      logical :: run_cc = .false.
         !! Follow the reference with coupled cluster. Set by the factory from
         !! the method name, not by a deck keyword.
      logical :: cc_triples = .false.
      integer :: cc_max_iter = 100
      real(dp) :: cc_tolerance = 1.0e-8_dp
      integer :: cc_diis_size = 8
      logical :: cc_spin_adapted = .true.
         !! Spatial-orbital coupled cluster rather than spin orbitals
   end type hf_options_t

   type, extends(qc_method_t) :: hf_method_t
      !! Hartree-Fock method implementation
      type(hf_options_t) :: options
   contains
      procedure :: calc_energy => hf_calc_energy
      procedure :: calc_gradient => hf_calc_gradient
      procedure :: calc_hessian => hf_calc_hessian
   end type hf_method_t

contains

   subroutine hf_calc_energy(this, fragment, result)
      !! Electronic energy of a fragment
      class(hf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call hf_run(this, fragment, result, want_gradient=.false.)
   end subroutine hf_calc_energy

   subroutine hf_run(this, fragment, result, want_gradient, want_hessian)
      !! Run the SCF through whichever backend `options%backend` resolves to
      class(hf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in) :: want_gradient
      logical, intent(in), optional :: want_hessian
         !! Only the CPU backend has an analytic Hessian, and only for some of
         !! what reaches it; the backend decides, not this routine.

      type(cuest_scf_settings_t) :: settings
      type(error_t) :: backend_error

      call apply_scf_settings(settings, this%options)
      call apply_properties_settings(settings, this%options%properties)
      settings%run_mp2 = this%options%run_mp2
      settings%corr_density_fitting = this%options%corr_density_fitting
      settings%run_cc = this%options%run_cc
      settings%cc_triples = this%options%cc_triples
      settings%cc_max_iter = this%options%cc_max_iter
      settings%cc_tolerance = this%options%cc_tolerance
      settings%cc_diis_size = this%options%cc_diis_size
      settings%cc_spin_adapted = this%options%cc_spin_adapted
      settings%scs_ss = this%options%scs_ss
      settings%scs_os = this%options%scs_os
      settings%functional = ""        ! empty selects pure Hartree-Fock
      call parse_backend_name(this%options%backend, settings%backend, backend_error)
      if (backend_error%has_error()) then
         call result%error%set(ERROR_VALIDATION, backend_error%get_message())
         result%has_error = .true.
         return
      end if

      ! Which backend, and refuse a request that cannot be honoured. Asking for
      ! cuEST on a CPU-only build reaches the stub `run_cuest_scf`, which
      ! reports the missing build rather than falling through to libcint.
      ! TODO(mqc): the MP2/CC refusal below guards BACKEND_CUEST only. On a
      ! cuEST build `backend: auto` falls into the default branch and runs
      ! `run_cuest_scf` for an MP2 or CC deck too -- which is exactly what the
      ! refusal's own message advises the user to do.
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
   end subroutine hf_run

   subroutine hf_calc_gradient(this, fragment, result)
      !! Energy and nuclear gradient of a fragment
      class(hf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call hf_run(this, fragment, result, want_gradient=.true.)
   end subroutine hf_calc_gradient

   subroutine hf_calc_hessian(this, fragment, result)
      !! The analytic Hessian where there is one, central differences otherwise
      !!
      !! **The backend decides, not this routine.** The request goes down and
      !! the answer comes back as `has_hessian`: true means it was computed,
      !! false with no error means it was declined and finite differences take
      !! over. Unlike a gradient, a declined Hessian is not a failure.
      class(hf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call hf_run(this, fragment, result, want_gradient=.true., want_hessian=.true.)
      if (result%has_error) return
      if (result%has_hessian) return

      ! Declined. The finite-difference path runs its own reference point, so
      ! nothing computed above is reusable.
      call result%destroy()
      call finite_difference_hessian(this, fragment, result, verbose=this%options%verbose, &
                                     displacement_in=this%options%hessian_displacement)
   end subroutine hf_calc_hessian

end module mqc_method_hf
