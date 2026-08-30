!! Hartree-Fock method implementation for metalquicha
module mqc_method_hf
   !! Closed-shell Hartree-Fock.
   !!
   !! With the cuEST backend (CMake, MQC_ENABLE_CUEST=ON) the integrals come
   !! the GPU and the SCF runs in `mqc_cuest_scf`. Because cuEST offers no
   !! conventional four-index ERI path, J and K are always density-fitted, so
   !! an auxiliary (JKFIT) basis is required alongside the orbital basis.
   !!
   !! Without that backend the method reports a clear error rather than a
   !! placeholder number -- a silently wrong energy inside a many-body
   !! expansion is far worse than a failed fragment.
   use pic_types, only: dp
   use mqc_config_types, only: guess_step_t
   use mqc_method_base, only: qc_method_t
   use mqc_method_config, only: scf_options_t, pcm_config_t, properties_config_t
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_semi_numerical_hessian, only: finite_difference_hessian
   use mqc_cuest_iface, only: apply_scf_settings, cuest_scf_settings_t, parse_backend_name, &
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
         !! from the method name rather than by a keyword: "mp2" is a method,
         !! not an option on Hartree-Fock.
      logical :: corr_density_fitting = .false.
      real(dp) :: scs_ss = 1.0_dp
      real(dp) :: scs_os = 1.0_dp
      logical :: run_cc = .false.
         !! Follow the reference with coupled cluster. Set by the factory from
         !! the method name, for the same reason `run_mp2` is: "ccsd" is a
         !! method, not an option on Hartree-Fock.
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
      !! Calculate electronic energy using Hartree-Fock method
      class(hf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call hf_run(this, fragment, result, want_gradient=.false.)
   end subroutine hf_calc_energy

   subroutine hf_run(this, fragment, result, want_gradient, want_hessian)
      !! Run a density-fitted RHF calculation through cuEST
      class(hf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in) :: want_gradient
      logical, intent(in), optional :: want_hessian
         !! Only the CPU backend has an analytic Hessian, and only for some of
         !! what reaches it. Passed along rather than decided here: which cases
         !! qualify is a fact about that backend.

      type(cuest_scf_settings_t) :: settings
      type(error_t) :: backend_error

      call apply_scf_settings(settings, this%options)
      settings%bonding_analysis = this%options%properties%bonding_analysis
      settings%bonding_threshold = this%options%properties%bonding_threshold
      settings%bonding_energy = this%options%properties%bonding_energy
      if (allocated(this%options%properties%fukui_population)) then
         settings%fukui_population = this%options%properties%fukui_population
      end if
      if (allocated(this%options%properties%charges_scheme)) then
         settings%charges_scheme = this%options%properties%charges_scheme
      end if
      settings%bonding_no_sharing = this%options%properties%bonding_no_sharing
      settings%bonding_restrict_localization = &
         this%options%properties%bonding_restrict_localization
      settings%bonding_no_sharing_ci = this%options%properties%bonding_no_sharing_ci
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
      ! Resolved here rather than carried as a string, so an unknown name fails
      ! once, before any integrals, instead of at each dispatch.
      call parse_backend_name(this%options%backend, settings%backend, backend_error)
      if (backend_error%has_error()) then
         call result%error%set(ERROR_VALIDATION, backend_error%get_message())
         result%has_error = .true.
         return
      end if

      ! cuEST when this build has it, because that is the production path and
      ! the only one with gradients. The CPU backend takes over when it does
      ! not, which is how a fragmented Hartree-Fock runs on a laptop at all --
      ! and, when both are built, gives the GPU one something to be compared
      ! against.
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
   end subroutine hf_run

   subroutine hf_calc_gradient(this, fragment, result)
      !! Calculate energy gradient using Hartree-Fock method
      class(hf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call hf_run(this, fragment, result, want_gradient=.true.)
   end subroutine hf_calc_gradient

   subroutine hf_calc_hessian(this, fragment, result)
      !! The analytic Hessian where there is one, central differences otherwise
      !!
      !! **The backend decides, not this routine.** Whether an analytic Hessian
      !! applies depends on things the backend knows and this does not -- what
      !! the SCF actually fitted, whether it converged restricted, whether a
      !! functional was built. So the request goes down and the answer comes
      !! back as `has_hessian`: true means it was computed, false with no error
      !! means it was declined, and finite differences take over.
      !!
      !! That is deliberately not the pattern the rest of this file uses. A
      !! backend that cannot honour a *gradient* request refuses loudly,
      !! because there is no other way to get one and a silent substitution
      !! would report a provenance that was not true. A Hessian has a second
      !! way to get one that is correct and merely slower, so declining is not
      !! a failure and should not read as one.
      class(hf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call hf_run(this, fragment, result, want_gradient=.true., want_hessian=.true.)
      if (result%has_error) return
      if (result%has_hessian) return

      ! Declined. Nothing computed above is reusable -- the finite-difference
      ! path runs its own reference point -- so this starts over rather than
      ! trying to keep the SCF that just ran.
      call result%destroy()
      call finite_difference_hessian(this, fragment, result, verbose=this%options%verbose)
   end subroutine hf_calc_hessian

end module mqc_method_hf
