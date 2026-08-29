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
   use mqc_method_config, only: pcm_config_t, properties_config_t
   use mqc_method_base, only: qc_method_t
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_semi_numerical_hessian, only: finite_difference_hessian
   use mqc_cuest_iface, only: cuest_scf_settings_t, parse_backend_name, &
                              BACKEND_CUEST, BACKEND_LIBCINT
   use mqc_cuest_bridge, only: run_cuest_scf
   use mqc_libcint_bridge, only: run_libcint_hf
   implicit none
   private

   public :: dft_method_t, dft_options_t

   type :: dft_options_t
      !! DFT calculation options
      character(len=32) :: basis_set = "sto-3g"
      character(len=32) :: ecp_set = ""
         !! Effective core potential set, empty for an all-electron run
         !! Basis set name
      character(len=32) :: functional = "b3lyp"
         !! Exchange-correlation functional
      integer :: max_iter = 100
         !! Maximum SCF iterations
      real(dp) :: energy_tol = 1.0e-8_dp
         !! Energy convergence threshold
      real(dp) :: density_tol = 1.0e-6_dp
      real(dp) :: linear_dependence = 0.0_dp
         !! Zero means the orthogonaliser's own cutoff. See `scf_config_t`.
         !! Density matrix convergence threshold
      real(dp) :: level_shift = 0.0_dp
         !! Hartree added to the virtual block before each diagonalisation.
         !! Zero is off. See `scf_config_t`.
      logical :: spherical = .true.
         !! Use spherical (true) or Cartesian (false) basis
      logical :: verbose = .false.
         !! Print SCF iterations
      integer :: device_rank = 0
         !! Node-local MPI rank, for spreading ranks across a node's GPUs
      logical :: unrestricted = .false.
      logical :: allow_crap_scf = .false.
         !! Keep a non-converged SCF instead of failing.
         !!
         !! Here because it was *not* here: `keywords.scf.allow_crap_scf` was
         !! read from the deck, passed the schema, reached the config, and was
         !! then dropped for every DFT run -- this type had no such field,
         !! `configure_dft` did not set it, and the block below did not copy
         !! it. Hartree-Fock had all three. The backend and the workflow were
         !! both correct and were simply never told, so the deck was obeyed for
         !! HF and silently ignored for Kohn-Sham.
         !!
         !! **This type and `hf_options_t` share twenty-two of their
         !! thirty-five fields**, hand-copied through three layers apiece, and
         !! have already drifted: `conv_tol` here is `energy_tol` there, and
         !! `density_fitting` is `use_density_fitting`. An SCF keyword added to
         !! one and not the other fails exactly this way, with no error
         !! anywhere. A shared `scf_options_t` would remove the class.
         !! Force UHF/UKS even for a closed shell
      character(len=32) :: guess = "auto"
         !! Initial guess: 'core', 'gwh', 'sac', 'sad', 'basis_set_projection',
         !! or 'auto'
      type(guess_step_t), allocatable :: guess_steps(:)
         !! The basis ladder for 'basis_set_projection', one entry per
         !! preliminary SCF in order.
         !!
         !! 'auto' means the backend picks, because the best starting point
         !! is a property of the backend rather than of the request: the CPU
         !! path resolves it to 'sad', and cuEST to 'gwh', each having
         !! measured its own. An explicit spelling always wins over both.

      ! Grid settings
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
      logical :: use_density_fitting = .false.
         !! Use RI-J approximation
      logical :: cartesian = .false.
         !! `model.cartesian`; see `mqc_config_t`.
      character(len=32) :: aux_basis_set = "def2-universal-jkfit"
         !! Auxiliary (JKFIT) basis. Required by the cuEST backend.
      logical :: aux_basis_named = .false.
         !! Whether the deck asked for it. See `scf_config_t`.

      ! Correlation, for the double hybrids. Nothing else on this path has a
      ! correlated term, but `b2plyp` and its relatives carry an MP2 and read
      ! `keywords.correlation` like any other method that does. Absent these,
      ! a deck asking to freeze the core got an all-electron answer with
      ! nothing in the output to say the request had been dropped.
      logical :: freeze_core = .false.
      integer :: n_frozen_core = -1
         !! -1 counts the core from the elements

      ! Dispersion correction
      logical :: use_dispersion = .false.
         !! Add empirical dispersion correction
      character(len=8) :: dispersion_type = "d3bj"
         !! Dispersion type: "d3", "d3bj", "d4"

      ! DIIS acceleration
      logical :: use_diis = .true.
         !! Use DIIS for SCF convergence
      integer :: diis_size = 8
         !! Number of Fock matrices in DIIS
      type(properties_config_t) :: properties
      type(pcm_config_t) :: pcm
         !! Continuum solvation. Only the cuEST path implements it; the CPU
         !! backend ignores it, which `run_cuest_scf`'s stub makes visible.
      character(len=16) :: backend = "auto"
         !! Integral backend request: "auto", "cuest"/"gpu", "libcint"/"cpu".
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

      settings%basis_set = this%options%basis_set
      settings%cartesian = this%options%cartesian
      settings%ecp_set = this%options%ecp_set
      settings%aux_basis_set = this%options%aux_basis_set
      settings%aux_basis_named = this%options%aux_basis_named
      settings%functional = this%options%functional
      settings%spherical = this%options%spherical
      settings%verbose = this%options%verbose
      settings%device_rank = this%options%device_rank
      ! Resolved here rather than carried as a string, so an unknown name fails
      ! once, before any integrals, instead of at each dispatch.
      call parse_backend_name(this%options%backend, settings%backend, backend_error)
      if (backend_error%has_error()) then
         call result%error%set(ERROR_VALIDATION, backend_error%get_message())
         result%has_error = .true.
         return
      end if
      settings%unrestricted = this%options%unrestricted
      settings%guess = this%options%guess
      if (allocated(this%options%guess_steps)) settings%guess_steps = this%options%guess_steps
      settings%max_iter = this%options%max_iter
      settings%allow_crap_scf = this%options%allow_crap_scf
      settings%energy_tol = this%options%energy_tol
      settings%density_tol = this%options%density_tol
      settings%level_shift = this%options%level_shift
      settings%linear_dependence = this%options%linear_dependence
      settings%use_diis = this%options%use_diis
      settings%diis_size = this%options%diis_size
      settings%radial_points = this%options%radial_points
      settings%angular_points = this%options%angular_points
      settings%grid_level = this%options%grid_level
      settings%nlc_grid_level = this%options%nlc_grid_level
      settings%screening_tolerance = this%options%screening_tolerance
      settings%block_size = this%options%block_size
      settings%pcm = this%options%pcm
      ! Where the molecule reacts. Absent here until now, so a DFT deck asking
      ! for `properties.fukui` got a normal DFT run and no analysis, with no
      ! error to say why: the bridge gates on this being allocated.
      if (allocated(this%options%properties%fukui_population)) then
         settings%fukui_population = this%options%properties%fukui_population
      end if
      if (allocated(this%options%properties%charges_scheme)) then
         settings%charges_scheme = this%options%properties%charges_scheme
      end if
      settings%bonding_analysis = this%options%properties%bonding_analysis
      settings%bonding_threshold = this%options%properties%bonding_threshold
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
      settings%density_fitting = this%options%use_density_fitting
      settings%freeze_core = this%options%freeze_core
      settings%n_frozen_core = this%options%n_frozen_core

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
