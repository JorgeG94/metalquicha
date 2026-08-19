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
   use mqc_method_config, only: pcm_config_t, properties_config_t
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

   public :: hf_method_t, hf_options_t

   type :: hf_options_t
      !! Hartree-Fock calculation options
      type(properties_config_t) :: properties
      type(pcm_config_t) :: pcm
         !! Continuum solvation. Carried here as well as on the Kohn-Sham
         !! options because a continuum is a property of the reference, not of
         !! the functional: Hartree-Fock, MP2 and coupled cluster all reach the
         !! backend through this type, and leaving it off meant `keywords.pcm`
         !! was read, validated, and then silently dropped for every one of them.
      character(len=32) :: basis_set = "sto-3g"
         !! Orbital basis set name
      character(len=32) :: aux_basis_set = "def2-universal-jkfit"
         !! Auxiliary (JKFIT) basis for the density-fitted J and K
      logical :: aux_basis_named = .false.
         !! Whether the deck asked for `aux_basis_set`. See `scf_config_t`.
      logical :: density_fitting = .false.
      logical :: run_mp2 = .false.
         !! Follow the reference with an MP2 correction. Set by the factory
         !! from the method name rather than by a keyword: "mp2" is a method,
         !! not an option on Hartree-Fock.
      logical :: freeze_core = .false.
      integer :: n_frozen_core = -1
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
         !! Fit J and K rather than computing exact integrals (CPU backend)
      logical :: spherical = .true.
         !! Use spherical (true) or Cartesian (false) basis
      logical :: verbose = .false.
         !! Print SCF iterations
      integer :: device_rank = 0
         !! Node-local MPI rank, for spreading ranks across a node's GPUs
      logical :: unrestricted = .false.
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

      ! SCF settings (from shared scf_config_t)
      logical :: allow_crap_scf = .false.  !! Keep a non-converged SCF instead of failing
      integer :: max_iter = 100
         !! Maximum SCF iterations
      real(dp) :: conv_tol = 1.0e-8_dp
         !! Energy convergence threshold
      real(dp) :: density_tol = 1.0e-6_dp
         !! Density matrix convergence threshold
      logical :: use_diis = .true.
         !! Use DIIS acceleration
      integer :: diis_size = 8
         !! Number of Fock matrices for DIIS
      character(len=16) :: backend = "auto"
         !! Integral backend request: "auto", "cuest"/"gpu", "libcint"/"cpu".
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

   subroutine hf_run(this, fragment, result, want_gradient)
      !! Run a density-fitted RHF calculation through cuEST
      class(hf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in) :: want_gradient

      type(cuest_scf_settings_t) :: settings
      type(error_t) :: backend_error

      settings%pcm = this%options%pcm
      settings%bonding_analysis = this%options%properties%bonding_analysis
      settings%bonding_threshold = this%options%properties%bonding_threshold
      settings%bonding_no_sharing = this%options%properties%bonding_no_sharing
      settings%basis_set = this%options%basis_set
      settings%aux_basis_set = this%options%aux_basis_set
      settings%aux_basis_named = this%options%aux_basis_named
      settings%density_fitting = this%options%density_fitting
      settings%run_mp2 = this%options%run_mp2
      settings%freeze_core = this%options%freeze_core
      settings%n_frozen_core = this%options%n_frozen_core
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
      settings%energy_tol = this%options%conv_tol
      settings%density_tol = this%options%density_tol
      settings%use_diis = this%options%use_diis
      settings%diis_size = this%options%diis_size

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
         call run_libcint_hf(settings, fragment, result, want_gradient)
      case default
#ifdef MQC_WITH_CUEST
         call run_cuest_scf(settings, fragment, result, want_gradient)
#else
         call run_libcint_hf(settings, fragment, result, want_gradient)
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
      !! Hessian by central differences of the analytic gradients
      class(hf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call finite_difference_hessian(this, fragment, result, verbose=this%options%verbose)
   end subroutine hf_calc_hessian

end module mqc_method_hf
