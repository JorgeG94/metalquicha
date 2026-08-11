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
   use mqc_method_base, only: qc_method_t
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_semi_numerical_hessian, only: finite_difference_hessian
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_cuest_bridge, only: run_cuest_scf
   use mqc_libcint_bridge, only: run_libcint_hf
   implicit none
   private

   public :: hf_method_t, hf_options_t

   type :: hf_options_t
      !! Hartree-Fock calculation options
      character(len=32) :: basis_set = "sto-3g"
         !! Orbital basis set name
      character(len=32) :: aux_basis_set = "def2-universal-jkfit"
         !! Auxiliary (JKFIT) basis for the density-fitted J and K
      logical :: density_fitting = .false.
      logical :: run_mp2 = .false.
         !! Follow the reference with an MP2 correction. Set by the factory
         !! from the method name rather than by a keyword: "mp2" is a method,
         !! not an option on Hartree-Fock.
      logical :: freeze_core = .false.
      integer :: n_frozen_core = -1
      real(dp) :: scs_ss = 1.0_dp
      real(dp) :: scs_os = 1.0_dp
         !! Fit J and K rather than computing exact integrals (CPU backend)
      logical :: spherical = .true.
         !! Use spherical (true) or Cartesian (false) basis
      logical :: verbose = .false.
         !! Print SCF iterations
      integer :: device_rank = 0
         !! Node-local MPI rank, for spreading ranks across a node's GPUs
      logical :: unrestricted = .false.
         !! Force UHF/UKS even for a closed shell
      character(len=16) :: guess = "gwh"
         !! Initial guess: 'core', 'gwh' or 'sac'

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

      settings%basis_set = this%options%basis_set
      settings%aux_basis_set = this%options%aux_basis_set
      settings%density_fitting = this%options%density_fitting
      settings%run_mp2 = this%options%run_mp2
      settings%freeze_core = this%options%freeze_core
      settings%n_frozen_core = this%options%n_frozen_core
      settings%scs_ss = this%options%scs_ss
      settings%scs_os = this%options%scs_os
      settings%functional = ""        ! empty selects pure Hartree-Fock
      settings%spherical = this%options%spherical
      settings%verbose = this%options%verbose
      settings%device_rank = this%options%device_rank
      settings%unrestricted = this%options%unrestricted
      settings%guess = this%options%guess
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
#ifdef MQC_WITH_CUEST
      call run_cuest_scf(settings, fragment, result, want_gradient)
#else
      ! Without cuEST this is the CPU backend's call, real or stubbed. Routing
      ! a build with neither to the cuEST stub would name only one of the two
      ! ways to fix it; the libcint stub names both.
      call run_libcint_hf(settings, fragment, result, want_gradient)
#endif
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
