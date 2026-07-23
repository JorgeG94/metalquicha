!! Hartree-Fock method implementation for metalquicha
module mqc_method_hf
   !! Closed-shell Hartree-Fock.
   !!
   !! When built with MQC_WITH_CUEST the integrals come from NVIDIA cuEST on
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
#ifdef MQC_WITH_CUEST
   use mqc_cuest_driver, only: cuest_scf_settings_t, run_cuest_scf
#endif
   implicit none
   private

   public :: hf_method_t, hf_options_t

   type :: hf_options_t
      !! Hartree-Fock calculation options
      character(len=32) :: basis_set = 'sto-3g'
         !! Orbital basis set name
      character(len=32) :: aux_basis_set = 'def2-universal-jkfit'
         !! Auxiliary (JKFIT) basis for the density-fitted J and K
      logical :: spherical = .true.
         !! Use spherical (true) or Cartesian (false) basis
      logical :: verbose = .false.
         !! Print SCF iterations

      ! SCF settings (from shared scf_config_t)
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
      procedure :: calc_hessian => null_hessian  !! Placeholder for Hessian calculation
   end type hf_method_t

contains

   subroutine hf_calc_energy(this, fragment, result)
      !! Calculate electronic energy using Hartree-Fock method
      class(hf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

#ifdef MQC_WITH_CUEST
      call hf_energy_cuest(this, fragment, result, want_gradient=.false.)
#else
      call result%error%set(ERROR_VALIDATION, &
                            "Hartree-Fock requires an integral backend; rebuild with "// &
                            "MQC_ENABLE_CUEST=ON")
      result%has_error = .true.
#endif
   end subroutine hf_calc_energy

#ifdef MQC_WITH_CUEST
   subroutine hf_energy_cuest(this, fragment, result, want_gradient)
      !! Run a density-fitted RHF calculation through cuEST
      class(hf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in) :: want_gradient

      type(cuest_scf_settings_t) :: settings

      settings%basis_set = this%options%basis_set
      settings%aux_basis_set = this%options%aux_basis_set
      settings%functional = ''        ! empty selects pure Hartree-Fock
      settings%spherical = this%options%spherical
      settings%verbose = this%options%verbose
      settings%max_iter = this%options%max_iter
      settings%energy_tol = this%options%conv_tol
      settings%density_tol = this%options%density_tol
      settings%use_diis = this%options%use_diis
      settings%diis_size = this%options%diis_size

      call run_cuest_scf(settings, fragment, result, want_gradient)
   end subroutine hf_energy_cuest
#endif

   subroutine hf_calc_gradient(this, fragment, result)
      !! Calculate energy gradient using Hartree-Fock method
      class(hf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

#ifdef MQC_WITH_CUEST
      call hf_energy_cuest(this, fragment, result, want_gradient=.true.)
#else
      call result%error%set(ERROR_VALIDATION, &
                            "Hartree-Fock requires an integral backend; rebuild with "// &
                            "MQC_ENABLE_CUEST=ON")
      result%has_error = .true.
#endif
   end subroutine hf_calc_gradient

   subroutine null_hessian(this, fragment, result)
      !! Placeholder for Hessian calculation
      class(hf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call result%error%set(ERROR_VALIDATION, &
                            "Hartree-Fock Hessians are not implemented yet")
      result%has_error = .true.
   end subroutine null_hessian

end module mqc_method_hf
