!! Density Functional Theory (DFT) method implementation for metalquicha
module mqc_method_dft
   !! Closed-shell Kohn-Sham DFT.
   !!
   !! With MQC_WITH_CUEST the Coulomb and exact-exchange matrices come from
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
   use mqc_method_base, only: qc_method_t
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_error, only: error_t, ERROR_VALIDATION
#ifdef MQC_WITH_CUEST
   use mqc_cuest_driver, only: cuest_scf_settings_t, run_cuest_scf
#endif
   implicit none
   private

   public :: dft_method_t, dft_options_t

   type :: dft_options_t
      !! DFT calculation options
      character(len=32) :: basis_set = 'sto-3g'
         !! Basis set name
      character(len=32) :: functional = 'b3lyp'
         !! Exchange-correlation functional
      integer :: max_iter = 100
         !! Maximum SCF iterations
      real(dp) :: energy_tol = 1.0e-8_dp
         !! Energy convergence threshold
      real(dp) :: density_tol = 1.0e-6_dp
         !! Density matrix convergence threshold
      logical :: spherical = .true.
         !! Use spherical (true) or Cartesian (false) basis
      logical :: verbose = .false.
         !! Print SCF iterations

      ! Grid settings
      character(len=16) :: grid_type = 'medium'
         !! Integration grid quality
      integer :: radial_points = 75
         !! Number of radial grid points per atom
      integer :: angular_points = 302
         !! Number of angular grid points (Lebedev)

      ! Density fitting
      logical :: use_density_fitting = .false.
         !! Use RI-J approximation
      character(len=32) :: aux_basis_set = 'def2-universal-jkfit'
         !! Auxiliary (JKFIT) basis. Required by the cuEST backend.

      ! Dispersion correction
      logical :: use_dispersion = .false.
         !! Add empirical dispersion correction
      character(len=8) :: dispersion_type = 'd3bj'
         !! Dispersion type: "d3", "d3bj", "d4"

      ! DIIS acceleration
      logical :: use_diis = .true.
         !! Use DIIS for SCF convergence
      integer :: diis_size = 8
         !! Number of Fock matrices in DIIS
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

#ifdef MQC_WITH_CUEST
      type(cuest_scf_settings_t) :: settings

      if (this%options%use_dispersion) then
         ! Silently dropping a dispersion correction would bias every fragment
         ! energy the same way, which is exactly the kind of error a many-body
         ! expansion amplifies rather than cancels.
         call result%error%set(ERROR_VALIDATION, &
                               "Empirical dispersion ("//trim(this%options%dispersion_type)// &
                               ") is not implemented for the cuEST backend")
         result%has_error = .true.
         return
      end if

      settings%basis_set = this%options%basis_set
      settings%aux_basis_set = this%options%aux_basis_set
      settings%functional = this%options%functional
      settings%spherical = this%options%spherical
      settings%verbose = this%options%verbose
      settings%max_iter = this%options%max_iter
      settings%energy_tol = this%options%energy_tol
      settings%density_tol = this%options%density_tol
      settings%use_diis = this%options%use_diis
      settings%diis_size = this%options%diis_size
      settings%radial_points = this%options%radial_points
      settings%angular_points = this%options%angular_points

      call run_cuest_scf(settings, fragment, result)
#else
      call result%error%set(ERROR_VALIDATION, &
                            "DFT requires an integral backend; rebuild with MQC_ENABLE_CUEST=ON")
      result%has_error = .true.
#endif
   end subroutine dft_calc_energy

   subroutine dft_calc_gradient(this, fragment, result)
      !! Calculate energy gradient using Kohn-Sham DFT
      class(dft_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      ! Analytic KS gradients need the cuEST *DerivativeCompute entry points,
      ! an energy-weighted density and the XC grid response; not wired up yet.
      call this%calc_energy(fragment, result)
      if (result%has_error) return

      call result%error%set(ERROR_VALIDATION, "DFT gradients are not implemented yet")
      result%has_error = .true.
   end subroutine dft_calc_gradient

   subroutine dft_calc_hessian(this, fragment, result)
      !! Placeholder for Hessian calculation
      class(dft_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call result%error%set(ERROR_VALIDATION, "DFT Hessians are not implemented yet")
      result%has_error = .true.
   end subroutine dft_calc_hessian

end module mqc_method_dft
