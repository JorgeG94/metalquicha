!! Calculation keyword types for structured configuration
module mqc_calculation_keywords
   !! Provides structured keyword types for calculation-specific settings
   !! These types are embedded in driver_config_t to organize keywords by category
   use pic_types, only: dp
   implicit none
   private

   public :: hessian_keywords_t, aimd_keywords_t, scf_keywords_t

   type :: hessian_keywords_t
      !! Hessian calculation keywords
      real(dp) :: displacement = 0.005_dp  !! Finite difference displacement (Bohr)
      real(dp) :: temperature = 298.15_dp  !! Temperature for thermochemistry (K)
      real(dp) :: pressure = 1.0_dp        !! Pressure for thermochemistry (atm)
   end type hessian_keywords_t

   type :: aimd_keywords_t
      !! Ab initio molecular dynamics keywords
      real(dp) :: dt = 1.0_dp                    !! Timestep (femtoseconds)
      integer :: nsteps = 0                      !! Number of MD steps (0 = no AIMD)
      real(dp) :: initial_temperature = 300.0_dp  !! Initial temperature for velocity init (K)
      integer :: output_frequency = 1            !! Write output every N steps
   end type aimd_keywords_t

   type :: scf_keywords_t
      !! SCF calculation keywords (placeholder for future use)
      logical :: use_diis = .true.               !! Use DIIS acceleration
      integer :: max_iterations = 100             !! Maximum SCF iterations
      real(dp) :: convergence_threshold = 1.0e-6_dp  !! Convergence threshold for SCF
   end type scf_keywords_t

end module mqc_calculation_keywords
