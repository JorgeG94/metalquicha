!! Calculation keyword types for structured configuration
module mqc_calculation_keywords
   !! Provides structured keyword types for calculation-specific settings
   !! These types are embedded in driver_config_t to organize keywords by category
   use pic_types, only: dp
   use mqc_calculation_defaults, only: DEFAULT_DISPLACEMENT, DEFAULT_TEMPERATURE, &
                                       DEFAULT_PRESSURE, DEFAULT_SCF_MAXITER, &
                                       DEFAULT_SCF_CONV, DEFAULT_USE_DIIS, &
                                       DEFAULT_AIMD_DT, DEFAULT_AIMD_NSTEPS, &
                                       DEFAULT_AIMD_TEMPERATURE, DEFAULT_AIMD_OUTPUT_FREQ
   implicit none
   private

   public :: hessian_keywords_t, aimd_keywords_t, scf_keywords_t

   type :: hessian_keywords_t
      !! Hessian calculation keywords
      real(dp) :: displacement = DEFAULT_DISPLACEMENT  !! Finite difference displacement (Bohr)
      real(dp) :: temperature = DEFAULT_TEMPERATURE    !! Temperature for thermochemistry (K)
      real(dp) :: pressure = DEFAULT_PRESSURE          !! Pressure for thermochemistry (atm)
   end type hessian_keywords_t

   type :: aimd_keywords_t
      !! Ab initio molecular dynamics keywords
      real(dp) :: dt = DEFAULT_AIMD_DT                       !! Timestep (femtoseconds)
      integer :: nsteps = DEFAULT_AIMD_NSTEPS                !! Number of MD steps (0 = no AIMD)
      real(dp) :: initial_temperature = DEFAULT_AIMD_TEMPERATURE  !! Initial temperature for velocity init (K)
      integer :: output_frequency = DEFAULT_AIMD_OUTPUT_FREQ !! Write output every N steps
   end type aimd_keywords_t

   type :: scf_keywords_t
      !! SCF calculation keywords (placeholder for future use)
      logical :: use_diis = DEFAULT_USE_DIIS                   !! Use DIIS acceleration
      integer :: max_iterations = DEFAULT_SCF_MAXITER          !! Maximum SCF iterations
      real(dp) :: convergence_threshold = DEFAULT_SCF_CONV     !! Convergence threshold for SCF
   end type scf_keywords_t

end module mqc_calculation_keywords
