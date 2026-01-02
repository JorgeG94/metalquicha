!! Calculation configuration types for keywords (hessian, aimd, etc.)
module mqc_calculation_config
   !! Provides structured configuration for calculation-specific keywords
   !! This separates calculation parameters (hessian, aimd) from driver parameters (fragmentation)
   use pic_types, only: dp
   implicit none
   private

   public :: hessian_keywords_t, aimd_keywords_t, calculation_config_t

   type :: hessian_keywords_t
      !! Hessian calculation keywords
      real(dp) :: displacement = 0.001_dp  !! Finite difference displacement (Bohr)
   end type hessian_keywords_t

   type :: aimd_keywords_t
      !! Ab initio molecular dynamics keywords
      real(dp) :: dt = 1.0_dp                    !! Timestep (femtoseconds)
      integer :: nsteps = 0                      !! Number of MD steps (0 = no AIMD)
      real(dp) :: initial_temperature = 300.0_dp  !! Initial temperature for velocity init (K)
      integer :: output_frequency = 1            !! Write output every N steps
   end type aimd_keywords_t

   type :: calculation_config_t
      !! Container for all calculation-specific configuration
      type(hessian_keywords_t) :: hessian
      type(aimd_keywords_t) :: aimd
   contains
      procedure :: get_hess_keywords => calc_config_get_hess_keywords
      procedure :: get_aimd_keywords => calc_config_get_aimd_keywords
   end type calculation_config_t

contains

   function calc_config_get_hess_keywords(this) result(hess_kw)
      !! Get Hessian keywords from calculation config
      class(calculation_config_t), intent(in) :: this
      type(hessian_keywords_t) :: hess_kw
      hess_kw = this%hessian
   end function calc_config_get_hess_keywords

   function calc_config_get_aimd_keywords(this) result(aimd_kw)
      !! Get AIMD keywords from calculation config
      class(calculation_config_t), intent(in) :: this
      type(aimd_keywords_t) :: aimd_kw
      aimd_kw = this%aimd
   end function calc_config_get_aimd_keywords

end module mqc_calculation_config
