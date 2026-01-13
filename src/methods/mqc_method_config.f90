!! Unified configuration type for quantum chemistry method creation
module mqc_method_config
   !! Provides a union-style configuration type that holds all possible settings
   !! for any quantum chemistry method. Each method's factory configuration routine
   !! extracts only the relevant fields - unused fields are ignored.
   use pic_types, only: int32, dp
   use mqc_method_types, only: METHOD_TYPE_UNKNOWN
   implicit none
   private

   public :: method_config_t

   type :: method_config_t
      !! Union-style configuration for all quantum chemistry methods
      !!
      !! This type contains all possible configuration options for any QC method.
      !! When creating a specific method, only the relevant fields are used.
      !! This design avoids inheritance complexity while keeping a clean API.

      !----- Core settings (all methods) -----
      integer(int32) :: method_type = METHOD_TYPE_UNKNOWN
         !! Method type constant (METHOD_TYPE_GFN1, METHOD_TYPE_GFN2, METHOD_TYPE_HF)
      logical :: verbose = .false.
         !! Enable verbose output during calculations

      !----- XTB-specific (GFN1, GFN2) -----
      real(dp) :: accuracy = 0.01_dp
         !! Numerical accuracy parameter for xTB
      real(dp) :: electronic_temp = 300.0_dp
         !! Electronic temperature in Kelvin (for Fermi smearing)

      ! Solvation settings
      character(len=32) :: solvent = ''
         !! Solvent name: "water", "ethanol", etc. Empty for gas phase
      character(len=16) :: solvation_model = ''
         !! Solvation model: "alpb" (default), "gbsa", or "cpcm"
      logical :: use_cds = .true.
         !! Include non-polar CDS (Cavity-Dispersion-Solvent) terms
      logical :: use_shift = .true.
         !! Include solution state shift correction
      real(dp) :: dielectric = -1.0_dp
         !! Direct dielectric constant (-1 = use solvent lookup table)
      integer :: cpcm_nang = 110
         !! Number of angular grid points for CPCM cavity
      real(dp) :: cpcm_rscale = 1.0_dp
         !! Radii scaling factor for CPCM cavity

      !----- HF/ab-initio specific -----
      character(len=32) :: basis_set = 'sto-3g'
         !! Basis set name for ab-initio methods
      integer :: max_scf_iter = 100
         !! Maximum SCF iterations
      real(dp) :: scf_convergence = 1.0e-8_dp
         !! SCF energy convergence threshold (Hartree)
      logical :: use_spherical = .true.
         !! Use spherical (true) or Cartesian (false) basis functions

   contains
      procedure :: has_solvation => config_has_solvation
      procedure :: reset => config_reset
   end type method_config_t

contains

   pure logical function config_has_solvation(this)
      !! Check if solvation is configured
      !!
      !! Returns true if either a solvent name is specified or
      !! a positive dielectric constant is set.
      class(method_config_t), intent(in) :: this

      config_has_solvation = len_trim(this%solvent) > 0 .or. this%dielectric > 0.0_dp
   end function config_has_solvation

   subroutine config_reset(this)
      !! Reset all configuration values to defaults
      class(method_config_t), intent(inout) :: this

      this%method_type = METHOD_TYPE_UNKNOWN
      this%verbose = .false.

      ! XTB defaults
      this%accuracy = 0.01_dp
      this%electronic_temp = 300.0_dp
      this%solvent = ''
      this%solvation_model = ''
      this%use_cds = .true.
      this%use_shift = .true.
      this%dielectric = -1.0_dp
      this%cpcm_nang = 110
      this%cpcm_rscale = 1.0_dp

      ! HF defaults
      this%basis_set = 'sto-3g'
      this%max_scf_iter = 100
      this%scf_convergence = 1.0e-8_dp
      this%use_spherical = .true.
   end subroutine config_reset

end module mqc_method_config
