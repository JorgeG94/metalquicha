!! Unified configuration type for quantum chemistry method creation
module mqc_method_config
   !! Provides configuration types for all quantum chemistry methods.
   !! Uses composition pattern: method_config_t contains nested config types
   !! for each method family. The factory reads from the appropriate nested type.
   use pic_types, only: int32, dp
   use mqc_method_types, only: METHOD_TYPE_UNKNOWN
   implicit none
   private

   public :: method_config_t
   public :: xtb_config_t, hf_config_t, dft_config_t, mcscf_config_t

   !============================================================================
   ! XTB Configuration (GFN1, GFN2)
   !============================================================================
   type :: xtb_config_t
      !! Configuration for semi-empirical xTB methods
      real(dp) :: accuracy = 0.01_dp
         !! Numerical accuracy parameter
      real(dp) :: electronic_temp = 300.0_dp
         !! Electronic temperature in Kelvin (Fermi smearing)

      ! Solvation
      character(len=32) :: solvent = ''
         !! Solvent name: "water", "ethanol", etc. Empty for gas phase
      character(len=16) :: solvation_model = ''
         !! Solvation model: "alpb", "gbsa", "cpcm"
      logical :: use_cds = .true.
         !! Include non-polar CDS terms
      logical :: use_shift = .true.
         !! Include solution state shift
      real(dp) :: dielectric = -1.0_dp
         !! Dielectric constant (-1 = use solvent table)
      integer :: cpcm_nang = 110
         !! Angular grid points for CPCM
      real(dp) :: cpcm_rscale = 1.0_dp
         !! Radii scaling for CPCM
   contains
      procedure :: has_solvation => xtb_has_solvation
   end type xtb_config_t

   !============================================================================
   ! Hartree-Fock Configuration
   !============================================================================
   type :: hf_config_t
      !! Configuration for Hartree-Fock method
      integer :: max_scf_iter = 100
         !! Maximum SCF iterations
      real(dp) :: energy_convergence = 1.0e-8_dp
         !! Energy convergence threshold (Hartree)
      real(dp) :: density_convergence = 1.0e-6_dp
         !! Density matrix convergence threshold
      logical :: use_diis = .true.
         !! Use DIIS acceleration
      integer :: diis_size = 8
         !! Number of Fock matrices for DIIS
   end type hf_config_t

   !============================================================================
   ! DFT Configuration
   !============================================================================
   type :: dft_config_t
      !! Configuration for Kohn-Sham DFT method
      character(len=32) :: functional = 'b3lyp'
         !! XC functional: "lda", "pbe", "b3lyp", "m06-2x", etc.

      ! Integration grid
      character(len=16) :: grid_type = 'medium'
         !! Grid quality: "coarse", "medium", "fine", "ultrafine"
      integer :: radial_points = 75
         !! Radial grid points per atom
      integer :: angular_points = 302
         !! Angular grid points (Lebedev)

      ! SCF settings (inherits some from HF)
      integer :: max_scf_iter = 100
         !! Maximum SCF iterations
      real(dp) :: energy_convergence = 1.0e-8_dp
         !! Energy convergence threshold
      real(dp) :: density_convergence = 1.0e-6_dp
         !! Density convergence threshold
      logical :: use_diis = .true.
         !! Use DIIS acceleration
      integer :: diis_size = 8
         !! DIIS history size

      ! Density fitting
      logical :: use_density_fitting = .false.
         !! Use RI-J approximation
      character(len=32) :: aux_basis_set = ''
         !! Auxiliary basis for density fitting

      ! Dispersion correction
      logical :: use_dispersion = .false.
         !! Add empirical dispersion
      character(len=8) :: dispersion_type = 'd3bj'
         !! Dispersion type: "d3", "d3bj", "d4"
   end type dft_config_t

   !============================================================================
   ! MCSCF/CASSCF Configuration
   !============================================================================
   type :: mcscf_config_t
      !! Configuration for MCSCF/CASSCF method

      ! Active space definition
      integer :: n_active_electrons = 0
         !! Number of active electrons
      integer :: n_active_orbitals = 0
         !! Number of active orbitals
      integer :: n_inactive_orbitals = -1
         !! Inactive orbitals (-1 = auto from nelec)

      ! State averaging
      integer :: n_states = 1
         !! Number of states for SA-CASSCF
      real(dp), allocatable :: state_weights(:)
         !! State weights (must sum to 1)

      ! Convergence
      integer :: max_macro_iter = 100
         !! Maximum orbital optimization iterations
      integer :: max_micro_iter = 50
         !! Maximum CI iterations per macro
      real(dp) :: orbital_convergence = 1.0e-6_dp
         !! Orbital gradient threshold
      real(dp) :: ci_convergence = 1.0e-8_dp
         !! CI energy threshold

      ! Perturbative corrections
      logical :: use_pt2 = .false.
         !! Apply CASPT2/NEVPT2 after CASSCF
      character(len=16) :: pt2_type = 'nevpt2'
         !! PT2 flavor: "caspt2", "nevpt2"
      real(dp) :: ipea_shift = 0.25_dp
         !! IPEA shift for CASPT2
      real(dp) :: imaginary_shift = 0.0_dp
         !! Imaginary shift for intruder states
   end type mcscf_config_t

   !============================================================================
   ! Main Configuration Type (Composition)
   !============================================================================
   type :: method_config_t
      !! Master configuration containing all method-specific configs
      !!
      !! Usage:
      !!   config%method_type = METHOD_TYPE_DFT
      !!   config%basis_set = 'cc-pvdz'
      !!   config%dft%functional = 'pbe0'
      !!   config%dft%use_dispersion = .true.

      !----- Common settings (all ab-initio methods) -----
      integer(int32) :: method_type = METHOD_TYPE_UNKNOWN
         !! Method type constant
      logical :: verbose = .false.
         !! Enable verbose output
      character(len=32) :: basis_set = 'sto-3g'
         !! Basis set name (HF, DFT, MCSCF)
      logical :: use_spherical = .true.
         !! Spherical vs Cartesian basis functions

      !----- Method-specific configurations -----
      type(xtb_config_t) :: xtb
         !! XTB settings (GFN1, GFN2)
      type(hf_config_t) :: hf
         !! Hartree-Fock settings
      type(dft_config_t) :: dft
         !! DFT settings
      type(mcscf_config_t) :: mcscf
         !! MCSCF/CASSCF settings

   contains
      procedure :: reset => config_reset
   end type method_config_t

contains

   pure logical function xtb_has_solvation(this)
      !! Check if solvation is configured for XTB
      class(xtb_config_t), intent(in) :: this
      xtb_has_solvation = len_trim(this%solvent) > 0 .or. this%dielectric > 0.0_dp
   end function xtb_has_solvation

   subroutine config_reset(this)
      !! Reset all configuration values to defaults
      class(method_config_t), intent(inout) :: this

      ! Common settings
      this%method_type = METHOD_TYPE_UNKNOWN
      this%verbose = .false.
      this%basis_set = 'sto-3g'
      this%use_spherical = .true.

      ! XTB defaults
      this%xtb%accuracy = 0.01_dp
      this%xtb%electronic_temp = 300.0_dp
      this%xtb%solvent = ''
      this%xtb%solvation_model = ''
      this%xtb%use_cds = .true.
      this%xtb%use_shift = .true.
      this%xtb%dielectric = -1.0_dp
      this%xtb%cpcm_nang = 110
      this%xtb%cpcm_rscale = 1.0_dp

      ! HF defaults
      this%hf%max_scf_iter = 100
      this%hf%energy_convergence = 1.0e-8_dp
      this%hf%density_convergence = 1.0e-6_dp
      this%hf%use_diis = .true.
      this%hf%diis_size = 8

      ! DFT defaults
      this%dft%functional = 'b3lyp'
      this%dft%grid_type = 'medium'
      this%dft%radial_points = 75
      this%dft%angular_points = 302
      this%dft%max_scf_iter = 100
      this%dft%energy_convergence = 1.0e-8_dp
      this%dft%density_convergence = 1.0e-6_dp
      this%dft%use_diis = .true.
      this%dft%diis_size = 8
      this%dft%use_density_fitting = .false.
      this%dft%aux_basis_set = ''
      this%dft%use_dispersion = .false.
      this%dft%dispersion_type = 'd3bj'

      ! MCSCF defaults
      this%mcscf%n_active_electrons = 0
      this%mcscf%n_active_orbitals = 0
      this%mcscf%n_inactive_orbitals = -1
      this%mcscf%n_states = 1
      if (allocated(this%mcscf%state_weights)) deallocate (this%mcscf%state_weights)
      this%mcscf%max_macro_iter = 100
      this%mcscf%max_micro_iter = 50
      this%mcscf%orbital_convergence = 1.0e-6_dp
      this%mcscf%ci_convergence = 1.0e-8_dp
      this%mcscf%use_pt2 = .false.
      this%mcscf%pt2_type = 'nevpt2'
      this%mcscf%ipea_shift = 0.25_dp
      this%mcscf%imaginary_shift = 0.0_dp
   end subroutine config_reset

end module mqc_method_config
