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
   public :: scf_config_t, xtb_config_t, dft_config_t, mcscf_config_t
   public :: correlation_config_t, cc_config_t, f12_config_t

   !============================================================================
   ! SCF Configuration (shared by HF and DFT)
   !============================================================================
   type :: scf_config_t
      !! Shared SCF settings for HF and DFT methods
      integer :: max_iter = 100
         !! Maximum SCF iterations
      real(dp) :: energy_convergence = 1.0e-8_dp
         !! Energy convergence threshold (Hartree)
      real(dp) :: density_convergence = 1.0e-6_dp
         !! Density matrix convergence threshold
      logical :: use_diis = .true.
         !! Use DIIS acceleration
      integer :: diis_size = 8
         !! Number of Fock matrices for DIIS
   end type scf_config_t

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
      procedure :: configure => xtb_configure
      procedure :: get_solvation_info => xtb_get_solvation_info
   end type xtb_config_t

   !============================================================================
   ! DFT Configuration (uses scf_config_t for SCF settings)
   !============================================================================
   type :: dft_config_t
      !! Configuration for Kohn-Sham DFT method
      !! Note: SCF settings (convergence, DIIS) come from scf_config_t
      character(len=32) :: functional = 'b3lyp'
         !! XC functional: "lda", "pbe", "b3lyp", "m06-2x", etc.

      ! Integration grid
      character(len=16) :: grid_type = 'medium'
         !! Grid quality: "coarse", "medium", "fine", "ultrafine"
      integer :: radial_points = 75
         !! Radial grid points per atom
      integer :: angular_points = 302
         !! Angular grid points (Lebedev)

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
   ! Correlation Configuration (shared by MP2, CC, etc.)
   !============================================================================
   type :: correlation_config_t
      !! Shared settings for all post-HF correlation methods
      real(dp) :: energy_convergence = 1.0e-8_dp
         !! Correlation energy convergence threshold

      ! Frozen core
      integer :: n_frozen_core = -1
         !! Number of frozen core orbitals (-1 = auto from elements)
      logical :: freeze_core = .true.
         !! Whether to freeze core orbitals

      ! Density fitting for correlation
      logical :: use_df = .true.
         !! Use density fitting (RI) for correlation integrals
      character(len=32) :: aux_basis = ''
         !! Auxiliary basis for RI (e.g., "cc-pvdz-ri", "cc-pvtz-ri")

      ! Local correlation
      logical :: use_local = .false.
         !! Use local correlation approximation
      character(len=16) :: local_type = 'dlpno'
         !! Local correlation type: "pno", "dlpno", "lmp2", "lno"
      real(dp) :: pno_threshold = 1.0e-7_dp
         !! PNO occupation threshold for truncation

      ! Spin-component scaling (for MP2)
      logical :: use_scs = .false.
         !! Use spin-component scaled MP2
      real(dp) :: scs_ss = 1.0_dp/3.0_dp
         !! Same-spin scaling factor (default: 1/3 for SCS-MP2)
      real(dp) :: scs_os = 1.2_dp
         !! Opposite-spin scaling factor (default: 6/5 for SCS-MP2)
   end type correlation_config_t

   !============================================================================
   ! Coupled Cluster Configuration
   !============================================================================
   type :: cc_config_t
      !! Coupled-cluster specific settings (CCSD, CCSD(T), CC2, CC3, etc.)
      integer :: max_iter = 100
         !! Maximum CC iterations
      real(dp) :: amplitude_convergence = 1.0e-7_dp
         !! T-amplitude convergence threshold

      ! Excitation level
      logical :: include_triples = .false.
         !! Include (T) triples correction
      logical :: perturbative_triples = .true.
         !! Use perturbative (T) vs full CCSDT

      ! DIIS for CC
      logical :: use_diis = .true.
         !! Use DIIS for amplitude equations
      integer :: diis_size = 8
         !! DIIS subspace size

      ! EOM-CC for excited states
      integer :: n_roots = 0
         !! Number of EOM-CC roots (0 = ground state only)
      character(len=8) :: eom_type = 'ee'
         !! EOM type: "ee" (excitation), "ip" (ionization), "ea" (attachment)
   end type cc_config_t

   !============================================================================
   ! F12 Explicitly Correlated Configuration
   !============================================================================
   type :: f12_config_t
      !! Settings for explicitly correlated F12 methods (MP2-F12, CCSD-F12, etc.)
      real(dp) :: geminal_exponent = 1.0_dp
         !! Slater-type geminal exponent (beta)
      character(len=8) :: ansatz = '3c'
         !! F12 ansatz: "3c", "3c(fix)", "2b", "2a"

      ! Auxiliary basis sets for F12
      character(len=32) :: cabs_basis = ''
         !! Complementary auxiliary basis (CABS) for RI
      character(len=32) :: optri_basis = ''
         !! Optional RI basis for F12 intermediates

      ! Approximations
      logical :: use_exponent_fit = .false.
         !! Fit geminal exponent to basis set
      logical :: scale_triples = .true.
         !! Apply F12 scaling to (T) correction
   end type f12_config_t

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

      !----- Shared configurations -----
      type(scf_config_t) :: scf
         !! Shared SCF settings (used by HF and DFT)
      type(correlation_config_t) :: corr
         !! Shared correlation settings (used by MP2, CC, etc.)

      !----- Method-specific configurations -----
      type(xtb_config_t) :: xtb
         !! XTB settings (GFN1, GFN2)
      type(dft_config_t) :: dft
         !! DFT-specific settings (functional, grid, dispersion)
      type(mcscf_config_t) :: mcscf
         !! MCSCF/CASSCF settings
      type(cc_config_t) :: cc
         !! Coupled-cluster specific settings (CCSD, CCSD(T), etc.)
      type(f12_config_t) :: f12
         !! F12 explicitly correlated settings

   contains
      procedure :: reset => config_reset
      procedure :: log_settings => config_log_settings
   end type method_config_t

contains

   pure logical function xtb_has_solvation(this)
      !! Check if solvation is configured for XTB
      class(xtb_config_t), intent(in) :: this
      xtb_has_solvation = len_trim(this%solvent) > 0 .or. this%dielectric > 0.0_dp
   end function xtb_has_solvation

   subroutine xtb_configure(this, use_cds, use_shift, dielectric, cpcm_nang, cpcm_rscale, solvent, solvation_model)
      !! Configure XTB solvation settings from driver configuration
      !!
      !! Sets up all XTB-specific parameters and applies default solvation model
      !! logic (alpb if solvent or dielectric is specified but no model given).
      class(xtb_config_t), intent(inout) :: this
      logical, intent(in) :: use_cds             !! Include CDS non-polar terms
      logical, intent(in) :: use_shift           !! Include solution state shift
      real(dp), intent(in) :: dielectric         !! Direct dielectric constant (-1 = use solvent lookup)
      integer, intent(in) :: cpcm_nang           !! Angular grid points for CPCM
      real(dp), intent(in) :: cpcm_rscale        !! Radii scaling for CPCM
      character(len=*), intent(in), optional :: solvent          !! Solvent name
      character(len=*), intent(in), optional :: solvation_model  !! Solvation model name

      this%use_cds = use_cds
      this%use_shift = use_shift
      this%dielectric = dielectric
      this%cpcm_nang = cpcm_nang
      this%cpcm_rscale = cpcm_rscale

      if (present(solvent)) this%solvent = solvent
      if (present(solvation_model)) then
         this%solvation_model = solvation_model
      else if (len_trim(this%solvent) > 0 .or. dielectric > 0.0_dp) then
         this%solvation_model = 'alpb'  ! Default solvation model
      end if
   end subroutine xtb_configure

   subroutine xtb_get_solvation_info(this, info_lines, n_lines)
      !! Get solvation configuration info for logging
      !!
      !! Returns array of info strings describing the solvation setup.
      !! If no solvation is configured, n_lines = 0.
      use pic_io, only: to_char
      class(xtb_config_t), intent(in) :: this
      character(len=128), intent(out) :: info_lines(4)  !! Up to 4 lines of info
      integer, intent(out) :: n_lines                    !! Number of lines populated

      n_lines = 0
      info_lines = ''

      if (.not. this%has_solvation()) return

      n_lines = 1
      if (trim(this%solvation_model) == 'cpcm') then
         if (this%dielectric > 0.0_dp) then
            info_lines(1) = "XTB solvation enabled: cpcm with dielectric = "//to_char(this%dielectric)
         else
            info_lines(1) = "XTB solvation enabled: cpcm with "//trim(this%solvent)
         end if
         n_lines = 3
         info_lines(2) = "  CPCM grid points (nang): "//to_char(this%cpcm_nang)
         info_lines(3) = "  CPCM radii scale: "//to_char(this%cpcm_rscale)
      else
         info_lines(1) = "XTB solvation enabled: "//trim(this%solvation_model)//" with "//trim(this%solvent)
         if (this%use_cds) then
            n_lines = n_lines + 1
            info_lines(n_lines) = "  CDS (non-polar) terms: enabled"
         end if
         if (this%use_shift) then
            n_lines = n_lines + 1
            info_lines(n_lines) = "  Solution state shift: enabled"
         end if
      end if
   end subroutine xtb_get_solvation_info

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

      ! SCF defaults (shared by HF and DFT)
      this%scf%max_iter = 100
      this%scf%energy_convergence = 1.0e-8_dp
      this%scf%density_convergence = 1.0e-6_dp
      this%scf%use_diis = .true.
      this%scf%diis_size = 8

      ! DFT-specific defaults
      this%dft%functional = 'b3lyp'
      this%dft%grid_type = 'medium'
      this%dft%radial_points = 75
      this%dft%angular_points = 302
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

      ! Correlation defaults (shared by MP2, CC, etc.)
      this%corr%energy_convergence = 1.0e-8_dp
      this%corr%n_frozen_core = -1
      this%corr%freeze_core = .true.
      this%corr%use_df = .true.
      this%corr%aux_basis = ''
      this%corr%use_local = .false.
      this%corr%local_type = 'dlpno'
      this%corr%pno_threshold = 1.0e-7_dp
      this%corr%use_scs = .false.
      this%corr%scs_ss = 1.0_dp/3.0_dp
      this%corr%scs_os = 1.2_dp

      ! Coupled-cluster defaults
      this%cc%max_iter = 100
      this%cc%amplitude_convergence = 1.0e-7_dp
      this%cc%include_triples = .false.
      this%cc%perturbative_triples = .true.
      this%cc%use_diis = .true.
      this%cc%diis_size = 8
      this%cc%n_roots = 0
      this%cc%eom_type = 'ee'

      ! F12 defaults
      this%f12%geminal_exponent = 1.0_dp
      this%f12%ansatz = '3c'
      this%f12%cabs_basis = ''
      this%f12%optri_basis = ''
      this%f12%use_exponent_fit = .false.
      this%f12%scale_triples = .true.
   end subroutine config_reset

   subroutine config_log_settings(this)
      !! Log method-specific settings based on method type
      !!
      !! Dispatches to appropriate logging based on the configured method.
      !! This allows the driver to log settings without knowing method details.
      use mqc_method_types, only: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2
      use pic_logger, only: logger => global_logger
      class(method_config_t), intent(in) :: this

      character(len=128) :: info_lines(4)
      integer :: n_lines, i

      ! Log XTB solvation settings if using XTB method
      select case (this%method_type)
      case (METHOD_TYPE_GFN1, METHOD_TYPE_GFN2)
         call this%xtb%get_solvation_info(info_lines, n_lines)
         do i = 1, n_lines
            call logger%info(trim(info_lines(i)))
         end do
      end select

   end subroutine config_log_settings

end module mqc_method_config
