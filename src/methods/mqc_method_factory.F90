!! Factory for creating quantum chemistry method instances
module mqc_method_factory
   !! Provides centralized creation of quantum chemistry method instances.
   !! The factory pattern encapsulates method instantiation and configuration,
   !! making it easy to add new methods without modifying calling code.
   use pic_types, only: int32, dp
   use mqc_method_types, only: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2, METHOD_TYPE_HF, &
                               METHOD_TYPE_DFT, METHOD_TYPE_MCSCF, method_type_to_string
   use mqc_method_config, only: method_config_t
   use mqc_method_base, only: qc_method_t
   use mqc_method_hf, only: hf_method_t
   use mqc_method_dft, only: dft_method_t
   use mqc_method_mcscf, only: mcscf_method_t
#ifndef MQC_WITHOUT_TBLITE
   use mqc_method_xtb, only: xtb_method_t
   use mctc_env, only: wp
#endif
   implicit none
   private

   public :: method_factory_t
   public :: create_method  !! Convenience function

   type :: method_factory_t
      !! Factory for creating quantum chemistry method instances
      !!
      !! Usage:
      !!   type(method_factory_t) :: factory
      !!   type(method_config_t) :: config
      !!   class(qc_method_t), allocatable :: method
      !!
      !!   config%method_type = METHOD_TYPE_DFT
      !!   config%basis_set = "cc-pvdz"
      !!   config%dft%functional = "pbe0"
      !!   method = factory%create(config)
   contains
      procedure :: create => factory_create
   end type method_factory_t

contains

   function factory_create(this, config) result(method)
      !! Create a quantum chemistry method instance from configuration
      !!
      !! Instantiates the appropriate concrete method type based on
      !! config%method_type and configures it from the nested config.
      class(method_factory_t), intent(in) :: this
      type(method_config_t), intent(in) :: config
      class(qc_method_t), allocatable :: method

      select case (config%method_type)
#ifndef MQC_WITHOUT_TBLITE
      case (METHOD_TYPE_GFN1, METHOD_TYPE_GFN2)
         allocate (xtb_method_t :: method)
         call configure_xtb(method, config)
#else
      case (METHOD_TYPE_GFN1, METHOD_TYPE_GFN2)
         error stop "XTB methods require tblite library (MQC_ENABLE_TBLITE)"
#endif

      case (METHOD_TYPE_HF)
         allocate (hf_method_t :: method)
         call configure_hf(method, config)

      case (METHOD_TYPE_DFT)
         allocate (dft_method_t :: method)
         call configure_dft(method, config)

      case (METHOD_TYPE_MCSCF)
         allocate (mcscf_method_t :: method)
         call configure_mcscf(method, config)

      case default
         error stop "Unknown method type in method_factory_t%create"
      end select
   end function factory_create

#ifndef MQC_WITHOUT_TBLITE
   subroutine configure_xtb(method, config)
      !! Configure an XTB method instance from config%xtb
      class(qc_method_t), intent(inout) :: method
      type(method_config_t), intent(in) :: config

      select type (m => method)
      type is (xtb_method_t)
         ! Core settings
         m%variant = method_type_to_string(config%method_type)
         m%verbose = config%verbose
         m%accuracy = real(config%xtb%accuracy, wp)

         ! Electronic temperature (convert K to Hartree)
         ! kt = T * k_B, where k_B = 3.166808578545117e-06 Hartree/K
         m%kt = real(config%xtb%electronic_temp, wp)*3.166808578545117e-06_wp

         ! Solvation settings from config%xtb
         if (config%xtb%has_solvation()) then
            m%solvent = trim(config%xtb%solvent)
            if (len_trim(config%xtb%solvation_model) > 0) then
               m%solvation_model = trim(config%xtb%solvation_model)
            else
               m%solvation_model = "alpb"  ! Default
            end if
            m%use_cds = config%xtb%use_cds
            m%use_shift = config%xtb%use_shift
            m%dielectric = real(config%xtb%dielectric, wp)
            m%cpcm_nang = config%xtb%cpcm_nang
            m%cpcm_rscale = real(config%xtb%cpcm_rscale, wp)
         end if
      end select
   end subroutine configure_xtb
#endif

   subroutine configure_hf(method, config)
      !! Configure a Hartree-Fock method instance from config%scf (shared SCF settings)
      class(qc_method_t), intent(inout) :: method
      type(method_config_t), intent(in) :: config

      select type (m => method)
      type is (hf_method_t)
         ! Common settings
         m%options%basis_set = config%basis_set
         m%options%spherical = config%use_spherical
         m%options%verbose = config%verbose

         ! SCF settings from shared config%scf
         m%options%max_iter = config%scf%max_iter
         m%options%conv_tol = config%scf%energy_convergence
         m%options%density_tol = config%scf%density_convergence
         m%options%use_diis = config%scf%use_diis
         m%options%diis_size = config%scf%diis_size
      end select
   end subroutine configure_hf

   subroutine configure_dft(method, config)
      !! Configure a DFT method instance from config%scf (shared) and config%dft (DFT-specific)
      class(qc_method_t), intent(inout) :: method
      type(method_config_t), intent(in) :: config

      select type (m => method)
      type is (dft_method_t)
         ! Common settings
         m%options%basis_set = config%basis_set
         m%options%spherical = config%use_spherical
         m%options%verbose = config%verbose

         ! SCF settings from shared config%scf
         m%options%max_iter = config%scf%max_iter
         m%options%energy_tol = config%scf%energy_convergence
         m%options%density_tol = config%scf%density_convergence
         m%options%use_diis = config%scf%use_diis
         m%options%diis_size = config%scf%diis_size

         ! DFT-specific from config%dft
         m%options%functional = config%dft%functional
         m%options%grid_type = config%dft%grid_type
         m%options%radial_points = config%dft%radial_points
         m%options%angular_points = config%dft%angular_points

         ! Density fitting
         m%options%use_density_fitting = config%dft%use_density_fitting
         m%options%aux_basis_set = config%dft%aux_basis_set

         ! Dispersion
         m%options%use_dispersion = config%dft%use_dispersion
         m%options%dispersion_type = config%dft%dispersion_type
      end select
   end subroutine configure_dft

   subroutine configure_mcscf(method, config)
      !! Configure a MCSCF method instance from config%mcscf
      class(qc_method_t), intent(inout) :: method
      type(method_config_t), intent(in) :: config

      select type (m => method)
      type is (mcscf_method_t)
         ! Common settings
         m%options%basis_set = config%basis_set
         m%options%spherical = config%use_spherical
         m%options%verbose = config%verbose

         ! Active space from config%mcscf
         m%options%n_active_electrons = config%mcscf%n_active_electrons
         m%options%n_active_orbitals = config%mcscf%n_active_orbitals
         m%options%n_inactive_orbitals = config%mcscf%n_inactive_orbitals

         ! State averaging
         m%options%n_states = config%mcscf%n_states
         if (allocated(config%mcscf%state_weights)) then
            if (allocated(m%options%state_weights)) deallocate (m%options%state_weights)
            allocate (m%options%state_weights(size(config%mcscf%state_weights)))
            m%options%state_weights = config%mcscf%state_weights
         end if

         ! Convergence
         m%options%max_macro_iter = config%mcscf%max_macro_iter
         m%options%max_micro_iter = config%mcscf%max_micro_iter
         m%options%orbital_tol = config%mcscf%orbital_convergence
         m%options%ci_tol = config%mcscf%ci_convergence

         ! PT2 corrections
         m%options%use_pt2 = config%mcscf%use_pt2
         m%options%pt2_type = config%mcscf%pt2_type
         m%options%ipea_shift = config%mcscf%ipea_shift
         m%options%imaginary_shift = config%mcscf%imaginary_shift
      end select
   end subroutine configure_mcscf

   function create_method(config) result(method)
      !! Convenience function to create a method without instantiating factory
      !!
      !! Usage:
      !!   use mqc_method_factory, only: create_method
      !!   method = create_method(config)
      type(method_config_t), intent(in) :: config
      class(qc_method_t), allocatable :: method

      type(method_factory_t) :: factory

      method = factory%create(config)
   end function create_method

end module mqc_method_factory
