!! Factory for creating quantum chemistry method instances
module mqc_method_factory
   !! Provides centralized creation of quantum chemistry method instances.
   !! The factory pattern encapsulates method instantiation and configuration,
   !! making it easy to add new methods without modifying calling code.
   use pic_types, only: int32, dp
   use mqc_method_types, only: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2, METHOD_TYPE_HF, &
                               method_type_to_string
   use mqc_method_config, only: method_config_t
   use mqc_method_base, only: qc_method_t
   use mqc_method_hf, only: hf_method_t
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
      !!   config%method_type = METHOD_TYPE_GFN2
      !!   config%solvent = "water"
      !!   method = factory%create(config)
   contains
      procedure :: create => factory_create
   end type method_factory_t

contains

   function factory_create(this, config) result(method)
      !! Create a quantum chemistry method instance from configuration
      !!
      !! Instantiates the appropriate concrete method type based on
      !! config%method_type and configures it with the relevant options.
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

      case default
         error stop "Unknown method type in method_factory_t%create"
      end select
   end function factory_create

#ifndef MQC_WITHOUT_TBLITE
   subroutine configure_xtb(method, config)
      !! Configure an XTB method instance from method_config_t
      class(qc_method_t), intent(inout) :: method
      type(method_config_t), intent(in) :: config

      select type (m => method)
      type is (xtb_method_t)
         ! Core settings
         m%variant = method_type_to_string(config%method_type)
         m%verbose = config%verbose
         m%accuracy = real(config%accuracy, wp)

         ! Electronic temperature (convert K to Hartree)
         ! kt = T * k_B, where k_B = 3.166808578545117e-06 Hartree/K
         m%kt = real(config%electronic_temp, wp)*3.166808578545117e-06_wp

         ! Solvation settings
         if (config%has_solvation()) then
            m%solvent = trim(config%solvent)
            if (len_trim(config%solvation_model) > 0) then
               m%solvation_model = trim(config%solvation_model)
            else
               m%solvation_model = "alpb"  ! Default solvation model
            end if
            m%use_cds = config%use_cds
            m%use_shift = config%use_shift
            m%dielectric = real(config%dielectric, wp)
            m%cpcm_nang = config%cpcm_nang
            m%cpcm_rscale = real(config%cpcm_rscale, wp)
         end if
      end select
   end subroutine configure_xtb
#endif

   subroutine configure_hf(method, config)
      !! Configure a Hartree-Fock method instance from method_config_t
      class(qc_method_t), intent(inout) :: method
      type(method_config_t), intent(in) :: config

      select type (m => method)
      type is (hf_method_t)
         m%options%max_iter = config%max_scf_iter
         m%options%conv_tol = config%scf_convergence
         m%options%spherical = config%use_spherical
         m%options%verbose = config%verbose
         ! Note: basis_set would be used when HF is actually implemented
      end select
   end subroutine configure_hf

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
