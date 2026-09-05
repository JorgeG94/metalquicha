!! Factory for creating quantum chemistry method instances
module mqc_method_factory
   !! One place that turns a `method_config_t` into an allocated `qc_method_t`.
   use pic_types, only: int32, dp
   use mqc_method_types, only: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2, METHOD_TYPE_HF, &
                               METHOD_TYPE_DFT, METHOD_TYPE_MCSCF, METHOD_TYPE_MP2, &
                               METHOD_TYPE_CCSD, METHOD_TYPE_CCSD_T, &
                               method_type_to_string
   use mqc_method_config, only: scf_options_t, method_config_t
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
   public :: method_backend_built  !! Whether this build can run a method at all

   type :: method_factory_t
      !! Factory for creating quantum chemistry method instances
   contains
      procedure :: create => factory_create
   end type method_factory_t

contains

   pure function method_backend_built(method_type) result(built)
      !! Whether the backend this method needs was compiled into this binary
      !!
      !! **Ask this before calling the factory.** The factory returns a
      !! polymorphic allocatable and has no error to set, so a method whose
      !! backend is absent can only `ERROR STOP`.
      !!
      !! Only tblite is asked about here. The cenzontle and cuEST paths each have
      !! a stub that reports the missing build on the result; tblite has none,
      !! because the method type it backs is not compiled at all without it.
      integer, intent(in) :: method_type
      logical :: built

      built = .true.
#ifdef MQC_WITHOUT_TBLITE
      if (method_type == METHOD_TYPE_GFN1 .or. method_type == METHOD_TYPE_GFN2) then
         built = .false.
      end if
#endif
   end function method_backend_built

   function factory_create(this, config) result(method)
      !! Create a quantum chemistry method instance from configuration
      class(method_factory_t), intent(in) :: this
      type(method_config_t), intent(in) :: config
      class(qc_method_t), allocatable :: method

      ! Each branch configures a CONCRETE local and then source-allocates the
      ! polymorphic result from it. The more natural
      !     allocate(<type> :: method); call configure_x(method, config)
      ! passes a polymorphic allocatable function result to a CLASS(..),
      ! INTENT(INOUT) dummy, which gfortran 13.2.0 miscompiles: the callee
      ! updates a temporary and every assignment is silently discarded.
      select case (config%method_type)
#ifndef MQC_WITHOUT_TBLITE
      case (METHOD_TYPE_GFN1, METHOD_TYPE_GFN2)
         block
            type(xtb_method_t) :: xtb
            call configure_xtb(xtb, config)
            allocate (method, source=xtb)
         end block
#else
      case (METHOD_TYPE_GFN1, METHOD_TYPE_GFN2)
         error stop "XTB methods require tblite library (MQC_ENABLE_TBLITE)"
#endif

      case (METHOD_TYPE_HF)
         block
            type(hf_method_t) :: hf
            call configure_hf(hf, config)
            allocate (method, source=hf)
         end block

      case (METHOD_TYPE_MP2)
         block
            type(hf_method_t) :: hf
            call configure_hf(hf, config, with_mp2=.true.)
            allocate (method, source=hf)
         end block

         ! Both spellings land on the same method object. Which of the two ran
         ! is carried by config%cc%include_triples, which the adapter set from
         ! the method type.
      case (METHOD_TYPE_CCSD, METHOD_TYPE_CCSD_T)
         block
            type(hf_method_t) :: hf
            call configure_hf(hf, config, with_cc=.true.)
            allocate (method, source=hf)
         end block

      case (METHOD_TYPE_DFT)
         block
            type(dft_method_t) :: dft
            call configure_dft(dft, config)
            allocate (method, source=dft)
         end block

      case (METHOD_TYPE_MCSCF)
         block
            type(mcscf_method_t) :: mcscf
            call configure_mcscf(mcscf, config)
            allocate (method, source=mcscf)
         end block

      case default
         error stop "Unknown method type in method_factory_t%create"
      end select
   end function factory_create

#ifndef MQC_WITHOUT_TBLITE
   subroutine configure_xtb(m, config)
      !! Configure an XTB method instance from config%xtb
      type(xtb_method_t), intent(inout) :: m
      type(method_config_t), intent(in) :: config

      ! Core settings
      m%variant = method_type_to_string(config%method_type)
      m%verbose = config%verbose
      m%accuracy = real(config%xtb%accuracy, wp)
      m%max_iter = config%scf%max_iter
      m%allow_crap_scf = config%scf%allow_crap_scf

      ! Electronic temperature, Kelvin to Hartree: kt = T * k_B.
      ! TODO(mqc): a second, different Boltzmann constant. `KB_HARTREE` in
      ! `mqc_physical_constants` is 3.1668115634556e-6 and this literal is
      ! 3.166808578545117e-6, so the xTB electronic temperature is about one
      ! part in a million off the value every other module uses.
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
   end subroutine configure_xtb
#endif

   subroutine configure_scf(options, config)
      !! Fill everything a self-consistent-field method shares, from the deck
      !!
      !! Every field declared on `scf_options_t` is assigned here, once. A
      !! caller adds only what is specific to its own reference.
      class(scf_options_t), intent(inout) :: options
      type(method_config_t), intent(in) :: config

      options%basis_set = config%basis_set
      options%ecp_set = config%ecp_set
      options%spherical = config%use_spherical
      options%verbose = config%verbose
      options%hessian_displacement = config%hessian_displacement
      options%hessian_response_tol = config%hessian_response_tol
      options%hessian_response_max_iter = config%hessian_response_max_iter
      options%hessian_response_batch = config%hessian_response_batch
      options%device_rank = config%device_rank
      options%backend = config%backend
      options%freeze_core = config%corr%freeze_core
      options%n_frozen_core = config%corr%n_frozen_core
      options%cartesian = config%scf%cartesian
      options%aux_basis_set = config%scf%aux_basis_set
      options%aux_basis_named = config%scf%aux_basis_named
      options%density_fitting = config%scf%density_fitting
      options%unrestricted = config%scf%unrestricted
      options%guess = config%scf%guess
      if (allocated(config%scf%guess_steps)) options%guess_steps = config%scf%guess_steps
      options%max_iter = config%scf%max_iter
      options%allow_crap_scf = config%scf%allow_crap_scf
      options%energy_tol = config%scf%energy_convergence
      options%density_tol = config%scf%density_convergence
      options%grad_tol = config%scf%gradient_convergence
      options%level_shift = config%scf%level_shift
      options%linear_dependence = config%scf%linear_dependence
      options%use_diis = config%scf%use_diis
      options%diis_size = config%scf%diis_size
      options%incremental_fock = config%scf%incremental_fock
      options%accelerator = config%scf%accelerator
      options%convergence_metric = config%scf%convergence_metric
      options%pcm = config%pcm
      options%properties = config%properties
   end subroutine configure_scf

   subroutine configure_hf(m, config, with_mp2, with_cc)
      !! Configure a Hartree-Fock method instance from config%scf (shared SCF settings)
      type(hf_method_t), intent(inout) :: m
      type(method_config_t), intent(in) :: config
      logical, intent(in), optional :: with_cc
         !! Follow the reference with coupled cluster.
      logical, intent(in), optional :: with_mp2
         !! Follow the reference with MP2. MP2 and coupled cluster are
         !! dispatched onto this same method object rather than classes of
         !! their own: each is a Hartree-Fock calculation plus a correction
         !! built from its orbitals.

      call configure_scf(m%options, config)
      if (present(with_mp2)) m%options%run_mp2 = with_mp2
      if (present(with_cc)) m%options%run_cc = with_cc
      m%options%cc_triples = config%cc%include_triples
      m%options%cc_max_iter = config%cc%max_iter
      m%options%cc_tolerance = config%cc%amplitude_convergence
      ! Zero turns DIIS off, which is how the SCF spells the same thing.
      m%options%cc_diis_size = config%cc%diis_size
      m%options%cc_spin_adapted = config%cc%spin_adapted
      if (.not. config%cc%use_diis) m%options%cc_diis_size = 0
      ! From config%corr, not config%scf: the reference and the correlation
      ! treatment are configured separately and may disagree.
      m%options%corr_density_fitting = config%corr%use_df
      ! One and one unless scaling was asked for, so an unscaled run cannot
      ! pick up factors that happen to be sitting in the config.
      if (config%corr%use_scs) then
         m%options%scs_ss = config%corr%scs_ss
         m%options%scs_os = config%corr%scs_os
      end if
   end subroutine configure_hf

   subroutine configure_dft(m, config)
      !! Configure a DFT method instance from config%scf (shared) and config%dft (DFT-specific)
      type(dft_method_t), intent(inout) :: m
      type(method_config_t), intent(in) :: config

      call configure_scf(m%options, config)
      m%options%functional = config%dft%functional
      m%options%grid_type = config%dft%grid_type
      m%options%grid_level = config%dft%grid_level
      m%options%nlc_grid_level = config%dft%nlc_grid_level
      m%options%screening_tolerance = config%dft%screening_tolerance
      m%options%block_size = config%dft%block_size
      m%options%radial_points = config%dft%radial_points
      m%options%angular_points = config%dft%angular_points

      ! The auxiliary basis comes from the shared SCF settings, which is what
      ! `model.aux_basis` fills; a DFT-specific override wins where it is given.
      ! `density_fitting` itself arrives through `configure_scf`.
      ! TODO(mqc): the `else` branch is empty -- dead since the removal of the
      ! OR against `config%dft%use_density_fitting`.
      if (len_trim(config%dft%aux_basis_set) > 0) then
         m%options%aux_basis_set = config%dft%aux_basis_set
         m%options%aux_basis_named = .true.
      else
      end if

      ! TODO(mqc): nothing copies `config%corr` on this path, and
      ! `dft_options_t` declares no correlation fields, so a double hybrid's
      ! perturbative term never sees `keywords.correlation.use_df` or the
      ! spin-component scaling that Hartree-Fock reads from the same block.
      m%options%use_dispersion = config%dft%use_dispersion
      m%options%dispersion_type = config%dft%dispersion_type
   end subroutine configure_dft

   subroutine configure_mcscf(m, config)
      !! Configure a MCSCF method instance from config%mcscf
      type(mcscf_method_t), intent(inout) :: m
      type(method_config_t), intent(in) :: config

      ! Everything the reference SCF shares: a CASSCF starts from a closed-shell
      ! SCF, and a deck that tightened `keywords.scf` meant that one too.
      call configure_scf(m%options, config)

      ! Active space from config%mcscf
      if (allocated(config%mcscf%avas_orbitals)) then
         m%options%avas_orbitals = config%mcscf%avas_orbitals
      end if
      m%options%avas_threshold = config%mcscf%avas_threshold
      m%options%full_valence = config%mcscf%full_valence
      if (allocated(config%mcscf%ormas_subspaces)) then
         m%options%ormas_subspaces = config%mcscf%ormas_subspaces
         m%options%ormas_min_electrons = config%mcscf%ormas_min_electrons
         m%options%ormas_max_electrons = config%mcscf%ormas_max_electrons
      end if
      m%options%n_active_electrons = config%mcscf%n_active_electrons
      m%options%n_active_orbitals = config%mcscf%n_active_orbitals
      m%options%n_inactive_orbitals = config%mcscf%n_inactive_orbitals
      m%options%optimize_orbitals = config%mcscf%optimize_orbitals

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
   end subroutine configure_mcscf

   function create_method(config) result(method)
      !! Create a method without instantiating a factory first
      type(method_config_t), intent(in) :: config
      class(qc_method_t), allocatable :: method

      type(method_factory_t) :: factory

      ! Intentionally allocate(..., source=) rather than the more natural
      ! `method = factory%create(config)`. gfortran 13.2.0 miscompiles intrinsic
      ! assignment from a CLASS(...), ALLOCATABLE *function result*: the
      ! result's _vptr is not set up, so the assignment jumps through a null
      ! pointer and segfaults. Do not "simplify" this back.
      allocate (method, source=factory%create(config))
   end function create_method

end module mqc_method_factory
