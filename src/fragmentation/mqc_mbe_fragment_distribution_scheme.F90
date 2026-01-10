!! Many-Body Expansion (MBE) calculation module
module mqc_mbe_fragment_distribution_scheme
   !! Implements hierarchical many-body expansion for fragment-based quantum chemistry
   !! calculations with MPI parallelization and energy/gradient computation.
   use pic_types, only: int32, int64, dp
   use pic_timer, only: timer_type
   use pic_blas_interfaces, only: pic_gemm, pic_dot
   use pic_mpi_lib, only: comm_t, send, recv, isend, irecv, wait, iprobe, MPI_Status, &
                          request_t, MPI_ANY_SOURCE, MPI_ANY_TAG, abort_comm
   use mqc_resources, only: resources_t
   use pic_logger, only: logger => global_logger, verbose_level, info_level
   use pic_io, only: to_char
   use mqc_mbe_io, only: print_fragment_xyz, print_unfragmented_json
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   use mqc_mbe, only: compute_mbe
   use mqc_mpi_tags, only: TAG_WORKER_REQUEST, TAG_WORKER_FRAGMENT, TAG_WORKER_FINISH, &
                           TAG_WORKER_SCALAR_RESULT, &
                           TAG_NODE_REQUEST, TAG_NODE_FRAGMENT, TAG_NODE_FINISH, &
                           TAG_NODE_SCALAR_RESULT
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, build_fragment_from_indices, &
                                    build_fragment_from_atom_list, to_angstrom, check_duplicate_atoms
   use mqc_method_types, only: method_type_to_string
   use mqc_calc_types, only: calc_type_to_string, CALC_TYPE_ENERGY, CALC_TYPE_GRADIENT, CALC_TYPE_HESSIAN
   use mqc_config_parser, only: bond_t
   use mqc_config_adapter, only: driver_config_t

   ! Method API imports
#ifndef MQC_WITHOUT_TBLITE
   use mqc_method_xtb, only: xtb_method_t
#endif
   use mqc_result_types, only: calculation_result_t, result_send, result_isend, result_recv, result_irecv
   implicit none
   private

   ! XTB method options (module-level, set once at startup)
   type :: xtb_options_t
      character(len=:), allocatable :: solvent  !! Solvent name or empty for gas phase
      character(len=:), allocatable :: solvation_model  !! "alpb" (default), "gbsa", or "cpcm"
      logical :: use_cds = .true.               !! Include CDS non-polar terms (not for CPCM)
      logical :: use_shift = .true.             !! Include solution state shift (not for CPCM)
      ! CPCM-specific settings
      real(dp) :: dielectric = -1.0_dp          !! Direct dielectric constant (-1 = use solvent lookup)
      integer :: cpcm_nang = 110                !! Number of angular grid points for CPCM
      real(dp) :: cpcm_rscale = 1.0_dp          !! Radii scaling factor for CPCM
   end type xtb_options_t

   type(xtb_options_t), save :: xtb_options  !! Module-level XTB options

   ! Public interface
   public :: do_fragment_work, global_coordinator, node_coordinator
   public :: serial_fragment_processor
   public :: node_worker, unfragmented_calculation, distributed_unfragmented_hessian
   public :: set_xtb_options  !! Set XTB solvation options

   interface
      module subroutine do_fragment_work(fragment_idx, result, method, phys_frag, calc_type, world_comm)
         implicit none
         integer(int64), intent(in) :: fragment_idx
         type(calculation_result_t), intent(out) :: result
         integer(int32), intent(in) :: method
         type(physical_fragment_t), intent(in), optional :: phys_frag
         integer(int32), intent(in) :: calc_type
         type(comm_t), intent(in), optional :: world_comm
      end subroutine do_fragment_work

      module subroutine global_coordinator(resources, total_fragments, polymers, max_level, &
                                           node_leader_ranks, num_nodes, sys_geom, calc_type, bonds)
         implicit none
         type(resources_t), intent(in) :: resources
         integer(int64), intent(in) :: total_fragments
         integer, intent(in) :: max_level, num_nodes
         integer, intent(in) :: polymers(:, :), node_leader_ranks(:)
         type(system_geometry_t), intent(in), optional :: sys_geom
         integer(int32), intent(in) :: calc_type
         type(bond_t), intent(in), optional :: bonds(:)
      end subroutine global_coordinator

      module subroutine node_coordinator(resources, calc_type)
         implicit none
         type(resources_t), intent(in) :: resources
         integer(int32), intent(in) :: calc_type
      end subroutine node_coordinator

     module subroutine serial_fragment_processor(total_fragments, polymers, max_level, sys_geom, method, calc_type, bonds)
         implicit none
         integer(int64), intent(in) :: total_fragments
         integer, intent(in) :: polymers(:, :), max_level
         type(system_geometry_t), intent(in) :: sys_geom
         integer(int32), intent(in) :: method
         integer(int32), intent(in) :: calc_type
         type(bond_t), intent(in), optional :: bonds(:)
      end subroutine serial_fragment_processor

      module subroutine node_worker(resources, sys_geom, method, calc_type, bonds)
         implicit none
         type(resources_t), intent(in) :: resources
         type(system_geometry_t), intent(in), optional :: sys_geom
         integer(int32), intent(in) :: method
         integer(int32), intent(in) :: calc_type
         type(bond_t), intent(in), optional :: bonds(:)
      end subroutine node_worker

      module subroutine unfragmented_calculation(sys_geom, method, calc_type, bonds, result_out)
         implicit none
         type(system_geometry_t), intent(in), optional :: sys_geom
         integer(int32), intent(in) :: method
         integer(int32), intent(in) :: calc_type
         type(bond_t), intent(in), optional :: bonds(:)
         type(calculation_result_t), intent(out), optional :: result_out
      end subroutine unfragmented_calculation

      module subroutine distributed_unfragmented_hessian(world_comm, sys_geom, method, driver_config)
         implicit none
         type(comm_t), intent(in) :: world_comm
         type(system_geometry_t), intent(in) :: sys_geom
         integer(int32), intent(in) :: method
         type(driver_config_t), intent(in), optional :: driver_config  !! Driver configuration
      end subroutine distributed_unfragmented_hessian

      module subroutine hessian_coordinator(world_comm, sys_geom, method, displacement)
         implicit none
         type(comm_t), intent(in) :: world_comm
         type(system_geometry_t), intent(in) :: sys_geom
         integer(int32), intent(in) :: method
         real(dp), intent(in) :: displacement  !! Finite difference displacement (Bohr)
      end subroutine hessian_coordinator

      module subroutine hessian_worker(world_comm, sys_geom, method, displacement)
         implicit none
         type(comm_t), intent(in) :: world_comm
         type(system_geometry_t), intent(in) :: sys_geom
         integer(int32), intent(in) :: method
         real(dp), intent(in) :: displacement  !! Finite difference displacement (Bohr)
      end subroutine hessian_worker

   end interface

contains

   subroutine set_xtb_options(solvent, solvation_model, use_cds, use_shift, dielectric, cpcm_nang, cpcm_rscale)
      !! Set module-level XTB solvation options
      !! Call this once at startup before running calculations
      character(len=*), intent(in), optional :: solvent  !! Solvent name or empty for gas phase
      character(len=*), intent(in), optional :: solvation_model  !! "alpb" (default), "gbsa", or "cpcm"
      logical, intent(in), optional :: use_cds           !! Include CDS non-polar terms (not for CPCM)
      logical, intent(in), optional :: use_shift         !! Include solution state shift (not for CPCM)
      real(dp), intent(in), optional :: dielectric       !! Direct dielectric constant for CPCM
      integer, intent(in), optional :: cpcm_nang         !! Angular grid points for CPCM
      real(dp), intent(in), optional :: cpcm_rscale      !! Radii scaling for CPCM

      if (present(solvent) .and. len_trim(solvent) > 0) then
         xtb_options%solvent = trim(solvent)
      end if
      if (present(solvation_model) .and. len_trim(solvation_model) > 0) then
         xtb_options%solvation_model = trim(solvation_model)
      end if
      if (present(use_cds)) then
         xtb_options%use_cds = use_cds
      end if
      if (present(use_shift)) then
         xtb_options%use_shift = use_shift
      end if
      if (present(dielectric)) then
         xtb_options%dielectric = dielectric
      end if
      if (present(cpcm_nang)) then
         xtb_options%cpcm_nang = cpcm_nang
      end if
      if (present(cpcm_rscale)) then
         xtb_options%cpcm_rscale = cpcm_rscale
      end if
   end subroutine set_xtb_options

end module mqc_mbe_fragment_distribution_scheme
