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
   use mqc_mbe_io, only: print_fragment_xyz
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
   use mqc_calculation_defaults, only: FRAGMENT_TYPE_MONOMERS, FRAGMENT_TYPE_ATOMS

   ! Method API imports
   use mqc_method_base, only: qc_method_t
   use mqc_method_config, only: method_config_t
   use mqc_method_factory, only: create_method
   use mqc_result_types, only: calculation_result_t, result_send, result_isend, result_recv, result_irecv
   use mqc_json_output_types, only: json_output_data_t, OUTPUT_MODE_UNFRAGMENTED, OUTPUT_MODE_MBE, OUTPUT_MODE_GMBE_PIE
   implicit none
   private

   ! Public interface - method_config is passed explicitly to all routines
   public :: do_fragment_work, global_coordinator, node_coordinator
   public :: serial_fragment_processor
   public :: node_worker, unfragmented_calculation, distributed_unfragmented_hessian

   interface
      module subroutine do_fragment_work(fragment_idx, result, method_config, phys_frag, calc_type, world_comm)
         implicit none
         integer(int64), intent(in) :: fragment_idx
         type(calculation_result_t), intent(out) :: result
         type(method_config_t), intent(in) :: method_config  !! Method configuration
         type(physical_fragment_t), intent(in), optional :: phys_frag
         integer(int32), intent(in) :: calc_type
         type(comm_t), intent(in), optional :: world_comm
      end subroutine do_fragment_work

      module subroutine global_coordinator(resources, total_fragments, polymers, max_level, &
                                       node_leader_ranks, num_nodes, sys_geom, method_config, calc_type, bonds, json_data)
         implicit none
         type(resources_t), intent(in) :: resources
         integer(int64), intent(in) :: total_fragments
         integer, intent(in) :: max_level, num_nodes
         integer, intent(in) :: polymers(:, :), node_leader_ranks(:)
         type(system_geometry_t), intent(in), optional :: sys_geom
         type(method_config_t), intent(in) :: method_config  !! Method configuration
         integer(int32), intent(in) :: calc_type
         type(bond_t), intent(in), optional :: bonds(:)
         type(json_output_data_t), intent(out), optional :: json_data  !! JSON output data
      end subroutine global_coordinator

      module subroutine node_coordinator(resources, method_config, calc_type)
         implicit none
         type(resources_t), intent(in) :: resources
         type(method_config_t), intent(in) :: method_config  !! Method configuration
         integer(int32), intent(in) :: calc_type
      end subroutine node_coordinator

      module subroutine serial_fragment_processor(total_fragments, polymers, max_level, sys_geom, &
                                                  method_config, calc_type, bonds, json_data)
         implicit none
         integer(int64), intent(in) :: total_fragments
         integer, intent(in) :: polymers(:, :), max_level
         type(system_geometry_t), intent(in) :: sys_geom
         type(method_config_t), intent(in) :: method_config  !! Method configuration
         integer(int32), intent(in) :: calc_type
         type(bond_t), intent(in), optional :: bonds(:)
         type(json_output_data_t), intent(out), optional :: json_data  !! JSON output data
      end subroutine serial_fragment_processor

      module subroutine node_worker(resources, sys_geom, method_config, calc_type, bonds)
         implicit none
         type(resources_t), intent(in) :: resources
         type(system_geometry_t), intent(in), optional :: sys_geom
         type(method_config_t), intent(in) :: method_config  !! Method configuration
         integer(int32), intent(in) :: calc_type
         type(bond_t), intent(in), optional :: bonds(:)
      end subroutine node_worker

      module subroutine unfragmented_calculation(sys_geom, method_config, calc_type, bonds, result_out, &
                                                 temperature, pressure, json_data)
         implicit none
         type(system_geometry_t), intent(in), optional :: sys_geom
         type(method_config_t), intent(in) :: method_config  !! Method configuration
         integer(int32), intent(in) :: calc_type
         type(bond_t), intent(in), optional :: bonds(:)
         type(calculation_result_t), intent(out), optional :: result_out
         real(dp), intent(in), optional :: temperature
         real(dp), intent(in), optional :: pressure
         type(json_output_data_t), intent(out), optional :: json_data
      end subroutine unfragmented_calculation

      module subroutine distributed_unfragmented_hessian(world_comm, sys_geom, method_config, driver_config, json_data)
         implicit none
         type(comm_t), intent(in) :: world_comm
         type(system_geometry_t), intent(in) :: sys_geom
         type(method_config_t), intent(in) :: method_config  !! Method configuration
         type(driver_config_t), intent(in), optional :: driver_config  !! Driver configuration
         type(json_output_data_t), intent(out), optional :: json_data  !! JSON output data
      end subroutine distributed_unfragmented_hessian

module subroutine hessian_coordinator(world_comm, sys_geom, method_config, displacement, temperature, pressure, json_data)
         implicit none
         type(comm_t), intent(in) :: world_comm
         type(system_geometry_t), intent(in) :: sys_geom
         type(method_config_t), intent(in) :: method_config  !! Method configuration
         real(dp), intent(in) :: displacement  !! Finite difference displacement (Bohr)
         real(dp), intent(in) :: temperature   !! Temperature for thermochemistry (K)
         real(dp), intent(in) :: pressure      !! Pressure for thermochemistry (atm)
         type(json_output_data_t), intent(out), optional :: json_data  !! JSON output data
      end subroutine hessian_coordinator

      module subroutine hessian_worker(world_comm, sys_geom, method_config, displacement)
         implicit none
         type(comm_t), intent(in) :: world_comm
         type(system_geometry_t), intent(in) :: sys_geom
         type(method_config_t), intent(in) :: method_config  !! Method configuration
         real(dp), intent(in) :: displacement  !! Finite difference displacement (Bohr)
      end subroutine hessian_worker

   end interface

end module mqc_mbe_fragment_distribution_scheme
