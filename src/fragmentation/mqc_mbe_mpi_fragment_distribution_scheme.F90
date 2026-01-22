submodule(mqc_mbe_fragment_distribution_scheme) mpi_fragment_work_smod
   use mqc_error, only: ERROR_VALIDATION, ERROR_GENERIC
   implicit none

contains

   module subroutine do_fragment_work(fragment_idx, result, method_config, phys_frag, calc_type, world_comm)
      !! Process a single fragment for quantum chemistry calculation
      !!
      !! Performs energy and gradient calculation on a molecular fragment using
      !! the factory pattern to create a calculator from the provided method_config.
      !! Verbosity is controlled by the global logger level.

      use pic_logger, only: verbose_level

      integer(int64), intent(in) :: fragment_idx        !! Fragment index for identification
      type(calculation_result_t), intent(out) :: result  !! Computation results
      type(method_config_t), intent(in) :: method_config  !! Method configuration
      type(physical_fragment_t), intent(in), optional :: phys_frag  !! Fragment geometry
      integer(int32), intent(in) :: calc_type  !! Calculation type
      type(comm_t), intent(in), optional :: world_comm  !! MPI communicator for abort

      integer :: current_log_level  !! Current logger verbosity level
      logical :: is_verbose  !! Whether verbose output is enabled
      integer(int32) :: calc_type_local  !! Local copy of calc_type
      type(method_config_t) :: local_config  !! Local copy for verbose override
      class(qc_method_t), allocatable :: calculator  !! Polymorphic calculator instance

      calc_type_local = calc_type

      ! Query logger to determine verbosity
      call logger%configuration(level=current_log_level)
      is_verbose = (current_log_level >= verbose_level)

      ! Print fragment geometry if provided and verbose mode is enabled
      if (present(phys_frag)) then
         if (is_verbose) then
            call print_fragment_xyz(fragment_idx, phys_frag)
         end if

         ! Copy config and override verbose based on logger level
         local_config = method_config
         local_config%verbose = is_verbose

         ! Create calculator using factory
         calculator = create_method(local_config)

         ! Run the calculation using polymorphic dispatch
         select case (calc_type_local)
         case (CALC_TYPE_ENERGY)
            call calculator%calc_energy(phys_frag, result)
         case (CALC_TYPE_GRADIENT)
            call calculator%calc_gradient(phys_frag, result)
         case (CALC_TYPE_HESSIAN)
            call calculator%calc_hessian(phys_frag, result)
         case default
            call result%error%set(ERROR_VALIDATION, "Unknown calc_type: "//calc_type_to_string(calc_type_local))
            result%has_error = .true.
            return
         end select

         ! Check for calculation errors
         if (result%has_error) then
            call result%error%add_context("do_fragment_work:fragment_"//to_char(fragment_idx))
            return
         end if

         ! Copy fragment distance to result for JSON output
         result%distance = phys_frag%distance

         ! Cleanup
         deallocate (calculator)
      else
         ! For empty fragments, set energy to zero
         call result%energy%reset()
         result%has_energy = .true.
      end if
   end subroutine do_fragment_work

   module subroutine global_coordinator(resources, total_fragments, polymers, max_level, &
                                       node_leader_ranks, num_nodes, sys_geom, method_config, calc_type, bonds, json_data)
      !! Global coordinator for distributing fragments to node coordinators
      !! will act as a node coordinator for a single node calculation
      !! Uses int64 for total_fragments to handle large fragment counts that overflow int32.
      use mqc_json_output_types, only: json_output_data_t
      type(resources_t), intent(in) :: resources
      integer(int64), intent(in) :: total_fragments
      integer, intent(in) :: max_level, num_nodes
      integer, intent(in) :: polymers(:, :), node_leader_ranks(:)
      type(system_geometry_t), intent(in), optional :: sys_geom
      type(method_config_t), intent(in) :: method_config  !! Method configuration
      integer(int32), intent(in) :: calc_type
      type(bond_t), intent(in), optional :: bonds(:)
      type(json_output_data_t), intent(out), optional :: json_data  !! JSON output data

      type(timer_type) :: coord_timer
      integer(int64) :: current_fragment, results_received
      integer :: finished_nodes
      integer :: request_source, dummy_msg
      integer(int64) :: fragment_idx
      type(MPI_Status) :: status, local_status
      logical :: handling_local_workers
      logical :: has_pending
      integer(int32) :: calc_type_local

      ! For local workers
      integer :: local_finished_workers, local_dummy

      ! Storage for results
      type(calculation_result_t), allocatable :: results(:)
      integer(int64) :: worker_fragment_map(resources%mpi_comms%node_comm%size())
      integer :: worker_source

      ! MPI request handles for non-blocking operations
      type(request_t) :: req

      calc_type_local = calc_type

      current_fragment = total_fragments
      finished_nodes = 0
      local_finished_workers = 0
      handling_local_workers = (resources%mpi_comms%node_comm%size() > 1)
      results_received = 0_int64

      ! Allocate storage for results
      allocate (results(total_fragments))
      worker_fragment_map = 0

      call logger%verbose("Global coordinator starting with "//to_char(total_fragments)// &
                          " fragments for "//to_char(num_nodes)//" nodes")

      call coord_timer%start()
      do while (finished_nodes < num_nodes)

         ! PRIORITY 1: Check for incoming results from local workers
         ! This MUST be checked before sending new work to avoid race conditions
         if (handling_local_workers) then
            ! Keep checking for results until there are none pending
            do
           call iprobe(resources%mpi_comms%node_comm, MPI_ANY_SOURCE, TAG_WORKER_SCALAR_RESULT, has_pending, local_status)
               if (.not. has_pending) exit

               worker_source = local_status%MPI_SOURCE

               ! Safety check: worker should have a fragment assigned
               if (worker_fragment_map(worker_source) == 0) then
                  call logger%error("Received result from worker "//to_char(worker_source)// &
                                    " but no fragment was assigned!")
                  call abort_comm(resources%mpi_comms%world_comm, 1)
               end if

               ! Receive result and store it using the fragment index for this worker
            call result_irecv(results(worker_fragment_map(worker_source)), resources%mpi_comms%node_comm, worker_source, &
                                 TAG_WORKER_SCALAR_RESULT, req)
               call wait(req)

               ! Check for calculation errors from worker
               if (results(worker_fragment_map(worker_source))%has_error) then
                  call logger%error("Fragment "//to_char(worker_fragment_map(worker_source))// &
                                    " calculation failed: "// &
                                    results(worker_fragment_map(worker_source))%error%get_message())
                  call abort_comm(resources%mpi_comms%world_comm, 1)
               end if

               ! Clear the mapping since we've received the result
               worker_fragment_map(worker_source) = 0
               results_received = results_received + 1
               if (mod(results_received, max(1_int64, total_fragments/10)) == 0 .or. &
                   results_received == total_fragments) then
                  call logger%info("  Processed "//to_char(results_received)//"/"// &
                                   to_char(total_fragments)//" fragments ["// &
                                   to_char(coord_timer%get_elapsed_time())//" s]")
               end if
            end do
         end if

         ! PRIORITY 1b: Check for incoming results from remote node coordinators
         do
            call iprobe(resources%mpi_comms%world_comm, MPI_ANY_SOURCE, TAG_NODE_SCALAR_RESULT, has_pending, status)
            if (.not. has_pending) exit

            ! Receive fragment index and result from node coordinator
            ! TODO: serialize the data for better performance
            call irecv(resources%mpi_comms%world_comm, fragment_idx, status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT, req)
            call wait(req)
  call result_irecv(results(fragment_idx), resources%mpi_comms%world_comm, status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT, req)
            call wait(req)

            ! Check for calculation errors from node coordinator
            if (results(fragment_idx)%has_error) then
               call logger%error("Fragment "//to_char(fragment_idx)//" calculation failed: "// &
                                 results(fragment_idx)%error%get_message())
               call abort_comm(resources%mpi_comms%world_comm, 1)
            end if

            results_received = results_received + 1
            if (mod(results_received, max(1_int64, total_fragments/10)) == 0 .or. &
                results_received == total_fragments) then
               call logger%info("  Processed "//to_char(results_received)//"/"// &
                                to_char(total_fragments)//" fragments ["// &
                                to_char(coord_timer%get_elapsed_time())//" s]")
            end if
         end do

         ! PRIORITY 2: Remote node coordinator requests
         call iprobe(resources%mpi_comms%world_comm, MPI_ANY_SOURCE, TAG_NODE_REQUEST, has_pending, status)
         if (has_pending) then
            call irecv(resources%mpi_comms%world_comm, dummy_msg, status%MPI_SOURCE, TAG_NODE_REQUEST, req)
            call wait(req)
            request_source = status%MPI_SOURCE

            if (current_fragment >= 1) then
               call send_fragment_to_node(resources%mpi_comms%world_comm, current_fragment, polymers, request_source)
               current_fragment = current_fragment - 1
            else
               call isend(resources%mpi_comms%world_comm, -1, request_source, TAG_NODE_FINISH, req)
               call wait(req)
               finished_nodes = finished_nodes + 1
            end if
         end if

         ! PRIORITY 3: Local workers (shared memory) - send new work
         if (handling_local_workers .and. local_finished_workers < resources%mpi_comms%node_comm%size() - 1) then
            call iprobe(resources%mpi_comms%node_comm, MPI_ANY_SOURCE, TAG_WORKER_REQUEST, has_pending, local_status)
            if (has_pending) then
               ! Only process work request if this worker doesn't have pending results
               if (worker_fragment_map(local_status%MPI_SOURCE) == 0) then
                  call irecv(resources%mpi_comms%node_comm, local_dummy, local_status%MPI_SOURCE, TAG_WORKER_REQUEST, req)
                  call wait(req)

                  if (current_fragment >= 1) then
                     call send_fragment_to_worker(resources%mpi_comms%node_comm, current_fragment, polymers, &
                                                  local_status%MPI_SOURCE)
                     ! Track which fragment was sent to this worker
                     worker_fragment_map(local_status%MPI_SOURCE) = current_fragment
                     current_fragment = current_fragment - 1
                  else
                     call isend(resources%mpi_comms%node_comm, -1, local_status%MPI_SOURCE, TAG_WORKER_FINISH, req)
                     call wait(req)
                     local_finished_workers = local_finished_workers + 1
                  end if
               end if
               ! If worker still has pending results, skip the work request
               ! It will be processed on the next iteration after results are received
            end if
         end if

         ! Finalize local worker completion
         if (handling_local_workers .and. local_finished_workers >= resources%mpi_comms%node_comm%size() - 1 &
             .and. results_received >= total_fragments) then
            handling_local_workers = .false.
            if (num_nodes == 1) then
               finished_nodes = finished_nodes + 1
               call logger%debug("Manually incremented finished_nodes for self")
            else
               finished_nodes = finished_nodes + 1
               call logger%verbose("Global coordinator finished local workers")
            end if
         end if
      end do

      call logger%verbose("Global coordinator finished all fragments")
      call coord_timer%stop()
      call logger%info("Time to evaluate all fragments "//to_char(coord_timer%get_elapsed_time())//" s")
      block
         use mqc_result_types, only: mbe_result_t
         type(mbe_result_t) :: mbe_result

         ! Compute the many-body expansion
         call logger%info(" ")
         call logger%info("Computing Many-Body Expansion (MBE)...")
         call coord_timer%start()

         ! Allocate mbe_result components based on calc_type
         call mbe_result%allocate_dipole()  ! Always compute dipole
         if (calc_type_local == CALC_TYPE_HESSIAN) then
            if (.not. present(sys_geom)) then
               call logger%error("sys_geom required for Hessian calculation in global_coordinator")
               call abort_comm(resources%mpi_comms%world_comm, 1)
            end if
            call mbe_result%allocate_gradient(sys_geom%total_atoms)
            call mbe_result%allocate_hessian(sys_geom%total_atoms)
         else if (calc_type_local == CALC_TYPE_GRADIENT) then
            if (.not. present(sys_geom)) then
               call logger%error("sys_geom required for gradient calculation in global_coordinator")
               call abort_comm(resources%mpi_comms%world_comm, 1)
            end if
            call mbe_result%allocate_gradient(sys_geom%total_atoms)
         end if

         call compute_mbe(polymers, total_fragments, max_level, results, mbe_result, &
                          sys_geom, bonds, resources%mpi_comms%world_comm, json_data)
         call mbe_result%destroy()

         call coord_timer%stop()
         call logger%info("Time to compute MBE "//to_char(coord_timer%get_elapsed_time())//" s")

      end block

      ! Cleanup
      deallocate (results)
   end subroutine global_coordinator

   subroutine send_fragment_to_node(world_comm, fragment_idx, polymers, dest_rank)
      !! Send fragment data to remote node coordinator
      !! Uses int64 for fragment_idx to handle large fragment indices that overflow int32.
      type(comm_t), intent(in) :: world_comm
      integer(int64), intent(in) :: fragment_idx
      integer, intent(in) :: dest_rank
      integer, intent(in) :: polymers(:, :)
      integer :: fragment_size
      integer(int32) :: fragment_type
      integer, allocatable :: fragment_indices(:)
      type(request_t) :: req(4)
      integer(int64) :: fragment_idx_int64

      fragment_size = count(polymers(fragment_idx, :) > 0)
      allocate (fragment_indices(fragment_size))
      fragment_indices = polymers(fragment_idx, 1:fragment_size)

      ! Standard MBE always uses monomer indices
      fragment_type = FRAGMENT_TYPE_MONOMERS

      ! TODO: serialize the data for better performance
      fragment_idx_int64 = int(fragment_idx, kind=int64)
      call isend(world_comm, fragment_idx_int64, dest_rank, TAG_NODE_FRAGMENT, req(1))
      call isend(world_comm, fragment_type, dest_rank, TAG_NODE_FRAGMENT, req(2))
      call isend(world_comm, fragment_size, dest_rank, TAG_NODE_FRAGMENT, req(3))
      call isend(world_comm, fragment_indices, dest_rank, TAG_NODE_FRAGMENT, req(4))

      ! Wait for all sends to complete
      call wait(req(1))
      call wait(req(2))
      call wait(req(3))
      call wait(req(4))

      deallocate (fragment_indices)
   end subroutine send_fragment_to_node

   subroutine send_fragment_to_worker(node_comm, fragment_idx, polymers, dest_rank)
      !! Send fragment data to local worker
      !! Uses int64 for fragment_idx to handle large fragment indices that overflow int32.
      type(comm_t), intent(in) :: node_comm
      integer(int64), intent(in) :: fragment_idx
      integer, intent(in) :: dest_rank
      integer, intent(in) :: polymers(:, :)
      integer :: fragment_size
      integer(int32) :: fragment_type
      integer, allocatable :: fragment_indices(:)
      type(request_t) :: req(4)
      integer(int64) :: fragment_idx_int64

      fragment_size = count(polymers(fragment_idx, :) > 0)
      allocate (fragment_indices(fragment_size))
      fragment_indices = polymers(fragment_idx, 1:fragment_size)

      ! Standard MBE always uses monomer indices
      fragment_type = FRAGMENT_TYPE_MONOMERS

      ! TODO: serialize the data for better performance
      fragment_idx_int64 = int(fragment_idx, kind=int64)
      call isend(node_comm, fragment_idx_int64, dest_rank, TAG_WORKER_FRAGMENT, req(1))
      call isend(node_comm, fragment_type, dest_rank, TAG_WORKER_FRAGMENT, req(2))
      call isend(node_comm, fragment_size, dest_rank, TAG_WORKER_FRAGMENT, req(3))
      call isend(node_comm, fragment_indices, dest_rank, TAG_WORKER_FRAGMENT, req(4))

      ! Wait for all sends to complete
      call wait(req(1))
      call wait(req(2))
      call wait(req(3))
      call wait(req(4))

      deallocate (fragment_indices)
   end subroutine send_fragment_to_worker

   module subroutine node_coordinator(resources, method_config, calc_type)
      !! Node coordinator for distributing fragments to local workers
      !! Handles work requests and result collection from local workers
      type(resources_t), intent(in) :: resources
      type(method_config_t), intent(in) :: method_config  !! Method configuration (passed through to workers)
      integer(int32), intent(in) :: calc_type

      integer(int64) :: fragment_idx
      integer(int32) :: fragment_size, fragment_type, dummy_msg
      integer(int32) :: finished_workers
      integer(int32), allocatable :: fragment_indices(:)
      type(MPI_Status) :: status, global_status
      logical :: local_message_pending, more_fragments, has_result
      integer(int32) :: local_dummy

      ! For tracking worker-fragment mapping and collecting results
      integer(int64) :: worker_fragment_map(resources%mpi_comms%node_comm%size())
      integer(int32) :: worker_source
      type(calculation_result_t) :: worker_result

      ! MPI request handles for non-blocking operations
      type(request_t) :: req

      finished_workers = 0
      more_fragments = .true.
      dummy_msg = 0
      worker_fragment_map = 0

      do while (finished_workers < resources%mpi_comms%node_comm%size() - 1)

         ! PRIORITY 1: Check for incoming results from local workers
         call iprobe(resources%mpi_comms%node_comm, MPI_ANY_SOURCE, TAG_WORKER_SCALAR_RESULT, has_result, status)
         if (has_result) then
            worker_source = status%MPI_SOURCE

            ! Safety check: worker should have a fragment assigned
            if (worker_fragment_map(worker_source) == 0) then
               call logger%error("Node coordinator received result from worker "//to_char(worker_source)// &
                                 " but no fragment was assigned!")
               call abort_comm(resources%mpi_comms%world_comm, 1)
            end if

            ! Receive result from worker
            call result_irecv(worker_result, resources%mpi_comms%node_comm, worker_source, TAG_WORKER_SCALAR_RESULT, req)
            call wait(req)

            ! Check for calculation errors before forwarding
            if (worker_result%has_error) then
               call logger%error("Fragment "//to_char(worker_fragment_map(worker_source))// &
                                 " calculation failed on worker "//to_char(worker_source)//": "// &
                                 worker_result%error%get_message())
               call abort_comm(resources%mpi_comms%world_comm, 1)
            end if

            ! Forward results to global coordinator with fragment index
call isend(resources%mpi_comms%world_comm, worker_fragment_map(worker_source), 0, TAG_NODE_SCALAR_RESULT, req)  ! fragment_idx
            call wait(req)
  call result_isend(worker_result, resources%mpi_comms%world_comm, 0, TAG_NODE_SCALAR_RESULT, req)                ! result
            call wait(req)

            ! Clear the mapping
            worker_fragment_map(worker_source) = 0
         end if

         ! PRIORITY 2: Check for work requests from local workers
         call iprobe(resources%mpi_comms%node_comm, MPI_ANY_SOURCE, TAG_WORKER_REQUEST, local_message_pending, status)

         if (local_message_pending) then
            ! Only process work request if this worker doesn't have pending results
            if (worker_fragment_map(status%MPI_SOURCE) == 0) then
               call irecv(resources%mpi_comms%node_comm, local_dummy, status%MPI_SOURCE, TAG_WORKER_REQUEST, req)
               call wait(req)

               if (more_fragments) then
                  call isend(resources%mpi_comms%world_comm, dummy_msg, 0, TAG_NODE_REQUEST, req)
                  call wait(req)
                  call irecv(resources%mpi_comms%world_comm, fragment_idx, 0, MPI_ANY_TAG, req)
                  call wait(req, global_status)

                  if (global_status%MPI_TAG == TAG_NODE_FRAGMENT) then
                     ! Receive fragment type (0 = monomer indices, 1 = intersection atom list)
                     call irecv(resources%mpi_comms%world_comm, fragment_type, 0, TAG_NODE_FRAGMENT, req)
                     call wait(req)
                     call irecv(resources%mpi_comms%world_comm, fragment_size, 0, TAG_NODE_FRAGMENT, req)
                     call wait(req)
                     ! Note: must use blocking recv for allocatable arrays since size is unknown
                     allocate (fragment_indices(fragment_size))
                     call recv(resources%mpi_comms%world_comm, fragment_indices, 0, TAG_NODE_FRAGMENT, global_status)

                     ! Forward to worker
                     call isend(resources%mpi_comms%node_comm, fragment_idx, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)
                     call isend(resources%mpi_comms%node_comm, fragment_type, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)
                     call isend(resources%mpi_comms%node_comm, fragment_size, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)
                  call isend(resources%mpi_comms%node_comm, fragment_indices, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)

                     ! Track which fragment was sent to this worker
                     worker_fragment_map(status%MPI_SOURCE) = fragment_idx

                     deallocate (fragment_indices)
                  else
                     call isend(resources%mpi_comms%node_comm, -1, status%MPI_SOURCE, TAG_WORKER_FINISH, req)
                     call wait(req)
                     finished_workers = finished_workers + 1
                     more_fragments = .false.
                  end if
               else
                  call isend(resources%mpi_comms%node_comm, -1, status%MPI_SOURCE, TAG_WORKER_FINISH, req)
                  call wait(req)
                  finished_workers = finished_workers + 1
               end if
            end if
         end if
      end do
   end subroutine node_coordinator

   module subroutine node_worker(resources, sys_geom, method_config, calc_type, bonds)
      !! Node worker for processing fragments assigned by node coordinator
      use mqc_error, only: error_t
      type(resources_t), intent(in) :: resources
      type(system_geometry_t), intent(in), optional :: sys_geom
      type(method_config_t), intent(in) :: method_config  !! Method configuration
      integer(int32), intent(in) :: calc_type
      type(bond_t), intent(in), optional :: bonds(:)

      integer(int64) :: fragment_idx
      integer(int32) :: fragment_size, dummy_msg
      integer(int32) :: fragment_type  !! 0 = monomer (indices), 1 = intersection (atom list)
      integer(int32), allocatable :: fragment_indices(:)
      type(calculation_result_t) :: result
      type(MPI_Status) :: status
      type(physical_fragment_t) :: phys_frag
      type(error_t) :: error

      ! MPI request handles for non-blocking operations
      type(request_t) :: req

      dummy_msg = 0

      do
         call isend(resources%mpi_comms%node_comm, dummy_msg, 0, TAG_WORKER_REQUEST, req)
         call wait(req)
         call irecv(resources%mpi_comms%node_comm, fragment_idx, 0, MPI_ANY_TAG, req)
         call wait(req, status)

         select case (status%MPI_TAG)
         case (TAG_WORKER_FRAGMENT)
            ! Receive fragment type (0 = monomer indices, 1 = intersection atom list)
            call irecv(resources%mpi_comms%node_comm, fragment_type, 0, TAG_WORKER_FRAGMENT, req)
            call wait(req)
            call irecv(resources%mpi_comms%node_comm, fragment_size, 0, TAG_WORKER_FRAGMENT, req)
            call wait(req)
            ! Note: must use blocking recv for allocatable arrays since size is unknown
            allocate (fragment_indices(fragment_size))
            call recv(resources%mpi_comms%node_comm, fragment_indices, 0, TAG_WORKER_FRAGMENT, status)

            ! Build physical fragment based on type
            if (present(sys_geom)) then
               if (fragment_type == 0) then
                  ! Monomer: fragment_indices are monomer indices
                  call build_fragment_from_indices(sys_geom, fragment_indices, phys_frag, error, bonds)
               else
                  ! Intersection: fragment_indices are atom indices
                  call build_fragment_from_atom_list(sys_geom, fragment_indices, fragment_size, phys_frag, error, bonds)
               end if

               if (error%has_error()) then
                  call logger%error(error%get_full_trace())
                  call abort_comm(resources%mpi_comms%world_comm, 1)
               end if

               ! Process the chemistry fragment with physical geometry
          call do_fragment_work(fragment_idx, result, method_config, phys_frag, calc_type, resources%mpi_comms%world_comm)

               call phys_frag%destroy()
            else
               ! Process without physical geometry (old behavior)
               call do_fragment_work(fragment_idx, result, method_config, &
                                     calc_type=calc_type, world_comm=resources%mpi_comms%world_comm)
            end if

            ! Send result back to coordinator
            call result_isend(result, resources%mpi_comms%node_comm, 0, TAG_WORKER_SCALAR_RESULT, req)
            call wait(req)

            ! Clean up result
            call result%destroy()
            deallocate (fragment_indices)
         case (TAG_WORKER_FINISH)
            exit
         case default
            ! Unexpected MPI tag - this should not happen in normal operation
            call logger%error("Worker received unexpected MPI tag: "//to_char(status%MPI_TAG))
            call logger%error("Expected TAG_WORKER_FRAGMENT or TAG_WORKER_FINISH")
            call abort_comm(resources%mpi_comms%world_comm, 1)
         end select
      end do
   end subroutine node_worker

end submodule mpi_fragment_work_smod
