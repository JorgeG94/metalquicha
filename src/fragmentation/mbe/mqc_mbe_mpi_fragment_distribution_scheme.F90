submodule(mqc_mbe_fragment_distribution_scheme) mpi_fragment_work_smod
   use mqc_error, only: ERROR_VALIDATION, ERROR_GENERIC
   implicit none

contains

   subroutine queue_init_from_list(queue, ids)
      !! Initialize queue from a list of fragment indices.
      type(fragment_queue_t), intent(out) :: queue
      integer(int64), intent(in) :: ids(:)

      queue%count = size(ids, kind=int64)
      if (queue%count > 0) then
         allocate (queue%ids(queue%count))
         queue%ids = ids
      end if
      queue%head = 1_int64
   end subroutine queue_init_from_list

   subroutine queue_pop(queue, fragment_idx, has_item)
      !! Pop the next fragment index from the queue.
      type(fragment_queue_t), intent(inout) :: queue
      integer(int64), intent(out) :: fragment_idx
      logical, intent(out) :: has_item

      if (queue%head > queue%count) then
         fragment_idx = -1_int64
         has_item = .false.
         return
      end if

      fragment_idx = queue%ids(queue%head)
      queue%head = queue%head + 1_int64
      has_item = .true.
   end subroutine queue_pop

   subroutine queue_destroy(queue)
      type(fragment_queue_t), intent(inout) :: queue
      if (allocated(queue%ids)) deallocate (queue%ids)
      queue%head = 1_int64
      queue%count = 0_int64
   end subroutine queue_destroy

   subroutine build_fragment_payload(fragment_idx, polymers, fragment_type, fragment_size, fragment_indices)
      !! Build fragment payload arrays for a given fragment index.
      !! Caller owns the returned fragment_indices and must deallocate it.
      integer(int64), intent(in) :: fragment_idx
      integer, intent(in) :: polymers(:, :)
      integer(int32), intent(out) :: fragment_type
      integer(int32), intent(out) :: fragment_size
      integer, allocatable, intent(out) :: fragment_indices(:)

      fragment_size = count(polymers(fragment_idx, :) > 0)
      allocate (fragment_indices(fragment_size))
      fragment_indices = polymers(fragment_idx, 1:fragment_size)

      ! Standard MBE always uses monomer indices
      fragment_type = FRAGMENT_TYPE_MONOMERS
   end subroutine build_fragment_payload

   subroutine send_fragment_payload(comm, tag, fragment_idx, polymers, dest_rank)
      !! Send fragment payload over the specified communicator/tag.
      !! Uses int64 for fragment_idx to handle large fragment indices.
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: tag
      integer(int64), intent(in) :: fragment_idx
      integer, intent(in) :: dest_rank
      integer, intent(in) :: polymers(:, :)
      integer(int32) :: fragment_size, fragment_type
      integer, allocatable :: fragment_indices(:)
      type(request_t) :: req(4)
      integer(int64) :: fragment_idx_int64

      call build_fragment_payload(fragment_idx, polymers, fragment_type, fragment_size, fragment_indices)

      fragment_idx_int64 = int(fragment_idx, kind=int64)
      call isend(comm, fragment_idx_int64, dest_rank, tag, req(1))
      call isend(comm, fragment_type, dest_rank, tag, req(2))
      call isend(comm, fragment_size, dest_rank, tag, req(3))
      call isend(comm, fragment_indices, dest_rank, tag, req(4))

      call wait(req(1))
      call wait(req(2))
      call wait(req(3))
      call wait(req(4))

      deallocate (fragment_indices)
   end subroutine send_fragment_payload

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

   module subroutine global_coordinator(ctx, json_data)
      !! Global coordinator for distributing fragments to node coordinators
      !! will act as a node coordinator for a single node calculation
      !! Uses int64 for total_fragments to handle large fragment counts that overflow int32.
      !! Bond connectivity is accessed via ctx%sys_geom%bonds
      use mqc_json_output_types, only: json_output_data_t
      use mqc_many_body_expansion, only: mbe_context_t
      class(*), intent(in) :: ctx
      type(json_output_data_t), intent(out), optional :: json_data  !! JSON output data

      ! Cast to mbe_context_t via select type
      select type (ctx)
      type is (mbe_context_t)
         call global_coordinator_impl(ctx, json_data)
      class default
         call logger%error("global_coordinator: expected mbe_context_t")
         call abort_comm(comm_world(), 1)
      end select
   end subroutine global_coordinator

   subroutine global_coordinator_impl(ctx, json_data)
      !! Internal implementation of global_coordinator with typed context
      use mqc_json_output_types, only: json_output_data_t
      use mqc_many_body_expansion, only: mbe_context_t
      type(mbe_context_t), intent(in) :: ctx
      type(json_output_data_t), intent(out), optional :: json_data  !! JSON output data

      type(timer_type) :: coord_timer
      integer(int64) :: results_received
      integer :: finished_nodes
      integer(int64) :: fragment_idx
      type(MPI_Status) :: status, local_status
      logical :: handling_local_workers
      logical :: has_pending
      integer(int32) :: calc_type_local

      ! For local workers
      integer :: local_finished_workers, local_dummy

      ! Storage for results
      type(calculation_result_t), allocatable :: results(:)
      integer(int64) :: worker_fragment_map(ctx%resources%mpi_comms%node_comm%size())
      integer :: worker_source
      type(fragment_queue_t) :: fragment_queue
      logical :: has_fragment

      ! MPI request handles for non-blocking operations
      type(request_t) :: req

      calc_type_local = ctx%calc_type

      finished_nodes = 0
      local_finished_workers = 0
      handling_local_workers = (ctx%resources%mpi_comms%node_comm%size() > 1)
      results_received = 0_int64

      ! Allocate storage for results
      allocate (results(ctx%total_fragments))
      worker_fragment_map = 0
      block
         integer(int64), allocatable :: temp_ids(:)
         integer(int64) :: i

         allocate (temp_ids(ctx%total_fragments))
         do i = 1_int64, ctx%total_fragments
            temp_ids(i) = ctx%total_fragments - i + 1_int64
         end do
         call queue_init_from_list(fragment_queue, temp_ids)
         deallocate (temp_ids)
      end block

      call logger%verbose("Global coordinator starting with "//to_char(ctx%total_fragments)// &
                          " fragments for "//to_char(ctx%num_nodes)//" nodes")

      call coord_timer%start()
      do while (finished_nodes < ctx%num_nodes)

         ! PRIORITY 1: Check for incoming results from local workers
         ! This MUST be checked before sending new work to avoid race conditions
         if (handling_local_workers) then
            call handle_local_worker_results(ctx, worker_fragment_map, results, results_received, coord_timer)
         end if

         ! PRIORITY 1b: Check for incoming results from remote node coordinators
         call handle_node_results(ctx, results, results_received, coord_timer)

         ! PRIORITY 2: Remote node coordinator requests
         call handle_node_requests(ctx, fragment_queue, finished_nodes)

         ! PRIORITY 3: Local workers (shared memory) - send new work
         if (handling_local_workers .and. local_finished_workers < ctx%resources%mpi_comms%node_comm%size() - 1) then
            call handle_local_worker_requests(ctx, fragment_queue, worker_fragment_map, local_finished_workers)
         end if

         ! Finalize local worker completion
         if (handling_local_workers .and. local_finished_workers >= ctx%resources%mpi_comms%node_comm%size() - 1 &
             .and. results_received >= ctx%total_fragments) then
            handling_local_workers = .false.
            if (ctx%num_nodes == 1) then
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
            if (.not. ctx%has_geometry()) then
               call logger%error("sys_geom required for Hessian calculation in global_coordinator")
               call abort_comm(ctx%resources%mpi_comms%world_comm, 1)
            end if
            call mbe_result%allocate_gradient(ctx%sys_geom%total_atoms)
            call mbe_result%allocate_hessian(ctx%sys_geom%total_atoms)
         else if (calc_type_local == CALC_TYPE_GRADIENT) then
            if (.not. ctx%has_geometry()) then
               call logger%error("sys_geom required for gradient calculation in global_coordinator")
               call abort_comm(ctx%resources%mpi_comms%world_comm, 1)
            end if
            call mbe_result%allocate_gradient(ctx%sys_geom%total_atoms)
         end if

         call compute_mbe(ctx%polymers, ctx%total_fragments, ctx%max_level, results, mbe_result, &
                          ctx%sys_geom, ctx%resources%mpi_comms%world_comm, json_data)
         call mbe_result%destroy()

         call coord_timer%stop()
         call logger%info("Time to compute MBE "//to_char(coord_timer%get_elapsed_time())//" s")

      end block

      ! Cleanup
      call queue_destroy(fragment_queue)
      deallocate (results)
   end subroutine global_coordinator_impl

   subroutine handle_node_requests(ctx, fragment_queue, finished_nodes)
      !! Handle a single pending node coordinator request, if any.
      use mqc_many_body_expansion, only: mbe_context_t
      type(mbe_context_t), intent(in) :: ctx
      type(fragment_queue_t), intent(inout) :: fragment_queue
      integer, intent(inout) :: finished_nodes

      integer :: request_source, dummy_msg
      integer(int64) :: fragment_idx
      type(MPI_Status) :: status
      logical :: has_pending, has_fragment
      type(request_t) :: req

      call iprobe(ctx%resources%mpi_comms%world_comm, MPI_ANY_SOURCE, TAG_NODE_REQUEST, has_pending, status)
      if (.not. has_pending) return

      call irecv(ctx%resources%mpi_comms%world_comm, dummy_msg, status%MPI_SOURCE, TAG_NODE_REQUEST, req)
      call wait(req)
      request_source = status%MPI_SOURCE

      call queue_pop(fragment_queue, fragment_idx, has_fragment)
      if (has_fragment) then
         call send_fragment_to_node(ctx%resources%mpi_comms%world_comm, fragment_idx, ctx%polymers, request_source)
      else
         call isend(ctx%resources%mpi_comms%world_comm, -1, request_source, TAG_NODE_FINISH, req)
         call wait(req)
         finished_nodes = finished_nodes + 1
      end if
   end subroutine handle_node_requests

   subroutine handle_local_worker_requests(ctx, fragment_queue, worker_fragment_map, local_finished_workers)
      !! Handle a single pending local worker request, if any.
      use mqc_many_body_expansion, only: mbe_context_t
      type(mbe_context_t), intent(in) :: ctx
      type(fragment_queue_t), intent(inout) :: fragment_queue
      integer(int64), intent(inout) :: worker_fragment_map(:)
      integer, intent(inout) :: local_finished_workers

      integer(int64) :: fragment_idx
      integer(int32) :: local_dummy
      type(MPI_Status) :: local_status
      logical :: has_pending, has_fragment
      type(request_t) :: req

      call iprobe(ctx%resources%mpi_comms%node_comm, MPI_ANY_SOURCE, TAG_WORKER_REQUEST, has_pending, local_status)
      if (.not. has_pending) return

      if (worker_fragment_map(local_status%MPI_SOURCE) /= 0) return

      call irecv(ctx%resources%mpi_comms%node_comm, local_dummy, local_status%MPI_SOURCE, TAG_WORKER_REQUEST, req)
      call wait(req)

      call queue_pop(fragment_queue, fragment_idx, has_fragment)
      if (has_fragment) then
      call send_fragment_to_worker(ctx%resources%mpi_comms%node_comm, fragment_idx, ctx%polymers, local_status%MPI_SOURCE)
         worker_fragment_map(local_status%MPI_SOURCE) = fragment_idx
      else
         call isend(ctx%resources%mpi_comms%node_comm, -1, local_status%MPI_SOURCE, TAG_WORKER_FINISH, req)
         call wait(req)
         local_finished_workers = local_finished_workers + 1
      end if
   end subroutine handle_local_worker_requests

   subroutine handle_local_worker_results(ctx, worker_fragment_map, results, results_received, coord_timer)
      !! Drain results from local workers and update tracking state.
      use mqc_many_body_expansion, only: mbe_context_t
      type(mbe_context_t), intent(in) :: ctx
      integer(int64), intent(inout) :: worker_fragment_map(:)
      type(calculation_result_t), intent(inout) :: results(:)
      integer(int64), intent(inout) :: results_received
      type(timer_type), intent(in) :: coord_timer

      type(MPI_Status) :: local_status
      logical :: has_pending
      integer :: worker_source
      type(request_t) :: req

      do
       call iprobe(ctx%resources%mpi_comms%node_comm, MPI_ANY_SOURCE, TAG_WORKER_SCALAR_RESULT, has_pending, local_status)
         if (.not. has_pending) exit

         worker_source = local_status%MPI_SOURCE

         if (worker_fragment_map(worker_source) == 0) then
            call logger%error("Received result from worker "//to_char(worker_source)// &
                              " but no fragment was assigned!")
            call abort_comm(ctx%resources%mpi_comms%world_comm, 1)
         end if

        call result_irecv(results(worker_fragment_map(worker_source)), ctx%resources%mpi_comms%node_comm, worker_source, &
                           TAG_WORKER_SCALAR_RESULT, req)
         call wait(req)

         if (results(worker_fragment_map(worker_source))%has_error) then
            call logger%error("Fragment "//to_char(worker_fragment_map(worker_source))// &
                              " calculation failed: "// &
                              results(worker_fragment_map(worker_source))%error%get_message())
            call abort_comm(ctx%resources%mpi_comms%world_comm, 1)
         end if

         worker_fragment_map(worker_source) = 0
         results_received = results_received + 1
         if (mod(results_received, max(1_int64, ctx%total_fragments/10)) == 0 .or. &
             results_received == ctx%total_fragments) then
            call logger%info("  Processed "//to_char(results_received)//"/"// &
                             to_char(ctx%total_fragments)//" fragments ["// &
                             to_char(coord_timer%get_elapsed_time())//" s]")
         end if
      end do
   end subroutine handle_local_worker_results

   subroutine handle_node_results(ctx, results, results_received, coord_timer)
      !! Drain results from remote node coordinators and update tracking state.
      use mqc_many_body_expansion, only: mbe_context_t
      type(mbe_context_t), intent(in) :: ctx
      type(calculation_result_t), intent(inout) :: results(:)
      integer(int64), intent(inout) :: results_received
      type(timer_type), intent(in) :: coord_timer

      integer(int64) :: fragment_idx
      type(MPI_Status) :: status
      logical :: has_pending
      type(request_t) :: req

      do
         call iprobe(ctx%resources%mpi_comms%world_comm, MPI_ANY_SOURCE, TAG_NODE_SCALAR_RESULT, has_pending, status)
         if (.not. has_pending) exit

         call irecv(ctx%resources%mpi_comms%world_comm, fragment_idx, status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT, req)
         call wait(req)
         call result_irecv(results(fragment_idx), ctx%resources%mpi_comms%world_comm, status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT, req)
         call wait(req)

         if (results(fragment_idx)%has_error) then
            call logger%error("Fragment "//to_char(fragment_idx)//" calculation failed: "// &
                              results(fragment_idx)%error%get_message())
            call abort_comm(ctx%resources%mpi_comms%world_comm, 1)
         end if

         results_received = results_received + 1
         if (mod(results_received, max(1_int64, ctx%total_fragments/10)) == 0 .or. &
             results_received == ctx%total_fragments) then
            call logger%info("  Processed "//to_char(results_received)//"/"// &
                             to_char(ctx%total_fragments)//" fragments ["// &
                             to_char(coord_timer%get_elapsed_time())//" s]")
         end if
      end do
   end subroutine handle_node_results

   subroutine send_fragment_to_node(world_comm, fragment_idx, polymers, dest_rank)
      !! Send fragment data to remote node coordinator
      !! Uses int64 for fragment_idx to handle large fragment indices that overflow int32.
      type(comm_t), intent(in) :: world_comm
      integer(int64), intent(in) :: fragment_idx
      integer, intent(in) :: dest_rank
      integer, intent(in) :: polymers(:, :)
      call send_fragment_payload(world_comm, TAG_NODE_FRAGMENT, fragment_idx, polymers, dest_rank)
   end subroutine send_fragment_to_node

   subroutine send_fragment_to_worker(node_comm, fragment_idx, polymers, dest_rank)
      !! Send fragment data to local worker
      !! Uses int64 for fragment_idx to handle large fragment indices that overflow int32.
      type(comm_t), intent(in) :: node_comm
      integer(int64), intent(in) :: fragment_idx
      integer, intent(in) :: dest_rank
      integer, intent(in) :: polymers(:, :)
      call send_fragment_payload(node_comm, TAG_WORKER_FRAGMENT, fragment_idx, polymers, dest_rank)
   end subroutine send_fragment_to_worker

   module subroutine node_coordinator(ctx)
      !! Node coordinator for distributing fragments to local workers
      !! Handles work requests and result collection from local workers
      use mqc_many_body_expansion, only: many_body_expansion_t
      class(*), intent(in) :: ctx

      ! Cast to many_body_expansion_t via select type
      select type (ctx)
      class is (many_body_expansion_t)
         call node_coordinator_impl(ctx)
      class default
         call logger%error("node_coordinator: expected many_body_expansion_t")
         call abort_comm(comm_world(), 1)
      end select
   end subroutine node_coordinator

   subroutine node_coordinator_impl(ctx)
      !! Internal implementation of node_coordinator with typed context
      use mqc_many_body_expansion, only: many_body_expansion_t
      class(many_body_expansion_t), intent(in) :: ctx

      integer(int64) :: fragment_idx
      integer(int32) :: fragment_size, fragment_type, dummy_msg
      integer(int32) :: finished_workers
      integer(int32), allocatable :: fragment_indices(:)
      type(MPI_Status) :: status, global_status
      logical :: local_message_pending, more_fragments, has_result
      integer(int32) :: local_dummy

      ! For tracking worker-fragment mapping and collecting results
      integer(int64) :: worker_fragment_map(ctx%resources%mpi_comms%node_comm%size())
      integer(int32) :: worker_source
      type(calculation_result_t) :: worker_result

      ! MPI request handles for non-blocking operations
      type(request_t) :: req

      finished_workers = 0
      more_fragments = .true.
      dummy_msg = 0
      worker_fragment_map = 0

      do while (finished_workers < ctx%resources%mpi_comms%node_comm%size() - 1)

         ! PRIORITY 1: Check for incoming results from local workers
         call iprobe(ctx%resources%mpi_comms%node_comm, MPI_ANY_SOURCE, TAG_WORKER_SCALAR_RESULT, has_result, status)
         if (has_result) then
            worker_source = status%MPI_SOURCE

            ! Safety check: worker should have a fragment assigned
            if (worker_fragment_map(worker_source) == 0) then
               call logger%error("Node coordinator received result from worker "//to_char(worker_source)// &
                                 " but no fragment was assigned!")
               call abort_comm(ctx%resources%mpi_comms%world_comm, 1)
            end if

            ! Receive result from worker
         call result_irecv(worker_result, ctx%resources%mpi_comms%node_comm, worker_source, TAG_WORKER_SCALAR_RESULT, req)
            call wait(req)

            ! Check for calculation errors before forwarding
            if (worker_result%has_error) then
               call logger%error("Fragment "//to_char(worker_fragment_map(worker_source))// &
                                 " calculation failed on worker "//to_char(worker_source)//": "// &
                                 worker_result%error%get_message())
               call abort_comm(ctx%resources%mpi_comms%world_comm, 1)
            end if

            ! Forward results to global coordinator with fragment index
call isend(ctx%resources%mpi_comms%world_comm, worker_fragment_map(worker_source), 0, TAG_NODE_SCALAR_RESULT, req)  ! fragment_idx
            call wait(req)
call result_isend(worker_result, ctx%resources%mpi_comms%world_comm, 0, TAG_NODE_SCALAR_RESULT, req)                ! result
            call wait(req)

            ! Clear the mapping
            worker_fragment_map(worker_source) = 0
         end if

         ! PRIORITY 2: Check for work requests from local workers
         call iprobe(ctx%resources%mpi_comms%node_comm, MPI_ANY_SOURCE, TAG_WORKER_REQUEST, local_message_pending, status)

         if (local_message_pending) then
            ! Only process work request if this worker doesn't have pending results
            if (worker_fragment_map(status%MPI_SOURCE) == 0) then
               call irecv(ctx%resources%mpi_comms%node_comm, local_dummy, status%MPI_SOURCE, TAG_WORKER_REQUEST, req)
               call wait(req)

               if (more_fragments) then
                  call isend(ctx%resources%mpi_comms%world_comm, dummy_msg, 0, TAG_NODE_REQUEST, req)
                  call wait(req)
                  call irecv(ctx%resources%mpi_comms%world_comm, fragment_idx, 0, MPI_ANY_TAG, req)
                  call wait(req, global_status)

                  if (global_status%MPI_TAG == TAG_NODE_FRAGMENT) then
                     ! Receive fragment type (0 = monomer indices, 1 = intersection atom list)
                     call irecv(ctx%resources%mpi_comms%world_comm, fragment_type, 0, TAG_NODE_FRAGMENT, req)
                     call wait(req)
                     call irecv(ctx%resources%mpi_comms%world_comm, fragment_size, 0, TAG_NODE_FRAGMENT, req)
                     call wait(req)
                     ! Note: must use blocking recv for allocatable arrays since size is unknown
                     allocate (fragment_indices(fragment_size))
                     call recv(ctx%resources%mpi_comms%world_comm, fragment_indices, 0, TAG_NODE_FRAGMENT, global_status)

                     ! Forward to worker
                  call isend(ctx%resources%mpi_comms%node_comm, fragment_idx, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)
                 call isend(ctx%resources%mpi_comms%node_comm, fragment_type, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)
                 call isend(ctx%resources%mpi_comms%node_comm, fragment_size, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)
              call isend(ctx%resources%mpi_comms%node_comm, fragment_indices, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)

                     ! Track which fragment was sent to this worker
                     worker_fragment_map(status%MPI_SOURCE) = fragment_idx

                     deallocate (fragment_indices)
                  else
                     call isend(ctx%resources%mpi_comms%node_comm, -1, status%MPI_SOURCE, TAG_WORKER_FINISH, req)
                     call wait(req)
                     finished_workers = finished_workers + 1
                     more_fragments = .false.
                  end if
               else
                  call isend(ctx%resources%mpi_comms%node_comm, -1, status%MPI_SOURCE, TAG_WORKER_FINISH, req)
                  call wait(req)
                  finished_workers = finished_workers + 1
               end if
            end if
         end if
      end do
   end subroutine node_coordinator_impl

   module subroutine node_worker(ctx)
      !! Node worker for processing fragments assigned by node coordinator
      !! Bond connectivity is accessed via ctx%sys_geom%bonds
      use mqc_many_body_expansion, only: many_body_expansion_t
      class(*), intent(in) :: ctx

      ! Cast to many_body_expansion_t via select type
      select type (ctx)
      class is (many_body_expansion_t)
         call node_worker_impl(ctx)
      class default
         call logger%error("node_worker: expected many_body_expansion_t")
         call abort_comm(comm_world(), 1)
      end select
   end subroutine node_worker

   subroutine node_worker_impl(ctx)
      !! Internal implementation of node_worker with typed context
      use mqc_error, only: error_t
      use mqc_many_body_expansion, only: many_body_expansion_t
      class(many_body_expansion_t), intent(in) :: ctx

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
         call isend(ctx%resources%mpi_comms%node_comm, dummy_msg, 0, TAG_WORKER_REQUEST, req)
         call wait(req)
         call irecv(ctx%resources%mpi_comms%node_comm, fragment_idx, 0, MPI_ANY_TAG, req)
         call wait(req, status)

         select case (status%MPI_TAG)
         case (TAG_WORKER_FRAGMENT)
            ! Receive fragment type (0 = monomer indices, 1 = intersection atom list)
            call irecv(ctx%resources%mpi_comms%node_comm, fragment_type, 0, TAG_WORKER_FRAGMENT, req)
            call wait(req)
            call irecv(ctx%resources%mpi_comms%node_comm, fragment_size, 0, TAG_WORKER_FRAGMENT, req)
            call wait(req)
            ! Note: must use blocking recv for allocatable arrays since size is unknown
            allocate (fragment_indices(fragment_size))
            call recv(ctx%resources%mpi_comms%node_comm, fragment_indices, 0, TAG_WORKER_FRAGMENT, status)

            ! Build physical fragment based on type
            if (ctx%has_geometry()) then
               if (fragment_type == 0) then
                  ! Monomer: fragment_indices are monomer indices
                  call build_fragment_from_indices(ctx%sys_geom, fragment_indices, phys_frag, error, ctx%sys_geom%bonds)
               else
                  ! Intersection: fragment_indices are atom indices
                  call build_fragment_from_atom_list(ctx%sys_geom, fragment_indices, fragment_size, &
                                                     phys_frag, error, ctx%sys_geom%bonds)
               end if

               if (error%has_error()) then
                  call logger%error(error%get_full_trace())
                  call abort_comm(ctx%resources%mpi_comms%world_comm, 1)
               end if

               ! Process the chemistry fragment with physical geometry
               call do_fragment_work(fragment_idx, result, ctx%method_config, phys_frag, ctx%calc_type, &
                                     ctx%resources%mpi_comms%world_comm)

               call phys_frag%destroy()
            else
               ! Process without physical geometry (old behavior)
               call do_fragment_work(fragment_idx, result, ctx%method_config, &
                                     calc_type=ctx%calc_type, world_comm=ctx%resources%mpi_comms%world_comm)
            end if

            ! Send result back to coordinator
            call result_isend(result, ctx%resources%mpi_comms%node_comm, 0, TAG_WORKER_SCALAR_RESULT, req)
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
            call abort_comm(ctx%resources%mpi_comms%world_comm, 1)
         end select
      end do
   end subroutine node_worker_impl

end submodule mpi_fragment_work_smod
