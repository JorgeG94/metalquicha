submodule(mqc_mbe_fragment_distribution_scheme) mpi_fragment_work_smod
   use mqc_error, only: ERROR_VALIDATION, ERROR_GENERIC
   use mqc_work_queue, only: queue_t, queue_init_from_list, queue_pop, queue_is_empty, queue_destroy
   use mqc_group_batching, only: flush_group_results, handle_local_worker_results_to_batch, &
                                 handle_node_results_to_batch, handle_group_results
   use mqc_group_shard_io, only: send_group_assignment_matrix, receive_group_assignment_matrix
   implicit none

contains

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

   subroutine build_fragment_payload_from_row(polymer_row, fragment_type, fragment_size, fragment_indices)
      !! Build fragment payload arrays for a given polymer row.
      !! Caller owns the returned fragment_indices and must deallocate it.
      integer, intent(in) :: polymer_row(:)
      integer(int32), intent(out) :: fragment_type
      integer(int32), intent(out) :: fragment_size
      integer, allocatable, intent(out) :: fragment_indices(:)

      fragment_size = count(polymer_row > 0)
      allocate (fragment_indices(fragment_size))
      fragment_indices = polymer_row(1:fragment_size)

      fragment_type = FRAGMENT_TYPE_MONOMERS
   end subroutine build_fragment_payload_from_row

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

   subroutine send_fragment_payload_from_row(comm, tag, fragment_idx, polymer_row, dest_rank)
      !! Send fragment payload over the specified communicator/tag using a polymer row.
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: tag
      integer(int64), intent(in) :: fragment_idx
      integer, intent(in) :: dest_rank
      integer, intent(in) :: polymer_row(:)
      integer(int32) :: fragment_size, fragment_type
      integer, allocatable :: fragment_indices(:)
      type(request_t) :: req(4)
      integer(int64) :: fragment_idx_int64

      call build_fragment_payload_from_row(polymer_row, fragment_type, fragment_size, fragment_indices)

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
   end subroutine send_fragment_payload_from_row

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

      type :: group_shard_t
         integer(int64), allocatable :: fragment_ids(:)
         integer, allocatable :: polymers(:, :)
      end type group_shard_t

      type(timer_type) :: coord_timer
      integer(int64) :: results_received
      integer :: group_done_count
      integer :: group0_node_count
      integer :: group0_finished_nodes
      integer :: group_id
      integer :: i
      integer :: local_finished_workers
      integer :: group0_done
      integer :: local_node_done
      integer(int32) :: calc_type_local

      ! Storage for results
      type(calculation_result_t), allocatable :: results(:)
      integer(int64) :: worker_fragment_map(ctx%resources%mpi_comms%node_comm%size())
      type(queue_t) :: group0_queue
      integer(int64), allocatable :: group0_fragment_ids(:)
      integer, allocatable :: group0_polymers(:, :)

      integer(int64) :: fragment_idx
      integer(int64) :: chunk_id, chunk_size
      integer(int64), allocatable :: group_counts(:)
      integer(int64), allocatable :: group_fill(:)
      integer, allocatable :: group_leader_by_group(:)
      integer, allocatable :: group_node_counts(:)
      integer :: n_cols
      type(group_shard_t), allocatable :: group_shards(:)

      ! MPI request handles for non-blocking operations
      type(request_t) :: req

      calc_type_local = ctx%calc_type

      results_received = 0_int64
      group_done_count = 0
      group0_finished_nodes = 0
      local_finished_workers = 0
      group0_done = 0
      local_node_done = 0

      ! Allocate storage for results
      allocate (results(ctx%total_fragments))
      worker_fragment_map = 0

      call logger%verbose("Super-global coordinator starting with "//to_char(ctx%total_fragments)// &
                          " fragments for "//to_char(ctx%num_nodes)//" nodes and "// &
                          to_char(ctx%global_groups)//" groups")

      ! Build group leader map and node counts
      allocate (group_leader_by_group(ctx%global_groups))
      group_leader_by_group = -1
      allocate (group_node_counts(ctx%global_groups))
      group_node_counts = 0
      do i = 1, size(ctx%node_leader_ranks)
         group_id = ctx%group_ids(i)
         group_node_counts(group_id) = group_node_counts(group_id) + 1
         if (group_leader_by_group(group_id) == -1) then
            group_leader_by_group(group_id) = ctx%node_leader_ranks(i)
         end if
      end do
      group0_node_count = group_node_counts(1)

      ! Partition fragments into group shards (chunked round-robin)
      allocate (group_counts(ctx%global_groups))
      group_counts = 0_int64
      if (ctx%total_fragments > 0) then
         chunk_size = max(1_int64, ctx%total_fragments/int(ctx%global_groups, int64))
         do fragment_idx = 1_int64, ctx%total_fragments
            chunk_id = (fragment_idx - 1_int64)/chunk_size + 1_int64
            group_id = int(mod(chunk_id - 1_int64, int(ctx%global_groups, int64)) + 1_int64)
            group_counts(group_id) = group_counts(group_id) + 1_int64
         end do
      end if

      allocate (group_shards(ctx%global_groups))
      allocate (group_fill(ctx%global_groups))
      group_fill = 0_int64
      n_cols = size(ctx%polymers, 2)
      do i = 1, ctx%global_groups
         if (group_counts(i) > 0_int64) then
            allocate (group_shards(i)%fragment_ids(group_counts(i)))
            allocate (group_shards(i)%polymers(group_counts(i), n_cols))
         end if
      end do

      if (ctx%total_fragments > 0) then
         do fragment_idx = 1_int64, ctx%total_fragments
            chunk_id = (fragment_idx - 1_int64)/chunk_size + 1_int64
            group_id = int(mod(chunk_id - 1_int64, int(ctx%global_groups, int64)) + 1_int64)
            group_fill(group_id) = group_fill(group_id) + 1_int64
            group_shards(group_id)%fragment_ids(group_fill(group_id)) = fragment_idx
            group_shards(group_id)%polymers(group_fill(group_id), :) = ctx%polymers(fragment_idx, :)
         end do
      end if

      ! Dispatch shards to group globals
      do i = 1, ctx%global_groups
         if (group_leader_by_group(i) == 0) then
            if (allocated(group_shards(i)%fragment_ids)) then
               call move_alloc(group_shards(i)%fragment_ids, group0_fragment_ids)
               call move_alloc(group_shards(i)%polymers, group0_polymers)
            else
               allocate (group0_fragment_ids(0))
               allocate (group0_polymers(0, n_cols))
            end if
         else if (group_leader_by_group(i) > 0) then
            call send_group_assignment_matrix(ctx%resources%mpi_comms%world_comm, group_leader_by_group(i), &
                                              group_shards(i)%fragment_ids, group_shards(i)%polymers)
         end if
         if (allocated(group_shards(i)%fragment_ids)) deallocate (group_shards(i)%fragment_ids)
         if (allocated(group_shards(i)%polymers)) deallocate (group_shards(i)%polymers)
      end do
      deallocate (group_shards)
      deallocate (group_counts)
      deallocate (group_fill)

      ! Initialize local group queue (group 0)
      if (.not. allocated(group0_fragment_ids)) then
         allocate (group0_fragment_ids(0))
         allocate (group0_polymers(0, n_cols))
      end if
      block
         integer(int64), allocatable :: temp_ids(:)
         integer(int64) :: idx

         if (size(group0_fragment_ids) > 0) then
            ! Queue stores local indices (1..N) into group0_fragment_ids/polymers.
            allocate (temp_ids(size(group0_fragment_ids)))
            do idx = 1_int64, size(group0_fragment_ids, kind=int64)
               temp_ids(idx) = idx
            end do
            call queue_init_from_list(group0_queue, temp_ids)
            deallocate (temp_ids)
         else
            group0_queue%count = 0_int64
            group0_queue%head = 1_int64
         end if
      end block

      call coord_timer%start()
      do while (group_done_count < ctx%global_groups)

         ! PRIORITY 1: Receive batched results from group globals
         call handle_group_results(ctx%resources%mpi_comms%world_comm, results, results_received, &
                                   ctx%total_fragments, coord_timer, group_done_count, "fragment")

         ! PRIORITY 2: Check for incoming results from local workers
         if (ctx%resources%mpi_comms%node_comm%size() > 1) then
            call handle_local_worker_results(ctx, worker_fragment_map, results, results_received, coord_timer)
         end if

         ! PRIORITY 3: Check for incoming results from node coordinators (group 0 only)
         call handle_node_results(ctx, results, results_received, coord_timer)

         ! PRIORITY 4: Remote node coordinator requests for group 0
         call handle_group_node_requests(ctx, group0_queue, group0_fragment_ids, group0_polymers, group0_finished_nodes)

         ! PRIORITY 5: Local workers (shared memory) - send new work for group 0
         if (ctx%resources%mpi_comms%node_comm%size() > 1 .and. &
             local_finished_workers < ctx%resources%mpi_comms%node_comm%size() - 1) then
            call handle_local_worker_requests_group(ctx, group0_queue, group0_fragment_ids, group0_polymers, &
                                                    worker_fragment_map, local_finished_workers)
         end if

         ! Mark local node completion once all local workers are finished and queue is empty
         if (local_node_done == 0) then
            if (queue_is_empty(group0_queue) .and. &
                (ctx%resources%mpi_comms%node_comm%size() == 1 .or. &
                 local_finished_workers >= ctx%resources%mpi_comms%node_comm%size() - 1)) then
               local_node_done = 1
               group0_finished_nodes = group0_finished_nodes + 1
            end if
         end if

         if (group0_done == 0) then
            if (group0_finished_nodes >= group0_node_count) then
               group0_done = 1
               group_done_count = group_done_count + 1
            end if
         end if
      end do

      call logger%verbose("Super-global coordinator finished all fragments")
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
      call queue_destroy(group0_queue)
      if (allocated(group0_fragment_ids)) deallocate (group0_fragment_ids)
      if (allocated(group0_polymers)) deallocate (group0_polymers)
      if (allocated(group_leader_by_group)) deallocate (group_leader_by_group)
      if (allocated(group_node_counts)) deallocate (group_node_counts)
      deallocate (results)
   end subroutine global_coordinator_impl

   subroutine handle_node_requests(ctx, fragment_queue, finished_nodes)
      !! Handle a single pending node coordinator request, if any.
      use mqc_many_body_expansion, only: mbe_context_t
      type(mbe_context_t), intent(in) :: ctx
      type(queue_t), intent(inout) :: fragment_queue
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

   subroutine handle_group_node_requests(ctx, fragment_queue, group_fragment_ids, group_polymers, finished_nodes)
      !! Handle a single pending node coordinator request for a group shard, if any.
      use mqc_many_body_expansion, only: many_body_expansion_t
      class(many_body_expansion_t), intent(in) :: ctx
      type(queue_t), intent(inout) :: fragment_queue
      integer(int64), intent(in) :: group_fragment_ids(:)
      integer, intent(in) :: group_polymers(:, :)
      integer, intent(inout) :: finished_nodes

      integer :: request_source, dummy_msg
      integer(int64) :: local_idx, fragment_idx
      type(MPI_Status) :: status
      logical :: has_pending, has_fragment
      type(request_t) :: req

      call iprobe(ctx%resources%mpi_comms%world_comm, MPI_ANY_SOURCE, TAG_NODE_REQUEST, has_pending, status)
      if (.not. has_pending) return

      call irecv(ctx%resources%mpi_comms%world_comm, dummy_msg, status%MPI_SOURCE, TAG_NODE_REQUEST, req)
      call wait(req)
      request_source = status%MPI_SOURCE

      call queue_pop(fragment_queue, local_idx, has_fragment)
      if (has_fragment) then
         fragment_idx = group_fragment_ids(local_idx)
         call send_fragment_payload_from_row(ctx%resources%mpi_comms%world_comm, TAG_NODE_FRAGMENT, fragment_idx, &
                                             group_polymers(local_idx, :), request_source)
      else
         call isend(ctx%resources%mpi_comms%world_comm, -1, request_source, TAG_NODE_FINISH, req)
         call wait(req)
         finished_nodes = finished_nodes + 1
      end if
   end subroutine handle_group_node_requests

   subroutine handle_local_worker_requests_group(ctx, fragment_queue, group_fragment_ids, group_polymers, &
                                                 worker_fragment_map, local_finished_workers)
      !! Handle a single pending local worker request for a group shard, if any.
      use mqc_many_body_expansion, only: many_body_expansion_t
      class(many_body_expansion_t), intent(in) :: ctx
      type(queue_t), intent(inout) :: fragment_queue
      integer(int64), intent(in) :: group_fragment_ids(:)
      integer, intent(in) :: group_polymers(:, :)
      integer(int64), intent(inout) :: worker_fragment_map(:)
      integer, intent(inout) :: local_finished_workers

      integer(int64) :: local_idx, fragment_idx
      integer(int32) :: local_dummy
      type(MPI_Status) :: local_status
      logical :: has_pending, has_fragment
      type(request_t) :: req

      call iprobe(ctx%resources%mpi_comms%node_comm, MPI_ANY_SOURCE, TAG_WORKER_REQUEST, has_pending, local_status)
      if (.not. has_pending) return

      if (worker_fragment_map(local_status%MPI_SOURCE) /= 0) return

      call irecv(ctx%resources%mpi_comms%node_comm, local_dummy, local_status%MPI_SOURCE, TAG_WORKER_REQUEST, req)
      call wait(req)

      call queue_pop(fragment_queue, local_idx, has_fragment)
      if (has_fragment) then
         fragment_idx = group_fragment_ids(local_idx)
         call send_fragment_payload_from_row(ctx%resources%mpi_comms%node_comm, TAG_WORKER_FRAGMENT, fragment_idx, &
                                             group_polymers(local_idx, :), local_status%MPI_SOURCE)
         worker_fragment_map(local_status%MPI_SOURCE) = fragment_idx
      else
         call isend(ctx%resources%mpi_comms%node_comm, -1, local_status%MPI_SOURCE, TAG_WORKER_FINISH, req)
         call wait(req)
         local_finished_workers = local_finished_workers + 1
      end if
   end subroutine handle_local_worker_requests_group

   subroutine handle_local_worker_requests(ctx, fragment_queue, worker_fragment_map, local_finished_workers)
      !! Handle a single pending local worker request, if any.
      use mqc_many_body_expansion, only: mbe_context_t
      type(mbe_context_t), intent(in) :: ctx
      type(queue_t), intent(inout) :: fragment_queue
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

   subroutine group_global_coordinator_impl(ctx)
      !! Group-global coordinator for distributing a fragment shard to node coordinators.
      use mqc_many_body_expansion, only: many_body_expansion_t
      class(many_body_expansion_t), intent(in) :: ctx

      integer(int64), allocatable :: group_fragment_ids(:)
      integer, allocatable :: group_polymers(:, :)
      type(queue_t) :: group_queue
      integer(int64), allocatable :: temp_ids(:)
      integer(int64) :: idx
      integer(int32) :: batch_count
      integer(int64), allocatable :: batch_ids(:)
      type(calculation_result_t), allocatable :: batch_results(:)
      integer :: finished_nodes
      integer :: local_finished_workers
      integer :: group_node_count
      integer :: group_leader_rank, group_id
      integer :: local_node_done
      integer(int64) :: worker_fragment_map(ctx%resources%mpi_comms%node_comm%size())
      type(request_t) :: req

      call get_group_leader_rank(ctx, ctx%resources%mpi_comms%world_comm%rank(), group_leader_rank, group_id)
      if (group_leader_rank /= ctx%resources%mpi_comms%world_comm%rank()) then
         call logger%error("group_global_coordinator_impl called on non-group leader rank")
         call abort_comm(ctx%resources%mpi_comms%world_comm, 1)
      end if
      group_node_count = count(ctx%group_ids == group_id)

      call receive_group_assignment_matrix(ctx%resources%mpi_comms%world_comm, group_fragment_ids, group_polymers)

      if (size(group_fragment_ids) > 0) then
         ! Queue stores local indices (1..N) into group_fragment_ids/group_polymers.
         allocate (temp_ids(size(group_fragment_ids)))
         do idx = 1_int64, size(group_fragment_ids, kind=int64)
            temp_ids(idx) = idx
         end do
         call queue_init_from_list(group_queue, temp_ids)
         deallocate (temp_ids)
      else
         group_queue%count = 0_int64
         group_queue%head = 1_int64
      end if

      batch_count = 0
      allocate (batch_ids(GROUP_RESULT_BATCH_SIZE))
      allocate (batch_results(GROUP_RESULT_BATCH_SIZE))
      finished_nodes = 0
      local_finished_workers = 0
      local_node_done = 0
      worker_fragment_map = 0

      do while (finished_nodes < group_node_count)

         call handle_local_worker_results_to_batch(ctx%resources%mpi_comms%node_comm, &
                                                   ctx%resources%mpi_comms%world_comm, &
                                                   worker_fragment_map, batch_count, batch_ids, batch_results)

         call handle_node_results_to_batch(ctx%resources%mpi_comms%world_comm, batch_count, batch_ids, batch_results)

         call handle_group_node_requests(ctx, group_queue, group_fragment_ids, group_polymers, finished_nodes)

         if (ctx%resources%mpi_comms%node_comm%size() > 1 .and. &
             local_finished_workers < ctx%resources%mpi_comms%node_comm%size() - 1) then
            call handle_local_worker_requests_group(ctx, group_queue, group_fragment_ids, group_polymers, &
                                                    worker_fragment_map, local_finished_workers)
         end if

         if (local_node_done == 0) then
            if (queue_is_empty(group_queue) .and. &
                (ctx%resources%mpi_comms%node_comm%size() == 1 .or. &
                 local_finished_workers >= ctx%resources%mpi_comms%node_comm%size() - 1)) then
               local_node_done = 1
               finished_nodes = finished_nodes + 1
            end if
         end if

         if (batch_count >= GROUP_RESULT_BATCH_SIZE) then
            call flush_group_results(ctx%resources%mpi_comms%world_comm, batch_count, batch_ids, batch_results)
         end if
      end do

      call flush_group_results(ctx%resources%mpi_comms%world_comm, batch_count, batch_ids, batch_results)

      call isend(ctx%resources%mpi_comms%world_comm, 0, 0, TAG_GROUP_DONE, req)
      call wait(req)

      call queue_destroy(group_queue)
      deallocate (group_fragment_ids)
      deallocate (group_polymers)
      deallocate (batch_ids)
      deallocate (batch_results)
   end subroutine group_global_coordinator_impl

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

   subroutine get_group_leader_rank(ctx, node_rank, group_leader_rank, group_id)
      !! Resolve group leader rank and group id for the given node leader rank.
      use mqc_many_body_expansion, only: many_body_expansion_t
      class(many_body_expansion_t), intent(in) :: ctx
      integer, intent(in) :: node_rank
      integer, intent(out) :: group_leader_rank
      integer, intent(out) :: group_id

      integer :: i

      group_leader_rank = 0
      group_id = 1
      if (.not. allocated(ctx%node_leader_ranks)) return
      if (.not. allocated(ctx%group_leader_ranks)) return
      if (.not. allocated(ctx%group_ids)) return

      do i = 1, size(ctx%node_leader_ranks)
         if (ctx%node_leader_ranks(i) == node_rank) then
            group_leader_rank = ctx%group_leader_ranks(i)
            group_id = ctx%group_ids(i)
            return
         end if
      end do
   end subroutine get_group_leader_rank

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

      integer :: group_leader_rank, group_id
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

      call get_group_leader_rank(ctx, ctx%resources%mpi_comms%world_comm%rank(), group_leader_rank, group_id)
      if (group_leader_rank == ctx%resources%mpi_comms%world_comm%rank()) then
         call group_global_coordinator_impl(ctx)
         return
      end if

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
            call isend(ctx%resources%mpi_comms%world_comm, worker_fragment_map(worker_source), &
                       group_leader_rank, TAG_NODE_SCALAR_RESULT, req)  ! fragment_idx
            call wait(req)
call result_isend(worker_result, ctx%resources%mpi_comms%world_comm, group_leader_rank, TAG_NODE_SCALAR_RESULT, req) ! result
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
                  call isend(ctx%resources%mpi_comms%world_comm, dummy_msg, group_leader_rank, TAG_NODE_REQUEST, req)
                  call wait(req)
                  call irecv(ctx%resources%mpi_comms%world_comm, fragment_idx, group_leader_rank, MPI_ANY_TAG, req)
                  call wait(req, global_status)

                  if (global_status%MPI_TAG == TAG_NODE_FRAGMENT) then
                     ! Receive fragment type (0 = monomer indices, 1 = intersection atom list)
                  call irecv(ctx%resources%mpi_comms%world_comm, fragment_type, group_leader_rank, TAG_NODE_FRAGMENT, req)
                     call wait(req)
                  call irecv(ctx%resources%mpi_comms%world_comm, fragment_size, group_leader_rank, TAG_NODE_FRAGMENT, req)
                     call wait(req)
                     ! Note: must use blocking recv for allocatable arrays since size is unknown
                     allocate (fragment_indices(fragment_size))
                     call recv(ctx%resources%mpi_comms%world_comm, fragment_indices, group_leader_rank, &
                               TAG_NODE_FRAGMENT, global_status)

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
