submodule(mqc_mbe_fragment_distribution_scheme) mpi_fragment_work_smod
   use mqc_error, only: ERROR_VALIDATION, ERROR_GENERIC
   use mqc_work_queue, only: queue_t, queue_init_from_list, queue_pop, queue_is_empty, queue_destroy
   use mqc_group_batching, only: flush_group_results, handle_local_worker_results_to_batch, &
                                 handle_node_results_to_batch, handle_group_results
   use mqc_group_shard_io, only: send_group_assignment_matrix, receive_group_assignment_matrix
   implicit none

contains

   subroutine build_hessian_task_table(polymers, total_fragments, sys_geom, expand, &
                                       task_rows, task_offset, frag_natoms, total_tasks, world_comm)
      !! Expand the screened fragment list into the flat list of independent evaluations
      !!
      !! A Hessian needs 6*n_atoms + 1 gradient evaluations per fragment and none of
      !! them depend on each other. Keeping the whole fragment as one work unit puts a
      !! hard floor of (6*n_max + 1) sequential gradients on the wall time and caps
      !! strong scaling at the fragment count -- which is brutal for a handful of large
      !! fragments. Making every displacement its own schedulable task removes both.
      !!
      !! A task row is the fragment's polymer row plus one trailing column holding the
      !! displacement code:
      !!    0   -> undisplaced reference (supplies the fragment's energy and gradient)
      !!   +k   -> atom (k-1)/3+1 displaced by +h along coordinate mod(k-1,3)+1
      !!   -k   -> the same coordinate displaced by -h
      !!
      !! Packing the code into the row means the existing shard machinery, which ships
      !! an arbitrary integer matrix, carries tasks with no change.
      !!
      !! Tasks are emitted fragment-major so the digest walks one contiguous block per
      !! fragment and can release it as soon as that Hessian is assembled.
      use mqc_error, only: error_t
      integer, intent(in) :: polymers(:, :)
      integer(int64), intent(in) :: total_fragments
      type(system_geometry_t), intent(in) :: sys_geom
      logical, intent(in) :: expand  !! .false. -> one task per fragment (energy/gradient runs)
      integer, allocatable, intent(out) :: task_rows(:, :)      !! (total_tasks, n_cols+1)
      integer(int64), allocatable, intent(out) :: task_offset(:)    !! first task index of each fragment
      integer, allocatable, intent(out) :: frag_natoms(:)       !! capped atom count of each fragment
      integer(int64), intent(out) :: total_tasks
      type(comm_t), intent(in) :: world_comm

      type(physical_fragment_t) :: frag
      type(error_t) :: error
      integer(int64) :: f, base
      integer :: fragment_size, k, n_disp, n_cols

      n_cols = size(polymers, 2)
      allocate (frag_natoms(total_fragments))
      allocate (task_offset(total_fragments))
      frag_natoms = 0

      total_tasks = 0_int64
      do f = 1_int64, total_fragments
         task_offset(f) = total_tasks + 1_int64
         if (expand) then
            ! Atom count must include hydrogen caps, so the fragment has to be built
            fragment_size = count(polymers(f, :) > 0)
            call build_fragment_from_indices(sys_geom, polymers(f, 1:fragment_size), frag, error, sys_geom%bonds)
            if (error%has_error()) then
               call logger%error("build_hessian_task_table: "//error%get_full_trace())
               call abort_comm(world_comm, 1)
            end if
            frag_natoms(f) = frag%n_atoms
            call frag%destroy()
            total_tasks = total_tasks + 1_int64 + 2_int64*3_int64*int(frag_natoms(f), int64)
         else
            total_tasks = total_tasks + 1_int64
         end if
      end do

      allocate (task_rows(total_tasks, n_cols + 1))
      task_rows = 0

      do f = 1_int64, total_fragments
         base = task_offset(f)

         ! Reference point first
         task_rows(base, 1:n_cols) = polymers(f, 1:n_cols)
         task_rows(base, n_cols + 1) = 0
         if (.not. expand) then
            ! Energy/gradient runs keep one task per fragment
            task_rows(base, n_cols + 1) = DISP_WHOLE_FRAGMENT
            cycle
         end if

         ! Then the forward/backward pair for each Cartesian coordinate
         n_disp = 3*frag_natoms(f)
         do k = 1, n_disp
            task_rows(base + 2*k - 1, 1:n_cols) = polymers(f, 1:n_cols)
            task_rows(base + 2*k - 1, n_cols + 1) = k

            task_rows(base + 2*k, 1:n_cols) = polymers(f, 1:n_cols)
            task_rows(base + 2*k, n_cols + 1) = -k
         end do
      end do

   end subroutine build_hessian_task_table

   subroutine build_task_payload_from_row(task_row, fragment_type, fragment_size, fragment_indices, disp_code)
      !! Split a task row into its fragment payload and displacement code.
      !! Caller owns the returned fragment_indices and must deallocate it.
      integer, intent(in) :: task_row(:)
      integer(int32), intent(out) :: fragment_type
      integer(int32), intent(out) :: fragment_size
      integer, allocatable, intent(out) :: fragment_indices(:)
      integer(int32), intent(out) :: disp_code

      integer :: n_cols

      n_cols = size(task_row) - 1
      disp_code = int(task_row(n_cols + 1), int32)

      fragment_size = count(task_row(1:n_cols) > 0)
      allocate (fragment_indices(fragment_size))
      fragment_indices = task_row(1:fragment_size)

      ! Standard MBE always uses monomer indices
      fragment_type = FRAGMENT_TYPE_MONOMERS
   end subroutine build_task_payload_from_row

   subroutine send_task_payload_from_row(comm, tag, task_idx, task_row, dest_rank)
      !! Send a task payload over the specified communicator/tag.
      !! task_idx indexes the flat task list, not the fragment list.
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: tag
      integer(int64), intent(in) :: task_idx
      integer, intent(in) :: dest_rank
      integer, intent(in) :: task_row(:)
      !! task index, fragment type, fragment size, monomer indices, displacement code
      integer, parameter :: N_TASK_PAYLOAD_MSGS = 5
      integer(int32) :: fragment_size, fragment_type, disp_code
      integer, allocatable :: fragment_indices(:)
      type(request_t) :: req(N_TASK_PAYLOAD_MSGS)
      integer(int64) :: task_idx_int64

      call build_task_payload_from_row(task_row, fragment_type, fragment_size, fragment_indices, disp_code)

      task_idx_int64 = int(task_idx, kind=int64)
      call isend(comm, task_idx_int64, dest_rank, tag, req(1))
      call isend(comm, fragment_type, dest_rank, tag, req(2))
      call isend(comm, fragment_size, dest_rank, tag, req(3))
      call isend(comm, fragment_indices, dest_rank, tag, req(4))
      call isend(comm, disp_code, dest_rank, tag, req(5))

      call wait(req(1))
      call wait(req(2))
      call wait(req(3))
      call wait(req(4))
      call wait(req(5))

      deallocate (fragment_indices)
   end subroutine send_task_payload_from_row

   module subroutine do_fragment_work(fragment_idx, result, method_config, phys_frag, calc_type, world_comm, &
                                      print_geometry)
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
      logical, intent(in), optional :: print_geometry  !! Dump the geometry when verbose (default .true.)

      integer :: current_log_level  !! Current logger verbosity level
      logical :: is_verbose  !! Whether verbose output is enabled
      logical :: show_geometry  !! Whether this task should dump its geometry
      integer(int32) :: calc_type_local  !! Local copy of calc_type
      type(method_config_t) :: local_config  !! Local copy for verbose override
      class(qc_method_t), allocatable :: calculator  !! Polymorphic calculator instance

      calc_type_local = calc_type
      show_geometry = .true.
      if (present(print_geometry)) show_geometry = print_geometry

      ! Query logger to determine verbosity
      call logger%configuration(level=current_log_level)
      is_verbose = (current_log_level >= verbose_level)

      ! Print fragment geometry if provided and verbose mode is enabled.
      ! Suppressed for displaced tasks: a flattened Hessian has 6N+1 tasks per
      ! fragment, and dumping a near-identical geometry for every one of them buries
      ! the log without adding anything the reference geometry does not already show.
      if (present(phys_frag)) then
         if (is_verbose .and. show_geometry) then
            call print_fragment_xyz(fragment_idx, phys_frag)
         end if

         ! Copy config and override verbose based on logger level
         local_config = method_config
         local_config%verbose = is_verbose

         ! Create calculator using factory
         ! allocate(..., source=) rather than plain assignment: see the note in
         ! mqc_method_factory -- gfortran 13.2.0 segfaults on intrinsic
         ! assignment from a polymorphic allocatable function result.
         allocate (calculator, source=create_method(local_config))

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

   subroutine assemble_one_fragment_hessian(f, polymers, sys_geom, task_results, &
                                            task_offset, frag_natoms, frag_result, world_comm)
      !! Fold one fragment's displacement results into its Hessian, then release them
      !!
      !! The fragment owns a contiguous block of task results laid out by
      !! build_hessian_task_table as [reference, +1, -1, +2, -2, ...]. Gathering that
      !! block gives exactly the forward/backward gradient arrays the serial path
      !! builds, so the same finite-difference routines produce a bit-identical
      !! Hessian.
      use mqc_error, only: error_t
      use mqc_finite_differences, only: finite_diff_hessian_from_gradients, &
                                        finite_diff_dipole_derivatives, DEFAULT_DISPLACEMENT
      integer(int64), intent(in) :: f
      integer, intent(in) :: polymers(:, :)
      type(system_geometry_t), intent(in) :: sys_geom
      type(calculation_result_t), intent(inout) :: task_results(:)
      integer(int64), intent(in) :: task_offset(:)
      integer, intent(in) :: frag_natoms(:)
      type(calculation_result_t), intent(inout) :: frag_result
      type(comm_t), intent(in) :: world_comm

      type(physical_fragment_t) :: frag
      type(error_t) :: error
      real(dp), allocatable :: fwd_grad(:, :, :), bwd_grad(:, :, :)
      real(dp), allocatable :: fwd_dip(:, :), bwd_dip(:, :)
      integer(int64) :: base
      integer :: fragment_size, k, n_disp, n_atoms
      logical :: have_dipoles

      n_atoms = frag_natoms(f)
      n_disp = 3*n_atoms
      base = task_offset(f)

      fragment_size = count(polymers(f, :) > 0)
      call build_fragment_from_indices(sys_geom, polymers(f, 1:fragment_size), frag, error, sys_geom%bonds)
      if (error%has_error()) then
         call logger%error("assemble_one_fragment_hessian: "//error%get_full_trace())
         call abort_comm(world_comm, 1)
      end if

      allocate (fwd_grad(n_disp, 3, n_atoms), bwd_grad(n_disp, 3, n_atoms))
      allocate (fwd_dip(n_disp, 3), bwd_dip(n_disp, 3))
      fwd_dip = 0.0_dp
      bwd_dip = 0.0_dp
      have_dipoles = .true.

      do k = 1, n_disp
         if (.not. task_results(base + 2*k - 1)%has_gradient .or. &
             .not. task_results(base + 2*k)%has_gradient) then
            call logger%error("Missing gradient for displacement "//to_char(k)// &
                              " of fragment "//to_char(f))
            call abort_comm(world_comm, 1)
         end if
         fwd_grad(k, :, :) = task_results(base + 2*k - 1)%gradient
         bwd_grad(k, :, :) = task_results(base + 2*k)%gradient

         if (task_results(base + 2*k - 1)%has_dipole .and. task_results(base + 2*k)%has_dipole) then
            fwd_dip(k, :) = task_results(base + 2*k - 1)%dipole
            bwd_dip(k, :) = task_results(base + 2*k)%dipole
         else
            have_dipoles = .false.
         end if
      end do

      call finite_diff_hessian_from_gradients(frag, fwd_grad, bwd_grad, DEFAULT_DISPLACEMENT, &
                                              frag_result%hessian)
      frag_result%has_hessian = .true.

      ! The undisplaced point carries the energy, gradient and dipole that belong
      ! beside this fragment's Hessian.
      frag_result%energy = task_results(base)%energy
      frag_result%has_energy = task_results(base)%has_energy
      frag_result%distance = task_results(base)%distance
      if (task_results(base)%has_gradient) then
         frag_result%gradient = task_results(base)%gradient
         frag_result%has_gradient = .true.
      end if
      if (task_results(base)%has_dipole) then
         frag_result%dipole = task_results(base)%dipole
         frag_result%has_dipole = .true.
      end if

      if (have_dipoles) then
         call finite_diff_dipole_derivatives(n_atoms, fwd_dip, bwd_dip, DEFAULT_DISPLACEMENT, &
                                             frag_result%dipole_derivatives)
         frag_result%has_dipole_derivatives = .true.
      end if

      deallocate (fwd_grad, bwd_grad, fwd_dip, bwd_dip)
      call frag%destroy()

      ! Release this fragment's task results now that they are folded in
      do k = 0, 2*n_disp
         call task_results(base + k)%destroy()
      end do

   end subroutine assemble_one_fragment_hessian

   subroutine fold_completed_fragments(polymers, total_fragments, sys_geom, task_results, &
                                       task_offset, frag_natoms, frag_results, &
                                       next_frag, scan_pos, world_comm)
      !! Fold every fragment whose displacement results have all landed
      !!
      !! Called from the coordinator loop so a fragment's Hessian is assembled and its
      !! gradients released while the rest of the queue is still running, instead of
      !! holding every task result until the end. Memory then tracks the number of
      !! fragments in flight rather than the total task count.
      !!
      !! Completion is detected by a monotone scan from the oldest unfolded fragment,
      !! so the bookkeeping is O(1) amortised per task and the drain routines need to
      !! report nothing. A straggler at the head stalls the scan, in which case this
      !! degrades to folding everything at the end -- never worse than that.
      integer, intent(in) :: polymers(:, :)
      integer(int64), intent(in) :: total_fragments
      type(system_geometry_t), intent(in) :: sys_geom
      type(calculation_result_t), intent(inout) :: task_results(:)
      integer(int64), intent(in) :: task_offset(:)
      integer, intent(in) :: frag_natoms(:)
      type(calculation_result_t), intent(inout) :: frag_results(:)
      integer(int64), intent(inout) :: next_frag  !! Oldest fragment not yet folded
      integer(int64), intent(inout) :: scan_pos   !! Next task of next_frag to verify
      type(comm_t), intent(in) :: world_comm

      integer(int64) :: last_task

      do while (next_frag <= total_fragments)
         last_task = task_offset(next_frag) + 2_int64*3_int64*int(frag_natoms(next_frag), int64)

         do while (scan_pos <= last_task)
            if (.not. task_results(scan_pos)%has_gradient) return
            scan_pos = scan_pos + 1_int64
         end do

         call assemble_one_fragment_hessian(next_frag, polymers, sys_geom, task_results, &
                                            task_offset, frag_natoms, frag_results(next_frag), &
                                            world_comm)
         next_frag = next_frag + 1_int64
         if (next_frag <= total_fragments) scan_pos = task_offset(next_frag)
      end do

   end subroutine fold_completed_fragments

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

      ! Storage for results, indexed by flat task index
      type(calculation_result_t), allocatable :: results(:)
      type(calculation_result_t), allocatable :: frag_results(:)

      ! Flat task list: one schedulable evaluation per entry
      integer, allocatable :: task_rows(:, :)
      integer(int64), allocatable :: task_offset(:)
      integer, allocatable :: frag_natoms(:)
      integer(int64) :: total_tasks
      integer(int64) :: next_frag, scan_pos
      logical :: expand_displacements
      integer(int64) :: worker_fragment_map(ctx%resources%mpi_comms%node_comm%size())
      type(queue_t) :: group0_queue
      integer(int64), allocatable :: group0_fragment_ids(:)
      integer, allocatable :: group0_polymers(:, :)

      integer(int64) :: task_idx
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

      ! A Hessian fans each fragment out into its 6N+1 independent gradient
      ! evaluations; energy and gradient runs stay one task per fragment.
      expand_displacements = (calc_type_local == CALC_TYPE_HESSIAN)
      call build_hessian_task_table(ctx%polymers, ctx%total_fragments, ctx%sys_geom, expand_displacements, &
                                    task_rows, task_offset, frag_natoms, total_tasks, &
                                    ctx%resources%mpi_comms%world_comm)

      ! Allocate storage for results
      allocate (results(total_tasks))
      worker_fragment_map = 0

      ! Fragment Hessians are folded in as their displacements land, so the targets
      ! have to exist before the coordinator loop starts.
      next_frag = 1_int64
      scan_pos = 1_int64
      if (expand_displacements) allocate (frag_results(ctx%total_fragments))

      if (expand_displacements) then
         call logger%info("Hessian displacements flattened into the fragment queue:")
         call logger%info("  Fragments:            "//to_char(ctx%total_fragments))
         call logger%info("  Schedulable tasks:    "//to_char(total_tasks))
         call logger%info("  Max useful ranks:     "//to_char(total_tasks)//" (was "// &
                          to_char(ctx%total_fragments)//")")
      end if

      call logger%verbose("Super-global coordinator starting with "//to_char(total_tasks)// &
                          " tasks for "//to_char(ctx%num_nodes)//" nodes and "// &
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

      ! Partition tasks into group shards (chunked round-robin)
      allocate (group_counts(ctx%global_groups))
      group_counts = 0_int64
      if (total_tasks > 0) then
         chunk_size = max(1_int64, total_tasks/int(ctx%global_groups, int64))
         do task_idx = 1_int64, total_tasks
            chunk_id = (task_idx - 1_int64)/chunk_size + 1_int64
            group_id = int(mod(chunk_id - 1_int64, int(ctx%global_groups, int64)) + 1_int64)
            group_counts(group_id) = group_counts(group_id) + 1_int64
         end do
      end if

      allocate (group_shards(ctx%global_groups))
      allocate (group_fill(ctx%global_groups))
      group_fill = 0_int64
      n_cols = size(task_rows, 2)
      do i = 1, ctx%global_groups
         if (group_counts(i) > 0_int64) then
            allocate (group_shards(i)%fragment_ids(group_counts(i)))
            allocate (group_shards(i)%polymers(group_counts(i), n_cols))
         end if
      end do

      if (total_tasks > 0) then
         do task_idx = 1_int64, total_tasks
            chunk_id = (task_idx - 1_int64)/chunk_size + 1_int64
            group_id = int(mod(chunk_id - 1_int64, int(ctx%global_groups, int64)) + 1_int64)
            group_fill(group_id) = group_fill(group_id) + 1_int64
            group_shards(group_id)%fragment_ids(group_fill(group_id)) = task_idx
            group_shards(group_id)%polymers(group_fill(group_id), :) = task_rows(task_idx, :)
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
      do while (group_done_count < ctx%global_groups .or. results_received < total_tasks)

         ! PRIORITY 1: Receive batched results from group globals
         call handle_group_results(ctx%resources%mpi_comms%world_comm, results, results_received, &
                                   total_tasks, coord_timer, group_done_count, "task")

         ! PRIORITY 2: Check for incoming results from local workers
         if (ctx%resources%mpi_comms%node_comm%size() > 1) then
            call handle_local_worker_results(ctx, worker_fragment_map, results, results_received, &
                                             total_tasks, coord_timer)
         end if

         ! PRIORITY 3: Check for incoming results from node coordinators (group 0 only)
         call handle_node_results(ctx, results, results_received, total_tasks, coord_timer)

         ! Fold any fragment whose displacements have all landed, releasing its
         ! task gradients while the rest of the queue is still running.
         if (expand_displacements) then
            call fold_completed_fragments(ctx%polymers, ctx%total_fragments, ctx%sys_geom, results, &
                                          task_offset, frag_natoms, frag_results, next_frag, scan_pos, &
                                          ctx%resources%mpi_comms%world_comm)
         end if

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

      call logger%verbose("Super-global coordinator finished all tasks")
      call coord_timer%stop()
      call logger%info("Time to evaluate all tasks "//to_char(coord_timer%get_elapsed_time())//" s")
      block
         use mqc_result_types, only: mbe_result_t
         type(mbe_result_t) :: mbe_result

         ! Most fragments were folded as their displacements landed; this only
         ! mops up whatever a straggler at the head of the scan held back.
         if (expand_displacements) then
            call coord_timer%start()
            if (next_frag <= ctx%total_fragments) then
               call logger%info("Assembling remaining fragment Hessians ("// &
                                to_char(ctx%total_fragments - next_frag + 1_int64)//" of "// &
                                to_char(ctx%total_fragments)//" deferred)...")
               call fold_completed_fragments(ctx%polymers, ctx%total_fragments, ctx%sys_geom, results, &
                                             task_offset, frag_natoms, frag_results, next_frag, scan_pos, &
                                             ctx%resources%mpi_comms%world_comm)
            end if
            if (next_frag <= ctx%total_fragments) then
               call logger%error("Fragment "//to_char(next_frag)//" never received all its displacements")
               call abort_comm(ctx%resources%mpi_comms%world_comm, 1)
            end if
            call move_alloc(frag_results, results)
            call coord_timer%stop()
            call logger%info("Fragment Hessians finalised in "//to_char(coord_timer%get_elapsed_time())//" s")
         end if

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
      if (allocated(task_rows)) deallocate (task_rows)
      if (allocated(task_offset)) deallocate (task_offset)
      if (allocated(frag_natoms)) deallocate (frag_natoms)
      deallocate (results)
   end subroutine global_coordinator_impl

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
         call send_task_payload_from_row(ctx%resources%mpi_comms%world_comm, TAG_NODE_FRAGMENT, fragment_idx, &
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
         call send_task_payload_from_row(ctx%resources%mpi_comms%node_comm, TAG_WORKER_FRAGMENT, fragment_idx, &
                                         group_polymers(local_idx, :), local_status%MPI_SOURCE)
         worker_fragment_map(local_status%MPI_SOURCE) = fragment_idx
      else
         call isend(ctx%resources%mpi_comms%node_comm, -1, local_status%MPI_SOURCE, TAG_WORKER_FINISH, req)
         call wait(req)
         local_finished_workers = local_finished_workers + 1
      end if
   end subroutine handle_local_worker_requests_group

   subroutine handle_local_worker_results(ctx, worker_fragment_map, results, results_received, total_items, coord_timer)
      !! Drain results from local workers and update tracking state.
      use mqc_many_body_expansion, only: mbe_context_t
      type(mbe_context_t), intent(in) :: ctx
      integer(int64), intent(inout) :: worker_fragment_map(:)
      type(calculation_result_t), intent(inout) :: results(:)
      integer(int64), intent(inout) :: results_received
      integer(int64), intent(in) :: total_items
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
         if (mod(results_received, max(1_int64, total_items/10)) == 0 .or. &
             results_received == total_items) then
            call logger%info("  Processed "//to_char(results_received)//"/"// &
                             to_char(total_items)//" tasks ["// &
                             to_char(coord_timer%get_elapsed_time())//" s]")
         end if
      end do
   end subroutine handle_local_worker_results

   subroutine handle_node_results(ctx, results, results_received, total_items, coord_timer)
      !! Drain results from remote node coordinators and update tracking state.
      use mqc_many_body_expansion, only: mbe_context_t
      type(mbe_context_t), intent(in) :: ctx
      type(calculation_result_t), intent(inout) :: results(:)
      integer(int64), intent(inout) :: results_received
      integer(int64), intent(in) :: total_items
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
         call result_irecv(results(fragment_idx), ctx%resources%mpi_comms%world_comm, status%MPI_SOURCE, &
                           TAG_NODE_SCALAR_RESULT, req)
         call wait(req)

         if (results(fragment_idx)%has_error) then
            call logger%error("Fragment "//to_char(fragment_idx)//" calculation failed: "// &
                              results(fragment_idx)%error%get_message())
            call abort_comm(ctx%resources%mpi_comms%world_comm, 1)
         end if

         results_received = results_received + 1
         if (mod(results_received, max(1_int64, total_items/10)) == 0 .or. &
             results_received == total_items) then
            call logger%info("  Processed "//to_char(results_received)//"/"// &
                             to_char(total_items)//" tasks ["// &
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
      integer(int64) :: results_received
      integer(int64) :: total_group_fragments
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
      results_received = 0_int64
      total_group_fragments = int(size(group_fragment_ids, kind=int64), int64)
      finished_nodes = 0
      local_finished_workers = 0
      local_node_done = 0
      worker_fragment_map = 0

      do while (finished_nodes < group_node_count .or. results_received < total_group_fragments)

         call handle_local_worker_results_to_batch(ctx%resources%mpi_comms%node_comm, &
                                                   ctx%resources%mpi_comms%world_comm, &
                                                   worker_fragment_map, batch_count, batch_ids, batch_results, &
                                                   results_received)

         call handle_node_results_to_batch(ctx%resources%mpi_comms%world_comm, batch_count, batch_ids, batch_results, &
                                           results_received)

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
      integer(int32) :: fragment_size, fragment_type, dummy_msg, disp_code
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
call result_isend(worker_result, ctx%resources%mpi_comms%world_comm, group_leader_rank, TAG_NODE_SCALAR_RESULT, req)  ! result
            call wait(req)

            ! Release before the next receive. pic-mpi's array recv only
            ! allocates when the target is unallocated -- it does not resize a
            ! buffer that is already there, but still receives the sender's full
            ! count into it. Reusing this one across fragments of different sizes
            ! therefore writes past the end of the smaller allocation and
            ! corrupts the heap, surfacing later as an abort in an unrelated
            ! deallocate.
            call worker_result%destroy()

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
                     call irecv(ctx%resources%mpi_comms%world_comm, disp_code, group_leader_rank, TAG_NODE_FRAGMENT, req)
                     call wait(req)

                     ! Forward to worker
                  call isend(ctx%resources%mpi_comms%node_comm, fragment_idx, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)
                 call isend(ctx%resources%mpi_comms%node_comm, fragment_type, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)
                 call isend(ctx%resources%mpi_comms%node_comm, fragment_size, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)
              call isend(ctx%resources%mpi_comms%node_comm, fragment_indices, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)
                     call isend(ctx%resources%mpi_comms%node_comm, disp_code, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
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
      use mqc_finite_differences, only: copy_and_displace_geometry, DEFAULT_DISPLACEMENT
      class(many_body_expansion_t), intent(in) :: ctx

      integer(int64) :: fragment_idx
      integer(int32) :: fragment_size, dummy_msg
      integer(int32) :: fragment_type  !! 0 = monomer (indices), 1 = intersection (atom list)
      integer(int32) :: disp_code      !! 0 = reference point, +/-k = displaced coordinate k
      integer(int32), allocatable :: fragment_indices(:)
      type(calculation_result_t) :: result
      type(MPI_Status) :: status
      type(physical_fragment_t) :: phys_frag, disp_frag
      type(error_t) :: error
      integer :: atom_idx, coord_idx
      real(dp) :: step

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
            call irecv(ctx%resources%mpi_comms%node_comm, disp_code, 0, TAG_WORKER_FRAGMENT, req)
            call wait(req)

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

               if (disp_code == DISP_WHOLE_FRAGMENT) then
                  ! One task per fragment: run the requested calculation as it stands.
                  ! This covers every energy/gradient run and the whole GMBE path.
                  call do_fragment_work(fragment_idx, result, ctx%method_config, phys_frag, ctx%calc_type, &
                                        ctx%resources%mpi_comms%world_comm)
               else if (disp_code == 0) then
                  ! Flattened Hessian, undisplaced reference point: supplies the
                  ! fragment's energy and gradient alongside its assembled Hessian.
                  call do_fragment_work(fragment_idx, result, ctx%method_config, phys_frag, &
                                        CALC_TYPE_GRADIENT, ctx%resources%mpi_comms%world_comm)
               else
                  ! Flattened Hessian, one displaced coordinate. The Hessian itself is
                  ! assembled by the coordinator once all of a fragment's tasks land.
                  atom_idx = (abs(disp_code) - 1)/3 + 1
                  coord_idx = mod(abs(disp_code) - 1, 3) + 1
                  step = sign(DEFAULT_DISPLACEMENT, real(disp_code, dp))
                  call copy_and_displace_geometry(phys_frag, atom_idx, coord_idx, step, disp_frag)
                  call do_fragment_work(fragment_idx, result, ctx%method_config, disp_frag, &
                                        CALC_TYPE_GRADIENT, ctx%resources%mpi_comms%world_comm, &
                                        print_geometry=.false.)
                  call disp_frag%destroy()
               end if

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
