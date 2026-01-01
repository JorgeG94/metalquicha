!! Many-Body Expansion (MBE) calculation module
module mqc_mbe_fragment_distribution_scheme
   !! Implements hierarchical many-body expansion for fragment-based quantum chemistry
   !! calculations with MPI parallelization and energy/gradient computation.
   use pic_types, only: int32, int64, dp
   use pic_timer, only: timer_type
   use pic_blas_interfaces, only: pic_gemm, pic_dot
 use pic_mpi_lib, only: comm_t, send, recv, isend, irecv, wait, iprobe, MPI_Status, request_t, MPI_ANY_SOURCE, MPI_ANY_TAG
   use pic_logger, only: logger => global_logger, verbose_level, info_level
   use pic_io, only: to_char
   use mqc_mbe_io, only: print_fragment_xyz, print_unfragmented_json
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   use mqc_mbe, only: compute_mbe_energy, compute_mbe_energy_gradient
   use mqc_mpi_tags, only: TAG_WORKER_REQUEST, TAG_WORKER_FRAGMENT, TAG_WORKER_FINISH, &
                           TAG_WORKER_SCALAR_RESULT, &
                           TAG_NODE_REQUEST, TAG_NODE_FRAGMENT, TAG_NODE_FINISH, &
                           TAG_NODE_SCALAR_RESULT
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, build_fragment_from_indices, &
                                    build_fragment_from_atom_list, to_angstrom, check_duplicate_atoms
   use mqc_method_types, only: method_type_to_string
   use mqc_calc_types, only: calc_type_to_string, CALC_TYPE_ENERGY, CALC_TYPE_GRADIENT
   use mqc_config_parser, only: bond_t

   ! Method API imports
#ifndef MQC_WITHOUT_TBLITE
   use mqc_method_xtb, only: xtb_method_t
#endif
   use mqc_result_types, only: calculation_result_t, result_send, result_isend, result_recv, result_irecv
   implicit none
   private

   ! Public interface
   public :: do_fragment_work, global_coordinator, node_coordinator
   public :: serial_fragment_processor
   public :: node_worker, unfragmented_calculation

contains

   subroutine do_fragment_work(fragment_idx, result, method, phys_frag, calc_type)
      !! Process a single fragment for quantum chemistry calculation
      !!
      !! Performs energy and gradient calculation on a molecular fragment using
      !! specified quantum chemistry method (GFN-xTB variants).
      !! Verbosity is controlled by the global logger level.

      use pic_logger, only: verbose_level

      integer, intent(in) :: fragment_idx        !! Fragment index for identification
      type(calculation_result_t), intent(out) :: result  !! Computation results
      integer(int32), intent(in) :: method       !! QC method
      type(physical_fragment_t), intent(in), optional :: phys_frag  !! Fragment geometry
      integer(int32), intent(in), optional :: calc_type  !! Calculation type

      integer :: current_log_level  !! Current logger verbosity level
      logical :: is_verbose  !! Whether verbose output is enabled
      integer(int32) :: calc_type_local  !! Local copy of calc_type with default
#ifndef MQC_WITHOUT_TBLITE
      type(xtb_method_t) :: xtb_calc  !! XTB calculator instance
#endif

      ! Set default calc_type if not provided
      if (present(calc_type)) then
         calc_type_local = calc_type
      else
         calc_type_local = CALC_TYPE_ENERGY  ! Default to energy-only calculation
      end if

      ! Query logger to determine verbosity
      call logger%configuration(level=current_log_level)
      is_verbose = (current_log_level >= verbose_level)

      ! Print fragment geometry if provided and verbose mode is enabled
      if (present(phys_frag)) then
         if (is_verbose) then
            call print_fragment_xyz(fragment_idx, phys_frag)
         end if

#ifndef MQC_WITHOUT_TBLITE
         ! Setup XTB method
         xtb_calc%variant = method_type_to_string(method)
         xtb_calc%verbose = is_verbose

         ! Run the calculation using the method API
         select case (calc_type_local)
         case (CALC_TYPE_ENERGY)
            call xtb_calc%calc_energy(phys_frag, result)
         case (CALC_TYPE_GRADIENT)
            call xtb_calc%calc_gradient(phys_frag, result)
         case default
            call logger%error("Unknown calc_type: "//calc_type_to_string(calc_type_local))
            error stop "Invalid calc_type in do_fragment_work"
         end select
#else
         call logger%error("XTB method requested but tblite support not compiled in")
         call logger%error("Please rebuild with -DMQC_ENABLE_TBLITE=ON")
         error stop "tblite support not available"
#endif
      else
         ! For empty fragments, set energy to zero
         call result%energy%reset()
         result%has_energy = .true.
      end if
   end subroutine do_fragment_work

   subroutine global_coordinator(world_comm, node_comm, total_fragments, polymers, max_level, &
                                 node_leader_ranks, num_nodes, sys_geom, calc_type, bonds)
      !! Global coordinator for distributing fragments to node coordinators
      !! will act as a node coordinator for a single node calculation
      !! Uses int64 for total_fragments to handle large fragment counts that overflow int32.
      type(comm_t), intent(in) :: world_comm, node_comm
      integer(int64), intent(in) :: total_fragments
      integer, intent(in) :: max_level, num_nodes
      integer, intent(in) :: polymers(:, :), node_leader_ranks(:)
      type(system_geometry_t), intent(in), optional :: sys_geom
      integer(int32), intent(in), optional :: calc_type
      type(bond_t), intent(in), optional :: bonds(:)

      type(timer_type) :: coord_timer
      integer(int64) :: current_fragment, results_received
      integer :: finished_nodes
      integer :: request_source, dummy_msg, fragment_idx
      type(MPI_Status) :: status, local_status
      logical :: handling_local_workers
      logical :: has_pending
      integer(int32) :: calc_type_local

      ! For local workers
      integer :: local_finished_workers, local_dummy

      ! Storage for results
      type(calculation_result_t), allocatable :: results(:)
      integer(int64) :: worker_fragment_map(node_comm%size())
      integer :: worker_source

      ! MPI request handles for non-blocking operations
      type(request_t) :: req

      ! Set default calc_type if not provided
      if (present(calc_type)) then
         calc_type_local = calc_type
      else
         calc_type_local = CALC_TYPE_ENERGY
      end if

      current_fragment = total_fragments
      finished_nodes = 0
      local_finished_workers = 0
      handling_local_workers = (node_comm%size() > 1)
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
               call iprobe(node_comm, MPI_ANY_SOURCE, TAG_WORKER_SCALAR_RESULT, has_pending, local_status)
               if (.not. has_pending) exit

               worker_source = local_status%MPI_SOURCE

               ! Safety check: worker should have a fragment assigned
               if (worker_fragment_map(worker_source) == 0) then
                  call logger%error("Received result from worker "//to_char(worker_source)// &
                                    " but no fragment was assigned!")
                  error stop "Invalid worker_fragment_map state"
               end if

               ! Receive result and store it using the fragment index for this worker
               call result_irecv(results(worker_fragment_map(worker_source)), node_comm, worker_source, &
                                 TAG_WORKER_SCALAR_RESULT, req)
               call wait(req)
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
            call iprobe(world_comm, MPI_ANY_SOURCE, TAG_NODE_SCALAR_RESULT, has_pending, status)
            if (.not. has_pending) exit

            ! Receive fragment index and result from node coordinator
            ! TODO: serialize the data for better performance
            call irecv(world_comm, fragment_idx, status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT, req)
            call wait(req)
            call result_irecv(results(fragment_idx), world_comm, status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT, req)
            call wait(req)
            results_received = results_received + 1
            if (mod(results_received, max(1_int64, total_fragments/10)) == 0 .or. &
                results_received == total_fragments) then
               call logger%info("  Processed "//to_char(results_received)//"/"// &
                                to_char(total_fragments)//" fragments ["// &
                                to_char(coord_timer%get_elapsed_time())//" s]")
            end if
         end do

         ! PRIORITY 2: Remote node coordinator requests
         call iprobe(world_comm, MPI_ANY_SOURCE, TAG_NODE_REQUEST, has_pending, status)
         if (has_pending) then
            call irecv(world_comm, dummy_msg, status%MPI_SOURCE, TAG_NODE_REQUEST, req)
            call wait(req)
            request_source = status%MPI_SOURCE

            if (current_fragment >= 1) then
               call send_fragment_to_node(world_comm, current_fragment, polymers, request_source)
               current_fragment = current_fragment - 1
            else
               call isend(world_comm, -1, request_source, TAG_NODE_FINISH, req)
               call wait(req)
               finished_nodes = finished_nodes + 1
            end if
         end if

         ! PRIORITY 3: Local workers (shared memory) - send new work
         if (handling_local_workers .and. local_finished_workers < node_comm%size() - 1) then
            call iprobe(node_comm, MPI_ANY_SOURCE, TAG_WORKER_REQUEST, has_pending, local_status)
            if (has_pending) then
               ! Only process work request if this worker doesn't have pending results
               if (worker_fragment_map(local_status%MPI_SOURCE) == 0) then
                  call irecv(node_comm, local_dummy, local_status%MPI_SOURCE, TAG_WORKER_REQUEST, req)
                  call wait(req)

                  if (current_fragment >= 1) then
                     call send_fragment_to_worker(node_comm, current_fragment, polymers, &
                                                  local_status%MPI_SOURCE)
                     ! Track which fragment was sent to this worker
                     worker_fragment_map(local_status%MPI_SOURCE) = current_fragment
                     current_fragment = current_fragment - 1
                  else
                     call isend(node_comm, -1, local_status%MPI_SOURCE, TAG_WORKER_FINISH, req)
                     call wait(req)
                     local_finished_workers = local_finished_workers + 1
                  end if
               end if
               ! If worker still has pending results, skip the work request
               ! It will be processed on the next iteration after results are received
            end if
         end if

         ! Finalize local worker completion
         if (handling_local_workers .and. local_finished_workers >= node_comm%size() - 1 &
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
         real(dp) :: mbe_total_energy
         real(dp), allocatable :: mbe_total_gradient(:, :)

         ! Compute the many-body expansion
         call logger%info(" ")
         call logger%info("Computing Many-Body Expansion (MBE)...")
         call coord_timer%start()

         ! Use combined function if computing gradients (more efficient)
         if (calc_type_local == CALC_TYPE_GRADIENT) then
            if (.not. present(sys_geom)) then
               call logger%error("sys_geom required for gradient calculation in global_coordinator")
               error stop "Missing sys_geom for gradient calculation"
            end if
            allocate (mbe_total_gradient(3, sys_geom%total_atoms))
            call compute_mbe_energy_gradient(polymers, total_fragments, max_level, results, sys_geom, &
                                             mbe_total_energy, mbe_total_gradient, bonds)
            deallocate (mbe_total_gradient)
         else
            call compute_mbe_energy(polymers, total_fragments, max_level, results, mbe_total_energy)
         end if

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
      integer(int32) :: fragment_idx_int32

      fragment_size = count(polymers(fragment_idx, :) > 0)
      allocate (fragment_indices(fragment_size))
      fragment_indices = polymers(fragment_idx, 1:fragment_size)

      ! Standard MBE always uses monomer indices (type 0)
      fragment_type = 0

      ! TODO: serialize the data for better performance
      fragment_idx_int32 = int(fragment_idx, kind=int32)
      call isend(world_comm, fragment_idx_int32, dest_rank, TAG_NODE_FRAGMENT, req(1))
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
      integer(int32) :: fragment_idx_int32

      fragment_size = count(polymers(fragment_idx, :) > 0)
      allocate (fragment_indices(fragment_size))
      fragment_indices = polymers(fragment_idx, 1:fragment_size)

      ! Standard MBE always uses monomer indices (type 0)
      fragment_type = 0

      ! TODO: serialize the data for better performance
      fragment_idx_int32 = int(fragment_idx, kind=int32)
      call isend(node_comm, fragment_idx_int32, dest_rank, TAG_WORKER_FRAGMENT, req(1))
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

   subroutine node_coordinator(world_comm, node_comm, calc_type)
      !! Node coordinator for distributing fragments to local workers
      !! Handles work requests and result collection from local workers
      class(comm_t), intent(in) :: world_comm, node_comm
      integer(int32), intent(in), optional :: calc_type

      integer(int32) :: fragment_idx, fragment_size, fragment_type, dummy_msg
      integer(int32) :: finished_workers
      integer(int32), allocatable :: fragment_indices(:)
      type(MPI_Status) :: status, global_status
      logical :: local_message_pending, more_fragments, has_result
      integer(int32) :: local_dummy

      ! For tracking worker-fragment mapping and collecting results
      integer(int32) :: worker_fragment_map(node_comm%size())
      integer(int32) :: worker_source
      type(calculation_result_t) :: worker_result

      ! MPI request handles for non-blocking operations
      type(request_t) :: req

      finished_workers = 0
      more_fragments = .true.
      dummy_msg = 0
      worker_fragment_map = 0

      do while (finished_workers < node_comm%size() - 1)

         ! PRIORITY 1: Check for incoming results from local workers
         call iprobe(node_comm, MPI_ANY_SOURCE, TAG_WORKER_SCALAR_RESULT, has_result, status)
         if (has_result) then
            worker_source = status%MPI_SOURCE

            ! Safety check: worker should have a fragment assigned
            if (worker_fragment_map(worker_source) == 0) then
               call logger%error("Node coordinator received result from worker "//to_char(worker_source)// &
                                 " but no fragment was assigned!")
               error stop "Invalid worker_fragment_map state in node coordinator"
            end if

            ! Receive result from worker
            call result_irecv(worker_result, node_comm, worker_source, TAG_WORKER_SCALAR_RESULT, req)
            call wait(req)

            ! Forward results to global coordinator with fragment index
            call isend(world_comm, worker_fragment_map(worker_source), 0, TAG_NODE_SCALAR_RESULT, req)  ! fragment_idx
            call wait(req)
            call result_isend(worker_result, world_comm, 0, TAG_NODE_SCALAR_RESULT, req)                ! result
            call wait(req)

            ! Clear the mapping
            worker_fragment_map(worker_source) = 0
         end if

         ! PRIORITY 2: Check for work requests from local workers
         call iprobe(node_comm, MPI_ANY_SOURCE, TAG_WORKER_REQUEST, local_message_pending, status)

         if (local_message_pending) then
            ! Only process work request if this worker doesn't have pending results
            if (worker_fragment_map(status%MPI_SOURCE) == 0) then
               call irecv(node_comm, local_dummy, status%MPI_SOURCE, TAG_WORKER_REQUEST, req)
               call wait(req)

               if (more_fragments) then
                  call isend(world_comm, dummy_msg, 0, TAG_NODE_REQUEST, req)
                  call wait(req)
                  call irecv(world_comm, fragment_idx, 0, MPI_ANY_TAG, req)
                  call wait(req, global_status)

                  if (global_status%MPI_TAG == TAG_NODE_FRAGMENT) then
                     ! Receive fragment type (0 = monomer indices, 1 = intersection atom list)
                     call irecv(world_comm, fragment_type, 0, TAG_NODE_FRAGMENT, req)
                     call wait(req)
                     call irecv(world_comm, fragment_size, 0, TAG_NODE_FRAGMENT, req)
                     call wait(req)
                     ! Note: must use blocking recv for allocatable arrays since size is unknown
                     allocate (fragment_indices(fragment_size))
                     call recv(world_comm, fragment_indices, 0, TAG_NODE_FRAGMENT, global_status)

                     ! Forward to worker
                     call isend(node_comm, fragment_idx, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)
                     call isend(node_comm, fragment_type, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)
                     call isend(node_comm, fragment_size, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)
                     call isend(node_comm, fragment_indices, status%MPI_SOURCE, TAG_WORKER_FRAGMENT, req)
                     call wait(req)

                     ! Track which fragment was sent to this worker
                     worker_fragment_map(status%MPI_SOURCE) = fragment_idx

                     deallocate (fragment_indices)
                  else
                     call isend(node_comm, -1, status%MPI_SOURCE, TAG_WORKER_FINISH, req)
                     call wait(req)
                     finished_workers = finished_workers + 1
                     more_fragments = .false.
                  end if
               else
                  call isend(node_comm, -1, status%MPI_SOURCE, TAG_WORKER_FINISH, req)
                  call wait(req)
                  finished_workers = finished_workers + 1
               end if
            end if
         end if
      end do
   end subroutine node_coordinator

   subroutine node_worker(world_comm, node_comm, sys_geom, method, calc_type, bonds)
      !! Node worker for processing fragments assigned by node coordinator
      class(comm_t), intent(in) :: world_comm, node_comm
      type(system_geometry_t), intent(in), optional :: sys_geom
      integer(int32), intent(in) :: method
      integer(int32), intent(in), optional :: calc_type
      type(bond_t), intent(in), optional :: bonds(:)

      integer(int32) :: fragment_idx, fragment_size, dummy_msg
      integer(int32) :: fragment_type  !! 0 = monomer (indices), 1 = intersection (atom list)
      integer(int32), allocatable :: fragment_indices(:)
      type(calculation_result_t) :: result
      type(MPI_Status) :: status
      type(physical_fragment_t) :: phys_frag

      ! MPI request handles for non-blocking operations
      type(request_t) :: req

      dummy_msg = 0

      do
         call isend(node_comm, dummy_msg, 0, TAG_WORKER_REQUEST, req)
         call wait(req)
         call irecv(node_comm, fragment_idx, 0, MPI_ANY_TAG, req)
         call wait(req, status)

         select case (status%MPI_TAG)
         case (TAG_WORKER_FRAGMENT)
            ! Receive fragment type (0 = monomer indices, 1 = intersection atom list)
            call irecv(node_comm, fragment_type, 0, TAG_WORKER_FRAGMENT, req)
            call wait(req)
            call irecv(node_comm, fragment_size, 0, TAG_WORKER_FRAGMENT, req)
            call wait(req)
            ! Note: must use blocking recv for allocatable arrays since size is unknown
            allocate (fragment_indices(fragment_size))
            call recv(node_comm, fragment_indices, 0, TAG_WORKER_FRAGMENT, status)

            ! Build physical fragment based on type
            if (present(sys_geom)) then
               if (fragment_type == 0) then
                  ! Monomer: fragment_indices are monomer indices
                  call build_fragment_from_indices(sys_geom, fragment_indices, phys_frag, bonds)
               else
                  ! Intersection: fragment_indices are atom indices
                  call build_fragment_from_atom_list(sys_geom, fragment_indices, fragment_size, phys_frag, bonds)
               end if

               ! Process the chemistry fragment with physical geometry
               call do_fragment_work(fragment_idx, result, method, phys_frag, calc_type)

               call phys_frag%destroy()
            else
               ! Process without physical geometry (old behavior)
               call do_fragment_work(fragment_idx, result, method, calc_type=calc_type)
            end if

            ! Send result back to coordinator
            call result_isend(result, node_comm, 0, TAG_WORKER_SCALAR_RESULT, req)
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
            error stop "MPI protocol error in node_worker"
         end select
      end do
   end subroutine node_worker

   subroutine unfragmented_calculation(sys_geom, method, calc_type, bonds)
      !! Run unfragmented calculation on the entire system (nlevel=0)
      !! This is a simple single-process calculation without MPI distribution
      type(system_geometry_t), intent(in), optional :: sys_geom
      integer(int32), intent(in) :: method
      integer(int32), intent(in), optional :: calc_type
      type(bond_t), intent(in), optional :: bonds(:)

      type(calculation_result_t) :: result
      integer :: total_atoms
      type(physical_fragment_t) :: full_system
      integer :: i

      if (.not. present(sys_geom)) then
         call logger%error("sys_geom required for unfragmented calculation")
         error stop "Missing geometry in unfragmented_calculation"
      end if

      total_atoms = sys_geom%total_atoms

      call logger%info("============================================")
      call logger%info("Running unfragmented calculation")
      call logger%info("  Total atoms: "//to_char(total_atoms))
      call logger%info("============================================")

      ! Build the full system as a single fragment
      ! For overlapping fragments, we use the full system directly (not concatenating fragments)
      full_system%n_atoms = total_atoms
      full_system%n_caps = 0
      allocate (full_system%element_numbers(total_atoms))
      allocate (full_system%coordinates(3, total_atoms))

      ! Copy all atoms from system geometry
      full_system%element_numbers = sys_geom%element_numbers
      full_system%coordinates = sys_geom%coordinates

      ! Set charge and multiplicity from system
      full_system%charge = sys_geom%charge
      full_system%multiplicity = sys_geom%multiplicity
      call full_system%compute_nelec()

      ! Validate geometry (check for spatially overlapping atoms)
      call check_duplicate_atoms(full_system)

      ! Process the full system
      call do_fragment_work(0_int32, result, method, phys_frag=full_system, calc_type=calc_type)

      call logger%info("============================================")
      call logger%info("Unfragmented calculation completed")
      block
         character(len=256) :: result_line
         integer :: current_log_level, iatom

         write (result_line, '(a,f25.15)') "  Final energy: ", result%energy%total()
         call logger%info(trim(result_line))

         if (result%has_gradient) then
            write (result_line, '(a,f25.15)') "  Gradient norm: ", sqrt(sum(result%gradient**2))
            call logger%info(trim(result_line))

            ! Print full gradient if verbose and system is small
            call logger%configuration(level=current_log_level)
            if (current_log_level >= verbose_level .and. total_atoms < 100) then
               call logger%info(" ")
               call logger%info("Gradient (Hartree/Bohr):")
               do iatom = 1, total_atoms
                  write (result_line, '(a,i5,a,3f20.12)') "  Atom ", iatom, ": ", &
                     result%gradient(1, iatom), result%gradient(2, iatom), result%gradient(3, iatom)
                  call logger%info(trim(result_line))
               end do
               call logger%info(" ")
            end if
         end if
      end block
      call logger%info("============================================")
      call print_unfragmented_json(result)

      call result%destroy()

   end subroutine unfragmented_calculation

   subroutine serial_fragment_processor(total_fragments, polymers, max_level, sys_geom, method, calc_type, bonds)
      !! Process all fragments serially in single-rank mode
      !! This is used when running with only 1 MPI rank
      integer(int64), intent(in) :: total_fragments
      integer, intent(in) :: polymers(:, :), max_level
      type(system_geometry_t), intent(in) :: sys_geom
      integer(int32), intent(in) :: method
      integer(int32), intent(in), optional :: calc_type
      type(bond_t), intent(in), optional :: bonds(:)

      integer(int64) :: frag_idx
      integer :: fragment_size, current_log_level, iatom
      integer, allocatable :: fragment_indices(:)
      type(calculation_result_t), allocatable :: results(:)
      real(dp) :: mbe_total_energy
      real(dp), allocatable :: mbe_total_gradient(:, :)
      type(physical_fragment_t) :: phys_frag
      type(timer_type) :: coord_timer
      integer(int32) :: calc_type_local

      ! Set default calc_type if not provided
      if (present(calc_type)) then
         calc_type_local = calc_type
      else
         calc_type_local = CALC_TYPE_ENERGY
      end if

      call logger%info("Processing "//to_char(total_fragments)//" fragments serially...")
      call logger%info("  Calculation type: "//calc_type_to_string(calc_type_local))

      allocate (results(total_fragments))

      call omp_set_num_threads(1)
      call coord_timer%start()
      do frag_idx = 1_int64, total_fragments
         fragment_size = count(polymers(frag_idx, :) > 0)
         allocate (fragment_indices(fragment_size))
         fragment_indices = polymers(frag_idx, 1:fragment_size)

         call build_fragment_from_indices(sys_geom, fragment_indices, phys_frag, bonds)

         call do_fragment_work(int(frag_idx), results(frag_idx), method, phys_frag, calc_type=calc_type_local)

         ! Debug output for gradients
         if (calc_type_local == CALC_TYPE_GRADIENT .and. results(frag_idx)%has_gradient) then
            call logger%configuration(level=current_log_level)
            if (current_log_level >= verbose_level) then
               block
                  character(len=512) :: debug_line
                  integer :: iatom_local
                  write (debug_line, '(a,i0,a,*(i0,1x))') "Fragment ", frag_idx, " monomers: ", fragment_indices
                  call logger%verbose(trim(debug_line))
                  write (debug_line, '(a,f25.15)') "  Energy: ", results(frag_idx)%energy%total()
                  call logger%verbose(trim(debug_line))
                  write (debug_line, '(a,f25.15)') "  Gradient norm: ", sqrt(sum(results(frag_idx)%gradient**2))
                  call logger%verbose(trim(debug_line))
                  if (size(results(frag_idx)%gradient, 2) <= 20) then
                     call logger%verbose("  Fragment gradient:")
                     do iatom_local = 1, size(results(frag_idx)%gradient, 2)
                        write (debug_line, '(a,i3,a,3f20.12)') "    Atom ", iatom_local, ": ", &
                           results(frag_idx)%gradient(1, iatom_local), &
                           results(frag_idx)%gradient(2, iatom_local), &
                           results(frag_idx)%gradient(3, iatom_local)
                        call logger%verbose(trim(debug_line))
                     end do
                  end if
               end block
            end if
         end if

         call phys_frag%destroy()
         deallocate (fragment_indices)

         if (mod(frag_idx, max(1_int64, total_fragments/10)) == 0 .or. frag_idx == total_fragments) then
            call logger%info("  Processed "//to_char(frag_idx)//"/"//to_char(total_fragments)// &
                             " fragments ["//to_char(coord_timer%get_elapsed_time())//" s]")
         end if
      end do
      call coord_timer%stop()
      call logger%info("Time to evaluate all fragments "//to_char(coord_timer%get_elapsed_time())//" s")
      call omp_set_num_threads(omp_get_max_threads())

      call logger%info("All fragments processed")

      call logger%info(" ")
      call logger%info("Computing Many-Body Expansion (MBE)...")
      call coord_timer%start()

      ! Use combined function if computing gradients (more efficient)
      if (calc_type_local == CALC_TYPE_GRADIENT) then
         allocate (mbe_total_gradient(3, sys_geom%total_atoms))
         call compute_mbe_energy_gradient(polymers, total_fragments, max_level, results, sys_geom, &
                                          mbe_total_energy, mbe_total_gradient, bonds)
         deallocate (mbe_total_gradient)
      else
         call compute_mbe_energy(polymers, total_fragments, max_level, results, mbe_total_energy)
      end if

      call coord_timer%stop()
      call logger%info("Time to compute MBE "//to_char(coord_timer%get_elapsed_time())//" s")

      deallocate (results)

   end subroutine serial_fragment_processor

end module mqc_mbe_fragment_distribution_scheme
