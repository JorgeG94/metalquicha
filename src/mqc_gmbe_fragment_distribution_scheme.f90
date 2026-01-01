!! Generalized Many-Body Expansion (GMBE) fragment distribution module
module mqc_gmbe_fragment_distribution_scheme
   !! Implements fragment distribution schemes for GMBE calculations with overlapping fragments
   !! Handles both serial and MPI-parallelized distribution of monomers and intersection fragments
   use pic_types, only: int32, int64, dp
   use pic_timer, only: timer_type
 use pic_mpi_lib, only: comm_t, send, recv, isend, irecv, wait, iprobe, MPI_Status, request_t, MPI_ANY_SOURCE, MPI_ANY_TAG
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_mpi_tags, only: TAG_WORKER_REQUEST, TAG_WORKER_FRAGMENT, TAG_WORKER_FINISH, &
                           TAG_WORKER_SCALAR_RESULT, &
                           TAG_NODE_REQUEST, TAG_NODE_FRAGMENT, TAG_NODE_FINISH, &
                           TAG_NODE_SCALAR_RESULT
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, build_fragment_from_indices, &
                                    build_fragment_from_atom_list
   use mqc_config_parser, only: bond_t
   use mqc_result_types, only: calculation_result_t, result_send, result_isend, result_recv, result_irecv
   use mqc_mbe_fragment_distribution_scheme, only: do_fragment_work
   use mqc_mbe_io, only: print_gmbe_json, print_gmbe_pie_json
   implicit none
   private

   ! Public interface
   public :: serial_gmbe_processor, gmbe_coordinator
   public :: serial_gmbe_pie_processor  !! New PIE-based processor

contains

   subroutine serial_gmbe_processor(n_monomers, polymers, intersections, intersection_sets, intersection_levels, &
                                    n_intersections, sys_geom, method, calc_type, bonds)
      !! Serial GMBE processor: builds all fragments (monomers + intersections) and computes GMBE energy
      use mqc_mbe, only: compute_gmbe_energy

      integer, intent(in) :: n_monomers  !! Number of monomers
      integer, intent(in) :: polymers(:, :)  !! Monomer indices (n_monomers, 1)
      integer, intent(in) :: intersections(:, :)  !! Intersection atom lists
      integer, intent(in) :: intersection_sets(:, :)  !! k-tuples that created intersections
      integer, intent(in) :: intersection_levels(:)  !! Level k of each intersection
      integer, intent(in) :: n_intersections  !! Number of intersection fragments
      type(system_geometry_t), intent(in) :: sys_geom
      integer(int32), intent(in) :: method, calc_type
      type(bond_t), intent(in), optional :: bonds(:)

      type(calculation_result_t), allocatable :: all_results(:)
      type(physical_fragment_t) :: phys_frag
      integer :: i, intersection_size
      integer, allocatable :: monomer_idx(:), monomer_indices(:)
      real(dp) :: total_energy

      ! Total results = monomers + intersections
      allocate (all_results(n_monomers + n_intersections))

      ! Build and calculate monomers
      call logger%info("Processing "//to_char(n_monomers)//" monomer fragments...")
      do i = 1, n_monomers
         allocate (monomer_idx(1))
         monomer_idx(1) = polymers(i, 1)
         ! Use build_fragment_from_indices for monomers to preserve fragment charges/multiplicities
         call build_fragment_from_indices(sys_geom, monomer_idx, phys_frag, bonds)
         call do_fragment_work(i, all_results(i), method, phys_frag, calc_type)
         call phys_frag%destroy()
         deallocate (monomer_idx)
      end do

      ! Build and calculate intersections
      if (n_intersections > 0) then
         call logger%info("Processing "//to_char(n_intersections)//" intersection fragments...")
         do i = 1, n_intersections
            ! Find size of this intersection (first zero marks end)
            intersection_size = 0
            do while (intersection_size < size(intersections, 1))
               if (intersections(intersection_size + 1, i) == 0) exit
               intersection_size = intersection_size + 1
            end do

            ! Build intersection fragment from atom list
            call build_fragment_from_atom_list(sys_geom, intersections(1:intersection_size, i), &
                                               intersection_size, phys_frag, bonds)
            call do_fragment_work(n_monomers + i, all_results(n_monomers + i), method, phys_frag, calc_type)
            call phys_frag%destroy()
         end do
      end if

      ! Compute GMBE energy
      allocate (monomer_indices(n_monomers))
      do i = 1, n_monomers
         monomer_indices(i) = polymers(i, 1)
      end do

      call compute_gmbe_energy(monomer_indices, n_monomers, all_results(1:n_monomers), &
                               n_intersections, all_results(n_monomers + 1:n_monomers + n_intersections), &
                               intersection_sets, intersection_levels, total_energy)

      call logger%info(" ")
      call logger%info("GMBE calculation completed successfully")
      call logger%info("Final GMBE energy: "//to_char(total_energy)//" Hartree")
      call logger%info(" ")

      ! Write JSON output
      if (n_intersections > 0) then
         call print_gmbe_json(n_monomers, monomer_indices, all_results(1:n_monomers), &
                              n_intersections, all_results(n_monomers + 1:n_monomers + n_intersections), &
                              intersection_sets, intersection_levels, total_energy)
      else
         ! No intersections - omit optional parameters
         call print_gmbe_json(n_monomers, monomer_indices, all_results(1:n_monomers), &
                              0, total_energy=total_energy)
      end if

      deallocate (all_results, monomer_indices)

   end subroutine serial_gmbe_processor

   subroutine gmbe_coordinator(world_comm, node_comm, n_monomers, polymers, intersections, intersection_sets, &
                   intersection_levels, n_intersections, node_leader_ranks, num_nodes, sys_geom, method, calc_type, bonds)
      !! MPI-enabled GMBE coordinator: distributes monomers and intersections across ranks
      !! Similar to global_coordinator but handles GMBE-specific fragment types
      use mqc_mbe, only: compute_gmbe_energy

      type(comm_t), intent(in) :: world_comm, node_comm
      integer, intent(in) :: n_monomers  !! Number of monomers
      integer, intent(in) :: polymers(:, :)  !! Monomer indices (n_monomers, 1)
      integer, intent(in) :: intersections(:, :)  !! Intersection atom lists
      integer, intent(in) :: intersection_sets(:, :)  !! k-tuples that created intersections
      integer, intent(in) :: intersection_levels(:)  !! Level k of each intersection
      integer, intent(in) :: n_intersections  !! Number of intersection fragments
      integer, intent(in) :: node_leader_ranks(:), num_nodes
      type(system_geometry_t), intent(in) :: sys_geom
      integer(int32), intent(in) :: method, calc_type
      type(bond_t), intent(in), optional :: bonds(:)

      type(timer_type) :: coord_timer
      integer :: current_fragment_idx, total_fragments, results_received
      integer :: finished_nodes
      integer :: request_source, dummy_msg, fragment_idx
      type(MPI_Status) :: status, local_status
      logical :: handling_local_workers
      logical :: has_pending
      logical :: is_monomer
      integer :: intersection_size

      ! For local workers
      integer :: local_finished_workers, local_dummy

      ! Storage for results
      type(calculation_result_t), allocatable :: monomer_results(:), intersection_results(:)
      integer :: worker_fragment_map(node_comm%size())
      integer :: worker_source
      integer, allocatable :: monomer_indices(:)
      real(dp) :: total_energy

      ! MPI request handles for non-blocking operations
      type(request_t) :: req

      ! Total fragments = monomers + intersections
      total_fragments = n_monomers + n_intersections
      current_fragment_idx = total_fragments
      finished_nodes = 0
      local_finished_workers = 0
      handling_local_workers = (node_comm%size() > 1)
      results_received = 0
      worker_fragment_map = 0

      ! Allocate storage for results
      allocate (monomer_results(n_monomers))
      if (n_intersections > 0) then
         allocate (intersection_results(n_intersections))
      end if

      call logger%verbose("GMBE coordinator starting with "//to_char(n_monomers)// &
                          " monomers and "//to_char(n_intersections)//" intersections for "// &
                          to_char(num_nodes)//" nodes")

      call coord_timer%start()
      do while (finished_nodes < num_nodes)

         ! PRIORITY 1: Check for incoming results from local workers
         if (handling_local_workers) then
            do
               call iprobe(node_comm, MPI_ANY_SOURCE, TAG_WORKER_SCALAR_RESULT, has_pending, local_status)
               if (.not. has_pending) exit

               worker_source = local_status%MPI_SOURCE

               if (worker_fragment_map(worker_source) == 0) then
                  call logger%error("Received result from worker "//to_char(worker_source)// &
                                    " but no fragment was assigned!")
                  error stop "Invalid worker_fragment_map state"
               end if

               ! Determine if this was a monomer or intersection
               fragment_idx = worker_fragment_map(worker_source)
               is_monomer = (fragment_idx <= n_monomers)

               ! Receive result and store it
               if (is_monomer) then
                  call result_irecv(monomer_results(fragment_idx), node_comm, worker_source, &
                                    TAG_WORKER_SCALAR_RESULT, req)
               else
                  call result_irecv(intersection_results(fragment_idx - n_monomers), node_comm, worker_source, &
                                    TAG_WORKER_SCALAR_RESULT, req)
               end if
               call wait(req)
               worker_fragment_map(worker_source) = 0
               results_received = results_received + 1
               if (mod(results_received, max(1, total_fragments/10)) == 0 .or. &
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

            ! Receive fragment index
            call irecv(world_comm, fragment_idx, status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT, req)
            call wait(req)

            ! Determine if this was a monomer or intersection
            is_monomer = (fragment_idx <= n_monomers)

            ! Receive result
            if (is_monomer) then
               call result_irecv(monomer_results(fragment_idx), world_comm, status%MPI_SOURCE, &
                                 TAG_NODE_SCALAR_RESULT, req)
            else
               call result_irecv(intersection_results(fragment_idx - n_monomers), world_comm, status%MPI_SOURCE, &
                                 TAG_NODE_SCALAR_RESULT, req)
            end if
            call wait(req)
            results_received = results_received + 1
            if (mod(results_received, max(1, total_fragments/10)) == 0 .or. &
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

            if (current_fragment_idx >= 1) then
               call send_gmbe_fragment_to_node(world_comm, current_fragment_idx, n_monomers, polymers, &
                                               intersections, request_source)
               current_fragment_idx = current_fragment_idx - 1
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
               if (worker_fragment_map(local_status%MPI_SOURCE) == 0) then
                  call irecv(node_comm, local_dummy, local_status%MPI_SOURCE, TAG_WORKER_REQUEST, req)
                  call wait(req)

                  if (current_fragment_idx >= 1) then
                     call send_gmbe_fragment_to_worker(node_comm, current_fragment_idx, n_monomers, polymers, &
                                                       intersections, local_status%MPI_SOURCE)
                     worker_fragment_map(local_status%MPI_SOURCE) = current_fragment_idx
                     current_fragment_idx = current_fragment_idx - 1
                  else
                     call isend(node_comm, -1, local_status%MPI_SOURCE, TAG_WORKER_FINISH, req)
                     call wait(req)
                     local_finished_workers = local_finished_workers + 1
                  end if
               end if
            end if
         end if

         ! Finalize local worker completion
         if (handling_local_workers .and. local_finished_workers >= node_comm%size() - 1 &
             .and. results_received >= total_fragments) then
            handling_local_workers = .false.
            finished_nodes = finished_nodes + 1
         end if
      end do

      call logger%verbose("GMBE coordinator finished all fragments")
      call coord_timer%stop()
      call logger%info("Time to evaluate all fragments "//to_char(coord_timer%get_elapsed_time())//" s")

      ! Compute GMBE energy
      call logger%info(" ")
      call logger%info("Computing Generalized Many-Body Expansion (GMBE)...")
      call coord_timer%start()

      allocate (monomer_indices(n_monomers))
      do fragment_idx = 1, n_monomers
         monomer_indices(fragment_idx) = polymers(fragment_idx, 1)
      end do

      if (n_intersections > 0) then
         call compute_gmbe_energy(monomer_indices, n_monomers, monomer_results, &
                                  n_intersections, intersection_results, &
                                  intersection_sets, intersection_levels, total_energy)
      else
         ! No intersections - omit optional parameters
         call compute_gmbe_energy(monomer_indices, n_monomers, monomer_results, &
                                  0, total_energy=total_energy)
      end if

      call coord_timer%stop()
      call logger%info("Time to compute GMBE "//to_char(coord_timer%get_elapsed_time())//" s")
      call logger%info(" ")
      call logger%info("GMBE calculation completed successfully")
      call logger%info("Final GMBE energy: "//to_char(total_energy)//" Hartree")
      call logger%info(" ")

      ! Write JSON output
      if (n_intersections > 0) then
         call print_gmbe_json(n_monomers, monomer_indices, monomer_results, &
                              n_intersections, intersection_results, &
                              intersection_sets, intersection_levels, total_energy)
      else
         ! No intersections - omit optional parameters
         call print_gmbe_json(n_monomers, monomer_indices, monomer_results, &
                              0, total_energy=total_energy)
      end if

      ! Cleanup
      deallocate (monomer_results, monomer_indices)
      if (allocated(intersection_results)) deallocate (intersection_results)

   end subroutine gmbe_coordinator

   subroutine send_gmbe_fragment_to_node(world_comm, fragment_idx, n_monomers, polymers, intersections, dest_rank)
      !! Send GMBE fragment (monomer or intersection) to remote node coordinator
      type(comm_t), intent(in) :: world_comm
      integer, intent(in) :: fragment_idx, n_monomers
      integer, intent(in) :: polymers(:, :)
      integer, intent(in) :: intersections(:, :)
      integer, intent(in) :: dest_rank

      integer :: fragment_size, intersection_size, i
      integer, allocatable :: fragment_indices(:)
      logical :: is_monomer
      type(request_t) :: req(3)

      is_monomer = (fragment_idx <= n_monomers)

      if (is_monomer) then
         ! Monomer fragment
         fragment_size = 1
         allocate (fragment_indices(fragment_size))
         fragment_indices(1) = polymers(fragment_idx, 1)
      else
         ! Intersection fragment
         intersection_size = 0
         do while (intersection_size < size(intersections, 1))
            if (intersections(intersection_size + 1, fragment_idx - n_monomers) == 0) exit
            intersection_size = intersection_size + 1
         end do
         fragment_size = intersection_size
         allocate (fragment_indices(fragment_size))
         fragment_indices = intersections(1:intersection_size, fragment_idx - n_monomers)
      end if

      call isend(world_comm, fragment_idx, dest_rank, TAG_NODE_FRAGMENT, req(1))
      call isend(world_comm, fragment_size, dest_rank, TAG_NODE_FRAGMENT, req(2))
      call isend(world_comm, fragment_indices, dest_rank, TAG_NODE_FRAGMENT, req(3))

      call wait(req(1))
      call wait(req(2))
      call wait(req(3))

      deallocate (fragment_indices)
   end subroutine send_gmbe_fragment_to_node

   subroutine send_gmbe_fragment_to_worker(node_comm, fragment_idx, n_monomers, polymers, intersections, dest_rank)
      !! Send GMBE fragment (monomer or intersection) to local worker
      type(comm_t), intent(in) :: node_comm
      integer, intent(in) :: fragment_idx, n_monomers
      integer, intent(in) :: polymers(:, :)
      integer, intent(in) :: intersections(:, :)
      integer, intent(in) :: dest_rank

      integer :: fragment_size, intersection_size
      integer, allocatable :: fragment_indices(:)
      logical :: is_monomer
      type(request_t) :: req(3)

      is_monomer = (fragment_idx <= n_monomers)

      if (is_monomer) then
         ! Monomer fragment
         fragment_size = 1
         allocate (fragment_indices(fragment_size))
         fragment_indices(1) = polymers(fragment_idx, 1)
      else
         ! Intersection fragment
         intersection_size = 0
         do while (intersection_size < size(intersections, 1))
            if (intersections(intersection_size + 1, fragment_idx - n_monomers) == 0) exit
            intersection_size = intersection_size + 1
         end do
         fragment_size = intersection_size
         allocate (fragment_indices(fragment_size))
         fragment_indices = intersections(1:intersection_size, fragment_idx - n_monomers)
      end if

      call isend(node_comm, fragment_idx, dest_rank, TAG_WORKER_FRAGMENT, req(1))
      call isend(node_comm, fragment_size, dest_rank, TAG_WORKER_FRAGMENT, req(2))
      call isend(node_comm, fragment_indices, dest_rank, TAG_WORKER_FRAGMENT, req(3))

      call wait(req(1))
      call wait(req(2))
      call wait(req(3))

      deallocate (fragment_indices)
   end subroutine send_gmbe_fragment_to_worker

   subroutine serial_gmbe_pie_processor(pie_atom_sets, pie_coefficients, n_pie_terms, sys_geom, method, calc_type, bonds)
      !! Serial GMBE processor using PIE coefficients
      !! Evaluates each unique atom set once and sums with PIE coefficients
      integer, intent(in) :: pie_atom_sets(:, :)  !! Unique atom sets (max_atoms, n_pie_terms)
      integer, intent(in) :: pie_coefficients(:)  !! PIE coefficient for each term
      integer, intent(in) :: n_pie_terms
      type(system_geometry_t), intent(in) :: sys_geom
      integer(int32), intent(in) :: method, calc_type
      type(bond_t), intent(in), optional :: bonds(:)

      type(physical_fragment_t) :: phys_frag
      type(calculation_result_t) :: result
      integer :: i, n_atoms, max_atoms
      integer, allocatable :: atom_list(:)
      real(dp) :: total_energy, term_energy
      real(dp), allocatable :: pie_energies(:)  !! Store individual energies for JSON output
      integer :: coeff

      call logger%info("Processing "//to_char(n_pie_terms)//" unique PIE terms...")

      total_energy = 0.0_dp
      max_atoms = size(pie_atom_sets, 1)
      allocate (pie_energies(n_pie_terms))

      do i = 1, n_pie_terms
         coeff = pie_coefficients(i)

         ! Skip terms with zero coefficient (shouldn't happen, but safety check)
         if (coeff == 0) then
            pie_energies(i) = 0.0_dp  ! Mark as skipped
            cycle
         end if

         ! Extract atom list for this term
         n_atoms = 0
         do while (n_atoms < max_atoms .and. pie_atom_sets(n_atoms + 1, i) >= 0)
            n_atoms = n_atoms + 1
         end do

         if (n_atoms == 0) then
            pie_energies(i) = 0.0_dp  ! Mark as skipped
            cycle
         end if

         allocate (atom_list(n_atoms))
         atom_list = pie_atom_sets(1:n_atoms, i)

         ! Build fragment from atom list
         call build_fragment_from_atom_list(sys_geom, atom_list, n_atoms, phys_frag, bonds)

         ! Compute energy
         call do_fragment_work(i, result, method, phys_frag, calc_type)
         term_energy = result%energy%total()

         ! Store energy for JSON output
         pie_energies(i) = term_energy

         ! Accumulate with PIE coefficient
         total_energy = total_energy + real(coeff, dp)*term_energy

         call logger%verbose("PIE term "//to_char(i)//"/"//to_char(n_pie_terms)// &
                             ": "//to_char(n_atoms)//" atoms, coeff="//to_char(coeff)// &
                             ", E="//to_char(term_energy))

         deallocate (atom_list)
         call phys_frag%destroy()
      end do

      call logger%info(" ")
      call logger%info("GMBE PIE calculation completed successfully")
      call logger%info("Final GMBE energy: "//to_char(total_energy)//" Hartree")
      call logger%info(" ")

      ! Write JSON output
      call print_gmbe_pie_json(pie_atom_sets, pie_coefficients, pie_energies, n_pie_terms, total_energy)

      deallocate (pie_energies)

   end subroutine serial_gmbe_pie_processor

end module mqc_gmbe_fragment_distribution_scheme
