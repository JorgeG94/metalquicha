!! Generalized Many-Body Expansion (GMBE) fragment distribution module
module mqc_gmbe_fragment_distribution_scheme
   !! Implements fragment distribution schemes for GMBE calculations with overlapping fragments
   !! Handles both serial and MPI-parallelized distribution of monomers and intersection fragments
   use pic_types, only: int32, int64, dp
   use pic_timer, only: timer_type
   use mqc_calc_types, only: CALC_TYPE_GRADIENT
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
   public :: serial_gmbe_pie_processor  !! PIE-based serial processor
   public :: gmbe_pie_coordinator  !! PIE-based MPI coordinator

contains

   subroutine serial_gmbe_pie_processor(pie_atom_sets, pie_coefficients, n_pie_terms, sys_geom, method, calc_type, bonds)
      !! Serial GMBE processor using PIE coefficients
      !! Evaluates each unique atom set once and sums with PIE coefficients
      !! Supports both energy-only and energy+gradient calculations
      use mqc_calc_types, only: CALC_TYPE_GRADIENT, CALC_TYPE_ENERGY, calc_type_to_string
      use mqc_physical_fragment, only: redistribute_cap_gradients
      use pic_logger, only: info_level
      integer, intent(in) :: pie_atom_sets(:, :)  !! Unique atom sets (max_atoms, n_pie_terms)
      integer, intent(in) :: pie_coefficients(:)  !! PIE coefficient for each term
      integer, intent(in) :: n_pie_terms
      type(system_geometry_t), intent(in) :: sys_geom
      integer(int32), intent(in) :: method, calc_type
      type(bond_t), intent(in), optional :: bonds(:)

      type(physical_fragment_t) :: phys_frag
      type(calculation_result_t), allocatable :: results(:)
      integer :: i, n_atoms, max_atoms, iatom, current_log_level
      integer, allocatable :: atom_list(:)
      real(dp) :: total_energy, term_energy
      real(dp), allocatable :: pie_energies(:)  !! Store individual energies for JSON output
      real(dp), allocatable :: total_gradient(:, :)  !! Total gradient (3, total_atoms)
      real(dp), allocatable :: term_gradient(:, :)  !! Temporary gradient for each term
      integer :: coeff

      call logger%info("Processing "//to_char(n_pie_terms)//" unique PIE terms...")
      call logger%info("  Calculation type: "//calc_type_to_string(calc_type))

      total_energy = 0.0_dp
      max_atoms = size(pie_atom_sets, 1)
      allocate (pie_energies(n_pie_terms))
      allocate (results(n_pie_terms))

      ! Allocate gradient arrays if needed
      if (calc_type == CALC_TYPE_GRADIENT) then
         allocate (total_gradient(3, sys_geom%total_atoms))
         allocate (term_gradient(3, sys_geom%total_atoms))
         total_gradient = 0.0_dp
      end if

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

         ! Compute energy (and gradient if requested)
         call do_fragment_work(i, results(i), method, phys_frag, calc_type)
         term_energy = results(i)%energy%total()

         ! Store energy for JSON output
         pie_energies(i) = term_energy

         ! Accumulate with PIE coefficient
         total_energy = total_energy + real(coeff, dp)*term_energy

         ! Accumulate gradient if present
         if (calc_type == CALC_TYPE_GRADIENT .and. results(i)%has_gradient) then
            ! Map fragment gradient to system coordinates with proper cap handling
            term_gradient = 0.0_dp
            call redistribute_cap_gradients(phys_frag, results(i)%gradient, term_gradient)

            ! Accumulate with PIE coefficient
            total_gradient = total_gradient + real(coeff, dp)*term_gradient
         end if

         call logger%verbose("PIE term "//to_char(i)//"/"//to_char(n_pie_terms)// &
                             ": "//to_char(n_atoms)//" atoms, coeff="//to_char(coeff)// &
                             ", E="//to_char(term_energy))

         deallocate (atom_list)
         call phys_frag%destroy()
      end do

      call logger%info(" ")
      call logger%info("GMBE PIE calculation completed successfully")
      call logger%info("Final GMBE energy: "//to_char(total_energy)//" Hartree")

      ! Print gradient info if computed
      if (calc_type == CALC_TYPE_GRADIENT) then
         call logger%info("GMBE PIE gradient computation completed")
         call logger%info("  Total gradient norm: "//to_char(sqrt(sum(total_gradient**2))))

         ! Print detailed gradient if info level and small system
         call logger%configuration(level=current_log_level)
         if (current_log_level >= info_level .and. sys_geom%total_atoms < 100) then
            call logger%info(" ")
            call logger%info("Total GMBE PIE Gradient (Hartree/Bohr):")
            do iatom = 1, sys_geom%total_atoms
               block
                  character(len=256) :: grad_line
                  write (grad_line, '(a,i5,a,3f20.12)') "  Atom ", iatom, ": ", &
                     total_gradient(1, iatom), total_gradient(2, iatom), total_gradient(3, iatom)
                  call logger%info(trim(grad_line))
               end block
            end do
            call logger%info(" ")
         end if
      end if

      call logger%info(" ")

      ! Write JSON output
      call print_gmbe_pie_json(pie_atom_sets, pie_coefficients, pie_energies, n_pie_terms, total_energy)

      deallocate (pie_energies, results)
      if (allocated(total_gradient)) deallocate (total_gradient)
      if (allocated(term_gradient)) deallocate (term_gradient)

   end subroutine serial_gmbe_pie_processor

   subroutine gmbe_pie_coordinator(world_comm, node_comm, pie_atom_sets, pie_coefficients, n_pie_terms, &
                                   node_leader_ranks, num_nodes, sys_geom, method, calc_type, bonds)
      !! MPI coordinator for PIE-based GMBE calculations
      !! Distributes PIE terms across MPI ranks and accumulates results
      use mqc_calc_types, only: CALC_TYPE_GRADIENT
      use mqc_physical_fragment, only: redistribute_cap_gradients

      type(comm_t), intent(in) :: world_comm, node_comm
      integer, intent(in) :: pie_atom_sets(:, :)  !! Unique atom sets (max_atoms, n_pie_terms)
      integer, intent(in) :: pie_coefficients(:)  !! PIE coefficient for each term
      integer, intent(in) :: n_pie_terms
      integer, intent(in) :: node_leader_ranks(:), num_nodes
      type(system_geometry_t), intent(in) :: sys_geom
      integer(int32), intent(in) :: method, calc_type
      type(bond_t), intent(in), optional :: bonds(:)

      type(timer_type) :: coord_timer
      integer :: current_term_idx, results_received, finished_nodes
      integer :: request_source, dummy_msg, term_idx
      type(MPI_Status) :: status, local_status
      logical :: handling_local_workers, has_pending
      integer :: local_finished_workers, local_dummy

      ! Storage for results
      type(calculation_result_t), allocatable :: results(:)
      integer :: worker_term_map(node_comm%size())
      integer :: worker_source
      real(dp) :: total_energy
      real(dp), allocatable :: total_gradient(:, :)

      ! MPI request handles
      type(request_t) :: req

      current_term_idx = n_pie_terms
      finished_nodes = 0
      local_finished_workers = 0
      handling_local_workers = (node_comm%size() > 1)
      results_received = 0
      worker_term_map = 0

      allocate (results(n_pie_terms))

      call logger%verbose("GMBE PIE coordinator starting with "//to_char(n_pie_terms)// &
                          " PIE terms for "//to_char(num_nodes)//" nodes")

      call coord_timer%start()
      do while (finished_nodes < num_nodes)

         ! PRIORITY 1: Check for incoming results from local workers
         if (handling_local_workers) then
            do
               call iprobe(node_comm, MPI_ANY_SOURCE, TAG_WORKER_SCALAR_RESULT, has_pending, local_status)
               if (.not. has_pending) exit

               worker_source = local_status%MPI_SOURCE
               if (worker_term_map(worker_source) == 0) then
                  call logger%error("Received result from worker "//to_char(worker_source)// &
                                    " but no term was assigned!")
                  error stop "Invalid worker_term_map state"
               end if

               call result_irecv(results(worker_term_map(worker_source)), node_comm, worker_source, &
                                 TAG_WORKER_SCALAR_RESULT, req)
               call wait(req)
               worker_term_map(worker_source) = 0
               results_received = results_received + 1
               if (mod(results_received, max(1, n_pie_terms/10)) == 0 .or. results_received == n_pie_terms) then
                  call logger%info("  Processed "//to_char(results_received)//"/"// &
                                   to_char(n_pie_terms)//" PIE terms ["// &
                                   to_char(coord_timer%get_elapsed_time())//" s]")
               end if
            end do
         end if

         ! PRIORITY 1b: Check for incoming results from remote node coordinators
         do
            call iprobe(world_comm, MPI_ANY_SOURCE, TAG_NODE_SCALAR_RESULT, has_pending, status)
            if (.not. has_pending) exit

            call irecv(world_comm, term_idx, status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT, req)
            call wait(req)
            call result_irecv(results(term_idx), world_comm, status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT, req)
            call wait(req)
            results_received = results_received + 1
            if (mod(results_received, max(1, n_pie_terms/10)) == 0 .or. results_received == n_pie_terms) then
               call logger%info("  Processed "//to_char(results_received)//"/"// &
                                to_char(n_pie_terms)//" PIE terms ["// &
                                to_char(coord_timer%get_elapsed_time())//" s]")
            end if
         end do

         ! PRIORITY 2: Remote node coordinator requests
         call iprobe(world_comm, MPI_ANY_SOURCE, TAG_NODE_REQUEST, has_pending, status)
         if (has_pending) then
            call irecv(world_comm, dummy_msg, status%MPI_SOURCE, TAG_NODE_REQUEST, req)
            call wait(req)
            request_source = status%MPI_SOURCE

            if (current_term_idx >= 1) then
               call send_pie_term_to_node(world_comm, current_term_idx, pie_atom_sets, request_source)
               current_term_idx = current_term_idx - 1
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
               if (worker_term_map(local_status%MPI_SOURCE) == 0) then
                  call irecv(node_comm, local_dummy, local_status%MPI_SOURCE, TAG_WORKER_REQUEST, req)
                  call wait(req)

                  if (current_term_idx >= 1) then
                     call send_pie_term_to_worker(node_comm, current_term_idx, pie_atom_sets, local_status%MPI_SOURCE)
                     worker_term_map(local_status%MPI_SOURCE) = current_term_idx
                     current_term_idx = current_term_idx - 1
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
             .and. results_received >= n_pie_terms) then
            handling_local_workers = .false.
            finished_nodes = finished_nodes + 1
         end if
      end do

      call logger%verbose("GMBE PIE coordinator finished all terms")
      call coord_timer%stop()
      call logger%info("Time to evaluate all PIE terms "//to_char(coord_timer%get_elapsed_time())//" s")

      ! Accumulate results with PIE coefficients
      call logger%info(" ")
      call logger%info("Computing GMBE PIE energy...")
      call coord_timer%start()

      total_energy = 0.0_dp
      do term_idx = 1, n_pie_terms
         total_energy = total_energy + real(pie_coefficients(term_idx), dp)*results(term_idx)%energy%total()
      end do

      ! Handle gradients if computed
      if (calc_type == CALC_TYPE_GRADIENT) then
         allocate (total_gradient(3, sys_geom%total_atoms))
         total_gradient = 0.0_dp

         do term_idx = 1, n_pie_terms
            if (results(term_idx)%has_gradient) then
               ! Map fragment gradient to system coordinates
               block
                  real(dp), allocatable :: term_gradient(:, :)
                  type(physical_fragment_t) :: phys_frag
                  integer :: n_atoms, max_atoms
                  integer, allocatable :: atom_list(:)

                  allocate (term_gradient(3, sys_geom%total_atoms))
                  term_gradient = 0.0_dp

                  ! Extract atom list for this term
                  max_atoms = size(pie_atom_sets, 1)
                  n_atoms = 0
                  do while (n_atoms < max_atoms .and. pie_atom_sets(n_atoms + 1, term_idx) >= 0)
                     n_atoms = n_atoms + 1
                  end do

                  if (n_atoms > 0) then
                     allocate (atom_list(n_atoms))
                     atom_list = pie_atom_sets(1:n_atoms, term_idx)

                     ! Build fragment to get proper mapping
                     call build_fragment_from_atom_list(sys_geom, atom_list, n_atoms, phys_frag, bonds)
                     call redistribute_cap_gradients(phys_frag, results(term_idx)%gradient, term_gradient)
                     call phys_frag%destroy()
                     deallocate (atom_list)
                  end if

                  ! Accumulate with PIE coefficient
                  total_gradient = total_gradient + real(pie_coefficients(term_idx), dp)*term_gradient
                  deallocate (term_gradient)
               end block
            end if
         end do

         ! Print gradient information
         call logger%info("GMBE PIE gradient computation completed")
         call logger%info("  Total gradient norm: "//to_char(sqrt(sum(total_gradient**2))))

         ! Print detailed gradient if info level and small system
         block
            use pic_logger, only: info_level
            integer :: iatom, current_log_level
            call logger%configuration(level=current_log_level)
            if (current_log_level >= info_level .and. sys_geom%total_atoms < 100) then
               call logger%info(" ")
               call logger%info("Total GMBE PIE Gradient (Hartree/Bohr):")
               do iatom = 1, sys_geom%total_atoms
                  block
                     character(len=256) :: grad_line
                     write (grad_line, '(a,i5,a,3f20.12)') "  Atom ", iatom, ": ", &
                        total_gradient(1, iatom), total_gradient(2, iatom), total_gradient(3, iatom)
                     call logger%info(trim(grad_line))
                  end block
               end do
               call logger%info(" ")
            end if
         end block

         deallocate (total_gradient)
      end if

      call coord_timer%stop()
      call logger%info("Time to compute GMBE PIE "//to_char(coord_timer%get_elapsed_time())//" s")
      call logger%info(" ")
      call logger%info("GMBE PIE calculation completed successfully")
      call logger%info("Final GMBE energy: "//to_char(total_energy)//" Hartree")
      call logger%info(" ")

      ! Write JSON output (reuse existing function)
      block
         real(dp), allocatable :: pie_energies(:)
         allocate (pie_energies(n_pie_terms))
         do term_idx = 1, n_pie_terms
            pie_energies(term_idx) = results(term_idx)%energy%total()
         end do
         call print_gmbe_pie_json(pie_atom_sets, pie_coefficients, pie_energies, n_pie_terms, total_energy)
         deallocate (pie_energies)
      end block

      deallocate (results)

   end subroutine gmbe_pie_coordinator

   subroutine send_pie_term_to_node(world_comm, term_idx, pie_atom_sets, dest_rank)
      !! Send PIE term (atom list) to remote node coordinator
      type(comm_t), intent(in) :: world_comm
      integer, intent(in) :: term_idx
      integer, intent(in) :: pie_atom_sets(:, :)
      integer, intent(in) :: dest_rank

      integer :: n_atoms, max_atoms
      integer, allocatable :: atom_list(:)
      integer(int32) :: fragment_type
      type(request_t) :: req(4)

      ! PIE terms always use atom lists (type 1)
      fragment_type = 1

      ! Extract atom list for this term
      max_atoms = size(pie_atom_sets, 1)
      n_atoms = 0
      do while (n_atoms < max_atoms .and. pie_atom_sets(n_atoms + 1, term_idx) >= 0)
         n_atoms = n_atoms + 1
      end do

      allocate (atom_list(n_atoms))
      atom_list = pie_atom_sets(1:n_atoms, term_idx)

      call isend(world_comm, term_idx, dest_rank, TAG_NODE_FRAGMENT, req(1))
      call isend(world_comm, fragment_type, dest_rank, TAG_NODE_FRAGMENT, req(2))
      call isend(world_comm, n_atoms, dest_rank, TAG_NODE_FRAGMENT, req(3))
      call isend(world_comm, atom_list, dest_rank, TAG_NODE_FRAGMENT, req(4))

      call wait(req(1))
      call wait(req(2))
      call wait(req(3))
      call wait(req(4))

      deallocate (atom_list)
   end subroutine send_pie_term_to_node

   subroutine send_pie_term_to_worker(node_comm, term_idx, pie_atom_sets, dest_rank)
      !! Send PIE term (atom list) to local worker
      type(comm_t), intent(in) :: node_comm
      integer, intent(in) :: term_idx
      integer, intent(in) :: pie_atom_sets(:, :)
      integer, intent(in) :: dest_rank

      integer :: n_atoms, max_atoms
      integer, allocatable :: atom_list(:)
      integer(int32) :: fragment_type
      type(request_t) :: req(4)

      ! PIE terms always use atom lists (type 1)
      fragment_type = 1

      ! Extract atom list for this term
      max_atoms = size(pie_atom_sets, 1)
      n_atoms = 0
      do while (n_atoms < max_atoms .and. pie_atom_sets(n_atoms + 1, term_idx) >= 0)
         n_atoms = n_atoms + 1
      end do

      allocate (atom_list(n_atoms))
      atom_list = pie_atom_sets(1:n_atoms, term_idx)

      call isend(node_comm, term_idx, dest_rank, TAG_WORKER_FRAGMENT, req(1))
      call isend(node_comm, fragment_type, dest_rank, TAG_WORKER_FRAGMENT, req(2))
      call isend(node_comm, n_atoms, dest_rank, TAG_WORKER_FRAGMENT, req(3))
      call isend(node_comm, atom_list, dest_rank, TAG_WORKER_FRAGMENT, req(4))

      call wait(req(1))
      call wait(req(2))
      call wait(req(3))
      call wait(req(4))

      deallocate (atom_list)
   end subroutine send_pie_term_to_worker

end module mqc_gmbe_fragment_distribution_scheme
