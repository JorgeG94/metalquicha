module mqc_chemistry_algorithms
   use pic_types, only: int32, dp
   use pic_timer, only: timer_type
   use pic_blas_interfaces, only: pic_gemm, pic_dot
   use pic_mpi_lib, only: comm_t, send, recv, iprobe, MPI_Status, MPI_ANY_SOURCE, MPI_ANY_TAG
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char

   use mqc_mpi_tags
   use mqc_physical_fragment, only: system_geometry_t, physical_fragment_t, build_fragment_from_indices, to_angstrom
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_elements, only: element_number_to_symbol

   use mctc_env, only: wp, error_type
   use mctc_io, only: structure_type, new

   use tblite_context_type, only: context_type
   use tblite_wavefunction, only: wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only: xtb_calculator
   use tblite_xtb_gfn2, only: new_gfn2_calculator
   use tblite_xtb_singlepoint, only: xtb_singlepoint
   implicit none

contains

   subroutine process_chemistry_fragment(fragment_idx, fragment_indices, fragment_size, matrix_size, &
                                         water_energy, C_flat, phys_frag, verbosity)

      integer, intent(in) :: fragment_idx, fragment_size, matrix_size
      integer, intent(in) :: fragment_indices(fragment_size)
      real(dp), intent(out) :: water_energy
      real(dp), allocatable, intent(out) :: C_flat(:)
      type(physical_fragment_t), intent(in), optional :: phys_frag
      integer, intent(in), optional :: verbosity
      integer :: i, verb_level

      ! GFN1 calculation variables
      type(error_type), allocatable :: error
      type(structure_type) :: mol
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: num(:)
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp) :: energy
      type(context_type) :: ctx
      real(wp), parameter :: acc = 0.01_wp
      real(wp), parameter :: kt = 300.0_wp*3.166808578545117e-06_wp

      ! Set verbosity level (default is 0 for silent operation)
      if (present(verbosity)) then
         verb_level = verbosity
      else
         verb_level = 0
      end if

      ! Print fragment geometry if provided
      if (present(phys_frag)) then
         if (verb_level == 1) then
            call print_fragment_xyz(fragment_idx, phys_frag)
         end if

         allocate (num(phys_frag%n_atoms))
         allocate (xyz(3, phys_frag%n_atoms))

         ! Coordinates are already in Bohr in sys_geom, just copy them
         num = phys_frag%element_numbers
         xyz = phys_frag%coordinates

         ! Create molecular structure and run GFN1 calculation
         call new(mol, num, xyz, charge=0.0_wp, uhf=0)
         call new_gfn2_calculator(calc, mol, error)

         if (.not. allocated(error)) then
            call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
            energy = 0.0_wp
            call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=verb_level)
            water_energy = real(energy, dp)
         end if

         deallocate (num, xyz)
      else
         water_energy = 0.0_dp
      end if

      ! Return empty vector for C_flat
      allocate (C_flat(1))
      C_flat(1) = 0.0_dp
   end subroutine process_chemistry_fragment

   subroutine print_fragment_xyz(fragment_idx, phys_frag)
      !! Print fragment geometry in XYZ format
      integer, intent(in) :: fragment_idx
      type(physical_fragment_t), intent(in) :: phys_frag
      integer :: i
      character(len=2) :: symbol
      character(len=256) :: coord_line

      call logger%info("=========================================")
      call logger%info(" Fragment "//to_char(fragment_idx))
      call logger%info(" Number of atoms: "//to_char(phys_frag%n_atoms))
      call logger%info(" Coordinates in Angstroms:")
      call logger%info("-----------------------------------------")
      do i = 1, phys_frag%n_atoms
         symbol = element_number_to_symbol(phys_frag%element_numbers(i))
         ! Convert from Bohr back to Angstroms for printing
         write (coord_line, '(a2,3f15.8)') symbol, to_angstrom(phys_frag%coordinates(1:3, i))
         call logger%info(trim(coord_line))
      end do
      call logger%info("=========================================")

   end subroutine print_fragment_xyz

   subroutine compute_mbe_energy(polymers, fragment_count, max_level, energies, total_energy)
      !! Compute the many-body expansion (MBE) energy
      !! Total = sum(E(i)) + sum(deltaE(ij)) + sum(deltaE(ijk)) + ...
      !! General n-body correction:
      !! deltaE(i1,i2,...,in) = E(i1,i2,...,in) - sum of all lower-order terms
      integer, intent(in) :: polymers(:, :), fragment_count, max_level
      real(dp), intent(in) :: energies(:)
      real(dp), intent(out) :: total_energy

      integer :: i, fragment_size, body_level
      real(dp), allocatable :: sum_by_level(:)
      real(dp) :: delta_E

      allocate (sum_by_level(max_level))
      sum_by_level = 0.0_dp

      ! Sum over all fragments by their size
      do i = 1, fragment_count
         fragment_size = count(polymers(i, :) > 0)

         if (fragment_size == 1) then
            ! 1-body terms: just the monomer energies
            sum_by_level(1) = sum_by_level(1) + energies(i)
         else if (fragment_size >= 2 .and. fragment_size <= max_level) then
            ! n-body corrections for n >= 2
            delta_E = compute_delta_nbody(polymers(i, 1:fragment_size), polymers, energies, &
                                          fragment_count, fragment_size)
            sum_by_level(fragment_size) = sum_by_level(fragment_size) + delta_E
         end if
      end do

      total_energy = sum(sum_by_level)

      call logger%info("MBE Energy breakdown:")
      do body_level = 1, max_level
         if (abs(sum_by_level(body_level)) > 1e-15_dp) then
            block
               character(len=256) :: energy_line
               write (energy_line, '(a,i0,a,f20.10)') "  ", body_level, "-body:  ", sum_by_level(body_level)
               call logger%info(trim(energy_line))
            end block
         end if
      end do
      block
         character(len=256) :: total_line
         write (total_line, '(a,f20.10)') "  Total:   ", total_energy
         call logger%info(trim(total_line))
      end block

      deallocate (sum_by_level)

   end subroutine compute_mbe_energy

   recursive function compute_delta_nbody(fragment, polymers, energies, fragment_count, n) result(delta_E)
      !! Compute general n-body correction using recursion
      !! deltaE(i1,i2,...,in) = E(i1,i2,...,in)
      !!                        - sum over all proper subsets of deltaE or E
      !! For n=2: deltaE(ij) = E(ij) - E(i) - E(j)
      !! For n=3: deltaE(ijk) = E(ijk) - deltaE(ij) - deltaE(ik) - deltaE(jk) - E(i) - E(j) - E(k)
      !! For n>=4: same pattern recursively
      integer, intent(in) :: fragment(:), polymers(:, :), fragment_count, n
      real(dp), intent(in) :: energies(:)
      real(dp) :: delta_E

      integer :: idx_n, subset_size, i, j, num_subsets
      integer, allocatable :: subset(:), indices(:)
      real(dp) :: E_n, subset_contribution

      ! Find the energy of this n-mer
      idx_n = find_fragment_index(fragment, polymers, fragment_count, n)
      E_n = energies(idx_n)

      ! Start with the full n-mer energy
      delta_E = E_n

      ! Subtract contributions from all proper subsets (size 1 to n-1)
      do subset_size = 1, n - 1
         ! Generate all subsets of this size and subtract their contributions
         call generate_and_subtract_subsets(fragment, subset_size, n, polymers, energies, &
                                            fragment_count, delta_E)
      end do

   end function compute_delta_nbody

   subroutine generate_and_subtract_subsets(fragment, subset_size, n, polymers, energies, &
                                            fragment_count, delta_E)
      !! Generate all subsets of given size and subtract their contribution
      integer, intent(in) :: fragment(:), subset_size, n, polymers(:, :), fragment_count
      real(dp), intent(in) :: energies(:)
      real(dp), intent(inout) :: delta_E

      integer, allocatable :: indices(:), subset(:)
      integer :: i

      allocate (indices(subset_size))
      allocate (subset(subset_size))

      ! Initialize indices for first combination
      do i = 1, subset_size
         indices(i) = i
      end do

      ! Loop through all combinations
      do
         ! Build the current subset
         do i = 1, subset_size
            subset(i) = fragment(indices(i))
         end do

         ! Subtract this subset's contribution
         if (subset_size == 1) then
            ! 1-body: just subtract the monomer energy
            delta_E = delta_E - energies(find_fragment_index(subset, polymers, fragment_count, 1))
         else
            ! n-body: recursively compute and subtract deltaE for this subset
            delta_E = delta_E - compute_delta_nbody(subset, polymers, energies, fragment_count, subset_size)
         end if

         ! Get next combination
         if (.not. next_combination(indices, subset_size, n)) exit
      end do

      deallocate (indices, subset)

   end subroutine generate_and_subtract_subsets

   function next_combination(indices, k, n) result(has_next)
      !! Generate next combination (updates indices in place)
      !! Returns .true. if there's a next combination, .false. if we're done
      integer, intent(inout) :: indices(:)
      integer, intent(in) :: k, n
      logical :: has_next
      integer :: i

      has_next = .true.

      ! Find rightmost index that can be incremented
      i = k
      do while (i >= 1)
         if (indices(i) < n - k + i) then
            indices(i) = indices(i) + 1
            ! Reset all indices to the right
            do while (i < k)
               i = i + 1
               indices(i) = indices(i - 1) + 1
            end do
            return
         end if
         i = i - 1
      end do

      ! No more combinations
      has_next = .false.

   end function next_combination

   function find_fragment_index(target_monomers, polymers, fragment_count, expected_size) result(idx)
      !! Find the fragment index that contains exactly the target monomers
      integer, intent(in) :: target_monomers(:), polymers(:, :), fragment_count, expected_size
      integer :: idx

      integer :: i, j, fragment_size
      logical :: match

      idx = -1

      do i = 1, fragment_count
         fragment_size = count(polymers(i, :) > 0)

         ! Check if this fragment has the right size
         if (fragment_size /= expected_size) cycle

         ! Check if all target monomers are in this fragment
         match = .true.
         do j = 1, expected_size
            if (.not. any(polymers(i, 1:fragment_size) == target_monomers(j))) then
               match = .false.
               exit
            end if
         end do

         if (match) then
            idx = i
            return
         end if
      end do

      ! If we get here, we didn't find the fragment
      block
         character(len=256) :: monomers_str
         integer :: k
         write (monomers_str, '(*(i0,1x))') (target_monomers(k), k=1, size(target_monomers))
         call logger%error("Could not find fragment with monomers: "//trim(monomers_str))
      end block
      error stop "Fragment not found in find_fragment_index"

   end function find_fragment_index

   subroutine calculate_exact_flops(polymers, fragment_count, max_level, matrix_size, total_flops)
      integer, intent(in) :: polymers(:, :), fragment_count, max_level, matrix_size
      real(dp), intent(out) :: total_flops

      integer :: i, fragment_size
      integer, allocatable :: n_mers(:)
      real(dp), allocatable :: mer_flops(:)
      real(dp) :: mer_size

      ! Allocate counters for each n-mer level (1 to max_level)
      allocate (n_mers(max_level))
      allocate (mer_flops(max_level))

      n_mers = 0
      mer_flops = 0.0_dp

      ! Count fragments by size
      do i = 1, fragment_count
         fragment_size = count(polymers(i, :) > 0)
         if (fragment_size >= 1 .and. fragment_size <= max_level) then
            n_mers(fragment_size) = n_mers(fragment_size) + 1
         end if
      end do

      ! Calculate FLOPs for each n-mer level (using JAXPY: C = alpha * H + S, which is 2*N^2 FLOPs)
      do i = 1, max_level
         mer_size = real(i*matrix_size, dp)
         mer_flops(i) = real(n_mers(i), dp)*2.0_dp*mer_size**2  ! 2*N^2 for JAXPY
      end do

      ! Total FLOPs
      total_flops = sum(mer_flops)

      ! Print breakdown
      call logger%info("Fragment breakdown:")
      do i = 1, max_level
         if (n_mers(i) > 0) then
            block
               character(len=256) :: flop_line
               write (flop_line, '(a,i0,a,i0,a,f12.3,a)') "  ", i, "-mers:  ", n_mers(i), &
                  " (", mer_flops(i)/1.0e9_dp, " GFLOP)"
               call logger%info(trim(flop_line))
            end block
         end if
      end do
      block
         character(len=256) :: total_line
         write (total_line, '(a,i0,a,f12.3,a)') "  Total:    ", fragment_count, &
            " (", total_flops/1.0e9_dp, " GFLOP)"
         call logger%info(trim(total_line))
      end block

      deallocate (n_mers, mer_flops)

   end subroutine calculate_exact_flops

   subroutine global_coordinator(world_comm, node_comm, total_fragments, polymers, max_level, &
                                 node_leader_ranks, num_nodes, matrix_size)
      type(comm_t), intent(in) :: world_comm, node_comm
      integer, intent(in) :: total_fragments, max_level, num_nodes, matrix_size
      integer, intent(in) :: polymers(:, :), node_leader_ranks(:)

      integer :: current_fragment, finished_nodes
      integer :: request_source, dummy_msg, fragment_idx
      type(MPI_Status) :: status, local_status
      logical :: handling_local_workers
      logical :: has_pending

      ! For local workers
      integer :: local_finished_workers, fragment_size, local_dummy
      integer, allocatable :: fragment_indices(:)

      ! Storage for results
      real(dp), allocatable :: scalar_results(:)
      real(dp), allocatable :: matrix_results(:, :)
      real(dp), allocatable :: temp_matrix(:)
      integer :: max_matrix_size
      integer :: worker_fragment_map(node_comm%size())
      integer :: worker_source
      integer :: results_received

      current_fragment = total_fragments
      finished_nodes = 0
      local_finished_workers = 0
      handling_local_workers = (node_comm%size() > 1)
      results_received = 0

      ! Allocate storage for results
      allocate (scalar_results(total_fragments))
      max_matrix_size = (max_level*matrix_size)**2
      allocate (matrix_results(max_matrix_size, total_fragments))
      scalar_results = 0.0_dp
      matrix_results = 0.0_dp
      worker_fragment_map = 0

      call logger%verbose("Global coordinator starting with "//to_char(total_fragments)// &
                          " fragments for "//to_char(num_nodes)//" nodes")

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

               ! Receive scalar result and store it using the fragment index for this worker
               call recv(node_comm, scalar_results(worker_fragment_map(worker_source)), worker_source, TAG_WORKER_SCALAR_RESULT)
               ! Receive matrix result into temporary array, then copy to storage
               call recv(node_comm, temp_matrix, worker_source, TAG_WORKER_MATRIX_RESULT, local_status)
               ! Copy only the received size, pad rest with zeros
               matrix_results(1:size(temp_matrix), worker_fragment_map(worker_source)) = temp_matrix
               if (size(temp_matrix) < max_matrix_size) then
                  matrix_results(size(temp_matrix) + 1:max_matrix_size, worker_fragment_map(worker_source)) = 0.0_dp
               end if
               ! Clear the mapping since we've received the result
               worker_fragment_map(worker_source) = 0
               if (allocated(temp_matrix)) deallocate (temp_matrix)
               results_received = results_received + 1
            end do
         end if

         ! PRIORITY 1b: Check for incoming results from remote node coordinators
         do
            call iprobe(world_comm, MPI_ANY_SOURCE, TAG_NODE_SCALAR_RESULT, has_pending, status)
            if (.not. has_pending) exit

            ! Receive fragment index, scalar result, and matrix result from node coordinator
            call recv(world_comm, fragment_idx, status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT)
            call recv(world_comm, scalar_results(fragment_idx), status%MPI_SOURCE, TAG_NODE_SCALAR_RESULT)
            call recv(world_comm, temp_matrix, status%MPI_SOURCE, TAG_NODE_MATRIX_RESULT, status)

            ! Copy matrix result into storage
            matrix_results(1:size(temp_matrix), fragment_idx) = temp_matrix
            if (size(temp_matrix) < max_matrix_size) then
               matrix_results(size(temp_matrix) + 1:max_matrix_size, fragment_idx) = 0.0_dp
            end if
            if (allocated(temp_matrix)) deallocate (temp_matrix)
            results_received = results_received + 1
         end do

         ! PRIORITY 2: Remote node coordinator requests
         call iprobe(world_comm, MPI_ANY_SOURCE, TAG_NODE_REQUEST, has_pending, status)
         if (has_pending) then
            call recv(world_comm, dummy_msg, status%MPI_SOURCE, TAG_NODE_REQUEST)
            request_source = status%MPI_SOURCE

            if (current_fragment >= 1) then
               call send_fragment_to_node(world_comm, current_fragment, polymers, max_level, request_source)
               current_fragment = current_fragment - 1
            else
               call send(world_comm, -1, request_source, TAG_NODE_FINISH)
               finished_nodes = finished_nodes + 1
            end if
         end if

         ! PRIORITY 3: Local workers (shared memory) - send new work
         if (handling_local_workers .and. local_finished_workers < node_comm%size() - 1) then
            call iprobe(node_comm, MPI_ANY_SOURCE, TAG_WORKER_REQUEST, has_pending, local_status)
            if (has_pending) then
               ! Only process work request if this worker doesn't have pending results
               if (worker_fragment_map(local_status%MPI_SOURCE) == 0) then
                  call recv(node_comm, local_dummy, local_status%MPI_SOURCE, TAG_WORKER_REQUEST)

                  if (current_fragment >= 1) then
                     call send_fragment_to_worker(node_comm, current_fragment, polymers, max_level, &
                                                  local_status%MPI_SOURCE, matrix_size)
                     ! Track which fragment was sent to this worker
                     worker_fragment_map(local_status%MPI_SOURCE) = current_fragment
                     current_fragment = current_fragment - 1
                  else
                     call send(node_comm, -1, local_status%MPI_SOURCE, TAG_WORKER_FINISH)
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
      block
         real(dp) :: mbe_total_energy

         ! Compute the many-body expansion energy
         call logger%info("")
         call logger%info("Computing Many-Body Expansion (MBE)...")
         call compute_mbe_energy(polymers, total_fragments, max_level, scalar_results, mbe_total_energy)

      end block

      ! Cleanup
      deallocate (scalar_results, matrix_results)
   end subroutine global_coordinator

   subroutine send_fragment_to_node(world_comm, fragment_idx, polymers, max_level, dest_rank)
      type(comm_t), intent(in) :: world_comm
      integer, intent(in) :: fragment_idx, max_level, dest_rank
      integer, intent(in) :: polymers(:, :)
      integer :: fragment_size
      integer, allocatable :: fragment_indices(:)

      fragment_size = count(polymers(fragment_idx, :) > 0)
      allocate (fragment_indices(fragment_size))
      fragment_indices = polymers(fragment_idx, 1:fragment_size)

      call send(world_comm, fragment_idx, dest_rank, TAG_NODE_FRAGMENT)
      call send(world_comm, fragment_size, dest_rank, TAG_NODE_FRAGMENT)
      call send(world_comm, fragment_indices, dest_rank, TAG_NODE_FRAGMENT)

      deallocate (fragment_indices)
   end subroutine send_fragment_to_node

   subroutine send_fragment_to_worker(node_comm, fragment_idx, polymers, max_level, dest_rank, matrix_size)
      type(comm_t), intent(in) :: node_comm
      integer, intent(in) :: fragment_idx, max_level, dest_rank, matrix_size
      integer, intent(in) :: polymers(:, :)
      integer :: fragment_size
      integer, allocatable :: fragment_indices(:)

      fragment_size = count(polymers(fragment_idx, :) > 0)
      allocate (fragment_indices(fragment_size))
      fragment_indices = polymers(fragment_idx, 1:fragment_size)

      call send(node_comm, fragment_idx, dest_rank, TAG_WORKER_FRAGMENT)
      call send(node_comm, fragment_size, dest_rank, TAG_WORKER_FRAGMENT)
      call send(node_comm, fragment_indices, dest_rank, TAG_WORKER_FRAGMENT)
      call send(node_comm, matrix_size, dest_rank, TAG_WORKER_FRAGMENT)

      deallocate (fragment_indices)
   end subroutine send_fragment_to_worker

   subroutine node_coordinator(world_comm, node_comm, max_level, matrix_size)
      class(comm_t), intent(in) :: world_comm, node_comm
      integer(int32), intent(in) :: max_level, matrix_size

      integer(int32) :: fragment_idx, fragment_size, dummy_msg
      integer(int32) :: finished_workers
      integer(int32), allocatable :: fragment_indices(:)
      type(MPI_Status) :: status, global_status
      logical :: local_message_pending, more_fragments, has_result
      integer(int32) :: local_dummy

      ! For tracking worker-fragment mapping and collecting results
      integer(int32) :: worker_fragment_map(node_comm%size())
      integer(int32) :: worker_source
      real(dp) :: scalar_result
      real(dp), allocatable :: matrix_result(:)

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

            ! Receive scalar and matrix results from worker
            call recv(node_comm, scalar_result, worker_source, TAG_WORKER_SCALAR_RESULT)
            call recv(node_comm, matrix_result, worker_source, TAG_WORKER_MATRIX_RESULT, status)

            ! Forward results to global coordinator with fragment index
            call send(world_comm, worker_fragment_map(worker_source), 0, TAG_NODE_SCALAR_RESULT)  ! fragment_idx
            call send(world_comm, scalar_result, 0, TAG_NODE_SCALAR_RESULT)                       ! scalar result
            call send(world_comm, matrix_result, 0, TAG_NODE_MATRIX_RESULT)                       ! matrix result

            ! Clear the mapping and deallocate
            worker_fragment_map(worker_source) = 0
            if (allocated(matrix_result)) deallocate (matrix_result)
         end if

         ! PRIORITY 2: Check for work requests from local workers
         call iprobe(node_comm, MPI_ANY_SOURCE, TAG_WORKER_REQUEST, local_message_pending, status)

         if (local_message_pending) then
            ! Only process work request if this worker doesn't have pending results
            if (worker_fragment_map(status%MPI_SOURCE) == 0) then
               call recv(node_comm, local_dummy, status%MPI_SOURCE, TAG_WORKER_REQUEST)

               if (more_fragments) then
                  call send(world_comm, dummy_msg, 0, TAG_NODE_REQUEST)
                  call recv(world_comm, fragment_idx, 0, MPI_ANY_TAG, global_status)

                  if (global_status%MPI_TAG == TAG_NODE_FRAGMENT) then
                     call recv(world_comm, fragment_size, 0, TAG_NODE_FRAGMENT, global_status)
                     allocate (fragment_indices(fragment_size))
                     call recv(world_comm, fragment_indices, 0, TAG_NODE_FRAGMENT, global_status)

                     call send(node_comm, fragment_idx, status%MPI_SOURCE, TAG_WORKER_FRAGMENT)
                     call send(node_comm, fragment_size, status%MPI_SOURCE, TAG_WORKER_FRAGMENT)
                     call send(node_comm, fragment_indices, status%MPI_SOURCE, TAG_WORKER_FRAGMENT)
                     call send(node_comm, matrix_size, status%MPI_SOURCE, TAG_WORKER_FRAGMENT)

                     ! Track which fragment was sent to this worker
                     worker_fragment_map(status%MPI_SOURCE) = fragment_idx

                     deallocate (fragment_indices)
                  else
                     call send(node_comm, -1, status%MPI_SOURCE, TAG_WORKER_FINISH)
                     finished_workers = finished_workers + 1
                     more_fragments = .false.
                  end if
               else
                  call send(node_comm, -1, status%MPI_SOURCE, TAG_WORKER_FINISH)
                  finished_workers = finished_workers + 1
               end if
            end if
         end if
      end do
   end subroutine node_coordinator

   subroutine node_worker(world_comm, node_comm, max_level, sys_geom)
      class(comm_t), intent(in) :: world_comm, node_comm
      integer, intent(in) :: max_level
      type(system_geometry_t), intent(in), optional :: sys_geom

      integer(int32) :: fragment_idx, fragment_size, dummy_msg, matrix_size
      integer(int32), allocatable :: fragment_indices(:)
      real(dp) :: dot_result
      real(dp), allocatable :: C_flat(:)
      type(MPI_Status) :: status
      type(physical_fragment_t) :: phys_frag

      dummy_msg = 0

      do
         call send(node_comm, dummy_msg, 0, TAG_WORKER_REQUEST)
         call recv(node_comm, fragment_idx, 0, MPI_ANY_TAG, status)

         select case (status%MPI_TAG)
         case (TAG_WORKER_FRAGMENT)
            call recv(node_comm, fragment_size, 0, TAG_WORKER_FRAGMENT, status)
            allocate (fragment_indices(fragment_size))
            call recv(node_comm, fragment_indices, 0, TAG_WORKER_FRAGMENT, status)
            call recv(node_comm, matrix_size, 0, TAG_WORKER_FRAGMENT, status)

            ! Build physical fragment from indices if sys_geom is available
            if (present(sys_geom)) then
               call build_fragment_from_indices(sys_geom, fragment_indices, phys_frag)

               ! Process the chemistry fragment with physical geometry
               call process_chemistry_fragment(fragment_idx, fragment_indices, fragment_size, matrix_size, &
                                               dot_result, C_flat, phys_frag)

               call phys_frag%destroy()
            else
               ! Process without physical geometry (old behavior)
               call process_chemistry_fragment(fragment_idx, fragment_indices, fragment_size, matrix_size, &
                                               dot_result, C_flat)
            end if

            ! Send results back to coordinator
            call send(node_comm, dot_result, 0, TAG_WORKER_SCALAR_RESULT)
            call send(node_comm, C_flat, 0, TAG_WORKER_MATRIX_RESULT)

            deallocate (fragment_indices, C_flat)
         case (TAG_WORKER_FINISH)
            exit
         end select
      end do
   end subroutine node_worker

   subroutine unfragmented_calculation(sys_geom)
      !! Run unfragmented calculation on the entire system (nlevel=0)
      !! This is a simple single-process calculation without MPI distribution
      type(system_geometry_t), intent(in), optional :: sys_geom

      real(dp) :: dot_result
      real(dp), allocatable :: C_flat(:)
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

      ! Build the full system as a single fragment (all monomers)
      block
         integer, allocatable :: all_monomer_indices(:)

         allocate (all_monomer_indices(sys_geom%n_monomers))
         do i = 1, sys_geom%n_monomers
            all_monomer_indices(i) = i
         end do

         call build_fragment_from_indices(sys_geom, all_monomer_indices, full_system)
         deallocate (all_monomer_indices)
      end block

      ! Process the full system with verbosity=1 for detailed output
      block
         integer, allocatable :: temp_indices(:)
         allocate (temp_indices(sys_geom%n_monomers))
         do i = 1, sys_geom%n_monomers
            temp_indices(i) = i
         end do

         call process_chemistry_fragment(0, temp_indices, sys_geom%n_monomers, &
                                         total_atoms, dot_result, C_flat, &
                                         phys_frag=full_system, verbosity=1)
         deallocate (temp_indices)
      end block

      call logger%info("============================================")
      call logger%info("Unfragmented calculation completed")
      block
         character(len=256) :: result_line
         write (result_line, '(a,es15.8)') "  Scalar result: ", dot_result
         call logger%info(trim(result_line))
      end block
      call logger%info("  Matrix size: "//to_char(size(C_flat)))
      call logger%info("============================================")

      if (allocated(C_flat)) deallocate (C_flat)

   end subroutine unfragmented_calculation

end module mqc_chemistry_algorithms
