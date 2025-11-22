module pic_chemistry_algorithms
   use mpi_f08
   use pic_types
   use pic_timer
   use mpi_comm_simple
   use pic_fragment, only: pic_fragment_block
   use pic_blas_interfaces, only: pic_gemm, pic_dot 
   implicit none

contains

   subroutine process_chemistry_fragment(fragment_idx, fragment_indices, fragment_size, matrix_size, &
                                          water_energy, C_flat)
      integer, intent(in) :: fragment_idx, fragment_size, matrix_size
      integer, intent(in) :: fragment_indices(fragment_size)
      real(dp), intent(out) :: water_energy
      real(dp), allocatable, intent(out) :: C_flat(:)
      real(dp) :: dot_result
      real(dp), allocatable :: H(:, :), S(:, :), C(:, :)
      real(dp), allocatable :: H_flat(:), S_flat(:)
      integer :: i, j, dims
      type(timer_type) :: compute_timer
      real(dp) :: elapsed_time
      real(dp), parameter :: alpha = 17.0_dp
      real(dp), parameter :: water_1 = -75.0_dp

      dims = fragment_size*matrix_size
      water_energy = water_1 * fragment_size

      ! Allocate matrices for Hamiltonian, Overlap, and Coefficients
      allocate (H(dims, dims), S(dims, dims), C(dims, dims))

      ! Initialize matrices (simulating chemistry computation)
      ! Scale values to get dot product in reasonable range
      do concurrent(j=1:dims, i=1:dims)
         H(i, j) = real(fragment_idx, dp) * 1e-4_dp  ! Hamiltonian matrix
         S(i, j) = real(fragment_idx, dp) * 1e-4_dp  ! Overlap matrix
         C(i, j) = 0.0_dp                             ! Coefficient matrix
      end do

      ! Simulate chemistry work (GEMM operation: C = H * S)
      call compute_timer%start()
      call pic_gemm(H,S,C)
      call compute_timer%stop()
      elapsed_time = compute_timer%get_elapsed_time()

      ! Collapse H and S into 1D arrays
      allocate(H_flat(dims*dims), S_flat(dims*dims))
      H_flat = reshape(H, [dims*dims])
      S_flat = reshape(S, [dims*dims])

      ! Compute dot product of H and S
      dot_result = pic_dot(H_flat, S_flat)

      ! With 1e-3 scaling on H and S, dot_result ~ dims² * (fragment_idx * 1e-3)²
      ! For dims=30, fragment_idx=1: dot_result ~ 900 * 1e-6 = 9e-4
      ! This gives us a contribution on the order of 1e-4 to water_energy
      block 
      real(dp) :: division 
      if(fragment_size > 2) then 
        division = 1e5_dp
      else 
        division = 1e3_dp
      end if
      water_energy = water_energy - dot_result/division
      end block

      ! Collapse C into 1D array for sending
      allocate(C_flat(dims*dims))
      C_flat = reshape(C, [dims*dims])

      deallocate (H, S, C, H_flat, S_flat)
   end subroutine process_chemistry_fragment

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

      allocate(sum_by_level(max_level))
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

      print *, "MBE Energy breakdown:"
      do body_level = 1, max_level
         if (abs(sum_by_level(body_level)) > 1e-15_dp) then
            print '(a,i0,a,f20.10)', "  ", body_level, "-body:  ", sum_by_level(body_level)
         end if
      end do
      print '(a,f20.10)', "  Total:   ", total_energy

      deallocate(sum_by_level)

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

      allocate(indices(subset_size))
      allocate(subset(subset_size))

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

      deallocate(indices, subset)

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
      print *, "ERROR: Could not find fragment with monomers:", target_monomers
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
      print *, "Fragment breakdown:"
      do i = 1, max_level
         if (n_mers(i) > 0) then
            print '(a,i0,a,i0,a,f12.3,a)', "  ", i, "-mers:  ", n_mers(i), &
               " (", mer_flops(i)/1.0e9_dp, " GFLOP)"
         end if
      end do
      print '(a,i0,a,f12.3,a)', "  Total:    ", fragment_count, &
         " (", total_flops/1.0e9_dp, " GFLOP)"

      deallocate (n_mers, mer_flops)

   end subroutine calculate_exact_flops

   subroutine test_global_coordinator(world_comm, node_comm, total_fragments, polymers, max_level, &
                                       node_leader_ranks, num_nodes, matrix_size)
      type(comm_t), intent(in) :: world_comm, node_comm
      integer, intent(in) :: total_fragments, max_level, num_nodes, matrix_size
      integer, intent(in) :: polymers(:, :), node_leader_ranks(:)

      integer :: current_fragment, finished_nodes
      integer :: request_source, dummy_msg
      type(MPI_Status) :: status, local_status
      logical :: handling_local_workers
      logical :: has_pending

      ! For local workers
      integer :: local_finished_workers, fragment_size, local_dummy
      integer, allocatable :: fragment_indices(:)

      ! Storage for results
      real(dp), allocatable :: scalar_results(:)
      real(dp), allocatable :: matrix_results(:,:)
      real(dp), allocatable :: temp_matrix(:)
      integer :: max_matrix_size
      integer :: worker_fragment_map(node_comm%size())
      integer :: worker_source

      current_fragment = total_fragments
      finished_nodes = 0
      local_finished_workers = 0
      handling_local_workers = (node_comm%size() > 1)

      ! Allocate storage for results
      allocate(scalar_results(total_fragments))
      max_matrix_size = (max_level * matrix_size)**2
      allocate(matrix_results(max_matrix_size, total_fragments))
      scalar_results = 0.0_dp
      matrix_results = 0.0_dp
      worker_fragment_map = 0

      print *, "Global coordinator starting with", total_fragments, "fragments for", num_nodes, "nodes"

      do while (finished_nodes < num_nodes)

         ! Check for incoming results from local workers (tag 203 = scalar, tag 204 = matrix)
         if (handling_local_workers) then
            call iprobe(node_comm, MPI_ANY_SOURCE, 203, has_pending, local_status)
            if (has_pending) then
               worker_source = local_status%MPI_SOURCE
               ! Receive scalar result and store it using the fragment index for this worker
               call recv(node_comm, scalar_results(worker_fragment_map(worker_source)), worker_source, 203)
               ! Receive matrix result into temporary array, then copy to storage
               call recv(node_comm, temp_matrix, worker_source, 204, local_status)
               matrix_results(:, worker_fragment_map(worker_source)) = temp_matrix
               deallocate(temp_matrix)
            end if
         end if

         ! Remote node coordinator requests
         call iprobe(world_comm, MPI_ANY_SOURCE, 300, has_pending, status)
         if (has_pending) then
            call recv(world_comm, dummy_msg, status%MPI_SOURCE, 300)
            request_source = status%MPI_SOURCE

            if (current_fragment >= 1) then
               call send_fragment_to_node(world_comm, current_fragment, polymers, max_level, request_source)
               current_fragment = current_fragment - 1
            else
               call send(world_comm, -1, request_source, 302)
               finished_nodes = finished_nodes + 1
            end if
         end if

         ! Local workers (shared memory)
         if (handling_local_workers .and. local_finished_workers < node_comm%size() - 1) then
            call iprobe(node_comm, MPI_ANY_SOURCE, 200, has_pending, local_status)
            if (has_pending) then
               call recv(node_comm, local_dummy, local_status%MPI_SOURCE, 200)

               if (current_fragment >= 1) then
                  call send_fragment_to_worker(node_comm, current_fragment, polymers, max_level, &
                                                local_status%MPI_SOURCE, matrix_size)
                  ! Track which fragment was sent to which worker
                  worker_fragment_map(local_status%MPI_SOURCE) = current_fragment
                  current_fragment = current_fragment - 1
               else
                  call send(node_comm, -1, local_status%MPI_SOURCE, 202)
                  local_finished_workers = local_finished_workers + 1
               end if
            end if
         end if

         ! Finalize local worker completion
         if (handling_local_workers .and. local_finished_workers >= node_comm%size() - 1) then
            handling_local_workers = .false.
            if (num_nodes == 1) then
               finished_nodes = finished_nodes + 1
               print *, "Manually incremented finished_nodes for self"
            else
               finished_nodes = finished_nodes + 1
               print *, "Global coordinator finished local workers"
            end if
         end if
      end do

      print *, "Global coordinator finished all fragments"
      print *, "Sample scalar result (fragment 1):", scalar_results(1)
      block
      integer :: i
      real(dp) :: mbe_total_energy
      !do i = 1, size(scalar_results,1)
      !  print *, "Fragment", i, "energy:", scalar_results(i)
      !end do

      ! Compute the many-body expansion energy
      print *
      print *, "Computing Many-Body Expansion (MBE)..."
      call compute_mbe_energy(polymers, total_fragments, max_level, scalar_results, mbe_total_energy)

      end block
      print *, "Sample matrix result (fragment 1, first element):", matrix_results(1, 1)

      ! Cleanup
      deallocate(scalar_results, matrix_results)
   end subroutine test_global_coordinator

   subroutine send_fragment_to_node(world_comm, fragment_idx, polymers, max_level, dest_rank)
      type(comm_t), intent(in) :: world_comm
      integer, intent(in) :: fragment_idx, max_level, dest_rank
      integer, intent(in) :: polymers(:, :)
      integer :: fragment_size
      integer, allocatable :: fragment_indices(:)

      fragment_size = count(polymers(fragment_idx, :) > 0)
      allocate (fragment_indices(fragment_size))
      fragment_indices = polymers(fragment_idx, 1:fragment_size)

      call send(world_comm, fragment_idx, dest_rank, 301)
      call send(world_comm, fragment_size, dest_rank, 301)
      call send(world_comm, fragment_indices, dest_rank, 301)

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

      call send(node_comm, fragment_idx, dest_rank, 201)
      call send(node_comm, fragment_size, dest_rank, 201)
      call send(node_comm, fragment_indices, dest_rank, 201)
      call send(node_comm, matrix_size, dest_rank, 201)

      deallocate (fragment_indices)
   end subroutine send_fragment_to_worker

   subroutine test_node_coordinator(world_comm, node_comm, max_level, matrix_size)
      class(comm_t), intent(in) :: world_comm, node_comm
      integer(int32), intent(in) :: max_level, matrix_size

      integer(int32) :: fragment_idx, fragment_size, dummy_msg
      integer(int32) :: finished_workers
      integer(int32), allocatable :: fragment_indices(:)
      type(MPI_Status) :: status, global_status
      logical :: local_message_pending, more_fragments
      integer(int32) :: local_dummy

      finished_workers = 0
      more_fragments = .true.
      dummy_msg = 0

      do while (finished_workers < node_comm%size() - 1)
         call iprobe(node_comm, MPI_ANY_SOURCE, 200, local_message_pending, status)

         if (local_message_pending) then
            call recv(node_comm, local_dummy, MPI_ANY_SOURCE, 200)

            if (more_fragments) then
               call send(world_comm, dummy_msg, 0, 300)
               call recv(world_comm, fragment_idx, 0, MPI_ANY_TAG, global_status)

               if (global_status%MPI_TAG == 301) then
                  call recv(world_comm, fragment_size, 0, 301, global_status)
                  allocate (fragment_indices(fragment_size))
                  call recv(world_comm, fragment_indices, 0, 301, global_status)

                  call send(node_comm, fragment_idx, status%MPI_SOURCE, 201)
                  call send(node_comm, fragment_size, status%MPI_SOURCE, 201)
                  call send(node_comm, fragment_indices, status%MPI_SOURCE, 201)
                  call send(node_comm, matrix_size, status%MPI_SOURCE, 201)

                  deallocate (fragment_indices)
               else
                  call send(node_comm, -1, status%MPI_SOURCE, 202)
                  finished_workers = finished_workers + 1
                  more_fragments = .false.
               end if
            else
               call send(node_comm, -1, status%MPI_SOURCE, 202)
               finished_workers = finished_workers + 1
            end if
         end if
      end do
   end subroutine test_node_coordinator

   subroutine test_node_worker(world_comm, node_comm, max_level)
      class(comm_t), intent(in) :: world_comm, node_comm
      integer, intent(in) :: max_level

      integer(int32) :: fragment_idx, fragment_size, dummy_msg, matrix_size
      integer(int32), allocatable :: fragment_indices(:)
      real(dp) :: dot_result
      real(dp), allocatable :: C_flat(:)
      type(MPI_Status) :: status

      dummy_msg = 0

      do
         call send(node_comm, dummy_msg, 0, 200)
         call recv(node_comm, fragment_idx, 0, MPI_ANY_TAG, status)

         select case (status%MPI_TAG)
         case (201)
            call recv(node_comm, fragment_size, 0, 201, status)
            allocate (fragment_indices(fragment_size))
            call recv(node_comm, fragment_indices, 0, 201, status)
            call recv(node_comm, matrix_size, 0, 201, status)

            ! Process the chemistry fragment
            call process_chemistry_fragment(fragment_idx, fragment_indices, fragment_size, matrix_size, &
                                           dot_result, C_flat)

            ! Send results back to coordinator
            call send(node_comm, dot_result, 0, 203)
            call send(node_comm, C_flat, 0, 204)

            deallocate (fragment_indices, C_flat)
         case (202)
            exit
         end select
      end do
   end subroutine test_node_worker

end module pic_chemistry_algorithms
