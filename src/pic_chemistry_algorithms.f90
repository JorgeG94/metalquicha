module pic_chemistry_algorithms
   use mpi_f08
   use pic_types
   use pic_timer
   use mpi_comm_simple
   use pic_fragment, only: pic_fragment_block
   use pic_blas_interfaces, only: pic_gemm
   implicit none

contains

   subroutine process_chemistry_fragment(fragment_idx, fragment_indices, fragment_size, matrix_size)
      integer, intent(in) :: fragment_idx, fragment_size, matrix_size
      integer, intent(in) :: fragment_indices(fragment_size)
      real(dp), allocatable :: H(:, :), S(:, :), C(:, :)
      integer :: i, j, dims
      type(timer_type) :: compute_timer
      real(dp) :: elapsed_time
      real(dp), parameter :: alpha = 17.0_dp

      dims = fragment_size*matrix_size

      ! Allocate matrices for Hamiltonian, Overlap, and Coefficients
      allocate (H(dims, dims), S(dims, dims), C(dims, dims))

      ! Initialize matrices (simulating chemistry computation)
      do concurrent(j=1:dims, i=1:dims)
         H(i, j) = real(fragment_idx, dp)  ! Hamiltonian matrix
         S(i, j) = real(fragment_idx, dp)  ! Overlap matrix
         C(i, j) = 0.0_dp                   ! Coefficient matrix
      end do

      ! Simulate chemistry work (JAXPY-like operation: C = alpha * H + S)
      call compute_timer%start()
      call pic_gemm(H,S,C)
      call compute_timer%stop()
      elapsed_time = compute_timer%get_elapsed_time()

      deallocate (H, S, C)
   end subroutine process_chemistry_fragment

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

      current_fragment = total_fragments
      finished_nodes = 0
      local_finished_workers = 0
      handling_local_workers = (node_comm%size() > 1)

      print *, "Global coordinator starting with", total_fragments, "fragments for", num_nodes, "nodes"

      do while (finished_nodes < num_nodes)

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
            call process_chemistry_fragment(fragment_idx, fragment_indices, fragment_size, matrix_size)

            deallocate (fragment_indices)
         case (202)
            exit
         end select
      end do
   end subroutine test_node_worker

end module pic_chemistry_algorithms
