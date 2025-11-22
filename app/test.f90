program hierarchical_mpi_mbe
   use mpi_f08
   use mpi_comm_simple
   use omp_lib
   !use pic_blas_interfaces, only: pic_gemm
   use pic_timer, only: timer_type
   use pic_types, only: dp, default_int
   use pic_io, only: to_char
   use pic_mbe, only: get_nfrags, create_monomer_list, generate_fragment_list
   use pic_fragment, only: pic_fragment_block, count_nonzeros
   use pic_chemistry_algorithms
   implicit none

   ! Fragment generation parameters
   integer(default_int) :: n_monomers
   integer(default_int) :: max_level
   integer(default_int) :: n
   integer :: narg
   character(len=32) :: arg

   ! MPI wrappers
   type(comm_t) :: world_comm, node_comm
   character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
   integer :: hostname_len

   ! Timing
   type(timer_type) :: timer
   real(dp) :: elapsed_time, flops

   ! Fragment data structures
   integer(default_int), allocatable :: monomers(:)
   integer(default_int), allocatable :: polymers(:, :)
   integer(default_int) :: n_expected_fragments, fragment_count
   type(pic_fragment_block), allocatable :: fragments(:)

   ! Node leader detection
   integer :: global_node_rank
   integer, allocatable :: all_node_leader_ranks(:), node_leader_ranks(:)
   integer :: i, j, num_nodes

   !==============================
   ! MPI Initialization
   !==============================
   call MPI_Init()
   world_comm = comm_world()
   node_comm = world_comm%split()   ! shared memory communicator
   narg = command_argument_count()
   n_monomers = 10
   max_level = 3
   n = 25  ! monomer matrix size

   if (narg >= 1) then
      call get_command_argument(1, arg)
      read (arg, *) n_monomers
   end if
   if (narg >= 2) then
      call get_command_argument(2, arg)
      read (arg, *) n
   end if
   if (narg >= 3) then
      call get_command_argument(3, arg)
      read (arg, *) max_level
   end if

   if (world_comm%size() < 2) then
      if (world_comm%leader()) print *, "This program requires at least 2 processes."
      call abort_comm(world_comm, 1)
   end if

   if (world_comm%leader()) then
      call timer%start()
   end if

   !==============================
   ! Determine node leaders
   !==============================
   global_node_rank = -1
   if (node_comm%leader()) global_node_rank = world_comm%rank()

   allocate (all_node_leader_ranks(world_comm%size()))
   call allgather(world_comm, global_node_rank, all_node_leader_ranks)

   num_nodes = count(all_node_leader_ranks /= -1)
   !print *, "RUNNING USING ", num_nodes, " NODES"
! After determining num_nodes and node sizes
   if (num_nodes > 1) then
      ! Multi-node: warn if any node has only a coordinator
      if (world_comm%size() < 3) then
         if (world_comm%leader()) print *, "This program requires at least 3 processes."
         call abort_comm(world_comm, 1)
      end if
   end if
   allocate (node_leader_ranks(num_nodes))
   i = 0
   do j = 1, world_comm%size()
      if (all_node_leader_ranks(j) /= -1) then
         i = i + 1
         node_leader_ranks(i) = all_node_leader_ranks(j)
      end if
   end do

   !==============================
   ! Generate fragments (on rank 0)
   !==============================
   if (world_comm%leader()) then
      n_expected_fragments = get_nfrags(n_monomers, max_level)
      allocate (monomers(n_monomers))
      allocate (polymers(n_expected_fragments, max_level))
      allocate (fragments(n_expected_fragments))

      print *, "Expected fragments = ", n_expected_fragments

      polymers = 0_default_int
      call create_monomer_list(monomers)
      polymers(:n_monomers, 1) = monomers(:n_monomers)
      fragment_count = n_monomers

      call generate_fragment_list(monomers, max_level, polymers, fragment_count)

      if (fragment_count /= n_expected_fragments) then
         print *, "Error: fragment_count does not match expected fragments!"
         call abort_comm(world_comm, 1)
      end if

      print *, "Generated ", fragment_count, " fragments for distribution"
   end if

   if (node_comm%leader() .eqv. .false.) then
      call get_processor_name(hostname, hostname_len)
      !print *, "Node:", trim(hostname), " node_comm rank:", node_comm%rank() - 1, " world rank:", world_comm%rank()
      !!call omp_set_default_device(node_comm%rank() - 1)
      !!!$acc set device_num(node_comm%rank() - 1)
   end if

   !==============================
   ! Role assignment
   !==============================
   if (world_comm%leader() .and. node_comm%leader()) then
      call test_global_coordinator(world_comm, node_comm, fragment_count, polymers, max_level, &
                                    node_leader_ranks, num_nodes, n)
   else if (node_comm%leader()) then
      call test_node_coordinator(world_comm, node_comm, max_level, n)
   else
      call test_node_worker(world_comm, node_comm, max_level)
   end if

   !==============================
   ! Final timing and flops
   !==============================
   call world_comm%barrier()
   if (world_comm%leader()) then
      call timer%stop()
      elapsed_time = timer%get_elapsed_time()
      flops = 0.0_dp
      print *, "N = ", n
      call calculate_exact_flops(polymers, fragment_count, max_level, n, flops)

      print *, "Total elapsed time for all fragments:", elapsed_time, "seconds"
      print *, "Total flops ", flops/1.0e9_dp, " GFLOPs"
      print *, "Estimated flop rate: ", flops/elapsed_time/1.0e9_dp, " GFLOP/s"
   end if

   !==============================
   ! Finalization
   !==============================
   call node_comm%finalize()
   call world_comm%finalize()
   call MPI_Finalize()

end program hierarchical_mpi_mbe
