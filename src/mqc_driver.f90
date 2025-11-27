module mqc_driver
   use pic_mpi_f08
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_chemistry_algorithms
   use mqc_mbe
   use mqc_physical_fragment
   use mqc_input_parser
   implicit none
   private

   public :: run_calculation

contains

   subroutine run_calculation(world_comm, node_comm, config, sys_geom)
      type(comm_t), intent(in) :: world_comm, node_comm
      type(input_config_t), intent(in) :: config
      type(system_geometry_t), intent(in) :: sys_geom

      ! Local variables
      integer :: world_rank, world_size, node_rank, node_size
      integer :: max_level, matrix_size
      integer :: total_fragments
      integer, allocatable :: polymers(:, :)
      integer :: num_nodes, i, j
      integer, allocatable :: node_leader_ranks(:)
      integer, allocatable :: monomers(:)
      integer :: n_expected_frags, n_rows
      integer :: global_node_rank
      integer, allocatable :: all_node_leader_ranks(:)

      ! Get MPI info
      world_rank = world_comm%rank()
      world_size = world_comm%size()
      node_rank = node_comm%rank()
      node_size = node_comm%size()

      ! Set max_level from config
      max_level = config%nlevel

      ! Set matrix_size based on atoms per monomer (natoms * 3 for gradient)
      matrix_size = sys_geom%atoms_per_monomer*3

      if (world_rank == 0) then
         call logger%info("============================================")
         call logger%info("Loaded geometry:")
         call logger%info("  Total monomers: "//to_char(sys_geom%n_monomers))
         call logger%info("  Atoms per monomer: "//to_char(sys_geom%atoms_per_monomer))
         call logger%info("  Total atoms: "//to_char(sys_geom%total_atoms))
         call logger%info("  Fragment level: "//to_char(max_level))
         call logger%info("  Matrix size (natoms*3): "//to_char(matrix_size))
         call logger%info("============================================")
      end if

      ! Check if this is an unfragmented calculation (nlevel=0)
      if (max_level == 0) then
         ! Validate that only a single rank is used for unfragmented calculation
         ! (parallelism comes from OpenMP threads, not MPI ranks)
         if (world_size > 1) then
            if (world_rank == 0) then
               call logger%error("")
               call logger%error("Unfragmented calculation (nlevel=0) requires exactly 1 MPI rank")
               call logger%error("  Parallelism is achieved through OpenMP threads, not MPI")
               call logger%error("  Current number of MPI ranks: "//to_char(world_size)//" (must be 1)")
               call logger%error("")
               call logger%error("Please run with a single MPI rank (e.g., mpirun -np 1 ...)")
               call logger%error("Use OMP_NUM_THREADS to control thread-level parallelism")
               call logger%error("")
            end if
            call abort_comm(world_comm, 1)
         end if

         if (world_rank == 0) then
            call logger%info("")
            call logger%info("nlevel=0 detected: Running unfragmented calculation")
            call logger%info("Parallelism provided by OpenMP threads")
            call logger%info("")
            call unfragmented_calculation(sys_geom)
         end if
         return
      end if

      ! Generate fragments (only rank 0 needs this for coordination)
      if (world_rank == 0) then
         ! Calculate expected number of fragments
         n_expected_frags = get_nfrags(sys_geom%n_monomers, max_level)
         n_rows = n_expected_frags

         ! Allocate monomer list and polymers array
         allocate (monomers(sys_geom%n_monomers))
         allocate (polymers(n_rows, max_level))
         polymers = 0

         ! Create monomer list [1, 2, 3, ..., n_monomers]
         call create_monomer_list(monomers)

         ! Generate all fragments (includes monomers in polymers array)
         total_fragments = 0

         ! First add monomers
         do i = 1, sys_geom%n_monomers
            total_fragments = total_fragments + 1
            polymers(total_fragments, 1) = i
         end do

         ! Then add n-mers for n >= 2
         call generate_fragment_list(monomers, max_level, polymers, total_fragments)

         deallocate (monomers)

         call logger%info("Generated fragments:")
         call logger%info("  Total fragments: "//to_char(total_fragments))
         call logger%info("  Max level: "//to_char(max_level))
      end if

      ! Broadcast total_fragments to all ranks
      call bcast(world_comm, total_fragments, 1, 0)

      ! Determine node leaders
      global_node_rank = -1
      if (node_comm%rank() == 0) global_node_rank = world_comm%rank()

      allocate (all_node_leader_ranks(world_comm%size()))
      call allgather(world_comm, global_node_rank, all_node_leader_ranks)

      num_nodes = count(all_node_leader_ranks /= -1)

      if (world_comm%rank() == 0) then
         call logger%info("Running with "//to_char(num_nodes)//" node(s)")
      end if

      allocate (node_leader_ranks(num_nodes))
      i = 0
      do j = 1, world_comm%size()
         if (all_node_leader_ranks(j) /= -1) then
            i = i + 1
            node_leader_ranks(i) = all_node_leader_ranks(j)
         end if
      end do
      deallocate (all_node_leader_ranks)

      ! Execute appropriate role
      if (world_comm%leader() .and. node_comm%leader()) then
         ! Global coordinator (rank 0, node leader on node 0)
         call logger%verbose("Rank 0: Acting as global coordinator")
         call global_coordinator(world_comm, node_comm, total_fragments, polymers, max_level, &
                                 node_leader_ranks, num_nodes, matrix_size)
      else if (node_comm%leader()) then
         ! Node coordinator (node leader on other nodes)
         call logger%verbose("Rank "//to_char(world_rank)//": Acting as node coordinator")
         call node_coordinator(world_comm, node_comm, max_level, matrix_size)
      else
         ! Worker
         call logger%verbose("Rank "//to_char(world_rank)//": Acting as worker")
         call node_worker(world_comm, node_comm, max_level, sys_geom)
      end if

      ! Cleanup
      if (world_rank == 0) then
         if (allocated(polymers)) deallocate (polymers)
         if (allocated(node_leader_ranks)) deallocate (node_leader_ranks)
      end if

   end subroutine run_calculation

end module mqc_driver
