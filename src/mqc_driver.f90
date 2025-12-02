!! Main calculation driver module for metalquicha
module mqc_driver
   !! Handles both fragmented (many-body expansion) and unfragmented calculations
   !! with MPI parallelization and node-based work distribution.
   use pic_mpi_lib, only: comm_t, abort_comm, bcast, allgather
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_mbe, only: global_coordinator, node_coordinator, node_worker, unfragmented_calculation
   use mqc_frag_utils, only: get_nfrags, create_monomer_list, generate_fragment_list
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_input_parser, only: input_config_t
   implicit none
   private

   public :: run_calculation  !! Main entry point for all calculations

contains

   subroutine run_calculation(world_comm, node_comm, config, sys_geom)
      !! Main calculation dispatcher - routes to fragmented or unfragmented calculation
      !!
      !! Determines calculation type based on nlevel and dispatches to appropriate
      !! calculation routine with proper MPI setup and validation.
      type(comm_t), intent(in) :: world_comm  !! Global MPI communicator
      type(comm_t), intent(in) :: node_comm   !! Node-local MPI communicator
      type(input_config_t), intent(in) :: config  !! Parsed input configuration
      type(system_geometry_t), intent(in) :: sys_geom  !! System geometry and fragment info

      ! Local variables
      integer :: max_level   !! Maximum fragment level (nlevel from config)
      integer :: matrix_size  !! Size of gradient matrix (natoms*3), tmp
      integer :: total_fragments  !! Total number of fragments generated
      integer, allocatable :: polymers(:, :)  !! Fragment indices array
      integer :: num_nodes   !! Number of compute nodes
      integer :: i, j        !! Loop counters
      integer, allocatable :: node_leader_ranks(:)  !! Ranks of node leaders
      integer, allocatable :: monomers(:)     !! Monomer indices for fragment generation
      integer :: n_expected_frags  !! Expected number of fragments
      integer :: n_rows      !! Number of rows for polymers array
      integer :: global_node_rank  !! Global rank if node leader, -1 otherwise
      integer, allocatable :: all_node_leader_ranks(:)  !! All node leader ranks

      ! Set max_level from config
      max_level = config%nlevel

      ! Set matrix_size based on atoms per monomer (natoms * 3 for gradient)
      matrix_size = sys_geom%atoms_per_monomer*3

      if (world_comm%rank() == 0) then
         call logger%info("============================================")
         call logger%info("Loaded geometry:")
         call logger%info("  Total monomers: "//to_char(sys_geom%n_monomers))
         call logger%info("  Atoms per monomer: "//to_char(sys_geom%atoms_per_monomer))
         call logger%info("  Total atoms: "//to_char(sys_geom%total_atoms))
         call logger%info("  Fragment level: "//to_char(max_level))
         call logger%info("  Matrix size (natoms*3): "//to_char(matrix_size))
         call logger%info("============================================")
      end if

      ! Choose calculation type based on nlevel
      if (max_level == 0) then
         call run_unfragmented_calculation(world_comm, sys_geom, config%method)
      else
         call run_fragmented_calculation(world_comm, node_comm, config%method, sys_geom, max_level, matrix_size)
      end if

   end subroutine run_calculation

   subroutine run_unfragmented_calculation(world_comm, sys_geom, method)
      !! Handle unfragmented calculation (nlevel=0)
      !!
      !! Validates single MPI rank requirement and runs direct calculation
      !! on entire system using OpenMP for parallelization.
      type(comm_t), intent(in) :: world_comm  !! Global MPI communicator
      type(system_geometry_t), intent(in) :: sys_geom  !! Complete system geometry
      character(len=*), intent(in) :: method  !! Quantum chemistry method (gfn1/gfn2)

      ! Validate that only a single rank is used for unfragmented calculation
      ! (parallelism comes from OpenMP threads, not MPI ranks)
      if (world_comm%size() > 1) then
         if (world_comm%rank() == 0) then
            call logger%error("")
            call logger%error("Unfragmented calculation (nlevel=0) requires exactly 1 MPI rank")
            call logger%error("  Parallelism is achieved through OpenMP threads, not MPI")
            call logger%error("  Current number of MPI ranks: "//to_char(world_comm%size())//" (must be 1)")
            call logger%error("")
            call logger%error("Please run with a single MPI rank (e.g., mpirun -np 1 ...)")
            call logger%error("Use OMP_NUM_THREADS to control thread-level parallelism")
            call logger%error("")
         end if
         call abort_comm(world_comm, 1)
      end if

      if (world_comm%rank() == 0) then
         call logger%info("")
         call logger%info("nlevel=0 detected: Running unfragmented calculation")
         call logger%info("Parallelism provided by OpenMP threads")
         call logger%info("")
         call unfragmented_calculation(sys_geom, method)
      end if

   end subroutine run_unfragmented_calculation

   subroutine run_fragmented_calculation(world_comm, node_comm, method, sys_geom, max_level, matrix_size)
      !! Handle fragmented calculation (nlevel > 0)
      !!
      !! Generates fragments, distributes work across MPI processes organized in nodes,
      !! and coordinates many-body expansion calculation using hierarchical parallelism.
      type(comm_t), intent(in) :: world_comm  !! Global MPI communicator
      type(comm_t), intent(in) :: node_comm   !! Node-local MPI communicator
      character(len=*), intent(in) :: method  !! Quantum chemistry method (gfn1/gfn2)
      type(system_geometry_t), intent(in) :: sys_geom  !! System geometry and fragment info
      integer, intent(in) :: max_level    !! Maximum fragment level for MBE
      integer, intent(in) :: matrix_size  !! Size of gradient matrix (natoms*3)

      integer :: total_fragments  !! Total number of fragments generated
      integer, allocatable :: polymers(:, :)  !! Fragment composition array (fragment, monomer_indices)
      integer :: num_nodes   !! Number of compute nodes
      integer :: i, j        !! Loop counters
      integer, allocatable :: node_leader_ranks(:)  !! Ranks of processes that lead each node
      integer, allocatable :: monomers(:)     !! Temporary monomer list for fragment generation
      integer :: n_expected_frags  !! Expected number of fragments based on combinatorics
      integer :: n_rows      !! Number of rows needed for polymers array
      integer :: global_node_rank  !! Global rank if this process leads a node, -1 otherwise
      integer, allocatable :: all_node_leader_ranks(:)  !! Node leader status for all ranks

      ! Generate fragments (only rank 0 needs this for coordination)
      if (world_comm%size() == 1) then
         if (world_comm%rank() == 0) then
            call logger%error("Fragmented calculation cannot be run with one rank")
         end if
         call abort_comm(world_comm, 1)
      end if
      if (world_comm%rank() == 0) then
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
         call logger%verbose("Rank "//to_char(world_comm%rank())//": Acting as node coordinator")
         call node_coordinator(world_comm, node_comm, max_level, matrix_size)
      else
         ! Worker
         call logger%verbose("Rank "//to_char(world_comm%rank())//": Acting as worker")
         call node_worker(world_comm, node_comm, max_level, sys_geom, method)
      end if

      ! Cleanup
      if (world_comm%rank() == 0) then
         if (allocated(polymers)) deallocate (polymers)
         if (allocated(node_leader_ranks)) deallocate (node_leader_ranks)
      end if

   end subroutine run_fragmented_calculation

end module mqc_driver
