program test_with_geometry
   use mpi_comm_simple
   use pic_chemistry_algorithms
   use pic_mbe
   use pic_physical_fragment
   use pic_input_parser
   implicit none

   type(comm_t) :: world_comm, node_comm
   integer :: world_rank, world_size, node_rank, node_size
   integer :: ierr

   ! Fragment data
   integer, parameter :: max_level = 3
   integer, parameter :: matrix_size = 10
   integer :: total_fragments
   integer, allocatable :: polymers(:, :)
   integer :: num_nodes, i
   integer, allocatable :: node_leader_ranks(:)
   integer, allocatable :: monomers(:)
   integer :: n_expected_frags, n_rows

   ! Geometry data
   type(input_config_t) :: config
   type(system_geometry_t) :: sys_geom
   integer :: stat
   character(len=:), allocatable :: errmsg
   character(len=256) :: input_file

   ! Initialize MPI
   call mpi_initialize()

   ! Create communicators
   world_comm = comm_world()
   node_comm = world_comm%split()

   world_rank = world_comm%rank()
   world_size = world_comm%size()
   node_rank = node_comm%rank()
   node_size = node_comm%size()

   ! All ranks read input file and load geometry
   input_file = "test.inp"

   call read_input_file(input_file, config, stat, errmsg)
   if (stat /= 0) then
      if (world_rank == 0) then
         print *, "ERROR reading input file: ", errmsg
      end if
      call abort_comm(world_comm, 1)
   end if

   call initialize_system_geometry(config%geom_file, config%monomer_file, sys_geom, stat, errmsg)
   if (stat /= 0) then
      if (world_rank == 0) then
         print *, "ERROR loading geometry: ", errmsg
      end if
      call abort_comm(world_comm, 1)
   end if

   if (world_rank == 0) then
      print *, "============================================"
      print *, "Loaded geometry:"
      print '(a,i0)', "  Total monomers: ", sys_geom%n_monomers
      print '(a,i0)', "  Atoms per monomer: ", sys_geom%atoms_per_monomer
      print '(a,i0)', "  Total atoms: ", sys_geom%total_atoms
      print *, "============================================"
   end if

   ! Generate fragments (only rank 0 needs this for coordination)
   if (world_rank == 0) then
      ! Calculate expected number of fragments
      n_expected_frags = get_nfrags(sys_geom%n_monomers, max_level)
      n_rows = n_expected_frags

      ! Allocate monomer list and polymers array
      allocate(monomers(sys_geom%n_monomers))
      allocate(polymers(n_rows, max_level))
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

      deallocate(monomers)

      print *, "Generated fragments:"
      print '(a,i0)', "  Total fragments: ", total_fragments
      print '(a,i0)', "  Max level: ", max_level
   end if

   ! Broadcast total_fragments to all ranks
   call bcast(world_comm, total_fragments, 1, 0)

   ! Determine number of nodes
   num_nodes = 1
   do i = 0, world_size - 1
      if (i == 0) cycle
      world_comm = comm_world()
      node_comm = world_comm%split()
      if (node_comm%rank() == 0) then
         num_nodes = num_nodes + 1
      end if
   end do

   ! Simplified: assume single node for this test
   num_nodes = 1
   allocate (node_leader_ranks(num_nodes))
   node_leader_ranks(1) = 0

   ! Reset communicators
   call world_comm%finalize()
   call node_comm%finalize()
   world_comm = comm_world()
   node_comm = world_comm%split()

   ! Run the test based on role
   if (world_rank == 0 .and. node_rank == 0) then
      ! Global coordinator
      print *, "Rank 0: Acting as global coordinator"
      call test_global_coordinator(world_comm, node_comm, total_fragments, polymers, max_level, &
                                    node_leader_ranks, num_nodes, matrix_size)
   else if (node_rank == 0) then
      ! Node coordinator (not used in single-node test)
      print *, "This should not happen in single-node test"
   else
      ! Worker
      print '(a,i0,a)', "Rank ", world_rank, ": Acting as worker"
      call test_node_worker(world_comm, node_comm, max_level, sys_geom)
   end if

   ! Cleanup
   if (world_rank == 0) then
      if (allocated(polymers)) deallocate (polymers)
      if (allocated(node_leader_ranks)) deallocate (node_leader_ranks)
   end if
   call config%destroy()
   call sys_geom%destroy()
   call world_comm%finalize()
   call node_comm%finalize()

   call mpi_finalize_wrapper()

end program test_with_geometry
