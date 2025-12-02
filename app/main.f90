!! This file simply contains the main program, look at it as the base of the calculation
!! and how everything sits together.
program main
   !! Orchestrates MPI initialization, input parsing, geometry loading,
   !! and dispatches to appropriate calculation routines (fragmented or unfragmented).
   use pic_logger, only: logger => global_logger, info_level
   use pic_io, only: to_char
   use pic_mpi_lib, only: pic_mpi_init, comm_world, comm_t, abort_comm, pic_mpi_finalize
   use mqc_driver, only: run_calculation
   use mqc_physical_fragment, only: initialize_system_geometry, system_geometry_t
   use mqc_input_parser, only: read_input_file, input_config_t
   use mqc_logo, only: print_logo
   use pic_timer, only: timer_type
   implicit none

   type(timer_type) :: my_timer     !! Execution timing
   type(comm_t) :: world_comm       !! Global MPI communicator
   type(comm_t) :: node_comm        !! Node-local MPI communicator
   type(input_config_t) :: config   !! Parsed input configuration
   type(system_geometry_t) :: sys_geom  !! Loaded molecular system
   integer :: stat                  !! Status code for error handling
   character(len=:), allocatable :: errmsg  !! Error messages
   character(len=256) :: input_file  !! Input file name

   ! Initialize MPI
   call pic_mpi_init()

   ! Create communicators
   world_comm = comm_world()
   node_comm = world_comm%split()

   ! Start timer on rank 0
   if (world_comm%rank() == 0) then
      call print_logo()
      call my_timer%start()
   end if

   ! Read input file
   input_file = "test.inp"
   call read_input_file(input_file, config, stat, errmsg)
   if (stat /= 0) then
      if (world_comm%rank() == 0) then
         call logger%error("Error reading input file "//errmsg)
      end if
      call abort_comm(world_comm, 1)
   end if

   ! Initialize system geometry
   call initialize_system_geometry(config%geom_file, config%monomer_file, sys_geom, stat, errmsg)
   if (stat /= 0) then
      if (world_comm%rank() == 0) then
         call logger%error("Error loeading geometry "//errmsg)
      end if
      call abort_comm(world_comm, 1)
   end if

   ! Run the calculation
   call run_calculation(world_comm, node_comm, config, sys_geom)

   ! Stop timer and report on rank 0
   if (world_comm%rank() == 0) then
      call my_timer%stop()
      call logger%info("Total processing time: "//to_char(my_timer%get_elapsed_time())//" s")
   end if

   ! Cleanup and finalize
   call config%destroy()
   call sys_geom%destroy()
   call world_comm%finalize()
   call node_comm%finalize()
   call pic_mpi_finalize()

end program main
