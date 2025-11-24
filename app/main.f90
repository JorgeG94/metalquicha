program main
   use mpi_comm_simple
   use pic_driver
   use pic_physical_fragment
   use pic_input_parser
   use pic_timer
   implicit none

   type(timer_type) :: my_timer
   type(comm_t) :: world_comm, node_comm
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

   ! Start timer on rank 0
   if (world_comm%rank() == 0) then
      call my_timer%start()
   end if

   ! Read input file
   input_file = "test.inp"
   call read_input_file(input_file, config, stat, errmsg)
   if (stat /= 0) then
      if (world_comm%rank() == 0) then
         print *, "ERROR reading input file: ", errmsg
      end if
      call abort_comm(world_comm, 1)
   end if

   ! Initialize system geometry
   call initialize_system_geometry(config%geom_file, config%monomer_file, sys_geom, stat, errmsg)
   if (stat /= 0) then
      if (world_comm%rank() == 0) then
         print *, "ERROR loading geometry: ", errmsg
      end if
      call abort_comm(world_comm, 1)
   end if

   ! Run the calculation
   call run_calculation(world_comm, node_comm, config, sys_geom)

   ! Stop timer and report on rank 0
   if (world_comm%rank() == 0) then
      call my_timer%stop()
      print *, "Total processing time ", my_timer%get_elapsed_time(), " s"
   end if

   ! Cleanup and finalize
   call config%destroy()
   call sys_geom%destroy()
   call world_comm%finalize()
   call node_comm%finalize()
   call mpi_finalize_wrapper()

end program main
