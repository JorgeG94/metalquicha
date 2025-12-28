!! Main program for metalquicha quantum chemistry calculations
!!
!! Input format: .mqc (section-based format)
!!
!! Usage: metalquicha input_file.mqc
program main
   !! Orchestrates MPI initialization, input parsing, geometry loading,
   !! and dispatches to appropriate calculation routines (fragmented or unfragmented).
   use pic_logger, only: logger => global_logger, info_level
   use pic_io, only: to_char
   use pic_mpi_lib, only: pic_mpi_init, comm_world, comm_t, abort_comm, pic_mpi_finalize
   use mqc_driver, only: run_calculation
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_config_parser, only: mqc_config_t, read_mqc_file
   use mqc_config_adapter, only: driver_config_t, config_to_driver, config_to_system_geometry, get_logger_level
   use mqc_logo, only: print_logo
   use pic_timer, only: timer_type
   implicit none

   type(timer_type) :: my_timer      !! Execution timing
   type(comm_t) :: world_comm        !! Global MPI communicator
   type(comm_t) :: node_comm         !! Node-local MPI communicator
   type(driver_config_t) :: config   !! Driver configuration
   type(mqc_config_t) :: mqc_config  !! Parsed .mqc file
   type(system_geometry_t) :: sys_geom  !! Loaded molecular system
   integer :: stat                  !! Status code for error handling
   character(len=:), allocatable :: errmsg  !! Error messages
   character(len=256) :: input_file  !! Input file name

   ! Initialize MPI
   ! pic-mpi will call mpi_init_thread when needed
   call pic_mpi_init()

   ! Create communicators
   world_comm = comm_world()
   node_comm = world_comm%split()

   if (world_comm%rank() == 0) then
      call print_logo()
      call my_timer%start()
   end if

   ! Parse command line arguments
   if (command_argument_count() == 0) then
      if (world_comm%rank() == 0) then
         call logger%error("No input file specified. Usage: metalquicha input_file.mqc")
      end if
      call abort_comm(world_comm, 1)
   else if (command_argument_count() == 1) then
      call get_command_argument(1, input_file, status=stat)
      if (stat /= 0) then
         if (world_comm%rank() == 0) then
            call logger%error("Error reading command line argument")
         end if
         call abort_comm(world_comm, 1)
      end if
      input_file = trim(input_file)

      ! Validate file extension
      if (.not. ends_with(input_file, '.mqc')) then
         if (world_comm%rank() == 0) then
            call logger%error("Invalid input file extension. Expected .mqc")
         end if
         call abort_comm(world_comm, 1)
      end if
   else
      if (world_comm%rank() == 0) then
         call logger%error("Too many arguments. Usage: metalquicha [input_file.mqc]")
      end if
      call abort_comm(world_comm, 1)
   end if

   ! Parse .mqc input file
   if (world_comm%rank() == 0) then
      call logger%info("Reading input file: "//trim(input_file))
   end if

   call read_mqc_file(input_file, mqc_config, stat, errmsg)
   if (stat /= 0) then
      if (world_comm%rank() == 0) then
         call logger%error("Error reading .mqc file: "//errmsg)
      end if
      call abort_comm(world_comm, 1)
   end if

   ! Convert to driver configuration
   call config_to_driver(mqc_config, config)

   ! Configure logger
   if (allocated(mqc_config%log_level)) then
      call logger%configure(get_logger_level(mqc_config%log_level))
      if (world_comm%rank() == 0) then
         call logger%info("Logger verbosity set to: "//trim(mqc_config%log_level))
      end if
   end if

   ! Convert geometry to system_geometry_t
   call config_to_system_geometry(mqc_config, sys_geom, stat, errmsg)
   if (stat /= 0) then
      if (world_comm%rank() == 0) then
         call logger%error("Error converting geometry: "//errmsg)
      end if
      call abort_comm(world_comm, 1)
   end if

   call run_calculation(world_comm, node_comm, config, sys_geom, mqc_config%bonds)

   if (world_comm%rank() == 0) then
      call my_timer%stop()
      call logger%info("Total processing time: "//to_char(my_timer%get_elapsed_time())//" s")
   end if

   call mqc_config%destroy()
   call sys_geom%destroy()
   call world_comm%finalize()
   call node_comm%finalize()
   call pic_mpi_finalize()

contains

   logical function ends_with(str, suffix)
      !! Check if string ends with suffix
      character(len=*), intent(in) :: str, suffix
      integer :: str_len, suffix_len

      str_len = len_trim(str)
      suffix_len = len_trim(suffix)

      if (suffix_len > str_len) then
         ends_with = .false.
      else
         ends_with = (str(str_len - suffix_len + 1:str_len) == suffix)
      end if
   end function ends_with

end program main
