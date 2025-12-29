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

   ! Configure logger
   if (allocated(mqc_config%log_level)) then
      call logger%configure(get_logger_level(mqc_config%log_level))
      if (world_comm%rank() == 0) then
         call logger%info("Logger verbosity set to: "//trim(mqc_config%log_level))
      end if
   end if

   ! Handle single vs multiple molecules
   if (mqc_config%nmol == 0) then
      ! Single molecule mode (backward compatible)
      call config_to_driver(mqc_config, config)
      call config_to_system_geometry(mqc_config, sys_geom, stat, errmsg)
      if (stat /= 0) then
         if (world_comm%rank() == 0) then
            call logger%error("Error converting geometry: "//errmsg)
         end if
         call abort_comm(world_comm, 1)
      end if

      call run_calculation(world_comm, node_comm, config, sys_geom, mqc_config%bonds)
      call sys_geom%destroy()
   else
      ! Multi-molecule mode: loop over all molecules
      call run_multi_molecule_calculations(world_comm, node_comm, mqc_config)
   end if

   if (world_comm%rank() == 0) then
      call my_timer%stop()
      call logger%info("Total processing time: "//to_char(my_timer%get_elapsed_time())//" s")
   end if

   call mqc_config%destroy()
   call world_comm%finalize()
   call node_comm%finalize()
   call pic_mpi_finalize()

contains

   subroutine run_multi_molecule_calculations(world_comm, node_comm, mqc_config)
      !! Run calculations for multiple molecules with MPI parallelization
      !! Each molecule is independent, so assign one molecule per rank
      type(comm_t), intent(in) :: world_comm
      type(comm_t), intent(in) :: node_comm
      type(mqc_config_t), intent(in) :: mqc_config

      type(driver_config_t) :: config
      type(system_geometry_t) :: sys_geom
      type(comm_t) :: mol_comm, mol_node_comm
      integer :: imol, stat, my_rank, num_ranks, color
      integer :: molecules_processed
      character(len=:), allocatable :: errmsg
      character(len=:), allocatable :: mol_name
      logical :: has_fragmented_molecules

      my_rank = world_comm%rank()
      num_ranks = world_comm%size()

      ! Check if any molecules have fragments (nlevel > 0)
      has_fragmented_molecules = .false.
      do imol = 1, mqc_config%nmol
         if (mqc_config%molecules(imol)%nfrag > 0) then
            has_fragmented_molecules = .true.
            exit
         end if
      end do

      if (my_rank == 0) then
         call logger%info(" ")
         call logger%info("============================================")
         call logger%info("Multi-molecule mode: "//to_char(mqc_config%nmol)//" molecules")
         call logger%info("MPI ranks: "//to_char(num_ranks))
         if (has_fragmented_molecules) then
            call logger%info("Mode: Sequential execution (fragmented molecules detected)")
            call logger%info("  Each molecule will use all "//to_char(num_ranks)//" rank(s) for its calculation")
         else if (num_ranks == 1) then
            call logger%info("Mode: Sequential execution (single rank)")
         else if (num_ranks > mqc_config%nmol) then
            call logger%info("Mode: Parallel execution (one molecule per rank)")
            call logger%info("Note: More ranks than molecules - ranks "//to_char(mqc_config%nmol)// &
                             " to "//to_char(num_ranks - 1)//" will be idle")
         else
            call logger%info("Mode: Parallel execution (one molecule per rank)")
         end if
         call logger%info("============================================")
         call logger%info(" ")
      end if

      ! Determine execution mode:
      ! 1. Sequential: Single rank OR fragmented molecules (each molecule needs all ranks)
      ! 2. Parallel: Multiple ranks AND unfragmented molecules (distribute molecules across ranks)
      molecules_processed = 0

      if (num_ranks == 1 .or. has_fragmented_molecules) then
         ! Sequential mode: process all molecules one after another
         ! Each molecule uses all available ranks for its calculation
         do imol = 1, mqc_config%nmol
            ! Determine molecule name for logging
            if (allocated(mqc_config%molecules(imol)%name)) then
               mol_name = mqc_config%molecules(imol)%name
            else
               mol_name = "molecule_"//to_char(imol)
            end if

            if (my_rank == 0) then
               call logger%info(" ")
               call logger%info("--------------------------------------------")
               call logger%info("Processing molecule "//to_char(imol)//"/"//to_char(mqc_config%nmol)//": "//mol_name)
               call logger%info("--------------------------------------------")
            end if

            ! Convert to driver configuration for this molecule
            call config_to_driver(mqc_config, config, molecule_index=imol)

            ! Convert geometry for this molecule
            call config_to_system_geometry(mqc_config, sys_geom, stat, errmsg, molecule_index=imol)
            if (stat /= 0) then
               if (my_rank == 0) then
                  call logger%error("Error converting geometry for "//mol_name//": "//errmsg)
               end if
               call abort_comm(world_comm, 1)
            end if

            ! Run calculation for this molecule
            call run_calculation(world_comm, node_comm, config, sys_geom, mqc_config%molecules(imol)%bonds)

            ! Clean up for this molecule
            call sys_geom%destroy()

            if (my_rank == 0) then
               call logger%info("Completed molecule "//to_char(imol)//"/"//to_char(mqc_config%nmol)//": "//mol_name)
            end if
            molecules_processed = molecules_processed + 1
         end do
      else
         ! Multiple ranks: distribute molecules across ranks (one per rank)
         if (my_rank < mqc_config%nmol) then
            imol = my_rank + 1
            ! This rank has a molecule to process
            ! Determine molecule name for logging
            if (allocated(mqc_config%molecules(imol)%name)) then
               mol_name = mqc_config%molecules(imol)%name
            else
               mol_name = "molecule_"//to_char(imol)
            end if

            call logger%info(" ")
            call logger%info("--------------------------------------------")
            call logger%info("Rank "//to_char(my_rank)//": Processing molecule "//to_char(imol)// &
                             "/"//to_char(mqc_config%nmol)//": "//mol_name)
            call logger%info("--------------------------------------------")

            ! Convert to driver configuration for this molecule
            call config_to_driver(mqc_config, config, molecule_index=imol)

            ! Convert geometry for this molecule
            call config_to_system_geometry(mqc_config, sys_geom, stat, errmsg, molecule_index=imol)
            if (stat /= 0) then
               call logger%error("Rank "//to_char(my_rank)//": Error converting geometry for "//mol_name//": "//errmsg)
               call abort_comm(world_comm, 1)
            end if

            ! Run calculation for this molecule
            call run_calculation(world_comm, node_comm, config, sys_geom, mqc_config%molecules(imol)%bonds)

            ! Clean up for this molecule
            call sys_geom%destroy()

            call logger%info("Rank "//to_char(my_rank)//": Completed molecule "//to_char(imol)// &
                             "/"//to_char(mqc_config%nmol)//": "//mol_name)

            molecules_processed = 1
         else
            ! Idle rank - no molecule assigned
            call logger%verbose("Rank "//to_char(my_rank)//": No molecule assigned (idle)")
         end if
      end if

      ! Synchronize all ranks
      call world_comm%barrier()

      if (my_rank == 0) then
         call logger%info(" ")
         call logger%info("============================================")
         call logger%info("All "//to_char(mqc_config%nmol)//" molecules completed")
         if (has_fragmented_molecules) then
            call logger%info("Execution: Sequential (each molecule used all ranks)")
         else if (num_ranks == 1) then
            call logger%info("Execution: Sequential (single rank)")
         else if (num_ranks > mqc_config%nmol) then
           call logger%info("Execution: Parallel (active ranks: "//to_char(mqc_config%nmol)//"/"//to_char(num_ranks)//")")
         else
            call logger%info("Execution: Parallel (all ranks active)")
         end if
         call logger%info("============================================")
      end if

   end subroutine run_multi_molecule_calculations

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
