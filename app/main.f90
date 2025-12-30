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
   use mqc_output_filename, only: set_output_json_filename, set_molecule_suffix, get_output_json_filename
   use mqc_logo, only: print_logo
   use pic_timer, only: timer_type
   use mqc_error, only: error_t
   implicit none

   type(timer_type) :: my_timer      !! Execution timing
   type(comm_t) :: world_comm        !! Global MPI communicator
   type(comm_t) :: node_comm         !! Node-local MPI communicator
   type(driver_config_t) :: config   !! Driver configuration
   type(mqc_config_t) :: mqc_config  !! Parsed .mqc file
   type(system_geometry_t) :: sys_geom  !! Loaded molecular system
   type(error_t) :: error            !! Error handling
   integer :: stat                   !! Status code for file I/O
   character(len=:), allocatable :: errmsg  !! Error messages for file I/O
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

      call set_output_json_filename(input_file)
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

   call read_mqc_file(input_file, mqc_config, error)
   if (error%has_error()) then
      if (world_comm%rank() == 0) then
         call logger%error("Error reading .mqc file: "//error%get_message())
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
      call config_to_system_geometry(mqc_config, sys_geom, error)
      if (error%has_error()) then
         if (world_comm%rank() == 0) then
            call logger%error("Error converting geometry: "//error%get_message())
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
      type(error_t) :: error
      integer :: imol, my_rank, num_ranks, color
      integer :: molecules_processed
      character(len=:), allocatable :: mol_name
      logical :: has_fragmented_molecules
      character(len=256), allocatable :: individual_json_files(:)

      my_rank = world_comm%rank()
      num_ranks = world_comm%size()

      ! Allocate array to track individual JSON files for merging
      allocate (individual_json_files(mqc_config%nmol))

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
            call config_to_system_geometry(mqc_config, sys_geom, error, molecule_index=imol)
            if (error%has_error()) then
               if (my_rank == 0) then
                  call logger%error("Error converting geometry for "//mol_name//": "//error%get_message())
               end if
               call abort_comm(world_comm, 1)
            end if

            ! Set output filename suffix for this molecule
            call set_molecule_suffix("_"//trim(mol_name))

            ! Run calculation for this molecule
            call run_calculation(world_comm, node_comm, config, sys_geom, mqc_config%molecules(imol)%bonds)

            ! Track the JSON filename for later merging
            individual_json_files(imol) = get_output_json_filename()

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
            call config_to_system_geometry(mqc_config, sys_geom, error, molecule_index=imol)
            if (error%has_error()) then
     call logger%error("Rank "//to_char(my_rank)//": Error converting geometry for "//mol_name//": "//error%get_message())
               call abort_comm(world_comm, 1)
            end if

            ! Set output filename suffix for this molecule
            call set_molecule_suffix("_"//trim(mol_name))

            ! Run calculation for this molecule
            call run_calculation(world_comm, node_comm, config, sys_geom, mqc_config%molecules(imol)%bonds)

            ! Track the JSON filename for later merging
            individual_json_files(imol) = get_output_json_filename()

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

      ! Merge individual JSON files into one combined file (rank 0 only)
      if (my_rank == 0) then
         call merge_multi_molecule_json(individual_json_files, mqc_config%nmol)
      end if

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

   subroutine merge_multi_molecule_json(individual_files, nmol)
      !! Merge individual molecule JSON files into a single combined file
      character(len=256), intent(in) :: individual_files(:)
      integer, intent(in) :: nmol

      integer :: imol, unit_in, unit_out, io_stat, slash_pos, dot_pos
      character(len=10000) :: line
      character(len=256) :: output_file, basename
      logical :: file_exists

      ! Determine combined output filename from first individual file
      ! Example: "output_multi_structure_molecule_1.json" -> "output_multi_structure.json"
      basename = individual_files(1)
      slash_pos = index(basename, '/', back=.true.)
      if (slash_pos > 0) then
         basename = basename(slash_pos + 1:)
      end if

      ! Remove "_molecule_1" or similar suffix
      dot_pos = index(basename, '_molecule_')
      if (dot_pos > 0) then
         output_file = basename(1:dot_pos - 1)//".json"
      else
         output_file = "output_combined.json"
      end if

      ! Open combined output file
      open (newunit=unit_out, file=trim(output_file), status='replace', action='write', iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         return
      end if

      call logger%info("Merging "//to_char(nmol)//" molecule JSON files into "//trim(output_file))

      ! Write opening brace and top-level key (basename without "output_" and ".json")
      dot_pos = index(output_file, '.json')
      if (dot_pos > 0) then
         basename = output_file(8:dot_pos - 1)  ! Skip "output_"
      else
         basename = "combined"
      end if

      write (unit_out, '(a)') "{"
      write (unit_out, '(a)') '  "'//trim(basename)//'": {'

      ! Process each individual JSON file
      do imol = 1, nmol
         inquire (file=trim(individual_files(imol)), exist=file_exists)
         if (.not. file_exists) cycle

         open (newunit=unit_in, file=trim(individual_files(imol)), status='old', action='read', iostat=io_stat)
         if (io_stat /= 0) cycle

         ! Read all lines from the individual JSON file
         call read_json_content(unit_in, imol, unit_out, individual_files(imol))

         close (unit_in)

         ! Delete individual file
         open (newunit=unit_in, file=trim(individual_files(imol)), status='old', action='readwrite')
         close (unit_in, status='delete')
      end do

      ! Close last molecule
      write (unit_out, '(a)') '    }'

      ! Close top-level key and file
      write (unit_out, '(a)') '  }'
      write (unit_out, '(a)') '}'

      close (unit_out)
      call logger%info("Combined JSON written to "//trim(output_file))

   end subroutine merge_multi_molecule_json

   subroutine read_json_content(unit_in, mol_index, unit_out, filename)
      !! Read and write JSON content from an individual molecule file
      !! Properly handles nested structures from fragmented calculations
      integer, intent(in) :: unit_in, mol_index, unit_out
      character(len=*), intent(in) :: filename

      character(len=10000), allocatable :: all_lines(:)
      character(len=10000) :: line
      integer :: io_stat, nlines, i, brace_count

      ! Read all lines into memory
      allocate (all_lines(1000))  ! Reasonable size for most JSON files
      nlines = 0

      do
         read (unit_in, '(a)', iostat=io_stat) line
         if (io_stat /= 0) exit
         nlines = nlines + 1
         if (nlines > size(all_lines)) then
            ! Reallocate if needed
            call logger%error("JSON file too large: "//trim(filename))
            return
         end if
         all_lines(nlines) = line
      end do

      ! Lines structure:
      ! 1: "{"
      ! 2: '  "molecule_name": {'
      ! 3..(n-2): content
      ! n-1: "  }"
      ! n: "}"

      if (nlines < 3) then
         call logger%error("Invalid JSON structure: "//trim(filename))
         return
      end if

      ! Write molecule key (extracted from filename)
      if (mol_index > 1) write (unit_out, '(a)') '    },'
      write (unit_out, '(a)') '    "'//trim(get_molecule_name(filename))//'" : {'

      ! Write all content lines (from line 3 to line n-2)
      do i = 3, nlines - 2
         write (unit_out, '(a)') '  '//trim(all_lines(i))  ! Add 2 spaces for proper indentation
      end do

      deallocate (all_lines)

   end subroutine read_json_content

   function get_molecule_name(filename) result(name)
      !! Extract molecule name from filename
      !! Example: "output_multi_structure_molecule_1.json" -> "molecule_1"
      character(len=*), intent(in) :: filename
      character(len=256) :: name
      integer :: start_pos, end_pos

      ! Find "_molecule_" or similar pattern
      start_pos = index(filename, '_molecule_')
      if (start_pos == 0) start_pos = index(filename, '_mol_')

      if (start_pos > 0) then
         start_pos = start_pos + 1  ! Skip leading underscore
         end_pos = index(filename, '.json') - 1
         if (end_pos > start_pos) then
            name = filename(start_pos:end_pos)
         else
            name = "unknown"
         end if
      else
         name = "unknown"
      end if
   end function get_molecule_name

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
