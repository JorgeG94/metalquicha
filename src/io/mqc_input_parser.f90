!! Input file parser for the hastily put together input file format
module mqc_input_parser
   !! Parses simple key=value input files to configure calculation parameters
   !! including geometry files, method selection, and fragment levels.
   implicit none
   private

   public :: input_config_t, read_input_file, get_logger_level  !! Main configuration type and parser

   type :: input_config_t
      !! Configuration data structure for calculation parameters
      !!
      !! Stores parsed input file contents including file paths,
      !! method selection, and fragmentation level.
      character(len=:), allocatable :: geom_file     !! Path to system geometry XYZ file
      character(len=:), allocatable :: monomer_file  !! Path to monomer template XYZ file
      character(len=:), allocatable :: method        !! QC method (gfn1, gfn2)
      character(len=:), allocatable :: calc_type     !! Calculation type (energy, gradient)
      character(len=:), allocatable :: log_level     !! Logger verbosity level (debug/verbose/info/warning/error)
      integer :: nlevel = 1  !! Fragmentation level (default: 1)
   contains
      procedure :: destroy => config_destroy  !! Cleanup allocated memory
   end type input_config_t

contains

   subroutine read_input_file(filename, config, stat, errmsg)
      !! Simple parser for key=value input files
      !! Looks for: geom="path/to/geometry.xyz"
      !!            monomer_symbols="path/to/monomer.xyz"
      !!            method="gfn1" or "gfn2" (defaults to gfn2)
      !!            calc_type="energy" or "gradient" (defaults to energy)
      !!            nlevel=N (fragmentation level, defaults to 1)
      !!            log_level="debug|verbose|info|performance|warning|error|knowledge" (defaults to info)
      character(len=*), intent(in) :: filename  !! Path to input file to parse
      type(input_config_t), intent(out) :: config  !! Parsed configuration data
      integer, intent(out) :: stat  !! Status code (0 = success, >0 = error)
      character(len=:), allocatable, intent(out) :: errmsg  !! Error message on failure

      integer :: unit     !! File unit number
      integer :: io_stat  !! I/O operation status
      character(len=512) :: line   !! Current line being parsed
      character(len=512) :: key    !! Parsed key name
      character(len=512) :: value  !! Parsed value string
      integer :: eq_pos   !! Position of '=' character in line
      logical :: file_exists  !! Whether input file exists

      stat = 0

      inquire (file=filename, exist=file_exists)
      if (.not. file_exists) then
         stat = 1
         errmsg = "Input file not found: "//trim(filename)
         return
      end if

      open (newunit=unit, file=filename, status='old', action='read', iostat=io_stat)
      if (io_stat /= 0) then
         stat = io_stat
         errmsg = "Error opening input file: "//trim(filename)
         return
      end if

      ! Read line by line
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) exit

         ! Skip empty lines and comments
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#' .or. line(1:1) == '!') cycle

         ! Find '=' sign
         eq_pos = index(line, '=')
         if (eq_pos == 0) cycle

         ! Extract key and value
         key = adjustl(line(1:eq_pos - 1))
         value = adjustl(line(eq_pos + 1:))

         ! Remove quotes from value if present
         value = remove_quotes(value)

         ! Parse recognized keys
         select case (trim(key))
         case ('geom')
            config%geom_file = trim(value)
         case ('monomer_symbols')
            config%monomer_file = trim(value)
         case ('method')
            ! Validate that method is gfn1 or gfn2
            select case (trim(value))
            case ('gfn1', 'gfn2')
               config%method = trim(value)
            case default
               stat = 1
               errmsg = "Invalid method: "//trim(value)//" (supported: gfn1, gfn2)"
               return
            end select
         case ('calc_type')
            ! Validate that calc_type is energy or gradient
            select case (trim(value))
            case ('energy', 'gradient')
               config%calc_type = trim(value)
            case default
               stat = 1
               errmsg = "Invalid calc_type: "//trim(value)//" (supported: energy, gradient)"
               return
            end select
         case ('nlevel')
            read (value, *, iostat=io_stat) config%nlevel
            if (io_stat /= 0) then
               stat = 1
               errmsg = "Invalid value for nlevel: "//trim(value)
               return
            end if
            if (config%nlevel < 0) then
               stat = 1
               errmsg = "nlevel must be >= 0 (0 for unfragmented calculation)"
               return
            end if
         case ('log_level')
            ! Validate log level
            select case (trim(value))
            case ('debug', 'Debug', 'DEBUG', &
                  'verbose', 'Verbose', 'VERBOSE', &
                  'info', 'Info', 'INFO', &
                  'performance', 'Performance', 'PERFORMANCE', &
                  'warning', 'Warning', 'WARNING', &
                  'error', 'Error', 'ERROR', &
                  'knowledge', 'Knowledge', 'KNOWLEDGE')
               config%log_level = trim(value)
            case default
               stat = 1
               errmsg = "Invalid log_level: "//trim(value)// &
                        " (supported: debug, verbose, info, performance, warning, error, knowledge)"
               return
            end select
         case default
            ! Ignore unrecognized keys
            continue
         end select
      end do

      close (unit)

      ! Validate required fields
      if (.not. allocated(config%geom_file)) then
         stat = 1
         errmsg = "Missing required field: geom"
         return
      end if

      if (.not. allocated(config%monomer_file)) then
         stat = 1
         errmsg = "Missing required field: monomer_symbols"
         return
      end if

      ! Set default method if not specified
      if (.not. allocated(config%method)) then
         config%method = "gfn2"  ! Default to GFN2-xTB
      end if

      ! Set default calc_type if not specified
      if (.not. allocated(config%calc_type)) then
         config%calc_type = "energy"  ! Default to energy-only calculation
      end if

      ! Set default log_level if not specified
      if (.not. allocated(config%log_level)) then
         config%log_level = "info"  ! Default to info level
      end if

   end subroutine read_input_file

   function remove_quotes(str) result(cleaned)
      !! Remove surrounding quotes from string
      character(len=*), intent(in) :: str
      character(len=:), allocatable :: cleaned
      character(len=len(str)) :: temp
      integer :: start_pos, end_pos

      temp = adjustl(str)
      start_pos = 1
      end_pos = len_trim(temp)

      ! Remove leading quote
      if (end_pos > 0) then
         if (temp(1:1) == '"' .or. temp(1:1) == "'") then
            start_pos = 2
         end if
      end if

      ! Remove trailing quote
      if (end_pos > 0) then
         if (temp(end_pos:end_pos) == '"' .or. temp(end_pos:end_pos) == "'") then
            end_pos = end_pos - 1
         end if
      end if

      if (end_pos >= start_pos) then
         cleaned = temp(start_pos:end_pos)
      else
         cleaned = ""
      end if

   end function remove_quotes

   function get_logger_level(level_string) result(level_int)
      !! Convert string log level to integer value
      !! This function uses the pic_logger constants
      use pic_logger, only: debug_level, verbose_level, info_level, performance_level, &
                            warning_level, error_level, knowledge_level
      character(len=*), intent(in) :: level_string
      integer :: level_int

      select case (trim(adjustl(level_string)))
      case ('debug', 'Debug', 'DEBUG')
         level_int = debug_level
      case ('verbose', 'Verbose', 'VERBOSE')
         level_int = verbose_level
      case ('info', 'Info', 'INFO')
         level_int = info_level
      case ('performance', 'Performance', 'PERFORMANCE')
         level_int = performance_level
      case ('warning', 'Warning', 'WARNING')
         level_int = warning_level
      case ('error', 'Error', 'ERROR')
         level_int = error_level
      case ('knowledge', 'Knowledge', 'KNOWLEDGE')
         level_int = knowledge_level
      case default
         ! Default to info level if unknown
         level_int = info_level
      end select
   end function get_logger_level

   subroutine config_destroy(this)
      !! Clean up allocated memory in input_config_t
      class(input_config_t), intent(inout) :: this
      if (allocated(this%geom_file)) deallocate (this%geom_file)
      if (allocated(this%monomer_file)) deallocate (this%monomer_file)
      if (allocated(this%method)) deallocate (this%method)
      if (allocated(this%calc_type)) deallocate (this%calc_type)
      if (allocated(this%log_level)) deallocate (this%log_level)
   end subroutine config_destroy

end module mqc_input_parser
