module mqc_input_parser
   !! Input file parser for PIC Chemistry configuration
   !!
   !! Parses simple key=value input files to configure calculation parameters
   !! including geometry files, method selection, and fragment levels.
   implicit none
   private

   public :: input_config_t, read_input_file  !! Main configuration type and parser

   type :: input_config_t
      !! Configuration data structure for calculation parameters
      !!
      !! Stores parsed input file contents including file paths,
      !! method selection, and fragmentation level.
      character(len=:), allocatable :: geom_file     !! Path to system geometry XYZ file
      character(len=:), allocatable :: monomer_file  !! Path to monomer template XYZ file
      character(len=:), allocatable :: method        !! QC method (gfn1, gfn2)
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
      !!            nlevel=N (fragmentation level, defaults to 1)
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

      ! Check if file exists
      inquire (file=filename, exist=file_exists)
      if (.not. file_exists) then
         stat = 1
         errmsg = "Input file not found: "//trim(filename)
         return
      end if

      ! Open input file
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

   subroutine config_destroy(this)
      class(input_config_t), intent(inout) :: this
      if (allocated(this%geom_file)) deallocate (this%geom_file)
      if (allocated(this%monomer_file)) deallocate (this%monomer_file)
      if (allocated(this%method)) deallocate (this%method)
   end subroutine config_destroy

end module mqc_input_parser
