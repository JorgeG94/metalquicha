!! Command line argument parsing for metalquicha
module mqc_cli_parser
   !! Handles parsing of command line options including geometry files,
   !! basis set specifications, and help/usage display.
   implicit none
   private

   public :: cli_args_type         !! Parsed command line arguments container
   public :: parse_command_line    !! Main argument parsing routine
   public :: print_usage           !! Display program usage information
   public :: normalize_basis_name  !! Standardize basis set names
   public :: find_basis_file       !! Locate basis set files

   type :: cli_args_type
      !! Container for parsed command line arguments
      !!
      !! Stores file paths and options extracted from command line,
      !! with automatic memory management for string allocations.
      character(len=:), allocatable :: xyz_file    !! Input XYZ geometry file path
      character(len=:), allocatable :: basis_name  !! Basis set name (e.g., "6-31G")
   contains
      procedure :: destroy => cli_args_destroy  !! Memory cleanup
   end type cli_args_type

contains

   subroutine parse_command_line(args, stat, errmsg)
      !! Parse command line arguments for geometry file and basis set
      !!
      !! Extracts XYZ file path and basis set name from command line,
      !! validates arguments, and handles help requests.
      type(cli_args_type), intent(out) :: args  !! Parsed argument container
      integer, intent(out) :: stat              !! Status (0=success, >0=error)
      character(len=:), allocatable, intent(out) :: errmsg  !! Error message

      integer :: nargs        !! Number of command line arguments
      character(len=256) :: arg_buffer  !! Temporary argument buffer
      integer :: arg_len      !! Length of current argument

      stat = 0

      ! Get number of command line arguments
      nargs = command_argument_count()

      ! Check for help flag
      if (nargs >= 1) then
         call get_command_argument(1, arg_buffer, arg_len, stat)
         if (stat /= 0) then
            errmsg = "Error reading command line argument 1"
            return
         end if
         arg_buffer = trim(arg_buffer)

         if (arg_buffer == "-h" .or. arg_buffer == "--help") then
            call print_usage()
            stat = -1  ! Special code to indicate help was requested
            return
         end if
      end if

      ! Validate number of arguments
      if (nargs < 2) then
         stat = 1
         errmsg = "Error: Insufficient arguments. Expected 2 arguments (geometry.xyz basis_name)"
         call print_usage()
         return
      end if

      if (nargs > 2) then
         stat = 1
         errmsg = "Error: Too many arguments. Expected 2 arguments (geometry.xyz basis_name)"
         call print_usage()
         return
      end if

      ! Parse argument 1: XYZ file
      call get_command_argument(1, arg_buffer, arg_len, stat)
      if (stat /= 0) then
         errmsg = "Error reading geometry file argument"
         return
      end if
      args%xyz_file = trim(arg_buffer)

      ! Parse argument 2: Basis set name
      call get_command_argument(2, arg_buffer, arg_len, stat)
      if (stat /= 0) then
         errmsg = "Error reading basis set name argument"
         return
      end if
      args%basis_name = trim(arg_buffer)

      ! Reset stat to success
      stat = 0

   end subroutine parse_command_line

   !> Print usage information
   subroutine print_usage()
      character(len=256) :: prog_name
      integer :: stat

      call get_command_argument(0, prog_name, status=stat)
      if (stat /= 0) prog_name = "pic_basis_reader"

      print *
      print *, "Usage: ", trim(prog_name), " <geometry.xyz> <basis_name>"
      print *
      print *, "Arguments:"
      print *, "  geometry.xyz   XYZ format molecular geometry file"
      print *, "  basis_name     Name of basis set (e.g., 6-31G, 6-311G**)"
      print *
      print *, "Options:"
      print *, "  -h, --help     Show this help message"
      print *
      print *, "Example:"
      print *, "  ", trim(prog_name), " water.xyz 6-31G"
      print *

   end subroutine print_usage

   !> Clean up CLI args
   subroutine cli_args_destroy(this)
      class(cli_args_type), intent(inout) :: this
      if (allocated(this%xyz_file)) deallocate (this%xyz_file)
      if (allocated(this%basis_name)) deallocate (this%basis_name)
   end subroutine cli_args_destroy

   !> Normalize basis set name to filename-safe format
  !! Rules: * -> s, + -> p, remove parentheses and their contents
  !! Examples:
  !!   6-31G*      -> 6-31Gs
  !!   6-31+G*     -> 6-31pGs
  !!   6-31G(d)    -> 6-31Gd
  !!   6-311G(d,p) -> 6-31Gdp
  !!   6-311++G**  -> 6-311ppGss
   function normalize_basis_name(basis_name) result(normalized)
      character(len=*), intent(in) :: basis_name
      character(len=:), allocatable :: normalized
      integer :: i, out_pos
      character(len=256) :: buffer
      logical :: in_parens

      buffer = ""
      out_pos = 0
      in_parens = .false.

      do i = 1, len_trim(basis_name)
         select case (basis_name(i:i))
         case ('*')
            ! Star becomes 's'
            out_pos = out_pos + 1
            buffer(out_pos:out_pos) = 's'

         case ('+')
            ! Plus becomes 'p'
            out_pos = out_pos + 1
            buffer(out_pos:out_pos) = 'p'

         case ('(')
            ! Start of parentheses - we'll extract contents
            in_parens = .true.

         case (')')
            ! End of parentheses
            in_parens = .false.

         case (',', ' ')
            ! Skip commas and spaces inside parentheses
            if (in_parens) cycle
            ! Keep them outside parentheses
            out_pos = out_pos + 1
            buffer(out_pos:out_pos) = basis_name(i:i)

         case default
            ! Copy character as-is
            out_pos = out_pos + 1
            buffer(out_pos:out_pos) = basis_name(i:i)
         end select
      end do

      normalized = trim(buffer(1:out_pos))
   end function normalize_basis_name

   !> Find basis file, trying multiple name variants
  !! Searches in:
  !!   1. basis_sets/ subdirectory
  !!   2. Current directory
  !! Tries name variants:
  !!   1. Exact name as given
  !!   2. Normalized name (with * -> s, + -> p, etc.)
  !!   3. Common synonyms
   subroutine find_basis_file(basis_name, filename, stat, errmsg)
      character(len=*), intent(in) :: basis_name
      character(len=:), allocatable, intent(out) :: filename
      integer, intent(out) :: stat
      character(len=:), allocatable, intent(out) :: errmsg

      character(len=:), allocatable :: normalized
      logical :: file_exists
      integer :: i, j
      character(len=256), dimension(10) :: variants
      integer :: nvariants
      character(len=256) :: test_path
      character(len=256), dimension(2) :: search_dirs

      stat = 0
      nvariants = 0

      ! Define search directories (in order of priority)
      search_dirs(1) = "basis_sets/"
      search_dirs(2) = ""  ! Current directory

      ! Variant 1: Exact name
      nvariants = nvariants + 1
      variants(nvariants) = trim(basis_name)//".txt"

      ! Variant 2: Normalized name
      normalized = normalize_basis_name(basis_name)
      if (trim(normalized) /= trim(basis_name)) then
         nvariants = nvariants + 1
         variants(nvariants) = trim(normalized)//".txt"
      end if

      ! Variant 3: Common synonyms
      ! 6-31G* = 6-31G(d)
      if (trim(basis_name) == "6-31G*" .or. trim(basis_name) == "6-31G(d)") then
         nvariants = nvariants + 1
         variants(nvariants) = "6-31Gs.txt"
         nvariants = nvariants + 1
         variants(nvariants) = "6-31Gd.txt"
      end if

      ! 6-31G** = 6-31G(d,p)
      if (trim(basis_name) == "6-31G**" .or. trim(basis_name) == "6-31G(d,p)") then
         nvariants = nvariants + 1
         variants(nvariants) = "6-31Gss.txt"
         nvariants = nvariants + 1
         variants(nvariants) = "6-31Gdp.txt"
      end if

      ! 6-311G* = 6-311G(d)
      if (trim(basis_name) == "6-311G*" .or. trim(basis_name) == "6-311G(d)") then
         nvariants = nvariants + 1
         variants(nvariants) = "6-311Gs.txt"
         nvariants = nvariants + 1
         variants(nvariants) = "6-311Gd.txt"
      end if

      ! Try each directory with each variant
      do j = 1, 2
         do i = 1, nvariants
            test_path = trim(search_dirs(j))//trim(variants(i))
            inquire (file=trim(test_path), exist=file_exists)
            if (file_exists) then
               filename = trim(test_path)
               return
            end if
         end do
      end do

      ! None found - return error with all paths tried
      stat = 1
      errmsg = "Basis set file not found. Tried: "
      do j = 1, 2
         do i = 1, nvariants
            test_path = trim(search_dirs(j))//trim(variants(i))
            if (i == 1 .and. j == 1) then
               errmsg = errmsg//trim(test_path)
            else
               errmsg = errmsg//", "//trim(test_path)
            end if
         end do
      end do

   end subroutine find_basis_file

end module mqc_cli_parser
