!! Command line argument parsing for metalquicha
module mqc_cli_parser
   !! Handles parsing of command line options including geometry files,
   !! basis set specifications, and help/usage display.
   use mqc_basis_utils, only: normalize_basis_name, find_basis_file
   use mqc_error, only: error_t, ERROR_PARSE, ERROR_IO
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

   subroutine parse_command_line(args, error)
      !! Parse command line arguments for geometry file and basis set
      !!
      !! Extracts XYZ file path and basis set name from command line,
      !! validates arguments, and handles help requests.
      type(cli_args_type), intent(out) :: args  !! Parsed argument container
      type(error_t), intent(out) :: error       !! Error object

      integer :: nargs        !! Number of command line arguments
      character(len=256) :: arg_buffer  !! Temporary argument buffer
      integer :: arg_len      !! Length of current argument
      integer :: stat         !! Local status for intrinsic calls

      ! Get number of command line arguments
      nargs = command_argument_count()

      ! Check for help flag
      if (nargs >= 1) then
         call get_command_argument(1, arg_buffer, arg_len, stat)
         if (stat /= 0) then
            call error%set(ERROR_PARSE, "Error reading command line argument 1")
            return
         end if
         arg_buffer = trim(arg_buffer)

         if (arg_buffer == "-h" .or. arg_buffer == "--help") then
            call print_usage()
            call error%set(ERROR_PARSE, "HELP_REQUESTED")  ! Special marker for help
            return
         end if
      end if

      ! Validate number of arguments
      if (nargs < 2) then
         call error%set(ERROR_PARSE, "Error: Insufficient arguments. Expected 2 arguments (geometry.xyz basis_name)")
         call print_usage()
         return
      end if

      if (nargs > 2) then
         call error%set(ERROR_PARSE, "Error: Too many arguments. Expected 2 arguments (geometry.xyz basis_name)")
         call print_usage()
         return
      end if

      ! Parse argument 1: XYZ file
      call get_command_argument(1, arg_buffer, arg_len, stat)
      if (stat /= 0) then
         call error%set(ERROR_PARSE, "Error reading geometry file argument")
         return
      end if
      args%xyz_file = trim(arg_buffer)

      ! Parse argument 2: Basis set name
      call get_command_argument(2, arg_buffer, arg_len, stat)
      if (stat /= 0) then
         call error%set(ERROR_PARSE, "Error reading basis set name argument")
         return
      end if
      args%basis_name = trim(arg_buffer)

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

end module mqc_cli_parser
