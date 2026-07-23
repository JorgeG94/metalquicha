!! Utilities for handling basis set names and files
module mqc_basis_utils
   !! Provides utilities for normalizing basis set names and locating basis set files
   !!
   !! Normalization rules:
   !!   * -> s   (e.g., 6-31G* -> 6-31Gs)
   !!   + -> p   (e.g., 6-31+G -> 6-31pG)
   !!   (d,p) -> dp (remove parentheses and commas)
   use mqc_error, only: error_t, ERROR_IO
   implicit none
   private

   public :: normalize_basis_name
   public :: find_basis_file

   integer, parameter, public :: BASIS_FORMAT_UNKNOWN = 0  !! No format determined
   integer, parameter, public :: BASIS_FORMAT_GAMESS = 1   !! GAMESS $DATA, `.txt`
   integer, parameter, public :: BASIS_FORMAT_GBS = 2      !! Gaussian94, `.gbs`
   integer, parameter, public :: BASIS_FORMAT_JSON = 3     !! Basis Set Exchange JSON, `.json`

contains

   pure function normalize_basis_name(basis_name) result(normalized)
      !! Normalize basis set name to filename-safe format
      !!
      !! Rules:
      !!   * -> s
      !!   + -> p
      !!   Remove parentheses and commas
      !!
      !! Examples:
      !!   6-31G*      -> 6-31Gs
      !!   6-31+G*     -> 6-31pGs
      !!   6-31G(d)    -> 6-31Gd
      !!   6-311G(d,p) -> 6-311Gdp
      !!   6-311++G**  -> 6-311ppGss
      !!   cc-pVDZ     -> cc-pVDZ (unchanged)
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
            ! Skip commas and spaces (inside or outside parentheses)
            continue

         case default
            ! Copy character as-is
            out_pos = out_pos + 1
            buffer(out_pos:out_pos) = basis_name(i:i)
         end select
      end do

      normalized = trim(buffer(1:out_pos))
   end function normalize_basis_name

   subroutine find_basis_file(basis_name, filename, error, basis_format)
      !! Find basis set file using normalized name
      !!
      !! Search strategy:
      !!   1. Normalize the basis name (e.g., 6-31G* -> 6-31Gs)
      !!   2. basis_sets/{normalized}.json  (Basis Set Exchange JSON)
      !!   3. basis_sets/{normalized}.gbs   (Gaussian94)
      !!   4. basis_sets/{normalized}.txt   (GAMESS $DATA)
      !!   5. If none exists, return error
      !!
      !! JSON is preferred: it keeps a combined shell's coefficient sets
      !! separate, so SP shells split cleanly instead of needing the
      !! `uncontract_spdf=1` download flag, and it carries ECP data too.
      character(len=*), intent(in) :: basis_name
      character(len=:), allocatable, intent(out) :: filename
      type(error_t), intent(out) :: error
      integer, intent(out), optional :: basis_format  !! Which BASIS_FORMAT_* was found

      character(len=:), allocatable :: normalized
      logical :: json_exists, gbs_exists, gamess_exists
      character(len=512) :: json_path, gbs_path, gamess_path

      normalized = normalize_basis_name(basis_name)

      json_path = "basis_sets/"//trim(normalized)//".json"
      gbs_path = "basis_sets/"//trim(normalized)//".gbs"
      gamess_path = "basis_sets/"//trim(normalized)//".txt"

      inquire (file=trim(json_path), exist=json_exists)
      inquire (file=trim(gbs_path), exist=gbs_exists)
      inquire (file=trim(gamess_path), exist=gamess_exists)

      if (present(basis_format)) basis_format = BASIS_FORMAT_UNKNOWN

      if (json_exists) then
         filename = trim(json_path)
         if (present(basis_format)) basis_format = BASIS_FORMAT_JSON
      else if (gbs_exists) then
         filename = trim(gbs_path)
         if (present(basis_format)) basis_format = BASIS_FORMAT_GBS
      else if (gamess_exists) then
         filename = trim(gamess_path)
         if (present(basis_format)) basis_format = BASIS_FORMAT_GAMESS
      else
         call error%set(ERROR_IO, "Basis set file not found: tried "// &
                        trim(json_path)//", "//trim(gbs_path)//" and "//trim(gamess_path)// &
                        " (from basis name: "//trim(basis_name)//")")
      end if

   end subroutine find_basis_file

end module mqc_basis_utils
