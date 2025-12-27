!! Utilities for handling basis set names and files
module mqc_basis_utils
   !! Provides utilities for normalizing basis set names and locating basis set files
   !!
   !! Normalization rules:
   !!   * -> s   (e.g., 6-31G* -> 6-31Gs)
   !!   + -> p   (e.g., 6-31+G -> 6-31pG)
   !!   (d,p) -> dp (remove parentheses and commas)
   implicit none
   private

   public :: normalize_basis_name
   public :: find_basis_file

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

   subroutine find_basis_file(basis_name, filename, stat, errmsg)
      !! Find basis set file using normalized name
      !!
      !! Search strategy:
      !!   1. Normalize the basis name (e.g., 6-31G* -> 6-31Gs)
      !!   2. Look for basis_sets/{normalized}.txt
      !!   3. If not found, return error
      !!
      !! This is a simple, straightforward approach that assumes
      !! the JSON/mqc input provides the correct basis set name.
      character(len=*), intent(in) :: basis_name
      character(len=:), allocatable, intent(out) :: filename
      integer, intent(out) :: stat
      character(len=:), allocatable, intent(out) :: errmsg

      character(len=:), allocatable :: normalized
      logical :: file_exists
      character(len=512) :: filepath

      stat = 0

      ! Normalize the basis name
      normalized = normalize_basis_name(basis_name)

      ! Construct file path: basis_sets/{normalized}.txt
      filepath = "basis_sets/"//trim(normalized)//".txt"

      ! Check if file exists
      inquire (file=trim(filepath), exist=file_exists)

      if (file_exists) then
         filename = trim(filepath)
      else
         stat = 1
         errmsg = "Basis set file not found: "//trim(filepath)// &
                  " (from basis name: "//trim(basis_name)//")"
      end if

   end subroutine find_basis_file

end module mqc_basis_utils
