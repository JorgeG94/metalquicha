!! this file contains the modules and routines to open and read a GAMESS formatted basis set file
module mqc_basis_file_reader
   !! Module for reading and parsing GAMESS formatted basis set files
   use pic_types, only: int32, dp
   use mqc_error, only: error_t, ERROR_IO, ERROR_VALIDATION
   implicit none

   private
   public :: basis_file_t, open_basis_file, extract_element, strings_equal

   type :: basis_file_t
      !! Container for basis set file contents
      character(len=:), allocatable :: full_content
      character(len=:), allocatable :: data_section
   end type basis_file_t

contains

   subroutine open_basis_file(basis_file, filename, error)
      !! Open and read a GAMESS formatted basis set file
      type(basis_file_t), intent(out) :: basis_file
      character(len=*), intent(in) :: filename
      type(error_t), intent(out) :: error

      integer :: unit, iostat, file_size
      logical :: file_exists
      integer :: data_start, data_end

      ! Check if file exists
      inquire (file=filename, exist=file_exists, size=file_size)
      if (.not. file_exists) then
         call error%set(ERROR_IO, "Basis set file not found: "//filename)
         return
      end if

      ! Allocate buffer for entire file
      allocate (character(len=file_size) :: basis_file%full_content)

      ! Open and read entire file
      open (newunit=unit, file=filename, status='old', action='read', &
            access='stream', form='unformatted', iostat=iostat)
      if (iostat /= 0) then
         call error%set(ERROR_IO, "Error opening file: "//filename)
         return
      end if

      read (unit, iostat=iostat) basis_file%full_content
      if (iostat /= 0) then
         close (unit)
         call error%set(ERROR_IO, "Error reading file: "//filename)
         return
      end if
      close (unit)

      ! Extract the $DATA section
      data_start = index(basis_file%full_content, "$DATA")
      if (data_start == 0) then
         call error%set(ERROR_VALIDATION, "Could not find $DATA section in basis set file: "//filename)
         return
      end if

      data_end = index(basis_file%full_content(data_start:), "$END")
      if (data_end == 0) then
         call error%set(ERROR_VALIDATION, "Could not find $END marker in basis set file: "//filename)
         return
      end if

      ! Store just the data section (between $DATA and $END)
      basis_file%data_section = basis_file%full_content(data_start + 5:data_start + data_end - 2)

   end subroutine open_basis_file

   subroutine extract_element(basis_file, element, element_content, error)
      !! Extract the basis set data for a specific element from the basis file
      type(basis_file_t), intent(in) :: basis_file
      character(len=*), intent(in) :: element
      character(len=:), allocatable, intent(out) :: element_content
      type(error_t), intent(out) :: error

      integer :: start_pos, end_pos, i
      character(len=:), allocatable :: search_element
      logical :: at_line_start

      ! Convert element to uppercase for searching
      search_element = uppercase(trim(element))

      ! Find the element name (it appears on its own line)
      start_pos = index(basis_file%data_section, new_line('a')//trim(search_element)//new_line('a'))

      if (start_pos == 0) then
         ! Try without leading newline (might be first element after $DATA)
         if (index(basis_file%data_section, trim(search_element)//new_line('a')) == 1) then
            start_pos = 1
         else
            call error%set(ERROR_VALIDATION, "Element not found in basis set file: "//element)
            return
         end if
      else
         start_pos = start_pos + 1  ! Skip the leading newline
      end if

      ! Find the next element by looking for a line that:
      ! - Starts with an uppercase letter
      ! - Has a second character that is also a letter (not a space or number)
      ! This distinguishes "CARBON" from "S   3"

      end_pos = len(basis_file%data_section)
      at_line_start = .false.

      i = start_pos + len(search_element) + 1
      do while (i < len(basis_file%data_section))
         if (basis_file%data_section(i:i) == new_line('a')) then
            at_line_start = .true.
            i = i + 1
            cycle
         end if

         if (at_line_start) then
            ! We're at the start of a new line
            if (is_uppercase_letter(basis_file%data_section(i:i))) then
               ! Check if next character is also a letter
               if (i + 1 <= len(basis_file%data_section)) then
                  if (is_letter(basis_file%data_section(i + 1:i + 1))) then
                     ! Found next element!
                     end_pos = i - 1
                     exit
                  end if
               end if
            end if
            at_line_start = .false.
         end if

         i = i + 1
      end do

      ! Extract the section
      element_content = basis_file%data_section(start_pos:end_pos)

   end subroutine extract_element

   pure function is_letter(c) result(is_alpha)
      !! Check if character is a letter (A-Z or a-z)
      character(len=1), intent(in) :: c
      logical :: is_alpha
      integer :: ic

      ic = iachar(c)
      is_alpha = (ic >= iachar('A') .and. ic <= iachar('Z')) .or. &
                 (ic >= iachar('a') .and. ic <= iachar('z'))
   end function is_letter

   pure function uppercase(str) result(upper)
      !! Convert a string to uppercase, should use pic_ascii!
      character(len=*), intent(in) :: str
      character(len=:), allocatable :: upper
      integer :: i, ic

      allocate (character(len=len(str)) :: upper)
      upper = str

      do i = 1, len(str)
         ic = iachar(str(i:i))
         if (ic >= iachar('a') .and. ic <= iachar('z')) then
            upper(i:i) = achar(ic - 32)
         end if
      end do
   end function uppercase

   pure function is_uppercase_letter(c) result(is_upper)
      !! Check if character is an uppercase letter (A-Z)
      character(len=1), intent(in) :: c
      logical :: is_upper
      integer :: ic

      ic = iachar(c)
      is_upper = (ic >= iachar('A') .and. ic <= iachar('Z'))
   end function is_uppercase_letter

   !> Compare two strings after trimming and adjusting (removing leading/trailing whitespace)
   pure function strings_equal(str1, str2) result(equal)
      !! Compare two strings for equality after trimming and adjusting (removing leading/trailing whitespace)
      character(len=*), intent(in) :: str1, str2
      logical :: equal
      equal = trim(adjustl(str1)) == trim(adjustl(str2))
   end function strings_equal

end module mqc_basis_file_reader
