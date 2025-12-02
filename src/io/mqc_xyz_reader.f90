!! XYZ molecular geometry file reader
module mqc_xyz_reader
   !! Provides functions to parse standard XYZ format files containing
   !! atomic coordinates and element symbols for molecular structures.
   use pic_types, only: dp
   use mqc_geometry, only: geometry_type
   implicit none
   private

   public :: read_xyz_file    !! Read XYZ file from disk
   public :: read_xyz_string  !! Parse XYZ data from string
   public :: split_lines      !! Split text into lines (for testing)

   ! Constants
   integer, parameter :: MAX_ELEMENT_SYMBOL_LEN = 4  !! Maximum element symbol length

contains

   subroutine read_xyz_file(filename, geom, stat, errmsg)
      !! Read molecular geometry from XYZ format file
      !!
      !! Parses standard XYZ files with format:
      !! Line 1: Number of atoms
      !! Line 2: Comment/title line
      !! Lines 3+: Element X Y Z (coordinates in Angstrom)
      character(len=*), intent(in) :: filename  !! Path to XYZ file
      type(geometry_type), intent(out) :: geom  !! Parsed molecular geometry
      integer, intent(out) :: stat              !! Status (0=success, >0=error)
      character(len=:), allocatable, intent(out) :: errmsg  !! Error message

      integer :: unit      !! File unit number
      integer :: io_stat   !! I/O operation status
      integer :: file_size  !! File size in bytes
      logical :: file_exists  !! Whether file exists on disk
      character(len=:), allocatable :: file_contents  !! Full file content buffer

      stat = 0

      ! Check if file exists
      inquire (file=filename, exist=file_exists, size=file_size)
      if (.not. file_exists) then
         stat = 1
         errmsg = "XYZ file not found: "//trim(filename)
         return
      end if

      ! Allocate buffer for entire file
      allocate (character(len=file_size) :: file_contents)

      ! Open and read entire file as stream
      open (newunit=unit, file=filename, status='old', action='read', &
            access='stream', form='unformatted', iostat=io_stat)
      if (io_stat /= 0) then
         stat = io_stat
         errmsg = "Error opening file: "//trim(filename)
         return
      end if

      read (unit, iostat=io_stat) file_contents
      close (unit)

      if (io_stat /= 0) then
         stat = io_stat
         errmsg = "Error reading file: "//trim(filename)
         return
      end if

      ! Parse the contents
      call read_xyz_string(file_contents, geom, stat, errmsg)

   end subroutine read_xyz_file

   pure subroutine read_xyz_string(xyz_string, geom, stat, errmsg)
      !! Parse molecular geometry from XYZ format string
      character(len=*), intent(in) :: xyz_string
      type(geometry_type), intent(out) :: geom
      integer, intent(out) :: stat
      character(len=:), allocatable, intent(out) :: errmsg

      character(len=:), allocatable :: lines(:)
      integer :: nlines, iatom, io_stat
      character(len=256) :: element
      real(dp) :: x, y, z

      stat = 0

      ! Split into lines
      call split_lines(xyz_string, lines, nlines)

      if (nlines < 2) then
         stat = 1
         errmsg = "XYZ file must have at least 2 lines (natoms + comment)"
         return
      end if

      ! Read number of atoms from first line
      read (lines(1), *, iostat=io_stat) geom%natoms
      if (io_stat /= 0) then
         stat = 1
         errmsg = "Failed to read number of atoms from first line"
         return
      end if

      if (geom%natoms < 0) then
         stat = 1
         errmsg = "Number of atoms must be non-negative"
         return
      end if

      ! Store comment line
      geom%comment = trim(adjustl(lines(2)))

      ! Check we have enough lines
      if (nlines < 2 + geom%natoms) then
         stat = 1
         errmsg = "XYZ file has insufficient lines: expected "// &
                  trim(int_to_string(2 + geom%natoms))//", got "// &
                  trim(int_to_string(nlines))
         return
      end if

      ! Allocate arrays
      allocate (character(len=MAX_ELEMENT_SYMBOL_LEN) :: geom%elements(geom%natoms))
      allocate (geom%coords(3, geom%natoms))

      ! Read atom data
      do iatom = 1, geom%natoms
         read (lines(2 + iatom), *, iostat=io_stat) element, x, y, z
         if (io_stat /= 0) then
            stat = 1
            errmsg = "Failed to parse atom data on line "// &
                     trim(int_to_string(2 + iatom))//": '"// &
                     trim(lines(2 + iatom))//"'"
            return
         end if

         geom%elements(iatom) = trim(adjustl(element))
         geom%coords(1, iatom) = x
         geom%coords(2, iatom) = y
         geom%coords(3, iatom) = z
      end do

   end subroutine read_xyz_string

   pure function int_to_string(i) result(str)
      !! Convert integer to string
      integer, intent(in) :: i
      character(len=:), allocatable :: str
      character(len=20) :: buffer

      write (buffer, '(I0)') i
      str = trim(adjustl(buffer))
   end function int_to_string

   pure subroutine split_lines(text, lines, nlines)
      !! Split input text into lines based on CR, LF, or CRLF line endings
      !! Trailing newlines do not create empty lines
      character(len=*), intent(in) :: text
      character(len=:), allocatable, intent(out) :: lines(:)
      integer, intent(out) :: nlines

      integer :: i, line_start, line_end, max_line_len
      character(len=:), allocatable :: temp_lines(:)

      if (len(text) == 0) then
         nlines = 0
         allocate (character(len=1) :: lines(0))
         return
      end if

      ! Pass 1: Count lines and find maximum line length
      nlines = 0
      max_line_len = 0
      line_start = 1
      i = 1

      do while (i <= len(text))
         ! Check for line ending
         if (text(i:i) == achar(13)) then  ! CR
            ! Check for CRLF
            if (i < len(text) .and. text(i + 1:i + 1) == achar(10)) then
               line_end = i - 1
               i = i + 2  ! Skip both CR and LF
            else
               line_end = i - 1
               i = i + 1
            end if
            nlines = nlines + 1
            max_line_len = max(max_line_len, line_end - line_start + 1)
            line_start = i
         else if (text(i:i) == achar(10)) then  ! LF
            line_end = i - 1
            nlines = nlines + 1
            max_line_len = max(max_line_len, line_end - line_start + 1)
            i = i + 1
            line_start = i
         else
            i = i + 1
         end if
      end do

      ! Handle last line if text doesn't end with newline
      if (line_start <= len(text)) then
         nlines = nlines + 1
         max_line_len = max(max_line_len, len(text) - line_start + 1)
      end if

      ! Handle empty text or ensure at least length 1
      if (max_line_len == 0) max_line_len = 1

      ! Allocate output array
      allocate (character(len=max_line_len) :: temp_lines(nlines))

      ! Pass 2: Extract lines
      nlines = 0
      line_start = 1
      i = 1

      do while (i <= len(text))
         ! Check for line ending
         if (text(i:i) == achar(13)) then  ! CR
            ! Check for CRLF
            if (i < len(text) .and. text(i + 1:i + 1) == achar(10)) then
               line_end = i - 1
               i = i + 2
            else
               line_end = i - 1
               i = i + 1
            end if
            nlines = nlines + 1
            temp_lines(nlines) = ""  ! Initialize line before copying
            if (line_end >= line_start) then
               ! Intel compiler workaround: use character-by-character copy
               block
                  integer :: j, line_len
                  line_len = line_end - line_start + 1
                  do j = 1, line_len
                     temp_lines(nlines) (j:j) = text(line_start + j - 1:line_start + j - 1)
                  end do
               end block
            end if
            line_start = i
         else if (text(i:i) == achar(10)) then  ! LF
            line_end = i - 1
            nlines = nlines + 1
            temp_lines(nlines) = ""  ! Initialize line before copying
            if (line_end >= line_start) then
               ! Intel compiler workaround: use character-by-character copy
               block
                  integer :: j, line_len
                  line_len = line_end - line_start + 1
                  do j = 1, line_len
                     temp_lines(nlines) (j:j) = text(line_start + j - 1:line_start + j - 1)
                  end do
               end block
            end if
            i = i + 1
            line_start = i
         else
            i = i + 1
         end if
      end do

      ! Handle last line if text doesn't end with newline
      if (line_start <= len(text)) then
         nlines = nlines + 1
         temp_lines(nlines) = ""  ! Initialize line before copying
         ! Intel compiler workaround: use character-by-character copy
         block
            integer :: j, line_len
            line_len = len(text) - line_start + 1
            do j = 1, line_len
               temp_lines(nlines) (j:j) = text(line_start + j - 1:line_start + j - 1)
            end do
         end block
      end if

      ! Copy to output (use explicit loop for Intel compiler compatibility)
      allocate (character(len=max_line_len) :: lines(nlines))
      block
         integer :: iline
         do iline = 1, nlines
            lines(iline) = temp_lines(iline)
         end do
      end block

   end subroutine split_lines

end module mqc_xyz_reader
