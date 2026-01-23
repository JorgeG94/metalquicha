!! Structure and geometry section parsers for MQC config files
!! Handles: structure, geometry sections and their generic helpers
submodule(mqc_config_parser) mqc_config_parser_structure
   implicit none

contains

   module subroutine parse_structure_section(unit, config, error)
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      call parse_structure_generic(unit, config%charge, config%multiplicity, error)

   end subroutine parse_structure_section

   module subroutine parse_geometry_section(unit, config, error)
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      call parse_geometry_generic(unit, config%geometry, error)

   end subroutine parse_geometry_section

   module subroutine parse_structure_generic(unit, charge, multiplicity, error)
      !! Generic parser for %structure section (works for both config and molecule)
      integer, intent(in) :: unit
      integer, intent(inout) :: charge, multiplicity
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value
      integer :: io_stat, eq_pos
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_IO, "Unexpected end of file in %structure section")
            return
         end if

         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#' .or. line(1:1) == '!') cycle

         if (trim(strip_comment(line)) == 'end') exit

         eq_pos = index(line, '=')
         if (eq_pos == 0) cycle

         key = adjustl(line(1:eq_pos - 1))
         value = adjustl(line(eq_pos + 1:))

         select case (trim(key))
         case ('charge')
            read (value, *, iostat=io_stat) charge
         case ('multiplicity')
            read (value, *, iostat=io_stat) multiplicity
         case default
            call error%set(ERROR_PARSE, "Unknown key in %structure section: "//trim(key))
            return
         end select
      end do

   end subroutine parse_structure_generic

   module subroutine parse_geometry_generic(unit, geom, error)
      !! Generic parser for %geometry section (works for both config and molecule)
      use mqc_geometry, only: geometry_type
      integer, intent(in) :: unit
      type(geometry_type), intent(inout) :: geom
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, elem
      integer :: io_stat, natoms, i
      real(dp) :: x, y, z
      ! Read number of atoms
      read (unit, '(A)', iostat=io_stat) line
      if (io_stat /= 0) then
         call error%set(ERROR_PARSE, "Error reading natoms in %geometry section")
         return
      end if

      read (line, *, iostat=io_stat) natoms
      if (io_stat /= 0) then
         call error%set(ERROR_PARSE, "Invalid natoms in %geometry section")
         return
      end if

      geom%natoms = natoms

      ! Read blank line (comment line in XYZ format)
      read (unit, '(A)', iostat=io_stat) line
      if (io_stat /= 0) then
         call error%set(ERROR_PARSE, "Error reading comment line in %geometry section")
         return
      end if

      geom%comment = trim(line)

      ! Allocate arrays
      allocate (character(len=4) :: geom%elements(natoms))
      allocate (geom%coords(3, natoms))

      ! Read coordinates
      do i = 1, natoms
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_PARSE, "Error reading geometry coordinates")
            return
         end if

         line = adjustl(line)
         if (trim(strip_comment(line)) == 'end') then
            call error%set(ERROR_PARSE, "Unexpected 'end' while reading geometry")
            return
         end if

         read (line, *, iostat=io_stat) elem, x, y, z
         if (io_stat /= 0) then
            call error%set(ERROR_PARSE, "Invalid coordinate format in %geometry section")
            return
         end if

         geom%elements(i) = trim(elem)
         geom%coords(1, i) = x
         geom%coords(2, i) = y
         geom%coords(3, i) = z
      end do

      ! Read 'end' marker
      read (unit, '(A)', iostat=io_stat) line
      if (io_stat /= 0) then
         call error%set(ERROR_VALIDATION, "Missing 'end' in %geometry section")
         return
      end if

      line = adjustl(line)
      if (trim(strip_comment(line)) /= 'end') then
         call error%set(ERROR_PARSE, "Expected 'end' after geometry coordinates")
         return
      end if

   end subroutine parse_geometry_generic

end submodule mqc_config_parser_structure
