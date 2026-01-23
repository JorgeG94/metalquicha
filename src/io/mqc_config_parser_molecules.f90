!! Molecules section parser for MQC config files
!! Handles: molecules section with multiple molecule blocks
submodule(mqc_config_parser) mqc_config_parser_molecules
   implicit none

contains

   module subroutine parse_molecules_section(unit, config, error)
      !! Parse %molecules section containing multiple %molecule blocks
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value
      character(len=256) :: msg
      integer :: io_stat, eq_pos, nmol, imol
      nmol = 0

      ! First pass: read nmol
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_IO, "Unexpected end of file in %molecules section")
            return
         end if

         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#' .or. line(1:1) == '!') cycle

         if (trim(strip_comment(line)) == 'end') exit

         eq_pos = index(line, '=')
         if (eq_pos > 0) then
            key = adjustl(line(1:eq_pos - 1))
            value = adjustl(line(eq_pos + 1:))

            if (trim(key) == 'nmol') then
               read (value, *, iostat=io_stat) nmol
               if (io_stat /= 0) then
                  call error%set(ERROR_PARSE, "Invalid nmol value")
                  return
               end if
               exit
            end if
         end if
      end do

      if (nmol == 0) then
         ! No molecules, just skip to end
         call skip_to_end(unit, error)
         return
      end if

      config%nmol = nmol
      allocate (config%molecules(nmol))

      ! Parse individual molecules
      imol = 0
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) exit

         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#' .or. line(1:1) == '!') cycle

         if (trim(strip_comment(line)) == 'end') exit

         if (trim(line) == '%molecule') then
            imol = imol + 1
            if (imol > nmol) then
               call error%set(ERROR_PARSE, "More molecules than declared nmol")
               return
            end if
            call parse_single_molecule(unit, config%molecules(imol), error)
            if (error%has_error()) then
               call error%add_context("mqc_config_parser:parse_molecules_section")
               return
            end if
         end if
      end do

      if (imol /= nmol) then
         write (msg, '(A,I0,A,I0)') "Expected ", nmol, " molecules, found ", imol
         call error%set(ERROR_PARSE, trim(msg))
         return
      end if

   end subroutine parse_molecules_section

   subroutine parse_single_molecule(unit, mol, error)
      !! Parse a single %molecule block with its sections
      integer, intent(in) :: unit
      type(molecule_t), intent(inout) :: mol
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value
      integer :: io_stat, eq_pos
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_IO, "Unexpected end of file in %molecule")
            return
         end if

         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#' .or. line(1:1) == '!') cycle

         if (trim(strip_comment(line)) == 'end') exit

         ! Check for key=value pairs (like name)
         eq_pos = index(line, '=')
         if (eq_pos > 0) then
            key = adjustl(line(1:eq_pos - 1))
            value = adjustl(line(eq_pos + 1:))
            if (trim(key) == 'name') then
               mol%name = trim(value)
               cycle
            end if
         end if

         ! Check for subsections
         if (line(1:1) == '%') then
            select case (trim(line))
            case ('%structure')
               call parse_molecule_structure(unit, mol, error)
            case ('%geometry')
               call parse_molecule_geometry(unit, mol, error)
            case ('%fragments')
               call parse_molecule_fragments(unit, mol, error)
            case ('%connectivity')
               call parse_molecule_connectivity(unit, mol, error)
            case default
               ! Skip unknown subsections
               call skip_to_end(unit, error)
            end select

            if (error%has_error()) then
               call error%add_context("mqc_config_parser:parse_single_molecule")
               return
            end if
         end if
      end do

   end subroutine parse_single_molecule

   subroutine parse_molecule_structure(unit, mol, error)
      !! Parse %structure section for a molecule
      integer, intent(in) :: unit
      type(molecule_t), intent(inout) :: mol
      type(error_t), intent(out) :: error

      call parse_structure_generic(unit, mol%charge, mol%multiplicity, error)

   end subroutine parse_molecule_structure

   subroutine parse_molecule_geometry(unit, mol, error)
      !! Parse %geometry section for a molecule
      integer, intent(in) :: unit
      type(molecule_t), intent(inout) :: mol
      type(error_t), intent(out) :: error

      call parse_geometry_generic(unit, mol%geometry, error)

   end subroutine parse_molecule_geometry

   subroutine parse_molecule_fragments(unit, mol, error)
      !! Parse %fragments section for a molecule
      integer, intent(in) :: unit
      type(molecule_t), intent(inout) :: mol
      type(error_t), intent(out) :: error

      call parse_fragments_generic(unit, mol%nfrag, mol%fragments, error)

   end subroutine parse_molecule_fragments

   subroutine parse_molecule_connectivity(unit, mol, error)
      !! Parse %connectivity section for a molecule
      integer, intent(in) :: unit
      type(molecule_t), intent(inout) :: mol
      type(error_t), intent(out) :: error

      call parse_connectivity_generic(unit, mol%nbonds, mol%nbroken, mol%bonds, error)

   end subroutine parse_molecule_connectivity

end submodule mqc_config_parser_molecules
