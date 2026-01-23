!! Fragments and connectivity section parsers for MQC config files
!! Handles: fragments, connectivity sections and their generic helpers
submodule(mqc_config_parser) mqc_config_parser_fragments
   implicit none

contains

   module subroutine parse_fragments_section(unit, config, error)
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      call parse_fragments_generic(unit, config%nfrag, config%fragments, error)

   end subroutine parse_fragments_section

   module subroutine parse_connectivity_section(unit, config, error)
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      call parse_connectivity_generic(unit, config%nbonds, config%nbroken, config%bonds, error)

   end subroutine parse_connectivity_section

   module subroutine parse_fragments_generic(unit, nfrag, fragments, error)
      !! Generic parser for %fragments section (works for both config and molecule)
      integer, intent(in) :: unit
      integer, intent(inout) :: nfrag
      type(input_fragment_t), allocatable, intent(inout) :: fragments(:)
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value
      character(len=256) :: msg
      integer :: io_stat, eq_pos, nfrag_local, ifrag
      nfrag_local = 0

      ! First pass: read nfrag
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_IO, "Unexpected end of file in %fragments section")
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

            if (trim(key) == 'nfrag') then
               read (value, *, iostat=io_stat) nfrag_local
               if (io_stat /= 0) then
                  call error%set(ERROR_PARSE, "Invalid nfrag value")
                  return
               end if
               exit
            end if
         end if
      end do

      if (nfrag_local == 0) then
         ! No fragments, just skip to end
         call skip_to_end(unit, error)
         return
      end if

      nfrag = nfrag_local
      allocate (fragments(nfrag))

      ! Parse individual fragments
      ifrag = 0
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) exit

         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#' .or. line(1:1) == '!') cycle

         if (trim(strip_comment(line)) == 'end') exit

         if (trim(line) == '%fragment') then
            ifrag = ifrag + 1
            if (ifrag > nfrag) then
               call error%set(ERROR_PARSE, "More fragments than declared nfrag")
               return
            end if
            call parse_fragment(unit, fragments(ifrag), error)
            if (error%has_error()) then
               call error%add_context("mqc_config_parser:parse_fragments_generic")
               return
            end if
         end if
      end do

      if (ifrag /= nfrag) then
         write (msg, '(A,I0,A,I0)') "Expected ", nfrag, " fragments, found ", ifrag
         call error%set(ERROR_PARSE, trim(msg))
         return
      end if

   end subroutine parse_fragments_generic

   subroutine parse_fragment(unit, fragment, error)
      integer, intent(in) :: unit
      type(input_fragment_t), intent(inout) :: fragment
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value
      integer :: io_stat, eq_pos
      logical :: in_indices
      in_indices = .false.

      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_IO, "Unexpected end of file in %fragment")
            return
         end if

         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#' .or. line(1:1) == '!') cycle

         if (trim(strip_comment(line)) == 'end') then
            if (in_indices) then
               in_indices = .false.
               cycle
            else
               exit
            end if
         end if

         if (trim(line) == '%indices') then
            in_indices = .true.
            cycle
         end if

         if (in_indices) then
            ! Read indices
            call parse_indices_line(line, fragment, error)
            if (error%has_error()) then
               call error%add_context("mqc_config_parser:parse_fragment")
               return
            end if
         else
            eq_pos = index(line, '=')
            if (eq_pos > 0) then
               key = adjustl(line(1:eq_pos - 1))
               value = adjustl(line(eq_pos + 1:))

               select case (trim(key))
               case ('charge')
                  read (value, *, iostat=io_stat) fragment%charge
               case ('multiplicity')
                  read (value, *, iostat=io_stat) fragment%multiplicity
               case default
                  call error%set(ERROR_PARSE, "Unknown key in fragment properties: "//trim(key))
                  return
               end select
            end if
         end if
      end do

   end subroutine parse_fragment

   subroutine parse_indices_line(line, fragment, error)
      character(len=*), intent(in) :: line
      type(input_fragment_t), intent(inout) :: fragment
      type(error_t), intent(out) :: error

      integer :: io_stat, pos, count, idx
      character(len=MAX_LINE_LEN) :: temp_line
      integer, allocatable :: temp_indices(:), new_indices(:)
      temp_line = line

      ! Count how many integers
      count = 0
      do
         read (temp_line, *, iostat=io_stat) idx
         if (io_stat /= 0) exit
         count = count + 1
         ! Remove the read integer from temp_line
         pos = scan(temp_line, ' ')
         if (pos == 0) exit
         temp_line = adjustl(temp_line(pos:))
      end do

      if (count == 0) return

      ! Allocate temporary array
      allocate (temp_indices(count))

      ! Read the integers
      read (line, *, iostat=io_stat) temp_indices
      if (io_stat /= 0) then
         call error%set(ERROR_PARSE, "Error reading fragment indices")
         deallocate (temp_indices)
         return
      end if

      ! Append to existing indices
      if (allocated(fragment%indices)) then
         allocate (new_indices(size(fragment%indices) + count))
         new_indices(1:size(fragment%indices)) = fragment%indices
         new_indices(size(fragment%indices) + 1:) = temp_indices
         call move_alloc(new_indices, fragment%indices)
      else
         call move_alloc(temp_indices, fragment%indices)
      end if

   end subroutine parse_indices_line

   module subroutine parse_connectivity_generic(unit, nbonds, nbroken, bonds, error)
      !! Generic parser for %connectivity section (works for both config and molecule)
      use mqc_physical_fragment, only: bond_t
      integer, intent(in) :: unit
      integer, intent(inout) :: nbonds, nbroken
      type(bond_t), allocatable, intent(inout) :: bonds(:)
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line, key, value, status_str
      integer :: io_stat, eq_pos, nbonds_local, ibond
      integer :: atom_i, atom_j, order
      nbonds_local = 0

      ! First pass: read nbonds
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) then
            call error%set(ERROR_IO, "Unexpected end of file in %connectivity section")
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

            if (trim(key) == 'nbonds') then
               read (value, *, iostat=io_stat) nbonds_local
               if (io_stat /= 0) then
                  call error%set(ERROR_PARSE, "Invalid nbonds value")
                  return
               end if
               exit
            end if
         end if
      end do

      if (nbonds_local == 0) then
         ! No bonds, just skip to end
         call skip_to_end(unit, error)
         return
      end if

      nbonds = nbonds_local
      allocate (bonds(nbonds))

      ! Read bonds
      ibond = 0
      do
         read (unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) exit

         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#' .or. line(1:1) == '!') cycle

         ! Check for key=value pairs (like nbroken=9)
         eq_pos = index(line, '=')
         if (eq_pos > 0) then
            key = adjustl(line(1:eq_pos - 1))
            value = adjustl(line(eq_pos + 1:))
            if (trim(key) == 'nbroken') then
               read (value, *, iostat=io_stat) nbroken
            end if
            cycle
         end if

         if (trim(strip_comment(line)) == 'end') exit

         ! Parse bond line: atom_i atom_j order broken/preserved
         read (line, *, iostat=io_stat) atom_i, atom_j, order, status_str
         if (io_stat /= 0) then
            call error%set(ERROR_PARSE, "Invalid bond format in %connectivity section")
            return
         end if

         ibond = ibond + 1
         if (ibond > nbonds) then
            call error%set(ERROR_PARSE, "More bonds than declared nbonds")
            return
         end if

         bonds(ibond)%atom_i = atom_i
         bonds(ibond)%atom_j = atom_j
         bonds(ibond)%order = order
         bonds(ibond)%is_broken = (trim(status_str) == 'broken')
      end do

   end subroutine parse_connectivity_generic

end submodule mqc_config_parser_fragments
