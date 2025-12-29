!! This file contains basis set reader routines for basis sets
module mqc_basis_reader
   !! Gaussian basis set parser and molecular basis construction
   !!
   !! Provides utilities for parsing Gaussian-type orbital basis sets
   !! from text files and building molecular basis sets for quantum calculations.
   use mqc_cgto, only: cgto_type, atomic_basis_type, molecular_basis_type
   use mqc_basis_file_reader, only: strings_equal
   use mqc_error, only: error_t, ERROR_PARSE
   use pic_types, only: dp
   implicit none
   private

   public :: classify_line        !! Determine basis file line type
   public :: parse_element_basis  !! Parse basis for single element
   public :: build_molecular_basis  !! Build complete molecular basis
   public :: ang_mom_char_to_int  !! Convert angular momentum character to integer
   public :: ang_mom_int_to_char  !! Convert angular momentum integer to character

   ! Basis file line classification constants
   integer, parameter, public :: LINE_UNKNOWN = 0   !! Unrecognized line type
   integer, parameter, public :: LINE_ATOM = 1      !! Element specification line
   integer, parameter, public :: LINE_SHELL = 2     !! Shell definition line
   integer, parameter, public :: LINE_FUNCTION = 3  !! Basis function coefficient line

contains

   pure function ang_mom_char_to_int(ang_mom_char) result(ang_mom)
      !! Convert angular momentum character to integer
      !!
      !! Standard mapping: S=0, P=1, D=2, F=3, G=4, H=5, I=6
      !! Special case: L=-1 (combined S+P shell, requires splitting)
      character(len=1), intent(in) :: ang_mom_char  !! Angular momentum symbol
      integer :: ang_mom  !! Corresponding integer value

      select case (ang_mom_char)
      case ('S')
         ang_mom = 0
      case ('P')
         ang_mom = 1
      case ('D')
         ang_mom = 2
      case ('F')
         ang_mom = 3
      case ('G')
         ang_mom = 4
      case ('H')
         ang_mom = 5
      case ('I')
         ang_mom = 6
      case ('L')
         ang_mom = -1  ! Special case: L shells are split into S+P
      case default
         ang_mom = -1
      end select
   end function ang_mom_char_to_int

   pure function ang_mom_int_to_char(ang_mom) result(ang_mom_char)
      !! Convert angular momentum integer to character
      !!
      !! Inverse mapping: 0=S, 1=P, 2=D, 3=F, 4=G, 5=H, 6=I
      !! Returns '?' for invalid input values.
      integer, intent(in) :: ang_mom  !! Angular momentum quantum number
      character(len=1) :: ang_mom_char  !! Corresponding symbol character

      select case (ang_mom)
      case (0)
         ang_mom_char = 'S'
      case (1)
         ang_mom_char = 'P'
      case (2)
         ang_mom_char = 'D'
      case (3)
         ang_mom_char = 'F'
      case (4)
         ang_mom_char = 'G'
      case (5)
         ang_mom_char = 'H'
      case (6)
         ang_mom_char = 'I'
      case default
         ang_mom_char = '?'
      end select
   end function ang_mom_int_to_char

   pure function classify_line(line) result(line_type)
      !! Classify a line from a gamess formatted basis set file
      character(len=*), intent(in) :: line
      integer :: line_type

      character(len=:), allocatable :: line_trim

      line_trim = trim(adjustl(line))

      if (is_blank_or_control(line_trim)) then
         line_type = LINE_UNKNOWN
      else if (is_function_line(line_trim)) then
         line_type = LINE_FUNCTION
      else if (is_shell_header(line_trim)) then
         line_type = LINE_SHELL
      else
         line_type = LINE_ATOM
      end if

   end function classify_line

   pure function is_blank_or_control(line) result(res)
      !! Check if a line is blank or a control line (starts with '$')
      character(len=*), intent(in) :: line
      logical :: res
      integer :: trimmed_len

      trimmed_len = len_trim(line)

      if (trimmed_len == 0) then
         res = .true.
      else
         res = (line(1:1) == '$')
      end if
   end function is_blank_or_control

   pure function is_function_line(line) result(res)
      !! Check if a line is a function coefficient line (starts with a number)
      character(len=*), intent(in) :: line
      logical :: res
      character(len=1) :: first_char

      if (len_trim(line) == 0) then
         res = .false.
         return
      end if

      first_char = line(1:1)
      res = (first_char >= '0' .and. first_char <= '9')
   end function is_function_line

   pure function is_shell_header(line) result(res)
      !! Check if a line is a shell header line (starts with S, P, D, F, G, H, I, or L)
      character(len=*), intent(in) :: line
      logical :: res
      character(len=1) :: first_char
      integer :: ios, dummy

      res = .false.
      if (len_trim(line) == 0) return

      first_char = line(1:1)

      if (.not. any(first_char == ['S', 'P', 'D', 'F', 'G', 'H', 'I', 'L'])) return

      read (line(2:), *, iostat=ios) dummy
      res = (ios == 0)

   end function is_shell_header

   pure subroutine parse_element_basis(basis_string, element_name, atom_basis, error)
      !! Parse basis set for a specific element from a GAMESS formatted basis string
      character(len=*), intent(in) :: basis_string
      character(len=*), intent(in) :: element_name
      type(atomic_basis_type), intent(out) :: atom_basis
      type(error_t), intent(out) :: error

      integer :: nshells

      ! Pass 1: Find the element and count its shells
      call count_shells_for_element(basis_string, element_name, nshells, error)
      if (error%has_error()) return

      if (nshells == 0) then
         call error%set(ERROR_PARSE, "Element "//trim(element_name)//" not found in basis file")
         return
      end if

      ! ! Allocate shells
      atom_basis%element = trim(element_name)
      call atom_basis%allocate_shells(nshells)

      ! ! Pass 2: Parse and fill shell data
      call fill_element_basis(basis_string, element_name, atom_basis, error)

   end subroutine parse_element_basis

   pure subroutine count_shells_for_element(basis_string, element_name, nshells, error)
      !! Count the number of shells for a specific element in a GAMESS formatted basis string,
      character(len=*), intent(in) :: basis_string
      character(len=*), intent(in) :: element_name
      integer, intent(out) :: nshells
      type(error_t), intent(out) :: error

      integer :: line_start, line_end, line_type
      character(len=256) :: line
      logical :: in_target_element, found_element
      character(len=1) :: ang_mom

      nshells = 0
      in_target_element = .false.
      found_element = .false.
      line_start = 1

      do while (line_start <= len(basis_string))
         call get_next_line(basis_string, line_start, line, line_end)
         if (line_end == 0) exit

         line = adjustl(line)
         line_type = classify_line(line)

         select case (line_type)
         case (LINE_ATOM)
            ! Check if this is our target element
            if (strings_equal(line, element_name)) then
               in_target_element = .true.
               found_element = .true.
            else
               ! Different element - stop counting if we were in target
               if (in_target_element) exit
               in_target_element = .false.
            end if

         case (LINE_SHELL)
            if (in_target_element) then
               ! Extract angular momentum
               line = adjustl(line)
               ang_mom = line(1:1)

               ! L shells become 2 shells (S + P)
               if (ang_mom == 'L') then
                  nshells = nshells + 2
               else
                  nshells = nshells + 1
               end if
            end if

         case (LINE_UNKNOWN)
            ! Skip blank lines and comments
            continue
         case default
            ! Skip any other line types (e.g., LINE_FUNCTION)
            continue
         end select

         line_start = line_end
      end do

      ! Check if we found the element at all
      if (.not. found_element) then
         call error%set(ERROR_PARSE, "Element not found in basis string: "//trim(element_name))
      end if

   end subroutine count_shells_for_element

   pure subroutine get_next_line(string, line_start, line, line_end)
      !! Extract the next line from a string starting at line_start
      character(len=*), intent(in) :: string
      integer, intent(in) :: line_start
      character(len=*), intent(out) :: line
      integer, intent(out) :: line_end

      integer :: newline_pos

      if (line_start > len(string)) then
         line = ''
         line_end = 0
         return
      end if

      newline_pos = index(string(line_start:), new_line('a'))

      if (newline_pos == 0) then
         ! Last line (no newline at end)
         line = string(line_start:)
         line_end = len(string) + 1
      else
         line = string(line_start:line_start + newline_pos - 2)
         line_end = line_start + newline_pos
      end if

   end subroutine get_next_line

   pure subroutine parse_shell_header(line, ang_mom, nfunc, stat)
      !! Parse shell header line (e.g., "S 2" or "L 3")
      character(len=*), intent(in) :: line
      character(len=1), intent(out) :: ang_mom
      integer, intent(out) :: nfunc
      integer, intent(out) :: stat

      character(len=256) :: line_trim

      line_trim = adjustl(line)
      ang_mom = line_trim(1:1)

      ! Read the number of functions
      read (line_trim(2:), *, iostat=stat) nfunc

   end subroutine parse_shell_header

   pure subroutine parse_function_line(line, func_num, exponent, coeff_s, coeff_p, has_p, stat)
      !! Parse function line (e.g., "1 1.0 2.0" or "1 1.0 2.0 3.0" for L shells)
      character(len=*), intent(in) :: line
      integer, intent(out) :: func_num
      real(dp), intent(out) :: exponent
      real(dp), intent(out) :: coeff_s
      real(dp), intent(out), optional :: coeff_p
      logical, intent(out) :: has_p
      integer, intent(out) :: stat

      real(dp) :: temp_p

      has_p = .false.

      ! Try to read 4 values (func_num, exponent, coeff_s, coeff_p)
      read (line, *, iostat=stat) func_num, exponent, coeff_s, temp_p

      if (stat == 0) then
         ! Successfully read 4 values - this is an L shell
         has_p = .true.
         if (present(coeff_p)) coeff_p = temp_p
      else
         ! Try reading just 3 values (func_num, exponent, coeff_s)
         read (line, *, iostat=stat) func_num, exponent, coeff_s
      end if

   end subroutine parse_function_line

   pure subroutine fill_element_basis(basis_string, element_name, atom_basis, error)
      !! Fill in the shell data for a specific element from a GAMESS formatted basis string
      character(len=*), intent(in) :: basis_string
      character(len=*), intent(in) :: element_name
      type(atomic_basis_type), intent(inout) :: atom_basis
      type(error_t), intent(out) :: error

      integer :: line_start, line_end, line_type
      character(len=256) :: line
      logical :: in_data_block, in_target_element
      character(len=1) :: ang_mom
      integer :: nfunc, func_num, ishell, ifunc
      real(dp) :: exponent, coeff_s, coeff_p
      logical :: has_p

      ! L shell handling: we split into two shells, need to track both
      logical :: reading_l_shell
      integer :: l_shell_s_idx, l_shell_p_idx
      integer :: stat

      in_data_block = .false.
      in_target_element = .false.
      ishell = 0
      reading_l_shell = .false.

      line_start = 1
      do while (line_start <= len(basis_string))
         call get_next_line(basis_string, line_start, line, line_end)
         if (line_end == 0) exit

         line = adjustl(line)
         line_type = classify_line(line)

         select case (line_type)
            ! case (LINE_UNKNOWN)
            !   if (index(line, '$DATA') > 0) then
            !     in_data_block = .true.
            !   else if (index(line, '$END') > 0) then
            !     exit
            !   end if

         case (LINE_ATOM)
            if (strings_equal(line, element_name)) then
               in_target_element = .true.
            else
               if (in_target_element) exit
               in_target_element = .false.
            end if

         case (LINE_SHELL)
            if (in_target_element) then
               ! Parse shell header
               call parse_shell_header(line, ang_mom, nfunc, stat)
               if (stat /= 0) then
                  call error%set(ERROR_PARSE, "Failed to parse shell header: "//trim(line))
                  return
               end if

               if (ang_mom == 'L') then
                  ! L shell: create two shells (S and P)
                  reading_l_shell = .true.

                  ishell = ishell + 1
                  l_shell_s_idx = ishell
                  atom_basis%shells(ishell)%ang_mom = 0  ! S
                  call atom_basis%shells(ishell)%allocate_arrays(nfunc)

                  ishell = ishell + 1
                  l_shell_p_idx = ishell
                  atom_basis%shells(ishell)%ang_mom = 1  ! P
                  call atom_basis%shells(ishell)%allocate_arrays(nfunc)

                  ifunc = 0  ! Reset function counter
               else
                  ! Regular shell
                  reading_l_shell = .false.
                  ishell = ishell + 1

                  ! Set angular momentum (S=0, P=1, D=2, F=3, G=4, H=5, I=6)
                  atom_basis%shells(ishell)%ang_mom = ang_mom_char_to_int(ang_mom)

                  call atom_basis%shells(ishell)%allocate_arrays(nfunc)
                  ifunc = 0
               end if
            end if

         case (LINE_FUNCTION)
            if (in_target_element) then
               call parse_function_line(line, func_num, exponent, coeff_s, coeff_p, has_p, stat)
               if (stat /= 0) then
                  call error%set(ERROR_PARSE, "Failed to parse function line: "//trim(line))
                  return
               end if

               ifunc = ifunc + 1

               if (reading_l_shell) then
               if (.not. has_p) then
                  call error%set(ERROR_PARSE, "L shell requires both S and P coefficients")
                  return
               end if
               ! Store in both S and P shells
               atom_basis%shells(l_shell_s_idx)%exponents(ifunc) = exponent
               atom_basis%shells(l_shell_s_idx)%coefficients(ifunc) = coeff_s

               atom_basis%shells(l_shell_p_idx)%exponents(ifunc) = exponent
               atom_basis%shells(l_shell_p_idx)%coefficients(ifunc) = coeff_p
               else
               ! Store in current shell
               atom_basis%shells(ishell)%exponents(ifunc) = exponent
               atom_basis%shells(ishell)%coefficients(ifunc) = coeff_s
               end if
            end if

         case default
            ! Skip unknown line types (e.g., LINE_UNKNOWN, blank lines, comments)
            continue
         end select

         line_start = line_end
      end do

   end subroutine fill_element_basis

   pure subroutine find_unique_strings(input_array, unique_array, nunique)
      !! Find unique strings in an array
      !! Returns array of unique strings and count
      character(len=*), intent(in) :: input_array(:)
      character(len=:), allocatable, intent(out) :: unique_array(:)
      integer, intent(out) :: nunique

      integer :: i, j, n
      logical :: is_unique
      character(len=len(input_array)), allocatable :: temp_unique(:)

      n = size(input_array)
      allocate (temp_unique(n))  ! Max possible size
      nunique = 0

      do i = 1, n
         is_unique = .true.

         ! Check if we've already seen this string
         do j = 1, nunique
            if (strings_equal(input_array(i), temp_unique(j))) then
               is_unique = .false.
               exit
            end if
         end do

         if (is_unique) then
            nunique = nunique + 1
            temp_unique(nunique) = input_array(i)
         end if
      end do

      ! Allocate output array with exact size and copy
      allocate (character(len=len(input_array)) :: unique_array(nunique))
      unique_array = temp_unique(1:nunique)

   end subroutine find_unique_strings

   pure subroutine copy_atomic_basis(source, dest)
      !! Deep copy of atomic basis data from source to dest
      type(atomic_basis_type), intent(in) :: source
      type(atomic_basis_type), intent(out) :: dest
      integer :: ishell

      dest%element = source%element
      call dest%allocate_shells(source%nshells)

      do ishell = 1, source%nshells
         dest%shells(ishell)%ang_mom = source%shells(ishell)%ang_mom
         call dest%shells(ishell)%allocate_arrays(source%shells(ishell)%nfunc)
         dest%shells(ishell)%exponents = source%shells(ishell)%exponents
         dest%shells(ishell)%coefficients = source%shells(ishell)%coefficients
      end do

   end subroutine copy_atomic_basis

   subroutine build_molecular_basis(basis_string, element_names, mol_basis, error)
      !! Build molecular basis from geometry and basis file
      !! Only parses unique elements, then copies basis data to atoms
      character(len=*), intent(in) :: basis_string
      character(len=*), intent(in) :: element_names(:)  !! Element for each atom in geometry order
      type(molecular_basis_type), intent(out) :: mol_basis
      type(error_t), intent(out) :: error

      integer :: iatom, natoms, iunique, nunique
      character(len=:), allocatable :: unique_elements(:)
      type(atomic_basis_type), allocatable :: unique_bases(:)
      integer :: match_idx

      match_idx = 0
      natoms = size(element_names)

      ! Find unique elements
      call find_unique_strings(element_names, unique_elements, nunique)

      print *, "Found ", nunique, " unique elements out of ", natoms, " atoms"

      ! Allocate for unique bases
      allocate (unique_bases(nunique))

      ! Parse basis for each unique element
      do iunique = 1, nunique
         print *, "Parsing basis for: ", trim(unique_elements(iunique))
         call parse_element_basis(basis_string, unique_elements(iunique), &
                                  unique_bases(iunique), error)
         if (error%has_error()) then
            ! Prepend context to error message
            call error%set(ERROR_PARSE, "Failed to parse basis for element "// &
                           trim(unique_elements(iunique))//": "//error%get_message())
            return
         end if
      end do

      ! Allocate molecular basis and assign to each atom
      call mol_basis%allocate_elements(natoms)

      do iatom = 1, natoms
         ! Find which unique element this atom corresponds to
         do iunique = 1, nunique
            if (strings_equal(element_names(iatom), unique_elements(iunique))) then
               match_idx = iunique
               exit
            end if
         end do

         ! Copy the basis data
         call copy_atomic_basis(unique_bases(match_idx), mol_basis%elements(iatom))
      end do

      ! Clean up
      do iunique = 1, nunique
         call unique_bases(iunique)%destroy()
      end do

   end subroutine build_molecular_basis

end module mqc_basis_reader
