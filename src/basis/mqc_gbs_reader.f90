!! Reader for Gaussian94-format (.gbs) basis set files
module mqc_gbs_reader
   !! Parses Gaussian94-format basis set files as served by the Basis Set
   !! Exchange, producing the same `mqc_cgto` types as the GAMESS reader.
   !!
   !! This is the format used for both orbital and JKFIT auxiliary basis sets,
   !! and is what the cuEST reference samples consume, so results computed
   !! through it are directly comparable with the C reference.
   !!
   !! Restrictions, matching the upstream cuEST helper:
   !!   * Combined SP/SPD shells are not handled. Download with
   !!     "Uncontract SPDF", or split the blocks into separate S and P.
   !!   * ECP definitions following the basis are not parsed.
   !!
   !! Format reminder: element blocks are delimited by `****`, each block opens
   !! with `<symbol> 0`, and each shell header is `<L> <nprim> <scale>` followed
   !! by `nprim` lines of `<exponent> <coefficient>`.
   use pic_types, only: dp
   use mqc_cgto, only: cgto_type, atomic_basis_type, molecular_basis_type
   use mqc_error, only: error_t, ERROR_IO, ERROR_PARSE
   implicit none
   private

   public :: read_gbs_element        !! Parse one element from a .gbs file
   public :: build_molecular_basis_gbs  !! Build a molecular basis from a .gbs file

   integer, parameter :: MAX_LINE = 512  !! Longest .gbs line we accept

   character(len=1), parameter :: SHELL_SYMBOLS(0:6) = &
                                  [character(len=1) :: 'S', 'P', 'D', 'F', 'G', 'H', 'I']
      !! Angular momentum symbols, indexed by L

contains

   pure function upcase(string) result(upper)
      !! Uppercase a string
      character(len=*), intent(in) :: string
      character(len=len(string)) :: upper
      integer :: i, code

      do i = 1, len(string)
         code = iachar(string(i:i))
         if (code >= iachar('a') .and. code <= iachar('z')) then
            upper(i:i) = achar(code - 32)
         else
            upper(i:i) = string(i:i)
         end if
      end do
   end function upcase

   pure function shell_symbol_to_ang_mom(symbol) result(ang_mom)
      !! Convert a shell symbol to an angular momentum quantum number
      !!
      !! Returns -1 for an unrecognized symbol (including the SP/SPD
      !! combined shells this reader deliberately does not support).
      character(len=1), intent(in) :: symbol
      integer :: ang_mom
      integer :: i

      ang_mom = -1
      do i = 0, ubound(SHELL_SYMBOLS, 1)
         if (SHELL_SYMBOLS(i) == upcase(symbol)) then
            ang_mom = i
            return
         end if
      end do
   end function shell_symbol_to_ang_mom

   pure function count_words(line) result(n_words)
      !! Count whitespace-separated words in a line
      character(len=*), intent(in) :: line
      integer :: n_words
      integer :: i
      logical :: in_word

      n_words = 0
      in_word = .false.
      do i = 1, len_trim(line)
         if (line(i:i) == ' ' .or. line(i:i) == achar(9)) then
            in_word = .false.
         else if (.not. in_word) then
            in_word = .true.
            n_words = n_words + 1
         end if
      end do
   end function count_words

   pure function first_token(line) result(token)
      !! First whitespace-separated token of a line, blank if there is none
      character(len=*), intent(in) :: line
      character(len=len(line)) :: token
      integer :: i, j

      token = ''
      i = 1
      do while (i <= len_trim(line))
         if (line(i:i) /= ' ' .and. line(i:i) /= achar(9)) exit
         i = i + 1
      end do
      if (i > len_trim(line)) return

      j = i
      do while (j <= len_trim(line))
         if (line(j:j) == ' ' .or. line(j:j) == achar(9)) exit
         j = j + 1
      end do
      token = line(i:j - 1)
   end function first_token

   pure function is_skippable(line) result(skip)
      !! True for blank lines and `!` comments
      character(len=*), intent(in) :: line
      logical :: skip
      character(len=len(line)) :: trimmed

      trimmed = adjustl(line)
      skip = (len_trim(trimmed) == 0)
      if (.not. skip) skip = (trimmed(1:1) == '!')
   end function is_skippable

   subroutine read_gbs_element(gbs_path, element_symbol, atom_basis, error)
      !! Parse the basis definition for one element from a .gbs file
      !!
      !! Two passes over the file: the first counts shells and primitives so
      !! the arrays can be sized exactly, the second fills them.
      character(len=*), intent(in) :: gbs_path        !! Path to the .gbs file
      character(len=*), intent(in) :: element_symbol  !! Element to extract, e.g. "O"
      type(atomic_basis_type), intent(out) :: atom_basis  !! Parsed shells
      type(error_t), intent(out) :: error

      integer :: unit, iostat, n_words, n_prim, n_to_skip
      integer :: n_shells, ishell, iprim, n_prim_to_parse, ang_mom
      logical :: in_block, found_block
      character(len=MAX_LINE) :: line
      character(len=8) :: wanted, block_symbol, shell_token
      real(dp) :: exponent_value, coefficient_value

      wanted = upcase(adjustl(element_symbol))

      open (newunit=unit, file=gbs_path, status='old', action='read', iostat=iostat)
      if (iostat /= 0) then
         call error%set(ERROR_IO, "Unable to open GBS basis file: "//trim(gbs_path))
         return
      end if

      ! ---- pass 1: count shells ------------------------------------------
      n_shells = 0
      in_block = .false.
      found_block = .false.
      n_to_skip = 0
      do
         read (unit, '(A)', iostat=iostat) line
         if (iostat /= 0) exit
         if (is_skippable(line)) cycle

         ! Primitive lines of a shell we have already counted
         if (n_to_skip > 0) then
            n_to_skip = n_to_skip - 1
            cycle
         end if

         if (index(line, "****") > 0) then
            if (found_block) exit
            in_block = .false.
            cycle
         end if

         n_words = count_words(line)
         if (n_words /= 2 .and. n_words /= 3) cycle

         if (.not. in_block) then
            block_symbol = upcase(first_token(line))
            if (block_symbol == wanted) found_block = .true.
            in_block = .true.
            cycle
         end if
         if (.not. found_block) cycle

         read (line, *, iostat=iostat) shell_token, n_prim
         if (iostat /= 0) then
            close (unit)
            call error%set(ERROR_PARSE, "Malformed shell header in "// &
                           trim(gbs_path)//": "//trim(line))
            return
         end if
         n_to_skip = n_prim
         n_shells = n_shells + 1
      end do

      if (.not. found_block) then
         close (unit)
         call error%set(ERROR_PARSE, "Element "//trim(element_symbol)// &
                        " not found in GBS basis file "//trim(gbs_path))
         return
      end if

      atom_basis%element = trim(adjustl(element_symbol))
      call atom_basis%allocate_shells(n_shells)

      ! ---- pass 2: fill shells -------------------------------------------
      rewind (unit)
      ishell = 0
      iprim = 0
      in_block = .false.
      found_block = .false.
      n_prim_to_parse = 0
      do
         read (unit, '(A)', iostat=iostat) line
         if (iostat /= 0) exit
         if (is_skippable(line)) cycle

         if (n_prim_to_parse > 0 .and. found_block) then
            ! List-directed READ accepts 1.23D-02 natively, so unlike the C
            ! helper there is no D -> e fixup to do here.
            read (line, *, iostat=iostat) exponent_value, coefficient_value
            if (iostat /= 0) then
               close (unit)
               call error%set(ERROR_PARSE, "Malformed primitive in "// &
                              trim(gbs_path)//": "//trim(line))
               return
            end if
            iprim = iprim + 1
            atom_basis%shells(ishell)%exponents(iprim) = exponent_value
            atom_basis%shells(ishell)%coefficients(iprim) = coefficient_value
            n_prim_to_parse = n_prim_to_parse - 1
            cycle
         end if

         if (index(line, "****") > 0) then
            if (found_block) exit
            in_block = .false.
            cycle
         end if

         n_words = count_words(line)
         if (n_words /= 2 .and. n_words /= 3) cycle

         if (.not. in_block) then
            block_symbol = upcase(first_token(line))
            if (block_symbol == wanted) found_block = .true.
            in_block = .true.
            cycle
         end if
         if (.not. found_block) cycle

         read (line, *, iostat=iostat) shell_token, n_prim
         if (iostat /= 0) then
            close (unit)
            call error%set(ERROR_PARSE, "Malformed shell header in "// &
                           trim(gbs_path)//": "//trim(line))
            return
         end if

         ang_mom = shell_symbol_to_ang_mom(shell_token(1:1))
         if (ang_mom < 0 .or. len_trim(shell_token) > 1) then
            close (unit)
            call error%set(ERROR_PARSE, "Unsupported shell type '"// &
                           trim(shell_token)//"' in "//trim(gbs_path)// &
                           " (combined SP/SPD shells must be split)")
            return
         end if

         ishell = ishell + 1
         atom_basis%shells(ishell)%ang_mom = ang_mom
         call atom_basis%shells(ishell)%allocate_arrays(n_prim)
         iprim = 0
         n_prim_to_parse = n_prim
      end do
      close (unit)

      if (ishell /= n_shells) then
         call error%set(ERROR_PARSE, "Inconsistent shell count parsing "// &
                        trim(element_symbol)//" from "//trim(gbs_path))
      end if
   end subroutine read_gbs_element

   subroutine build_molecular_basis_gbs(gbs_path, element_symbols, mol_basis, error)
      !! Build a molecular basis from a .gbs file for a list of atoms
      !!
      !! Each unique element is parsed once and its shells copied to every
      !! atom of that element, in geometry order.
      character(len=*), intent(in) :: gbs_path            !! Path to the .gbs file
      character(len=*), intent(in) :: element_symbols(:)  !! Element per atom, geometry order
      type(molecular_basis_type), intent(out) :: mol_basis
      type(error_t), intent(out) :: error

      type(atomic_basis_type), allocatable :: unique_bases(:)
      character(len=len(element_symbols)), allocatable :: unique_symbols(:)
      integer :: n_atoms, n_unique, iatom, iunique, match_idx
      logical :: is_new

      n_atoms = size(element_symbols)

      ! ---- collect unique elements ----------------------------------------
      allocate (unique_symbols(n_atoms))
      n_unique = 0
      do iatom = 1, n_atoms
         is_new = .true.
         do iunique = 1, n_unique
            if (upcase(adjustl(element_symbols(iatom))) == unique_symbols(iunique)) then
               is_new = .false.
               exit
            end if
         end do
         if (is_new) then
            n_unique = n_unique + 1
            unique_symbols(n_unique) = upcase(adjustl(element_symbols(iatom)))
         end if
      end do

      ! ---- parse each unique element once ---------------------------------
      allocate (unique_bases(n_unique))
      do iunique = 1, n_unique
         call read_gbs_element(gbs_path, unique_symbols(iunique), &
                               unique_bases(iunique), error)
         if (error%has_error()) then
            call error%add_context("mqc_gbs_reader:build_molecular_basis_gbs")
            return
         end if
      end do

      ! ---- assign to every atom -------------------------------------------
      call mol_basis%allocate_elements(n_atoms)
      do iatom = 1, n_atoms
         match_idx = 1
         do iunique = 1, n_unique
            if (upcase(adjustl(element_symbols(iatom))) == unique_symbols(iunique)) then
               match_idx = iunique
               exit
            end if
         end do
         call copy_atomic_basis(unique_bases(match_idx), mol_basis%elements(iatom))
      end do

      do iunique = 1, n_unique
         call unique_bases(iunique)%destroy()
      end do
   end subroutine build_molecular_basis_gbs

   pure subroutine copy_atomic_basis(source, dest)
      !! Deep copy of an atomic basis
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

end module mqc_gbs_reader
