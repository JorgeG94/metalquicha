!! Reader for Basis Set Exchange JSON basis sets
module mqc_json_basis_reader
   !! Parses the JSON format served by the
   !! [Basis Set Exchange](https://www.basissetexchange.org), producing the
   !! same `mqc_cgto` types as the Gaussian94 and GAMESS readers.
   !!
   !! This is the format worth preferring. Unlike Gaussian94, it keeps a
   !! combined shell's coefficient sets separate:
   !!
   !!     "angular_momentum": [0, 1],
   !!     "exponents":    ["...", "...", "..."],
   !!     "coefficients": [[s...], [p...]]
   !!
   !! so an SP shell splits into an S and a P sharing exponents, with no need
   !! for the `uncontract_spdf=1` download flag the `.gbs` reader requires.
   !! It also carries ECP data in the same file, for when that gets wired up.
   !!
   !! Elements are keyed by atomic number as a string, so a symbol is
   !! converted before lookup. Exponents and coefficients are stored as
   !! strings to preserve every digit; they are read with list-directed input,
   !! which accepts the `0.1307093214E+03` form BSE emits.
   use pic_types, only: dp
   use mqc_cgto, only: cgto_type, atomic_basis_type, molecular_basis_type
   use mqc_elements, only: element_symbol_to_number
   use mqc_error, only: error_t, ERROR_IO, ERROR_PARSE
   use json_module, only: json_file
   implicit none
   private

   public :: read_json_basis_element    !! Parse one element from a BSE JSON file
   public :: build_molecular_basis_json  !! Build a molecular basis from a BSE JSON file

contains

   pure function integer_to_key(value) result(key)
      !! Atomic number as the string BSE keys its elements by
      integer, intent(in) :: value
      character(len=:), allocatable :: key
      character(len=12) :: buffer

      write (buffer, "(I0)") value
      key = trim(buffer)
   end function integer_to_key

   subroutine read_json_basis_element(json_path, element_symbol, atom_basis, error)
      !! Parse the basis definition for one element from a BSE JSON file
      !!
      !! Shells listing several angular momenta are split into one shell per
      !! angular momentum, all sharing the exponent list.
      character(len=*), intent(in) :: json_path
      character(len=*), intent(in) :: element_symbol
      type(atomic_basis_type), intent(out) :: atom_basis
      type(error_t), intent(out) :: error

      type(json_file) :: json
      character(len=:), allocatable :: element_key, shell_path
      character(len=:), allocatable :: value_text
      integer, allocatable :: angular_momenta(:)
      integer :: atomic_number, n_shells_json, n_shells_total
      integer :: ishell, imom, iprim, n_prim, shell_index, read_status
      logical :: found
      real(dp) :: value
      logical :: file_exists

      inquire (file=json_path, exist=file_exists)
      if (.not. file_exists) then
         call error%set(ERROR_IO, "JSON basis file not found: "//trim(json_path))
         return
      end if

      atomic_number = element_symbol_to_number(trim(adjustl(element_symbol)))
      if (atomic_number <= 0) then
         call error%set(ERROR_PARSE, "Unrecognized element symbol: "//trim(element_symbol))
         return
      end if
      element_key = integer_to_key(atomic_number)

      call json%initialize()
      call json%load(filename=json_path)
      if (json%failed()) then
         call error%set(ERROR_PARSE, "Could not parse JSON basis file: "//trim(json_path))
         call json%destroy()
         return
      end if

      ! ---- how many shells, and how many after splitting? -------------------
      call json%info("elements."//element_key//".electron_shells", found=found, n_children=n_shells_json)
      if (.not. found .or. n_shells_json <= 0) then
         call error%set(ERROR_PARSE, "Element "//trim(element_symbol)// &
                        " not found in JSON basis file "//trim(json_path))
         call json%destroy()
         return
      end if

      n_shells_total = 0
      do ishell = 1, n_shells_json
         call shell_angular_momenta(json, element_key, ishell, angular_momenta, found)
         if (.not. found) cycle
         n_shells_total = n_shells_total + size(angular_momenta)
      end do

      if (n_shells_total == 0) then
         call error%set(ERROR_PARSE, "No usable shells for "//trim(element_symbol)// &
                        " in "//trim(json_path))
         call json%destroy()
         return
      end if

      atom_basis%element = trim(adjustl(element_symbol))
      call atom_basis%allocate_shells(n_shells_total)

      ! ---- fill, splitting combined shells ----------------------------------
      shell_index = 0
      do ishell = 1, n_shells_json
         call shell_angular_momenta(json, element_key, ishell, angular_momenta, found)
         if (.not. found) cycle

         shell_path = "elements."//element_key//".electron_shells("//integer_to_key(ishell)//")"

         ! Values are read one at a time by indexed path. BSE stores them as
         ! strings to preserve every digit, and json-fortran's scalar string
         ! accessor is the portable way to reach them; basis sets are small
         ! enough that the extra lookups do not matter.
         call json%info(shell_path//".exponents", found=found, n_children=n_prim)
         if (.not. found .or. n_prim <= 0) then
            call error%set(ERROR_PARSE, "Shell without exponents in "//trim(json_path))
            call json%destroy()
            return
         end if

         do imom = 1, size(angular_momenta)
            shell_index = shell_index + 1
            atom_basis%shells(shell_index)%ang_mom = angular_momenta(imom)
            call atom_basis%shells(shell_index)%allocate_arrays(n_prim)

            do iprim = 1, n_prim
               call json%get(shell_path//".exponents("//integer_to_key(iprim)//")", &
                             value_text, found)
               if (.not. found) then
                  call error%set(ERROR_PARSE, "Missing exponent in "//trim(json_path))
                  call json%destroy()
                  return
               end if
               read (value_text, *, iostat=read_status) value
               if (read_status /= 0) then
                  call error%set(ERROR_PARSE, "Malformed exponent '"//trim(value_text)// &
                                 "' in "//trim(json_path))
                  call json%destroy()
                  return
               end if
               atom_basis%shells(shell_index)%exponents(iprim) = value

               call json%get(shell_path//".coefficients("//integer_to_key(imom)//")("// &
                             integer_to_key(iprim)//")", value_text, found)
               if (.not. found) then
                  call error%set(ERROR_PARSE, "Coefficient set does not match the "// &
                                 "exponent count in "//trim(json_path))
                  call json%destroy()
                  return
               end if
               read (value_text, *, iostat=read_status) value
               if (read_status /= 0) then
                  call error%set(ERROR_PARSE, "Malformed coefficient '"//trim(value_text)// &
                                 "' in "//trim(json_path))
                  call json%destroy()
                  return
               end if
               atom_basis%shells(shell_index)%coefficients(iprim) = value
            end do
         end do
      end do

      call json%destroy()
   end subroutine read_json_basis_element

   subroutine shell_angular_momenta(json, element_key, ishell, angular_momenta, found)
      !! Angular momenta carried by one shell entry
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: element_key
      integer, intent(in) :: ishell
      integer, allocatable, intent(out) :: angular_momenta(:)
      logical, intent(out) :: found

      character(len=:), allocatable :: path
      character(len=12) :: index_text

      write (index_text, "(I0)") ishell
      path = "elements."//element_key//".electron_shells("//trim(index_text)//").angular_momentum"
      call json%get(path, angular_momenta, found)
   end subroutine shell_angular_momenta

   subroutine build_molecular_basis_json(json_path, element_symbols, mol_basis, error)
      !! Build a molecular basis from a BSE JSON file for a list of atoms
      character(len=*), intent(in) :: json_path
      character(len=*), intent(in) :: element_symbols(:)
      type(molecular_basis_type), intent(out) :: mol_basis
      type(error_t), intent(out) :: error

      type(atomic_basis_type), allocatable :: unique_bases(:)
      integer, allocatable :: unique_z(:)
      integer :: n_atoms, n_unique, iatom, iunique, match, z
      logical :: is_new

      n_atoms = size(element_symbols)
      allocate (unique_z(n_atoms), unique_bases(n_atoms))

      n_unique = 0
      do iatom = 1, n_atoms
         z = element_symbol_to_number(trim(adjustl(element_symbols(iatom))))
         is_new = .true.
         do iunique = 1, n_unique
            if (unique_z(iunique) == z) then
               is_new = .false.
               exit
            end if
         end do
         if (is_new) then
            n_unique = n_unique + 1
            unique_z(n_unique) = z
            call read_json_basis_element(json_path, element_symbols(iatom), &
                                         unique_bases(n_unique), error)
            if (error%has_error()) then
               call error%add_context("mqc_json_basis_reader:build_molecular_basis_json")
               return
            end if
         end if
      end do

      call mol_basis%allocate_elements(n_atoms)
      do iatom = 1, n_atoms
         z = element_symbol_to_number(trim(adjustl(element_symbols(iatom))))
         match = 1
         do iunique = 1, n_unique
            if (unique_z(iunique) == z) then
               match = iunique
               exit
            end if
         end do
         call copy_atomic_basis(unique_bases(match), mol_basis%elements(iatom))
      end do

      do iunique = 1, n_unique
         call unique_bases(iunique)%destroy()
      end do
      deallocate (unique_bases, unique_z)
   end subroutine build_molecular_basis_json

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

end module mqc_json_basis_reader
