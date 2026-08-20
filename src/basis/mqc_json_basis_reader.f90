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
   !!
   !! Every shell also carries a `function_type` -- `gto`, `gto_spherical` or
   !! `gto_cartesian` -- and it is not decoration. 6-31G* specifies Cartesian d
   !! functions and BSE marks them so; a Cartesian d shell carries six
   !! functions where a spherical one carries five. Reading one as the other
   !! computes a different basis under the same name, worth 1.4 mHartree on
   !! water, with nothing to say it happened. The reader records the convention
   !! on the basis and the integral layer routes on it.
   use pic_types, only: dp
   use mqc_cgto, only: cgto_type, atomic_basis_type, molecular_basis_type, &
                       ANGULAR_FORM_UNSET, ANGULAR_FORM_SPHERICAL, ANGULAR_FORM_CARTESIAN
   use mqc_elements, only: element_symbol_to_number
   use mqc_error, only: error_t, ERROR_IO, ERROR_PARSE
   use json_module, only: json_file
   use mqc_string_utils, only: int_to_text
   implicit none
   private

   public :: read_json_basis_element    !! Parse one element from a BSE JSON file
   public :: build_molecular_basis_json  !! Build a molecular basis from a BSE JSON file

   !> A basis file, parsed once and kept
   type :: cached_basis_t
      character(len=:), allocatable :: path
      type(json_file), pointer :: json => null()
         !! Heap-allocated, and a pointer on purpose. `json_file` has a
         !! finaliser, so holding one by value would mean that growing the cache
         !! array -- copy into a bigger array, `move_alloc`, old array finalised
         !! -- destroys the parsed tree that the surviving copy still points at.
         !! Every basis loaded before the second one would go stale the moment
         !! the second arrived. Storing pointers means growth copies addresses
         !! and nothing is finalised.
   end type cached_basis_t

   !> Every basis file read so far, for the life of the process
   !>
   !> Parsing was happening once per *element*, not once per file, so a water
   !> molecule re-read and re-parsed the whole basis set twice and a five-element
   !> system five times. Anything that builds molecules in a loop -- a fragment
   !> method especially, which builds one per fragment per iteration -- paid that
   !> over and over, and on a shared filesystem with many ranks it is the same
   !> file being opened by everyone repeatedly.
   !>
   !> A run touches one or two basis sets, so this never holds much. Entries are
   !> never evicted: the whole point is that the second request is free, and a
   !> parsed basis set is small next to anything it is used to compute.
   !>
   !> Not thread safe, and does not need to be -- basis sets are read while
   !> setting a molecule up, not from inside a parallel region.
   type(cached_basis_t), allocatable, target :: basis_cache(:)

contains

   subroutine cached_basis_json(json_path, json, error)
      !! The parsed basis file, loading it the first time it is asked for
      character(len=*), intent(in) :: json_path
      type(json_file), pointer, intent(out) :: json
      type(error_t), intent(out) :: error

      type(cached_basis_t), allocatable :: bigger(:)
      logical :: file_exists
      integer :: i, n

      nullify (json)
      if (allocated(basis_cache)) then
         do i = 1, size(basis_cache)
            if (basis_cache(i)%path == json_path) then
               json => basis_cache(i)%json
               ! A failed lookup earlier leaves the exception flag set, which
               ! would make every later query on this file look like a failure.
               call json%clear_exceptions()
               return
            end if
         end do
      end if

      inquire (file=json_path, exist=file_exists)
      if (.not. file_exists) then
         call error%set(ERROR_IO, "JSON basis file not found: "//trim(json_path))
         return
      end if

      n = 0
      if (allocated(basis_cache)) n = size(basis_cache)
      allocate (bigger(n + 1))
      if (n > 0) bigger(1:n) = basis_cache
      call move_alloc(bigger, basis_cache)

      basis_cache(n + 1)%path = json_path
      allocate (basis_cache(n + 1)%json)
      call basis_cache(n + 1)%json%initialize()
      call basis_cache(n + 1)%json%load(filename=json_path)
      if (basis_cache(n + 1)%json%failed()) then
         call error%set(ERROR_PARSE, "Could not parse JSON basis file: "//trim(json_path))
         return
      end if
      json => basis_cache(n + 1)%json
   end subroutine cached_basis_json

   subroutine read_json_basis_element(json_path, element_symbol, atom_basis, error)
      !! Parse the basis definition for one element from a BSE JSON file
      !!
      !! Shells listing several angular momenta are split into one shell per
      !! angular momentum, all sharing the exponent list.
      character(len=*), intent(in) :: json_path
      character(len=*), intent(in) :: element_symbol
      type(atomic_basis_type), intent(out) :: atom_basis
      type(error_t), intent(out) :: error

      type(json_file), pointer :: json
      character(len=:), allocatable :: element_key, shell_path
      character(len=:), allocatable :: value_text
      integer, allocatable :: angular_momenta(:)
      integer :: atomic_number, n_shells_json, n_shells_total, n_contract
      integer :: ishell, imom, iprim, n_prim, shell_index, read_status
      integer :: shell_form
      logical :: found
      real(dp) :: value

      logical :: saw_cartesian, saw_spherical

      atomic_number = element_symbol_to_number(trim(adjustl(element_symbol)))
      if (atomic_number <= 0) then
         call error%set(ERROR_PARSE, "Unrecognized element symbol: "//trim(element_symbol))
         return
      end if
      element_key = int_to_text(atomic_number)

      call cached_basis_json(json_path, json, error)
      if (error%has_error()) return

      ! ---- how many shells, and how many after splitting? -------------------
      call json%info("elements."//element_key//".electron_shells", found=found, n_children=n_shells_json)
      if (.not. found .or. n_shells_json <= 0) then
         call error%set(ERROR_PARSE, "Element "//trim(element_symbol)// &
                        " not found in JSON basis file "//trim(json_path))
         return
      end if

      n_shells_total = 0
      saw_cartesian = .false.
      saw_spherical = .false.
      do ishell = 1, n_shells_json
         call shell_angular_momenta(json, element_key, ishell, angular_momenta, found)
         if (.not. found) cycle
         n_shells_total = n_shells_total + shell_contractions(json, element_key, ishell, &
                                                              size(angular_momenta))

         shell_form = shell_angular_form(json, element_key, ishell, angular_momenta)
         if (shell_form == ANGULAR_FORM_CARTESIAN) saw_cartesian = .true.
         if (shell_form == ANGULAR_FORM_SPHERICAL) saw_spherical = .true.
      end do

      if (n_shells_total == 0) then
         call error%set(ERROR_PARSE, "No usable shells for "//trim(element_symbol)// &
                        " in "//trim(json_path))
         return
      end if

      ! One element, two conventions. libcint chooses spherical or Cartesian
      ! per *call*, not per shell, so there is no way to honour both -- and
      ! silently picking one would be the bug this reader was fixed for, only
      ! quieter. 6-31G* does this on Sc through Zn, where the d shell is marked
      ! Cartesian and the f shell spherical.
      if (saw_cartesian .and. saw_spherical) then
         call error%set(ERROR_PARSE, "Element "//trim(element_symbol)//" in "// &
                        trim(json_path)//" mixes Cartesian and spherical shells "// &
                        "(function_type 'gto_cartesian' on one shell above p and "// &
                        "'gto_spherical' on another). The integrals are built in one "// &
                        "form or the other for the whole molecule, so this basis "// &
                        "cannot be used for this element; choose a set that is "// &
                        "consistently one or the other.")
         return
      end if

      atom_basis%element = trim(adjustl(element_symbol))
      call atom_basis%allocate_shells(n_shells_total)
      if (saw_cartesian) then
         atom_basis%angular_form = ANGULAR_FORM_CARTESIAN
      else if (saw_spherical) then
         atom_basis%angular_form = ANGULAR_FORM_SPHERICAL
      else
         atom_basis%angular_form = ANGULAR_FORM_UNSET
      end if

      ! ---- fill, splitting combined shells ----------------------------------
      shell_index = 0
      do ishell = 1, n_shells_json
         call shell_angular_momenta(json, element_key, ishell, angular_momenta, found)
         if (.not. found) cycle

         shell_path = "elements."//element_key//".electron_shells("//int_to_text(ishell)//")"

         ! Values are read one at a time by indexed path. BSE stores them as
         ! strings to preserve every digit, and json-fortran's scalar string
         ! accessor is the portable way to reach them; basis sets are small
         ! enough that the extra lookups do not matter.
         call json%info(shell_path//".exponents", found=found, n_children=n_prim)
         if (.not. found .or. n_prim <= 0) then
            call error%set(ERROR_PARSE, "Shell without exponents in "//trim(json_path))
            return
         end if

         ! One shell per *coefficient column*, not per angular momentum. The
         ! two coincide for an SP shell, where each column carries its own l --
         ! which is why 6-31G was right. They do not for a general contraction:
         ! cc-pVDZ oxygen is nine s primitives with three columns and a single
         ! l, so looping over angular momenta read column one and discarded the
         ! other two. That gave oxygen nine basis functions instead of
         ! fourteen, and an SCF that converged a full Hartree high.
         n_contract = shell_contractions(json, element_key, ishell, size(angular_momenta))
         do imom = 1, n_contract
            shell_index = shell_index + 1
            ! Each column has its own l when the file lists one per column;
            ! a general contraction lists one l for all of them.
            if (size(angular_momenta) == n_contract) then
               atom_basis%shells(shell_index)%ang_mom = angular_momenta(imom)
            else
               atom_basis%shells(shell_index)%ang_mom = angular_momenta(1)
            end if
            call atom_basis%shells(shell_index)%allocate_arrays(n_prim)

            do iprim = 1, n_prim
               call json%get(shell_path//".exponents("//int_to_text(iprim)//")", &
                             value_text, found)
               if (.not. found) then
                  call error%set(ERROR_PARSE, "Missing exponent in "//trim(json_path))
                  return
               end if
               read (value_text, *, iostat=read_status) value
               if (read_status /= 0) then
                  call error%set(ERROR_PARSE, "Malformed exponent '"//trim(value_text)// &
                                 "' in "//trim(json_path))
                  return
               end if
               atom_basis%shells(shell_index)%exponents(iprim) = value

               call json%get(shell_path//".coefficients("//int_to_text(imom)//")("// &
                             int_to_text(iprim)//")", value_text, found)
               if (.not. found) then
                  call error%set(ERROR_PARSE, "Coefficient set does not match the "// &
                                 "exponent count in "//trim(json_path))
                  return
               end if
               read (value_text, *, iostat=read_status) value
               if (read_status /= 0) then
                  call error%set(ERROR_PARSE, "Malformed coefficient '"//trim(value_text)// &
                                 "' in "//trim(json_path))
                  return
               end if
               atom_basis%shells(shell_index)%coefficients(iprim) = value
            end do
         end do
      end do

   end subroutine read_json_basis_element

   function shell_angular_form(json, element_key, ishell, angular_momenta) result(form)
      !! What one shell entry's `function_type` says about the angular form
      !!
      !! Returns `ANGULAR_FORM_UNSET` for anything at or below p, whatever the
      !! file calls it. An s shell is one function and a p shell is three in
      !! both conventions, so those entries carry no information -- and BSE
      !! labels them the neutral `gto` for exactly that reason. Only d and
      !! above distinguish 6d from 5d, and only those are evidence.
      !!
      !! Anything that is not spelled `gto_cartesian` reads as spherical. That
      !! is the same default the reader had before it looked at the key at all,
      !! so an unfamiliar spelling behaves as it always did rather than
      !! flipping the convention on a typo.
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: element_key
      integer, intent(in) :: ishell
      integer, intent(in) :: angular_momenta(:)
      integer :: form

      character(len=:), allocatable :: function_type
      logical :: found

      form = ANGULAR_FORM_UNSET
      if (size(angular_momenta) == 0) return
      if (maxval(angular_momenta) < 2) return

      call json%get("elements."//element_key//".electron_shells("// &
                    int_to_text(ishell)//").function_type", function_type, found)
      if (.not. found) then
         form = ANGULAR_FORM_SPHERICAL
      else if (function_type == "gto_cartesian") then
         form = ANGULAR_FORM_CARTESIAN
      else
         form = ANGULAR_FORM_SPHERICAL
      end if
   end function shell_angular_form

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
      integer :: molecule_form
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

      ! One convention for the molecule, from the atoms that are actually in
      ! it. Hydrogen in 6-31G* has no shell above s and so votes for nothing;
      ! oxygen has a Cartesian d and settles it. An element that genuinely
      ! disagrees with another is refused rather than resolved -- the integrals
      ! are built one way for every shell at once.
      molecule_form = ANGULAR_FORM_UNSET
      do iunique = 1, n_unique
         if (unique_bases(iunique)%angular_form == ANGULAR_FORM_UNSET) cycle
         if (molecule_form == ANGULAR_FORM_UNSET) then
            molecule_form = unique_bases(iunique)%angular_form
         else if (molecule_form /= unique_bases(iunique)%angular_form) then
            call error%set(ERROR_PARSE, "the basis in "//trim(json_path)// &
                           " is Cartesian for some of this molecule's elements and "// &
                           "spherical for others, and the integrals must be built in "// &
                           "one form for the whole molecule. Choose a basis that is "// &
                           "consistently one or the other.")
            call error%add_context("mqc_json_basis_reader:build_molecular_basis_json")
            return
         end if
      end do

      call mol_basis%allocate_elements(n_atoms)
      mol_basis%angular_form = molecule_form
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

   function shell_contractions(json, element_key, ishell, n_angular) result(n)
      !! How many contracted functions one shell entry defines
      !!
      !! The length of the `coefficients` array, which is the number of
      !! columns -- one per contracted function sharing this shell's
      !! exponents. Falls back to the angular momentum count if the key
      !! cannot be read, which keeps a malformed file behaving as it did
      !! rather than silently producing no shells at all.
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: element_key
      integer, intent(in) :: ishell, n_angular
      integer :: n

      logical :: found
      integer :: n_children

      call json%info("elements."//element_key//".electron_shells("// &
                     int_to_text(ishell)//").coefficients", &
                     found=found, n_children=n_children)
      if (found .and. n_children > 0) then
         n = n_children
      else
         n = n_angular
      end if
   end function shell_contractions

   pure subroutine copy_atomic_basis(source, dest)
      !! Deep copy of an atomic basis
      type(atomic_basis_type), intent(in) :: source
      type(atomic_basis_type), intent(out) :: dest
      integer :: ishell

      dest%element = source%element
      dest%angular_form = source%angular_form
      call dest%allocate_shells(source%nshells)
      do ishell = 1, source%nshells
         dest%shells(ishell)%ang_mom = source%shells(ishell)%ang_mom
         call dest%shells(ishell)%allocate_arrays(source%shells(ishell)%nfunc)
         dest%shells(ishell)%exponents = source%shells(ishell)%exponents
         dest%shells(ishell)%coefficients = source%shells(ishell)%coefficients
      end do
   end subroutine copy_atomic_basis

end module mqc_json_basis_reader
