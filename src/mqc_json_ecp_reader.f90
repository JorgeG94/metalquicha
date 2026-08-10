!! Reader for Basis Set Exchange JSON effective core potentials
module mqc_json_ecp_reader
   !! Parses the `ecp_potentials` blocks the
   !! [Basis Set Exchange](https://www.basissetexchange.org) ships alongside
   !! orbital basis sets, producing the `mqc_ecp` types.
   !!
   !! The BSE layout for one element is
   !!
   !!     "ecp_electrons": 28,
   !!     "ecp_potentials": [
   !!       { "angular_momentum": [3], "r_exponents": [2],
   !!         "gaussian_exponents": ["3.8431140"],
   !!         "coefficients": [["-12.3169000"]] },
   !!       { "angular_momentum": [0], ... }
   !!     ]
   !!
   !! **The local channel is the one with the highest angular momentum**, not
   !! the one listed first. For every def2 element those happen to coincide --
   !! the order is 3, 0, 1, 2 -- but the physics is defined by the angular
   !! momentum, so that is what this selects on. A tabulation that listed its
   !! channels in ascending order would otherwise silently produce a potential
   !! with the local and projected parts exchanged, which is a wrong answer
   !! rather than an error.
   !!
   !! Exponents and coefficients are stored as strings to preserve every digit,
   !! exactly as in the orbital basis files, and are read back with
   !! list-directed input.
   !!
   !! Not every element in a file has an ECP, and not every basis set has an
   !! ECP file at all. Both are ordinary: a missing element yields an entry
   !! with `has_ecp` false rather than an error, so a molecule of light atoms
   !! and heavy ones needs no special casing.
   use pic_types, only: dp
   use mqc_ecp, only: ecp_shell_type, atomic_ecp_type, molecular_ecp_type
   use mqc_elements, only: element_symbol_to_number
   use mqc_error, only: error_t, ERROR_IO, ERROR_PARSE
   use json_module, only: json_file
   implicit none
   private

   public :: read_json_ecp_element    !! Parse one element's ECP from a BSE file
   public :: build_molecular_ecp_json  !! Build a per-atom ECP set for a molecule

contains

   pure function integer_to_key(value) result(key)
      !! An integer as the string BSE keys its elements by
      integer, intent(in) :: value
      character(len=:), allocatable :: key
      character(len=12) :: buffer

      write (buffer, "(I0)") value
      key = trim(buffer)
   end function integer_to_key

   subroutine read_json_ecp_element(json_path, element_symbol, atom_ecp, error)
      !! Parse the ECP for one element, if the file defines one
      !!
      !! An element the file does not carry is not an error: `atom_ecp` comes
      !! back with `has_ecp` false.
      character(len=*), intent(in) :: json_path
      character(len=*), intent(in) :: element_symbol
      type(atomic_ecp_type), intent(out) :: atom_ecp
      type(error_t), intent(out) :: error

      type(json_file) :: json
      character(len=:), allocatable :: element_key, base_path
      integer :: atomic_number, n_potentials, ipot, local_index, iproj
      integer, allocatable :: channel_l(:)
      logical :: found, file_exists

      atom_ecp%element = trim(adjustl(element_symbol))

      inquire (file=json_path, exist=file_exists)
      if (.not. file_exists) then
         call error%set(ERROR_IO, "JSON ECP file not found: "//trim(json_path))
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
         call error%set(ERROR_PARSE, "Could not parse JSON ECP file: "//trim(json_path))
         call json%destroy()
         return
      end if

      base_path = "elements."//element_key
      call json%info(base_path//".ecp_potentials", found=found, n_children=n_potentials)
      if (.not. found .or. n_potentials <= 0) then
         ! A light element in a file that only covers the heavy ones.
         call json%destroy()
         return
      end if

      call json%get(base_path//".ecp_electrons", atom_ecp%core_electrons, found)
      if (.not. found) then
         call error%set(ERROR_PARSE, "Element "//trim(element_symbol)// &
                        " has ecp_potentials but no ecp_electrons in "//trim(json_path))
         call json%destroy()
         return
      end if

      ! ---- which channel is local? -----------------------------------------
      allocate (channel_l(n_potentials))
      do ipot = 1, n_potentials
         call channel_angular_momentum(json, base_path, ipot, channel_l(ipot), found)
         if (.not. found) then
            call error%set(ERROR_PARSE, "ECP channel without angular_momentum for "// &
                           trim(element_symbol)//" in "//trim(json_path))
            deallocate (channel_l)
            call json%destroy()
            return
         end if
      end do
      local_index = maxloc(channel_l, dim=1)

      call atom_ecp%allocate_projected(n_potentials - 1)
      call read_channel(json, base_path, local_index, channel_l(local_index), &
                        atom_ecp%local, json_path, error)
      if (error%has_error()) then
         deallocate (channel_l)
         call json%destroy()
         return
      end if

      iproj = 0
      do ipot = 1, n_potentials
         if (ipot == local_index) cycle
         iproj = iproj + 1
         call read_channel(json, base_path, ipot, channel_l(ipot), &
                           atom_ecp%projected(iproj), json_path, error)
         if (error%has_error()) then
            deallocate (channel_l)
            call json%destroy()
            return
         end if
      end do

      deallocate (channel_l)
      atom_ecp%has_ecp = .true.
      call json%destroy()
   end subroutine read_json_ecp_element

   subroutine channel_angular_momentum(json, base_path, ipot, ang_mom, found)
      !! The single angular momentum one ECP channel carries
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: base_path
      integer, intent(in) :: ipot
      integer, intent(out) :: ang_mom
      logical, intent(out) :: found

      integer, allocatable :: momenta(:)

      ang_mom = -1
      call json%get(base_path//".ecp_potentials("//integer_to_key(ipot)// &
                    ").angular_momentum", momenta, found)
      if (.not. found .or. .not. allocated(momenta)) then
         found = .false.
         return
      end if
      ! Unlike an orbital shell, an ECP channel carries exactly one l.
      if (size(momenta) < 1) then
         found = .false.
         return
      end if
      ang_mom = momenta(1)
   end subroutine channel_angular_momentum

   subroutine read_channel(json, base_path, ipot, ang_mom, shell, json_path, error)
      !! Read one channel's primitives
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: base_path
      integer, intent(in) :: ipot, ang_mom
      type(ecp_shell_type), intent(inout) :: shell
      character(len=*), intent(in) :: json_path
      type(error_t), intent(out) :: error

      character(len=:), allocatable :: path, value_text
      integer, allocatable :: powers(:)
      integer :: nprim, iprim, read_status
      real(dp) :: value
      logical :: found

      path = base_path//".ecp_potentials("//integer_to_key(ipot)//")"

      call json%get(path//".r_exponents", powers, found)
      if (.not. found .or. .not. allocated(powers)) then
         call error%set(ERROR_PARSE, "ECP channel without r_exponents in "//trim(json_path))
         return
      end if
      nprim = size(powers)

      call json%info(path//".gaussian_exponents", found=found, n_children=iprim)
      if (.not. found .or. iprim /= nprim) then
         call error%set(ERROR_PARSE, "ECP channel with mismatched exponent count in "// &
                        trim(json_path))
         return
      end if

      call shell%allocate_arrays(nprim)
      shell%ang_mom = ang_mom
      shell%radial_powers = powers

      do iprim = 1, nprim
         call json%get(path//".gaussian_exponents("//integer_to_key(iprim)//")", &
                       value_text, found)
         if (.not. found) then
            call error%set(ERROR_PARSE, "Missing ECP exponent in "//trim(json_path))
            return
         end if
         read (value_text, *, iostat=read_status) value
         if (read_status /= 0) then
            call error%set(ERROR_PARSE, "Malformed ECP exponent '"//trim(value_text)// &
                           "' in "//trim(json_path))
            return
         end if
         shell%exponents(iprim) = value

         ! Coefficients are nested one level deeper: a channel has a single
         ! coefficient set, stored as a list of lists to match the orbital
         ! shells, which can have several.
         call json%get(path//".coefficients(1)("//integer_to_key(iprim)//")", &
                       value_text, found)
         if (.not. found) then
            call error%set(ERROR_PARSE, "Missing ECP coefficient in "//trim(json_path))
            return
         end if
         read (value_text, *, iostat=read_status) value
         if (read_status /= 0) then
            call error%set(ERROR_PARSE, "Malformed ECP coefficient '"//trim(value_text)// &
                           "' in "//trim(json_path))
            return
         end if
         shell%coefficients(iprim) = value
      end do
   end subroutine read_channel

   subroutine build_molecular_ecp_json(json_path, element_symbols, mol_ecp, error)
      !! Build a per-atom ECP set for a list of atoms
      !!
      !! Each distinct element is read once and copied to every atom of that
      !! element, as the orbital basis reader does -- parsing the same file per
      !! atom would dominate the setup for anything but a tiny molecule.
      character(len=*), intent(in) :: json_path
      character(len=*), intent(in) :: element_symbols(:)
      type(molecular_ecp_type), intent(out) :: mol_ecp
      type(error_t), intent(out) :: error

      type(atomic_ecp_type), allocatable :: unique_ecps(:)
      integer, allocatable :: unique_z(:)
      integer :: n_atoms, n_unique, iatom, iunique, match, z
      logical :: is_new

      n_atoms = size(element_symbols)
      call mol_ecp%allocate_atoms(n_atoms)
      if (n_atoms == 0) return

      allocate (unique_z(n_atoms), unique_ecps(n_atoms))

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
            call read_json_ecp_element(json_path, element_symbols(iatom), &
                                       unique_ecps(n_unique), error)
            if (error%has_error()) then
               call error%add_context("mqc_json_ecp_reader:build_molecular_ecp_json")
               deallocate (unique_ecps, unique_z)
               return
            end if
         end if
      end do

      mol_ecp%n_ecp_atoms = 0
      do iatom = 1, n_atoms
         z = element_symbol_to_number(trim(adjustl(element_symbols(iatom))))
         match = 1
         do iunique = 1, n_unique
            if (unique_z(iunique) == z) then
               match = iunique
               exit
            end if
         end do
         call copy_atomic_ecp(unique_ecps(match), mol_ecp%atoms(iatom))
         if (mol_ecp%atoms(iatom)%has_ecp) mol_ecp%n_ecp_atoms = mol_ecp%n_ecp_atoms + 1
      end do

      do iunique = 1, n_unique
         call unique_ecps(iunique)%destroy()
      end do
      deallocate (unique_ecps, unique_z)
   end subroutine build_molecular_ecp_json

   subroutine copy_atomic_ecp(source, dest)
      !! Deep copy of one atom's ECP
      type(atomic_ecp_type), intent(in) :: source
      type(atomic_ecp_type), intent(out) :: dest

      integer :: i

      if (allocated(source%element)) dest%element = source%element
      dest%core_electrons = source%core_electrons
      dest%has_ecp = source%has_ecp
      if (.not. source%has_ecp) return

      call copy_ecp_shell(source%local, dest%local)
      call dest%allocate_projected(source%n_projected)
      do i = 1, source%n_projected
         call copy_ecp_shell(source%projected(i), dest%projected(i))
      end do
   end subroutine copy_atomic_ecp

   subroutine copy_ecp_shell(source, dest)
      !! Deep copy of one channel
      type(ecp_shell_type), intent(in) :: source
      type(ecp_shell_type), intent(out) :: dest

      dest%ang_mom = source%ang_mom
      if (source%nprim <= 0) return
      call dest%allocate_arrays(source%nprim)
      dest%radial_powers = source%radial_powers
      dest%exponents = source%exponents
      dest%coefficients = source%coefficients
   end subroutine copy_ecp_shell

end module mqc_json_ecp_reader
