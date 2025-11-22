module pic_physical_fragment
   use pic_types, only: dp, default_int
   use pic_xyz_reader, only: geometry_type, read_xyz_file
   implicit none
   private

   public :: physical_fragment_t, system_geometry_t
   public :: initialize_system_geometry, build_fragment_from_indices
   public :: element_symbol_to_number, element_number_to_symbol

   !! Physical fragment with actual atomic coordinates
   type :: physical_fragment_t
      integer :: n_atoms
      integer, allocatable :: element_numbers(:)     ! Atomic numbers (e.g., 8 for O, 1 for H)
      real(dp), allocatable :: coordinates(:, :)     ! xyz coords (3, n_atoms)
   contains
      procedure :: destroy => fragment_destroy
   end type physical_fragment_t

   !! System geometry holding the full molecular cluster
   type :: system_geometry_t
      integer :: n_monomers                          ! Number of monomers
      integer :: atoms_per_monomer                   ! Atoms in each monomer
      integer :: total_atoms                         ! Total atoms in system
      integer, allocatable :: element_numbers(:)     ! Atomic numbers for all atoms
      real(dp), allocatable :: coordinates(:, :)     ! All coordinates (3, total_atoms)
   contains
      procedure :: destroy => system_destroy
   end type system_geometry_t

contains

   subroutine initialize_system_geometry(full_geom_file, monomer_file, sys_geom, stat, errmsg)
      !! Read full geometry and monomer template, initialize system_geometry_t
      character(len=*), intent(in) :: full_geom_file, monomer_file
      type(system_geometry_t), intent(out) :: sys_geom
      integer, intent(out) :: stat
      character(len=:), allocatable, intent(out) :: errmsg

      type(geometry_type) :: full_geom, monomer_geom
      integer :: i

      ! Read full system geometry
      call read_xyz_file(full_geom_file, full_geom, stat, errmsg)
      if (stat /= 0) return

      ! Read monomer template
      call read_xyz_file(monomer_file, monomer_geom, stat, errmsg)
      if (stat /= 0) then
         call full_geom%destroy()
         return
      end if

      ! Validate that full geometry is a multiple of monomer size
      sys_geom%atoms_per_monomer = monomer_geom%natoms
      sys_geom%total_atoms = full_geom%natoms

      if (mod(sys_geom%total_atoms, sys_geom%atoms_per_monomer) /= 0) then
         stat = 1
         errmsg = "Full geometry atoms not a multiple of monomer atoms"
         call full_geom%destroy()
         call monomer_geom%destroy()
         return
      end if

      sys_geom%n_monomers = sys_geom%total_atoms / sys_geom%atoms_per_monomer

      ! Allocate and copy data
      allocate (sys_geom%element_numbers(sys_geom%total_atoms))
      allocate (sys_geom%coordinates(3, sys_geom%total_atoms))

      ! Convert element symbols to atomic numbers
      do i = 1, sys_geom%total_atoms
         sys_geom%element_numbers(i) = element_symbol_to_number(full_geom%elements(i))
      end do

      sys_geom%coordinates = full_geom%coords

      ! Cleanup
      call full_geom%destroy()
      call monomer_geom%destroy()

      stat = 0

   end subroutine initialize_system_geometry

   subroutine build_fragment_from_indices(sys_geom, monomer_indices, fragment)
      !! Build a fragment on-the-fly from monomer indices
      !! e.g., monomer_indices = [1, 3, 5] extracts waters 1, 3, and 5
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: monomer_indices(:)
      type(physical_fragment_t), intent(out) :: fragment

      integer :: n_monomers_in_frag, atoms_per_monomer
      integer :: i, mono_idx, atom_start, atom_end, frag_atom_idx

      n_monomers_in_frag = size(monomer_indices)
      atoms_per_monomer = sys_geom%atoms_per_monomer

      fragment%n_atoms = n_monomers_in_frag * atoms_per_monomer

      allocate (fragment%element_numbers(fragment%n_atoms))
      allocate (fragment%coordinates(3, fragment%n_atoms))

      frag_atom_idx = 0

      ! Loop over requested monomers and extract their atoms
      do i = 1, n_monomers_in_frag
         mono_idx = monomer_indices(i)
         atom_start = (mono_idx - 1) * atoms_per_monomer + 1
         atom_end = mono_idx * atoms_per_monomer

         ! Copy atoms from this monomer
         fragment%element_numbers(frag_atom_idx + 1:frag_atom_idx + atoms_per_monomer) = &
            sys_geom%element_numbers(atom_start:atom_end)

         fragment%coordinates(:, frag_atom_idx + 1:frag_atom_idx + atoms_per_monomer) = &
            sys_geom%coordinates(:, atom_start:atom_end)

         frag_atom_idx = frag_atom_idx + atoms_per_monomer
      end do

   end subroutine build_fragment_from_indices

   function element_symbol_to_number(symbol) result(atomic_number)
      !! Convert element symbol to atomic number
      character(len=*), intent(in) :: symbol
      integer :: atomic_number

      character(len=2) :: sym

      ! Normalize: uppercase first letter, lowercase second
      sym = adjustl(symbol)
      if (len_trim(sym) >= 1) sym(1:1) = to_upper(sym(1:1))
      if (len_trim(sym) >= 2) sym(2:2) = to_lower(sym(2:2))

      select case (trim(sym))
      case ('H');  atomic_number = 1
      case ('He'); atomic_number = 2
      case ('Li'); atomic_number = 3
      case ('Be'); atomic_number = 4
      case ('B');  atomic_number = 5
      case ('C');  atomic_number = 6
      case ('N');  atomic_number = 7
      case ('O');  atomic_number = 8
      case ('F');  atomic_number = 9
      case ('Ne'); atomic_number = 10
      case ('Na'); atomic_number = 11
      case ('Mg'); atomic_number = 12
      case ('Al'); atomic_number = 13
      case ('Si'); atomic_number = 14
      case ('P');  atomic_number = 15
      case ('S');  atomic_number = 16
      case ('Cl'); atomic_number = 17
      case ('Ar'); atomic_number = 18
      case ('K');  atomic_number = 19
      case ('Ca'); atomic_number = 20
      ! Add more as needed
      case default
         atomic_number = 0  ! Unknown element
      end select

   end function element_symbol_to_number

   function element_number_to_symbol(atomic_number) result(symbol)
      !! Convert atomic number to element symbol
      integer, intent(in) :: atomic_number
      character(len=2) :: symbol

      select case (atomic_number)
      case (1);  symbol = 'H '
      case (2);  symbol = 'He'
      case (3);  symbol = 'Li'
      case (4);  symbol = 'Be'
      case (5);  symbol = 'B '
      case (6);  symbol = 'C '
      case (7);  symbol = 'N '
      case (8);  symbol = 'O '
      case (9);  symbol = 'F '
      case (10); symbol = 'Ne'
      case (11); symbol = 'Na'
      case (12); symbol = 'Mg'
      case (13); symbol = 'Al'
      case (14); symbol = 'Si'
      case (15); symbol = 'P '
      case (16); symbol = 'S '
      case (17); symbol = 'Cl'
      case (18); symbol = 'Ar'
      case (19); symbol = 'K '
      case (20); symbol = 'Ca'
      case default
         symbol = 'Xx'  ! Unknown
      end select

   end function element_number_to_symbol

   pure function to_upper(c) result(uc)
      character(len=1), intent(in) :: c
      character(len=1) :: uc
      integer :: i

      i = iachar(c)
      if (i >= iachar('a') .and. i <= iachar('z')) then
         uc = achar(i - 32)
      else
         uc = c
      end if
   end function to_upper

   pure function to_lower(c) result(lc)
      character(len=1), intent(in) :: c
      character(len=1) :: lc
      integer :: i

      i = iachar(c)
      if (i >= iachar('A') .and. i <= iachar('Z')) then
         lc = achar(i + 32)
      else
         lc = c
      end if
   end function to_lower

   subroutine fragment_destroy(this)
      class(physical_fragment_t), intent(inout) :: this
      if (allocated(this%element_numbers)) deallocate (this%element_numbers)
      if (allocated(this%coordinates)) deallocate (this%coordinates)
      this%n_atoms = 0
   end subroutine fragment_destroy

   subroutine system_destroy(this)
      class(system_geometry_t), intent(inout) :: this
      if (allocated(this%element_numbers)) deallocate (this%element_numbers)
      if (allocated(this%coordinates)) deallocate (this%coordinates)
      this%n_monomers = 0
      this%atoms_per_monomer = 0
      this%total_atoms = 0
   end subroutine system_destroy

end module pic_physical_fragment
