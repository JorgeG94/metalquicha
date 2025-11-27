module mqc_physical_fragment
   use pic_types, only: dp, default_int
   use mqc_geometry, only: geometry_type
   use mqc_xyz_reader, only: read_xyz_file
   use mqc_elements, only: element_symbol_to_number, element_number_to_symbol
   implicit none
   private

   public :: physical_fragment_t, system_geometry_t
   public :: initialize_system_geometry, build_fragment_from_indices
   public :: to_angstrom, to_bohr

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

   ! Bohr radius constant
   real(dp), parameter :: bohr_radius = 0.52917721092_dp

contains

   pure elemental function to_angstrom(bohr_value) result(angstrom_value)
      !! Convert coordinate from Bohr to Angstrom
      real(dp), intent(in) :: bohr_value
      real(dp) :: angstrom_value
      angstrom_value = bohr_value*bohr_radius
   end function to_angstrom

   pure elemental function to_bohr(angstrom_value) result(bohr_value)
      !! Convert coordinate from Angstrom to Bohr
      real(dp), intent(in) :: angstrom_value
      real(dp) :: bohr_value
      bohr_value = angstrom_value/bohr_radius
   end function to_bohr

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

      sys_geom%n_monomers = sys_geom%total_atoms/sys_geom%atoms_per_monomer

      ! Allocate and copy data
      allocate (sys_geom%element_numbers(sys_geom%total_atoms))
      allocate (sys_geom%coordinates(3, sys_geom%total_atoms))

      ! Convert element symbols to atomic numbers
      do i = 1, sys_geom%total_atoms
         sys_geom%element_numbers(i) = element_symbol_to_number(full_geom%elements(i))
      end do

      ! Store coordinates in Bohr (convert from Angstroms)
      sys_geom%coordinates = to_bohr(full_geom%coords)

      ! Cleanup
      call full_geom%destroy()
      call monomer_geom%destroy()

      stat = 0

   end subroutine initialize_system_geometry

   pure subroutine build_fragment_from_indices(sys_geom, monomer_indices, fragment)
      !! Build a fragment on-the-fly from monomer indices
      !! e.g., monomer_indices = [1, 3, 5] extracts waters 1, 3, and 5
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: monomer_indices(:)
      type(physical_fragment_t), intent(out) :: fragment

      integer :: n_monomers_in_frag, atoms_per_monomer
      integer :: i, mono_idx, atom_start, atom_end, frag_atom_idx

      n_monomers_in_frag = size(monomer_indices)
      atoms_per_monomer = sys_geom%atoms_per_monomer

      fragment%n_atoms = n_monomers_in_frag*atoms_per_monomer

      allocate (fragment%element_numbers(fragment%n_atoms))
      allocate (fragment%coordinates(3, fragment%n_atoms))

      frag_atom_idx = 0

      ! Loop over requested monomers and extract their atoms
      do i = 1, n_monomers_in_frag
         mono_idx = monomer_indices(i)
         atom_start = (mono_idx - 1)*atoms_per_monomer + 1
         atom_end = mono_idx*atoms_per_monomer

         ! Copy atoms from this monomer
         fragment%element_numbers(frag_atom_idx + 1:frag_atom_idx + atoms_per_monomer) = &
            sys_geom%element_numbers(atom_start:atom_end)

         fragment%coordinates(:, frag_atom_idx + 1:frag_atom_idx + atoms_per_monomer) = &
            sys_geom%coordinates(:, atom_start:atom_end)

         frag_atom_idx = frag_atom_idx + atoms_per_monomer
      end do

   end subroutine build_fragment_from_indices

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

end module mqc_physical_fragment
