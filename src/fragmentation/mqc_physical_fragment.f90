module mqc_physical_fragment
   use pic_types, only: dp, default_int
   use mqc_geometry, only: geometry_type
   use mqc_xyz_reader, only: read_xyz_file
   use mqc_elements, only: element_symbol_to_number, element_number_to_symbol, element_mass
   implicit none
   private

   public :: physical_fragment_t, system_geometry_t
   public :: initialize_system_geometry, build_fragment_from_indices
   public :: to_angstrom, to_bohr
   public :: fragment_centroid, fragment_center_of_mass
   public :: distance_between_points, distance_between_fragments
   public :: minimal_distance_between_fragments

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

   pure function fragment_centroid(fragment) result(centroid)
      !! Calculate the geometric centroid (center of geometry) of a fragment
      !! This is the simple average of all atomic coordinates
      !! Returns coordinates in the same units as the fragment (typically Bohr)
      type(physical_fragment_t), intent(in) :: fragment
      real(dp) :: centroid(3)

      integer :: i

      centroid = 0.0_dp

      ! Sum all atomic positions
      do i = 1, fragment%n_atoms
         centroid = centroid + fragment%coordinates(:, i)
      end do

      ! Divide by number of atoms to get average
      centroid = centroid/real(fragment%n_atoms, dp)

   end function fragment_centroid

   pure function fragment_center_of_mass(fragment) result(com)
      !! Calculate the center of mass of a fragment
      !! Weights each atomic position by its atomic mass
      !! Returns coordinates in the same units as the fragment (typically Bohr)
      type(physical_fragment_t), intent(in) :: fragment
      real(dp) :: com(3)

      real(dp) :: total_mass, atom_mass
      integer :: i

      com = 0.0_dp
      total_mass = 0.0_dp

      ! Sum mass-weighted positions and total mass
      do i = 1, fragment%n_atoms
         atom_mass = element_mass(fragment%element_numbers(i))
         com = com + atom_mass*fragment%coordinates(:, i)
         total_mass = total_mass + atom_mass
      end do

      ! Divide by total mass to get center of mass
      com = com/total_mass

   end function fragment_center_of_mass

   pure function distance_between_points(point1, point2) result(distance)
      !! Calculate Euclidean distance between two 3D points
      !! Points should be in the same units (typically Bohr)
      real(dp), intent(in) :: point1(3), point2(3)
      real(dp) :: distance

      real(dp) :: diff(3)

      diff = point2 - point1
      distance = sqrt(dot_product(diff, diff))

   end function distance_between_points

   pure function distance_between_fragments(frag1, frag2, use_com) result(distance)
      !! Calculate distance between two fragments
      !! If use_com is .true., uses center of mass; otherwise uses centroid
      !! Distance is in the same units as the fragment coordinates (typically Bohr)
      type(physical_fragment_t), intent(in) :: frag1, frag2
      logical, intent(in) :: use_com
      real(dp) :: distance

      real(dp) :: point1(3), point2(3)

      if (use_com) then
         point1 = fragment_center_of_mass(frag1)
         point2 = fragment_center_of_mass(frag2)
      else
         point1 = fragment_centroid(frag1)
         point2 = fragment_centroid(frag2)
      end if

      distance = distance_between_points(point1, point2)

   end function distance_between_fragments

   pure function minimal_distance_between_fragments(frag1, frag2) result(min_distance)
      !! Calculate the minimal distance between any two atoms in two fragments
      !! This iterates over all atom pairs and finds the closest pair
      !! Distance is in the same units as the fragment coordinates (typically Bohr)
      type(physical_fragment_t), intent(in) :: frag1, frag2
      real(dp) :: min_distance

      real(dp) :: current_distance
      integer :: i, j

      ! Initialize with a very large value
      min_distance = huge(1.0_dp)

      ! Loop over all atoms in fragment 1
      do i = 1, frag1%n_atoms
         ! Loop over all atoms in fragment 2
         do j = 1, frag2%n_atoms
            ! Calculate distance between atom i in frag1 and atom j in frag2
            current_distance = distance_between_points(frag1%coordinates(:, i), &
                                                       frag2%coordinates(:, j))

            ! Update minimum if this distance is smaller
            if (current_distance < min_distance) then
               min_distance = current_distance
            end if
         end do
      end do

   end function minimal_distance_between_fragments

end module mqc_physical_fragment
