!! This file contains all routines and types to represent a "physical" fragment or molecule
!! i.e., with atomic coordinates, element types, electronic properties, etc.
module mqc_physical_fragment
   !! Physical molecular fragment representation and geometry handling
   !!
   !! Provides data structures and utilities for managing molecular fragments
   !! with atomic coordinates, electronic properties, and geometric operations.
   use pic_types, only: dp, default_int
   use mqc_geometry, only: geometry_type
   use mqc_xyz_reader, only: read_xyz_file
   use mqc_elements, only: element_symbol_to_number, element_number_to_symbol, element_mass
   use mqc_cgto, only: molecular_basis_type
   use mqc_config_parser, only: bond_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: physical_fragment_t         !! Single molecular fragment type
   public :: system_geometry_t          !! Complete system geometry type
   public :: initialize_system_geometry  !! System geometry initialization
   public :: build_fragment_from_indices  !! Extract fragment from system
   public :: build_fragment_from_atom_list  !! Build fragment from explicit atom indices (for intersections)
   public :: to_angstrom, to_bohr       !! Unit conversion utilities
   public :: fragment_centroid          !! Geometric centroid calculation
   public :: fragment_center_of_mass    !! Mass-weighted center calculation
   public :: distance_between_points    !! Point-to-point distance
   public :: distance_between_fragments  !! Inter-fragment distance
   public :: minimal_distance_between_fragments  !! Closest approach distance

   type :: physical_fragment_t
      !! Physical molecular fragment with atomic coordinates and properties
      !!
      !! Represents a molecular fragment containing atomic positions, element types,
      !! electronic structure information, and basis set data for quantum calculations.
      integer :: n_atoms  !! Number of atoms in this fragment
      integer, allocatable :: element_numbers(:)     !! Atomic numbers (Z values)
      real(dp), allocatable :: coordinates(:, :)     !! Cartesian coordinates (3, n_atoms) in Bohr

      ! Electronic structure properties
      integer :: charge = 0        !! Net molecular charge (electrons)
      integer :: multiplicity = 1  !! Spin multiplicity (2S+1)
      integer :: nelec = 0         !! Total number of electrons

      ! Hydrogen capping for broken bonds
      integer :: n_caps = 0  !! Number of hydrogen caps added (always at end of atom list)
      integer, allocatable :: cap_replaces_atom(:)  !! Original atom index that each cap replaces (size: n_caps)

      ! Quantum chemistry basis set
      type(molecular_basis_type), allocatable :: basis  !! Gaussian basis functions
   contains
      procedure :: destroy => fragment_destroy          !! Memory cleanup
      procedure :: compute_nelec => fragment_compute_nelec  !! Calculate electron count
      procedure :: set_basis => fragment_set_basis      !! Assign basis set
   end type physical_fragment_t

   type :: system_geometry_t
      !! Complete molecular system geometry for fragment-based calculations
      !!
      !! Contains the full atomic structure of a molecular cluster organized
      !! by monomers for efficient fragment generation and MBE calculations.
      integer :: n_monomers        !! Number of monomer units in system
      integer :: atoms_per_monomer  !! Atoms in each monomer (0 if variable-sized)
      integer :: total_atoms       !! Total number of atoms
      integer, allocatable :: element_numbers(:)  !! Atomic numbers for all atoms
      real(dp), allocatable :: coordinates(:, :)  !! All coordinates (3, total_atoms) in Bohr

      ! Electronic structure properties
      integer :: charge         !! Net molecular charge (electrons)
      integer :: multiplicity   !! Spin multiplicity (2S+1)

      ! For variable-sized fragments (explicit fragment definitions)
      integer, allocatable :: fragment_sizes(:)      !! Number of atoms in each fragment (n_monomers)
      integer, allocatable :: fragment_atoms(:, :)   !! Atom indices for each fragment (max_frag_size, n_monomers), 0-indexed
      integer, allocatable :: fragment_charges(:)    !! Charge for each fragment (n_monomers)
      integer, allocatable :: fragment_multiplicities(:)  !! Multiplicity for each fragment (n_monomers)
   contains
      procedure :: destroy => system_destroy  !! Memory cleanup
   end type system_geometry_t

   ! Physical constants
   real(dp), parameter :: bohr_radius = 0.52917721092_dp  !! Bohr radius in Ångström

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

   subroutine initialize_system_geometry(full_geom_file, monomer_file, sys_geom, error)
      !! Read full geometry and monomer template, initialize system_geometry_t
      character(len=*), intent(in) :: full_geom_file, monomer_file
      type(system_geometry_t), intent(out) :: sys_geom
      type(error_t), intent(out) :: error

      type(geometry_type) :: full_geom, monomer_geom
      integer :: i

      call read_xyz_file(full_geom_file, full_geom, error)
      if (error%has_error()) return

      ! Read monomer template
      ! this will be changed once we have a proper input file parsing
      call read_xyz_file(monomer_file, monomer_geom, error)
      if (error%has_error()) then
         call full_geom%destroy()
         return
      end if

      ! Validate that full geometry is a multiple of monomer size
      sys_geom%atoms_per_monomer = monomer_geom%natoms
      sys_geom%total_atoms = full_geom%natoms

      if (mod(sys_geom%total_atoms, sys_geom%atoms_per_monomer) /= 0) then
         call error%set(ERROR_VALIDATION, "Full geometry atoms not a multiple of monomer atoms")
         call full_geom%destroy()
         call monomer_geom%destroy()
         return
      end if

      sys_geom%n_monomers = sys_geom%total_atoms/sys_geom%atoms_per_monomer

      ! TODO JORGE: this can be a sys_geom%allocate()
      allocate (sys_geom%element_numbers(sys_geom%total_atoms))
      allocate (sys_geom%coordinates(3, sys_geom%total_atoms))

      do i = 1, sys_geom%total_atoms
         sys_geom%element_numbers(i) = element_symbol_to_number(full_geom%elements(i))
      end do

      ! Store coordinates in Bohr (convert from Angstroms)
      ! TODO JORGE: need a way to handle units
      sys_geom%coordinates = to_bohr(full_geom%coords)

      call full_geom%destroy()
      call monomer_geom%destroy()

   end subroutine initialize_system_geometry

   subroutine build_fragment_from_indices(sys_geom, monomer_indices, fragment, bonds)
      !! Build a fragment on-the-fly from monomer indices with hydrogen capping for broken bonds
      !!
      !! Extracts atoms from specified monomers and adds hydrogen caps where bonds are broken.
      !! Caps are always added at the end of the atom list.
      !! Supports both fixed-size (identical monomers) and variable-sized fragments.
      !!
      !! Example: monomer_indices = [1, 3, 5] extracts waters 1, 3, and 5
      !!          If connectivity shows broken bonds, hydrogens are capped at positions of missing atoms
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: monomer_indices(:)
      type(physical_fragment_t), intent(out) :: fragment
      type(bond_t), intent(in), optional :: bonds(:)  !! Connectivity information for capping

      integer :: n_monomers_in_frag, atoms_per_monomer, n_atoms_no_caps
      integer :: i, j, mono_idx, atom_start, atom_end, frag_atom_idx
      integer :: ibond, atom_i, atom_j, n_caps, cap_idx
      logical :: atom_i_in_frag, atom_j_in_frag
      integer, allocatable :: atoms_in_fragment(:)  !! List of all atom indices in this fragment
      integer :: iatom, atom_global_idx
      logical :: use_explicit_fragments

      n_monomers_in_frag = size(monomer_indices)

      ! Determine if we're using explicit fragment definitions or regular monomer-based
      use_explicit_fragments = allocated(sys_geom%fragment_atoms)

      if (use_explicit_fragments) then
         ! Variable-sized fragments: count total atoms from fragment definitions
         n_atoms_no_caps = 0
         do i = 1, n_monomers_in_frag
            mono_idx = monomer_indices(i)
            n_atoms_no_caps = n_atoms_no_caps + sys_geom%fragment_sizes(mono_idx)
         end do

         ! Build list of atom indices (0-indexed) from explicit fragment definitions
         allocate (atoms_in_fragment(n_atoms_no_caps))
         iatom = 0
         do i = 1, n_monomers_in_frag
            mono_idx = monomer_indices(i)
            do j = 1, sys_geom%fragment_sizes(mono_idx)
               iatom = iatom + 1
               atoms_in_fragment(iatom) = sys_geom%fragment_atoms(j, mono_idx)
            end do
         end do
      else
         ! Fixed-size monomers: use atoms_per_monomer
         atoms_per_monomer = sys_geom%atoms_per_monomer
         n_atoms_no_caps = n_monomers_in_frag*atoms_per_monomer

         ! Build list of atom indices in this fragment (0-indexed to match bond indices)
         allocate (atoms_in_fragment(n_atoms_no_caps))
         iatom = 0
         do i = 1, n_monomers_in_frag
            mono_idx = monomer_indices(i)
            atom_start = (mono_idx - 1)*atoms_per_monomer
            do atom_i = 0, atoms_per_monomer - 1
               iatom = iatom + 1
               atoms_in_fragment(iatom) = atom_start + atom_i
            end do
         end do
      end if

      ! Count how many caps we need
      n_caps = 0
      if (present(bonds)) then
         do ibond = 1, size(bonds)
            if (.not. bonds(ibond)%is_broken) cycle

            ! Check if exactly one atom of this bond is in the fragment
            atom_i_in_frag = any(atoms_in_fragment == bonds(ibond)%atom_i)
            atom_j_in_frag = any(atoms_in_fragment == bonds(ibond)%atom_j)

            ! Add cap only if one atom in fragment, other not (XOR condition)
            if ((atom_i_in_frag .and. .not. atom_j_in_frag) .or. &
                (.not. atom_i_in_frag .and. atom_j_in_frag)) then
               n_caps = n_caps + 1
            end if
         end do
      end if

      ! Allocate arrays with space for original atoms + caps
      fragment%n_atoms = n_atoms_no_caps + n_caps
      fragment%n_caps = n_caps
      allocate (fragment%element_numbers(fragment%n_atoms))
      allocate (fragment%coordinates(3, fragment%n_atoms))
      if (n_caps > 0) allocate (fragment%cap_replaces_atom(n_caps))

      ! Copy original atoms
      frag_atom_idx = 0
      if (use_explicit_fragments) then
         ! Variable-sized: copy atoms based on explicit fragment definitions
         do i = 1, n_monomers_in_frag
            mono_idx = monomer_indices(i)
            do j = 1, sys_geom%fragment_sizes(mono_idx)
               frag_atom_idx = frag_atom_idx + 1
               ! fragment_atoms is 0-indexed, so +1 for Fortran arrays
               atom_global_idx = sys_geom%fragment_atoms(j, mono_idx) + 1
               fragment%element_numbers(frag_atom_idx) = sys_geom%element_numbers(atom_global_idx)
               fragment%coordinates(:, frag_atom_idx) = sys_geom%coordinates(:, atom_global_idx)
            end do
         end do
      else
         ! Fixed-size: use atoms_per_monomer
         do i = 1, n_monomers_in_frag
            mono_idx = monomer_indices(i)
            atom_start = (mono_idx - 1)*atoms_per_monomer + 1
            atom_end = mono_idx*atoms_per_monomer

            fragment%element_numbers(frag_atom_idx + 1:frag_atom_idx + atoms_per_monomer) = &
               sys_geom%element_numbers(atom_start:atom_end)

            fragment%coordinates(:, frag_atom_idx + 1:frag_atom_idx + atoms_per_monomer) = &
               sys_geom%coordinates(:, atom_start:atom_end)

            frag_atom_idx = frag_atom_idx + atoms_per_monomer
         end do
      end if

      ! Add hydrogen caps at end (if any)
      if (present(bonds) .and. n_caps > 0) then
         cap_idx = 0
         do ibond = 1, size(bonds)
            if (.not. bonds(ibond)%is_broken) cycle

            atom_i_in_frag = any(atoms_in_fragment == bonds(ibond)%atom_i)
            atom_j_in_frag = any(atoms_in_fragment == bonds(ibond)%atom_j)

            if (atom_i_in_frag .and. .not. atom_j_in_frag) then
               ! atom_i is in fragment, atom_j is not → cap at position of atom_j
               cap_idx = cap_idx + 1
               fragment%element_numbers(n_atoms_no_caps + cap_idx) = 1  ! Hydrogen
               ! Place H at position of atom_j (1-indexed for coordinates array)
               fragment%coordinates(:, n_atoms_no_caps + cap_idx) = &
                  sys_geom%coordinates(:, bonds(ibond)%atom_j + 1)
               fragment%cap_replaces_atom(cap_idx) = bonds(ibond)%atom_j

            else if (atom_j_in_frag .and. .not. atom_i_in_frag) then
               ! atom_j is in fragment, atom_i is not → cap at position of atom_i
               cap_idx = cap_idx + 1
               fragment%element_numbers(n_atoms_no_caps + cap_idx) = 1  ! Hydrogen
               ! Place H at position of atom_i (1-indexed for coordinates array)
               fragment%coordinates(:, n_atoms_no_caps + cap_idx) = &
                  sys_geom%coordinates(:, bonds(ibond)%atom_i + 1)
               fragment%cap_replaces_atom(cap_idx) = bonds(ibond)%atom_i
            end if
         end do
      end if

      ! Set electronic structure properties from system geometry
      if (use_explicit_fragments .and. allocated(sys_geom%fragment_charges) .and. &
          allocated(sys_geom%fragment_multiplicities)) then
         ! Explicit fragments: sum charges and multiplicities from constituent fragments
         fragment%charge = 0
         fragment%multiplicity = 1  ! Start with singlet assumption

         do i = 1, n_monomers_in_frag
            mono_idx = monomer_indices(i)
            fragment%charge = fragment%charge + sys_geom%fragment_charges(mono_idx)
         end do

         ! For single fragment, use its specific multiplicity
         if (n_monomers_in_frag == 1) then
            fragment%multiplicity = sys_geom%fragment_multiplicities(monomer_indices(1))
         else
            ! For multi-fragment composites, multiplicity needs careful treatment
            ! For now, default to system multiplicity (this may need refinement)
            fragment%multiplicity = sys_geom%multiplicity
         end if
      else
         ! Fixed-size monomers: use system defaults
         fragment%charge = sys_geom%charge
         fragment%multiplicity = sys_geom%multiplicity
      end if
      call fragment%compute_nelec()

      deallocate (atoms_in_fragment)

   end subroutine build_fragment_from_indices

   subroutine build_fragment_from_atom_list(sys_geom, atom_indices, n_atoms, fragment, bonds)
      !! Build a fragment from explicit atom list (for GMBE intersection fragments)
      !!
      !! Similar to build_fragment_from_indices but takes atom indices directly instead of
      !! monomer indices. Used for building intersection fragments in GMBE calculations.
      !! Intersection fragments are ALWAYS NEUTRAL (charge=0, multiplicity=1).
      !!
      !! Example: atom_indices = [3, 4, 5] builds fragment from atoms 3, 4, 5 of the system
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: atom_indices(:)  !! 0-indexed atom indices
      integer, intent(in) :: n_atoms          !! Number of atoms in list
      type(physical_fragment_t), intent(out) :: fragment
      type(bond_t), intent(in), optional :: bonds(:)  !! Connectivity for capping

      integer :: i, frag_atom_idx, atom_global_idx
      integer :: ibond, n_caps, cap_idx
      logical :: atom_i_in_frag, atom_j_in_frag

      ! Count how many caps we need
      n_caps = 0
      if (present(bonds)) then
         do ibond = 1, size(bonds)
            if (.not. bonds(ibond)%is_broken) cycle

            ! Check if exactly one atom of this bond is in the fragment
            atom_i_in_frag = any(atom_indices(1:n_atoms) == bonds(ibond)%atom_i)
            atom_j_in_frag = any(atom_indices(1:n_atoms) == bonds(ibond)%atom_j)

            ! Add cap only if one atom in fragment, other not (XOR condition)
            if ((atom_i_in_frag .and. .not. atom_j_in_frag) .or. &
                (.not. atom_i_in_frag .and. atom_j_in_frag)) then
               n_caps = n_caps + 1
            end if
         end do
      end if

      ! Allocate arrays with space for original atoms + caps
      fragment%n_atoms = n_atoms + n_caps
      fragment%n_caps = n_caps
      allocate (fragment%element_numbers(fragment%n_atoms))
      allocate (fragment%coordinates(3, fragment%n_atoms))
      if (n_caps > 0) allocate (fragment%cap_replaces_atom(n_caps))

      ! Copy original atoms (atom_indices are 0-indexed, add 1 for Fortran arrays)
      do i = 1, n_atoms
         atom_global_idx = atom_indices(i) + 1  ! Convert to 1-indexed
         fragment%element_numbers(i) = sys_geom%element_numbers(atom_global_idx)
         fragment%coordinates(:, i) = sys_geom%coordinates(:, atom_global_idx)
      end do

      ! Add hydrogen caps at end (if any)
      if (present(bonds) .and. n_caps > 0) then
         cap_idx = 0
         do ibond = 1, size(bonds)
            if (.not. bonds(ibond)%is_broken) cycle

            atom_i_in_frag = any(atom_indices(1:n_atoms) == bonds(ibond)%atom_i)
            atom_j_in_frag = any(atom_indices(1:n_atoms) == bonds(ibond)%atom_j)

            if (atom_i_in_frag .and. .not. atom_j_in_frag) then
               ! atom_i is in fragment, atom_j is not → cap at position of atom_j
               cap_idx = cap_idx + 1
               fragment%element_numbers(n_atoms + cap_idx) = 1  ! Hydrogen
               fragment%coordinates(:, n_atoms + cap_idx) = &
                  sys_geom%coordinates(:, bonds(ibond)%atom_j + 1)
               fragment%cap_replaces_atom(cap_idx) = bonds(ibond)%atom_j

            else if (atom_j_in_frag .and. .not. atom_i_in_frag) then
               ! atom_j is in fragment, atom_i is not → cap at position of atom_i
               cap_idx = cap_idx + 1
               fragment%element_numbers(n_atoms + cap_idx) = 1  ! Hydrogen
               fragment%coordinates(:, n_atoms + cap_idx) = &
                  sys_geom%coordinates(:, bonds(ibond)%atom_i + 1)
               fragment%cap_replaces_atom(cap_idx) = bonds(ibond)%atom_i
            end if
         end do
      end if

      ! Intersection fragments are ALWAYS NEUTRAL
      ! Rationale: For polypeptides, intersections are backbone atoms;
      ! charged side chains are in non-overlapping regions
      fragment%charge = 0
      fragment%multiplicity = 1
      call fragment%compute_nelec()

   end subroutine build_fragment_from_atom_list

   subroutine fragment_destroy(this)
      !! Clean up allocated memory in physical_fragment_t
      class(physical_fragment_t), intent(inout) :: this
      if (allocated(this%element_numbers)) deallocate (this%element_numbers)
      if (allocated(this%coordinates)) deallocate (this%coordinates)
      if (allocated(this%cap_replaces_atom)) deallocate (this%cap_replaces_atom)
      if (allocated(this%basis)) then
         call this%basis%destroy()
         deallocate (this%basis)
      end if
      this%n_atoms = 0
      this%charge = 0
      this%multiplicity = 1
      this%nelec = 0
      this%n_caps = 0
   end subroutine fragment_destroy

   subroutine fragment_compute_nelec(this)
      !! Compute number of electrons from atomic numbers and charge
      class(physical_fragment_t), intent(inout) :: this
      integer :: nuclear_charge

      nuclear_charge = sum(this%element_numbers)
      this%nelec = nuclear_charge - this%charge
   end subroutine fragment_compute_nelec

   subroutine fragment_set_basis(this, basis)
      !! Set the basis set for this fragment
      class(physical_fragment_t), intent(inout) :: this
      type(molecular_basis_type), intent(in) :: basis

      if (allocated(this%basis)) then
         call this%basis%destroy()
         deallocate (this%basis)
      end if

      allocate (this%basis)
      this%basis = basis
   end subroutine fragment_set_basis

   subroutine system_destroy(this)
      !! Clean up allocated memory in system_geometry_t
      class(system_geometry_t), intent(inout) :: this
      if (allocated(this%element_numbers)) deallocate (this%element_numbers)
      if (allocated(this%coordinates)) deallocate (this%coordinates)
      if (allocated(this%fragment_sizes)) deallocate (this%fragment_sizes)
      if (allocated(this%fragment_atoms)) deallocate (this%fragment_atoms)
      if (allocated(this%fragment_charges)) deallocate (this%fragment_charges)
      if (allocated(this%fragment_multiplicities)) deallocate (this%fragment_multiplicities)
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

      do i = 1, fragment%n_atoms
         centroid = centroid + fragment%coordinates(:, i)
      end do

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

      do i = 1, fragment%n_atoms
         atom_mass = element_mass(fragment%element_numbers(i))
         com = com + atom_mass*fragment%coordinates(:, i)
         total_mass = total_mass + atom_mass
      end do

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

      do i = 1, frag1%n_atoms
         do j = 1, frag2%n_atoms
            current_distance = distance_between_points(frag1%coordinates(:, i), &
                                                       frag2%coordinates(:, j))

            if (current_distance < min_distance) then
               min_distance = current_distance
            end if
         end do
      end do

   end function minimal_distance_between_fragments

end module mqc_physical_fragment
