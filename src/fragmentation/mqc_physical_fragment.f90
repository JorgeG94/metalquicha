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
   public :: check_duplicate_atoms      !! Validate fragment has no overlapping atoms
   ! TODO: in theory there should be a nice way to redistribute for a general matrix of any shape, need to think about this!
   public :: redistribute_cap_gradients  !! Redistribute hydrogen cap gradients to original atoms
   public :: redistribute_cap_hessian    !! Redistribute hydrogen cap Hessian to original atoms
   public :: to_angstrom, to_bohr       !! Unit conversion utilities
   public :: calculate_monomer_distance  !! Calculate minimal distance between monomers in a fragment

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

      ! Gradient redistribution support
      integer, allocatable :: local_to_global(:)  !! Map fragment atom index to system atom index (size: n_atoms - n_caps)

      ! Fragment distance (for screening)
      real(dp) :: distance = 0.0_dp  !! Minimal atomic distance between monomers in fragment (Angstrom, 0 for monomers)

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
      if (error%has_error()) then
         call error%add_context("mqc_physical_fragment:initialize_system_geometry")
         return
      end if

      ! Read monomer template
      ! this will be changed once we have a proper input file parsing
      call read_xyz_file(monomer_file, monomer_geom, error)
      if (error%has_error()) then
         call error%add_context("mqc_physical_fragment:initialize_system_geometry")
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

   subroutine count_hydrogen_caps(atoms_in_fragment, bonds, n_caps)
      !! Count how many hydrogen caps are needed for a fragment
      !! A cap is needed when exactly one atom of a broken bond is in the fragment
      integer, intent(in) :: atoms_in_fragment(:)  !! 0-indexed atom indices in fragment
      type(bond_t), intent(in), optional :: bonds(:)
      integer, intent(out) :: n_caps

      integer :: ibond
      logical :: atom_i_in_frag, atom_j_in_frag

      n_caps = 0
      if (.not. present(bonds)) return

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

   end subroutine count_hydrogen_caps

   subroutine add_hydrogen_caps(atoms_in_fragment, bonds, sys_geom, fragment, base_atom_count)
      !! Add hydrogen caps to fragment for broken bonds
      !! Caps are placed at the position of the atom outside the fragment
      integer, intent(in) :: atoms_in_fragment(:)  !! 0-indexed atom indices in fragment
      type(bond_t), intent(in) :: bonds(:)
      type(system_geometry_t), intent(in) :: sys_geom
      type(physical_fragment_t), intent(inout) :: fragment
      integer, intent(in) :: base_atom_count  !! Number of non-cap atoms

      integer :: ibond, cap_idx
      logical :: atom_i_in_frag, atom_j_in_frag

      if (fragment%n_caps == 0) return

      cap_idx = 0
      do ibond = 1, size(bonds)
         if (.not. bonds(ibond)%is_broken) cycle

         atom_i_in_frag = any(atoms_in_fragment == bonds(ibond)%atom_i)
         atom_j_in_frag = any(atoms_in_fragment == bonds(ibond)%atom_j)

         if (atom_i_in_frag .and. .not. atom_j_in_frag) then
            ! atom_i is in fragment, atom_j is not → cap at position of atom_j
            cap_idx = cap_idx + 1
            fragment%element_numbers(base_atom_count + cap_idx) = 1  ! Hydrogen
            ! Place H at position of atom_j (1-indexed for coordinates array)
            fragment%coordinates(:, base_atom_count + cap_idx) = &
               sys_geom%coordinates(:, bonds(ibond)%atom_j + 1)
            fragment%cap_replaces_atom(cap_idx) = bonds(ibond)%atom_j

         else if (atom_j_in_frag .and. .not. atom_i_in_frag) then
            ! atom_j is in fragment, atom_i is not → cap at position of atom_i
            cap_idx = cap_idx + 1
            fragment%element_numbers(base_atom_count + cap_idx) = 1  ! Hydrogen
            ! Place H at position of atom_i (1-indexed for coordinates array)
            fragment%coordinates(:, base_atom_count + cap_idx) = &
               sys_geom%coordinates(:, bonds(ibond)%atom_i + 1)
            fragment%cap_replaces_atom(cap_idx) = bonds(ibond)%atom_i
         end if
      end do

   end subroutine add_hydrogen_caps

   subroutine build_fragment_from_indices(sys_geom, monomer_indices, fragment, error, bonds)
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
      type(error_t), intent(out) :: error
      type(bond_t), intent(in), optional :: bonds(:)  !! Connectivity information for capping

      integer :: n_monomers_in_frag, atoms_per_monomer, n_atoms_no_caps
      integer :: i, j, mono_idx, atom_start, atom_end, frag_atom_idx
      integer :: atom_i, atom_j, n_caps
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
      call count_hydrogen_caps(atoms_in_fragment, bonds, n_caps)

      ! Allocate arrays with space for original atoms + caps
      fragment%n_atoms = n_atoms_no_caps + n_caps
      fragment%n_caps = n_caps
      allocate (fragment%element_numbers(fragment%n_atoms))
      allocate (fragment%coordinates(3, fragment%n_atoms))
      if (n_caps > 0) allocate (fragment%cap_replaces_atom(n_caps))
      allocate (fragment%local_to_global(n_atoms_no_caps))  ! Only non-cap atoms

      ! Copy original atoms and build local→global mapping
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
               fragment%local_to_global(frag_atom_idx) = atom_global_idx  ! Store 1-indexed global position
            end do
         end do
      else
         ! Fixed-size: use atoms_per_monomer
         do i = 1, n_monomers_in_frag
            mono_idx = monomer_indices(i)
            atom_start = (mono_idx - 1)*atoms_per_monomer + 1
            atom_end = mono_idx*atoms_per_monomer

            ! Copy coordinates and elements
            do atom_i = atom_start, atom_end
               frag_atom_idx = frag_atom_idx + 1
               fragment%element_numbers(frag_atom_idx) = sys_geom%element_numbers(atom_i)
               fragment%coordinates(:, frag_atom_idx) = sys_geom%coordinates(:, atom_i)
               fragment%local_to_global(frag_atom_idx) = atom_i  ! Store 1-indexed global position
            end do
         end do
      end if

      ! Add hydrogen caps at end (if any)
      if (present(bonds) .and. n_caps > 0) then
         call add_hydrogen_caps(atoms_in_fragment, bonds, sys_geom, fragment, n_atoms_no_caps)
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

      ! Validate: check for spatially overlapping atoms
      call check_duplicate_atoms(fragment, error)
      if (error%has_error()) then
         call error%add_context("mqc_physical_fragment:build_fragment_from_indices")
         return
      end if

      ! Calculate minimal distance between monomers in this fragment
      fragment%distance = calculate_monomer_distance(sys_geom, monomer_indices)

      deallocate (atoms_in_fragment)

   end subroutine build_fragment_from_indices

   subroutine build_fragment_from_atom_list(sys_geom, atom_indices, n_atoms, fragment, error, bonds)
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
      type(error_t), intent(out) :: error
      type(bond_t), intent(in), optional :: bonds(:)  !! Connectivity for capping

      integer :: i, frag_atom_idx, atom_global_idx
      integer :: n_caps

      ! Count how many caps we need
      call count_hydrogen_caps(atom_indices(1:n_atoms), bonds, n_caps)

      ! Allocate arrays with space for original atoms + caps
      fragment%n_atoms = n_atoms + n_caps
      fragment%n_caps = n_caps
      allocate (fragment%element_numbers(fragment%n_atoms))
      allocate (fragment%coordinates(3, fragment%n_atoms))
      if (n_caps > 0) allocate (fragment%cap_replaces_atom(n_caps))
      allocate (fragment%local_to_global(n_atoms))  ! Only non-cap atoms

      ! Copy original atoms and build local→global mapping (atom_indices are 0-indexed, add 1 for Fortran arrays)
      do i = 1, n_atoms
         atom_global_idx = atom_indices(i) + 1  ! Convert to 1-indexed
         fragment%element_numbers(i) = sys_geom%element_numbers(atom_global_idx)
         fragment%coordinates(:, i) = sys_geom%coordinates(:, atom_global_idx)
         fragment%local_to_global(i) = atom_global_idx  ! Store 1-indexed global position
      end do

      ! Add hydrogen caps at end (if any)
      if (present(bonds) .and. n_caps > 0) then
         call add_hydrogen_caps(atom_indices(1:n_atoms), bonds, sys_geom, fragment, n_atoms)
      end if

      ! Intersection fragments are ALWAYS NEUTRAL
      ! Rationale: For polypeptides, intersections are backbone atoms;
      ! charged side chains are in non-overlapping regions
      fragment%charge = 0
      fragment%multiplicity = 1
      call fragment%compute_nelec()

      ! Validate: check for spatially overlapping atoms
      call check_duplicate_atoms(fragment, error)
      if (error%has_error()) then
         call error%add_context("mqc_physical_fragment:build_fragment_from_atom_list")
         return
      end if

   end subroutine build_fragment_from_atom_list

   subroutine redistribute_cap_gradients(fragment, fragment_gradient, system_gradient)
      !! Redistribute hydrogen cap gradients to original atoms
      !!
      !! This subroutine handles gradient redistribution for fragments with hydrogen caps.
      !! Hydrogen caps are virtual atoms added at broken bonds - their gradients represent
      !! forces at the bond breakpoint and must be transferred to the original atoms they replace.
      !!
      !! Algorithm:
      !!   1. For real atoms (indices 1 to n_atoms - n_caps):
      !!      Accumulate gradient to system using local_to_global mapping
      !!   2. For hydrogen caps (indices n_atoms - n_caps + 1 to n_atoms):
      !!      Add cap gradient to the original atom it replaces (from cap_replaces_atom)
      !!
      !! Example:
      !!   Fragment: [C, C, H_cap] where H_cap replaces atom 5 in system
      !!   Fragment gradient: [(3,1), (3,2), (3,3)]
      !!   - Atoms 1,2: accumulate to system using local_to_global
      !!   - Atom 3 (cap): add gradient to system atom 5 (cap_replaces_atom(1) + 1)
      type(physical_fragment_t), intent(in) :: fragment
      real(dp), intent(in) :: fragment_gradient(:, :)   !! (3, n_atoms_fragment)
      real(dp), intent(inout) :: system_gradient(:, :)  !! (3, n_atoms_system)

      integer :: i, local_idx, global_idx
      integer :: i_cap, local_cap_idx, global_original_idx
      integer :: n_real_atoms

      n_real_atoms = fragment%n_atoms - fragment%n_caps

      ! Accumulate gradients for real atoms using local→global mapping
      do i = 1, n_real_atoms
         global_idx = fragment%local_to_global(i)
         system_gradient(:, global_idx) = system_gradient(:, global_idx) + fragment_gradient(:, i)
      end do

      ! Redistribute cap gradients to original atoms they replace
      if (fragment%n_caps > 0) then
         do i_cap = 1, fragment%n_caps
            local_cap_idx = n_real_atoms + i_cap
            ! cap_replaces_atom is 0-indexed, add 1 for Fortran arrays
            global_original_idx = fragment%cap_replaces_atom(i_cap) + 1

            ! Add cap gradient to the atom it replaces
            system_gradient(:, global_original_idx) = system_gradient(:, global_original_idx) + &
                                                      fragment_gradient(:, local_cap_idx)
         end do
      end if

   end subroutine redistribute_cap_gradients

   subroutine redistribute_cap_hessian(fragment, fragment_hessian, system_hessian)
      !! Redistribute hydrogen cap Hessian to original atoms
      !!
      !! This subroutine handles Hessian redistribution for fragments with hydrogen caps.
      !! The Hessian is a rank-2 tensor (3N × 3N) representing second derivatives of energy
      !! with respect to atomic coordinates. Similar to gradient redistribution, cap contributions
      !! must be transferred to the original atoms they replace.
      !!
      !! Algorithm:
      !!   1. For real atoms (indices 1 to n_atoms - n_caps):
      !!      Accumulate Hessian blocks to system using local_to_global mapping for both dimensions
      !!   2. For hydrogen caps (indices n_atoms - n_caps + 1 to n_atoms):
      !!      Add cap Hessian blocks (row and column) to the original atom it replaces
      !!
      !! Note: Hessian is stored as a flattened 2D array (3*n_atoms, 3*n_atoms)
      !!       where rows and columns are grouped by atoms (x,y,z for atom 1, then x,y,z for atom 2, etc.)
      type(physical_fragment_t), intent(in) :: fragment
      real(dp), intent(in) :: fragment_hessian(:, :)   !! (3*n_atoms_fragment, 3*n_atoms_fragment)
      real(dp), intent(inout) :: system_hessian(:, :)  !! (3*n_atoms_system, 3*n_atoms_system)

      integer :: i, j, local_i, local_j, global_i, global_j
      integer :: icart, jcart
      integer :: i_cap, local_cap_idx, global_original_idx
      integer :: n_real_atoms
      integer :: i_cap_2, local_cap_idx_2, global_original_idx_2

      n_real_atoms = fragment%n_atoms - fragment%n_caps

      ! Accumulate Hessian blocks for real atoms using local→global mapping
      ! Both row (i) and column (j) dimensions need mapping
      do i = 1, n_real_atoms
         global_i = fragment%local_to_global(i)
         do j = 1, n_real_atoms
            global_j = fragment%local_to_global(j)

            ! Copy 3×3 block for atom pair (i,j)
            do icart = 0, 2  ! x, y, z for atom i
               do jcart = 0, 2  ! x, y, z for atom j
                  system_hessian(3*(global_i - 1) + icart + 1, 3*(global_j - 1) + jcart + 1) = &
                     system_hessian(3*(global_i - 1) + icart + 1, 3*(global_j - 1) + jcart + 1) + &
                     fragment_hessian(3*(i - 1) + icart + 1, 3*(j - 1) + jcart + 1)
               end do
            end do
         end do
      end do

      ! Redistribute cap Hessian blocks to original atoms they replace
      if (fragment%n_caps > 0) then
         do i_cap = 1, fragment%n_caps
            local_cap_idx = n_real_atoms + i_cap
            global_original_idx = fragment%cap_replaces_atom(i_cap) + 1

            ! Cap rows: redistribute to original atom (cap derivatives w.r.t. all other atoms)
            do j = 1, n_real_atoms
               global_j = fragment%local_to_global(j)
               do icart = 0, 2
                  do jcart = 0, 2
                     system_hessian(3*(global_original_idx - 1) + icart + 1, 3*(global_j - 1) + jcart + 1) = &
                        system_hessian(3*(global_original_idx - 1) + icart + 1, 3*(global_j - 1) + jcart + 1) + &
                        fragment_hessian(3*(local_cap_idx - 1) + icart + 1, 3*(j - 1) + jcart + 1)
                  end do
               end do
            end do

            ! Cap columns: redistribute to original atom (all other atoms' derivatives w.r.t. cap)
            do i = 1, n_real_atoms
               global_i = fragment%local_to_global(i)
               do icart = 0, 2
                  do jcart = 0, 2
                     system_hessian(3*(global_i - 1) + icart + 1, 3*(global_original_idx - 1) + jcart + 1) = &
                        system_hessian(3*(global_i - 1) + icart + 1, 3*(global_original_idx - 1) + jcart + 1) + &
                        fragment_hessian(3*(i - 1) + icart + 1, 3*(local_cap_idx - 1) + jcart + 1)
                  end do
               end do
            end do

            ! Cap-cap blocks: redistribute to original atom diagonal block
            do i_cap_2 = 1, fragment%n_caps
               local_cap_idx_2 = n_real_atoms + i_cap_2
               global_original_idx_2 = fragment%cap_replaces_atom(i_cap_2) + 1

               do icart = 0, 2
                  do jcart = 0, 2
                    system_hessian(3*(global_original_idx - 1) + icart + 1, 3*(global_original_idx_2 - 1) + jcart + 1) = &
                    system_hessian(3*(global_original_idx - 1) + icart + 1, 3*(global_original_idx_2 - 1) + jcart + 1) + &
                        fragment_hessian(3*(local_cap_idx - 1) + icart + 1, 3*(local_cap_idx_2 - 1) + jcart + 1)
                  end do
               end do
            end do
         end do
      end if

   end subroutine redistribute_cap_hessian

   subroutine check_duplicate_atoms(fragment, error)
      !! Validate that fragment has no spatially overlapping atoms
      !! Checks if any two atoms are too close together (< 0.01 Bohr)
      !! This catches bugs in geometry construction or fragment building
      use pic_logger, only: logger => global_logger
      use pic_io, only: to_char

      type(physical_fragment_t), intent(in) :: fragment
      type(error_t), intent(out) :: error

      integer :: i, j, n_atoms
      real(dp) :: distance, dx, dy, dz
      real(dp), parameter :: MIN_ATOM_DISTANCE = 0.01_dp  !! Bohr - atoms closer than this are overlapping

      ! Only check non-cap atoms (caps can be close to replaced atoms)
      n_atoms = fragment%n_atoms - fragment%n_caps

      if (n_atoms < 2) return

      do i = 1, n_atoms - 1
         do j = i + 1, n_atoms
            dx = fragment%coordinates(1, i) - fragment%coordinates(1, j)
            dy = fragment%coordinates(2, i) - fragment%coordinates(2, j)
            dz = fragment%coordinates(3, i) - fragment%coordinates(3, j)
            distance = sqrt(dx*dx + dy*dy + dz*dz)

            if (distance < MIN_ATOM_DISTANCE) then
               ! Build detailed error message
               call error%set(ERROR_VALIDATION, &
                              "Fragment contains overlapping atoms "//to_char(i)//" and "//to_char(j)// &
                              " (distance: "//to_char(distance)//" Bohr). "// &
                              "This indicates bad input geometry or a bug in fragment construction.")

               ! Log detailed information for debugging
               call logger%error("ERROR: Fragment contains overlapping atoms!")
               call logger%error("  Atoms "//to_char(i)//" and "//to_char(j)//" are too close together")
               call logger%error("  Distance: "//to_char(distance)//" Bohr ("// &
                                 to_char(distance*0.529177_dp)//" Angstrom)")
               call logger%error("  Atom "//to_char(i)//": "// &
                                 element_number_to_symbol(fragment%element_numbers(i))// &
                                 " at ("//to_char(fragment%coordinates(1, i))//", "// &
                                 to_char(fragment%coordinates(2, i))//", "// &
                                 to_char(fragment%coordinates(3, i))//") Bohr")
               call logger%error("  Atom "//to_char(j)//": "// &
                                 element_number_to_symbol(fragment%element_numbers(j))// &
                                 " at ("//to_char(fragment%coordinates(1, j))//", "// &
                                 to_char(fragment%coordinates(2, j))//", "// &
                                 to_char(fragment%coordinates(3, j))//") Bohr")
               return
            end if
         end do
      end do
   end subroutine check_duplicate_atoms

   subroutine fragment_destroy(this)
      !! Clean up allocated memory in physical_fragment_t
      class(physical_fragment_t), intent(inout) :: this
      if (allocated(this%element_numbers)) deallocate (this%element_numbers)
      if (allocated(this%coordinates)) deallocate (this%coordinates)
      if (allocated(this%cap_replaces_atom)) deallocate (this%cap_replaces_atom)
      if (allocated(this%local_to_global)) deallocate (this%local_to_global)
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

   pure function calculate_monomer_distance(sys_geom, monomer_indices) result(min_distance)
      !! Calculate minimal atomic distance between monomers in a fragment
      !! For single monomer (size 1), returns 0.0
      !! For multi-monomer fragments, returns minimal distance between atoms in different monomers
      !! Result is in Angstrom
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: monomer_indices(:)
      real(dp) :: min_distance

      integer :: n_monomers, i, j, iatom, jatom
      integer :: mon_i, mon_j
      integer :: atom_start_i, atom_end_i, atom_start_j, atom_end_j
      real(dp) :: dist, dx, dy, dz
      logical :: is_variable_size

      n_monomers = size(monomer_indices)

      ! Monomers have distance 0
      if (n_monomers == 1) then
         min_distance = 0.0_dp
         return
      end if

      ! Check if we have variable-sized fragments
      is_variable_size = allocated(sys_geom%fragment_sizes)

      ! Initialize with huge value
      min_distance = huge(1.0_dp)

      ! Loop over all pairs of monomers
      do i = 1, n_monomers - 1
         mon_i = monomer_indices(i)
         do j = i + 1, n_monomers
            mon_j = monomer_indices(j)

            if (is_variable_size) then
               ! Variable-sized fragments
               do iatom = 1, sys_geom%fragment_sizes(mon_i)
                  atom_start_i = sys_geom%fragment_atoms(iatom, mon_i) + 1  ! Convert to 1-indexed
                  do jatom = 1, sys_geom%fragment_sizes(mon_j)
                     atom_start_j = sys_geom%fragment_atoms(jatom, mon_j) + 1  ! Convert to 1-indexed

                     ! Calculate distance (coordinates in Bohr)
                     dx = sys_geom%coordinates(1, atom_start_i) - sys_geom%coordinates(1, atom_start_j)
                     dy = sys_geom%coordinates(2, atom_start_i) - sys_geom%coordinates(2, atom_start_j)
                     dz = sys_geom%coordinates(3, atom_start_i) - sys_geom%coordinates(3, atom_start_j)
                     dist = sqrt(dx*dx + dy*dy + dz*dz)

                     if (dist < min_distance) min_distance = dist
                  end do
               end do
            else
               ! Fixed-sized monomers
               atom_start_i = (mon_i - 1)*sys_geom%atoms_per_monomer + 1
               atom_end_i = mon_i*sys_geom%atoms_per_monomer

               atom_start_j = (mon_j - 1)*sys_geom%atoms_per_monomer + 1
               atom_end_j = mon_j*sys_geom%atoms_per_monomer

               ! Loop over all atom pairs
               do iatom = atom_start_i, atom_end_i
                  do jatom = atom_start_j, atom_end_j
                     ! Calculate distance (coordinates in Bohr)
                     dx = sys_geom%coordinates(1, iatom) - sys_geom%coordinates(1, jatom)
                     dy = sys_geom%coordinates(2, iatom) - sys_geom%coordinates(2, jatom)
                     dz = sys_geom%coordinates(3, iatom) - sys_geom%coordinates(3, jatom)
                     dist = sqrt(dx*dx + dy*dy + dz*dz)

                     if (dist < min_distance) min_distance = dist
                  end do
               end do
            end if
         end do
      end do

      ! Convert from Bohr to Angstrom
      min_distance = to_angstrom(min_distance)

   end function calculate_monomer_distance

end module mqc_physical_fragment
