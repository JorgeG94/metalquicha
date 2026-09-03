!! Physical fragment and whole-system geometry types
module mqc_physical_fragment
   !! Fragments as the QC methods see them, the system they are cut from, and
   !! the redistribution of cap derivatives back onto the atoms a cap replaces.
   use pic_types, only: dp, default_int
   use mqc_geometry, only: geometry_type
   use mqc_xyz_reader, only: read_xyz_file
   use mqc_elements, only: element_symbol_to_number, element_number_to_symbol, element_mass
   use mqc_cgto, only: molecular_basis_type
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_physical_constants, only: BOHR_TO_ANGSTROM, ANGSTROM_TO_BOHR
   use mqc_program_limits, only: MIN_ATOM_DISTANCE
   implicit none
   private

   public :: bond_t                      !! Bond connectivity type
   public :: physical_fragment_t         !! Single molecular fragment type
   public :: system_geometry_t          !! Complete system geometry type
   public :: initialize_system_geometry  !! System geometry initialization
   public :: build_fragment_from_indices
   public :: fragment_charge_multiplicity  !! Charge/multiplicity a set of monomers forms
   public :: build_fragment_from_atom_list  !! Build fragment from explicit atom indices (for intersections)
   public :: check_duplicate_atoms      !! Validate fragment has no overlapping atoms
   public :: check_system_geometry      !! Validate the whole system before any work starts
   ! TODO(mqc): one redistribution for an array of any rank, rather than the
   ! three near-identical routines below.
   public :: redistribute_cap_gradients  !! Redistribute hydrogen cap gradients to original atoms
   public :: redistribute_cap_hessian    !! Redistribute hydrogen cap Hessian to original atoms
   public :: redistribute_cap_dipole_derivatives  !! Redistribute hydrogen cap dipole derivatives to original atoms
   public :: to_angstrom, to_bohr       !! Unit conversion utilities
   public :: calculate_monomer_distance  !! Calculate minimal distance between monomers in a fragment

   type :: bond_t
      !! A bond between two system atoms, and whether a partition cut it
      integer :: atom_i = 0        !! First atom index (0-indexed)
      integer :: atom_j = 0        !! Second atom index (0-indexed)
      integer :: order = 1         !! Bond order (1=single, 2=double, 3=triple)
      logical :: is_broken = .false.  !! Whether this bond crosses fragment boundaries
   end type bond_t

   type :: physical_fragment_t
      !! One fragment as a QC method sees it: geometry, electron count, basis
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
      integer, allocatable :: cap_bonded_to(:)
         !! The in-fragment atom at the other end of each cut bond (size: n_caps),
         !! 0-based like `cap_replaces_atom`. Needed only to differentiate: the
         !! cap sits on the line between the two, so both ends move it.
      real(dp) :: cap_scale = 1.0_dp
         !! Where each cap sits along its cut bond,
         !! `R_H = R_kept + s (R_gone - R_kept)`. One value for the fragment.

      ! Gradient redistribution support
      integer, allocatable :: local_to_global(:)  !! Map fragment atom index to system atom index (size: n_atoms - n_caps)
      logical, allocatable :: is_ghost(:)
         !! Which atoms carry basis functions but no nucleus and no electrons
         !! (size: n_atoms). Unallocated means none, which is every fragment a
         !! plain expansion builds. `charge`, `nelec` and `multiplicity` count
         !! the real atoms only.

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
      !! The whole system, laid out by monomer, that fragments are cut from
      integer :: n_monomers        !! Number of monomer units in system
      integer :: atoms_per_monomer  !! Atoms in each monomer (0 if variable-sized)
      integer :: total_atoms       !! Total number of atoms
      integer, allocatable :: element_numbers(:)  !! Atomic numbers for all atoms
      real(dp), allocatable :: coordinates(:, :)  !! All coordinates (3, total_atoms) in Bohr

      ! Electronic structure properties
      integer :: charge = 0        !! Net molecular charge (electrons)
      integer :: multiplicity = 1  !! Spin multiplicity (2S+1)

      ! For variable-sized fragments (explicit fragment definitions)
      integer, allocatable :: fragment_sizes(:)      !! Number of atoms in each fragment (n_monomers)
      integer, allocatable :: fragment_atoms(:, :)   !! Atom indices for each fragment (max_frag_size, n_monomers), 0-indexed
      integer, allocatable :: fragment_charges(:)    !! Charge for each fragment (n_monomers)
      integer, allocatable :: fragment_multiplicities(:)  !! Multiplicity for each fragment (n_monomers)

      ! Connectivity information for hydrogen capping
      type(bond_t), allocatable :: bonds(:)  !! Bond connectivity (for H-capping broken bonds)
      real(dp) :: cap_scale = 1.0_dp
         !! Where a cap sits along the bond it closes; copied to every fragment
         !! built from this geometry, for the same reason `bonds` is carried here.

      ! 256 characters a path, matching `driver_config_t%fragment_potentials`,
      ! which is where these are copied to.
      character(len=256), allocatable :: fragment_potentials(:)
         !! One effective fragment potential file per fragment, in fragment
         !! order, or unallocated when this is an ordinary quantum system.
         !! For the callers with no deck -- the C and Python interfaces send no
         !! `molecules` block for a potential to travel in.
   contains
      procedure :: destroy => system_destroy  !! Memory cleanup
   end type system_geometry_t

contains

   pure elemental function to_angstrom(bohr_value) result(angstrom_value)
      !! Convert coordinate from Bohr to Angstrom
      real(dp), intent(in) :: bohr_value
      real(dp) :: angstrom_value
      angstrom_value = bohr_value*BOHR_TO_ANGSTROM
   end function to_angstrom

   pure elemental function to_bohr(angstrom_value) result(bohr_value)
      !! Convert coordinate from Angstrom to Bohr
      real(dp), intent(in) :: angstrom_value
      real(dp) :: bohr_value
      bohr_value = angstrom_value*ANGSTROM_TO_BOHR
   end function to_bohr

   subroutine initialize_system_geometry(full_geom_file, monomer_file, sys_geom, error)
      !! Read full geometry and monomer template into a `system_geometry_t`
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

      ! TODO(mqc): predates the JSON deck; a monomer template read from a
      ! second xyz is not how any current input describes a partition.
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

      ! TODO(mqc): this can be a sys_geom%allocate()
      allocate (sys_geom%element_numbers(sys_geom%total_atoms))
      allocate (sys_geom%coordinates(3, sys_geom%total_atoms))

      do i = 1, sys_geom%total_atoms
         sys_geom%element_numbers(i) = element_symbol_to_number(full_geom%elements(i))
      end do

      ! Store coordinates in Bohr (convert from Angstroms)
      ! TODO(mqc): need a way to handle units
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
      integer :: kept, gone
      real(dp) :: s

      if (fragment%n_caps == 0) return

      ! `R_H = R_kept + s (R_gone - R_kept)`. At s = 1 the cap lands exactly on
      ! the removed atom, which is what this did before the scale existed and is
      ! still the default; the redistribution below then puts the whole cap
      ! gradient on that atom, as it always has.
      s = fragment%cap_scale

      cap_idx = 0
      do ibond = 1, size(bonds)
         if (.not. bonds(ibond)%is_broken) cycle

         atom_i_in_frag = any(atoms_in_fragment == bonds(ibond)%atom_i)
         atom_j_in_frag = any(atoms_in_fragment == bonds(ibond)%atom_j)

         if (atom_i_in_frag .neqv. atom_j_in_frag) then
            if (atom_i_in_frag) then
               kept = bonds(ibond)%atom_i
               gone = bonds(ibond)%atom_j
            else
               kept = bonds(ibond)%atom_j
               gone = bonds(ibond)%atom_i
            end if

            cap_idx = cap_idx + 1
            fragment%element_numbers(base_atom_count + cap_idx) = 1  ! Hydrogen
            ! 1-indexed for the coordinates array; the stored indices stay 0-based
            fragment%coordinates(:, base_atom_count + cap_idx) = &
               sys_geom%coordinates(:, kept + 1) + &
               s*(sys_geom%coordinates(:, gone + 1) - sys_geom%coordinates(:, kept + 1))
            fragment%cap_replaces_atom(cap_idx) = gone
            fragment%cap_bonded_to(cap_idx) = kept
         end if
      end do

   end subroutine add_hydrogen_caps

   subroutine build_fragment_core(sys_geom, monomer_indices, fragment, error, bonds)
      !! Build a fragment from monomer indices, capping any bond it cut
      !!
      !! Caps go at the end of the atom list. Handles both fixed-size monomers
      !! and an explicit variable-sized partition.
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
      fragment%cap_scale = sys_geom%cap_scale
      allocate (fragment%element_numbers(fragment%n_atoms))
      allocate (fragment%coordinates(3, fragment%n_atoms))
      if (n_caps > 0) allocate (fragment%cap_replaces_atom(n_caps))
      if (n_caps > 0) allocate (fragment%cap_bonded_to(n_caps))
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

      ! Charge and multiplicity of the composite, by the same rule the
      ! per-fragment breakdown reports -- computed in one place so the number
      ! the SCF runs and the number the table prints cannot drift apart.
      call fragment_charge_multiplicity(sys_geom, monomer_indices, &
                                        fragment%charge, fragment%multiplicity)
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

   end subroutine build_fragment_core

   subroutine build_fragment_from_indices(sys_geom, signed_indices, fragment, error, bonds)
      !! Build a fragment from monomer indices, ghosting any that are negative
      !!
      !! A positive entry is a real monomer, a negative one contributes ghost
      !! centres; the geometry is the union either way, which is what a
      !! counterpoise-corrected term needs. All-positive input is exactly
      !! `build_fragment_core`.
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: signed_indices(:)
      type(physical_fragment_t), intent(out) :: fragment
      type(error_t), intent(out) :: error
      type(bond_t), intent(in), optional :: bonds(:)

      integer, allocatable :: union_indices(:), real_indices(:)
      integer :: i, j, n_real, first, n_here

      allocate (union_indices(size(signed_indices)))
      union_indices = abs(signed_indices)

      call build_fragment_core(sys_geom, union_indices, fragment, error, bonds)
      if (error%has_error()) return

      n_real = count(signed_indices > 0)
      if (n_real == size(signed_indices)) return   ! nothing ghosted

      ! Mark the ghosted monomers' atoms. `build_fragment_core` lays the atoms
      ! out monomer by monomer in the order given and appends caps at the end,
      ! so walking the same order finds them. Caps stay real.
      allocate (fragment%is_ghost(fragment%n_atoms))
      fragment%is_ghost = .false.
      first = 1
      do i = 1, size(signed_indices)
         if (allocated(sys_geom%fragment_atoms)) then
            n_here = sys_geom%fragment_sizes(union_indices(i))
         else
            n_here = sys_geom%atoms_per_monomer
         end if
         if (signed_indices(i) < 0) then
            do j = first, first + n_here - 1
               fragment%is_ghost(j) = .true.
            end do
         end if
         first = first + n_here
      end do

      ! Charge and multiplicity belong to the real monomers alone -- a ghost has
      ! no charge to contribute and no electrons to pair. Recomputed rather than
      ! adjusted, so the rule stays in one place.
      allocate (real_indices(n_real))
      j = 0
      do i = 1, size(signed_indices)
         if (signed_indices(i) > 0) then
            j = j + 1
            real_indices(j) = signed_indices(i)
         end if
      end do
      call fragment_charge_multiplicity(sys_geom, real_indices, &
                                        fragment%charge, fragment%multiplicity)
      call fragment%compute_nelec()
   end subroutine build_fragment_from_indices

   pure subroutine fragment_charge_multiplicity(sys_geom, monomer_indices, charge, multiplicity)
      !! Total charge and spin multiplicity of the fragment these monomers form
      !!
      !! With an explicit partition, monomer charges sum and unpaired-electron
      !! counts sum with them -- high-spin coupling, which a low-spin state
      !! overrides by giving the fragment its own multiplicity. With fixed-size
      !! monomers it is the whole-system charge and multiplicity.
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: monomer_indices(:)
      integer, intent(out) :: charge, multiplicity

      integer :: i, mono_idx

      if (allocated(sys_geom%fragment_atoms) .and. allocated(sys_geom%fragment_charges) &
          .and. allocated(sys_geom%fragment_multiplicities)) then
         charge = 0
         multiplicity = 1
         do i = 1, size(monomer_indices)
            mono_idx = monomer_indices(i)
            charge = charge + sys_geom%fragment_charges(mono_idx)
            multiplicity = multiplicity + (sys_geom%fragment_multiplicities(mono_idx) - 1)
         end do
      else
         charge = sys_geom%charge
         multiplicity = sys_geom%multiplicity
      end if
   end subroutine fragment_charge_multiplicity

   subroutine build_fragment_from_atom_list(sys_geom, atom_indices, n_atoms, fragment, error, bonds)
      !! Build a fragment from explicit atom indices, for GMBE intersections
      !!
      !! Intersection fragments are always neutral singlets.
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
      fragment%cap_scale = sys_geom%cap_scale
      allocate (fragment%element_numbers(fragment%n_atoms))
      allocate (fragment%coordinates(3, fragment%n_atoms))
      if (n_caps > 0) allocate (fragment%cap_replaces_atom(n_caps))
      if (n_caps > 0) allocate (fragment%cap_bonded_to(n_caps))
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

      ! Always neutral: in a polypeptide the intersections are backbone atoms,
      ! and charged side chains fall in the non-overlapping regions.
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

   pure subroutine cap_targets(fragment, i_cap, idx, wgt)
      !! Which atoms a cap's derivative belongs to, and in what proportion
      !!
      !! The cap is a function of two real atoms and nothing else,
      !!
      !!     R_H = R_kept + s (R_gone - R_kept)
      !!
      !! so `dR_H/dR_gone = s` and `dR_H/dR_kept = 1 - s` -- the chain rule
      !! exactly, and the same two weights serve the gradient and the Hessian.
      !! At `s = 1` all the weight is on the removed atom.
      type(physical_fragment_t), intent(in) :: fragment
      integer, intent(in) :: i_cap
      integer, intent(out) :: idx(2)   !! 1-based system atom indices
      real(dp), intent(out) :: wgt(2)

      idx(1) = fragment%cap_replaces_atom(i_cap) + 1
      wgt(1) = fragment%cap_scale
      if (allocated(fragment%cap_bonded_to)) then
         idx(2) = fragment%cap_bonded_to(i_cap) + 1
      else
         idx(2) = idx(1)      ! a fragment built before this existed: s is 1 and
      end if                  ! the second weight is zero, so the index is unused
      wgt(2) = 1.0_dp - fragment%cap_scale
   end subroutine cap_targets

   subroutine redistribute_cap_gradients(fragment, fragment_gradient, system_gradient, scale)
      !! Accumulate a fragment gradient into the system gradient
      !!
      !! Real atoms go through `local_to_global`; each cap is split between the
      !! two atoms that move it, by the weights `cap_targets` returns.
      type(physical_fragment_t), intent(in) :: fragment
      real(dp), intent(in) :: fragment_gradient(:, :)   !! (3, n_atoms_fragment)
      real(dp), intent(inout) :: system_gradient(:, :)  !! (3, n_atoms_system)
      real(dp), intent(in), optional :: scale  !! Weight applied to this fragment (default 1)

      integer :: i, global_idx
      integer :: i_cap, local_cap_idx
      integer :: n_real_atoms
      real(dp) :: w
      integer :: tgt(2)
      integer :: t
      real(dp) :: tw(2)

      w = 1.0_dp
      if (present(scale)) w = scale

      n_real_atoms = fragment%n_atoms - fragment%n_caps

      ! Accumulate gradients for real atoms using local→global mapping
      do i = 1, n_real_atoms
         global_idx = fragment%local_to_global(i)
         system_gradient(:, global_idx) = system_gradient(:, global_idx) + w*fragment_gradient(:, i)
      end do

      ! Redistribute cap gradients to original atoms they replace
      if (fragment%n_caps > 0) then
         do i_cap = 1, fragment%n_caps
            local_cap_idx = n_real_atoms + i_cap
            call cap_targets(fragment, i_cap, tgt, tw)

            ! Split the cap's gradient between the two atoms that move it
            do t = 1, 2
               if (tw(t) == 0.0_dp) cycle
               system_gradient(:, tgt(t)) = system_gradient(:, tgt(t)) + &
                                            w*tw(t)*fragment_gradient(:, local_cap_idx)
            end do
         end do
      end if

   end subroutine redistribute_cap_gradients

   subroutine redistribute_cap_hessian(fragment, fragment_hessian, system_hessian, scale)
      !! Accumulate a fragment Hessian into the system Hessian
      !!
      !! As `redistribute_cap_gradients`, but the chain rule applies to both
      !! indices, so a cap-cap block spreads over four atom pairs. Rows and
      !! columns are grouped by atom: x, y, z for atom 1, then atom 2, and so on.
      type(physical_fragment_t), intent(in) :: fragment
      real(dp), intent(in) :: fragment_hessian(:, :)   !! (3*n_atoms_fragment, 3*n_atoms_fragment)
      real(dp), intent(inout) :: system_hessian(:, :)  !! (3*n_atoms_system, 3*n_atoms_system)
      real(dp), intent(in), optional :: scale  !! Weight applied to this fragment (default 1)

      integer :: i, j, global_i, global_j
      integer :: icart, jcart
      integer :: i_cap, local_cap_idx, global_original_idx
      integer :: n_real_atoms
      integer :: i_cap_2, local_cap_idx_2
      real(dp) :: w
      integer :: tgt(2)
      integer :: tgt2(2)
      integer :: t, t2
      real(dp) :: tw(2)
      real(dp) :: tw2(2)

      w = 1.0_dp
      if (present(scale)) w = scale

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
                     w*fragment_hessian(3*(i - 1) + icart + 1, 3*(j - 1) + jcart + 1)
               end do
            end do
         end do
      end do

      ! Redistribute cap Hessian blocks to original atoms they replace
      if (fragment%n_caps > 0) then
         do i_cap = 1, fragment%n_caps
            local_cap_idx = n_real_atoms + i_cap
            call cap_targets(fragment, i_cap, tgt, tw)

            ! Cap rows: one row of second derivatives, shared by the two atoms
            ! that move the cap, in the same proportion the gradient used
            do t = 1, 2
               if (tw(t) == 0.0_dp) cycle
               global_original_idx = tgt(t)
               do j = 1, n_real_atoms
                  global_j = fragment%local_to_global(j)
                  do icart = 0, 2
                     do jcart = 0, 2
                        system_hessian(3*(global_original_idx - 1) + icart + 1, 3*(global_j - 1) + jcart + 1) = &
                           system_hessian(3*(global_original_idx - 1) + icart + 1, 3*(global_j - 1) + jcart + 1) + &
                           w*tw(t)*fragment_hessian(3*(local_cap_idx - 1) + icart + 1, 3*(j - 1) + jcart + 1)
                     end do
                  end do
               end do
            end do

            ! Cap columns: the mirror of the rows, same two weights
            do t = 1, 2
               if (tw(t) == 0.0_dp) cycle
               do i = 1, n_real_atoms
                  global_i = fragment%local_to_global(i)
                  do icart = 0, 2
                     do jcart = 0, 2
                        system_hessian(3*(global_i - 1) + icart + 1, 3*(tgt(t) - 1) + jcart + 1) = &
                           system_hessian(3*(global_i - 1) + icart + 1, 3*(tgt(t) - 1) + jcart + 1) + &
                           w*tw(t)*fragment_hessian(3*(i - 1) + icart + 1, 3*(local_cap_idx - 1) + jcart + 1)
                     end do
                  end do
               end do
            end do

            ! Cap-cap blocks: both indices are caps, so the chain rule applies
            ! twice and the block spreads over four atom pairs with the outer
            ! product of the weights. At s = 1 three of the four weigh nothing
            ! and this is the single block it always was.
            do i_cap_2 = 1, fragment%n_caps
               local_cap_idx_2 = n_real_atoms + i_cap_2
               call cap_targets(fragment, i_cap_2, tgt2, tw2)

               do t = 1, 2
                  if (tw(t) == 0.0_dp) cycle
                  do t2 = 1, 2
                     if (tw2(t2) == 0.0_dp) cycle
                     do icart = 0, 2
                        do jcart = 0, 2
                           system_hessian(3*(tgt(t) - 1) + icart + 1, 3*(tgt2(t2) - 1) + jcart + 1) = &
                              system_hessian(3*(tgt(t) - 1) + icart + 1, 3*(tgt2(t2) - 1) + jcart + 1) + &
                              w*tw(t)*tw2(t2)* &
                              fragment_hessian(3*(local_cap_idx - 1) + icart + 1, 3*(local_cap_idx_2 - 1) + jcart + 1)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end if

   end subroutine redistribute_cap_hessian

   subroutine redistribute_cap_dipole_derivatives(fragment, fragment_dipole_derivs, system_dipole_derivs, scale)
      !! Accumulate fragment dipole derivatives into the system's
      !!
      !! Shape (3, 3*n_atoms): one column per Cartesian coordinate, mapped like
      !! the column dimension of the Hessian.
      type(physical_fragment_t), intent(in) :: fragment
      real(dp), intent(in) :: fragment_dipole_derivs(:, :)   !! (3, 3*n_atoms_fragment)
      real(dp), intent(inout) :: system_dipole_derivs(:, :)  !! (3, 3*n_atoms_system)
      real(dp), intent(in), optional :: scale  !! Weight applied to this fragment (default 1)

      integer :: tgt(2)
      integer :: i, local_i, global_i, icart
      integer :: i_cap, local_cap_idx, t
      integer :: n_real_atoms
      real(dp) :: tw(2)
      real(dp) :: w

      w = 1.0_dp
      if (present(scale)) w = scale

      n_real_atoms = fragment%n_atoms - fragment%n_caps

      ! Accumulate dipole derivative columns for real atoms
      do i = 1, n_real_atoms
         local_i = i
         global_i = fragment%local_to_global(local_i)
         if (global_i <= 0) cycle

         do icart = 1, 3
            system_dipole_derivs(:, (global_i - 1)*3 + icart) = &
               system_dipole_derivs(:, (global_i - 1)*3 + icart) + &
               w*fragment_dipole_derivs(:, (local_i - 1)*3 + icart)
         end do
      end do

      ! Redistribute cap contributions to their original atoms
      ! Split the cap between the two atoms that move it, through the same
      ! `cap_targets` the gradient and the Hessian use.
      !
      ! This used to put the whole contribution on `cap_replaces_atom` with
      ! weight one, which is right only at `cap_scale = 1` -- the default, and
      ! why it went unnoticed. At any other scale the intensities were
      ! redistributed by different weights than the forces they belong to:
      ! silent, and wrong in exactly the quantity a capped fragment's infrared
      ! spectrum is read from. The chain rule is `dR_H/dR_gone = s` and
      ! `dR_H/dR_kept = 1 - s`, and a dipole derivative obeys it for the same
      ! reason a gradient does.
      if (fragment%n_caps > 0 .and. allocated(fragment%cap_replaces_atom)) then
         do i_cap = 1, fragment%n_caps
            local_cap_idx = n_real_atoms + i_cap
            call cap_targets(fragment, i_cap, tgt, tw)

            do t = 1, 2
               if (tw(t) == 0.0_dp) cycle
               do icart = 1, 3
                  system_dipole_derivs(:, (tgt(t) - 1)*3 + icart) = &
                     system_dipole_derivs(:, (tgt(t) - 1)*3 + icart) + &
                     w*tw(t)*fragment_dipole_derivs(:, (local_cap_idx - 1)*3 + icart)
               end do
            end do
         end do
      end if

   end subroutine redistribute_cap_dipole_derivatives

   subroutine check_system_geometry(sys_geom, error)
      !! Refuse a geometry with two atoms in the same place, before anything runs
      !!
      !! The per-fragment check finds this only once some fragment holding both
      !! copies is built, and distance screening can drop that fragment
      !! entirely. A duplicated atom also changes the electron count, so every
      !! number afterwards is wrong rather than merely unconverged.
      !!
      !! Sorted along the widest axis and swept, `O(N log N)` where the pairwise
      !! comparison would be 5e11 pairs at a million atoms. The axis is chosen
      !! by extent so a planar molecule does not fall back to quadratic; that is
      !! a cost choice, not a correctness one.
      use pic_sorting, only: sort_index
      use pic_types, only: int_index
      use pic_io, only: to_char

      type(system_geometry_t), intent(in) :: sys_geom
      type(error_t), intent(out) :: error

      real(dp), allocatable :: key(:)
      ! Pre-allocated and of `int_index`: `sort_index` writes into it rather
      ! than allocating, and handing it an unallocated array is a segmentation
      ! fault rather than an error. The permutation comes back 1-indexed.
      integer(int_index), allocatable :: order(:)
      real(dp) :: extent(3)
      real(dp) :: lo, hi, separation
      integer :: n, axis, widest, i, j, a, b
      ! TODO(mqc): the same threshold as `MIN_ATOM_DISTANCE`, which this module
      ! already imports and `check_duplicate_atoms` uses. Two spellings of one
      ! constant, which is what MQC002 exists to prevent.
      real(dp), parameter :: MIN_SEPARATION = 0.01_dp   !! Bohr, as the fragment check uses

      n = sys_geom%total_atoms
      if (n < 2) return

      do axis = 1, 3
         lo = minval(sys_geom%coordinates(axis, 1:n))
         hi = maxval(sys_geom%coordinates(axis, 1:n))
         extent(axis) = hi - lo
      end do
      widest = maxloc(extent, dim=1)

      allocate (key(n), order(n))
      key = sys_geom%coordinates(widest, 1:n)
      call sort_index(key, order)

      do i = 1, n - 1
         a = int(order(i))
         do j = i + 1, n
            b = int(order(j))
            ! Sorted, so once the gap along this axis exceeds the threshold no
            ! later atom can be within it either.
            if (sys_geom%coordinates(widest, b) - sys_geom%coordinates(widest, a) &
                > MIN_SEPARATION) exit
            separation = norm2(sys_geom%coordinates(:, b) - sys_geom%coordinates(:, a))
            if (separation < MIN_SEPARATION) then
               ! Set and not logged. This runs on every rank, while the
               ! caller prints the message on rank zero alone -- logging here
               ! would put the same two lines in front of that one once per
               ! rank, which on any real job buries the readable copy.
               call error%set(ERROR_VALIDATION, "atoms "//to_char(min(a, b))//" and "// &
                              to_char(max(a, b))//" are "//to_char(separation)// &
                              " Bohr apart, which is closer than any two nuclei can "// &
                              "be. The usual cause is an atom repeated in the input "// &
                              "geometry, which also makes the electron count wrong, "// &
                              "so this is refused before anything is computed rather "// &
                              "than left to surface as a fragment that will not "// &
                              "converge.")
               deallocate (key, order)
               return
            end if
         end do
      end do

      deallocate (key, order)
   end subroutine check_system_geometry

   subroutine check_duplicate_atoms(fragment, error)
      !! Refuse a fragment with two atoms closer than `MIN_ATOM_DISTANCE`
      !!
      !! Non-cap atoms only: a cap sits on the bond it closes, so it is
      !! legitimately close to the atom it replaces.
      ! TODO(mqc): pairwise, where `check_system_geometry` sorts and sweeps for
      ! the same test. Called once per fragment, so the quadratic cost lands on
      ! every term of the expansion.
      ! TODO(mqc): sets the error *and* logs it, on every rank. The message the
      ! caller prints on rank zero gets buried under one copy per rank -- the
      ! reason `check_system_geometry` deliberately sets without logging.
      use pic_logger, only: logger => global_logger
      use pic_io, only: to_char

      type(physical_fragment_t), intent(in) :: fragment
      type(error_t), intent(out) :: error

      integer :: i, j, n_atoms
      real(dp) :: distance, dx, dy, dz

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
                                 to_char(to_angstrom(distance))//" Angstrom)")
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
      if (allocated(this%cap_bonded_to)) deallocate (this%cap_bonded_to)
      if (allocated(this%local_to_global)) deallocate (this%local_to_global)
      if (allocated(this%is_ghost)) deallocate (this%is_ghost)
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

      ! Ghost centres carry basis functions and no nucleus, so they contribute
      ! no electrons either. Counted here rather than at each call site, since
      ! this is the only place the number is formed.
      if (allocated(this%is_ghost)) then
         nuclear_charge = sum(this%element_numbers, mask=.not. this%is_ghost)
      else
         nuclear_charge = sum(this%element_numbers)
      end if
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
      if (allocated(this%bonds)) deallocate (this%bonds)
      if (allocated(this%fragment_potentials)) deallocate (this%fragment_potentials)
      this%n_monomers = 0
      this%atoms_per_monomer = 0
      this%total_atoms = 0
   end subroutine system_destroy

   pure function calculate_monomer_distance(sys_geom, monomer_indices) result(min_distance)
      !! Closest approach between any two different monomers, in Angstrom
      !!
      !! Zero for a single monomer. This is what distance screening cuts on.
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
