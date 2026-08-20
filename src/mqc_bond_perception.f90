!! Distance-based bond perception
module mqc_bond_perception
   !! Works out a molecule's connectivity from its geometry, so that a caller
   !! fragmenting a covalent system does not have to enumerate every bond by
   !! hand -- and, more usefully, so that a caller who did can be told what
   !! they left out.
   !!
   !! The rule is the usual one: two atoms are bonded when they are closer than
   !! the sum of their covalent radii, allowing a tolerance for the fact that
   !! real bonds stretch and tabulated radii are averages.
   !!
   !!     d(i,j) < TOLERANCE * (r_i + r_j)
   !!
   !! **This is a heuristic and is wrong sometimes.** It will invent a bond
   !! across a short hydrogen bond, miss a long dative one, and has no idea
   !! what a bond order is -- everything it finds is reported as single. It is
   !! a good default and a bad authority, which is why the important entry
   !! point here is `missing_broken_bonds`: it audits what a caller declared
   !! rather than replacing their judgement with its own.
   !!
   !! Radii are in Angstrom and coordinates in Bohr, so the comparison converts.
   !! An element with no tabulated radius bonds to nothing rather than being
   !! given a guessed one.
   use pic_types, only: dp
   use mqc_elements, only: element_covalent_radius
   use mqc_physical_fragment, only: system_geometry_t, bond_t, to_angstrom
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: perceive_bonds        !! Connectivity implied by the geometry
   public :: missing_broken_bonds  !! Cut bonds a caller's list does not mention
   public :: auto_monomers         !! Partition into covalently connected molecules
   public :: connected_components  !! Which covalently connected group each atom is in
   public :: monomer_of            !! Which monomer holds a 0-based atom
   public :: DEFAULT_BOND_TOLERANCE

   real(dp), parameter :: DEFAULT_BOND_TOLERANCE = 1.2_dp
      !! Slack on the radius sum. 1.2 is the common choice: loose enough for a
      !! stretched bond in an unrelaxed structure, tight enough that a hydrogen
      !! bond at 1.8 Angstrom is not mistaken for a covalent one.

contains

   pure function bonded(sys_geom, iatom, jatom, tolerance) result(is_bonded)
      !! Whether two atoms are within bonding distance
      !!
      !! Indices are 1-based here, into the geometry's own arrays; the bond
      !! lists this module produces are 0-based, like everything else.
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: iatom, jatom
      real(dp), intent(in) :: tolerance
      logical :: is_bonded

      real(dp) :: r_i, r_j, distance

      is_bonded = .false.
      r_i = element_covalent_radius(sys_geom%element_numbers(iatom))
      r_j = element_covalent_radius(sys_geom%element_numbers(jatom))
      ! No radius, no opinion.
      if (r_i <= 0.0_dp .or. r_j <= 0.0_dp) return

      distance = to_angstrom(norm2(sys_geom%coordinates(:, iatom) - &
                                   sys_geom%coordinates(:, jatom)))
      is_bonded = (distance < tolerance*(r_i + r_j))
   end function bonded

   subroutine perceive_bonds(sys_geom, bonds, n_bonds, tolerance)
      !! Every bond the geometry implies, with `is_broken` set from the partition
      !!
      !! A bond whose atoms fall in different monomers is one fragmentation
      !! cuts, and is marked so -- the same rule the JSON reader applies, which
      !! is what makes a perceived list usable directly.
      !!
      !! Order is always reported as 1. Nothing here can tell a double bond
      !! from a short single one, and capping only cares whether a bond breaks.
      type(system_geometry_t), intent(in) :: sys_geom
      type(bond_t), allocatable, intent(out) :: bonds(:)
      integer, intent(out) :: n_bonds
      real(dp), intent(in), optional :: tolerance

      real(dp) :: tol
      integer :: iatom, jatom, count

      tol = DEFAULT_BOND_TOLERANCE
      if (present(tolerance)) tol = tolerance

      ! Counted before allocating: the alternative is growing an array inside
      ! an O(N^2) loop, and this is the loop a large cluster spends real time in.
      count = 0
      do iatom = 1, sys_geom%total_atoms
         do jatom = iatom + 1, sys_geom%total_atoms
            if (bonded(sys_geom, iatom, jatom, tol)) count = count + 1
         end do
      end do

      n_bonds = count
      allocate (bonds(max(count, 1)))
      count = 0
      do iatom = 1, sys_geom%total_atoms
         do jatom = iatom + 1, sys_geom%total_atoms
            if (.not. bonded(sys_geom, iatom, jatom, tol)) cycle
            count = count + 1
            bonds(count)%atom_i = iatom - 1
            bonds(count)%atom_j = jatom - 1
            bonds(count)%order = 1
            bonds(count)%is_broken = &
               (monomer_of(sys_geom, iatom - 1) /= monomer_of(sys_geom, jatom - 1))
         end do
      end do
   end subroutine perceive_bonds

   subroutine missing_broken_bonds(sys_geom, declared, n_declared, missing_i, missing_j, &
                                   n_missing, tolerance)
      !! Bonds the geometry says are cut that the caller never mentioned
      !!
      !! The audit the declaration check cannot do on its own. Verifying that a
      !! declared bond is marked broken catches a caller who mis-labelled one;
      !! it cannot catch a caller who left one out entirely, because nothing in
      !! the declared list refers to it. This looks at the geometry instead.
      !!
      !! Reported rather than acted on. Perception is a heuristic, and a caller
      !! may have good reason to omit something it thinks is a bond -- so this
      !! hands back a list and lets the caller decide whether it is a mistake.
      type(system_geometry_t), intent(in) :: sys_geom
      type(bond_t), intent(in) :: declared(:)
      integer, intent(in) :: n_declared
      integer, allocatable, intent(out) :: missing_i(:), missing_j(:)
         !! 0-based atom indices, one entry per unmentioned cut
      integer, intent(out) :: n_missing
      real(dp), intent(in), optional :: tolerance

      type(bond_t), allocatable :: perceived(:)
      integer :: n_perceived, ibond, count

      call perceive_bonds(sys_geom, perceived, n_perceived, tolerance)

      count = 0
      do ibond = 1, n_perceived
         if (.not. perceived(ibond)%is_broken) cycle
         if (is_declared(declared, n_declared, perceived(ibond))) cycle
         count = count + 1
      end do

      n_missing = count
      allocate (missing_i(max(count, 1)), missing_j(max(count, 1)))
      missing_i = 0
      missing_j = 0

      count = 0
      do ibond = 1, n_perceived
         if (.not. perceived(ibond)%is_broken) cycle
         if (is_declared(declared, n_declared, perceived(ibond))) cycle
         count = count + 1
         missing_i(count) = perceived(ibond)%atom_i
         missing_j(count) = perceived(ibond)%atom_j
      end do

      deallocate (perceived)
   end subroutine missing_broken_bonds

   pure function is_declared(declared, n_declared, candidate) result(found)
      !! Whether a bond appears in the caller's list, either way round
      type(bond_t), intent(in) :: declared(:)
      integer, intent(in) :: n_declared
      type(bond_t), intent(in) :: candidate
      logical :: found

      integer :: ibond

      found = .false.
      do ibond = 1, n_declared
         if ((declared(ibond)%atom_i == candidate%atom_i .and. &
              declared(ibond)%atom_j == candidate%atom_j) .or. &
             (declared(ibond)%atom_i == candidate%atom_j .and. &
              declared(ibond)%atom_j == candidate%atom_i)) then
            found = .true.
            return
         end if
      end do
   end function is_declared

   subroutine connected_components(sys_geom, component, n_components, tolerance)
      !! Label each atom with the covalently connected group it belongs to
      !!
      !! Components are numbered 1..n in first-atom order, so a group's number
      !! reflects where it appears in the geometry rather than where the
      !! union-find happened to root it.
      !!
      !! Separate from `auto_monomers` because two callers want the labelling and
      !! only one of them wants a partition: an internal-coordinate scheme needs
      !! to know which atoms are bonded to which so it can build coordinates per
      !! group, and that is true whether or not the *energy* is being fragmented.
      !! Those are different questions about the same geometry.
      type(system_geometry_t), intent(in) :: sys_geom
      integer, allocatable, intent(out) :: component(:)
      integer, intent(out) :: n_components
      real(dp), intent(in), optional :: tolerance

      type(bond_t), allocatable :: bonds(:)
      integer :: n_bonds, iatom, ibond, root_i, root_j

      n_components = 0
      if (sys_geom%total_atoms <= 0) then
         allocate (component(0))
         return
      end if

      call perceive_bonds(sys_geom, bonds, n_bonds, tolerance)

      ! Union-find over the bond graph. Every atom starts as its own component;
      ! a bond merges two.
      allocate (component(sys_geom%total_atoms))
      do iatom = 1, sys_geom%total_atoms
         component(iatom) = iatom
      end do
      do ibond = 1, n_bonds
         root_i = find_root(component, bonds(ibond)%atom_i + 1)
         root_j = find_root(component, bonds(ibond)%atom_j + 1)
         if (root_i /= root_j) component(root_j) = root_i
      end do
      deallocate (bonds)

      ! Flatten to component ids, then renumber them 1..n in first-atom order.
      do iatom = 1, sys_geom%total_atoms
         component(iatom) = find_root(component, iatom)
      end do
      do iatom = 1, sys_geom%total_atoms
         if (component(iatom) /= iatom) cycle
         n_components = n_components + 1
         ! Temporarily negative so a renumbered id is not mistaken for a root.
         component(iatom) = -n_components
      end do
      do iatom = 1, sys_geom%total_atoms
         if (component(iatom) > 0) component(iatom) = component(component(iatom))
      end do
      component = -component
   end subroutine connected_components

   subroutine auto_monomers(sys_geom, error, tolerance)
      !! Make each covalently connected molecule a monomer
      !!
      !! The obvious partition for a cluster: perceive the bonds, take the
      !! connected components, and every water or every ligand becomes its own
      !! monomer. Hydrogen bonds are not perceived as covalent, so a hydrogen
      !! bonded cluster comes apart the way a chemist would expect.
      !!
      !! **It refuses on a single connected molecule**, which is the covalent
      !! case. There the components are a partition of one, which cannot be
      !! fragmented -- and the reason is not that the algorithm is too weak but
      !! that the question is genuinely the user's: where a molecule should be
      !! cut is chemistry, not connectivity. A peptide can be split per residue,
      !! per secondary structure, or somewhere else entirely, and nothing in a
      !! bond graph prefers one. So it says so rather than guessing, and rather
      !! than returning a single monomer that quietly makes fragmentation a no-op.
      !!
      !! Two connected molecules are a fine partition even if each is large --
      !! that is an ordinary MBE(2) on a dimer, and not this routine's business
      !! to second-guess.
      type(system_geometry_t), intent(inout) :: sys_geom
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: tolerance

      type(bond_t), allocatable :: bonds(:)
      integer, allocatable :: component(:), sizes(:)
      integer :: n_bonds, n_components, iatom, ibond, imon, largest
      integer :: root_i, root_j, slot

      if (sys_geom%total_atoms <= 0) then
         call error%set(ERROR_VALIDATION, "auto monomers: the system has no atoms")
         return
      end if

      call connected_components(sys_geom, component, n_components, tolerance)

      if (n_components < 2) then
         deallocate (component)
         call error%set(ERROR_VALIDATION, &
                        "auto monomers: the system is one covalently connected "// &
                        "molecule, so there is nothing to partition. Where to cut it "// &
                        "is a chemical choice -- supply the monomers explicitly, and "// &
                        "the bonds they break")
         return
      end if

      ! ---- fill the partition ----------------------------------------------
      allocate (sizes(n_components))
      sizes = 0
      do iatom = 1, sys_geom%total_atoms
         sizes(component(iatom)) = sizes(component(iatom)) + 1
      end do
      largest = maxval(sizes)

      if (allocated(sys_geom%fragment_sizes)) deallocate (sys_geom%fragment_sizes)
      if (allocated(sys_geom%fragment_atoms)) deallocate (sys_geom%fragment_atoms)
      if (allocated(sys_geom%fragment_charges)) deallocate (sys_geom%fragment_charges)
      if (allocated(sys_geom%fragment_multiplicities)) deallocate (sys_geom%fragment_multiplicities)
      allocate (sys_geom%fragment_sizes(n_components))
      allocate (sys_geom%fragment_atoms(largest, n_components))
      allocate (sys_geom%fragment_charges(n_components))
      allocate (sys_geom%fragment_multiplicities(n_components))

      sys_geom%fragment_atoms = 0
      sys_geom%fragment_sizes = 0
      ! Neutral singlets: connectivity says nothing about charge, and a caller
      ! who needs otherwise has to say so.
      sys_geom%fragment_charges = 0
      sys_geom%fragment_multiplicities = 1

      do iatom = 1, sys_geom%total_atoms
         imon = component(iatom)
         slot = sys_geom%fragment_sizes(imon) + 1
         sys_geom%fragment_sizes(imon) = slot
         sys_geom%fragment_atoms(slot, imon) = iatom - 1
      end do

      sys_geom%n_monomers = n_components
      if (all(sizes == sizes(1))) then
         sys_geom%atoms_per_monomer = sizes(1)
      else
         sys_geom%atoms_per_monomer = 0
      end if

      deallocate (component, sizes)
   end subroutine auto_monomers

   pure recursive function find_root(component, node) result(root)
      !! Representative of a node's component
      integer, intent(in) :: component(:)
      integer, intent(in) :: node
      integer :: root

      if (component(node) == node) then
         root = node
      else
         root = find_root(component, component(node))
      end if
   end function find_root

   pure function monomer_of(sys_geom, atom) result(imon)
      !! Which monomer holds a 0-based atom, or 0 if the partition omits it
      !!
      !! Zero for an unpartitioned atom is deliberate: two atoms nobody claimed
      !! compare equal and their bond is not called a cut, which is the right
      !! answer for a partition that does not cover the molecule.
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: atom
      integer :: imon

      integer :: jmon

      imon = 0
      if (.not. allocated(sys_geom%fragment_atoms)) return
      do jmon = 1, sys_geom%n_monomers
         if (any(sys_geom%fragment_atoms(1:sys_geom%fragment_sizes(jmon), jmon) == atom)) then
            imon = jmon
            return
         end if
      end do
   end function monomer_of

end module mqc_bond_perception
