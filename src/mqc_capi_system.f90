!! C entry points for describing the system to be fragmented
module mqc_capi_system
   !! The molecule, its monomer partition and its bonds, reachable from C.
   !!
   !! Wraps `system_geometry_t`, which already carries all three -- the atoms,
   !! the partition fragments are built from, and the connectivity hydrogen
   !! capping needs. Nothing new is modelled here; this is the way in.
   !!
   !! Three conventions that are invisible in the signatures and wrong by
   !! default if a caller assumes otherwise:
   !!
   !!   * **Coordinates arrive in Angstrom and are stored in Bohr.** Everything
   !!     inside metalquicha is atomic units; every input path converts. A
   !!     caller passing Bohr gets a molecule 1.89 times too large, which
   !!     converges to a perfectly plausible wrong energy.
   !!   * **Atom indices are 0-based**, in the monomer partition and in the
   !!     bonds alike, matching the JSON schema and `bond_t`. Monomer indices
   !!     in a *fragment term list* are 1-based -- different things, counted
   !!     differently, and the one trap in this interface.
   !!   * **Bonds must be declared, and broken ones must be declared broken.**
   !!     `set_bonds` is not optional: a system that never calls it is refused
   !!     rather than assumed to have none, because those two states are
   !!     indistinguishable from outside and only one of them is safe. Passing
   !!     `n_bonds = 0` is the way to say a system has nothing to cut.
   !!
   !!     Where the connectivity is given, every bond crossing a monomer
   !!     boundary must carry `is_broken`. That is checked here rather than
   !!     trusted, because the alternative is a fragment with a dangling
   !!     valence: an uncapped radical treated as a closed shell, which
   !!     converges to a plausible number and warns about nothing.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_char, c_null_char, &
                                                                             c_f_pointer, c_loc, c_associated
   use pic_types, only: dp
   use mqc_physical_fragment, only: system_geometry_t, bond_t, to_bohr
   use mqc_bond_perception, only: perceive_bonds, missing_broken_bonds, &
                                  auto_monomers, DEFAULT_BOND_TOLERANCE
   use mqc_error, only: error_t
   implicit none
   private

   public :: mqc_system_new, mqc_system_free
   public :: mqc_system_set_geometry, mqc_system_set_monomers, mqc_system_set_bonds
   public :: mqc_system_n_atoms, mqc_system_n_monomers, mqc_system_n_bonds
   public :: mqc_system_bonds_declared
   public :: mqc_system_perceive_bonds, mqc_system_count_missing_bonds
   public :: mqc_system_auto_monomers
   public :: mqc_system_last_error
   public :: system_handle_t
      !! Exposed for the sibling C-API modules that must open a system handle,
      !! not for general use. Nothing outside `src/interface` should hold one:
      !! the rest of the code works in `system_geometry_t`.

   integer(c_int), parameter :: MQC_OK = 0
   integer(c_int), parameter :: MQC_FAIL = 1
   integer(c_int), parameter :: MQC_BAD_HANDLE = 2

   integer, parameter :: MESSAGE_LEN = 512
   character(len=MESSAGE_LEN), save :: last_message = ""

   type :: system_handle_t
      !! The system, plus what the caller has actually told us about it
      !!
      !! `bonds_declared` is not the same as having bonds. A caller that has
      !! said "no bonds" and a caller that has not got round to it look
      !! identical in `system_geometry_t`, and one of them is about to
      !! fragment a covalent molecule into radicals.
      type(system_geometry_t) :: geom
      logical :: bonds_declared = .false.
   end type system_handle_t

contains

   function mqc_system_new() result(handle) bind(C, name="mqc_system_new")
      !! Allocate an empty system and return its handle
      type(c_ptr) :: handle

      type(system_handle_t), pointer :: h

      allocate (h)
      h%geom%n_monomers = 0
      h%geom%atoms_per_monomer = 0
      h%geom%total_atoms = 0
      h%geom%charge = 0
      h%geom%multiplicity = 1
      h%bonds_declared = .false.
      handle = c_loc(h)
   end function mqc_system_new

   subroutine mqc_system_free(handle) bind(C, name="mqc_system_free")
      !! Release a system. Safe on a null handle.
      type(c_ptr), value :: handle

      type(system_handle_t), pointer :: h

      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, h)
      call release(h%geom)
      deallocate (h)
   end subroutine mqc_system_free

   function mqc_system_set_geometry(handle, n_atoms, atomic_numbers, coordinates, &
                                    charge, multiplicity) result(status) &
      bind(C, name="mqc_system_set_geometry")
      !! The atoms, in Angstrom
      type(c_ptr), value :: handle
      integer(c_int), value :: n_atoms
      integer(c_int), intent(in) :: atomic_numbers(n_atoms)
      real(c_double), intent(in) :: coordinates(3*n_atoms)
         !! x,y,z per atom, contiguous, ANGSTROM
      integer(c_int), value :: charge
      integer(c_int), value :: multiplicity
      integer(c_int) :: status

      type(system_handle_t), pointer :: h
      integer :: iatom

      status = MQC_BAD_HANDLE
      if (.not. c_associated(handle)) then
         last_message = "null system handle"
         return
      end if
      call c_f_pointer(handle, h)

      if (n_atoms <= 0) then
         last_message = "mqc_system_set_geometry: a system needs at least one atom"
         status = MQC_FAIL
         return
      end if
      if (multiplicity < 1) then
         last_message = "mqc_system_set_geometry: multiplicity must be at least 1"
         status = MQC_FAIL
         return
      end if

      if (allocated(h%geom%element_numbers)) deallocate (h%geom%element_numbers)
      if (allocated(h%geom%coordinates)) deallocate (h%geom%coordinates)
      allocate (h%geom%element_numbers(n_atoms))
      allocate (h%geom%coordinates(3, n_atoms))

      h%geom%total_atoms = n_atoms
      h%geom%charge = charge
      h%geom%multiplicity = multiplicity
      do iatom = 1, n_atoms
         h%geom%element_numbers(iatom) = atomic_numbers(iatom)
         ! The one conversion. Internal units are Bohr throughout.
         h%geom%coordinates(1, iatom) = to_bohr(coordinates(3*(iatom - 1) + 1))
         h%geom%coordinates(2, iatom) = to_bohr(coordinates(3*(iatom - 1) + 2))
         h%geom%coordinates(3, iatom) = to_bohr(coordinates(3*(iatom - 1) + 3))
      end do

      status = MQC_OK
   end function mqc_system_set_geometry

   function mqc_system_set_monomers(handle, n_monomers, max_size, sizes, atoms, &
                                    charges, multiplicities) result(status) &
      bind(C, name="mqc_system_set_monomers")
      !! The partition fragments are built from
      !!
      !! `atoms` is (max_size, n_monomers) with one monomer per COLUMN, matching
      !! `fragment_atoms` -- the transpose of the term list's row-per-term
      !! layout, because that is what the two underlying types already use.
      type(c_ptr), value :: handle
      integer(c_int), value :: n_monomers
      integer(c_int), value :: max_size
      integer(c_int), intent(in) :: sizes(n_monomers)
      integer(c_int), intent(in) :: atoms(max_size*n_monomers)
         !! 0-based atom indices, column per monomer, unused slots ignored
      integer(c_int), intent(in) :: charges(n_monomers)
      integer(c_int), intent(in) :: multiplicities(n_monomers)
      integer(c_int) :: status

      type(system_handle_t), pointer :: h
      integer :: imon, iatom, base

      status = MQC_BAD_HANDLE
      if (.not. c_associated(handle)) then
         last_message = "null system handle"
         return
      end if
      call c_f_pointer(handle, h)

      if (h%geom%total_atoms <= 0) then
         last_message = "mqc_system_set_monomers: set the geometry first"
         status = MQC_FAIL
         return
      end if
      if (n_monomers <= 0 .or. max_size <= 0) then
         last_message = "mqc_system_set_monomers: monomer count and size must be positive"
         status = MQC_FAIL
         return
      end if
      do imon = 1, n_monomers
         if (sizes(imon) < 1 .or. sizes(imon) > max_size) then
            last_message = "mqc_system_set_monomers: a monomer size is outside 1..max_size"
            status = MQC_FAIL
            return
         end if
      end do

      if (allocated(h%geom%fragment_sizes)) deallocate (h%geom%fragment_sizes)
      if (allocated(h%geom%fragment_atoms)) deallocate (h%geom%fragment_atoms)
      if (allocated(h%geom%fragment_charges)) deallocate (h%geom%fragment_charges)
      if (allocated(h%geom%fragment_multiplicities)) deallocate (h%geom%fragment_multiplicities)
      allocate (h%geom%fragment_sizes(n_monomers))
      allocate (h%geom%fragment_atoms(max_size, n_monomers))
      allocate (h%geom%fragment_charges(n_monomers))
      allocate (h%geom%fragment_multiplicities(n_monomers))
      h%geom%fragment_atoms = 0

      do imon = 1, n_monomers
         base = (imon - 1)*max_size
         do iatom = 1, sizes(imon)
            if (atoms(base + iatom) < 0 .or. atoms(base + iatom) >= h%geom%total_atoms) then
               last_message = "mqc_system_set_monomers: an atom index is out of range "// &
                              "(indices are 0-based)"
               status = MQC_FAIL
               return
            end if
            h%geom%fragment_atoms(iatom, imon) = atoms(base + iatom)
         end do
         h%geom%fragment_sizes(imon) = sizes(imon)
         h%geom%fragment_charges(imon) = charges(imon)
         h%geom%fragment_multiplicities(imon) = multiplicities(imon)
      end do

      h%geom%n_monomers = n_monomers
      ! Zero means "variable-sized", which is what the rest of the code reads it
      ! as; only a uniform partition may claim a fixed size.
      if (all(sizes == sizes(1))) then
         h%geom%atoms_per_monomer = sizes(1)
      else
         h%geom%atoms_per_monomer = 0
      end if

      status = MQC_OK
   end function mqc_system_set_monomers

   function mqc_system_set_bonds(handle, n_bonds, atom_i, atom_j, order, is_broken) &
      result(status) bind(C, name="mqc_system_set_bonds")
      !! The connectivity hydrogen capping works from
      !!
      !! `is_broken` is non-zero for a bond fragmentation severs. Unlike the
      !! JSON path, nothing is derived: a caller supplying an arbitrary term
      !! list has no single partition to derive it from, so it says. Passing
      !! every entry as broken is the plain "here is my list of broken bonds".
      type(c_ptr), value :: handle
      integer(c_int), value :: n_bonds
      integer(c_int), intent(in) :: atom_i(n_bonds)  !! 0-based
      integer(c_int), intent(in) :: atom_j(n_bonds)  !! 0-based
      integer(c_int), intent(in) :: order(n_bonds)
      integer(c_int), intent(in) :: is_broken(n_bonds)  !! non-zero = broken
      integer(c_int) :: status

      type(system_handle_t), pointer :: h
      integer :: ibond

      status = MQC_BAD_HANDLE
      if (.not. c_associated(handle)) then
         last_message = "null system handle"
         return
      end if
      call c_f_pointer(handle, h)

      if (h%geom%total_atoms <= 0) then
         last_message = "mqc_system_set_bonds: set the geometry first"
         status = MQC_FAIL
         return
      end if
      if (h%geom%n_monomers <= 0) then
         ! The partition is what decides which bonds are cut, so it has to
         ! exist before the claim about them can be checked.
         last_message = "mqc_system_set_bonds: set the monomers first"
         status = MQC_FAIL
         return
      end if
      if (n_bonds < 0) then
         last_message = "mqc_system_set_bonds: bond count cannot be negative"
         status = MQC_FAIL
         return
      end if

      do ibond = 1, n_bonds
         if (atom_i(ibond) < 0 .or. atom_i(ibond) >= h%geom%total_atoms .or. &
             atom_j(ibond) < 0 .or. atom_j(ibond) >= h%geom%total_atoms) then
            last_message = "mqc_system_set_bonds: an atom index is out of range "// &
                           "(indices are 0-based)"
            status = MQC_FAIL
            return
         end if
         if (atom_i(ibond) == atom_j(ibond)) then
            last_message = "mqc_system_set_bonds: a bond joins an atom to itself"
            status = MQC_FAIL
            return
         end if
         if (order(ibond) < 1) then
            last_message = "mqc_system_set_bonds: bond order must be positive"
            status = MQC_FAIL
            return
         end if
      end do

      ! Every bond whose ends fall in different monomers is cut by the
      ! fragmentation, and must say so. Left unmarked, the fragment keeps a
      ! dangling valence: an uncapped radical run as a closed shell, which
      ! converges to a plausible number and warns about nothing. This is the
      ! one thing here worth failing hard over.
      do ibond = 1, n_bonds
         if (is_broken(ibond) /= 0) cycle
         if (monomer_of(h%geom, atom_i(ibond)) /= monomer_of(h%geom, atom_j(ibond))) then
            last_message = "mqc_system_set_bonds: the bond between atoms "// &
                           int_text(atom_i(ibond))//" and "//int_text(atom_j(ibond))// &
                           " crosses a monomer boundary but is not marked broken; "// &
                           "fragmenting it would leave an uncapped valence"
            status = MQC_FAIL
            return
         end if
      end do

      if (allocated(h%geom%bonds)) deallocate (h%geom%bonds)
      allocate (h%geom%bonds(max(n_bonds, 1)))
      do ibond = 1, n_bonds
         h%geom%bonds(ibond)%atom_i = atom_i(ibond)
         h%geom%bonds(ibond)%atom_j = atom_j(ibond)
         h%geom%bonds(ibond)%order = order(ibond)
         h%geom%bonds(ibond)%is_broken = (is_broken(ibond) /= 0)
      end do

      ! Even zero bonds is a statement, and the one this records.
      h%bonds_declared = .true.
      status = MQC_OK
   end function mqc_system_set_bonds

   function mqc_system_bonds_declared(handle) result(declared) &
      bind(C, name="mqc_system_bonds_declared")
      !! Whether the caller has said anything about bonds at all
      !!
      !! A run refuses to start without this. Saying "no bonds" is a fine
      !! answer; not having been asked is not, because the two are
      !! indistinguishable afterwards and only one of them is safe.
      type(c_ptr), value :: handle
      integer(c_int) :: declared

      type(system_handle_t), pointer :: h

      declared = 0_c_int
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, h)
      if (h%bonds_declared) declared = 1_c_int
   end function mqc_system_bonds_declared

   function mqc_system_auto_monomers(handle, tolerance) result(status) &
      bind(C, name="mqc_system_auto_monomers")
      !! Make each covalently connected molecule a monomer
      !!
      !! The right default for a cluster and refused for anything else. A
      !! single connected molecule has no automatic partition -- where to cut
      !! it is a chemical choice, and a bond graph has no opinion between per
      !! residue, per functional group, or something else. It fails rather than
      !! returning one monomer and quietly making fragmentation a no-op.
      !!
      !! So for a covalent system the monomers are mandatory, which is the
      !! contract this enforces rather than documents.
      type(c_ptr), value :: handle
      real(c_double), value :: tolerance
      integer(c_int) :: status

      type(system_handle_t), pointer :: h
      type(error_t) :: error
      real(dp) :: tol

      status = MQC_BAD_HANDLE
      if (.not. c_associated(handle)) then
         last_message = "null system handle"
         return
      end if
      call c_f_pointer(handle, h)

      if (h%geom%total_atoms <= 0) then
         last_message = "mqc_system_auto_monomers: set the geometry first"
         status = MQC_FAIL
         return
      end if

      tol = DEFAULT_BOND_TOLERANCE
      if (tolerance > 0.0_dp) tol = tolerance

      call auto_monomers(h%geom, error, tol)
      if (error%has_error()) then
         last_message = error%get_message()
         status = MQC_FAIL
         return
      end if
      status = MQC_OK
   end function mqc_system_auto_monomers

   function mqc_system_perceive_bonds(handle, tolerance) result(status) &
      bind(C, name="mqc_system_perceive_bonds")
      !! Work the connectivity out from the geometry and declare it
      !!
      !! Convenience, not authority. It uses covalent radii, so it will invent
      !! a bond across a short contact and miss a long one, and it calls every
      !! bond single. Good enough to save typing on an ordinary organic
      !! molecule; not a substitute for saying, which is why declaring bonds by
      !! hand remains the other way in.
      !!
      !! `tolerance` <= 0 asks for the default.
      type(c_ptr), value :: handle
      real(c_double), value :: tolerance
      integer(c_int) :: status

      type(system_handle_t), pointer :: h
      type(bond_t), allocatable :: found(:)
      integer :: n_found
      real(dp) :: tol

      status = MQC_BAD_HANDLE
      if (.not. c_associated(handle)) then
         last_message = "null system handle"
         return
      end if
      call c_f_pointer(handle, h)

      if (h%geom%total_atoms <= 0) then
         last_message = "mqc_system_perceive_bonds: set the geometry first"
         status = MQC_FAIL
         return
      end if
      if (h%geom%n_monomers <= 0) then
         last_message = "mqc_system_perceive_bonds: set the monomers first"
         status = MQC_FAIL
         return
      end if

      tol = DEFAULT_BOND_TOLERANCE
      if (tolerance > 0.0_dp) tol = tolerance

      call perceive_bonds(h%geom, found, n_found, tol)
      if (allocated(h%geom%bonds)) deallocate (h%geom%bonds)
      allocate (h%geom%bonds(max(n_found, 1)))
      if (n_found > 0) h%geom%bonds(1:n_found) = found(1:n_found)
      deallocate (found)

      ! Perception marks cuts from the partition itself, so what it produces is
      ! consistent by construction and needs no re-checking.
      h%bonds_declared = .true.
      status = MQC_OK
   end function mqc_system_perceive_bonds

   function mqc_system_count_missing_bonds(handle, tolerance) result(n_missing) &
      bind(C, name="mqc_system_count_missing_bonds")
      !! How many cuts the geometry implies that the declared bonds omit
      !!
      !! The audit `set_bonds` cannot perform. Checking that a declared bond is
      !! marked broken catches a mislabelled one; it cannot catch one left out
      !! entirely, because nothing in the list refers to it.
      !!
      !! Counted rather than enforced: perception is a heuristic and a caller
      !! may have a reason. Non-zero means look, not stop.
      !!
      !! A system with no declared bonds is answered rather than refused, and
      !! is the case most worth answering -- it is the caller who has not
      !! thought about bonds at all, so every cut the geometry implies is one
      !! they do not know about. `-1` therefore means only that the question
      !! could not be asked: a bad handle, no atoms, or no partition.
      type(c_ptr), value :: handle
      real(c_double), value :: tolerance
      integer(c_int) :: n_missing

      type(system_handle_t), pointer :: h
      type(bond_t), allocatable :: declared(:)
      integer, allocatable :: missing_i(:), missing_j(:)
      integer :: count
      real(dp) :: tol

      n_missing = -1_c_int
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, h)
      if (h%geom%total_atoms <= 0 .or. h%geom%n_monomers <= 0) return

      tol = DEFAULT_BOND_TOLERANCE
      if (tolerance > 0.0_dp) tol = tolerance

      if (allocated(h%geom%bonds)) then
         declared = h%geom%bonds
      else
         allocate (declared(0))
      end if
      call missing_broken_bonds(h%geom, declared, size(declared), &
                                missing_i, missing_j, count, tol)
      n_missing = int(count, c_int)
   end function mqc_system_count_missing_bonds

   pure function monomer_of(sys, atom) result(imon)
      !! Which monomer holds a 0-based atom, or 0 if the partition omits it
      !!
      !! Zero for an unpartitioned atom is deliberate: two atoms nobody claimed
      !! compare equal and their bond is not called a cut, which is the right
      !! answer for a partition that does not cover the molecule.
      type(system_geometry_t), intent(in) :: sys
      integer(c_int), intent(in) :: atom
      integer :: imon

      integer :: jmon

      imon = 0
      if (.not. allocated(sys%fragment_atoms)) return
      do jmon = 1, sys%n_monomers
         if (any(sys%fragment_atoms(1:sys%fragment_sizes(jmon), jmon) == int(atom))) then
            imon = jmon
            return
         end if
      end do
   end function monomer_of

   pure function int_text(value) result(text)
      integer(c_int), intent(in) :: value
      character(len=:), allocatable :: text
      character(len=16) :: buffer

      write (buffer, "(I0)") value
      text = trim(adjustl(buffer))
   end function int_text

   function mqc_system_n_atoms(handle) result(n) bind(C, name="mqc_system_n_atoms")
      !! Atoms in the system, or -1 for a bad handle
      type(c_ptr), value :: handle
      integer(c_int) :: n

      type(system_handle_t), pointer :: h

      n = -1_c_int
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, h)
      n = int(h%geom%total_atoms, c_int)
   end function mqc_system_n_atoms

   function mqc_system_n_monomers(handle) result(n) bind(C, name="mqc_system_n_monomers")
      !! Monomers in the partition, or -1 for a bad handle
      !!
      !! This is what a term list is generated over, so a caller reads it
      !! rather than tracking it separately.
      type(c_ptr), value :: handle
      integer(c_int) :: n

      type(system_handle_t), pointer :: h

      n = -1_c_int
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, h)
      n = int(h%geom%n_monomers, c_int)
   end function mqc_system_n_monomers

   function mqc_system_n_bonds(handle) result(n) bind(C, name="mqc_system_n_bonds")
      !! Bonds recorded, or -1 for a bad handle
      type(c_ptr), value :: handle
      integer(c_int) :: n

      type(system_handle_t), pointer :: h

      n = -1_c_int
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, h)
      n = 0_c_int
      if (allocated(h%geom%bonds)) n = int(size(h%geom%bonds), c_int)
   end function mqc_system_n_bonds

   subroutine mqc_system_last_error(buffer_len, buffer) &
      bind(C, name="mqc_system_last_error")
      !! Copy the most recent failure message out as a C string
      integer(c_int), value :: buffer_len
      character(kind=c_char), intent(inout) :: buffer(buffer_len)

      integer :: n, i

      if (buffer_len <= 0) return
      n = min(len_trim(last_message), int(buffer_len) - 1)
      do i = 1, n
         buffer(i) = last_message(i:i)
      end do
      buffer(n + 1) = c_null_char
   end subroutine mqc_system_last_error

   subroutine release(sys)
      !! Free everything the system owns
      type(system_geometry_t), intent(inout) :: sys

      if (allocated(sys%element_numbers)) deallocate (sys%element_numbers)
      if (allocated(sys%coordinates)) deallocate (sys%coordinates)
      if (allocated(sys%fragment_sizes)) deallocate (sys%fragment_sizes)
      if (allocated(sys%fragment_atoms)) deallocate (sys%fragment_atoms)
      if (allocated(sys%fragment_charges)) deallocate (sys%fragment_charges)
      if (allocated(sys%fragment_multiplicities)) deallocate (sys%fragment_multiplicities)
      if (allocated(sys%bonds)) deallocate (sys%bonds)
      sys%total_atoms = 0
      sys%n_monomers = 0
      sys%atoms_per_monomer = 0
   end subroutine release

end module mqc_capi_system
