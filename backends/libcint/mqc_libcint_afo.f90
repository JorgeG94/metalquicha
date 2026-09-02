!! Adjusted frozen orbitals: the small system a cut bond's orbital comes from
module mqc_libcint_afo
   !! A bond that fragmentation cuts has to be represented to both sides, and
   !! the way FMO does it is to freeze an orbital rather than cap the bond. This
   !! module builds the thing that orbital is derived from: a small model system
   !! around the cut, closed off with hydrogens, which is solved and localized so
   !! that the orbital sitting on the bond can be lifted out of it.
   !!
   !! **Why a model system rather than the fragment.** Capping a *fragment* is
   !! what was tried first and withdrawn: a cap belongs to a group, so a dimer
   !! spanning a cut bond must carry no cap there, and n-mers assembled from
   !! their members arrive holding caps that stand in the middle of a bond the
   !! dimer restored. A model system has none of that difficulty. It is built
   !! once per bond, is never assembled into anything, and contributes no energy
   !! -- only an orbital. The capping rule is the same one `cap_targets` already
   !! implements, applied where it is unambiguous.
   !!
   !! **The scale is derived rather than declared.** A deck's `cap_scale` is a
   !! choice about a fragment; here the cap only has to make a chemically
   !! sensible small molecule, so the hydrogen goes at the standard bond length
   !! for the atom it hangs off -- covalent radii summed -- as a fraction of the
   !! real internuclear distance. This is what GAMESS's automatic model does,
   !! with its own table of X-H lengths in place of the radius sum.
   use pic_types, only: dp
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_covalent_radius, element_number_to_symbol
   use mqc_physical_fragment, only: to_bohr
   use mqc_bond_perception, only: severed_bond_t, perceive_bonds, DEFAULT_BOND_TOLERANCE
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule, atom_ao_blocks
   use mqc_libcint_rhf, only: run_libcint_rhf, rhf_result_t
   use mqc_libcint_localize, only: boys_localize
   use mqc_scf_types, only: scf_numerics_t
   implicit none
   private

   public :: afo_model_t
   public :: afo_options_t
   public :: build_afo_model
   public :: bond_hybrid
   public :: afo_hybrid_t
   public :: build_group_frozen
   public :: cuts_outside_group
   public :: group_electron_shift
   public :: DEFAULT_MODEL_RADIUS
   public :: BOND_ORBITAL_REACH

   !! How far from either end of the cut bond the model reaches, in Angstrom.
   !!
   !! Wide enough to carry the bond's immediate chemical environment -- the
   !! substituents that decide what the hybrid looks like -- and no wider, since
   !! every atom past that is an SCF cost paid once per cut bond for an orbital
   !! that has stopped moving. GAMESS calls this RAFO and defaults it to a
   !! similar distance.
   real(dp), parameter :: DEFAULT_MODEL_RADIUS = 2.5_dp

   !! How close to the bond midpoint a localized orbital's centroid must sit to
   !! count as being *on* that bond, as a fraction of the bond length.
   !!
   !! **It must be well under a half, and that is structural rather than
   !! empirical.** Anything atom-centred -- a core orbital, a lone pair -- has
   !! its centroid on a nucleus, which is at exactly half the bond length from
   !! the midpoint. So a half admits every core on both atoms: measured on
   !! ethane the two carbon 1s orbitals come in at 1.4513013 Bohr against a
   !! half-bond of 1.4513096, inside it by eight decimal places, and a single
   !! sigma bond reports three orbitals on it and looks like a triple one.
   !!
   !! The measured spectrum for ethane, in Bohr: the C-C orbital at 8e-11, the
   !! two cores at 1.4513, the six C-H orbitals at 2.3546. A third of the bond
   !! length sits in the wide gap between the first two, with room for a polar
   !! bond whose orbital leans towards the electronegative end.
   real(dp), parameter :: BOND_ORBITAL_REACH = 0.35_dp

   !! What to solve the model system with
   type :: afo_options_t
      character(len=64) :: basis = "6-31g"
      type(scf_numerics_t) :: scf
         !! How the model system's SCF is driven. Its own type carried three
         !! fields long before this, and the caller filled none of them -- only
         !! `basis` was ever set -- so a model system converged on this type's
         !! defaults whatever the deck or the enclosing FMO run said.
         !! **Only the drive settings here are read** -- the accelerator, DIIS
         !! subspace, level shift, linear-dependence threshold and incremental
         !! Fock switch. Its `max_iter`, `energy_tol` and `density_tol` are
         !! NOT: the three bare fields below are, and they are passed
         !! positionally so they win. Two declarations of one concept in one
         !! type is a trap, and this comment is the guard rail until the bare
         !! three are folded in -- which cannot happen until their deliberately
         !! tighter values have somewhere else to live.
      integer :: scf_max_iter = 100
      real(dp) :: scf_energy_tol = 1.0e-10_dp
      real(dp) :: scf_density_tol = 1.0e-8_dp
   end type afo_options_t

   !! One cut bond's frozen orbital, over its bond-detached atom's functions
   type :: afo_hybrid_t
      !! Held one per cut bond rather than in a rectangular array, because two
      !! bonds can be detached at atoms of different elements and their hybrids
      !! are then different lengths.
      real(dp), allocatable :: coeff(:)
   end type afo_hybrid_t

   !! The small molecule a frozen orbital is derived from
   type :: afo_model_t
      !! Caps are the last `n_caps` entries of `z`, `sym` and `xyz`, and have no
      !! entry in `from_system`. That is the same layout `fragment_t` uses and
      !! for the same reason: an index into the system is meaningless for an
      !! atom the system does not contain.
      integer :: n_atoms = 0
      integer :: n_caps = 0
      integer, allocatable :: z(:)
      character(len=2), allocatable :: sym(:)
      real(dp), allocatable :: xyz(:, :)     !! (3, n_atoms), Bohr
      integer, allocatable :: from_system(:)  !! (n_atoms - n_caps), 1-based
      integer :: bda_local = 0               !! Where the bond's first atom sits
      integer :: baa_local = 0               !! Where its second atom sits
      integer :: nelec = 0
   end type afo_model_t

contains

   subroutine build_afo_model(z, coords, cut, model, error, radius, tolerance)
      !! The model system for one cut bond
      !!
      !! Both ends of the bond, everything within `radius` of either of them,
      !! every hydrogen hanging off what that selected, and a hydrogen cap for
      !! each bond that leaves the set.
      !!
      !! Hydrogens are pulled in wholesale rather than by distance because the
      !! alternative is absurd: a hydrogen just outside the radius would be
      !! replaced by a cap hydrogen a few hundredths of an Angstrom away, which
      !! is a different molecule for no reason and a discontinuity in anything
      !! that moves the geometry.
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      type(severed_bond_t), intent(in) :: cut
      type(afo_model_t), intent(out) :: model
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: radius
      real(dp), intent(in), optional :: tolerance

      logical, allocatable :: chosen(:)
      integer, allocatable :: order(:)
      real(dp) :: reach, tol
      integer :: n_atoms, i, j, n_real, n_caps, slot

      if (error%has_error()) return
      n_atoms = size(z)
      reach = to_bohr(DEFAULT_MODEL_RADIUS)
      if (present(radius)) reach = to_bohr(radius)
      tol = DEFAULT_BOND_TOLERANCE
      if (present(tolerance)) tol = tolerance

      if (cut%atom_a < 1 .or. cut%atom_a > n_atoms .or. &
          cut%atom_b < 1 .or. cut%atom_b > n_atoms) then
         call error%set(ERROR_VALIDATION, "afo model: the cut bond names an atom the "// &
                        "system does not have")
         return
      end if

      allocate (chosen(n_atoms), source=.false.)
      chosen(cut%atom_a) = .true.
      chosen(cut%atom_b) = .true.

      do i = 1, n_atoms
         if (chosen(i)) cycle
         if (near(coords, i, cut%atom_a, reach) .or. near(coords, i, cut%atom_b, reach)) then
            chosen(i) = .true.
         end if
      end do

      ! Hydrogens on anything already taken. One pass: a hydrogen brings nothing
      ! else in, so this cannot cascade.
      do i = 1, n_atoms
         if (chosen(i)) cycle
         if (z(i) /= 1) cycle
         do j = 1, n_atoms
            if (.not. chosen(j)) cycle
            if (bonded_pair(z, coords, i, j, tol)) then
               chosen(i) = .true.
               exit
            end if
         end do
      end do

      n_real = count(chosen)
      allocate (order(n_real))
      slot = 0
      do i = 1, n_atoms
         if (.not. chosen(i)) cycle
         slot = slot + 1
         order(slot) = i
      end do

      n_caps = 0
      do slot = 1, n_real
         i = order(slot)
         do j = 1, n_atoms
            if (chosen(j)) cycle
            if (bonded_pair(z, coords, i, j, tol)) n_caps = n_caps + 1
         end do
      end do

      model%n_atoms = n_real + n_caps
      model%n_caps = n_caps
      allocate (model%z(model%n_atoms), model%sym(model%n_atoms))
      allocate (model%xyz(3, model%n_atoms))
      allocate (model%from_system(n_real), source=order)

      do slot = 1, n_real
         i = order(slot)
         model%z(slot) = z(i)
         model%sym(slot) = element_number_to_symbol(z(i))
         model%xyz(:, slot) = coords(:, i)
         if (i == cut%atom_a) model%bda_local = slot
         if (i == cut%atom_b) model%baa_local = slot
      end do

      slot = n_real
      do i = 1, n_real
         do j = 1, n_atoms
            if (chosen(j)) cycle
            if (.not. bonded_pair(z, coords, order(i), j, tol)) cycle
            slot = slot + 1
            model%z(slot) = 1
            model%sym(slot) = element_number_to_symbol(1)
            model%xyz(:, slot) = cap_position(z, coords, order(i), j)
         end do
      end do

      model%nelec = sum(model%z)
      if (mod(model%nelec, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "afo model: the model system for the bond "// &
                        "between atoms "//to_char(cut%atom_a)//" and "// &
                        to_char(cut%atom_b)//" has an odd electron count, so the "// &
                        "capping did not close every valence it opened")
         return
      end if

      deallocate (chosen, order)
   end subroutine build_afo_model

   subroutine build_group_frozen(mol, bda_slot, occupied, hybrids, frozen, n_frozen_occ, error)
      !! Place each boundary's hybrid into this group's own basis
      !!
      !! A hybrid is stored over its bond-detached atom's functions alone, which
      !! is what makes it transferable. Putting it to work means writing it into
      !! that atom's block of whichever molecule is being solved and leaving the
      !! rest zero -- an index map, because the atom's functions are the same
      !! functions wherever it appears.
      !!
      !! **Occupied first, and that is not a formatting choice.** The constraint
      !! names its blocks by index range and `build_frozen_basis` orthonormalises
      !! them separately, so an occupied orbital sitting after a virtual one
      !! would be held at the level shift. Which of the two a boundary is comes
      !! from the assignment `$FMOBND` states: the fragment that gets all of the
      !! bond holds its hybrid occupied, the one that gets nothing of it holds
      !! the same hybrid empty.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: bda_slot(:)
         !! Where each boundary's bond-detached atom sits in `mol`, 1-based. It
         !! may be a ghost: a group holding the attached end carries the
         !! detached atom's functions without its nucleus.
      logical, intent(in) :: occupied(:)
      type(afo_hybrid_t), intent(in) :: hybrids(:)
      real(dp), allocatable, intent(out) :: frozen(:, :)
      integer, intent(out) :: n_frozen_occ
      type(error_t), intent(inout) :: error

      integer, allocatable :: offsets(:), counts(:)
      integer :: n, i, col

      if (error%has_error()) return
      n = size(bda_slot)
      if (size(occupied) /= n .or. size(hybrids) /= n) then
         call error%set(ERROR_VALIDATION, "afo: the boundaries of this group are "// &
                        "described by lists of different lengths")
         return
      end if

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)
      allocate (frozen(mol%nao, max(n, 1)), source=0.0_dp)
      n_frozen_occ = count(occupied)

      col = 0
      do i = 1, n
         if (.not. occupied(i)) cycle
         col = col + 1
         call place_hybrid(mol, offsets, counts, bda_slot(i), hybrids(i), &
                           frozen(:, col), error)
         if (error%has_error()) return
      end do
      do i = 1, n
         if (occupied(i)) cycle
         col = col + 1
         call place_hybrid(mol, offsets, counts, bda_slot(i), hybrids(i), &
                           frozen(:, col), error)
         if (error%has_error()) return
      end do
   end subroutine build_group_frozen

   subroutine place_hybrid(mol, offsets, counts, slot, hybrid, column, error)
      !! Write one hybrid into one atom's block of a group's basis
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: offsets(:), counts(:)
      integer, intent(in) :: slot
      type(afo_hybrid_t), intent(in) :: hybrid
      real(dp), intent(out) :: column(:)
      type(error_t), intent(inout) :: error

      integer :: first, last

      if (slot < 1 .or. slot > mol%natm) then
         call error%set(ERROR_VALIDATION, "afo: a boundary names an atom this group "// &
                        "does not contain")
         return
      end if
      if (.not. allocated(hybrid%coeff)) then
         call error%set(ERROR_VALIDATION, "afo: a boundary has no hybrid orbital to "// &
                        "freeze")
         return
      end if
      if (size(hybrid%coeff) /= counts(slot)) then
         call error%set(ERROR_VALIDATION, "afo: a hybrid has "// &
                        to_char(size(hybrid%coeff))//" coefficients but the atom it "// &
                        "belongs to has "//to_char(counts(slot))//" basis functions, "// &
                        "so it was built against a different basis set")
         return
      end if

      first = offsets(slot) + 1
      last = first + counts(slot) - 1
      column = 0.0_dp
      column(first:last) = hybrid%coeff
   end subroutine place_hybrid

   subroutine cuts_outside_group(cuts, n_cuts, members, outside, n_outside)
      !! Which cut bonds this n-mer is still cut across
      !!
      !! **The rule that the hydrogen caps got wrong.** A cut belongs to a
      !! *group*, not to a fragment. A bond severed between monomers I and J is
      !! whole again inside the dimer IJ, so that dimer must carry nothing
      !! standing in for it -- no cap, and no frozen orbital -- while still
      !! being cut against every fragment outside itself. Decide this per
      !! fragment and inherit it, and an n-mer arrives holding a stand-in for a
      !! bond it contains: 11.25 Hartree on a tripeptide, which is what
      !! `a605ee50d` withdrew.
      !!
      !! So it is computed from the group's own member list, every time, and
      !! never passed down from the members.
      type(severed_bond_t), intent(in) :: cuts(:)
      integer, intent(in) :: n_cuts
      integer, intent(in) :: members(:)   !! Fragment indices making up this n-mer
      integer, allocatable, intent(out) :: outside(:)  !! Indices into `cuts`
      integer, intent(out) :: n_outside

      integer :: i, count
      logical :: has_a, has_b

      count = 0
      do i = 1, n_cuts
         if (spans(cuts(i), members)) count = count + 1
      end do

      n_outside = count
      allocate (outside(max(count, 1)))
      count = 0
      do i = 1, n_cuts
         if (.not. spans(cuts(i), members)) cycle
         count = count + 1
         outside(count) = i
      end do
   end subroutine cuts_outside_group

   pure function spans(cut, members) result(crosses)
      !! Does this bond have exactly one end inside the group?
      !!
      !! Both ends inside means the group restored it. Neither means it is
      !! somebody else's boundary. Only one means the group's own edge.
      type(severed_bond_t), intent(in) :: cut
      integer, intent(in) :: members(:)
      logical :: crosses

      crosses = any(members == cut%frag_a) .neqv. any(members == cut%frag_b)
   end function spans

   pure function group_electron_shift(cuts, n_cuts, members) result(shift)
      !! How many electrons this group gains or loses at its boundaries
      !!
      !! GAMESS's `$FMOBND` documentation states the convention outright: of the
      !! pair `I J`, "the I-atom will get nothing of it, effectively remove one
      !! electron (1/2 of a single covalent bond) from its fragment. The J-atom
      !! will get all of the bond and thus adds one electron to its fragment".
      !!
      !! So a group holding the bond-detached end of a bond it does not contain
      !! is one electron short of the naive sum over its atoms, and one holding
      !! the attached end is one electron over. Summed over every fragment the
      !! shifts cancel, which is the check worth writing down: the electrons
      !! are moved between fragments, never created.
      type(severed_bond_t), intent(in) :: cuts(:)
      integer, intent(in) :: n_cuts
      integer, intent(in) :: members(:)
      integer :: shift

      integer :: i

      shift = 0
      do i = 1, n_cuts
         if (.not. spans(cuts(i), members)) cycle
         if (any(members == cuts(i)%frag_a)) then
            shift = shift - 1     ! the detached end gets nothing of the bond
         else
            shift = shift + 1     ! the attached end gets all of it
         end if
      end do
   end function group_electron_shift

   subroutine bond_hybrid(model, opts, hybrid, n_on_bond, error, centroid_distance)
      !! The orbital on the cut bond, expressed over the BDA's own basis functions
      !!
      !! Solve the model, localize its occupied orbitals, take the one sitting on
      !! the cut bond, and keep the part of it that lives on the atom the
      !! fragment will own.
      !!
      !! **Why only that atom's coefficients.** The frozen orbital has to be
      !! transferable into a fragment that does not contain the model system, and
      !! the one thing the two certainly share is the bond-detached atom and its
      !! basis functions. Restricting to that block turns the transfer into an
      !! index map rather than a projection between different molecules.
      !!
      !! **Normalisation carries across for free.** An atom's diagonal block of
      !! `S` is an integral over that atom's own functions and does not know what
      !! else is in the molecule, so it is bit-identical in the model and in the
      !! fragment. Normalising here is therefore normalising there, and nothing
      !! has to be renormalised on arrival.
      !!
      !! `n_on_bond` is the multiplicity question that geometry could not answer.
      !! A single sigma bond puts one localized orbital on the midpoint; a double
      !! bond puts two. Returned rather than acted on, because what to do about
      !! it -- refuse, or freeze both -- belongs to the caller. GAMESS makes the
      !! same choice available as `MODAFO`.
      type(afo_model_t), intent(in) :: model
      type(afo_options_t), intent(in) :: opts
      real(dp), allocatable, intent(out) :: hybrid(:)
      integer, intent(out) :: n_on_bond
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: centroid_distance(:)
         !! Every localized orbital's distance from the bond midpoint, in Bohr,
         !! so the cut between "on this bond" and "not" can be looked at rather
         !! than trusted.

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(scf_numerics_t) :: scf_numerics
      real(dp), allocatable :: localized(:, :), centroids(:, :), s(:, :), distance(:)
      integer, allocatable :: offsets(:), counts(:)
      real(dp) :: midpoint(3)
      real(dp) :: reach, bond_length, norm
      integer :: n_occ, i, best, first, last, n_bda

      if (error%has_error()) return
      if (model%bda_local < 1 .or. model%baa_local < 1) then
         call error%set(ERROR_VALIDATION, "afo: the model does not say where the cut "// &
                        "bond sits in it")
         return
      end if

      call build_libcint_molecule(model%z, model%sym, model%xyz, trim(opts%basis), &
                                  mol, error)
      if (error%has_error()) return

      ! `afo_options_t` carries the iteration count and the two tolerances twice
      ! -- once bare and once inside its `scf_numerics_t` -- and
      ! `run_libcint_rhf` reads only the positional ones, so setting
      ! `opts%scf%energy_tol` did nothing and said nothing. The bare fields are
      ! the ones callers and tests actually set, so they win; copying them over
      ! the numerics before the call means the two halves cannot disagree about
      ! what this SCF was asked for.
      scf_numerics = opts%scf
      scf_numerics%max_iter = opts%scf_max_iter
      scf_numerics%energy_tol = opts%scf_energy_tol
      scf_numerics%density_tol = opts%scf_density_tol
      call run_libcint_rhf(mol, model%nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                           opts%scf_density_tol, .false., scf, error, scf=scf_numerics)
      if (error%has_error()) return
      if (.not. scf%converged) then
         call error%set(ERROR_VALIDATION, "afo: the model system's SCF did not converge, "// &
                        "so there is no orbital to freeze")
         return
      end if

      n_occ = scf%n_occupied
      call boys_localize(mol, scf%orbitals, n_occ, localized, centroids, error)
      if (error%has_error()) return

      midpoint = 0.5_dp*(model%xyz(:, model%bda_local) + model%xyz(:, model%baa_local))
      bond_length = sqrt(sum((model%xyz(:, model%bda_local) &
                              - model%xyz(:, model%baa_local))**2))
      reach = BOND_ORBITAL_REACH*bond_length

      allocate (distance(n_occ))
      do i = 1, n_occ
         distance(i) = sqrt(sum((centroids(:, i) - midpoint)**2))
      end do

      best = minloc(distance, dim=1)
      n_on_bond = count(distance < reach)
      if (n_on_bond == 0) then
         call error%set(ERROR_VALIDATION, "afo: no localized orbital sits on the cut "// &
                        "bond, so the model system is not describing the bond it was "// &
                        "built for")
         return
      end if

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)
      first = offsets(model%bda_local) + 1
      n_bda = counts(model%bda_local)
      last = first + n_bda - 1

      call mol%overlap(s)
      allocate (hybrid(n_bda))
      hybrid = localized(first:last, best)

      norm = dot_product(hybrid, matmul(s(first:last, first:last), hybrid))
      if (norm <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "afo: the bond orbital has no weight on the "// &
                        "bond-detached atom, so there is no hybrid to take from it")
         return
      end if
      hybrid = hybrid/sqrt(norm)

      if (present(centroid_distance)) centroid_distance = distance
   end subroutine bond_hybrid

   pure function cap_position(z, coords, kept, gone) result(r)
      !! `R_H = R_kept + s (R_gone - R_kept)`, with `s` the standard bond length
      !!
      !! The same placement `cap_targets` uses, so a gradient through this cap
      !! has the weights that routine already returns. What differs is only
      !! where `s` comes from: there it is a declared fraction, here it is the
      !! covalent-radius sum over the real internuclear distance, which puts the
      !! hydrogen at a sensible bond length whatever the bond was.
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: kept, gone
      real(dp) :: r(3)

      real(dp) :: along(3)
      real(dp) :: distance, wanted, s

      along = coords(:, gone) - coords(:, kept)
      distance = sqrt(sum(along**2))
      wanted = to_bohr(element_covalent_radius(z(kept)) + element_covalent_radius(1))

      ! An element with no tabulated radius, or a contact so short the standard
      ! length is longer than it: leave the hydrogen where the atom was, which
      ! is the old unscaled behaviour and never worse than a negative scale.
      s = 1.0_dp
      if (wanted > 0.0_dp .and. distance > 0.0_dp) s = min(wanted/distance, 1.0_dp)

      r = coords(:, kept) + s*along
   end function cap_position

   pure function near(coords, i, j, reach) result(within)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: i, j
      real(dp), intent(in) :: reach
      logical :: within

      within = sum((coords(:, i) - coords(:, j))**2) <= reach*reach
   end function near

   pure function bonded_pair(z, coords, i, j, tol) result(is_bond)
      !! The same criterion `mqc_bond_perception` uses, on two atoms
      !!
      !! Restated here rather than reached for because that module's test takes
      !! a whole `system_geometry_t` and this is asked once per candidate pair
      !! inside two nested loops. The rule is one line and the tolerance is
      !! imported, so the two cannot drift on anything that matters.
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: i, j
      real(dp), intent(in) :: tol
      logical :: is_bond

      real(dp) :: radii, distance

      is_bond = .false.
      radii = element_covalent_radius(z(i)) + element_covalent_radius(z(j))
      if (radii <= 0.0_dp) return
      distance = sqrt(sum((coords(:, i) - coords(:, j))**2))
      is_bond = distance < tol*to_bohr(radii)
   end function bonded_pair

end module mqc_libcint_afo
