!! Adjusted frozen orbitals: the small system a cut bond's orbital comes from
module mqc_czt_afo
   !! A bond that fragmentation cuts has to be represented to both sides, and
   !! FMO does it by freezing an orbital rather than capping the bond. This
   !! module builds the thing that orbital is derived from: a small model system
   !! around the cut, closed off with hydrogens, which is solved and localized so
   !! that the orbital sitting on the bond can be lifted out of it.
   !!
   !! A model system is built once per bond, is never assembled into anything,
   !! and contributes no energy -- only an orbital.
   !!
   !! A deck's `cap_scale` does not apply to these caps. The hydrogen goes at
   !! the standard bond length for the atom it hangs off -- covalent radii
   !! summed -- as a fraction of the real internuclear distance.
   use pic_types, only: dp
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_covalent_radius, element_number_to_symbol
   use mqc_physical_fragment, only: to_bohr
   use mqc_physical_constants, only: BOHR_TO_ANGSTROM
   use mqc_bond_perception, only: severed_bond_t, perceive_bonds, DEFAULT_BOND_TOLERANCE
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, atom_ao_blocks
   use mqc_czt_rhf, only: run_czt_rhf, rhf_result_t
   use mqc_czt_localize, only: boys_localize, er_localize, LOCALIZER_BOYS, LOCALIZER_ER
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
   public :: peptide_bond_advice
   public :: afo_lmo_set_t
   public :: orient_cut
   public :: build_bonded_model
   public :: bond_lmo_set
   public :: build_group_frozen_set
   public :: lmo_set_pack_size
   public :: lmo_set_pack
   public :: lmo_set_unpack
   public :: DEFAULT_MODEL_RADIUS
   public :: BOND_ORBITAL_REACH

   real(dp), parameter :: DEFAULT_MODEL_RADIUS = 2.5_dp
      !! How far from either end of the cut bond the model reaches, in Angstrom.
      !! Wide enough to carry the bond's immediate chemical environment; GAMESS
      !! calls this RAFO and defaults it to a similar distance.

   real(dp), parameter :: BOND_ORBITAL_REACH = 0.35_dp
      !! How close to the bond midpoint a localized orbital's centroid must sit
      !! to count as being *on* that bond, as a fraction of the bond length.
      !!
      !! **It must be well under a half.** Anything atom-centred -- a core
      !! orbital, a lone pair -- has its centroid on a nucleus, which is at
      !! exactly half the bond length from the midpoint, so a half admits every
      !! core on both atoms and a single sigma bond looks like a triple one.

   type :: afo_options_t
      !! What to solve the model system with
      character(len=64) :: basis = ""
         !! **Empty on purpose, and refused rather than defaulted.** This field
         !! used to start at "6-31g", which no run ever saw: every caller
         !! overwrites it from the deck, and a deck that omits `model.basis`
         !! gets "sto-3g" from `mqc_method_config`. So the initialiser named a
         !! basis nothing was ever computed in, which is worse than no default
         !! at all -- a plumbing bug that lost the deck's basis would have
         !! silently produced 6-31G numbers.
      type(scf_numerics_t) :: scf
         !! How the model system's SCF is driven. **Only the drive settings are
         !! read** -- the accelerator, DIIS subspace, level shift,
         !! linear-dependence threshold and incremental Fock switch.
      ! TODO(mqc): `scf%max_iter`, `scf%energy_tol` and `scf%density_tol` are
      ! ignored -- the three bare fields below are passed positionally and win.
      ! Two declarations of one concept in one type, so setting the `scf` ones
      ! has no effect and no complaint.
      integer :: scf_max_iter = 200
      real(dp) :: scf_energy_tol = 1.0e-11_dp
      real(dp) :: scf_density_tol = 1.0e-9_dp
      real(dp) :: scf_grad_tol = 1.0e-9_dp
         !! The model's own convergence, deliberately not the fragments'. What
         !! leaves the model is its orbitals, which are first order in the
         !! commutator, so that is bounded outright rather than derived as
         !! `sqrt(scf_energy_tol)`; and a model is a dozen atoms, so the extra
         !! iterations cost nothing. Bounding it at 1e-5, which is what the
         !! fragment tolerance used to hand down, left the glycine
         !! tripeptide's ER monomers 2e-7 short.
      logical :: cartesian = .false.
         !! Build the model system Cartesian whatever the basis declares. The
         !! hybrid is transferred through the detached atom's block of
         !! functions, so the model has to be built in the angular form the
         !! fragment is: EFMO's are Cartesian throughout, FMO's follow the basis.
      character(len=8) :: localization = LOCALIZER_ER
         !! How the model's occupied orbitals are localized: "er"
         !! (Edmiston-Ruedenberg, GAMESS's default for the model) or "boys".
      logical :: show_scf = .false.
         !! Print the model system's SCF table, at the caller's verbose level
   end type afo_options_t

   integer, parameter :: GAMESS_MAX_Z = 86
      !! Elements GAMESS's `PAIRBOND` tabulates
   integer, parameter :: GAMESS_MAX_CAP_Z = 17
      !! Elements GAMESS's cap length table reaches

   real(dp), parameter :: GAMESS_RCOV(GAMESS_MAX_Z) = [ &
                          0.30_dp, 1.22_dp, &
                          1.23_dp, 0.89_dp, 0.88_dp, 0.77_dp, 0.70_dp, 0.66_dp, 0.58_dp, 1.60_dp, &
                          1.66_dp, 1.36_dp, 1.25_dp, 1.17_dp, 1.10_dp, 1.04_dp, 0.99_dp, 1.91_dp, &
                          2.03_dp, 1.74_dp, 1.44_dp, 1.32_dp, 1.22_dp, 1.19_dp, 1.17_dp, 1.165_dp, &
                          1.16_dp, 1.15_dp, 1.17_dp, 1.25_dp, 1.25_dp, 1.22_dp, 1.21_dp, 1.17_dp, &
                          1.14_dp, 1.98_dp, &
                          2.22_dp, 1.92_dp, 1.62_dp, 1.45_dp, 1.34_dp, 1.29_dp, 1.27_dp, 1.24_dp, &
                          1.25_dp, 1.28_dp, 1.34_dp, 1.41_dp, 1.50_dp, 1.40_dp, 1.41_dp, 1.37_dp, &
                          1.33_dp, 2.09_dp, &
                          2.35_dp, 1.98_dp, 1.69_dp, 1.65_dp, 1.65_dp, 1.64_dp, 1.65_dp, 1.66_dp, &
                          1.65_dp, 1.61_dp, 1.59_dp, 1.59_dp, 1.58_dp, 1.57_dp, 1.56_dp, 1.56_dp, &
                          1.56_dp, 1.44_dp, 1.34_dp, 1.30_dp, 1.28_dp, 1.26_dp, 1.26_dp, 1.29_dp, &
                          1.34_dp, 1.44_dp, 1.55_dp, 1.54_dp, 1.52_dp, 1.53_dp, 1.50_dp, 2.20_dp]
      !! GAMESS's `PAIRBOND` radii, Angstrom (Emsley, with its own guesses).
      !! Not `covalent_radius_emsley`, which is a different Emsley table: this
      !! one decides which atoms a GAMESS model system holds, and matching
      !! GAMESS's frozen orbitals means matching its model.

   real(dp), parameter :: GAMESS_CAP_LENGTH(GAMESS_MAX_CAP_Z) = [ &
                          0.74_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.19_dp, 1.09_dp, 1.01_dp, 0.96_dp, &
                          0.92_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.48_dp, 1.44_dp, 1.34_dp, 1.27_dp]
      !! X-H length a cap hydrogen is put at, Angstrom, by the element it hangs
      !! off; GAMESS's `RH`. Zero means untabulated, and 1.6 is used, as there.

   type :: afo_hybrid_t
      !! One cut bond's frozen orbital, over its bond-detached atom's functions
      !!
      !! Held one per cut bond rather than in a rectangular array: two bonds can
      !! be detached at atoms of different elements, and their hybrids are then
      !! different lengths.
      real(dp), allocatable :: coeff(:)
   end type afo_hybrid_t

   type :: afo_lmo_set_t
      !! One cut bond's frozen orbitals: the model's localized orbitals that
      !! belong to the bond-detached atom, over the real atoms near the bond
      !!
      !! GAMESS's adjusted frozen orbitals (`fmolib.src` AFO construction and
      !! `flmovec`). Orbital 1 is the bond's own, frozen occupied in the
      !! fragment carrying the detached atom as a ghost and frozen empty in the
      !! fragment that owns it. Orbitals 2 on -- the atom's core and its other
      !! bonds -- are frozen empty in the ghost-holding fragment only, which is
      !! what keeps that fragment from using the borrowed functions for
      !! anything but the bond pair. A carbon has five: 1 + 4.
      !!
      !! Each orbital is kept on every real model atom bonded to either end of
      !! the bond, and truncated to whichever of those a group holds when it
      !! is placed. It is not renormalised: the group orthonormalises the set.
      integer :: n_lmo = 0
      integer :: n_at = 0
      integer, allocatable :: atoms(:)       !! (n_at), system indices
      integer, allocatable :: n_func(:)      !! (n_at), functions on each
      real(dp), allocatable :: coeff(:, :)   !! (sum n_func, n_lmo), atom by atom
      real(dp), allocatable :: pop(:, :)
         !! (n_at, n_lmo), each orbital's population on each atom, `c_A S_AA c_A`
   end type afo_lmo_set_t

   type :: afo_model_t
      !! The small molecule a frozen orbital is derived from
      !!
      !! Caps are the last `n_caps` entries of `z`, `sym` and `xyz`, and have no
      !! entry in `from_system`.
      integer :: n_atoms = 0
      integer :: n_caps = 0
      integer, allocatable :: z(:)
      character(len=2), allocatable :: sym(:)
      real(dp), allocatable :: xyz(:, :)     !! (3, n_atoms), Bohr
      integer, allocatable :: from_system(:)  !! (n_atoms - n_caps), 1-based
      integer :: bda_local = 0               !! Where the bond's first atom sits
      integer :: baa_local = 0               !! Where its second atom sits
      integer :: charge = 0
         !! Net charge of the charged groups the sphere took in whole; see
         !! `model_formal_charge`
      integer :: nelec = 0
   end type afo_model_t

contains

   subroutine build_afo_model(z, coords, cut, model, error, radius, tolerance)
      !! The model system for one cut bond
      !!
      !! Both ends of the bond, everything within `radius` of either of them,
      !! every singly bonded atom hanging off what that selected, and a
      !! hydrogen cap for each bond that leaves the set. `radius` is in
      !! Angstrom and defaults to `DEFAULT_MODEL_RADIUS`.
      !!
      !! A singly bonded neighbour comes in wholesale rather than by distance,
      !! so that one just outside the radius is not replaced by a cap hydrogen
      !! a few hundredths of an Angstrom away -- a discontinuity in anything
      !! that moves the geometry. For a heavy one there is a second reason: a
      !! cap hydrogen closes one electron pair, and a carbonyl oxygen is held
      !! by two, so capping it would leave the model a radical.
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

      ! Terminal heavy atoms on anything already taken -- a carbonyl oxygen
      ! above all. A cap hydrogen stands in for exactly one electron pair, so
      ! replacing a doubly bonded oxygen by one leaves a valence open and the
      ! model comes back a radical. Perception here is distance-based and
      ! cannot report a bond order, so the order is not guessed: the atom comes
      ! in whole instead. One pass, for the reason the hydrogens are one pass --
      ! an atom with a single bond brings nothing else with it -- and it keeps
      ! the model continuous in the geometry for the same reason too.
      do i = 1, n_atoms
         if (chosen(i)) cycle
         if (z(i) == 1) cycle
         if (count_bonds(z, coords, i, tol) /= 1) cycle
         do j = 1, n_atoms
            if (.not. chosen(j)) cycle
            if (bonded_pair(z, coords, i, j, tol)) then
               chosen(i) = .true.
               exit
            end if
         end do
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

      model%charge = model_formal_charge(z, coords, chosen, tol)
      model%nelec = sum(model%z) - model%charge
      if (mod(model%nelec, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "afo model: the model system for the bond "// &
                        "between atoms "//to_char(cut%atom_a)//" and "// &
                        to_char(cut%atom_b)//" has an odd electron count, so the "// &
                        "capping did not close every valence it opened. A cap "// &
                        "hydrogen closes one electron pair, so the sphere around "// &
                        "the bond has clipped a multiple bond at an atom that is "// &
                        "not terminal, has taken in a charged group this builder "// &
                        "does not recognise, or the system is not a closed shell "// &
                        "to begin with")
         return
      end if

      deallocate (chosen, order)
   end subroutine build_afo_model

   subroutine build_group_frozen(mol, bda_slot, occupied, hybrids, frozen, n_frozen_occ, error)
      !! Place each boundary's hybrid into this group's own basis
      !!
      !! A hybrid is stored over its bond-detached atom's functions alone, so
      !! putting it to work is an index map: write it into that atom's block of
      !! whichever molecule is being solved and leave the rest zero.
      !!
      !! **Columns come back occupied first.** The constraint names its blocks
      !! by index range and `build_frozen_basis` orthonormalises them
      !! separately, so an occupied orbital sitting after a virtual one would be
      !! held at the level shift. The fragment that gets all of the bond holds
      !! its hybrid occupied, the one that gets nothing of it holds the same
      !! hybrid empty.
      type(czt_molecule_t), intent(in) :: mol
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
      type(czt_molecule_t), intent(in) :: mol
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
      !! **A cut belongs to a group, not to a fragment.** A bond severed between
      !! monomers I and J is whole again inside the dimer IJ, so that dimer
      !! carries nothing standing in for it -- no cap, and no frozen orbital --
      !! while still being cut against every fragment outside itself.
      !!
      !! Computed from the group's own member list every time, and never passed
      !! down from the members.
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
      !! By the `$FMOBND` convention, of the pair `I J` the I-atom gets nothing
      !! of the bond and the J-atom gets all of it. So a group holding the
      !! bond-detached end of a bond it does not contain is one electron short
      !! of the naive sum over its atoms, and one holding the attached end is
      !! one electron over. Summed over every fragment the shifts cancel.
      !!
      !! **The electron moves either way; the nucleus may or may not follow.**
      !! Where it does, a unit of nuclear charge crosses the same boundary and
      !! the fragment comes out neutral rather than charged -- but that half
      !! is applied per atom rather than per group, and only under a field;
      !! see `nuc_charge` on `group_t` and `splits_nucleus` in
      !! [[mqc_czt_fmo]]. What is counted here is the electron alone, which is
      !! the `$FMOBND` assignment and not a convention.
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
      !! Solve the model, localize its occupied orbitals, take the one whose
      !! centroid sits on the cut bond, and keep the part of it that lives on
      !! the bond-detached atom -- the one block the model and the fragment
      !! certainly share, which is what makes the hybrid transferable.
      !!
      !! Normalised against that atom's diagonal block of `S`, which is an
      !! integral over the atom's own functions and so is identical in the model
      !! and in the fragment. Nothing is renormalised on arrival.
      !!
      !! `n_on_bond` counts the orbitals within `BOND_ORBITAL_REACH` of the
      !! midpoint: one for a single sigma bond, two for a double. Returned
      !! rather than acted on -- whether to refuse or to freeze both belongs to
      !! the caller.
      type(afo_model_t), intent(in) :: model
      type(afo_options_t), intent(in) :: opts
      real(dp), allocatable, intent(out) :: hybrid(:)
      integer, intent(out) :: n_on_bond
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: centroid_distance(:)
         !! Every localized orbital's distance from the bond midpoint, in Bohr,
         !! so the cut between "on this bond" and "not" can be looked at rather
         !! than trusted.

      type(czt_molecule_t) :: mol
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
      if (len_trim(opts%basis) == 0) then
         call error%set(ERROR_VALIDATION, "afo: no orbital basis was named for the "// &
                        "model system. It has to be the one the fragment is solved "// &
                        "in -- the hybrid is transferred between the two through a "// &
                        "shared atomic block -- so there is no basis to default to "// &
                        "here.")
         return
      end if

      call build_czt_molecule(model%z, model%sym, model%xyz, trim(opts%basis), &
                              mol, error, force_cartesian=opts%cartesian)
      if (error%has_error()) return

      ! `afo_options_t` carries the iteration count and the two tolerances twice
      ! -- once bare and once inside its `scf_numerics_t` -- and
      ! `run_czt_rhf` reads only the positional ones, so setting
      ! `opts%scf%energy_tol` did nothing and said nothing. The bare fields are
      ! the ones callers and tests actually set, so they win; copying them over
      ! the numerics before the call means the two halves cannot disagree about
      ! what this SCF was asked for.
      scf_numerics = opts%scf
      scf_numerics%max_iter = opts%scf_max_iter
      scf_numerics%energy_tol = opts%scf_energy_tol
      scf_numerics%density_tol = opts%scf_density_tol
      scf_numerics%grad_tol = opts%scf_grad_tol
      call run_czt_rhf(mol, model%nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                       opts%scf_density_tol, opts%show_scf, scf, error, scf=scf_numerics, &
                       grad_tol=opts%scf_grad_tol)
      if (error%has_error()) return
      if (.not. scf%converged) then
         call error%set(ERROR_VALIDATION, "afo: the model system's SCF did not converge "// &
                        "in "//to_char(scf%iterations)//" iterations, orbital gradient "// &
                        to_char(scf%commutator)//" at the last, so there is no orbital "// &
                        "to freeze")
         return
      end if

      n_occ = scf%n_occupied
      call localize_model(mol, scf%orbitals, n_occ, opts, localized, centroids, error)
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

   subroutine localize_model(mol, orbitals, n_occ, opts, localized, centroids, error)
      !! The model's occupied orbitals, localized as `opts%localization` asks
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_occ
      type(afo_options_t), intent(in) :: opts
      real(dp), allocatable, intent(out) :: localized(:, :)   !! (n_ao, n_occ)
      real(dp), allocatable, intent(out) :: centroids(:, :)   !! (3, n_occ), Bohr
      type(error_t), intent(inout) :: error

      select case (trim(opts%localization))
      case (LOCALIZER_ER)
         call er_localize(mol, orbitals, n_occ, localized, centroids, error)
      case (LOCALIZER_BOYS)
         call boys_localize(mol, orbitals, n_occ, localized, centroids, error)
      case default
         call error%set(ERROR_VALIDATION, "afo: unknown localization '"// &
                        trim(opts%localization)//"'; expected 'er' or 'boys'")
      end select
   end subroutine localize_model

   pure function cap_position(z, coords, kept, gone) result(r)
      !! `R_H = R_kept + s (R_gone - R_kept)`, with `s` the standard bond length
      !!
      !! The same placement `cap_targets` uses, so a gradient through this cap
      !! has the weights that routine returns. Only `s` differs: here it is the
      !! covalent-radius sum over the real internuclear distance, capped at 1.
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
      ! length is longer than it: leave the hydrogen where the atom was.
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

   function peptide_bond_advice(z, coords, atom_a, atom_b, tolerance) result(advice)
      !! What to cut instead, when the bond named is a backbone peptide bond
      !!
      !! Empty unless the two atoms are a carbon and a nitrogen and that carbon
      !! also carries a terminal oxygen, which on a protein is the amide of the
      !! backbone. FMO's convention there is to leave the peptide bond whole
      !! and cut the C-alpha--C(=O) bond one place along instead, so that the
      !! detached bond is a nonpolar single one outside the carbonyl's
      !! conjugation. Named by index where the carbonyl has exactly one carbon
      !! neighbour to name.
      !!
      !! Atoms are numbered from one, as everywhere in this backend and one
      !! more than the deck's own numbering.
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: atom_a, atom_b
      real(dp), intent(in), optional :: tolerance
      character(len=:), allocatable :: advice

      real(dp) :: tol
      integer :: carbon, nitrogen, j, alpha, n_alpha
      logical :: carbonyl

      advice = ""
      tol = DEFAULT_BOND_TOLERANCE
      if (present(tolerance)) tol = tolerance

      if (.not. bonded_pair(z, coords, atom_a, atom_b, tol)) return
      if (z(atom_a) == 6 .and. z(atom_b) == 7) then
         carbon = atom_a
         nitrogen = atom_b
      else if (z(atom_a) == 7 .and. z(atom_b) == 6) then
         carbon = atom_b
         nitrogen = atom_a
      else
         return
      end if

      carbonyl = .false.
      do j = 1, size(z)
         if (j == carbon) cycle
         if (z(j) /= 8) cycle
         if (.not. bonded_pair(z, coords, carbon, j, tol)) cycle
         if (count_bonds(z, coords, j, tol) == 1) carbonyl = .true.
      end do
      if (.not. carbonyl) return

      n_alpha = 0
      alpha = 0
      do j = 1, size(z)
         if (j == carbon) cycle
         if (z(j) /= 6) cycle
         if (.not. bonded_pair(z, coords, carbon, j, tol)) cycle
         n_alpha = n_alpha + 1
         alpha = j
      end do

      advice = " Atoms "//to_char(min(nitrogen, carbon))//" and "// &
               to_char(max(nitrogen, carbon))//" are a "// &
               "peptide bond -- an amide C(=O)-N, conjugated with the carbonyl "// &
               "and not a plain single bond. The FMO convention for a protein is "// &
               "to leave it whole and cut the C-alpha--C(=O) bond one place along "// &
               "the backbone instead"
      if (n_alpha == 1) then
         advice = advice//", which here is the bond between atoms "// &
                  to_char(min(alpha, carbon))//" and "// &
                  to_char(max(alpha, carbon))//"."
      else
         advice = advice//"."
      end if
   end function peptide_bond_advice

   pure function model_formal_charge(z, coords, chosen, tol) result(q)
      !! The net charge of the ionised groups among the chosen atoms
      !!
      !! A model system is closed with neutral caps, so any charge it holds is
      !! a group it took in whole: a peptide's N-terminal ammonium sits inside
      !! the sphere around the first C-alpha--C(=O) bond, and a C-terminal
      !! carboxylate inside the last one. Without the charge such a model has
      !! an odd electron count and no hybrid can be taken off it.
      !!
      !! Recognised, from neighbour counts alone since perception here reports
      !! no bond orders: a nitrogen with four neighbours (ammonium, +1); a
      !! carbon whose three neighbours are nitrogens with three neighbours each
      !! (guanidinium, +1); and a carbon with three neighbours, two of them
      !! terminal oxygens (carboxylate, -1). A neutral carboxylic acid has a
      !! hydrogen on one oxygen, so that oxygen is not terminal and the group
      !! is not counted.
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      logical, intent(in) :: chosen(:)
      real(dp), intent(in) :: tol
      integer :: q

      integer :: i, j, n_n3, n_o1

      q = 0
      do i = 1, size(z)
         if (.not. chosen(i)) cycle
         select case (z(i))
         case (7)
            if (count_bonds(z, coords, i, tol) == 4) q = q + 1
         case (6)
            if (count_bonds(z, coords, i, tol) /= 3) cycle
            n_n3 = 0
            n_o1 = 0
            do j = 1, size(z)
               if (j == i) cycle
               if (.not. bonded_pair(z, coords, i, j, tol)) cycle
               if (z(j) == 7 .and. count_bonds(z, coords, j, tol) == 3) n_n3 = n_n3 + 1
               if (z(j) == 8 .and. chosen(j) .and. count_bonds(z, coords, j, tol) == 1) then
                  n_o1 = n_o1 + 1
               end if
            end do
            if (n_n3 == 3) q = q + 1
            if (n_o1 == 2) q = q - 1
         case default
            ! Nothing else is recognised; a charged group of another element
            ! leaves the count odd and the model is refused.
         end select
      end do
   end function model_formal_charge

   pure function count_bonds(z, coords, i, tol) result(n)
      !! How many atoms `i` is bonded to, by the same distance criterion
      !!
      !! One means terminal: a hydrogen, or the oxygen of a carbonyl. Since
      !! perception is distance-based the count is of *neighbours*, not of
      !! electron pairs, so a doubly bonded oxygen counts one.
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: i
      real(dp), intent(in) :: tol
      integer :: n

      integer :: j

      n = 0
      do j = 1, size(z)
         if (j == i) cycle
         if (bonded_pair(z, coords, i, j, tol)) n = n + 1
      end do
   end function count_bonds

   pure function bonded_pair(z, coords, i, j, tol) result(is_bond)
      !! The same criterion `mqc_bond_perception` uses, on two atoms
      !!
      !! Restated here because that module's test takes a whole
      !! `system_geometry_t`; the tolerance is imported from it.
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

   pure function gamess_bonded(z1, xyz1, z2, xyz2) result(bonded)
      !! GAMESS's `PAIRBOND` at a scale of one
      !!
      !! `|r1 - r2| <= 1.2 r1 + 1.2 r2` in Angstrom, except that hydrogen's
      !! radius is not scaled and an untabulated element counts 1.6. Two atoms
      !! at one point are bonded, which is what makes an atom its own
      !! neighbour below.
      integer, intent(in) :: z1, z2
      real(dp), intent(in) :: xyz1(3), xyz2(3)   !! Bohr
      logical :: bonded

      bonded = norm2(xyz1 - xyz2)*BOHR_TO_ANGSTROM <= pair_radius(z1) + pair_radius(z2)
   end function gamess_bonded

   pure function pair_radius(z) result(r)
      integer, intent(in) :: z
      real(dp) :: r

      r = 1.6_dp
      if (z == 1) then
         r = GAMESS_RCOV(1)
      else if (z > 1 .and. z <= size(GAMESS_RCOV)) then
         r = 1.2_dp*GAMESS_RCOV(z)
      end if
   end function pair_radius

   subroutine build_bonded_model(z, coords, cut, model, error)
      !! The model system GAMESS builds for one cut bond
      !!
      !! Both ends of the bond and every atom bonded to either (`gamess_bonded`),
      !! then, for each atom left out that is bonded to one taken: a hydrogen
      !! comes in where it is, and anything else is replaced by a cap hydrogen
      !! at `GAMESS_CAP_LENGTH` from the atom it hangs off. That is GAMESS's
      !! construction at `RAFO(1)=1,1,1` (`fmolib.src`, the AFO model loop).
      !!
      !! One addition GAMESS does not make: a *terminal* heavy atom bonded to a
      !! taken one comes in whole, since a cap hydrogen closes one electron pair
      !! and a carbonyl oxygen is held by two. It never applies where GAMESS's
      !! own model would be closed-shell to begin with. Charged groups taken in
      !! are counted as `build_afo_model` counts them.
      !!
      !! Ends first: the bond-detached atom is `bda_local = 1` and the attached
      !! one `baa_local = 2`.
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      type(severed_bond_t), intent(in) :: cut
      type(afo_model_t), intent(out) :: model
      type(error_t), intent(inout) :: error

      logical, allocatable :: chosen(:)
      integer, allocatable :: order(:)
      integer :: n_atoms, i, j, n_real, n_caps, slot, a, b, n_nb

      if (error%has_error()) return
      n_atoms = size(z)
      a = cut%atom_a
      b = cut%atom_b
      if (a < 1 .or. a > n_atoms .or. b < 1 .or. b > n_atoms) then
         call error%set(ERROR_VALIDATION, "afo model: the cut bond names an atom the "// &
                        "system does not have")
         return
      end if

      allocate (chosen(n_atoms), source=.false.)
      chosen(a) = .true.
      chosen(b) = .true.
      do i = 1, n_atoms
         if (gamess_bonded(z(a), coords(:, a), z(i), coords(:, i)) .or. &
             gamess_bonded(z(b), coords(:, b), z(i), coords(:, i))) chosen(i) = .true.
      end do

      ! What hangs off the neighbours: hydrogens as they are, and terminal
      ! heavy atoms whole. Decided against the set above only, so it cannot
      ! cascade.
      do i = 1, n_atoms
         if (chosen(i)) cycle
         if (.not. hangs_off(i)) cycle
         if (z(i) == 1) then
            chosen(i) = .true.
            cycle
         end if
         n_nb = 0
         do j = 1, n_atoms
            if (j == i) cycle
            if (gamess_bonded(z(i), coords(:, i), z(j), coords(:, j))) n_nb = n_nb + 1
         end do
         if (n_nb == 1) chosen(i) = .true.
      end do

      n_real = count(chosen)
      allocate (order(n_real))
      order(1) = a
      order(2) = b
      slot = 2
      do i = 1, n_atoms
         if (.not. chosen(i) .or. i == a .or. i == b) cycle
         slot = slot + 1
         order(slot) = i
      end do

      n_caps = 0
      do slot = 1, n_real
         do j = 1, n_atoms
            if (chosen(j)) cycle
            if (z(order(slot)) == 1) cycle
            if (gamess_bonded(z(order(slot)), coords(:, order(slot)), z(j), coords(:, j))) then
               n_caps = n_caps + 1
            end if
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
      end do
      model%bda_local = 1
      model%baa_local = 2

      slot = n_real
      do i = 1, n_real
         if (z(order(i)) == 1) cycle
         do j = 1, n_atoms
            if (chosen(j)) cycle
            if (.not. gamess_bonded(z(order(i)), coords(:, order(i)), z(j), coords(:, j))) cycle
            slot = slot + 1
            model%z(slot) = 1
            model%sym(slot) = element_number_to_symbol(1)
            model%xyz(:, slot) = coords(:, order(i)) + cap_length(z(order(i))) &
                                 *(coords(:, j) - coords(:, order(i))) &
                                 /norm2(coords(:, j) - coords(:, order(i)))
         end do
      end do

      model%charge = model_formal_charge(z, coords, chosen, DEFAULT_BOND_TOLERANCE)
      model%nelec = sum(model%z) - model%charge
      if (mod(model%nelec, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "afo model: the model system for the bond "// &
                        "between atoms "//to_char(a)//" and "//to_char(b)//" has an "// &
                        "odd electron count, so the capping did not close every "// &
                        "valence it opened")
         return
      end if

   contains

      pure function hangs_off(k) result(hangs)
         integer, intent(in) :: k
         logical :: hangs
         integer :: m

         hangs = .false.
         do m = 1, n_atoms
            if (.not. chosen(m)) cycle
            if (.not. (m == a .or. m == b .or. &
                       gamess_bonded(z(a), coords(:, a), z(m), coords(:, m)) .or. &
                       gamess_bonded(z(b), coords(:, b), z(m), coords(:, m)))) cycle
            if (gamess_bonded(z(m), coords(:, m), z(k), coords(:, k))) then
               hangs = .true.
               return
            end if
         end do
      end function hangs_off

   end subroutine build_bonded_model

   pure function cap_length(z) result(r)
      !! Where a cap hydrogen goes from atom `z`, Bohr
      integer, intent(in) :: z
      real(dp) :: r

      r = 1.6_dp
      if (z >= 1 .and. z <= size(GAMESS_CAP_LENGTH)) then
         if (GAMESS_CAP_LENGTH(z) > 0.0_dp) r = GAMESS_CAP_LENGTH(z)
      end if
      r = to_bohr(r)
   end function cap_length

   pure function atom_lmo_count(z) result(n)
      !! How many occupied orbitals an atom has when fully bonded: C is 1s and
      !! four sp3, five. GAMESS's `LMOATOM` for the main groups.
      integer, intent(in) :: z
      integer :: n

      if (z <= 2) then
         n = 1
      else if (z <= 4) then
         n = 2
      else if (z <= 10) then
         n = 5
         if (z == 5) n = 4
      else if (z <= 12) then
         n = 6
      else if (z <= 18) then
         n = 9
      else if (z <= 20) then
         n = 10
      else if (z <= 30) then
         n = 15
      else
         n = 18
      end if
   end function atom_lmo_count

   subroutine bond_lmo_set(model, opts, set, n_on_bond, error)
      !! The frozen orbitals a cut bond contributes, from its model system
      !!
      !! Solve the model and localize every occupied orbital as
      !! `opts%localization` says. The detached atom's own are the
      !! `atom_lmo_count` with the largest population on it; of those, the one
      !! with the largest population on the attached atom is the bond's and
      !! goes first. Each is kept on every real atom of the model bonded to
      !! either end of the bond. All of it is GAMESS's selection (`fmolib.src`, after `LMOX`, with `CRITLOC` the
      !! on-atom population and the default `RAFO`).
      !!
      !! `n_on_bond` counts localized orbitals whose centroid sits within
      !! `BOND_ORBITAL_REACH` of the midpoint, as `bond_hybrid` does, for the
      !! caller to refuse a bond that is not single.
      type(afo_model_t), intent(in) :: model
      type(afo_options_t), intent(in) :: opts
      type(afo_lmo_set_t), intent(out) :: set
      integer, intent(out) :: n_on_bond
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(scf_numerics_t) :: scf_numerics
      real(dp), allocatable :: localized(:, :), centroids(:, :), s(:, :), pop_bda(:)
      integer, allocatable :: offsets(:), counts(:), keep(:), pick(:)
      logical, allocatable :: taken(:)
      real(dp) :: midpoint(3)
      real(dp) :: reach, best
      integer :: n_occ, i, k, n_real, n_keep, at, row, first, n, special

      if (error%has_error()) return
      n_on_bond = 0
      if (len_trim(opts%basis) == 0) then
         call error%set(ERROR_VALIDATION, "afo: no orbital basis was named for the "// &
                        "model system")
         return
      end if

      call build_czt_molecule(model%z, model%sym, model%xyz, trim(opts%basis), &
                              mol, error, force_cartesian=opts%cartesian)
      if (error%has_error()) return
      scf_numerics = opts%scf
      scf_numerics%max_iter = opts%scf_max_iter
      scf_numerics%energy_tol = opts%scf_energy_tol
      scf_numerics%density_tol = opts%scf_density_tol
      scf_numerics%grad_tol = opts%scf_grad_tol
      call run_czt_rhf(mol, model%nelec, opts%scf_max_iter, opts%scf_energy_tol, &
                       opts%scf_density_tol, opts%show_scf, scf, error, scf=scf_numerics, &
                       grad_tol=opts%scf_grad_tol)
      if (error%has_error()) return
      if (.not. scf%converged) then
         call error%set(ERROR_VALIDATION, "afo: the model system's SCF did not converge "// &
                        "in "//to_char(scf%iterations)//" iterations, orbital gradient "// &
                        to_char(scf%commutator)//" at the last, so there is no orbital "// &
                        "to freeze")
         return
      end if
      n_occ = scf%n_occupied
      call localize_model(mol, scf%orbitals, n_occ, opts, localized, centroids, error)
      if (error%has_error()) return

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)
      call mol%overlap(s)

      midpoint = 0.5_dp*(model%xyz(:, 1) + model%xyz(:, 2))
      reach = BOND_ORBITAL_REACH*norm2(model%xyz(:, 1) - model%xyz(:, 2))
      do i = 1, n_occ
         if (norm2(centroids(:, i) - midpoint) < reach) n_on_bond = n_on_bond + 1
      end do

      ! The detached atom's orbitals: the largest populations on it.
      allocate (pop_bda(n_occ))
      do i = 1, n_occ
         pop_bda(i) = atom_population(localized(:, i), 1)
      end do
      n = min(atom_lmo_count(model%z(1)), n_occ)
      allocate (pick(n), source=0)
      allocate (taken(n_occ), source=.false.)
      do k = 1, n
         best = -huge(1.0_dp)
         do i = 1, n_occ
            if (taken(i)) cycle
            if (pop_bda(i) > best) then
               best = pop_bda(i)
               pick(k) = i
            end if
         end do
         taken(pick(k)) = .true.
      end do
      ! Of those, the bond's: the largest population on the attached atom.
      special = 1
      do k = 2, n
         if (atom_population(localized(:, pick(k)), 2) > &
             atom_population(localized(:, pick(special)), 2)) special = k
      end do
      if (special /= 1) pick([1, special]) = pick([special, 1])

      ! Kept on every real atom bonded to either end.
      n_real = model%n_atoms - model%n_caps
      allocate (keep(n_real))
      n_keep = 0
      do at = 1, n_real
         if (gamess_bonded(model%z(1), model%xyz(:, 1), model%z(at), model%xyz(:, at)) .or. &
             gamess_bonded(model%z(2), model%xyz(:, 2), model%z(at), model%xyz(:, at))) then
            n_keep = n_keep + 1
            keep(n_keep) = at
         end if
      end do

      set%n_lmo = n
      set%n_at = n_keep
      allocate (set%atoms(n_keep), set%n_func(n_keep), set%pop(n_keep, n))
      do k = 1, n_keep
         set%atoms(k) = model%from_system(keep(k))
         set%n_func(k) = counts(keep(k))
      end do
      allocate (set%coeff(sum(set%n_func), n))
      do i = 1, n
         row = 0
         do k = 1, n_keep
            first = offsets(keep(k)) + 1
            set%coeff(row + 1:row + counts(keep(k)), i) = &
               localized(first:first + counts(keep(k)) - 1, pick(i))
            set%pop(k, i) = atom_population(localized(:, pick(i)), keep(k))
            row = row + counts(keep(k))
         end do
      end do
      call mol%destroy()

   contains

      pure function atom_population(c, atom) result(p)
         real(dp), intent(in) :: c(:)
         integer, intent(in) :: atom
         real(dp) :: p

         integer :: f, l

         f = offsets(atom) + 1
         l = offsets(atom) + counts(atom)
         p = dot_product(c(f:l), matmul(s(f:l, f:l), c(f:l)))
      end function atom_population

   end subroutine bond_lmo_set

   subroutine build_group_frozen_set(mol, atom_of, cut_bda, cut_baa, occupied, sets, &
                                     frozen, n_frozen_occ, error)
      !! A group's frozen orbitals, placed into its own basis
      !!
      !! For each boundary the group is cut across: the bond's orbital,
      !! occupied if the group holds the detached atom as a ghost and empty if
      !! it owns it; and, when it holds the ghost, the detached atom's other
      !! orbitals as well, empty -- GAMESS's `flmovec`, "1 occupied and 4
      !! virtual" for a carbon. Each orbital goes onto whichever of its atoms
      !! the group has, the rest dropped.
      !!
      !! Where a second boundary of the group runs through the same ghost atom
      !! its bond orbital is already frozen, so this one's copy of it is left
      !! out of the empties. GAMESS does not do this, and on propane cut at
      !! both bonds its SCF for the two end fragments together fails to
      !! converge; see the comment in the loop.
      !!
      !! Occupied columns first, as `build_frozen_basis` wants them.
      type(czt_molecule_t), intent(in) :: mol
      integer, intent(in) :: atom_of(:)       !! (mol%natm), the system atom of each slot
      integer, intent(in) :: cut_bda(:), cut_baa(:)   !! per boundary, system atoms
      logical, intent(in) :: occupied(:)      !! per boundary, the group holds the ghost
      type(afo_lmo_set_t), intent(in) :: sets(:)      !! per boundary
      real(dp), allocatable, intent(out) :: frozen(:, :)
      integer, intent(out) :: n_frozen_occ
      type(error_t), intent(inout) :: error

      integer, allocatable :: offsets(:), counts(:), slot_of(:)
      integer :: nb, i, j, k, col, n_col, pass, skip, at_baa, other
      real(dp) :: best

      if (error%has_error()) return
      nb = size(cut_bda)
      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)
      allocate (slot_of(max(maxval(atom_of), maxval(cut_bda), 1)), source=0)
      do i = 1, mol%natm
         slot_of(atom_of(i)) = i
      end do

      n_frozen_occ = count(occupied)
      n_col = 0
      do i = 1, nb
         n_col = n_col + 1
         if (occupied(i)) n_col = n_col + sets(i)%n_lmo - 1
      end do
      allocate (frozen(mol%nao, max(n_col, 1)), source=0.0_dp)

      col = 0
      ! Pass one: the occupied bond orbitals. Pass two: every empty one.
      do pass = 1, 2
         do i = 1, nb
            if (pass == 1) then
               if (.not. occupied(i)) cycle
               col = col + 1
               call place(sets(i), 1, frozen(:, col))
               if (error%has_error()) return
               cycle
            end if
            if (.not. occupied(i)) then
               col = col + 1
               call place(sets(i), 1, frozen(:, col))
               if (error%has_error()) return
               cycle
            end if
            ! Two boundaries on one ghost -- an atom detached from two
            ! neighbours, both held here -- each carry the atom's core and its
            ! other bonds, from two different model systems. One copy is
            ! frozen, the first boundary's; a second, slightly different copy
            ! would leave meaningless directions after orthogonalisation.
            if (any(occupied(:i - 1) .and. cut_bda(:i - 1) == cut_bda(i))) cycle
            do k = 2, sets(i)%n_lmo
               ! Another boundary of this group runs through the same ghost
               ! atom, so its bond orbital is frozen already -- occupied where
               ! the group holds the ghost for it too, empty where the group
               ! owns the atom at its other end. Freezing this set's copy of
               ! that bond as well would put a second, slightly different
               ! truncation of one orbital into the frozen space, and what
               ! survives orthogonalisation is then a meaningless direction.
               ! The copy is the orbital with the largest population on the
               ! bond's other atom.
               skip = 0
               do j = 1, nb
                  if (j == i) cycle
                  other = 0
                  if (cut_bda(j) == cut_bda(i) .and. occupied(j)) other = cut_baa(j)
                  if (cut_baa(j) == cut_bda(i) .and. .not. occupied(j)) other = cut_bda(j)
                  if (other == 0) cycle
                  at_baa = findloc(sets(i)%atoms, other, dim=1)
                  if (at_baa == 0) cycle
                  best = maxval(sets(i)%pop(at_baa, 2:sets(i)%n_lmo))
                  if (sets(i)%pop(at_baa, k) >= best) skip = 1
               end do
               if (skip == 1) then
                  n_col = n_col - 1
                  cycle
               end if
               col = col + 1
               call place(sets(i), k, frozen(:, col))
               if (error%has_error()) return
            end do
         end do
      end do
      if (col < size(frozen, 2)) frozen = frozen(:, :max(col, 1))

   contains

      subroutine place(set, k, column)
         type(afo_lmo_set_t), intent(in) :: set
         integer, intent(in) :: k
         real(dp), intent(out) :: column(:)

         integer :: a, row, slot

         column = 0.0_dp
         row = 0
         do a = 1, set%n_at
            slot = 0
            if (set%atoms(a) <= size(slot_of)) slot = slot_of(set%atoms(a))
            if (slot > 0) then
               if (counts(slot) /= set%n_func(a)) then
                  call error%set(ERROR_VALIDATION, "afo: a frozen orbital's block on "// &
                                 "atom "//to_char(set%atoms(a))//" has "// &
                                 to_char(set%n_func(a))//" coefficients where the "// &
                                 "group has "//to_char(counts(slot))//" functions, so "// &
                                 "it was built against a different basis")
                  return
               end if
               column(offsets(slot) + 1:offsets(slot) + counts(slot)) = &
                  set%coeff(row + 1:row + set%n_func(a), k)
            end if
            row = row + set%n_func(a)
         end do
         if (all(column == 0.0_dp)) then
            call error%set(ERROR_VALIDATION, "afo: a frozen orbital has none of its "// &
                           "atoms in the group it is frozen in")
         end if
      end subroutine place

   end subroutine build_group_frozen_set

   pure function lmo_set_pack_size(set) result(n)
      !! Reals `lmo_set_pack` writes for one set
      type(afo_lmo_set_t), intent(in) :: set
      integer :: n

      n = 2 + 2*set%n_at + size(set%coeff) + size(set%pop)
   end function lmo_set_pack_size

   pure subroutine lmo_set_pack(set, buf)
      !! Flatten a set into reals, for sharing across ranks; counts are exact
      !! integers in a double
      type(afo_lmo_set_t), intent(in) :: set
      real(dp), intent(out) :: buf(:)

      integer :: at

      buf(1) = real(set%n_lmo, dp)
      buf(2) = real(set%n_at, dp)
      at = 2
      buf(at + 1:at + set%n_at) = real(set%atoms, dp)
      at = at + set%n_at
      buf(at + 1:at + set%n_at) = real(set%n_func, dp)
      at = at + set%n_at
      buf(at + 1:at + size(set%coeff)) = reshape(set%coeff, [size(set%coeff)])
      at = at + size(set%coeff)
      buf(at + 1:at + size(set%pop)) = reshape(set%pop, [size(set%pop)])
   end subroutine lmo_set_pack

   subroutine lmo_set_unpack(buf, set)
      !! The inverse of `lmo_set_pack`
      real(dp), intent(in) :: buf(:)
      type(afo_lmo_set_t), intent(out) :: set

      integer :: at, nf

      set%n_lmo = nint(buf(1))
      set%n_at = nint(buf(2))
      at = 2
      set%atoms = nint(buf(at + 1:at + set%n_at))
      at = at + set%n_at
      set%n_func = nint(buf(at + 1:at + set%n_at))
      at = at + set%n_at
      nf = sum(set%n_func)
      set%coeff = reshape(buf(at + 1:at + nf*set%n_lmo), [nf, set%n_lmo])
      at = at + nf*set%n_lmo
      set%pop = reshape(buf(at + 1:at + set%n_at*set%n_lmo), [set%n_at, set%n_lmo])
   end subroutine lmo_set_unpack

   subroutine orient_cut(z, coords, cut, error, detached)
      !! Decide which end of a cut bond is the bond-detached atom
      !!
      !! `atom_a` comes out as the detached end and `frag_a` as its owner. A
      !! deck that names detached atoms is followed: exactly one end of the
      !! bond has to be named, and naming both is refused. Otherwise the sp3
      !! end is detached -- the one with four neighbours when the other has
      !! fewer, which on a protein backbone cut at C-alpha--C(=O) is the
      !! C-alpha, the choice GAMESS's documentation and every FMO protein
      !! study make -- and where both or neither are, the lower-numbered one,
      !! which is how `find_severed_bonds` hands the bond over.
      integer, intent(in) :: z(:)
      real(dp), intent(in) :: coords(:, :)
      type(severed_bond_t), intent(inout) :: cut
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: detached(:)
         !! System atoms, 1-based, a deck named as detached ends

      logical :: named_a, named_b, flip
      integer :: n_a, n_b, swap

      if (error%has_error()) return
      named_a = .false.
      named_b = .false.
      if (present(detached)) then
         named_a = any(detached == cut%atom_a)
         named_b = any(detached == cut%atom_b)
      end if
      if (named_a .and. named_b) then
         call error%set(ERROR_VALIDATION, "afo: both atoms of the cut bond between "// &
                        to_char(cut%atom_a)//" and "//to_char(cut%atom_b)//" are named "// &
                        "as detached; a bond has one detached end")
         return
      end if

      if (named_a .or. named_b) then
         flip = named_b
      else
         n_a = neighbour_count(cut%atom_a)
         n_b = neighbour_count(cut%atom_b)
         flip = n_b == 4 .and. n_a /= 4
      end if
      if (flip) then
         swap = cut%atom_a
         cut%atom_a = cut%atom_b
         cut%atom_b = swap
         swap = cut%frag_a
         cut%frag_a = cut%frag_b
         cut%frag_b = swap
      end if

   contains

      pure function neighbour_count(i) result(n)
         integer, intent(in) :: i
         integer :: n
         integer :: j

         n = 0
         do j = 1, size(z)
            if (j == i) cycle
            if (gamess_bonded(z(i), coords(:, i), z(j), coords(:, j))) n = n + 1
         end do
      end function neighbour_count

   end subroutine orient_cut

end module mqc_czt_afo
