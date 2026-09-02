!! Semantic validation of a system and its fragment terms
module mqc_validate
   !! Checks that a system and a term list describe a calculation that means
   !! something, before anything spends a node-hour finding out it did not.
   !!
   !! This is not schema validation. `mqc_json_schema` checks that a *document*
   !! is well formed -- keys spelled right, required ones present. These checks
   !! are about the thing the document describes, and apply equally to a system
   !! assembled through the C interface, which never sees a document at all.
   !!
   !! **None of these failures crashes.** A term list that is not closed under
   !! subsets, fragment charges that do not sum, an atom in no monomer -- each
   !! yields a converged number that is wrong.
   !!
   !! `strict` is the default and refuses. The permissive mode warns instead,
   !! for the cases that are conventions rather than laws -- a partition that
   !! deliberately omits atoms, a hand-built term list with a scheme of its
   !! own. It does not stop checking.
   use pic_types, only: dp, default_int, int64
   use mqc_string_utils, only: int_to_text
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_fraglist, only: fraglist_t
   use mqc_bond_perception, only: missing_broken_bonds
   use mqc_error, only: error_t, ERROR_VALIDATION
   use pic_logger, only: logger => global_logger
   implicit none
   private

   public :: validate_system   !! Check a system describes a fragmentable molecule
   public :: validate_terms    !! Check a term list can actually be assembled

contains

   subroutine validate_system(sys_geom, strict, error, check_bonds)
      !! Check a system before it is fragmented
      !!
      !! `check_bonds` runs the perception audit, which is `O(N^2)` in the atom
      !! count and so is asked for rather than assumed.
      ! TODO(mqc): reads `fragment_sizes` and `fragment_atoms` without checking
      ! that either is allocated, where the rest of the code branches on
      ! `allocated(sys_geom%fragment_atoms)`. A fixed-size-monomer geometry --
      ! what `initialize_system_geometry` builds -- reaches here from the
      ! driver and dereferences both.
      type(system_geometry_t), intent(in) :: sys_geom
      logical, intent(in) :: strict
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: check_bonds

      integer, allocatable :: owner_count(:)
      integer, allocatable :: missing_i(:), missing_j(:)
      integer :: iatom, imon, islot, atom, total_charge, n_missing
      logical :: audit

      if (sys_geom%total_atoms <= 0) then
         call complain(strict, error, "the system has no atoms")
         return
      end if
      if (sys_geom%multiplicity < 1) then
         call complain(strict, error, "multiplicity must be at least 1")
         if (strict) return
      end if
      if (sys_geom%n_monomers <= 0) then
         call complain(strict, error, "the system has no monomers to fragment")
         return
      end if

      ! ---- the partition ----------------------------------------------------
      allocate (owner_count(sys_geom%total_atoms))
      owner_count = 0
      do imon = 1, sys_geom%n_monomers
         do islot = 1, sys_geom%fragment_sizes(imon)
            atom = sys_geom%fragment_atoms(islot, imon)
            if (atom < 0 .or. atom >= sys_geom%total_atoms) then
               call complain(strict, error, "monomer "//int_to_text(imon)// &
                             " names atom "//int_to_text(atom)//", which does not exist "// &
                             "(atom indices are 0-based)")
               if (strict) return
               cycle
            end if
            owner_count(atom + 1) = owner_count(atom + 1) + 1
         end do
      end do

      ! An atom in no monomer is in no fragment either, so it contributes to
      ! nothing and its electrons are simply absent from the expansion.
      do iatom = 1, sys_geom%total_atoms
         if (owner_count(iatom) == 0) then
            call complain(strict, error, "atom "//int_to_text(iatom - 1)// &
                          " belongs to no monomer, so no fragment will ever contain it")
            if (strict) return
            exit
         end if
      end do

      ! ---- charges ----------------------------------------------------------
      ! Overlapping partitions legitimately count an atom more than once, so
      ! this only holds where each atom has exactly one owner.
      if (allocated(sys_geom%fragment_charges) .and. all(owner_count == 1)) then
         total_charge = sum(sys_geom%fragment_charges(1:sys_geom%n_monomers))
         if (total_charge /= sys_geom%charge) then
            call complain(strict, error, "monomer charges sum to "//int_to_text(total_charge)// &
                          " but the system charge is "//int_to_text(sys_geom%charge)// &
                          "; the fragments describe a different number of electrons "// &
                          "than the molecule")
            if (strict) return
         end if
      end if

      ! ---- bonds ------------------------------------------------------------
      audit = .false.
      if (present(check_bonds)) audit = check_bonds
      if (audit .and. allocated(sys_geom%bonds)) then
         call missing_broken_bonds(sys_geom, sys_geom%bonds, size(sys_geom%bonds), &
                                   missing_i, missing_j, n_missing)
         if (n_missing > 0) then
            call complain(strict, error, "the geometry implies "//int_to_text(n_missing)// &
                          " bond(s) crossing monomer boundaries that were never "// &
                          "declared, starting with atoms "//int_to_text(missing_i(1))//" and "// &
                          int_to_text(missing_j(1))//"; those fragments would have uncapped "// &
                          "valences")
            if (strict) return
         end if
      end if

      deallocate (owner_count)
   end subroutine validate_system

   subroutine validate_terms(terms, sys_geom, strict, error)
      !! Check a term list can be assembled into an expansion
      !!
      !! The subset check is the one worth the cost; everything else here is a
      !! typo detector.
      type(fraglist_t), intent(in) :: terms
      type(system_geometry_t), intent(in) :: sys_geom
      logical, intent(in) :: strict
      type(error_t), intent(inout) :: error

      integer(int64) :: iterm
      integer(default_int) :: level, islot, jslot

      if (terms%n_terms <= 0) then
         call complain(strict, error, "the term list is empty; there is nothing to compute")
         return
      end if

      do iterm = 1, terms%n_terms
         level = terms%level_of(iterm)
         if (level <= 0) then
            call complain(strict, error, "term "//text64(iterm)//" names no monomers")
            if (strict) return
            cycle
         end if

         do islot = 1, level
            ! Monomer indices are 1-based here, unlike the atom indices above.
            if (terms%terms(iterm, islot) < 1 .or. &
                terms%terms(iterm, islot) > sys_geom%n_monomers) then
               call complain(strict, error, "term "//text64(iterm)//" names monomer "// &
                             int_to_text(int(terms%terms(iterm, islot)))//", but the system has "// &
                             int_to_text(sys_geom%n_monomers)//" (monomer indices are 1-based)")
               if (strict) return
               exit
            end if
            do jslot = islot + 1, level
               if (terms%terms(iterm, islot) == terms%terms(iterm, jslot)) then
                  call complain(strict, error, "term "//text64(iterm)// &
                                " names monomer "//int_to_text(int(terms%terms(iterm, islot)))// &
                                " twice")
                  if (strict) return
                  exit
               end if
            end do
         end do
      end do

      call check_subset_closure(terms, strict, error)
   end subroutine validate_terms

   subroutine check_subset_closure(terms, strict, error)
      !! Every proper subset of every term must itself be a term
      !!
      !! An n-body contribution is the term's energy less the delta of each
      !! proper subset, so a missing one is not an approximation -- the
      !! assembly looks it up and fails. `fraglist_t%close_subsets` is the fix.
      ! TODO(mqc): `contains_term` scans the whole list for each subset of each
      ! term, so this is quadratic in the term count on a check every run pays.
      ! It also compares slots in order, so a supplied term whose monomers are
      ! not in ascending order is reported missing when it is present.
      type(fraglist_t), intent(in) :: terms
      logical, intent(in) :: strict
      type(error_t), intent(inout) :: error

      integer(int64) :: iterm
      integer(default_int) :: level, missing_slot

      do iterm = 1, terms%n_terms
         level = terms%level_of(iterm)
         if (level < 2) cycle

         ! Only the (n-1)-subsets are checked: if each is present, its own
         ! subsets follow from the same check applied to it, so the closure
         ! holds without enumerating 2^n sets per term.
         missing_slot = first_missing_subset(terms, iterm, level)
         if (missing_slot > 0) then
            call complain(strict, error, "term "//text64(iterm)// &
                          " is present but the subset without its monomer "// &
                          int_to_text(int(terms%terms(iterm, missing_slot)))//" is not; "// &
                          "the expansion cannot be assembled from this list "// &
                          "(see close_subsets)")
            if (strict) return
            return
         end if
      end do
   end subroutine check_subset_closure

   pure function first_missing_subset(terms, iterm, level) result(slot)
      !! Which dropped monomer leaves a subset that is not in the list, or 0
      type(fraglist_t), intent(in) :: terms
      integer(int64), intent(in) :: iterm
      integer(default_int), intent(in) :: level
      integer(default_int) :: slot

      integer(default_int) :: candidate(level - 1)
      integer(default_int) :: drop, islot, kept

      do drop = 1, level
         kept = 0
         do islot = 1, level
            if (islot == drop) cycle
            kept = kept + 1
            candidate(kept) = terms%terms(iterm, islot)
         end do
         if (.not. contains_term(terms, candidate, level - 1)) then
            slot = drop
            return
         end if
      end do
      slot = 0
   end function first_missing_subset

   pure function contains_term(terms, wanted, level) result(found)
      !! Whether a term with exactly these monomers is in the list
      type(fraglist_t), intent(in) :: terms
      integer(default_int), intent(in) :: wanted(:)
      integer(default_int), intent(in) :: level
      logical :: found

      integer(int64) :: iterm

      found = .false.
      do iterm = 1, terms%n_terms
         if (terms%level_of(iterm) /= level) cycle
         if (all(terms%terms(iterm, 1:level) == wanted(1:level))) then
            found = .true.
            return
         end if
      end do
   end function contains_term

   subroutine complain(strict, error, message)
      !! Record a problem as a refusal or a warning
      !!
      !! The first refusal wins and later ones are not overwritten, so the
      !! message a caller sees is the first thing that was wrong rather than
      !! the last thing checked.
      logical, intent(in) :: strict
      type(error_t), intent(inout) :: error
      character(len=*), intent(in) :: message

      if (strict) then
         if (.not. error%has_error()) call error%set(ERROR_VALIDATION, message)
      else
         call logger%warning("unchecked input: "//message)
      end if
   end subroutine complain

   pure function text64(value) result(out)
      integer(int64), intent(in) :: value
      character(len=:), allocatable :: out
      character(len=24) :: buffer
      write (buffer, "(I0)") value
      out = trim(adjustl(buffer))
   end function text64

end module mqc_validate
