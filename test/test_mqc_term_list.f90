module test_mqc_term_list
   !! `generate_mbe_term_list`, which a frozen optimization rests on
   !!
   !! The driver and `mqc_geometry_optimizer` both call this, and a geometry
   !! optimization freezes whatever it returns and then uses it for every step
   !! of the run. Two properties therefore have to hold, and neither is
   !! self-evident from reading the code:
   !!
   !!   * the same geometry gives the same list, or the frozen list is not the
   !!     list the run would otherwise have used
   !!   * a moved geometry can give a *different* list, which is the whole
   !!     reason freezing exists
   !!
   !! No DL-FIND needed: this is the fragment layer, so it runs in CI whatever
   !! the optimizer backend is set to.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_frag_utils, only: generate_mbe_term_list, get_nfrags, binomial, fragment_lookup_t
   use mqc_combinatorics, only: is_auxiliary_row, real_count_of, rows_match_counterpoise, &
                                COUNTERPOISE_VMFC, COUNTERPOISE_SSFC
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_config_adapter, only: driver_config_t
   use pic_types, only: dp, int64
   implicit none
   private
   public :: collect_mqc_term_list_tests

contains

   subroutine collect_mqc_term_list_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("all_terms_without_cutoffs", test_all_terms_without_cutoffs), &
                  new_unittest("deterministic", test_deterministic), &
                  new_unittest("screening_removes_terms", test_screening_removes_terms), &
                  new_unittest("screened_list_keeps_monomers", test_screened_keeps_monomers), &
                  new_unittest("list_moves_with_geometry", test_list_moves_with_geometry), &
                  new_unittest("counterpoise_gives_every_n_mer_its_subsets", test_vmfc_rows), &
                  new_unittest("counterpoise_follows_screening", test_vmfc_after_screening), &
                  new_unittest("full_basis_counterpoise_rows_ghost_the_whole_system", test_ssfc_rows), &
                  new_unittest("full_basis_counterpoise_follows_screening", test_ssfc_after_screening), &
                  new_unittest("full_basis_counterpoise_ghosts_the_system_under_a_reference", &
                               test_ssfc_with_reference), &
                  new_unittest("full_basis_counterpoise_of_two_monomers_is_vmfc", test_ssfc_two_monomers), &
                  new_unittest("full_basis_counterpoise_is_ordered_by_real_monomers", test_ssfc_order), &
                  new_unittest("reference_closure_has_the_derived_length", test_reference_count), &
                  new_unittest("reference_closure_is_closed_complete_minimal", test_reference_closure), &
                  new_unittest("reference_closure_follows_screening", test_reference_after_screening) &
                  ]
   end subroutine collect_mqc_term_list_tests

   subroutine make_chain(sys_geom, spacing)
      !! Four one-atom monomers in a line, `spacing` Bohr apart
      !!
      !! One atom per monomer keeps the inter-monomer distance equal to the
      !! atom separation, so a cutoff in the test means exactly what it says.
      type(system_geometry_t), intent(out) :: sys_geom
      real(dp), intent(in) :: spacing

      integer :: i

      sys_geom%n_monomers = 4
      sys_geom%atoms_per_monomer = 1
      sys_geom%total_atoms = 4
      sys_geom%charge = 0
      sys_geom%multiplicity = 1
      allocate (sys_geom%element_numbers(4))
      allocate (sys_geom%coordinates(3, 4))
      sys_geom%element_numbers = 10  ! neon: closed shell, never bonded
      sys_geom%coordinates = 0.0_dp
      do i = 1, 4
         sys_geom%coordinates(1, i) = real(i - 1, dp)*spacing
      end do
   end subroutine make_chain

   subroutine test_all_terms_without_cutoffs(error)
      !! With no cutoffs the list is the full combinatorial one
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_terms

      call make_chain(sys_geom, 4.0_dp)
      config%nlevel = 2

      call generate_mbe_term_list(sys_geom, config, 2, polymers, n_terms)

      ! 4 monomers + 6 pairs
      call check(error, n_terms == 10_int64, &
                 "MBE(2) over 4 monomers is 4 monomers and 6 dimers")
      if (allocated(error)) return

      call check(error, n_terms == get_nfrags(4, 2), &
                 "the count should match get_nfrags, which sizes the array")
   end subroutine test_all_terms_without_cutoffs

   subroutine test_deterministic(error)
      !! The same geometry gives the same list, term for term
      !!
      !! This is what makes freezing sound: the list captured once at the start
      !! is the list every later step would have generated for that geometry.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: first(:, :), second(:, :)
      integer(int64) :: n_first, n_second

      call make_chain(sys_geom, 4.0_dp)
      config%nlevel = 2
      allocate (config%fragment_cutoffs(2))
      config%fragment_cutoffs = [0.0_dp, 5.0_dp]

      call generate_mbe_term_list(sys_geom, config, 2, first, n_first)
      call generate_mbe_term_list(sys_geom, config, 2, second, n_second)

      call check(error, n_first == n_second, "the term count should not change between calls")
      if (allocated(error)) return

      call check(error, all(first(1:n_first, :) == second(1:n_second, :)), &
                 "the terms should be identical, in the same order")
   end subroutine test_deterministic

   subroutine test_screening_removes_terms(error)
      !! A cutoff shorter than the chain drops the distant pairs
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_terms

      ! Neighbours 4 Bohr apart, so the pairs are at 4, 8 and 12 Bohr. The
      ! cutoff is in Angstrom, and 4 Bohr is about 2.117 Angstrom, so a cutoff
      ! of 3 Angstrom keeps only the three adjacent pairs.
      call make_chain(sys_geom, 4.0_dp)
      config%nlevel = 2
      allocate (config%fragment_cutoffs(2))
      config%fragment_cutoffs = [0.0_dp, 3.0_dp]

      call generate_mbe_term_list(sys_geom, config, 2, polymers, n_terms)

      call check(error, n_terms == 7_int64, &
                 "4 monomers and the 3 adjacent dimers should survive a 3 Angstrom cutoff")
   end subroutine test_screening_removes_terms

   subroutine test_screened_keeps_monomers(error)
      !! Screening never removes a monomer
      !!
      !! Monomers carry the leading term of the expansion, so a screen that
      !! dropped one would not be an approximation but a different system.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_terms, iterm
      integer :: n_monomer_terms

      call make_chain(sys_geom, 20.0_dp)  ! far enough that every pair is screened
      config%nlevel = 2
      allocate (config%fragment_cutoffs(2))
      config%fragment_cutoffs = [0.0_dp, 2.0_dp]

      call generate_mbe_term_list(sys_geom, config, 2, polymers, n_terms)

      n_monomer_terms = 0
      do iterm = 1, n_terms
         if (count(polymers(iterm, :) > 0) == 1) n_monomer_terms = n_monomer_terms + 1
      end do

      call check(error, n_monomer_terms == 4, &
                 "every monomer should survive however tight the cutoff")
      if (allocated(error)) return

      call check(error, n_terms == 4_int64, &
                 "and with every pair beyond the cutoff, nothing else should")
   end subroutine test_screened_keeps_monomers

   subroutine test_list_moves_with_geometry(error)
      !! The same system at a different geometry can give a different list
      !!
      !! The reason `keywords.optimization.freeze_terms` exists. Move the
      !! monomers closer and a pair that was screened out comes back, which
      !! mid-optimization is a step in the energy the optimizer reads as real.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: spread_out, close_up
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_spread, n_close

      config%nlevel = 2
      allocate (config%fragment_cutoffs(2))
      config%fragment_cutoffs = [0.0_dp, 3.0_dp]

      call make_chain(spread_out, 4.0_dp)
      call generate_mbe_term_list(spread_out, config, 2, polymers, n_spread)

      call make_chain(close_up, 2.0_dp)
      call generate_mbe_term_list(close_up, config, 2, polymers, n_close)

      call check(error, n_close > n_spread, &
                 "bringing the monomers together should bring screened pairs back")
   end subroutine test_list_moves_with_geometry

   subroutine test_vmfc_rows(error)
      !! Counterpoise at level 3, where the recursion first has depth
      !!
      !! Level 2 is the easy case: a pair has two subsets and both are
      !! monomers. Level 3 is where the rule has to be a rule -- a trimer
      !! contributes six subsets, three of them pairs, and each of those pairs
      !! must be ghosted against the *trimer* rather than against itself. Get
      !! that wrong and the pair rows collide with the level-2 ones, which is a
      !! wrong answer and not a crash.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_terms, i
      integer :: n_aux, n_real

      call make_chain(sys_geom, 4.0_dp)
      config%nlevel = 3
      config%counterpoise = "vmfc"

      call generate_mbe_term_list(sys_geom, config, 3, polymers, n_terms)

      ! 4 monomers, 6 pairs and 4 trimers is 14 ordinary rows. Each pair adds
      ! its 2 subsets and each trimer its 6, so 14 + 12 + 24.
      call check(error, n_terms == 50_int64, &
                 "MBE(3) over 4 monomers under counterpoise is 50 rows, not "// &
                 int_str(n_terms))
      if (allocated(error)) return

      n_aux = 0
      n_real = 0
      do i = 1, n_terms
         if (is_auxiliary_row(polymers(i, :), COUNTERPOISE_VMFC)) then
            n_aux = n_aux + 1
         else
            n_real = n_real + 1
         end if
      end do

      call check(error, n_real == 14, "the ordinary expansion must survive intact")
      if (allocated(error)) return
      call check(error, n_aux == 36, "every n-mer owes 2**n - 2 ghosted subsets")
      if (allocated(error)) return

      ! Every auxiliary row must belong to a parent that is actually being
      ! computed, or it is subtracted from nothing.
      do i = 1, n_terms
         if (.not. is_auxiliary_row(polymers(i, :), COUNTERPOISE_VMFC)) cycle
         call check(error, has_parent(polymers, n_terms, abs(polymers(i, :))), &
                    "a ghosted row has no parent n-mer in the list")
         if (allocated(error)) return
      end do

      ! And no auxiliary row may be all-ghost: something has to be real.
      do i = 1, n_terms
         if (.not. is_auxiliary_row(polymers(i, :), COUNTERPOISE_VMFC)) cycle
         call check(error, real_count_of(polymers(i, :)) >= 1, &
                    "a ghosted row with no real monomer is not a term")
         if (allocated(error)) return
      end do
   end subroutine test_vmfc_rows

   subroutine test_vmfc_after_screening(error)
      !! A screened-out pair does not bring ghosted monomers along with it
      !!
      !! The rows are added after screening for exactly this reason. Added
      !! before, a distant pair would be dropped and its two ghosted monomers
      !! would stay -- auxiliary rows subtracted by nothing, paid for in full.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_terms, i
      integer :: n_pairs, n_aux

      call make_chain(sys_geom, 4.0_dp)
      config%nlevel = 2
      config%counterpoise = "vmfc"
      allocate (config%fragment_cutoffs(2))
      ! Cutoffs are in Angstrom and the chain is 4 Bohr (~2.117 Angstrom) per
      ! step, so 3 Angstrom keeps the adjacent pairs and drops the rest.
      config%fragment_cutoffs = [0.0_dp, 3.0_dp]

      call generate_mbe_term_list(sys_geom, config, 2, polymers, n_terms)

      n_pairs = 0
      n_aux = 0
      do i = 1, n_terms
         if (is_auxiliary_row(polymers(i, :), COUNTERPOISE_VMFC)) then
            n_aux = n_aux + 1
         else if (real_count_of(polymers(i, :)) == 2) then
            n_pairs = n_pairs + 1
         end if
      end do

      call check(error, n_pairs == 3, &
                 "a 3 Angstrom cutoff on a 4 Bohr chain leaves the 3 adjacent pairs")
      if (allocated(error)) return
      call check(error, n_aux == 2*n_pairs, &
                 "each surviving pair owes exactly two ghosted monomers, and a "// &
                 "screened one owes none")
   end subroutine test_vmfc_after_screening

   subroutine test_ssfc_rows(error)
      !! The full-basis list: one own-basis and one full-basis row per monomer,
      !! one full-basis row per term of two or more
      !!
      !! Levels 2 to 4 over four monomers, so the whole system is a term at the
      !! last, where its complement is empty and it stays all-real. Level 1 adds
      !! nothing under either scheme.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_terms
      integer :: level

      call make_chain(sys_geom, 4.0_dp)
      config%counterpoise = "ssfc"

      config%nlevel = 1
      call generate_mbe_term_list(sys_geom, config, 1, polymers, n_terms)
      call check(error, n_terms == 4_int64, "level 1 adds no rows under the full-basis scheme")
      if (allocated(error)) return
      call check(error,.not. any(polymers(1:n_terms, :) < 0), "level 1 ghosts nothing")
      if (allocated(error)) return

      do level = 2, 4
         config%nlevel = level
         call generate_mbe_term_list(sys_geom, config, level, polymers, n_terms)

         ! 4 own-basis monomers, 4 full-basis ones, and a row for every term
         ! of 2 to `level` monomers: 6, then 4 more, then 1 more.
         call check(error, n_terms == 8_int64 + get_nfrags(4, level) - 4_int64, &
                    "SSFC("//int_str(int(level, int64))//") over 4 monomers is 8 monomer rows "// &
                    "and one per term of order 2 and above, not "//int_str(n_terms))
         if (allocated(error)) return

         call check_full_basis_list(error, polymers, n_terms, 4, level)
         if (allocated(error)) return
         call check(error,.not. rows_match_counterpoise(polymers, n_terms, 4, COUNTERPOISE_VMFC), &
                    "a full-basis list of four monomers is not a VMFC list")
         if (allocated(error)) return
      end do
   end subroutine test_ssfc_rows

   subroutine test_ssfc_after_screening(error)
      !! Screening decides the terms, and the full-basis rows follow it
      !!
      !! A 3 Angstrom cutoff on the 4 Bohr chain keeps the three adjacent pairs;
      !! a 1 Angstrom one keeps none, leaving the monomers alone -- two rows
      !! each and no summed term. A kept pair still ghosts the monomers on both
      !! sides of it, kept terms or not. In a chain with one monomer out of
      !! reach of the rest, that monomer keeps both its rows and the kept pairs
      !! are the only summed terms.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_terms

      call make_chain(sys_geom, 4.0_dp)
      config%nlevel = 2
      config%counterpoise = "ssfc"
      allocate (config%fragment_cutoffs(2))

      config%fragment_cutoffs = [0.0_dp, 3.0_dp]
      call generate_mbe_term_list(sys_geom, config, 2, polymers, n_terms)
      call check(error, n_terms == 11_int64, &
                 "4 monomers twice and the 3 adjacent pairs once is 11 rows, not "//int_str(n_terms))
      if (allocated(error)) return
      call check_full_basis_list(error, polymers, n_terms, 4, 2)
      if (allocated(error)) return

      config%fragment_cutoffs = [0.0_dp, 1.0_dp]
      call generate_mbe_term_list(sys_geom, config, 2, polymers, n_terms)
      call check(error, n_terms == 8_int64, &
                 "with every pair screened only the monomers remain, twice each: "//int_str(n_terms))
      if (allocated(error)) return
      call check_full_basis_list(error, polymers, n_terms, 4, 2)
      if (allocated(error)) return

      ! Mixed: the fourth monomer is moved far from the others, so a 3 Angstrom
      ! cutoff keeps the pairs 1-2 and 2-3 and leaves monomer 4 with no pair.
      ! It still has both its own-basis and its full-basis row.
      sys_geom%coordinates(1, 4) = 40.0_dp
      config%fragment_cutoffs = [0.0_dp, 3.0_dp]
      call generate_mbe_term_list(sys_geom, config, 2, polymers, n_terms)
      call check(error, n_terms == 10_int64, &
                 "4 monomers twice and the 2 kept pairs once is 10 rows, not "//int_str(n_terms))
      if (allocated(error)) return
      call check_full_basis_list(error, polymers, n_terms, 4, 2)
      if (allocated(error)) return
      call check(error, count(polymers(1:n_terms, 1) == 4 .and. all(polymers(1:n_terms, 2:) == 0, dim=2)) == 1, &
                 "the isolated monomer keeps its own-basis row")
      if (allocated(error)) return
      call check(error, count(polymers(1:n_terms, 1) == 4 .and. polymers(1:n_terms, 2) < 0) == 1, &
                 "the isolated monomer keeps its full-basis row")
   end subroutine test_ssfc_after_screening

   subroutine test_ssfc_with_reference(error)
      !! Reducing to one fragment's interactions keeps the whole system's basis
      !!
      !! The reduced list drops far monomers and the terms only they need, but
      !! the basis is still the cluster's: every summed row ghosts the monomers
      !! the list no longer has as well as the ones it does. Full-basis monomer
      !! rows are built for the monomers kept, and the count is the reduced
      !! ordinary list's plus one row for each of those. `n_full` is what the
      !! same deck gives with no reference.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config, plain
      integer, allocatable :: polymers(:, :), reduced(:, :), everything(:, :)
      integer(int64) :: n_terms, n_reduced, n_all, n_full
      integer :: ref, level

      call make_line(sys_geom, 6, 4.0_dp)

      ! Unscreened at levels 2 and 3, then screened so that some monomers are
      ! dropped by the reduction.
      do level = 2, 3
         do ref = 1, 6
            config%nlevel = level
            config%reference_fragment = ref
            config%counterpoise = "ssfc"
            call generate_mbe_term_list(sys_geom, config, level, polymers, n_terms, n_full=n_full)

            plain = config
            plain%counterpoise = "none"
            call generate_mbe_term_list(sys_geom, plain, level, reduced, n_reduced)
            call check(error, n_terms == n_reduced + count_rows_of_size(reduced, n_reduced, 1), &
                       "the reduced full-basis list is the reduced ordinary list and one row "// &
                       "for each monomer")
            if (allocated(error)) return

            plain%reference_fragment = 0
            plain%counterpoise = "ssfc"
            call generate_mbe_term_list(sys_geom, plain, level, everything, n_all)
            call check(error, n_full == n_all, "n_full should be the unreduced full-basis count")
            if (allocated(error)) return

            call check_full_basis_list(error, polymers, n_terms, 6, level)
            if (allocated(error)) return
         end do
      end do

      allocate (config%fragment_cutoffs(3))
      config%fragment_cutoffs = [0.0_dp, 5.0_dp, 3.0_dp]
      do ref = 1, 6
         config%nlevel = 3
         config%reference_fragment = ref
         call generate_mbe_term_list(sys_geom, config, 3, polymers, n_terms)
         call check_full_basis_list(error, polymers, n_terms, 6, 3)
         if (allocated(error)) return

         plain = config
         plain%counterpoise = "none"
         call generate_mbe_term_list(sys_geom, plain, 3, reduced, n_reduced)
         call check(error, n_terms == n_reduced + count_rows_of_size(reduced, n_reduced, 1), &
                    "screened and reduced: the ordinary list and one row for each kept monomer")
         if (allocated(error)) return
      end do
   end subroutine test_ssfc_with_reference

   subroutine test_ssfc_two_monomers(error)
      !! For two monomers the two schemes are one list
      !!
      !! `[1]`, `[2]`, `[1,-2]`, `[-1,2]` and `[1,2]`: the pair's complement is
      !! empty, so it is all-real, and each monomer in the basis of both is the
      !! VMFC subset. Asserted so that nobody "fixes" it, and so that either
      !! scheme's check accepts the other's list.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: vmfc, ssfc
      integer, allocatable :: from_vmfc(:, :), from_ssfc(:, :)
      integer(int64) :: n_vmfc, n_ssfc, i, j
      integer :: a(2), b(2)
      logical :: found

      call make_line(sys_geom, 2, 4.0_dp)
      vmfc%nlevel = 2
      vmfc%counterpoise = "vmfc"
      ssfc%nlevel = 2
      ssfc%counterpoise = "ssfc"
      call generate_mbe_term_list(sys_geom, vmfc, 2, from_vmfc, n_vmfc)
      call generate_mbe_term_list(sys_geom, ssfc, 2, from_ssfc, n_ssfc)

      call check(error, n_ssfc == 5_int64 .and. n_vmfc == 5_int64, &
                 "two monomers give 5 rows under either scheme")
      if (allocated(error)) return

      ! As sets: the order within a list is the scheme's own.
      do i = 1, n_ssfc
         found = .false.
         a = from_ssfc(i, 1:2)
         call ascending(a)
         do j = 1, n_vmfc
            b = from_vmfc(j, 1:2)
            call ascending(b)
            if (all(a == b)) found = .true.
         end do
         call check(error, found, "an SSFC row of two monomers is missing from the VMFC list")
         if (allocated(error)) return
      end do

      call check(error, rows_match_counterpoise(from_ssfc, n_ssfc, 2, COUNTERPOISE_VMFC), &
                 "the VMFC check should accept the SSFC list of two monomers")
      if (allocated(error)) return
      call check(error, rows_match_counterpoise(from_vmfc, n_vmfc, 2, COUNTERPOISE_SSFC), &
                 "the SSFC check should accept the VMFC list of two monomers")
   end subroutine test_ssfc_two_monomers

   subroutine test_ssfc_order(error)
      !! The list is ordered by real monomers, largest first
      !!
      !! Every full-basis row names the whole system, so the size of a row says
      !! nothing about its cost; the real monomers set the electron count. The
      !! summation does not depend on the order -- it groups by real monomers
      !! itself -- so this is about the load balance, not the answer.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_terms, i

      call make_line(sys_geom, 5, 4.0_dp)
      config%nlevel = 3
      config%counterpoise = "ssfc"
      call generate_mbe_term_list(sys_geom, config, 3, polymers, n_terms)

      do i = 2, n_terms
         call check(error, real_count_of(polymers(i, :)) <= real_count_of(polymers(i - 1, :)), &
                    "a row with more real monomers follows one with fewer")
         if (allocated(error)) return
      end do
      call check(error, all(polymers(n_terms, 2:) == 0) .and. polymers(n_terms, 1) > 0, &
                 "the cheapest row is an own-basis monomer, last of all")
   end subroutine test_ssfc_order

   subroutine check_full_basis_list(error, polymers, n_terms, n_monomers, level)
      !! What the full-basis scheme promises of a list, row by row
      !!
      !! Every ghosted row names each monomer of the system once, so its
      !! ghosts are exactly the system's complement of its real monomers. No
      !! all-real row has two to `n_monomers - 1` monomers, and the whole
      !! system is a row only at full level, all-real. Each own-basis monomer
      !! has its full-basis row. The scheme's own check accepts the list.
      type(error_type), allocatable, intent(inout) :: error
      integer, intent(in) :: polymers(:, :)
      integer(int64), intent(in) :: n_terms
      integer, intent(in) :: n_monomers, level

      integer(int64) :: i, n_whole, n_own, n_full_basis_monomers
      integer :: n_real, k
      integer :: support(n_monomers)

      call check(error, size(polymers, 2) >= n_monomers, "the rows are as wide as the system")
      if (allocated(error)) return

      n_whole = 0_int64
      n_own = 0_int64
      n_full_basis_monomers = 0_int64
      do i = 1, n_terms
         n_real = real_count_of(polymers(i, :))
         if (any(polymers(i, :) < 0)) then
            call check(error, count(polymers(i, :) /= 0) == n_monomers, &
                       "a ghosted row names every monomer of the system")
            if (allocated(error)) return
            support = abs(polymers(i, 1:n_monomers))
            call ascending(support)
            call check(error, all(support == [(k, k=1, n_monomers)]), &
                       "a ghosted row ghosts exactly the system's complement of its real monomers")
            if (allocated(error)) return
            if (n_real == 1) n_full_basis_monomers = n_full_basis_monomers + 1_int64
         else if (n_real == 1) then
            n_own = n_own + 1_int64
         else
            call check(error, n_real == n_monomers, &
                       "no all-real row of 2 to N-1 monomers remains")
            if (allocated(error)) return
            n_whole = n_whole + 1_int64
         end if
      end do

      call check(error, n_whole == merge(1_int64, 0_int64, level == n_monomers), &
                 "the whole system is a row at full level and only then")
      if (allocated(error)) return
      call check(error, n_own == n_full_basis_monomers, &
                 "each own-basis monomer row has a full-basis row")
      if (allocated(error)) return
      call check(error, rows_match_counterpoise(polymers, n_terms, n_monomers, COUNTERPOISE_SSFC), &
                 "the full-basis check should accept the list the generator made")
   end subroutine check_full_basis_list

   pure function count_rows_of_size(polymers, n_terms, n) result(n_rows)
      !! How many of the first `n_terms` rows name exactly `n` monomers
      integer, intent(in) :: polymers(:, :)
      integer(int64), intent(in) :: n_terms
      integer, intent(in) :: n
      integer(int64) :: n_rows

      integer(int64) :: i

      n_rows = 0_int64
      do i = 1, n_terms
         if (count(polymers(i, :) /= 0) == n) n_rows = n_rows + 1_int64
      end do
   end function count_rows_of_size

   logical function has_parent(polymers, n_terms, support)
      !! Is this row's full support -- real and ghosted alike -- a real term
      !!
      !! Compared as a set, not term for term: a ghosted key lists its real
      !! monomers first and the ghosted ones after, so `abs()` of it is the
      !! parent's monomers in a different order.
      integer, intent(in) :: polymers(:, :)
      integer(int64), intent(in) :: n_terms
      integer, intent(in) :: support(:)

      integer(int64) :: j
      integer :: want(size(support)), have(size(support))

      want = support
      call ascending(want)

      has_parent = .false.
      do j = 1, n_terms
         if (is_auxiliary_row(polymers(j, :), COUNTERPOISE_VMFC)) cycle
         have = polymers(j, :)
         call ascending(have)
         if (all(have == want)) then
            has_parent = .true.
            return
         end if
      end do
   end function has_parent

   subroutine ascending(a)
      !! Insertion sort; the arrays here are four wide
      integer, intent(inout) :: a(:)

      integer :: i, j, key

      do i = 2, size(a)
         key = a(i)
         j = i - 1
         do while (j >= 1)
            if (a(j) <= key) exit
            a(j + 1) = a(j)
            j = j - 1
         end do
         a(j + 1) = key
      end do
   end subroutine ascending

   subroutine make_line(sys_geom, n, spacing)
      !! `n` one-atom monomers in a line, `spacing` Bohr apart
      type(system_geometry_t), intent(out) :: sys_geom
      integer, intent(in) :: n
      real(dp), intent(in) :: spacing

      integer :: i

      sys_geom%n_monomers = n
      sys_geom%atoms_per_monomer = 1
      sys_geom%total_atoms = n
      sys_geom%charge = 0
      sys_geom%multiplicity = 1
      allocate (sys_geom%element_numbers(n))
      allocate (sys_geom%coordinates(3, n))
      sys_geom%element_numbers = 10
      sys_geom%coordinates = 0.0_dp
      do i = 1, n
         sys_geom%coordinates(1, i) = real(i - 1, dp)*spacing
      end do
   end subroutine make_line

   pure function reference_count(n, level) result(kept)
      !! The reduced list's length with no screening: every term of up to
      !! `level` monomers holding the reference, and every term of up to
      !! `level - 1` without it
      !!
      !!     sum_{k=0}^{L-1} C(n-1, k)  +  sum_{k=1}^{L-1} C(n-1, k)
      integer, intent(in) :: n, level
      integer(int64) :: kept

      integer :: k

      kept = 1_int64
      do k = 1, level - 1
         kept = kept + 2_int64*binomial(n - 1, k)
      end do
   end function reference_count

   subroutine test_reference_count(error)
      !! The reduced list has the length the subset closure predicts
      !!
      !! Every system size from 2 to 7, every level up to 4 that fits, and
      !! every choice of reference -- the first, the last and all between --
      !! so a rule that works only for fragment 1, or only at level 2, fails.
      !! `n_full` must be the ordinary expansion's count, which `get_nfrags`
      !! is.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_terms, n_full
      integer :: n, level, ref

      do n = 2, 7
         call make_line(sys_geom, n, 4.0_dp)
         do level = 2, min(4, n)
            do ref = 1, n
               config%nlevel = level
               config%reference_fragment = ref
               call generate_mbe_term_list(sys_geom, config, level, polymers, n_terms, n_full=n_full)
               call check(error, n_terms == reference_count(n, level), &
                          "reference "//trim(int_str(int(ref, int64)))//" of "//trim(int_str(int(n, int64)))// &
                          " at level "//trim(int_str(int(level, int64)))//" kept "//trim(int_str(n_terms))// &
                          " terms, not "//trim(int_str(reference_count(n, level))))
               if (allocated(error)) return
               call check(error, n_full == get_nfrags(n, level), &
                          "n_full should be the ordinary expansion's count")
               if (allocated(error)) return
            end do
         end do
      end do
   end subroutine test_reference_count

   subroutine test_reference_closure(error)
      !! The reduced list is closed under subsets, holds every term with the
      !! reference in it, and nothing that none of those needs
      !!
      !! Closure is what `compute_mbe_delta` relies on: a missing subset aborts
      !! the run. Every term holding the reference is what the answer is made
      !! of. Minimality is the saving: a term not holding the reference is kept
      !! only when the same term with the reference added is kept, since that
      !! is the term whose correction needs it.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config, plain
      integer, allocatable :: polymers(:, :), full(:, :)
      integer(int64) :: n_terms, n_all
      integer :: ref, level

      call make_line(sys_geom, 6, 4.0_dp)
      do level = 2, 4
         plain%nlevel = level
         call generate_mbe_term_list(sys_geom, plain, level, full, n_all)
         do ref = 1, 6
            config%nlevel = level
            config%reference_fragment = ref
            call generate_mbe_term_list(sys_geom, config, level, polymers, n_terms)
            call check_reduced_list(error, polymers, n_terms, full, n_all, ref, level)
            if (allocated(error)) return
         end do
      end do
   end subroutine test_reference_closure

   subroutine test_reference_after_screening(error)
      !! The reduction on a screened list is still closed and still minimal
      !!
      !! Screening removes terms first, so the reduced list is not the
      !! combinatorial one: a subset is needed only by the terms holding the
      !! reference that survived the screen. Checked against the screened list
      !! the ordinary run would have used.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config, plain
      integer, allocatable :: polymers(:, :), full(:, :)
      integer(int64) :: n_terms, n_all, n_full
      integer :: ref

      ! 4 Bohr is 2.117 Angstrom: a 5 Angstrom dimer cutoff keeps neighbours
      ! up to two apart, and a 3 Angstrom trimer cutoff only adjacent triples.
      call make_line(sys_geom, 6, 4.0_dp)
      allocate (plain%fragment_cutoffs(3))
      plain%fragment_cutoffs = [0.0_dp, 5.0_dp, 3.0_dp]
      plain%nlevel = 3
      call generate_mbe_term_list(sys_geom, plain, 3, full, n_all)

      do ref = 1, 6
         config = plain
         config%reference_fragment = ref
         call generate_mbe_term_list(sys_geom, config, 3, polymers, n_terms, n_full=n_full)
         call check(error, n_full == n_all, "n_full should be the screened ordinary count")
         if (allocated(error)) return
         call check(error, n_terms < n_all, "the reduction should skip something here")
         if (allocated(error)) return
         call check_reduced_list(error, polymers, n_terms, full, n_all, ref, 3)
         if (allocated(error)) return
      end do
   end subroutine test_reference_after_screening

   subroutine check_reduced_list(error, polymers, n_terms, full, n_all, ref, level)
      !! The three properties `apply_reference_closure` promises
      type(error_type), allocatable, intent(inout) :: error
      integer, intent(in) :: polymers(:, :), full(:, :)
      integer(int64), intent(in) :: n_terms, n_all
      integer, intent(in) :: ref, level

      type(fragment_lookup_t) :: kept, everything
      integer(int64) :: i
      integer :: n, mask, j, k
      integer :: sub(level), with_ref(level + 1)

      call kept%init(n_terms)
      do i = 1, n_terms
         call kept%insert(polymers(i, :), count(polymers(i, :) /= 0), i)
      end do
      call everything%init(n_all)
      do i = 1, n_all
         call everything%insert(full(i, :), count(full(i, :) /= 0), i)
      end do

      ! Closed: every proper, non-empty subset of every kept term is kept.
      do i = 1, n_terms
         n = count(polymers(i, :) /= 0)
         do mask = 1, 2**n - 2
            k = 0
            do j = 1, n
               if (btest(mask, j - 1)) then
                  k = k + 1
                  sub(k) = polymers(i, j)
               end if
            end do
            call check(error, kept%find(sub(1:k), k) > 0, &
                       "a subset of a kept term is missing, which compute_mbe_delta would abort on")
            if (allocated(error)) return
         end do
      end do

      ! Complete: every term of the ordinary list that holds the reference.
      do i = 1, n_all
         n = count(full(i, :) /= 0)
         if (.not. any(full(i, 1:n) == ref)) cycle
         call check(error, kept%find(full(i, 1:n), n) > 0, &
                    "a term holding the reference was skipped")
         if (allocated(error)) return
      end do

      ! Minimal: a kept term without the reference is a subset of a kept one with it.
      do i = 1, n_terms
         n = count(polymers(i, :) /= 0)
         if (any(polymers(i, 1:n) == ref)) cycle
         call check(error, n < level, "a term of the full level without the reference was kept")
         if (allocated(error)) return
         with_ref(1:n) = polymers(i, 1:n)
         with_ref(n + 1) = ref
         call check(error, kept%find(with_ref(1:n + 1), n + 1) > 0, &
                    "a term was kept that no term holding the reference needs")
         if (allocated(error)) return
      end do

      call kept%destroy()
      call everything%destroy()
   end subroutine check_reduced_list

   function int_str(n) result(s)
      !! The count, for a failure message that says what it actually found
      integer(int64), intent(in) :: n
      character(len=32) :: s

      write (s, "(i0)") n
      s = adjustl(s)
   end function int_str

end module test_mqc_term_list

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_term_list, only: collect_mqc_term_list_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_term_list", collect_mqc_term_list_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if

end program tester
