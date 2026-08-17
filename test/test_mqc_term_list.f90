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
   use mqc_frag_utils, only: generate_mbe_term_list, get_nfrags
   use mqc_combinatorics, only: is_auxiliary_row, real_count_of
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
                  new_unittest("counterpoise_follows_screening", test_vmfc_after_screening) &
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
         if (is_auxiliary_row(polymers(i, :))) then
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
         if (.not. is_auxiliary_row(polymers(i, :))) cycle
         call check(error, has_parent(polymers, n_terms, abs(polymers(i, :))), &
                    "a ghosted row has no parent n-mer in the list")
         if (allocated(error)) return
      end do

      ! And no auxiliary row may be all-ghost: something has to be real.
      do i = 1, n_terms
         if (.not. is_auxiliary_row(polymers(i, :))) cycle
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
         if (is_auxiliary_row(polymers(i, :))) then
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
         if (is_auxiliary_row(polymers(j, :))) cycle
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
