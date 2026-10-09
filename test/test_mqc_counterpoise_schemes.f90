module test_mqc_counterpoise_schemes
   !! Which counterpoise scheme a term list belongs to, and what each scheme sums
   !!
   !! A term list carries the scheme in its rows, and two schemes can share a
   !! row: at L = N the SSFC top term is the whole system, an all-real row
   !! exactly like VMFC's. So the scheme is named by the caller and the rows are
   !! held to it by `rows_match_counterpoise`, and which rows are summed follows
   !! the scheme through `is_auxiliary_row`.
   !!
   !! Every SSFC list below is built by hand, by `ssfc_rows` in this file, from
   !! the definition of the scheme rather than by production code: one own-basis
   !! row per monomer; at level 2 and above one full-basis row per monomer, with
   !! every other monomer of the system ghosted; and one full-basis row per term
   !! of two or more monomers, the system's complement ghosted, the whole system
   !! being an all-real row when the level reaches it. VMFC and plain lists come
   !! from `generate_mbe_term_list`, the production generator.
   !!
   !! The energies are made up, symmetric in the row and different for every
   !! ghost pattern, so a term looked up in the wrong basis moves the total.
   !! The SSFC totals are compared with inclusion-exclusion done directly over
   !! full-basis energies, which uses none of the recursion under test, and at
   !! L = N with the Boys-Bernardi identity
   !!
   !!     E = E_cluster + sum_i [ E_i(i) - E_i(FB) ]
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp, int64
   use mqc_mbe, only: compute_mbe
   use mqc_result_types, only: calculation_result_t, mbe_result_t
   use mqc_frag_utils, only: generate_mbe_term_list
   use mqc_combinatorics, only: is_auxiliary_row, rows_match_counterpoise, &
                                counterpoise_scheme_of, counterpoise_scheme_name, &
                                COUNTERPOISE_NONE, COUNTERPOISE_VMFC, COUNTERPOISE_SSFC
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_config_adapter, only: driver_config_t
   implicit none
   private
   public :: collect_mqc_counterpoise_schemes_tests

   real(dp), parameter :: TOLERANCE = 1.0e-12_dp

contains

   subroutine collect_mqc_counterpoise_schemes_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("the_scheme_is_named_by_its_spelling", test_scheme_names), &
                  new_unittest("only_one_real_monomer_is_auxiliary_under_ssfc", test_predicate_table), &
                  new_unittest("each_scheme_accepts_its_own_lists_only", test_list_matrix), &
                  new_unittest("screened_lists_keep_their_scheme", test_screened_lists), &
                  new_unittest("a_damaged_list_matches_no_scheme", test_damaged_lists), &
                  new_unittest("ssfc_at_the_top_level_sums_the_whole_system", test_ssfc_top_level), &
                  new_unittest("ssfc_truncated_follows_full_basis_inclusion_exclusion", test_ssfc_truncated), &
                  new_unittest("ssfc_with_no_pairs_is_the_own_basis_monomers", test_ssfc_no_pairs), &
                  new_unittest("level_one_adds_no_rows_in_either_scheme", test_level_one), &
                  new_unittest("two_monomers_vmfc_and_ssfc_coincide", test_two_monomers), &
                  new_unittest("ssfc_reference_with_screening", test_ssfc_reference) &
                  ]
   end subroutine collect_mqc_counterpoise_schemes_tests

   !---------------------------------------------------------------------------
   ! Fixtures
   !---------------------------------------------------------------------------

   subroutine make_chain(sys_geom, n_monomers)
      !! `n_monomers` one-atom monomers in a line, far enough apart for no screening
      type(system_geometry_t), intent(out) :: sys_geom
      integer, intent(in) :: n_monomers

      integer :: i

      sys_geom%n_monomers = n_monomers
      sys_geom%atoms_per_monomer = 1
      sys_geom%total_atoms = n_monomers
      sys_geom%charge = 0
      sys_geom%multiplicity = 1
      allocate (sys_geom%element_numbers(n_monomers))
      allocate (sys_geom%coordinates(3, n_monomers))
      sys_geom%element_numbers = 10  ! neon: closed shell, never bonded
      sys_geom%coordinates = 0.0_dp
      do i = 1, n_monomers
         sys_geom%coordinates(1, i) = real(i - 1, dp)*4.0_dp
      end do
   end subroutine make_chain

   subroutine vmfc_or_plain_rows(n_monomers, level, scheme_name, polymers, count, cutoff, reference)
      !! The production term list: ordinary, or with VMFC's ghosted rows added
      integer, intent(in) :: n_monomers, level
      character(len=*), intent(in) :: scheme_name
      integer, allocatable, intent(out) :: polymers(:, :)
      integer(int64), intent(out) :: count
      real(dp), intent(in), optional :: cutoff
         !! Pair cutoff in Angstrom; the chain steps are 4 Bohr, about 2.12
         !! Angstrom, so 3.0 keeps the adjacent pairs and nothing further
      integer, intent(in), optional :: reference
         !! Reduce the list to this monomer's terms, 1-based, as an
         !! interaction energy does

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config

      call make_chain(sys_geom, n_monomers)
      config%nlevel = level
      config%counterpoise = scheme_name
      if (present(cutoff)) then
         allocate (config%fragment_cutoffs(level))
         config%fragment_cutoffs = 0.0_dp
         config%fragment_cutoffs(2:) = cutoff
      end if
      if (present(reference)) config%reference_fragment = reference
      call generate_mbe_term_list(sys_geom, config, level, polymers, count)
   end subroutine vmfc_or_plain_rows

   function row_of(mask, n_monomers) result(row)
      !! The full-basis row of the monomers in `mask`: those real, the rest of
      !! the system ghosted. Width `n_monomers`; the whole system is all real.
      integer, intent(in) :: mask, n_monomers
      integer :: row(n_monomers)

      integer :: m, k

      row = 0
      k = 0
      do m = 1, n_monomers
         if (btest(mask, m - 1)) then
            k = k + 1
            row(k) = m
         end if
      end do
      do m = 1, n_monomers
         if (.not. btest(mask, m - 1)) then
            k = k + 1
            row(k) = -m
         end if
      end do
   end function row_of

   subroutine ssfc_rows(n_monomers, level, with_terms, polymers, count)
      !! An SSFC term list, built by hand from its definition
      !!
      !! Own-basis monomers; at `level` >= 2 a full-basis row for every
      !! monomer, and, unless `with_terms` is false (screening kept no pairs), a
      !! full-basis row for every term of 2 to `level` monomers.
      integer, intent(in) :: n_monomers, level
      logical, intent(in) :: with_terms
      integer, allocatable, intent(out) :: polymers(:, :)
      integer(int64), intent(out) :: count

      integer :: mask, size_of_term, m

      allocate (polymers(2*2**n_monomers, n_monomers))
      polymers = 0
      count = 0_int64
      do m = 1, n_monomers
         count = count + 1_int64
         polymers(count, 1) = m
      end do
      if (level < 2) return

      do mask = 1, 2**n_monomers - 1
         size_of_term = popcnt(mask)
         if (size_of_term > level) cycle
         if (size_of_term >= 2 .and. .not. with_terms) cycle
         count = count + 1_int64
         polymers(count, :) = row_of(mask, n_monomers)
      end do
   end subroutine ssfc_rows

   function made_up(row) result(e)
      !! An energy for a row, symmetric in it and different for every ghost pattern
      integer, intent(in) :: row(:)
      real(dp) :: e

      integer :: a, b, c

      e = 0.0_dp
      do a = 1, size(row)
         if (row(a) == 0) cycle
         if (row(a) > 0) then
            e = e - 10.0_dp - 0.1_dp*real(row(a), dp)
         else
            e = e - 1.3_dp + 0.07_dp*real(row(a), dp)
         end if
         do b = a + 1, size(row)
            if (row(b) == 0) cycle
            e = e + 0.01_dp*real(mod(abs(row(a)*row(b)), 7) + 1, dp) &
                *(1.0_dp + 0.5_dp*real(count([row(a), row(b)] < 0), dp))
            do c = b + 1, size(row)
               if (row(c) == 0) cycle
               e = e + 1.0e-3_dp*real(abs(row(a)*row(b)*row(c)), dp)/7.0_dp
            end do
         end do
      end do
   end function made_up

   subroutine fill_results(polymers, count, results)
      !! Made-up energies for every row, and dipoles that are linear in them
      integer, intent(in) :: polymers(:, :)
      integer(int64), intent(in) :: count
      type(calculation_result_t), allocatable, intent(out) :: results(:)

      integer(int64) :: i

      allocate (results(count))
      do i = 1_int64, count
         results(i)%has_energy = .true.
         results(i)%energy%scf = made_up(polymers(i, :))
         allocate (results(i)%dipole(3))
         results(i)%dipole = [1.0_dp, 2.0_dp, -0.5_dp]*results(i)%energy%scf
         results(i)%has_dipole = .true.
      end do
   end subroutine fill_results

   subroutine run_expansion(polymers, count, level, scheme, mbe_result, reference, sys_geom)
      !! `compute_mbe` over the made-up energies of a list
      integer, intent(in) :: polymers(:, :)
      integer(int64), intent(in) :: count
      integer, intent(in) :: level, scheme
      type(mbe_result_t), intent(inout) :: mbe_result
      integer, intent(in), optional :: reference
      type(system_geometry_t), intent(in), optional :: sys_geom

      type(calculation_result_t), allocatable :: results(:)

      call fill_results(polymers, count, results)
      call mbe_result%allocate_dipole()
      call compute_mbe(polymers, count, level, results, mbe_result, sys_geom=sys_geom, &
                       reference=reference, counterpoise_scheme=scheme)
   end subroutine run_expansion

   function ssfc_expected(n_monomers, level) result(total)
      !! The SSFC total by inclusion-exclusion over full-basis energies
      !!
      !! Own-basis monomers, plus for every term T of 2 to `level` monomers the
      !! alternating sum of the full-basis energies of its subsets. No
      !! recursion and no lookup, so it shares nothing with `compute_mbe`.
      integer, intent(in) :: n_monomers, level
      real(dp) :: total

      integer :: term, sub, m
      real(dp) :: increment

      total = 0.0_dp
      do m = 1, n_monomers
         total = total + made_up([m])
      end do
      do term = 1, 2**n_monomers - 1
         if (popcnt(term) < 2 .or. popcnt(term) > level) cycle
         increment = 0.0_dp
         do sub = 1, 2**n_monomers - 1
            if (iand(sub, term) /= sub) cycle
            increment = increment + merge(1.0_dp, -1.0_dp, mod(popcnt(term) - popcnt(sub), 2) == 0) &
                        *made_up(row_of(sub, n_monomers))
         end do
         total = total + increment
      end do
   end function ssfc_expected

   function same_row_set(a, na, b, nb) result(same)
      !! Whether two lists hold the same rows, as sets of signed monomers
      integer, intent(in) :: a(:, :), b(:, :)
      integer(int64), intent(in) :: na, nb
      logical :: same

      integer(int64) :: i, j
      logical :: found

      same = na == nb
      if (.not. same) return
      do i = 1_int64, na
         found = .false.
         do j = 1_int64, nb
            if (as_set(a(i, :), b(j, :))) then
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            same = .false.
            return
         end if
      end do
   end function same_row_set

   function as_set(x, y) result(equal)
      !! Whether two zero-padded rows name the same signed monomers
      integer, intent(in) :: x(:), y(:)
      logical :: equal

      integer :: i

      equal = count(x /= 0) == count(y /= 0)
      if (.not. equal) return
      do i = 1, size(x)
         if (x(i) == 0) cycle
         if (.not. any(y == x(i))) then
            equal = .false.
            return
         end if
      end do
   end function as_set

   subroutine expect_matches(error, polymers, count, n_monomers, accepted_by, what)
      !! Each scheme accepts the list exactly when `accepted_by` says it does
      type(error_type), allocatable, intent(inout) :: error
      integer, intent(in) :: polymers(:, :)
      integer(int64), intent(in) :: count
      integer, intent(in) :: n_monomers
      logical, intent(in) :: accepted_by(0:2)
         !! Indexed by the `COUNTERPOISE_*` constant
      character(len=*), intent(in) :: what

      integer :: scheme

      do scheme = COUNTERPOISE_NONE, COUNTERPOISE_SSFC
         call check(error, rows_match_counterpoise(polymers, count, n_monomers, scheme) .eqv. &
                    accepted_by(scheme), &
                    what//": the "//counterpoise_scheme_name(scheme)//" check should "// &
                    merge("accept", "reject", accepted_by(scheme))//" it")
         if (allocated(error)) return
      end do
   end subroutine expect_matches

   !---------------------------------------------------------------------------
   ! Names and the predicate
   !---------------------------------------------------------------------------

   subroutine test_scheme_names(error)
      !! The deck's spelling maps to the constant and back
      type(error_type), allocatable, intent(out) :: error

      call check(error, counterpoise_scheme_of("vmfc") == COUNTERPOISE_VMFC, "vmfc names VMFC")
      if (allocated(error)) return
      call check(error, counterpoise_scheme_of("none") == COUNTERPOISE_NONE, "none names no scheme")
      if (allocated(error)) return
      call check(error, counterpoise_scheme_of("vmfc   ") == COUNTERPOISE_VMFC, &
                 "the padding of a fixed-length string is not part of the name")
      if (allocated(error)) return
      call check(error, counterpoise_scheme_of("ssfc") == COUNTERPOISE_SSFC, "ssfc names SSFC")
      if (allocated(error)) return
      call check(error, counterpoise_scheme_name(COUNTERPOISE_SSFC), "ssfc", &
                 "the SSFC constant prints as ssfc")
   end subroutine test_scheme_names

   subroutine test_predicate_table(error)
      !! The summation table: which rows are subtracted and never summed
      !!
      !!     row                     VMFC        SSFC
      !!     all real                summed      summed
      !!     ghosted, 1 real         auxiliary   auxiliary
      !!     ghosted, 2+ real        auxiliary   summed
      type(error_type), allocatable, intent(out) :: error

      ! All real, padded or not: summed by both
      call check(error,.not. is_auxiliary_row([1, 0, 0, 0], COUNTERPOISE_VMFC), "VMFC sums a monomer")
      if (allocated(error)) return
      call check(error,.not. is_auxiliary_row([1, 0, 0, 0], COUNTERPOISE_SSFC), "SSFC sums a monomer")
      if (allocated(error)) return
      call check(error,.not. is_auxiliary_row([1, 2, 3, 4], COUNTERPOISE_VMFC), "VMFC sums an all-real n-mer")
      if (allocated(error)) return
      call check(error,.not. is_auxiliary_row([1, 2, 3, 4], COUNTERPOISE_SSFC), &
                 "SSFC sums the whole system at L = N")
      if (allocated(error)) return

      ! Ghosted, one real monomer: auxiliary in both, whichever side the ghosts sit
      call check(error, is_auxiliary_row([1, -2, -3, -4], COUNTERPOISE_VMFC), "VMFC: full-basis monomer")
      if (allocated(error)) return
      call check(error, is_auxiliary_row([1, -2, -3, -4], COUNTERPOISE_SSFC), "SSFC: full-basis monomer")
      if (allocated(error)) return
      call check(error, is_auxiliary_row([-1, 2, -3, -4], COUNTERPOISE_SSFC), &
                 "SSFC: the real monomer need not come first")
      if (allocated(error)) return
      call check(error, is_auxiliary_row([1, -2, 0, 0], COUNTERPOISE_VMFC), "VMFC: monomer in a pair basis")
      if (allocated(error)) return

      ! Ghosted, two or more real: the scheme decides
      call check(error, is_auxiliary_row([1, 2, -3, 0], COUNTERPOISE_VMFC), &
                 "VMFC: a pair in the trimer basis is subtracted by the trimer")
      if (allocated(error)) return
      call check(error,.not. is_auxiliary_row([1, 2, -3, -4], COUNTERPOISE_SSFC), &
                 "SSFC: a pair in the cluster basis is a term")
      if (allocated(error)) return
      call check(error,.not. is_auxiliary_row([1, 2, 3, -4], COUNTERPOISE_SSFC), &
                 "SSFC: a trimer in the cluster basis is a term")
      if (allocated(error)) return

      ! No scheme has nothing to subtract
      call check(error,.not. is_auxiliary_row([1, 2, 0, 0], COUNTERPOISE_NONE), "NONE sums a pair")
      if (allocated(error)) return
      call check(error,.not. is_auxiliary_row([1, 0, 0, 0], COUNTERPOISE_NONE), "NONE sums a monomer")
   end subroutine test_predicate_table

   !---------------------------------------------------------------------------
   ! The assertion
   !---------------------------------------------------------------------------

   subroutine test_list_matrix(error)
      !! Each scheme accepts its own list and rejects the other two
      !!
      !! Four monomers, levels 2 to 4. Level 4 is the case inference cannot
      !! settle: the SSFC list holds the whole system as an all-real row, and
      !! its ghosted rows are the ghosted subsets of that row, which is also
      !! what the VMFC check looks for beneath an all-real n-mer.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N = 4
      integer, allocatable :: plain(:, :), vmfc(:, :), ssfc(:, :)
      integer(int64) :: n_plain, n_vmfc, n_ssfc
      integer :: level

      do level = 2, N
         call vmfc_or_plain_rows(N, level, "none", plain, n_plain)
         call vmfc_or_plain_rows(N, level, "vmfc", vmfc, n_vmfc)
         call ssfc_rows(N, level, .true., ssfc, n_ssfc)

         call expect_matches(error, plain, n_plain, N, [.true., .false., .false.], &
                             "an ordinary list at level "//digit(level))
         if (allocated(error)) return
         call expect_matches(error, vmfc, n_vmfc, N, [.false., .true., .false.], &
                             "a VMFC list at level "//digit(level))
         if (allocated(error)) return
         call expect_matches(error, ssfc, n_ssfc, N, [.false., .false., .true.], &
                             "an SSFC list at level "//digit(level))
         if (allocated(error)) return
      end do
   end subroutine test_list_matrix

   subroutine test_screened_lists(error)
      !! Screening changes which terms a list holds, not which scheme it is
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N = 4
      integer, allocatable :: vmfc(:, :), ssfc(:, :), plain(:, :)
      integer(int64) :: n_vmfc, n_ssfc, n_plain

      ! The 3 Angstrom cutoff keeps the three adjacent pairs of the chain.
      call vmfc_or_plain_rows(N, 3, "vmfc", vmfc, n_vmfc, cutoff=3.0_dp)
      call vmfc_or_plain_rows(N, 3, "none", plain, n_plain, cutoff=3.0_dp)
      call expect_matches(error, vmfc, n_vmfc, N, [.false., .true., .false.], &
                          "a screened VMFC list")
      if (allocated(error)) return
      call expect_matches(error, plain, n_plain, N, [.true., .false., .false.], &
                          "a screened ordinary list")
      if (allocated(error)) return

      ! Screening that keeps no pair leaves the monomers and, for SSFC, their
      ! full-basis rows. The VMFC generator adds nothing for it.
      call ssfc_rows(N, 3, .false., ssfc, n_ssfc)
      call expect_matches(error, ssfc, n_ssfc, N, [.false., .false., .true.], &
                          "an SSFC list that screening left with no pairs")
      if (allocated(error)) return
      call check(error, n_ssfc == 2*N, "monomers own-basis and full-basis, and nothing else")
      if (allocated(error)) return

      ! A reference fragment reduces the list to its own terms and their
      ! subsets. With the 3 Angstrom cutoff, reference 2 keeps monomers 1 to 3
      ! and the pairs (1,2) and (2,3): pair (3,4) and monomer 4 are gone. The
      ! reduced list has the ghosted subsets of the pairs it kept and no others.
      call vmfc_or_plain_rows(N, 3, "vmfc", vmfc, n_vmfc, cutoff=3.0_dp, reference=2)
      call vmfc_or_plain_rows(N, 3, "none", plain, n_plain, cutoff=3.0_dp, reference=2)
      call check(error, n_plain == 5_int64, "reference 2 with screening keeps three monomers and two pairs")
      if (allocated(error)) return
      call check(error, n_vmfc == 9_int64, "and VMFC adds the two ghosted subsets of each pair")
      if (allocated(error)) return
      call expect_matches(error, vmfc, n_vmfc, N, [.false., .true., .false.], &
                          "a VMFC list reduced to a reference fragment and screened")
      if (allocated(error)) return
      call expect_matches(error, plain, n_plain, N, [.true., .false., .false.], &
                          "an ordinary list reduced to a reference fragment and screened")
      if (allocated(error)) return

      ! Without screening the reduction drops only what no kept term needs:
      ! reference 1 at level 3 keeps every term holding 1 and their subsets,
      ! which is all of the 14 terms but the triple (2,3,4).
      call vmfc_or_plain_rows(N, 3, "vmfc", vmfc, n_vmfc, reference=1)
      call vmfc_or_plain_rows(N, 3, "none", plain, n_plain, reference=1)
      call check(error, n_plain == 13_int64, "reference 1 at level 3 drops the triple (2,3,4)")
      if (allocated(error)) return
      call expect_matches(error, vmfc, n_vmfc, N, [.false., .true., .false.], &
                          "a VMFC list reduced to a reference fragment")
      if (allocated(error)) return

      call vmfc_or_plain_rows(N, 3, "vmfc", vmfc, n_vmfc, cutoff=0.5_dp)
      call check(error, n_vmfc == N, "a cutoff that drops every pair leaves only the monomers")
      if (allocated(error)) return
      call expect_matches(error, vmfc, n_vmfc, N, [.true., .true., .true.], &
                          "monomers alone are the same list under every scheme")
   end subroutine test_screened_lists

   subroutine test_damaged_lists(error)
      !! A list that is missing a row, repeats one or holds a stray one is no scheme's
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N = 4
      integer, allocatable :: list(:, :)
      integer(int64) :: n_rows
      integer :: row

      ! VMFC: drop each row in turn, then repeat each in turn. Dropping a
      ! monomer or an n-mer, a ghosted row or a parent, breaks the count or the
      ! closure either way.
      call vmfc_or_plain_rows(N, 3, "vmfc", list, n_rows)
      call check(error, rows_match_counterpoise(list, n_rows, N, COUNTERPOISE_VMFC), &
                 "the undamaged VMFC list must match")
      if (allocated(error)) return
      do row = 1, int(n_rows)
         call check(error,.not. rows_match_counterpoise(drop_row(list, n_rows, row), n_rows - 1_int64, &
                                                        N, COUNTERPOISE_VMFC), &
                    "a VMFC list with row "//digit(row)//" missing was accepted")
         if (allocated(error)) return
         call check(error,.not. rows_match_counterpoise(repeat_row(list, n_rows, row), n_rows + 1_int64, &
                                                        N, COUNTERPOISE_VMFC), &
                    "a VMFC list with row "//digit(row)//" twice was accepted")
         if (allocated(error)) return
      end do

      ! SSFC, with each row dropped or repeated, and for a system of the wrong size.
      call ssfc_rows(N, 3, .true., list, n_rows)
      call check(error, rows_match_counterpoise(list, n_rows, N, COUNTERPOISE_SSFC), &
                 "the undamaged SSFC list must match")
      if (allocated(error)) return
      do row = 1, int(n_rows)
         ! Dropping a term of the top level leaves a list that screening could
         ! have produced; dropping anything beneath one leaves a hole.
         call check(error, rows_match_counterpoise(drop_row(list, n_rows, row), n_rows - 1_int64, &
                                                   N, COUNTERPOISE_SSFC) .eqv. &
                    (count(list(row, :) > 0) == 3), &
                    "an SSFC list with row "//digit(row)//" missing was judged wrongly")
         if (allocated(error)) return
         call check(error,.not. rows_match_counterpoise(repeat_row(list, n_rows, row), n_rows + 1_int64, &
                                                        N, COUNTERPOISE_SSFC), &
                    "an SSFC list with row "//digit(row)//" twice was accepted")
         if (allocated(error)) return
      end do
      call check(error,.not. rows_match_counterpoise(list, n_rows, N + 1, COUNTERPOISE_SSFC), &
                 "an SSFC list for four monomers was accepted for a system of five")
      if (allocated(error)) return

      ! A scheme the code does not know matches nothing.
      call check(error,.not. rows_match_counterpoise(list, n_rows, N, 99), "an unknown scheme matched")
   end subroutine test_damaged_lists

   function drop_row(list, n_rows, row) result(shorter)
      integer, intent(in) :: list(:, :)
      integer(int64), intent(in) :: n_rows
      integer, intent(in) :: row
      integer, allocatable :: shorter(:, :)

      allocate (shorter(n_rows - 1_int64, size(list, 2)))
      shorter(1:row - 1, :) = list(1:row - 1, :)
      shorter(row:, :) = list(row + 1:n_rows, :)
   end function drop_row

   function repeat_row(list, n_rows, row) result(longer)
      integer, intent(in) :: list(:, :)
      integer(int64), intent(in) :: n_rows
      integer, intent(in) :: row
      integer, allocatable :: longer(:, :)

      allocate (longer(n_rows + 1_int64, size(list, 2)))
      longer(1:n_rows, :) = list(1:n_rows, :)
      longer(n_rows + 1_int64, :) = list(row, :)
   end function repeat_row

   function digit(n) result(text)
      integer, intent(in) :: n
      character(len=:), allocatable :: text

      character(len=12) :: buffer

      write (buffer, "(i0)") n
      text = trim(buffer)
   end function digit

   !---------------------------------------------------------------------------
   ! The sums
   !---------------------------------------------------------------------------

   subroutine test_ssfc_top_level(error)
      !! SSFC at L = N: the whole-system row is summed, and so is every ghosted
      !! row with two or more real monomers
      !!
      !! Three monomers, all of it. The total must be the Boys-Bernardi
      !! counterpoise-corrected cluster energy,
      !! `E_cluster + sum_i [E_i(i) - E_i(FB)]`, which the telescoping of the
      !! many-body sum guarantees. A rule that left `[1,2,-3]` unsummed, as VMFC
      !! does, returns the sum of the own-basis monomers and the whole system
      !! alone, and misses the dimer terms.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N = 3
      integer, allocatable :: polymers(:, :)
      integer(int64) :: count
      type(mbe_result_t) :: mbe_result
      type(system_geometry_t) :: sys_geom
      real(dp) :: boys_bernardi
      integer :: i

      call make_chain(sys_geom, N)
      call ssfc_rows(N, N, .true., polymers, count)
      ! Three own-basis monomers, three full-basis monomers, three full-basis
      ! pairs, and the whole system.
      call check(error, count == 10_int64, "SSFC(3) over three monomers is ten rows")
      if (allocated(error)) return

      call run_expansion(polymers, count, N, COUNTERPOISE_SSFC, mbe_result, sys_geom=sys_geom)

      boys_bernardi = made_up([1, 2, 3])
      do i = 1, N
         boys_bernardi = boys_bernardi + made_up([i]) - made_up(row_of(2**(i - 1), N))
      end do
      call check(error, mbe_result%total_energy, boys_bernardi, thr=TOLERANCE, &
                 message="SSFC(N) is not the Boys-Bernardi counterpoise-corrected total")
      if (allocated(error)) return
      call check(error, mbe_result%total_energy, ssfc_expected(N, N), thr=TOLERANCE, &
                 message="SSFC(N) differs from inclusion-exclusion over full-basis energies")
      if (allocated(error)) return
      call check(error, mbe_result%dipole(2), 2.0_dp*boys_bernardi, thr=TOLERANCE, &
                 message="the SSFC dipole does not follow the energy's recursion")
      if (allocated(error)) return

      call mbe_result%destroy()
   end subroutine test_ssfc_top_level

   subroutine test_ssfc_truncated(error)
      !! SSFC(2) and SSFC(3) over four monomers, against full-basis inclusion-exclusion
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N = 4
      integer, allocatable :: polymers(:, :)
      integer(int64) :: count
      type(mbe_result_t) :: mbe_result
      integer :: level

      do level = 2, N
         call ssfc_rows(N, level, .true., polymers, count)
         call run_expansion(polymers, count, level, COUNTERPOISE_SSFC, mbe_result)
         call check(error, mbe_result%total_energy, ssfc_expected(N, level), thr=TOLERANCE, &
                    message="SSFC("//digit(level)//") over four monomers is not the "// &
                    "full-basis inclusion-exclusion total")
         if (allocated(error)) return
         call check(error, mbe_result%dipole(1), ssfc_expected(N, level), thr=TOLERANCE, &
                    message="SSFC("//digit(level)//") dipole differs from the energy's recursion")
         if (allocated(error)) return
         call mbe_result%destroy()
      end do
   end subroutine test_ssfc_truncated

   subroutine test_ssfc_no_pairs(error)
      !! Screening that keeps no pairs: the total is the own-basis monomers
      !!
      !! The list is the monomers and their full-basis rows. Every row with one
      !! real monomer is auxiliary, so no level above the first is summed.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N = 3
      integer, allocatable :: polymers(:, :)
      integer(int64) :: count
      type(mbe_result_t) :: mbe_result

      call ssfc_rows(N, 2, .false., polymers, count)
      call run_expansion(polymers, count, 2, COUNTERPOISE_SSFC, mbe_result)

      call check(error, mbe_result%total_energy, made_up([1]) + made_up([2]) + made_up([3]), &
                 thr=TOLERANCE, message="with no pairs kept SSFC is the own-basis monomers, "// &
                 "and the full-basis rows must not be summed")
      if (allocated(error)) return
      call check(error, mbe_result%total_energy, ssfc_expected(N, 1), thr=TOLERANCE, &
                 message="and so is inclusion-exclusion with nothing above order one")
      call mbe_result%destroy()
   end subroutine test_ssfc_no_pairs

   subroutine test_level_one(error)
      !! Level 1 adds no rows in either scheme, and is the sum of the monomers
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N = 3
      integer, allocatable :: plain(:, :), vmfc(:, :), ssfc(:, :)
      integer(int64) :: n_plain, n_vmfc, n_ssfc
      type(mbe_result_t) :: mbe_result
      real(dp) :: monomers
      integer :: scheme

      call vmfc_or_plain_rows(N, 1, "none", plain, n_plain)
      call vmfc_or_plain_rows(N, 1, "vmfc", vmfc, n_vmfc)
      call ssfc_rows(N, 1, .true., ssfc, n_ssfc)
      call check(error, n_vmfc == n_plain .and. n_plain == N, &
                 "VMFC at level 1 must add no ghosted rows")
      if (allocated(error)) return
      call check(error, n_ssfc == N, "SSFC at level 1 must add no rows")
      if (allocated(error)) return

      monomers = made_up([1]) + made_up([2]) + made_up([3])
      do scheme = COUNTERPOISE_NONE, COUNTERPOISE_SSFC
         call check(error, rows_match_counterpoise(plain, n_plain, N, scheme), &
                    "a level-1 list is every scheme's: "//counterpoise_scheme_name(scheme))
         if (allocated(error)) return
         call run_expansion(plain, n_plain, 1, scheme, mbe_result)
         call check(error, mbe_result%total_energy, monomers, thr=TOLERANCE, &
                    message="level 1 under "//counterpoise_scheme_name(scheme)// &
                    " is not the sum of the monomers")
         if (allocated(error)) return
         call mbe_result%destroy()
      end do
   end subroutine test_level_one

   subroutine test_two_monomers(error)
      !! N = 2: VMFC(2) and SSFC(2) are the same rows and the same energy
      !!
      !! Asserted so that nobody "fixes" the coincidence. The pair's subsets in
      !! the pair's basis are its subsets in the cluster's basis, because the
      !! pair is the cluster.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N = 2
      integer, allocatable :: vmfc(:, :), ssfc(:, :)
      integer(int64) :: n_vmfc, n_ssfc
      type(mbe_result_t) :: by_vmfc, by_ssfc
      real(dp) :: expected

      call vmfc_or_plain_rows(N, 2, "vmfc", vmfc, n_vmfc)
      call ssfc_rows(N, 2, .true., ssfc, n_ssfc)

      call check(error, n_vmfc == 5_int64, "VMFC(2) over two monomers is [1] [2] [1,-2] [-1,2] [1,2]")
      if (allocated(error)) return
      call check(error, same_row_set(vmfc, n_vmfc, ssfc, n_ssfc), &
                 "SSFC(2) over two monomers must be the same five rows as VMFC(2)")
      if (allocated(error)) return

      call expect_matches(error, vmfc, n_vmfc, N, [.false., .true., .true.], "the VMFC(2) list")
      if (allocated(error)) return
      call expect_matches(error, ssfc, n_ssfc, N, [.false., .true., .true.], "the SSFC(2) list")
      if (allocated(error)) return

      call run_expansion(vmfc, n_vmfc, 2, COUNTERPOISE_VMFC, by_vmfc)
      call run_expansion(ssfc, n_ssfc, 2, COUNTERPOISE_SSFC, by_ssfc)
      expected = made_up([1]) + made_up([2]) + made_up([1, 2]) - made_up([1, -2]) - made_up([-1, 2])
      call check(error, by_vmfc%total_energy, expected, thr=TOLERANCE, &
                 message="VMFC(2) is not the counterpoise-corrected dimer")
      if (allocated(error)) return
      call check(error, by_ssfc%total_energy, by_vmfc%total_energy, thr=TOLERANCE, &
                 message="SSFC(2) and VMFC(2) differ over two monomers")
      if (allocated(error)) return
      call check(error, by_ssfc%dipole(3), by_vmfc%dipole(3), thr=TOLERANCE, &
                 message="SSFC(2) and VMFC(2) dipoles differ over two monomers")
      call by_vmfc%destroy()
      call by_ssfc%destroy()
   end subroutine test_two_monomers

   subroutine test_ssfc_reference(error)
      !! InteractionEnergy under SSFC with screening and a reference fragment
      !!
      !! Four monomers, the reference is monomer 1, and screening keeps the
      !! pairs (1,2) and (1,3) only. Monomer 4 is then in no kept term, so the
      !! reduced list holds full-basis rows for monomers 1 to 3 and none for 4,
      !! yet every row still ghosts all four monomers: the basis is the
      !! cluster's, not the reduced list's. The reference's share is the sum of
      !! its pair terms in the full basis, and its energy is the own-basis one.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: N = 4
      integer :: polymers(8, N)
      integer, allocatable :: damaged(:, :)
      integer(int64), parameter :: COUNT = 8_int64
      type(mbe_result_t) :: mbe_result
      real(dp) :: pair_12, pair_13

      polymers = 0
      ! Own basis, kept monomers only
      polymers(1, 1) = 1
      polymers(2, 1) = 2
      polymers(3, 1) = 3
      ! Full basis, kept monomers only
      polymers(4, :) = [1, -2, -3, -4]
      polymers(5, :) = [2, -1, -3, -4]
      polymers(6, :) = [3, -1, -2, -4]
      ! Kept pairs, the complement of the whole system ghosted
      polymers(7, :) = [1, 2, -3, -4]
      polymers(8, :) = [1, 3, -2, -4]

      call check(error, rows_match_counterpoise(polymers, COUNT, N, COUNTERPOISE_SSFC), &
                 "the reduced SSFC list must match: three monomers twice and two pairs")
      if (allocated(error)) return
      call check(error,.not. rows_match_counterpoise(polymers, COUNT, N, COUNTERPOISE_VMFC), &
                 "the reduced SSFC list must not read as VMFC")
      if (allocated(error)) return

      ! A full-basis row for the monomer no term holds, without its own-basis
      ! row, is a row the list should not have.
      allocate (damaged(COUNT + 1_int64, N))
      damaged(1:COUNT, :) = polymers
      damaged(COUNT + 1_int64, :) = [4, -1, -2, -3]
      call check(error,.not. rows_match_counterpoise(damaged, COUNT + 1_int64, N, COUNTERPOISE_SSFC), &
                 "a full-basis row for a monomer that was dropped was accepted")
      if (allocated(error)) return

      ! A kept monomer without its full-basis row has lost one.
      call check(error,.not. rows_match_counterpoise(drop_row(polymers, COUNT, 6), COUNT - 1_int64, N, &
                                                     COUNTERPOISE_SSFC), &
                 "a list missing a kept monomer's full-basis row was accepted")
      if (allocated(error)) return

      ! Order is not part of the list.
      damaged(1:COUNT, :) = polymers
      damaged(5, :) = polymers(8, :)
      damaged(8, :) = polymers(5, :)
      call check(error, rows_match_counterpoise(damaged, COUNT, N, COUNTERPOISE_SSFC), &
                 "reordering rows must not change the verdict")
      if (allocated(error)) return

      call run_expansion(polymers, COUNT, 2, COUNTERPOISE_SSFC, mbe_result, reference=1)

      pair_12 = made_up([1, 2, -3, -4]) - made_up([1, -2, -3, -4]) - made_up([2, -1, -3, -4])
      pair_13 = made_up([1, 3, -2, -4]) - made_up([1, -2, -3, -4]) - made_up([3, -1, -2, -4])
      call check(error, mbe_result%has_interaction, "a reference run reported no interaction energy")
      if (allocated(error)) return
      call check(error, mbe_result%reference_energy, made_up([1]), thr=TOLERANCE, &
                 message="the reference energy is the own-basis monomer's, not the full-basis one")
      if (allocated(error)) return
      call check(error, mbe_result%interaction_by_level(2), pair_12 + pair_13, thr=TOLERANCE, &
                 message="the reference's share is not its pair terms in the full basis")
      if (allocated(error)) return
      call check(error, mbe_result%interaction_count_by_level(2) == 2_int64, &
                 "two kept terms hold the reference")
      call mbe_result%destroy()
   end subroutine test_ssfc_reference

end module test_mqc_counterpoise_schemes

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_counterpoise_schemes, only: collect_mqc_counterpoise_schemes_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_counterpoise_schemes", collect_mqc_counterpoise_schemes_tests) &
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
