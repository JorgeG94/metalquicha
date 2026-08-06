module test_mqc_basis_utils
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_basis_utils, only: normalize_basis_name, find_basis_file, basis_search_path
   use mqc_error, only: error_t
   implicit none
   private
   public :: collect_mqc_basis_utils_tests

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_basis_utils_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("normalize_basis_stars", test_normalize_stars), &
                  new_unittest("normalize_basis_plus", test_normalize_plus), &
                  new_unittest("normalize_basis_parens", test_normalize_parens), &
                  new_unittest("normalize_basis_complex", test_normalize_complex), &
                  new_unittest("normalize_basis_case", test_normalize_case), &
                  new_unittest("parenthesized_names_stay_distinct", &
                               test_parenthesized_names_stay_distinct), &
                  new_unittest("find_basis_file_exists", test_find_basis_exists), &
                  new_unittest("find_basis_file_not_found", test_find_basis_not_found), &
                  new_unittest("basis_search_path_order", test_search_path_order) &
                  ]
   end subroutine collect_mqc_basis_utils_tests

   subroutine test_normalize_stars(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: result

      ! Single star
      result = normalize_basis_name("6-31G*")
      call check(error, trim(result), "6-31g_st_", "6-31G* should become 6-31g_st_")
      if (allocated(error)) return

      ! Double star
      result = normalize_basis_name("6-31G**")
      call check(error, trim(result), "6-31g_st__st_", "6-31G** should become 6-31g_st__st_")
   end subroutine test_normalize_stars

   subroutine test_normalize_plus(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: result

      ! Single plus
      result = normalize_basis_name("6-31+G")
      call check(error, trim(result), "6-31+g", "6-31+G keeps its plus")
      if (allocated(error)) return

      ! Double plus
      result = normalize_basis_name("6-311++G")
      call check(error, trim(result), "6-311++g", "6-311++G keeps both pluses")
   end subroutine test_normalize_plus

   subroutine test_normalize_parens(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: result

      ! Single parentheses
      result = normalize_basis_name("6-31G(d)")
      call check(error, trim(result), "6-31g(d)", "6-31G(d) keeps its parentheses")
      if (allocated(error)) return

      ! Multiple in parentheses
      result = normalize_basis_name("6-311G(d,p)")
      call check(error, trim(result), "6-311g(d,p)", "6-311G(d,p) keeps parentheses and comma")
      if (allocated(error)) return

      ! Triple in parentheses
      result = normalize_basis_name("6-311G(2d,2p)")
      call check(error, trim(result), "6-311g(2d,2p)", "6-311G(2d,2p) keeps parentheses and comma")
   end subroutine test_normalize_parens

   subroutine test_normalize_complex(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: result

      ! Plus and star combined
      result = normalize_basis_name("6-31+G*")
      call check(error, trim(result), "6-31+g_st_", "6-31+G* should become 6-31+g_st_")
      if (allocated(error)) return

      ! Double plus and double star
      result = normalize_basis_name("6-311++G**")
      call check(error, trim(result), "6-311++g_st__st_", "6-311++G** should become 6-311++g_st__st_")
      if (allocated(error)) return

      ! Plus with parentheses
      result = normalize_basis_name("6-31+G(d)")
      call check(error, trim(result), "6-31+g(d)", "6-31+G(d) should become 6-31+g(d)")
   end subroutine test_normalize_complex

   subroutine test_normalize_case(error)
      !! Names with no punctuation are only lowercased
      !!
      !! Lowercasing is not cosmetic: the Basis Set Exchange bundle the files
      !! are extracted from is entirely lowercase, so "cc-pVDZ" and "cc-pvdz"
      !! have to land on the same file or half the validation inputs stop
      !! resolving.
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: result

      result = normalize_basis_name("6-31G")
      call check(error, trim(result), "6-31g", "6-31G should become 6-31g")
      if (allocated(error)) return

      result = normalize_basis_name("STO-3G")
      call check(error, trim(result), "sto-3g", "STO-3G should become sto-3g")
      if (allocated(error)) return

      result = normalize_basis_name("cc-pVDZ")
      call check(error, trim(result), "cc-pvdz", "cc-pVDZ should become cc-pvdz")
      if (allocated(error)) return

      ! Already lowercase, and the spelling most inputs use
      result = normalize_basis_name("def2-universal-jkfit")
      call check(error, trim(result), "def2-universal-jkfit", &
                 "an already-normalized name should be unchanged")
   end subroutine test_normalize_case

   subroutine test_parenthesized_names_stay_distinct(error)
      !! def2-SV(P) and def2-SVP must not collapse onto one file
      !!
      !! They are different basis sets -- SV(P) polarizes only the heavy atoms
      !! -- and an earlier normalization dropped parentheses, so both resolved
      !! to def2-svp.json. Whichever was extracted last won, and every
      !! calculation after that ran in a basis nobody asked for. Same story for
      !! their RIFIT auxiliaries.
      type(error_type), allocatable, intent(out) :: error

      call check(error, normalize_basis_name("def2-SV(P)"), "def2-sv(p)", &
                 "def2-SV(P) should keep its parentheses")
      if (allocated(error)) return

      call check(error, normalize_basis_name("def2-SVP"), "def2-svp", &
                 "def2-SVP should normalize to def2-svp")
      if (allocated(error)) return

      call check(error, normalize_basis_name("def2-SV(P)") /= &
                 normalize_basis_name("def2-SVP"), &
                 "def2-SV(P) and def2-SVP must resolve to different files")
      if (allocated(error)) return

      call check(error, normalize_basis_name("def2-SV(P)-RIFIT") /= &
                 normalize_basis_name("def2-SVP-RIFIT"), &
                 "the SV(P) and SVP RIFIT auxiliaries must stay distinct too")
   end subroutine test_parenthesized_names_stay_distinct

   subroutine test_find_basis_exists(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: filename
      type(error_t) :: parse_error
      integer :: stat
      logical :: file_exists

      ! Create basis_sets directory if it doesn't exist
      call execute_command_line("mkdir -p basis_sets", exitstat=stat)

      ! Create a test basis file in basis_sets directory. Mixed case on the
      ! way in, lowercase on disk -- that is the mapping being tested.
      call create_test_basis_file("basis_sets/test_basis.json")

      call find_basis_file("Test_Basis", filename, parse_error)

      if (parse_error%has_error()) then
         call check(error, .false., "Should find test_basis.json: "//parse_error%get_message())
         call delete_file("basis_sets/test_basis.json")
         return
      end if

      ! Check that the returned filename actually exists
      inquire (file=trim(filename), exist=file_exists)
      call check(error, file_exists, &
                 "Returned filename should exist: "//trim(filename))
      if (allocated(error)) then
         call delete_file("basis_sets/test_basis.json")
         return
      end if

      ! Check that the filename is correct
      call check(error, trim(filename), "basis_sets/test_basis.json", &
                 "Filename should be basis_sets/test_basis.json")

      call delete_file("basis_sets/test_basis.json")
   end subroutine test_find_basis_exists

   subroutine test_search_path_order(error)
      !! ./basis_sets must be searched, and searched before the built-in default
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: directories(:)

      directories = basis_search_path()

      call check(error, size(directories) >= 1, "the search path is never empty")
      if (allocated(error)) return

      ! With MQC_BASIS_PATH unset -- as it is under ctest -- the relative
      ! directory leads, which is what the validation scripts depend on.
      call check(error, trim(directories(1)), "basis_sets", &
                 "./basis_sets should lead when MQC_BASIS_PATH is unset")
   end subroutine test_search_path_order

   subroutine test_find_basis_not_found(error)
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: filename
      type(error_t) :: parse_error

      ! Try to find non-existent basis
      call find_basis_file("NONEXISTENT-BASIS", filename, parse_error)

      call check(error, parse_error%has_error(), "Should fail to find non-existent basis")
   end subroutine test_find_basis_not_found

   ! Helper subroutines

   subroutine create_test_basis_file(filename)
      character(len=*), intent(in) :: filename
      integer :: unit
      open (newunit=unit, file=filename, status="replace", action="write")
      write (unit, "(a)") '{"elements": {}}'
      close (unit)
   end subroutine create_test_basis_file

   subroutine delete_file(filename)
      character(len=*), intent(in) :: filename
      integer :: unit, stat
      open (newunit=unit, file=filename, status="old", action="read", iostat=stat)
      if (stat == 0) close (unit, status="delete")
   end subroutine delete_file

end module test_mqc_basis_utils

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_basis_utils, only: collect_mqc_basis_utils_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_basis_utils", collect_mqc_basis_utils_tests) &
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
