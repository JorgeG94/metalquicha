module test_mqc_json_basis_reader
   !! Cross-checks the Basis Set Exchange JSON reader against the Gaussian94
   !! reader.
   !!
   !! Both parse the same basis from different files written in different
   !! formats by different tooling, so agreement to machine precision is strong
   !! evidence that neither is mangling the data. This matters because
   !! find_basis_file now prefers `.json`, meaning the JSON reader silently
   !! took over every calculation that was previously validated through `.gbs`.
   !!
   !! It also pins the behaviour that motivated the format change: a combined
   !! SP shell must split into separate S and P shells sharing exponents,
   !! which Gaussian94 can only express by being downloaded pre-split.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_json_basis_reader, only: read_json_basis_element
   use mqc_gbs_reader, only: read_gbs_element
   use mqc_cgto, only: atomic_basis_type
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_json_basis_reader_tests

contains

   subroutine collect_mqc_json_basis_reader_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("json_matches_gbs_oxygen", test_matches_gbs_oxygen), &
                  new_unittest("json_matches_gbs_hydrogen", test_matches_gbs_hydrogen), &
                  new_unittest("sp_shell_is_split", test_sp_shell_split), &
                  new_unittest("missing_element_errors", test_missing_element) &
                  ]
   end subroutine collect_mqc_json_basis_reader_tests

   subroutine compare_readers(element, error)
      !! Parse one element both ways and require identical shells
      character(len=*), intent(in) :: element
      type(error_type), allocatable, intent(out) :: error

      type(atomic_basis_type) :: from_json, from_gbs
      type(error_t) :: json_error, gbs_error
      integer :: ishell, iprim
      logical :: have_both

      call read_json_basis_element("../basis_sets/def2-svp.json", element, from_json, json_error)
      call read_gbs_element("../basis_sets/def2-svp.gbs", element, from_gbs, gbs_error)

      ! Both fixtures are checked in, so a missing one is a real failure --
      ! silently skipping would turn this into a test that passes by doing
      ! nothing, which is worse than not having it.
      have_both = .not. json_error%has_error() .and. .not. gbs_error%has_error()
      call check(error, have_both, "both def2-svp fixtures must parse: "// &
                 trim(json_error%get_message())//" "//trim(gbs_error%get_message()))
      if (allocated(error)) return

      call check(error, from_json%nshells > 0, "the comparison must actually see shells")
      if (allocated(error)) return

      call check(error, from_json%nshells, from_gbs%nshells, &
                 "both readers must find the same number of shells")
      if (allocated(error)) return

      do ishell = 1, from_json%nshells
         call check(error, from_json%shells(ishell)%ang_mom, from_gbs%shells(ishell)%ang_mom, &
                    "angular momentum must agree shell by shell")
         if (allocated(error)) return
         call check(error, from_json%shells(ishell)%nfunc, from_gbs%shells(ishell)%nfunc, &
                    "primitive count must agree shell by shell")
         if (allocated(error)) return

         do iprim = 1, from_json%shells(ishell)%nfunc
            call check(error, from_json%shells(ishell)%exponents(iprim), &
                       from_gbs%shells(ishell)%exponents(iprim), &
                       "exponents must agree", thr=1.0e-10_dp)
            if (allocated(error)) return
            call check(error, from_json%shells(ishell)%coefficients(iprim), &
                       from_gbs%shells(ishell)%coefficients(iprim), &
                       "coefficients must agree", thr=1.0e-10_dp)
            if (allocated(error)) return
         end do
      end do
   end subroutine compare_readers

   subroutine test_matches_gbs_oxygen(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_readers("O", error)
   end subroutine test_matches_gbs_oxygen

   subroutine test_matches_gbs_hydrogen(error)
      type(error_type), allocatable, intent(out) :: error
      call compare_readers("H", error)
   end subroutine test_matches_gbs_hydrogen

   subroutine test_sp_shell_split(error)
      !! STO-3G oxygen is one S shell plus one SP shell in the JSON. The reader
      !! must present that as three shells -- S, S, P -- with the S and P of
      !! the combined entry sharing exponents.
      type(error_type), allocatable, intent(out) :: error
      type(atomic_basis_type) :: basis
      type(error_t) :: parse_error
      integer :: iprim

      call read_json_basis_element("../basis_sets/sto-3g.json", "O", basis, parse_error)
      call check(error,.not. parse_error%has_error(), "sto-3g.json must parse: "// &
                 trim(parse_error%get_message()))
      if (allocated(error)) return

      call check(error, basis%nshells, 3, "STO-3G oxygen must split into three shells")
      if (allocated(error)) return
      call check(error, basis%shells(1)%ang_mom, 0, "shell 1 is S")
      if (allocated(error)) return
      call check(error, basis%shells(2)%ang_mom, 0, "shell 2 is the S of the SP pair")
      if (allocated(error)) return
      call check(error, basis%shells(3)%ang_mom, 1, "shell 3 is the P of the SP pair")
      if (allocated(error)) return

      do iprim = 1, basis%shells(2)%nfunc
         call check(error, basis%shells(2)%exponents(iprim), basis%shells(3)%exponents(iprim), &
                    "a split SP shell must share its exponents", thr=1.0e-12_dp)
         if (allocated(error)) return
      end do
   end subroutine test_sp_shell_split

   subroutine test_missing_element(error)
      !! An element absent from the file must be an error, not an empty basis
      type(error_type), allocatable, intent(out) :: error
      type(atomic_basis_type) :: basis
      type(error_t) :: parse_error

      call read_json_basis_element("../basis_sets/def2-svp.json", "Fe", basis, parse_error)
      call check(error, parse_error%has_error(), "a missing element must report an error")
   end subroutine test_missing_element

end module test_mqc_json_basis_reader

program tester_mqc_json_basis_reader
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_json_basis_reader, only: collect_mqc_json_basis_reader_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_json_basis_reader", collect_mqc_json_basis_reader_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_json_basis_reader
