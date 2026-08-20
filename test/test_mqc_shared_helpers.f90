!! Unit tests for the small helpers that were deduplicated out of the backends
module test_mqc_shared_helpers
   !! `hund_multiplicity`, `int_to_text` and `ends_with` each existed in two or
   !! more copies before being hoisted. These tests pin the behaviour the
   !! callers depended on, including the case-sensitivity change `ends_with`
   !! deliberately made.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_atomic_guess_common, only: hund_multiplicity
   use mqc_string_utils, only: int_to_text
   use mqc_io_helpers, only: ends_with
   implicit none
   private

   public :: collect_mqc_shared_helpers

contains

   subroutine collect_mqc_shared_helpers(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("hund_closed_shells_are_singlets", test_hund_closed), &
                  new_unittest("hund_first_row_follows_the_periodic_table", test_hund_first_row), &
                  new_unittest("hund_is_madelung_not_experiment_for_cr_and_cu", test_hund_madelung), &
                  new_unittest("hund_multiplicity_is_always_at_least_one", test_hund_positive), &
                  new_unittest("int_to_text_has_no_padding", test_int_to_text), &
                  new_unittest("ends_with_ignores_case", test_ends_with_case), &
                  new_unittest("ends_with_rejects_a_longer_suffix", test_ends_with_edges) &
                  ]
   end subroutine collect_mqc_shared_helpers

   subroutine test_hund_closed(error)
      !! Noble gases and the alkaline earths close a subshell exactly
      type(error_type), allocatable, intent(out) :: error

      call check(error, hund_multiplicity(2), 1)    ! He, 1s2
      if (allocated(error)) return
      call check(error, hund_multiplicity(4), 1)    ! Be, 2s2
      if (allocated(error)) return
      call check(error, hund_multiplicity(10), 1)   ! Ne, 2p6
      if (allocated(error)) return
      call check(error, hund_multiplicity(18), 1)   ! Ar, 3p6
      if (allocated(error)) return
      call check(error, hund_multiplicity(36), 1)   ! Kr, 4p6
   end subroutine test_hund_closed

   subroutine test_hund_first_row(error)
      !! The p block fills singly first, so the multiplicity peaks at nitrogen
      type(error_type), allocatable, intent(out) :: error

      call check(error, hund_multiplicity(1), 2)    ! H,  1s1
      if (allocated(error)) return
      call check(error, hund_multiplicity(5), 2)    ! B,  2p1
      if (allocated(error)) return
      call check(error, hund_multiplicity(6), 3)    ! C,  2p2
      if (allocated(error)) return
      call check(error, hund_multiplicity(7), 4)    ! N,  2p3, half filled
      if (allocated(error)) return
      call check(error, hund_multiplicity(8), 3)    ! O,  2p4
      if (allocated(error)) return
      call check(error, hund_multiplicity(9), 2)    ! F,  2p5
   end subroutine test_hund_first_row

   subroutine test_hund_madelung(error)
      !! Chromium and copper are the documented Madelung failures
      !!
      !! Experimentally Cr is 3d5 4s1 (septet) and Cu is 3d10 4s1 (doublet).
      !! Strict Madelung filling gives 3d4 4s2 and 3d9 4s2 instead. The guess is
      !! allowed to be wrong here -- it costs iterations, not correctness -- and
      !! this test records that it *is* wrong, so the behaviour is deliberate
      !! rather than discovered later as a surprise.
      type(error_type), allocatable, intent(out) :: error

      call check(error, hund_multiplicity(24), 5)   ! Cr as 3d4, not the 7 observed
      if (allocated(error)) return
      call check(error, hund_multiplicity(29), 2)   ! Cu as 3d9, which matches by luck
   end subroutine test_hund_madelung

   subroutine test_hund_positive(error)
      !! No element may come back with a nonsensical multiplicity
      type(error_type), allocatable, intent(out) :: error
      integer :: z

      do z = 1, 118
         call check(error, hund_multiplicity(z) >= 1, "multiplicity must be at least one")
         if (allocated(error)) return
      end do
      ! A ghost centre has no electrons and so is a singlet.
      call check(error, hund_multiplicity(0), 1)
   end subroutine test_hund_positive

   subroutine test_int_to_text(error)
      !! The callers concatenate this straight into paths and messages
      type(error_type), allocatable, intent(out) :: error

      call check(error, int_to_text(0), "0")
      if (allocated(error)) return
      call check(error, int_to_text(7), "7")
      if (allocated(error)) return
      call check(error, int_to_text(118), "118")
      if (allocated(error)) return
      call check(error, int_to_text(-42), "-42")
      if (allocated(error)) return
      ! No leading or trailing blanks, which is what makes it safe to build a
      ! JSON path like "elements.6.electron_shells(2)" by concatenation.
      call check(error, len(int_to_text(6)), 1)
   end subroutine test_int_to_text

   subroutine test_ends_with_case(error)
      !! Insensitive, so an uppercased extension is still recognised
      !!
      !! This is the behaviour change the merge made: the copy in the deck
      !! reader used to be case-sensitive, which refused `WATER.JSON`.
      type(error_type), allocatable, intent(out) :: error

      call check(error, ends_with("water.json", ".json"), "the plain case must match")
      if (allocated(error)) return
      call check(error, ends_with("WATER.JSON", ".json"), "an uppercase name must match")
      if (allocated(error)) return
      call check(error, ends_with("water.JsOn", ".json"), "a mixed-case name must match")
      if (allocated(error)) return
      call check(error, ends_with("run.h5", ".H5"), "an uppercase suffix must match")
      if (allocated(error)) return
      call check(error,.not. ends_with("water.xyz", ".json"), "a different extension must not match")
   end subroutine test_ends_with_case

   subroutine test_ends_with_edges(error)
      !! A suffix longer than the text cannot match, and must not read past it
      type(error_type), allocatable, intent(out) :: error

      call check(error,.not. ends_with("a", ".json"), "a longer suffix cannot match")
      if (allocated(error)) return
      call check(error,.not. ends_with("", ".json"), "an empty name cannot match")
      if (allocated(error)) return
      call check(error,.not. ends_with("water.json", ""), "an empty suffix must not match")
      if (allocated(error)) return
      call check(error, ends_with("water.json", "water.json"), "the whole string is a valid suffix")
   end subroutine test_ends_with_edges

end module test_mqc_shared_helpers

program tester_mqc_shared_helpers
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_shared_helpers, only: collect_mqc_shared_helpers
   implicit none

   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_shared_helpers", collect_mqc_shared_helpers)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_shared_helpers
