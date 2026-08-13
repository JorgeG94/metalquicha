module test_mqc_capi_fraglist
   !! Tests the C entry points for the fragment term list.
   !!
   !! `mqc_fraglist` itself is exercised by `test_mqc_fraglist`; this suite is
   !! about the boundary the C wrapper adds -- the opaque handle, the two-call
   !! read-out (ask the count, then fill a buffer), the row-major packing a
   !! Python caller indexes, and the status codes that stand in for an
   !! `error_t` nobody on the far side could read. The round trip that matters
   !! is set-then-get: what a screen hands back must come out in the order and
   !! shape it went in.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_int64_t, &
                                                                             c_char, c_null_char
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_capi_fraglist, only: mqc_fraglist_new, mqc_fraglist_free, &
                                mqc_fraglist_generate, mqc_fraglist_count, &
                                mqc_fraglist_max_level, mqc_fraglist_get, &
                                mqc_fraglist_set, mqc_fraglist_close_subsets, &
                                mqc_fraglist_last_error
   implicit none
   private
   public :: collect_mqc_capi_fraglist_tests

   integer(c_int), parameter :: MQC_OK = 0
   integer(c_int), parameter :: MQC_FAIL = 1
   integer(c_int), parameter :: MQC_BAD_HANDLE = 2

contains

   subroutine collect_mqc_capi_fraglist_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("new_list_is_empty", test_new_empty), &
                  new_unittest("generate_counts_combinations", test_generate), &
                  new_unittest("get_reads_the_terms_out", test_get), &
                  new_unittest("get_refuses_a_small_buffer", test_get_small_buffer), &
                  new_unittest("set_then_get_round_trips", test_round_trip), &
                  new_unittest("close_subsets_grows_the_list", test_close_subsets), &
                  new_unittest("bad_handle_everywhere", test_bad_handle), &
                  new_unittest("generate_rejects_bad_level", test_generate_bad_level), &
                  new_unittest("set_rejects_bad_shape", test_set_bad_shape), &
                  new_unittest("last_error_is_retrievable", test_last_error) &
                  ]
   end subroutine collect_mqc_capi_fraglist_tests

   subroutine test_new_empty(error)
      !! A fresh list holds nothing
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle

      handle = mqc_fraglist_new()
      call check(error, int(mqc_fraglist_count(handle)), 0, "a new list has no terms")

      call mqc_fraglist_free(handle)
   end subroutine test_new_empty

   subroutine test_generate(error)
      !! Five monomers at level 2 is C(5,2) = 10 pairs, rows of width 2
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status

      handle = mqc_fraglist_new()
      status = mqc_fraglist_generate(handle, 5_c_int, 2_c_int)
      call check(error, int(status), int(MQC_OK), "generation succeeds")
      if (allocated(error)) return
      call check(error, int(mqc_fraglist_count(handle)), 10, "C(5,2) = 10 dimers")
      if (allocated(error)) return
      call check(error, int(mqc_fraglist_max_level(handle)), 2, "rows are two wide")

      call mqc_fraglist_free(handle)
   end subroutine test_generate

   subroutine test_get(error)
      !! Three monomers give exactly the pairs {1,2},{1,3},{2,3}, in some order
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: terms(6)
      logical :: seen_12, seen_13, seen_23
      integer :: i, a, b

      handle = mqc_fraglist_new()
      status = mqc_fraglist_generate(handle, 3_c_int, 2_c_int)
      call check(error, int(mqc_fraglist_count(handle)), 3, "C(3,2) = 3 pairs")
      if (allocated(error)) return

      terms = -1
      status = mqc_fraglist_get(handle, terms, 3_c_int64_t, 2_c_int)
      call check(error, int(status), int(MQC_OK), "read-out succeeds")
      if (allocated(error)) return

      seen_12 = .false.; seen_13 = .false.; seen_23 = .false.
      do i = 1, 3
         a = terms(2*(i - 1) + 1)
         b = terms(2*(i - 1) + 2)
         if (min(a, b) == 1 .and. max(a, b) == 2) seen_12 = .true.
         if (min(a, b) == 1 .and. max(a, b) == 3) seen_13 = .true.
         if (min(a, b) == 2 .and. max(a, b) == 3) seen_23 = .true.
      end do
      call check(error, seen_12 .and. seen_13 .and. seen_23, &
                 "every pair of three monomers is present exactly once")

      call mqc_fraglist_free(handle)
   end subroutine test_get

   subroutine test_get_small_buffer(error)
      !! A buffer sized from a stale count is refused, not overrun
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: terms(10)   ! room for 5 rows, but the list has 10

      handle = mqc_fraglist_new()
      status = mqc_fraglist_generate(handle, 5_c_int, 2_c_int)
      status = mqc_fraglist_get(handle, terms, 5_c_int64_t, 2_c_int)
      call check(error, int(status), int(MQC_FAIL), "a short buffer is refused")

      call mqc_fraglist_free(handle)
   end subroutine test_get_small_buffer

   subroutine test_round_trip(error)
      !! What a caller sets is what a caller gets, in order
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: given(6), got(6)

      ! Three dimers of a caller's own choosing and order, row-major.
      given = [4, 7, 1, 2, 3, 9]

      handle = mqc_fraglist_new()
      status = mqc_fraglist_set(handle, given, 3_c_int64_t, 2_c_int)
      call check(error, int(status), int(MQC_OK), "set succeeds")
      if (allocated(error)) return
      call check(error, int(mqc_fraglist_count(handle)), 3, "three terms went in")
      if (allocated(error)) return
      call check(error, int(mqc_fraglist_max_level(handle)), 2, "width preserved")
      if (allocated(error)) return

      got = -1
      status = mqc_fraglist_get(handle, got, 3_c_int64_t, 2_c_int)
      call check(error, int(status), int(MQC_OK), "get succeeds")
      if (allocated(error)) return
      call check(error, all(got == given), &
                 "the caller's terms come back unchanged and unreordered")

      call mqc_fraglist_free(handle)
   end subroutine test_round_trip

   subroutine test_close_subsets(error)
      !! A lone trimer needs its three dimers and three monomers back: seven
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: given(3)

      given = [1, 2, 3]

      handle = mqc_fraglist_new()
      status = mqc_fraglist_set(handle, given, 1_c_int64_t, 3_c_int)
      call check(error, int(status), int(MQC_OK), "set the trimer")
      if (allocated(error)) return

      status = mqc_fraglist_close_subsets(handle)
      call check(error, int(status), int(MQC_OK), "closure succeeds")
      if (allocated(error)) return
      call check(error, int(mqc_fraglist_count(handle)), 7, &
                 "3 monomers + 3 dimers + the trimer itself")

      call mqc_fraglist_free(handle)
   end subroutine test_close_subsets

   subroutine test_bad_handle(error)
      !! Every entry point copes with a null handle
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: status
      integer(c_int) :: terms(2)

      call check(error, int(mqc_fraglist_count(c_null_ptr)), -1, "null -> -1 count")
      if (allocated(error)) return
      call check(error, int(mqc_fraglist_max_level(c_null_ptr)), -1, "null -> -1 level")
      if (allocated(error)) return

      status = mqc_fraglist_generate(c_null_ptr, 4_c_int, 2_c_int)
      call check(error, int(status), int(MQC_BAD_HANDLE), "generate -> bad handle")
      if (allocated(error)) return

      terms = 0
      status = mqc_fraglist_get(c_null_ptr, terms, 1_c_int64_t, 2_c_int)
      call check(error, int(status), int(MQC_BAD_HANDLE), "get -> bad handle")
      if (allocated(error)) return
      status = mqc_fraglist_set(c_null_ptr, terms, 1_c_int64_t, 2_c_int)
      call check(error, int(status), int(MQC_BAD_HANDLE), "set -> bad handle")
      if (allocated(error)) return
      status = mqc_fraglist_close_subsets(c_null_ptr)
      call check(error, int(status), int(MQC_BAD_HANDLE), "close_subsets -> bad handle")
   end subroutine test_bad_handle

   subroutine test_generate_bad_level(error)
      !! Level 1 has no n-mer terms and is rejected, not silently emptied
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status

      handle = mqc_fraglist_new()
      status = mqc_fraglist_generate(handle, 5_c_int, 1_c_int)
      call check(error, int(status), int(MQC_FAIL), "level 1 is rejected")

      call mqc_fraglist_free(handle)
   end subroutine test_generate_bad_level

   subroutine test_set_bad_shape(error)
      !! A non-positive width cannot describe a term
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: terms(2)

      handle = mqc_fraglist_new()
      terms = 0
      status = mqc_fraglist_set(handle, terms, 1_c_int64_t, 0_c_int)
      call check(error, int(status), int(MQC_FAIL), "max_level below 1 is refused")

      call mqc_fraglist_free(handle)
   end subroutine test_set_bad_shape

   subroutine test_last_error(error)
      !! A failed call leaves a message the C side can read back
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      character(len=:), allocatable :: message

      handle = mqc_fraglist_new()
      status = mqc_fraglist_generate(handle, 5_c_int, 1_c_int)
      call check(error, int(status), int(MQC_FAIL), "provoke a failure")
      if (allocated(error)) return

      message = last_fraglist_error()
      call check(error, len_trim(message) > 0, "the failure left a message behind")
      if (allocated(error)) return
      call check(error, index(message, "max_level") > 0, &
                 "the message names the actual problem, got: "//message)

      call mqc_fraglist_free(handle)
   end subroutine test_last_error

   function last_fraglist_error() result(message)
      !! Pull the most recent error out through the C string boundary
      !!
      !! Note the argument order here is (buffer, buffer_len), the reverse of
      !! the system and session variants -- a wart the wrapper carries, and a
      !! reason to call it through a helper rather than inline.
      character(len=:), allocatable :: message
      integer, parameter :: BUFLEN = 512
      character(kind=c_char) :: buffer(BUFLEN)
      integer :: i

      buffer = c_null_char
      call mqc_fraglist_last_error(buffer, int(BUFLEN, c_int))
      message = ""
      do i = 1, BUFLEN
         if (buffer(i) == c_null_char) exit
         message = message//buffer(i)
      end do
   end function last_fraglist_error

end module test_mqc_capi_fraglist

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_capi_fraglist, only: collect_mqc_capi_fraglist_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_capi_fraglist", collect_mqc_capi_fraglist_tests) &
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
