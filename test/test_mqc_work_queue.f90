module test_mqc_work_queue
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_work_queue, only: queue_t, queue_init_from_list, queue_pop, queue_is_empty, queue_destroy
   use pic_types, only: int64
   implicit none
   private
   public :: collect_mqc_work_queue_tests

contains

   subroutine collect_mqc_work_queue_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("empty_queue", test_empty_queue), &
                  new_unittest("queue_order", test_queue_order), &
                  new_unittest("queue_destroy", test_queue_destroy) &
                  ]
   end subroutine collect_mqc_work_queue_tests

   subroutine test_empty_queue(error)
      type(error_type), allocatable, intent(out) :: error
      type(queue_t) :: queue
      integer(int64), allocatable :: ids(:)
      integer(int64) :: item
      logical :: has_item

      allocate (ids(0))
      call queue_init_from_list(queue, ids)

      call check(error, queue_is_empty(queue), "Empty queue should be empty")
      if (allocated(error)) return

      call queue_pop(queue, item, has_item)
      call check(error,.not. has_item, "Pop from empty queue should not return item")
      if (allocated(error)) return

      call check(error, item, -1_int64, "Pop from empty queue should return -1")
      if (allocated(error)) return

      call queue_destroy(queue)
      deallocate (ids)
   end subroutine test_empty_queue

   subroutine test_queue_order(error)
      type(error_type), allocatable, intent(out) :: error
      type(queue_t) :: queue
      integer(int64) :: item
      logical :: has_item
      integer(int64) :: ids(3)

      ids = [3_int64, 1_int64, 5_int64]
      call queue_init_from_list(queue, ids)

      call queue_pop(queue, item, has_item)
      call check(error, has_item, "First pop should return item")
      if (allocated(error)) return
      call check(error, item, 3_int64, "First item should be 3")
      if (allocated(error)) return

      call queue_pop(queue, item, has_item)
      call check(error, has_item, "Second pop should return item")
      if (allocated(error)) return
      call check(error, item, 1_int64, "Second item should be 1")
      if (allocated(error)) return

      call queue_pop(queue, item, has_item)
      call check(error, has_item, "Third pop should return item")
      if (allocated(error)) return
      call check(error, item, 5_int64, "Third item should be 5")
      if (allocated(error)) return

      call check(error, queue_is_empty(queue), "Queue should be empty after pops")
      if (allocated(error)) return

      call queue_pop(queue, item, has_item)
      call check(error,.not. has_item, "Pop after empty should not return item")
      if (allocated(error)) return

      call queue_destroy(queue)
   end subroutine test_queue_order

   subroutine test_queue_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(queue_t) :: queue
      integer(int64) :: ids(1)

      ids = [7_int64]
      call queue_init_from_list(queue, ids)

      call queue_destroy(queue)

      call check(error, queue_is_empty(queue), "Destroyed queue should be empty")
      if (allocated(error)) return
      call check(error, queue%count, 0_int64, "Destroyed queue should reset count to 0")
   end subroutine test_queue_destroy

end module test_mqc_work_queue

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_work_queue, only: collect_mqc_work_queue_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_work_queue", collect_mqc_work_queue_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program tester
