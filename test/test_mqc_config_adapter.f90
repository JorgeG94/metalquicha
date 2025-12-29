module test_mqc_config_adapter
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_config_parser, only: input_fragment_t
   use mqc_config_adapter, only: check_fragment_overlap
   use mqc_error, only: error_t
   implicit none
   private
   public :: collect_mqc_config_adapter_tests

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_config_adapter_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("no_overlap_detection", test_no_overlap), &
                  new_unittest("overlap_detection", test_overlap_detected), &
                  new_unittest("single_fragment_no_overlap", test_single_fragment) &
                  ]
   end subroutine collect_mqc_config_adapter_tests

   subroutine test_no_overlap(error)
      !! Test that non-overlapping fragments pass validation
      type(error_type), allocatable, intent(out) :: error
      type(input_fragment_t), allocatable :: fragments(:)
      type(error_t) :: parse_error

      ! Create two non-overlapping fragments
      allocate (fragments(2))

      ! Fragment 1: atoms 0, 1, 2
      allocate (fragments(1)%indices(3))
      fragments(1)%indices = [0, 1, 2]
      fragments(1)%charge = 0
      fragments(1)%multiplicity = 1

      ! Fragment 2: atoms 3, 4, 5
      allocate (fragments(2)%indices(3))
      fragments(2)%indices = [3, 4, 5]
      fragments(2)%charge = 0
      fragments(2)%multiplicity = 1

      ! Check for overlap - should find none
      call check_fragment_overlap(fragments, 2, parse_error)

      call check(error,.not. parse_error%has_error(), "Non-overlapping fragments should pass validation")
      if (allocated(error)) return

      ! Clean up
      call fragments(1)%destroy()
      call fragments(2)%destroy()
      deallocate (fragments)
   end subroutine test_no_overlap

   subroutine test_overlap_detected(error)
      !! Test that overlapping fragments are detected
      type(error_type), allocatable, intent(out) :: error
      type(input_fragment_t), allocatable :: fragments(:)
      type(error_t) :: parse_error

      ! Create two overlapping fragments
      allocate (fragments(2))

      ! Fragment 1: atoms 0, 1, 2
      allocate (fragments(1)%indices(3))
      fragments(1)%indices = [0, 1, 2]
      fragments(1)%charge = 0
      fragments(1)%multiplicity = 1

      ! Fragment 2: atoms 2, 3, 4 (overlaps on atom 2)
      allocate (fragments(2)%indices(3))
      fragments(2)%indices = [2, 3, 4]
      fragments(2)%charge = 0
      fragments(2)%multiplicity = 1

      ! Check for overlap - should detect it
      call check_fragment_overlap(fragments, 2, parse_error)

      call check(error, parse_error%has_error(), "Overlapping fragments should be detected")
      if (allocated(error)) return

      ! Check that error message mentions the overlapping atom
      call check(error, index(parse_error%get_message(), "atom 2") > 0, &
                 "Error message should mention overlapping atom")
      if (allocated(error)) return

      ! Clean up
      call fragments(1)%destroy()
      call fragments(2)%destroy()
      deallocate (fragments)
   end subroutine test_overlap_detected

   subroutine test_single_fragment(error)
      !! Test that a single fragment has no overlap (edge case)
      type(error_type), allocatable, intent(out) :: error
      type(input_fragment_t), allocatable :: fragments(:)
      type(error_t) :: parse_error

      ! Create a single fragment
      allocate (fragments(1))

      ! Fragment 1: atoms 0, 1, 2
      allocate (fragments(1)%indices(3))
      fragments(1)%indices = [0, 1, 2]
      fragments(1)%charge = 0
      fragments(1)%multiplicity = 1

      ! Check for overlap - single fragment should pass
      call check_fragment_overlap(fragments, 1, parse_error)

      call check(error,.not. parse_error%has_error(), "Single fragment should pass validation")
      if (allocated(error)) return

      ! Clean up
      call fragments(1)%destroy()
      deallocate (fragments)
   end subroutine test_single_fragment

end module test_mqc_config_adapter

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_config_adapter, only: collect_mqc_config_adapter_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_config_adapter", collect_mqc_config_adapter_tests) &
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
