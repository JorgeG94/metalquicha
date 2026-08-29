module test_mqc_ecp_refusals
   !! Tests the two guards that keep an ECP out of code that cannot carry it.
   !!
   !! Neither guard protects against a crash. Both protect against a number: a
   !! nuclear derivative missing the potential's contribution, and a frozen core
   !! counted as though the potential had not already removed one. In each case
   !! the calculation converges, the output looks ordinary, and the answer is
   !! wrong -- so what is being tested here is that the refusal happens at all,
   !! and that it says enough for the user to act on.
   !!
   !! These take the per-atom core-electron count rather than a molecule, which
   !! is what makes them testable without building an integral environment: an
   !! array of zeros is an all-electron system and a non-zero entry anywhere is
   !! an ECP one. That is exactly the distinction the callers make.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_libcint_ecp, only: ecp_refuses_derivatives, ecp_refuses_auto_frozen_core
   use mqc_error, only: error_t
   implicit none
   private
   public :: collect_mqc_ecp_refusals_tests

contains

   subroutine collect_mqc_ecp_refusals_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("all_electron_gradient_is_allowed", test_grad_all_electron), &
                  new_unittest("ecp_gradient_is_refused", test_grad_refused), &
                  new_unittest("one_heavy_atom_is_enough", test_grad_mixed), &
                  new_unittest("gradient_message_names_the_derivative", test_grad_names_what), &
                  new_unittest("unallocated_count_is_all_electron", test_grad_unallocated), &
                  new_unittest("all_electron_frozen_core_is_allowed", test_frozen_all_electron), &
                  new_unittest("ecp_auto_frozen_core_is_refused", test_frozen_refused), &
                  new_unittest("frozen_core_message_offers_the_way_out", test_frozen_message) &
                  ]
   end subroutine collect_mqc_ecp_refusals_tests

   ! ---- nuclear derivatives -------------------------------------------------

   subroutine test_grad_all_electron(error)
      !! An all-electron system must pass through untouched. A guard that
      !! refused everything would be caught by the whole gradient suite, but
      !! only after it had been built -- this states it directly.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer, allocatable :: core(:)

      core = [0, 0, 0]
      call check(error,.not. ecp_refuses_derivatives(core, "nuclear gradient", err), &
                 "an all-electron gradient should not be refused")
      if (allocated(error)) return
      call check(error,.not. err%has_error(), &
                 "an all-electron gradient should not set an error")
   end subroutine test_grad_all_electron

   subroutine test_grad_refused(error)
      !! Iodine's def2-ECP: 28 electrons replaced, so the derivative is refused.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer, allocatable :: core(:)

      core = [28, 0]
      call check(error, ecp_refuses_derivatives(core, "nuclear gradient", err), &
                 "a gradient over an ECP should be refused")
      if (allocated(error)) return
      call check(error, err%has_error(), "the refusal should set an error")
   end subroutine test_grad_refused

   subroutine test_grad_mixed(error)
      !! The count is per atom, and one heavy atom among light ones is the
      !! ordinary case -- HgCl2 has a potential on one centre of three. Testing
      !! `any` rather than `all` matters: the opposite would let every mixed
      !! molecule through, which is most of them.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer, allocatable :: core(:)

      core = [0, 0, 60, 0, 0]
      call check(error, ecp_refuses_derivatives(core, "nuclear gradient", err), &
                 "one atom carrying a potential should refuse the whole gradient")
   end subroutine test_grad_mixed

   subroutine test_grad_names_what(error)
      !! Four derivative paths share this function, so the message has to say
      !! which one refused or the user cannot tell an MP2 gradient from a
      !! Hessian in a log.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer, allocatable :: core(:)
      character(len=:), allocatable :: message
      logical :: refused

      core = [28]
      refused = ecp_refuses_derivatives(core, "MP2 gradient", err)
      call check(error, refused, "should be refused")
      if (allocated(error)) return
      message = err%get_message()
      call check(error, index(message, "MP2 gradient") > 0, &
                 "the message should name the derivative but said: "//message)
   end subroutine test_grad_names_what

   subroutine test_grad_unallocated(error)
      !! An unallocated count is an all-electron run rather than a bug: the
      !! array is only filled in when a potential was read.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer, allocatable :: core(:)

      call check(error,.not. ecp_refuses_derivatives(core, "nuclear gradient", err), &
                 "an unallocated core count should not refuse")
   end subroutine test_grad_unallocated

   ! ---- automatically counted frozen cores ----------------------------------

   subroutine test_frozen_all_electron(error)
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer, allocatable :: core(:)

      core = [0, 0, 0]
      call check(error,.not. ecp_refuses_auto_frozen_core(core, err), &
                 "an all-electron frozen core should not be refused")
   end subroutine test_frozen_all_electron

   subroutine test_frozen_refused(error)
      !! The count would come from the atomic number, and the potential has
      !! already taken the orbitals it would freeze -- so it would freeze
      !! valence instead, by an amount that reads as a small energy difference.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer, allocatable :: core(:)

      core = [28, 0]
      call check(error, ecp_refuses_auto_frozen_core(core, err), &
                 "an automatically counted frozen core over an ECP should be refused")
      if (allocated(error)) return
      call check(error, err%has_error(), "the refusal should set an error")
   end subroutine test_frozen_refused

   subroutine test_frozen_message(error)
      !! This one is recoverable -- the user sets the count themselves -- so the
      !! message has to name the keyword that recovers it. A refusal with no way
      !! forward is why the automatic count was tempting in the first place.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer, allocatable :: core(:)
      character(len=:), allocatable :: message
      logical :: refused

      core = [60]
      refused = ecp_refuses_auto_frozen_core(core, err)
      call check(error, refused, "should be refused")
      if (allocated(error)) return
      message = err%get_message()
      call check(error, index(message, "n_frozen_core") > 0, &
                 "the message should name the keyword that fixes it but said: "//message)
   end subroutine test_frozen_message

end module test_mqc_ecp_refusals

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_ecp_refusals, only: collect_mqc_ecp_refusals_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_ecp_refusals", collect_mqc_ecp_refusals_tests) &
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
