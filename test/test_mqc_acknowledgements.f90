module test_mqc_acknowledgements
   !! That the backend a run credits is the backend it will actually use
   !!
   !! The banner is a claim about provenance, so the thing worth testing is not
   !! that it prints -- it is that `acknowledged_backend` resolves the same way
   !! the dispatch in `mqc_method_hf`/`mqc_method_dft` does. Thanking libcint
   !! for a run that happened on the GPU would be a lie in the output the user
   !! keeps, and the two resolutions can only drift apart silently.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_acknowledgements, only: acknowledged_backend, &
                                   ACK_CUEST, ACK_LIBCINT, ACK_TBLITE
   use mqc_method_types, only: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2, METHOD_TYPE_HF, &
                               METHOD_TYPE_DFT, METHOD_TYPE_MP2, METHOD_TYPE_CCSD_T
   use mqc_cuest_bridge, only: cuest_backend_available
   implicit none
   private
   public :: collect_mqc_acknowledgements_tests

contains

   subroutine collect_mqc_acknowledgements_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("semi_empirical_credits_tblite", test_tblite), &
                  new_unittest("explicit_request_is_honoured", test_explicit), &
                  new_unittest("auto_follows_the_build", test_auto), &
                  new_unittest("auto_keeps_cpu_only_methods_on_cpu", test_auto_cpu_only) &
                  ]
   end subroutine collect_mqc_acknowledgements_tests

   subroutine test_tblite(error)
      !! The semi-empirical methods are tblite's own and take no backend
      type(error_type), allocatable, intent(out) :: error

      call check(error, acknowledged_backend(METHOD_TYPE_GFN1, "auto") == ACK_TBLITE, &
                 "GFN1 runs through tblite")
      if (allocated(error)) return
      call check(error, acknowledged_backend(METHOD_TYPE_GFN2, "auto") == ACK_TBLITE, &
                 "GFN2 runs through tblite")
      if (allocated(error)) return
      ! Even spelled with a backend, which they ignore -- the credit must follow
      ! what ran, not what the deck happened to say.
      call check(error, acknowledged_backend(METHOD_TYPE_GFN2, "cpu") == ACK_TBLITE, &
                 "a backend name must not redirect the credit for xTB")
   end subroutine test_tblite

   subroutine test_explicit(error)
      !! A named backend is credited as named, on any build
      type(error_type), allocatable, intent(out) :: error

      call check(error, acknowledged_backend(METHOD_TYPE_HF, "cpu") == ACK_LIBCINT, &
                 "'cpu' credits libcint")
      if (allocated(error)) return
      call check(error, acknowledged_backend(METHOD_TYPE_HF, "libcint") == ACK_LIBCINT, &
                 "'libcint' credits libcint")
      if (allocated(error)) return
      call check(error, acknowledged_backend(METHOD_TYPE_DFT, "gpu") == ACK_CUEST, &
                 "'gpu' credits cuEST")
      if (allocated(error)) return
      call check(error, acknowledged_backend(METHOD_TYPE_HF, "cuest") == ACK_CUEST, &
                 "'cuest' credits cuEST")
   end subroutine test_explicit

   subroutine test_auto(error)
      !! `auto` resolves to whatever this build would actually dispatch to
      type(error_type), allocatable, intent(out) :: error

      integer :: expected

      if (cuest_backend_available()) then
         expected = ACK_CUEST
      else
         expected = ACK_LIBCINT
      end if

      call check(error, acknowledged_backend(METHOD_TYPE_HF, "auto") == expected, &
                 "auto HF must follow the build")
      if (allocated(error)) return
      call check(error, acknowledged_backend(METHOD_TYPE_DFT, "") == expected, &
                 "an empty backend is auto, and must follow the build too")
   end subroutine test_auto

   subroutine test_auto_cpu_only(error)
      !! A correlated method is CPU-only, so `auto` credits libcint either way
      !!
      !! This is the case a build with cuEST gets wrong if the resolution
      !! forgets to ask whether cuEST implements the method at all.
      type(error_type), allocatable, intent(out) :: error

      call check(error, acknowledged_backend(METHOD_TYPE_MP2, "auto") == ACK_LIBCINT, &
                 "MP2 has no GPU implementation, so auto is libcint")
      if (allocated(error)) return
      call check(error, acknowledged_backend(METHOD_TYPE_CCSD_T, "auto") == ACK_LIBCINT, &
                 "CCSD(T) has no GPU implementation, so auto is libcint")
   end subroutine test_auto_cpu_only

end module test_mqc_acknowledgements

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_acknowledgements, only: collect_mqc_acknowledgements_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_acknowledgements", collect_mqc_acknowledgements_tests) &
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
