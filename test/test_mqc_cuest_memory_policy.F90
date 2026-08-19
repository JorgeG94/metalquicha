#ifdef MQC_WITH_CUEST
module test_mqc_cuest_memory_policy
   !! Ties the hand-copied memory-policy numbers to the real cuEST bindings.
   !!
   !! `parse_df_derivative_memory` lives in `mqc_cuest_iface`, which is compiled
   !! on the fpm path where the cuEST bindings do not exist, so it cannot import
   !! the enumerators it returns -- it repeats their values as literals. That is
   !! a copy, and copies drift: the day cuEST renumbers or inserts a policy, a
   !! deck asking for `recompute` would quietly get `hostcache`, the gradient
   !! would still be right, and the only symptom would be running out of memory
   !! anyway. Nothing else in the build would notice.
   !!
   !! This test notices. It only builds where the bindings are present, which is
   !! exactly where the mismatch would matter.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_error, only: error_t
   use mqc_cuest_iface, only: parse_df_derivative_memory
   use cuest, only: CUEST_DFSYMMETRICDERIVATIVECOMPUTE_MEMORY_POLICY_DEVICECACHE, &
                    CUEST_DFSYMMETRICDERIVATIVECOMPUTE_MEMORY_POLICY_HOSTCACHE, &
                    CUEST_DFSYMMETRICDERIVATIVECOMPUTE_MEMORY_POLICY_OVERWRITE, &
                    CUEST_DFSYMMETRICDERIVATIVECOMPUTE_MEMORY_POLICY_RECOMPUTE
   implicit none
   private
   public :: collect_mqc_cuest_memory_policy_tests

contains

   subroutine collect_mqc_cuest_memory_policy_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("policy names match the bindings", test_policies_match) &
                  ]
   end subroutine collect_mqc_cuest_memory_policy_tests

   subroutine test_policies_match(error)
      !! Every name resolves to the enumerator cuEST actually defines
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: parse_error
      integer :: policy

      call parse_df_derivative_memory("devicecache", policy, parse_error)
      call check(error, policy, int(CUEST_DFSYMMETRICDERIVATIVECOMPUTE_MEMORY_POLICY_DEVICECACHE))
      if (allocated(error)) return

      call parse_df_derivative_memory("hostcache", policy, parse_error)
      call check(error, policy, int(CUEST_DFSYMMETRICDERIVATIVECOMPUTE_MEMORY_POLICY_HOSTCACHE))
      if (allocated(error)) return

      call parse_df_derivative_memory("overwrite", policy, parse_error)
      call check(error, policy, int(CUEST_DFSYMMETRICDERIVATIVECOMPUTE_MEMORY_POLICY_OVERWRITE))
      if (allocated(error)) return

      call parse_df_derivative_memory("recompute", policy, parse_error)
      call check(error, policy, int(CUEST_DFSYMMETRICDERIVATIVECOMPUTE_MEMORY_POLICY_RECOMPUTE))
      if (allocated(error)) return

      ! 'auto' must stay outside the enumeration. If cuEST ever numbered a
      ! policy -1 this would fire, which is the point: the sentinel and the
      ! vendor's values share one integer.
      call parse_df_derivative_memory("auto", policy, parse_error)
      call check(error, policy < 0, "'auto' is no longer distinguishable from a real policy")
   end subroutine test_policies_match

end module test_mqc_cuest_memory_policy

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_cuest_memory_policy, only: collect_mqc_cuest_memory_policy_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_cuest_memory_policy", &
                               collect_mqc_cuest_memory_policy_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
#else
program tester
   !! cuEST is not in this build, so there are no bindings to check against.
   implicit none
end program tester
#endif
