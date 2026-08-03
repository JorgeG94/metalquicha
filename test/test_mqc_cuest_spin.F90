#ifdef MQC_WITH_CUEST
module test_mqc_cuest_spin
   !! Unit tests for the open-shell occupation arithmetic.
   !!
   !! This is the only part of the unrestricted path that can be checked
   !! without a GPU, and it is exactly the part where an off-by-one would be
   !! invisible: the SCF would converge happily on the wrong number of
   !! electrons in each spin channel.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_cuest_scf, only: spin_occupations
   implicit none
   private
   public :: collect_mqc_cuest_spin_tests

contains

   subroutine collect_mqc_cuest_spin_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("closed_shell", test_closed_shell), &
                  new_unittest("doublet_radical", test_doublet), &
                  new_unittest("triplet_o2", test_triplet), &
                  new_unittest("inconsistent_rejected", test_inconsistent) &
                  ]
   end subroutine collect_mqc_cuest_spin_tests

   subroutine test_closed_shell(error)
      !! Water: 10 electrons, singlet -> 5 alpha, 5 beta
      type(error_type), allocatable, intent(out) :: error
      integer :: n_alpha, n_beta
      logical :: ok

      call spin_occupations(10, 1, n_alpha, n_beta, ok)
      call check(error, ok, "10 electrons as a singlet is consistent")
      if (allocated(error)) return
      call check(error, n_alpha, 5, "singlet alpha count")
      if (allocated(error)) return
      call check(error, n_beta, 5, "singlet beta count")
   end subroutine test_closed_shell

   subroutine test_doublet(error)
      !! OH radical: 9 electrons, doublet -> 5 alpha, 4 beta
      type(error_type), allocatable, intent(out) :: error
      integer :: n_alpha, n_beta
      logical :: ok

      call spin_occupations(9, 2, n_alpha, n_beta, ok)
      call check(error, ok, "9 electrons as a doublet is consistent")
      if (allocated(error)) return
      call check(error, n_alpha, 5, "doublet alpha count")
      if (allocated(error)) return
      call check(error, n_beta, 4, "doublet beta count")
      if (allocated(error)) return
      call check(error, n_alpha + n_beta, 9, "occupations must sum to the electron count")
   end subroutine test_doublet

   subroutine test_triplet(error)
      !! O2: 16 electrons, triplet -> 9 alpha, 7 beta
      type(error_type), allocatable, intent(out) :: error
      integer :: n_alpha, n_beta
      logical :: ok

      call spin_occupations(16, 3, n_alpha, n_beta, ok)
      call check(error, ok, "16 electrons as a triplet is consistent")
      if (allocated(error)) return
      call check(error, n_alpha, 9, "triplet alpha count")
      if (allocated(error)) return
      call check(error, n_beta, 7, "triplet beta count")
      if (allocated(error)) return
      call check(error, n_alpha - n_beta, 2, "2S must equal multiplicity - 1")
   end subroutine test_triplet

   subroutine test_inconsistent(error)
      !! An even electron count cannot be a doublet, and the unpaired count
      !! cannot exceed the electrons present. Both must be rejected rather
      !! than silently producing a fractional or negative occupation.
      type(error_type), allocatable, intent(out) :: error
      integer :: n_alpha, n_beta
      logical :: ok

      call spin_occupations(10, 2, n_alpha, n_beta, ok)
      call check(error,.not. ok, "10 electrons cannot be a doublet")
      if (allocated(error)) return

      call spin_occupations(9, 1, n_alpha, n_beta, ok)
      call check(error,.not. ok, "9 electrons cannot be a singlet")
      if (allocated(error)) return

      call spin_occupations(2, 5, n_alpha, n_beta, ok)
      call check(error,.not. ok, "more unpaired electrons than electrons must be rejected")
   end subroutine test_inconsistent

end module test_mqc_cuest_spin

program tester_mqc_cuest_spin
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_cuest_spin, only: collect_mqc_cuest_spin_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_cuest_spin", collect_mqc_cuest_spin_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_cuest_spin

#else
program tester_mqc_cuest_spin
   !! The cuEST backend is not part of the fpm build, so the open-shell
   !! occupation test it exercises is a no-op here. It runs in full under the
   !! CMake build, which is where the cuEST backend is compiled.
   implicit none
   write (*, '(A)') "# mqc_cuest_spin: skipped (cuEST backend not built)"
end program tester_mqc_cuest_spin
#endif
