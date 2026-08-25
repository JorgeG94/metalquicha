!! Coupled-cluster amplitudes surviving a restart, and refusing not to
!!
!! What is checked here is the set of failures that would otherwise be
!! *silent* -- a resumed run that iterates smoothly to a number nobody can
!! tell is wrong:
!!
!!   * `t2` and `l2` come back with the right **shape**, not merely the right
!!     element count. A rank-4 block read as the wrong extents is the same
!!     total size and a different tensor
!!   * a file whose last write was interrupted reports `resumable = .false.`
!!     rather than handing back half an iteration. This is the one that has no
!!     loud failure mode at all: the arrays are present and readable, they
!!     just do not solve anything
!!   * amplitudes from a **different calculation** are refused, not loaded.
!!     Resuming from the wrong fingerprint converges -- to a stationary point
!!     of another set of equations
!!   * asking for a block by the wrong rank is refused rather than reshaped
!!
!! None of this is reachable from a chemical deck. A CC run either resumes or
!! it does not; it cannot be killed at a chosen instant and then asked which
!! numbers came back. That is why it belongs in the unit suite.
module test_mqc_hdf5_amplitudes
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_hdf5_amplitudes, only: hdf5_amplitudes_t, hdf5_amplitudes_available
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_hdf5_amplitudes

   character(len=*), parameter :: PATH = "test_amps.h5"
   character(len=*), parameter :: FP = "a41c07de9b3f5521"
   character(len=*), parameter :: OTHER_FP = "0000deadbeef0000"

   integer, parameter :: NOCC = 3, NVIR = 5

contains

   subroutine collect_hdf5_amplitudes(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("round trip preserves shape and values", test_round_trip), &
                  new_unittest("an interrupted write is not resumable", test_torn_write), &
                  new_unittest("a foreign fingerprint is refused", test_wrong_fingerprint), &
                  new_unittest("a missing block is refused", test_missing_block), &
                  new_unittest("the wrong rank is refused", test_wrong_rank), &
                  new_unittest("a second commit overwrites in place", test_overwrite) &
                  ]
   end subroutine collect_hdf5_amplitudes

   subroutine test_round_trip(error)
      type(error_type), allocatable, intent(out) :: error

      type(hdf5_amplitudes_t) :: store
      type(error_t) :: err
      real(dp), allocatable :: t1(:, :), t2(:, :, :, :)
      real(dp), allocatable :: back1(:, :), back2(:, :, :, :)

      if (.not. hdf5_amplitudes_available()) return
      call fresh()
      call make_t1(t1)
      call make_t2(t2)

      call store%open(PATH, FP, err)
      call check(error,.not. err%has_error(), "opening a new store should succeed")
      if (allocated(error)) return
      call check(error,.not. store%resumable, "a new store has nothing to resume")
      if (allocated(error)) return

      call store%put("t1", t1)
      call store%put("t2", t2)
      call store%commit(7, -0.25_dp)
      call store%close()

      call store%open(PATH, FP, err)
      call check(error,.not. err%has_error(), "reopening should succeed")
      if (allocated(error)) return
      call check(error, store%resumable, "a committed store is resumable")
      if (allocated(error)) return
      call check(error, store%iteration, 7)
      if (allocated(error)) return
      call check(error, store%energy, -0.25_dp, thr=0.0_dp)
      if (allocated(error)) return

      call store%get("t1", back1, err)
      call check(error,.not. err%has_error(), "t1 should read back")
      if (allocated(error)) return
      call store%get("t2", back2, err)
      call check(error,.not. err%has_error(), "t2 should read back")
      if (allocated(error)) return
      call store%close()

      ! Shape before values: a wrong extent that happens to hold the right
      ! numbers in the wrong places is the failure this is here to catch.
      call check(error, all(shape(back1) == [NOCC, NVIR]), "t1 shape")
      if (allocated(error)) return
      call check(error, all(shape(back2) == [NOCC, NVIR, NOCC, NVIR]), "t2 shape")
      if (allocated(error)) return
      call check(error, maxval(abs(back1 - t1)), 0.0_dp, thr=0.0_dp)
      if (allocated(error)) return
      call check(error, maxval(abs(back2 - t2)), 0.0_dp, thr=0.0_dp)
   end subroutine test_round_trip

   subroutine test_torn_write(error)
      !! Write, then write again without committing, and reopen
      !!
      !! This is what a kill looks like from the next run's side: the datasets
      !! exist and are readable, and there is nothing in their contents to say
      !! they are half of one iteration and half of another. Only `complete`
      !! knows, which is why `put` clears it before touching the data.
      type(error_type), allocatable, intent(out) :: error

      type(hdf5_amplitudes_t) :: store
      type(error_t) :: err
      real(dp), allocatable :: t1(:, :), t2(:, :, :, :)

      if (.not. hdf5_amplitudes_available()) return
      call fresh()
      call make_t1(t1)
      call make_t2(t2)

      call store%open(PATH, FP, err)
      call store%put("t1", t1)
      call store%put("t2", t2)
      call store%commit(1, -0.1_dp)
      ! A second iteration that never reaches its commit.
      call store%put("t1", 2.0_dp*t1)
      call store%close()

      call store%open(PATH, FP, err)
      call check(error,.not. err%has_error(), "the file should still open")
      if (allocated(error)) return
      call check(error,.not. store%resumable, &
                 "an uncommitted write must leave the store unresumable")
      call store%close()
   end subroutine test_torn_write

   subroutine test_wrong_fingerprint(error)
      type(error_type), allocatable, intent(out) :: error

      type(hdf5_amplitudes_t) :: store
      type(error_t) :: err
      real(dp), allocatable :: t1(:, :)

      if (.not. hdf5_amplitudes_available()) return
      call fresh()
      call make_t1(t1)

      call store%open(PATH, FP, err)
      call store%put("t1", t1)
      call store%commit(1, -0.1_dp)
      call store%close()

      call store%open(PATH, OTHER_FP, err)
      call check(error, err%has_error(), &
                 "amplitudes from another calculation must be refused")
      call store%close()
   end subroutine test_wrong_fingerprint

   subroutine test_missing_block(error)
      type(error_type), allocatable, intent(out) :: error

      type(hdf5_amplitudes_t) :: store
      type(error_t) :: err
      real(dp), allocatable :: t1(:, :), back(:, :, :, :)

      if (.not. hdf5_amplitudes_available()) return
      call fresh()
      call make_t1(t1)

      call store%open(PATH, FP, err)
      call store%put("t1", t1)
      call store%commit(1, -0.1_dp)
      call store%get("l2", back, err)
      call check(error, err%has_error(), "a block never written must be refused")
      call store%close()
   end subroutine test_missing_block

   subroutine test_wrong_rank(error)
      !! Ask for the rank-2 block as rank-4
      type(error_type), allocatable, intent(out) :: error

      type(hdf5_amplitudes_t) :: store
      type(error_t) :: err
      real(dp), allocatable :: t1(:, :), back(:, :, :, :)

      if (.not. hdf5_amplitudes_available()) return
      call fresh()
      call make_t1(t1)

      call store%open(PATH, FP, err)
      call store%put("t1", t1)
      call store%commit(1, -0.1_dp)
      call store%get("t1", back, err)
      call check(error, err%has_error(), "reading a rank-2 block as rank-4 must be refused")
      call store%close()
   end subroutine test_wrong_rank

   subroutine test_overwrite(error)
      !! The file must not grow a second copy when an iteration is rewritten
      type(error_type), allocatable, intent(out) :: error

      type(hdf5_amplitudes_t) :: store
      type(error_t) :: err
      real(dp), allocatable :: t1(:, :), back(:, :)

      if (.not. hdf5_amplitudes_available()) return
      call fresh()
      call make_t1(t1)

      call store%open(PATH, FP, err)
      call store%put("t1", t1)
      call store%commit(1, -0.1_dp)
      call store%put("t1", 3.0_dp*t1)
      call store%commit(2, -0.2_dp)
      call store%close()

      call store%open(PATH, FP, err)
      call check(error, store%iteration, 2)
      if (allocated(error)) return
      call store%get("t1", back, err)
      call store%close()
      call check(error,.not. err%has_error(), "the rewritten block should read back")
      if (allocated(error)) return
      call check(error, all(shape(back) == [NOCC, NVIR]), "shape after rewrite")
      if (allocated(error)) return
      call check(error, maxval(abs(back - 3.0_dp*t1)), 0.0_dp, thr=0.0_dp)
      if (allocated(error)) return
      call fresh()
   end subroutine test_overwrite

   ! -- fixtures ------------------------------------------------------------
   ! Every element encodes its own indices, so a block that came back
   ! transposed or offset is caught by its values and not only by its size.

   subroutine make_t1(t1)
      real(dp), allocatable, intent(out) :: t1(:, :)
      integer :: i, a

      allocate (t1(NOCC, NVIR))
      do a = 1, NVIR
         do i = 1, NOCC
            t1(i, a) = real(100*i + a, dp)
         end do
      end do
   end subroutine make_t1

   subroutine make_t2(t2)
      real(dp), allocatable, intent(out) :: t2(:, :, :, :)
      integer :: i, a, j, b

      allocate (t2(NOCC, NVIR, NOCC, NVIR))
      do b = 1, NVIR
         do j = 1, NOCC
            do a = 1, NVIR
               do i = 1, NOCC
                  t2(i, a, j, b) = real(1000*i + 100*a + 10*j + b, dp)
               end do
            end do
         end do
      end do
   end subroutine make_t2

   subroutine fresh()
      !! Remove the scratch file, so a test starts from nothing
      integer :: unit
      logical :: exists

      inquire (file=PATH, exist=exists)
      if (exists) then
         open (newunit=unit, file=PATH, status="old", action="readwrite")
         close (unit, status="delete")
      end if
   end subroutine fresh

end module test_mqc_hdf5_amplitudes

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_hdf5_amplitudes, only: collect_hdf5_amplitudes
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_hdf5_amplitudes", collect_hdf5_amplitudes)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
