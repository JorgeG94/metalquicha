!! Gradients and Hessians surviving a checkpoint round trip
!!
!! The text checkpoint held one energy per line, so a gradient run resuming
!! from one came back without its gradient and died. This is the format that
!! fixes that, and what is checked here is the part that would fail *silently*
!! rather than loudly:
!!
!!   * a gradient comes back with the right **shape**, not merely the right
!!     count. Fragments have different atom counts, so the payload is ragged
!!     and packed; get the offsets wrong by one record and every fragment
!!     reads its neighbour's numbers, all correctly sized
!!   * records are found by monomer tuple after being **sorted**, and the
!!     slices have to follow their records. This is the one that broke first:
!!     cumulative offsets cannot be permuted
!!   * a checkpoint from a different calculation is refused rather than
!!     half-believed
!!
!! None of that is reachable from a chemical system. A deck either resumes or
!! it does not; it cannot ask for three fragments of different sizes written
!! out of order and then check which numbers came back where. That is why this
!! belongs in the unit suite rather than in a validation deck.
module test_mqc_hdf5_checkpoint
   use pic_types, only: dp, int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_hdf5_checkpoint, only: hdf5_checkpoint_t
   use mqc_error, only: error_t
   implicit none
   private

   public :: collect_hdf5_checkpoint

   character(len=*), parameter :: PATH = "test_ckpt.h5"
   character(len=*), parameter :: FP = "f5f16b57fa5252ab"

   !> Every element of a block says which block it came from, so a slice that
   !> travelled to the wrong record is caught by its values and not only by
   !> its size.
   real(dp), parameter :: BASE_MONOMER_G = 100.0_dp
   real(dp), parameter :: BASE_MONOMER_H = 10.0_dp
   real(dp), parameter :: BASE_DIMER_23 = 200.0_dp
   real(dp), parameter :: BASE_DIMER_12 = 300.0_dp

contains

   subroutine collect_hdf5_checkpoint(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ragged_records_survive_a_round_trip", test_round_trip), &
                  new_unittest("another_calculation_is_refused", test_wrong_fingerprint) &
                  ]
   end subroutine collect_hdf5_checkpoint

   subroutine test_round_trip(error)
      !! Three fragments of different sizes, written out of order, read back
      !!
      !! Out of order on purpose: the reader sorts, and the slices have to be
      !! carried along by that sort rather than left behind.
      type(error_type), allocatable, intent(out) :: error

      type(hdf5_checkpoint_t) :: ck, back
      type(error_t) :: err
      logical :: found
      real(dp) :: energy
      integer :: status, natoms
      real(dp), allocatable :: g(:, :), h(:, :)

      call fresh()

      call ck%open(PATH, FP, 2, err)
      call check(error, .not. err%has_error(), "open for writing: "//err%get_message())
      if (allocated(error)) return

      call ck%record([2, 3], -2.5_dp, 1, 6, gradient=ramp(3, 6, BASE_DIMER_23))
      call ck%record([1, 0], -0.5_dp, 1, 3, gradient=ramp(3, 3, BASE_MONOMER_G), &
                     hessian=ramp(9, 9, BASE_MONOMER_H))
      call ck%record([1, 2], -1.5_dp, 2, 9, gradient=ramp(3, 9, BASE_DIMER_12))
      call ck%close()

      call back%open(PATH, FP, 2, err)
      call check(error, .not. err%has_error(), "reopen: "//err%get_message())
      if (allocated(error)) return
      call check(error, back%n_loaded == 3_int64, "three records came back")
      if (allocated(error)) return

      ! The middle fragment by size, and the first by sort order.
      call back%lookup([1, 0], found, energy, status, natoms, gradient=g, hessian=h)
      call check(error, found, "the monomer should be found")
      if (allocated(error)) return
      call check(error, energy, -0.5_dp, "monomer energy", thr=1.0e-14_dp)
      if (allocated(error)) return
      call check(error, natoms, 3, "monomer atom count")
      if (allocated(error)) return
      call check(error, allocated(g), "the monomer gradient should be present")
      if (allocated(error)) return
      call check(error, size(g, 1) == 3 .and. size(g, 2) == 3, "monomer gradient shape")
      if (allocated(error)) return
      call check(error, matches(g, BASE_MONOMER_G), "monomer gradient values")
      if (allocated(error)) return
      call check(error, allocated(h), "the monomer hessian should be present")
      if (allocated(error)) return
      call check(error, size(h, 1) == 9 .and. size(h, 2) == 9, "monomer hessian shape")
      if (allocated(error)) return
      call check(error, matches(h, BASE_MONOMER_H), "monomer hessian values")
      if (allocated(error)) return

      ! Written first, moved last by the sort. If the offsets did not travel
      ! with it this reads someone else's numbers -- at the right size, which
      ! is exactly what makes it worth checking.
      deallocate (g)
      call back%lookup([2, 3], found, energy, status, natoms, gradient=g)
      call check(error, found, "dimer (2,3) should be found")
      if (allocated(error)) return
      call check(error, natoms, 6, "dimer (2,3) atom count")
      if (allocated(error)) return
      call check(error, size(g, 2), 6, "dimer (2,3) gradient shape")
      if (allocated(error)) return
      call check(error, matches(g, BASE_DIMER_23), &
                 "dimer (2,3) read its neighbour's gradient values")
      if (allocated(error)) return

      deallocate (g)
      call back%lookup([1, 2], found, energy, status, natoms, gradient=g)
      call check(error, found, "dimer (1,2) should be found")
      if (allocated(error)) return
      call check(error, status, 2, "the scf status should survive")
      if (allocated(error)) return
      call check(error, size(g, 2), 9, "dimer (1,2) gradient shape")
      if (allocated(error)) return
      call check(error, matches(g, BASE_DIMER_12), "dimer (1,2) gradient values")
      if (allocated(error)) return

      call back%lookup([4, 5], found, energy, status, natoms)
      call check(error, .not. found, "a term nobody wrote must not be found")

      call back%close()
      call fresh()
   end subroutine test_round_trip

   subroutine test_wrong_fingerprint(error)
      !! A checkpoint from a different calculation is refused
      !!
      !! The fingerprint covers the geometry, the partition, the method and the
      !! thresholds. Resuming across a change in any of them would mix results
      !! from two different calculations into one total, which is worse than
      !! starting again.
      type(error_type), allocatable, intent(out) :: error

      type(hdf5_checkpoint_t) :: ck, other
      type(error_t) :: err, other_err

      call fresh()
      call ck%open(PATH, FP, 2, err)
      call check(error, .not. err%has_error(), "open: "//err%get_message())
      if (allocated(error)) return
      call ck%record([1, 0], -0.5_dp, 1, 3)
      call ck%close()

      call other%open(PATH, "0000000000000000", 2, other_err)
      call check(error, other_err%has_error(), &
                 "a checkpoint from another calculation must be refused")
      call other%close()
      call fresh()
   end subroutine test_wrong_fingerprint

   function ramp(rows, cols, base) result(a)
      !! A block whose every element says which block it came from
      integer, intent(in) :: rows, cols
      real(dp), intent(in) :: base
      real(dp), allocatable :: a(:, :)

      integer :: i, j

      allocate (a(rows, cols))
      do j = 1, cols
         do i = 1, rows
            a(i, j) = base + real(i, dp) + 0.001_dp*real(j, dp)
         end do
      end do
   end function ramp

   function matches(a, base) result(same)
      !! Does this block still say it came from `base`
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(in) :: base
      logical :: same

      integer :: i, j

      same = .true.
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            if (abs(a(i, j) - (base + real(i, dp) + 0.001_dp*real(j, dp))) > 1.0e-12_dp) then
               same = .false.
               return
            end if
         end do
      end do
   end function matches

   subroutine fresh()
      !! Remove the scratch checkpoint, so a test starts from nothing
      integer :: unit
      logical :: exists

      inquire (file=PATH, exist=exists)
      if (exists) then
         open (newunit=unit, file=PATH, status="old", action="readwrite")
         close (unit, status="delete")
      end if
   end subroutine fresh

end module test_mqc_hdf5_checkpoint

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_hdf5_checkpoint, only: collect_hdf5_checkpoint
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_hdf5_checkpoint", collect_hdf5_checkpoint)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
