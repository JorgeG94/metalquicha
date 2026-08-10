!! Manual check that gradients and Hessians survive a checkpoint round trip
!!
!!     cmake -B build -DMQC_ENABLE_HDF5=ON && ./build/check_hdf5_checkpoint
!!
!! The text checkpoint holds one energy per line, so a gradient run resuming
!! from one came back without its gradient and died. This is the format that
!! fixes that, and the properties worth checking are the ones that would fail
!! silently rather than loudly:
!!
!!   * A gradient comes back with the right *shape*, not merely the right
!!     count. Fragments have different atom counts, so the payload is ragged
!!     and packed; get the offsets wrong by one record and every fragment
!!     reads its neighbour's numbers, all correctly sized.
!!   * Records are found by monomer tuple after being **sorted**, and the
!!     slices must follow their records. This is the one that broke first:
!!     cumulative offsets cannot be permuted.
!!   * `n_valid` bounds what a kill costs. Records written after the last
!!     commit must be ignored on reload -- not read, not half-read.
program check_hdf5_checkpoint
   use pic_types, only: dp, int64
   use mqc_hdf5_checkpoint, only: hdf5_checkpoint_t
   use mqc_error, only: error_t
   implicit none

   character(len=*), parameter :: PATH = "check_ckpt.h5"
   character(len=*), parameter :: FP = "f5f16b57fa5252ab"

   type(hdf5_checkpoint_t) :: ck
   type(error_t) :: err
   integer :: failures

   failures = 0
   call fresh()

   call phase_write()
   call phase_read()
   call phase_wrong_fingerprint()

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[ckpt] all ok -- gradients and Hessians survive the round trip"
   else
      write (*, "(A,I0,A)") "[ckpt] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine phase_write()
      !! Three fragments of different sizes, written out of order
      !!
      !! Out of order on purpose: the reader sorts, and the slices have to be
      !! carried along by that sort rather than left behind.
      call ck%open(PATH, FP, 2, err)
      call expect(.not. err%has_error(), "open for writing: "//err%get_message())

      call ck%record([2, 3], -2.5_dp, 1, 6, gradient=ramp(3, 6, 200.0_dp))
      call ck%record([1, 0], -0.5_dp, 1, 3, gradient=ramp(3, 3, 100.0_dp), &
                     hessian=ramp(9, 9, 10.0_dp))
      call ck%record([1, 2], -1.5_dp, 2, 9, gradient=ramp(3, 9, 300.0_dp))
      call ck%close()
   end subroutine phase_write

   subroutine phase_read()
      type(hdf5_checkpoint_t) :: back
      logical :: found
      real(dp) :: energy
      integer :: status, natoms
      real(dp), allocatable :: g(:, :), h(:, :)

      call back%open(PATH, FP, 2, err)
      call expect(.not. err%has_error(), "reopen: "//err%get_message())
      call expect(back%n_loaded == 3_int64, "three records came back")

      ! The middle fragment by size, and the first by sort order.
      call back%lookup([1, 0], found, energy, status, natoms, gradient=g, hessian=h)
      call expect(found, "monomer found")
      call expect(abs(energy + 0.5_dp) < 1.0e-14_dp, "monomer energy")
      call expect(natoms == 3, "monomer atom count")
      call expect(allocated(g), "monomer gradient present")
      if (allocated(g)) then
         call expect(size(g, 1) == 3 .and. size(g, 2) == 3, "monomer gradient shape")
         call expect(matches(g, 100.0_dp), "monomer gradient values")
      end if
      call expect(allocated(h), "monomer hessian present")
      if (allocated(h)) then
         call expect(size(h, 1) == 9 .and. size(h, 2) == 9, "monomer hessian shape")
         call expect(matches(h, 10.0_dp), "monomer hessian values")
      end if

      ! The one written first, which sorting moved last. If the offsets did
      ! not travel with it, this reads someone else's numbers -- at the right
      ! size, which is what makes it worth checking.
      if (allocated(g)) deallocate (g)
      call back%lookup([2, 3], found, energy, status, natoms, gradient=g)
      call expect(found, "dimer (2,3) found")
      call expect(natoms == 6, "dimer (2,3) atom count")
      if (allocated(g)) then
         call expect(size(g, 2) == 6, "dimer (2,3) gradient shape")
         call expect(matches(g, 200.0_dp), "dimer (2,3) gradient values -- not its neighbour's")
      end if

      if (allocated(g)) deallocate (g)
      call back%lookup([1, 2], found, energy, status, natoms, gradient=g)
      call expect(found, "dimer (1,2) found")
      call expect(status == 2, "scf status survived")
      if (allocated(g)) then
         call expect(size(g, 2) == 9, "dimer (1,2) gradient shape")
         call expect(matches(g, 300.0_dp), "dimer (1,2) gradient values")
      end if

      call back%lookup([4, 5], found, energy, status, natoms)
      call expect(.not. found, "a term nobody wrote is not found")

      call back%close()
   end subroutine phase_read

   subroutine phase_wrong_fingerprint()
      type(hdf5_checkpoint_t) :: other
      type(error_t) :: other_err

      call other%open(PATH, "0000000000000000", 2, other_err)
      call expect(other_err%has_error(), &
                  "a checkpoint from another calculation must be refused")
      call other%close()
   end subroutine phase_wrong_fingerprint

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

   subroutine expect(condition, what)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: what

      if (.not. condition) then
         write (*, "(A)") "  FAILED: "//trim(what)
         failures = failures + 1
      end if
   end subroutine expect

   subroutine fresh()
      integer :: unit
      logical :: exists
      inquire (file=PATH, exist=exists)
      if (exists) then
         open (newunit=unit, file=PATH, status="old", action="readwrite")
         close (unit, status="delete")
      end if
   end subroutine fresh

end program check_hdf5_checkpoint
