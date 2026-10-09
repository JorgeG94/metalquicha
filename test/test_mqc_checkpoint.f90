!! Tests for the append-as-you-go checkpoint
module test_mqc_checkpoint
   !! Three properties, and one of them is the reason the file exists.
   !!
   !!   * What was written comes back, keyed on the monomers rather than the
   !!     position, so a differently-ordered rerun still finds it.
   !!   * A half-written last line is dropped, because that is what a killed
   !!     job leaves and recovering from a killed job is the point.
   !!   * A checkpoint from a different calculation is **refused**. Silently
   !!     reusing those energies would produce a converged total belonging to
   !!     neither run, which is the failure mode this whole mechanism is built
   !!     to make impossible.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp, int64
   use mqc_checkpoint, only: checkpoint_t
   use mqc_hdf5_checkpoint, only: hdf5_checkpoint_available
   use mqc_error, only: error_t
   use mqc_result_types, only: SCF_CONVERGED, SCF_NOT_CONVERGED
   implicit none
   private

   public :: collect_mqc_checkpoint

   character(len=*), parameter :: PATH = "test_checkpoint.tmp"
   character(len=*), parameter :: FP = "0123456789abcdef"

contains

   subroutine collect_mqc_checkpoint(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("roundtrip", test_roundtrip), &
                  new_unittest("absent_term_not_found", test_absent), &
                  new_unittest("status_survives", test_status), &
                  new_unittest("truncated_line_dropped", test_truncated), &
                  new_unittest("wrong_fingerprint_refused", test_wrong_fingerprint), &
                  new_unittest("not_a_checkpoint_refused", test_not_a_checkpoint), &
                  new_unittest("derivative_run_refused", test_derivative_refused), &
                  new_unittest("a_term_wider_than_the_level_round_trips", test_wide_terms) &
                  ]
   end subroutine collect_mqc_checkpoint

   subroutine test_roundtrip(error)
      !! Written, closed, reopened, found -- and found by monomers
      !!
      !! The lookup is deliberately made in a different order from the writes.
      !! A checkpoint keyed on position would pass a same-order test and pair
      !! every energy with the wrong fragment the first time a screened list
      !! was resumed.
      type(error_type), allocatable, intent(out) :: error
      type(checkpoint_t) :: ck
      type(error_t) :: err
      logical :: found
      real(dp) :: energy
      integer :: status

      call fresh()
      call ck%open(PATH, FP, 2, .true., err)
      call ck%record([1, 2], -1.5_dp, SCF_CONVERGED)
      call ck%record([1, 0], -0.5_dp, SCF_CONVERGED)
      call ck%record([2, 3], -2.5_dp, SCF_CONVERGED)
      call ck%close()

      call reopen(ck, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, int(ck%n_loaded), 3)
      if (allocated(error)) return

      call ck%lookup([2, 3], found, energy, status)
      call check(error, found, "a recorded term must be found")
      if (allocated(error)) return
      call check(error, abs(energy + 2.5_dp) < 1.0e-12_dp)
      if (allocated(error)) return

      call ck%lookup([1, 0], found, energy, status)
      call check(error, found, "a monomer must be found")
      if (allocated(error)) return
      call check(error, abs(energy + 0.5_dp) < 1.0e-12_dp)
      call ck%close()
      call cleanup()
   end subroutine test_roundtrip

   subroutine test_absent(error)
      !! A term nobody computed is not quietly reported as done
      type(error_type), allocatable, intent(out) :: error
      type(checkpoint_t) :: ck
      type(error_t) :: err
      logical :: found
      real(dp) :: energy
      integer :: status

      call fresh()
      call ck%open(PATH, FP, 2, .true., err)
      call ck%record([1, 2], -1.5_dp, SCF_CONVERGED)
      call ck%close()

      call reopen(ck, err)
      call ck%lookup([1, 3], found, energy, status)
      call check(error,.not. found, "an unrecorded term must not be found")
      call ck%close()
      call cleanup()
   end subroutine test_absent

   subroutine test_status(error)
      !! Convergence survives the round trip
      !!
      !! A resumed run that forgot this would report every reused fragment as
      !! converged, which is the silence the scf column was added to break.
      type(error_type), allocatable, intent(out) :: error
      type(checkpoint_t) :: ck
      type(error_t) :: err
      logical :: found
      real(dp) :: energy
      integer :: status

      call fresh()
      call ck%open(PATH, FP, 2, .true., err)
      call ck%record([1, 2], -1.5_dp, SCF_NOT_CONVERGED)
      call ck%close()

      call reopen(ck, err)
      call ck%lookup([1, 2], found, energy, status)
      call check(error, found)
      if (allocated(error)) return
      call check(error, status, SCF_NOT_CONVERGED)
      call ck%close()
      call cleanup()
   end subroutine test_status

   subroutine test_truncated(error)
      !! The job died mid-write; keep everything before that
      type(error_type), allocatable, intent(out) :: error
      type(checkpoint_t) :: ck
      type(error_t) :: err
      integer :: unit

      call fresh()
      call ck%open(PATH, FP, 2, .true., err)
      call ck%record([1, 2], -1.5_dp, SCF_CONVERGED)
      call ck%record([1, 3], -2.5_dp, SCF_CONVERGED)
      call ck%close()

      ! Half a line, as a kill leaves it.
      open (newunit=unit, file=PATH, status="old", position="append", action="write")
      write (unit, "(a)") "2 2 3   -1.15"
      close (unit)

      call reopen(ck, err)
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call check(error, int(ck%n_loaded), 2)
      call ck%close()
      call cleanup()
   end subroutine test_truncated

   subroutine test_wrong_fingerprint(error)
      !! The one that matters
      type(error_type), allocatable, intent(out) :: error
      type(checkpoint_t) :: ck
      type(error_t) :: err

      call fresh()
      call ck%open(PATH, FP, 2, .true., err)
      call ck%record([1, 2], -1.5_dp, SCF_CONVERGED)
      call ck%close()

      call ck%open(PATH, "fedcba9876543210", 2, .true., err)
      call check(error, err%has_error(), &
                 "a checkpoint from another calculation must be refused, not reused")
      call ck%close()
      call cleanup()
   end subroutine test_wrong_fingerprint

   subroutine test_not_a_checkpoint(error)
      !! Pointed at some other file, refuse rather than truncate it
      type(error_type), allocatable, intent(out) :: error
      type(checkpoint_t) :: ck
      type(error_t) :: err
      integer :: unit

      call fresh()
      open (newunit=unit, file=PATH, status="new", action="write")
      write (unit, "(a)") "something else entirely"
      close (unit)

      call ck%open(PATH, FP, 2, .true., err)
      call check(error, err%has_error(), "a non-checkpoint file must be refused")
      call ck%close()
      call cleanup()
   end subroutine test_not_a_checkpoint

   subroutine test_derivative_refused(error)
      !! A derivative run gets a checkpoint only if HDF5 can hold one
      !!
      !! A text checkpoint line holds an energy and nothing else, so a gradient
      !! run reusing one comes back without the gradient and dies assembling
      !! the total -- which it did, before this was refused.
      !!
      !! With HDF5 compiled in that refusal is wrong: storing derivatives is
      !! the entire reason the backend exists, and open() routes a derivative
      !! run to it. So which outcome is correct depends on the build, and the
      !! test asks rather than assumes. Asserting refusal unconditionally
      !! passed for as long as nothing built with HDF5 enabled, and failed the
      !! moment something did.
      type(error_type), allocatable, intent(out) :: error
      type(checkpoint_t) :: ck
      type(error_t) :: err

      call fresh()
      call ck%open(PATH, FP, 2, .false., err)

      if (hdf5_checkpoint_available()) then
         call check(error,.not. err%has_error(), &
                    "with HDF5 built in, a derivative run must get a checkpoint: "// &
                    err%get_message())
      else
         call check(error, err%has_error(), &
                    "without HDF5, a run needing derivatives must be refused a checkpoint")
      end if

      call ck%close()
      call cleanup()
   end subroutine test_derivative_refused

   subroutine test_wide_terms(error)
      !! A ghosted term as wide as the system is kept whole, and found by all of it
      !!
      !! The file used to be opened at the MBE level plus one entry, which cut
      !! a full-cluster-basis row to its first few monomers. Two rows that share
      !! that prefix then stood for each other on a restart. These two agree on
      !! their first three entries and differ only past them, and the width is
      !! thirteen: twelve monomers and the displacement code.
      type(error_type), allocatable, intent(out) :: error
      type(checkpoint_t) :: ck, back
      type(error_t) :: err
      logical :: found
      real(dp) :: energy
      integer :: status
      integer, parameter :: WIDTH = 13
      integer :: all_ghosted(WIDTH), few_ghosted(WIDTH), absent(WIDTH)

      all_ghosted = [1, 2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, 0]
      few_ghosted = [1, 2, -3, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0]
      absent = [1, 2, -3, -4, -5, -6, -7, -8, -9, -10, -11, 0, 0]

      call fresh()
      call ck%open(PATH, FP, WIDTH, .true., err)
      call check(error,.not. err%has_error(), "could not open a checkpoint thirteen wide")
      if (allocated(error)) return
      call ck%record(all_ghosted, -1.5_dp, SCF_CONVERGED)
      call ck%record(few_ghosted, -2.5_dp, SCF_CONVERGED)
      call ck%close()

      call back%open(PATH, FP, WIDTH, .true., err)
      call check(error,.not. err%has_error(), "could not reopen it")
      if (allocated(error)) return
      call check(error, back%n_loaded == 2, "both wide terms should have been loaded")
      if (allocated(error)) return

      call back%lookup(all_ghosted, found, energy, status)
      call check(error, found, "the term ghosting ten monomers was not found")
      if (allocated(error)) return
      call check(error, abs(energy + 1.5_dp) < 1.0e-14_dp, "the term ghosting ten monomers came back changed")
      if (allocated(error)) return

      call back%lookup(few_ghosted, found, energy, status)
      call check(error, found, "the term ghosting two monomers was not found")
      if (allocated(error)) return
      call check(error, abs(energy + 2.5_dp) < 1.0e-14_dp, &
                 "the term ghosting two monomers was matched to its neighbour's energy")
      if (allocated(error)) return

      call back%lookup(absent, found, energy, status)
      call check(error,.not. found, "a term that was never recorded matched on its first entries")

      call back%close()
      call cleanup()
      if (allocated(error)) return

      ! Past a few hundred characters a line outgrows any fixed read buffer
      call round_trip_at_width(error, 150)
      if (allocated(error)) return
      call round_trip_at_width(error, 250)
   end subroutine test_wide_terms

   subroutine round_trip_at_width(error, width)
      !! A ghosted term `width` entries long comes back with every field exact
      !!
      !! The line holds the term, the energy, the status and the frontier
      !! orbitals, so it is several times `width` characters long. A read that
      !! cuts it drops the record or, cut inside the orbital fields, returns a
      !! frontier energy with its exponent changed.
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: width

      type(checkpoint_t) :: ck, back
      type(error_t) :: err
      integer :: term(width)
      integer :: k
      logical :: found, has_orbitals
      real(dp) :: energy, homo, lumo
      integer :: status
      real(dp), parameter :: ENERGY_IN = -1234.5678901234567_dp
      real(dp), parameter :: HOMO_IN = -4.0331768786648070e-01_dp
      real(dp), parameter :: LUMO_IN = 5.5312345678901234e+01_dp

      term = 0
      term(1:2) = [1, 2]
      do k = 3, width - 1
         term(k) = -k
      end do

      call fresh()
      call ck%open(PATH, FP, width, .true., err)
      call check(error,.not. err%has_error(), "could not open a checkpoint of width "//to_str(width))
      if (allocated(error)) return
      call ck%record(term, ENERGY_IN, SCF_CONVERGED, homo=HOMO_IN, lumo=LUMO_IN, has_orbitals=.true.)
      call ck%close()

      call back%open(PATH, FP, width, .true., err)
      call check(error,.not. err%has_error(), "could not reopen it at width "//to_str(width))
      if (allocated(error)) return
      call check(error, back%n_loaded == 1, "the record of width "//to_str(width)//" was dropped on loading")
      if (allocated(error)) return

      call back%lookup(term, found, energy, status, homo=homo, lumo=lumo, has_orbitals=has_orbitals)
      call check(error, found, "the term of width "//to_str(width)//" was not found")
      if (allocated(error)) return
      call check(error, energy == ENERGY_IN, "the energy of width "//to_str(width)//" came back changed")
      if (allocated(error)) return
      call check(error, has_orbitals, "the frontier orbitals of width "//to_str(width)//" were lost")
      if (allocated(error)) return
      call check(error, homo == HOMO_IN .and. lumo == LUMO_IN, &
                 "the frontier energies of width "//to_str(width)//" came back changed")
      call back%close()
      call cleanup()
   end subroutine round_trip_at_width

   function to_str(n) result(text)
      integer, intent(in) :: n
      character(len=:), allocatable :: text

      character(len=16) :: buffer

      write (buffer, "(i0)") n
      text = trim(buffer)
   end function to_str

   ! -- helpers ---------------------------------------------------------------

   subroutine reopen(ck, err)
      !! Load from disk into a fresh object
      type(checkpoint_t), intent(out) :: ck
      type(error_t), intent(out) :: err
      call ck%open(PATH, FP, 2, .true., err)
   end subroutine reopen

   subroutine fresh()
      call cleanup()
   end subroutine fresh

   subroutine cleanup()
      integer :: unit
      logical :: exists
      inquire (file=PATH, exist=exists)
      if (exists) then
         open (newunit=unit, file=PATH, status="old", action="readwrite")
         close (unit, status="delete")
      end if
   end subroutine cleanup

end module test_mqc_checkpoint

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_checkpoint, only: collect_mqc_checkpoint
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_checkpoint", collect_mqc_checkpoint)]
   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
