module test_mqc_mbe
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_mbe, only: compute_mbe, collect_unconverged, score_unconverged
   use mqc_result_types, only: SCF_CONVERGED, SCF_NOT_CONVERGED, SCF_UNKNOWN
   use mqc_result_types, only: calculation_result_t, mbe_result_t
   use pic_types, only: dp, int64
   implicit none
   private
   public :: collect_mqc_mbe_tests

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_mbe_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("mbe_monomers_only", test_mbe_monomers_only), &
                  new_unittest("vmfc_dipole_excludes_ghosted_rows", test_vmfc_dipole), &
                  new_unittest("mbe_simple_dimer", test_mbe_simple_dimer), &
                  new_unittest("mbe_sorted_order", test_mbe_sorted_order), &
                  new_unittest("mbe_reverse_order", test_mbe_reverse_order), &
                  new_unittest("mbe_random_order", test_mbe_random_order), &
                  new_unittest("unconverged_carry_their_monomers", test_unconverged_carry_their_monomers), &
                  new_unittest("unconverged_ignores_silent_methods", test_unconverged_ignores_silent_methods), &
                  new_unittest("converged_run_says_so", test_converged_run_says_so), &
                  new_unittest("failures_name_their_culprit", test_failures_name_their_culprit) &
                  ]
   end subroutine collect_mqc_mbe_tests

   subroutine test_vmfc_dipole(error)
      !! The counterpoise dipole uses the ghosted rows without summing them
      !!
      !! An auxiliary row -- one carrying a negative entry, such as `[1,-2]` --
      !! is a monomer solved in its parent's basis. It exists so the parent has
      !! something to subtract, and the energy path says so directly: "never to
      !! be summed". The dipole path did sum it, so every ghosted monomer went
      !! back into the total.
      !!
      !! It also sliced its row to the count of *real* monomers, which drops the
      !! ghost entirely, so the subset was looked up in the AB basis rather than
      !! the parent ABC one -- the wrong number at VMFC(2), and the wrong basis
      !! recursed at VMFC(3) and above.
      !!
      !! With per-fragment dipoles set by hand, the expansion has one right
      !! answer and both faults move it:
      !!
      !!     d = d_1 + d_2 + (d_AB - d_A(b) - d_B(a))
      type(error_type), allocatable, intent(out) :: error

      type(calculation_result_t), allocatable :: results(:)
      type(mbe_result_t) :: mbe_result
      integer, allocatable :: polymers(:, :)
      integer(int64) :: fragment_count
      integer :: max_level, i
      real(dp) :: expected(3)

      fragment_count = 5
      max_level = 2
      allocate (polymers(fragment_count, max_level))
      allocate (results(fragment_count))

      ! Two monomers in their own basis, the same two in the pair's basis, and
      ! the pair. The negative entries are what mark the ghosted rows.
      polymers(1, :) = [1, 0]
      polymers(2, :) = [2, 0]
      polymers(3, :) = [1, -2]
      polymers(4, :) = [-1, 2]
      polymers(5, :) = [1, 2]

      do i = 1, int(fragment_count)
         allocate (results(i)%dipole(3))
         results(i)%has_dipole = .true.
         results(i)%has_energy = .true.
      end do

      results(1)%energy%scf = -10.0_dp; results(1)%dipole = [1.0_dp, 0.0_dp, 0.0_dp]
      results(2)%energy%scf = -15.0_dp; results(2)%dipole = [0.0_dp, 1.0_dp, 0.0_dp]
      results(3)%energy%scf = -10.1_dp; results(3)%dipole = [1.1_dp, 0.0_dp, 0.0_dp]
      results(4)%energy%scf = -15.2_dp; results(4)%dipole = [0.0_dp, 1.2_dp, 0.0_dp]
      results(5)%energy%scf = -25.5_dp; results(5)%dipole = [2.0_dp, 2.0_dp, 0.0_dp]

      allocate (mbe_result%dipole(3))
      mbe_result%dipole = 0.0_dp

      call compute_mbe(polymers, fragment_count, max_level, results, mbe_result)

      call check(error, mbe_result%has_dipole, "the MBE dipole was not computed")
      if (allocated(error)) return

      expected = results(1)%dipole + results(2)%dipole &
                 + (results(5)%dipole - results(3)%dipole - results(4)%dipole)

      do i = 1, 3
         call check(error, mbe_result%dipole(i), expected(i), thr=1.0e-12_dp, &
                    message="the counterpoise dipole does not match the expansion; "// &
                    "a ghosted row was summed, or a subset was looked up in the wrong basis")
         if (allocated(error)) return
      end do

      ! The energy is assembled by the path this mirrors, so it pins the same
      ! two rules independently.
      call check(error, mbe_result%total_energy, &
                 -10.0_dp - 15.0_dp + (-25.5_dp + 10.1_dp + 15.2_dp), thr=1.0e-12_dp, &
                 message="the counterpoise energy changed, so the fixture is wrong "// &
                 "rather than the dipole")
   end subroutine test_vmfc_dipole

   subroutine test_mbe_monomers_only(error)
      !! Test MBE energy with monomers only (nlevel=1)
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), allocatable :: results(:)
      type(mbe_result_t) :: mbe_result
      integer, allocatable :: polymers(:, :)
      integer(int64) :: fragment_count
      integer :: max_level

      ! Three monomers with energies -10.0, -15.0, -20.0
      fragment_count = 3
      max_level = 1
      allocate (polymers(fragment_count, max_level))
      allocate (results(fragment_count))

      ! Monomers in order: [1], [2], [3]
      polymers(1, 1) = 1
      polymers(2, 1) = 2
      polymers(3, 1) = 3

      ! Set monomer energies
      results(1)%energy%scf = -10.0_dp
      results(1)%has_energy = .true.
      results(2)%energy%scf = -15.0_dp
      results(2)%has_energy = .true.
      results(3)%energy%scf = -20.0_dp
      results(3)%has_energy = .true.

      ! Compute MBE energy
      call compute_mbe(polymers, fragment_count, max_level, results, mbe_result)

      ! Total should be -10 + (-15) + (-20) = -45
      call check(error, mbe_result%total_energy, -45.0_dp, thr=1.0e-10_dp, &
                 message="MBE monomers only should equal sum of monomer energies")
      call mbe_result%destroy()
      if (allocated(error)) return

      if (allocated(results)) deallocate (results)
      if (allocated(polymers)) deallocate (polymers)

   end subroutine test_mbe_monomers_only

   subroutine test_mbe_simple_dimer(error)
      !! Test MBE with 2 monomers and 1 dimer (nlevel=2)
      !! E_total = E(1) + E(2) + deltaE(1,2)
      !! where deltaE(1,2) = E(1,2) - E(1) - E(2)
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), allocatable :: results(:)
      type(mbe_result_t) :: mbe_result
      integer, allocatable :: polymers(:, :)
      integer(int64) :: fragment_count
      integer :: max_level

      ! 2 monomers + 1 dimer
      fragment_count = 3
      max_level = 2
      allocate (polymers(fragment_count, max_level))
      allocate (results(fragment_count))
      polymers = 0

      ! Fragment 1: monomer [1]
      polymers(1, 1) = 1
      results(1)%energy%scf = -10.0_dp
      results(1)%has_energy = .true.

      ! Fragment 2: monomer [2]
      polymers(2, 1) = 2
      results(2)%energy%scf = -15.0_dp
      results(2)%has_energy = .true.

      ! Fragment 3: dimer [1, 2]
      polymers(3, 1) = 1
      polymers(3, 2) = 2
      results(3)%energy%scf = -26.0_dp  ! Slightly more stable than sum
      results(3)%has_energy = .true.

      ! Compute MBE energy
      call compute_mbe(polymers, fragment_count, max_level, results, mbe_result)

      ! deltaE(1,2) = E(1,2) - E(1) - E(2) = -26 - (-10) - (-15) = -26 + 10 + 15 = -1
      ! Total = E(1) + E(2) + deltaE(1,2) = -10 + (-15) + (-1) = -26
      call check(error, mbe_result%total_energy, -26.0_dp, thr=1.0e-10_dp, &
                 message="MBE with simple dimer")
      call mbe_result%destroy()
      if (allocated(error)) return

      if (allocated(results)) deallocate (results)
      if (allocated(polymers)) deallocate (polymers)

   end subroutine test_mbe_simple_dimer

   subroutine test_mbe_sorted_order(error)
      !! Test MBE with fragments in size order (monomers first)
      !! 3 monomers + 3 dimers + 1 trimer
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), allocatable :: results(:)
      type(mbe_result_t) :: mbe_result
      integer, allocatable :: polymers(:, :)
      integer(int64) :: fragment_count
      integer :: max_level

      fragment_count = 7  ! 3 monomers + 3 dimers + 1 trimer
      max_level = 3
      allocate (polymers(fragment_count, max_level))
      allocate (results(fragment_count))
      polymers = 0

      ! Monomers: [1], [2], [3]
      polymers(1, 1) = 1
      polymers(2, 1) = 2
      polymers(3, 1) = 3
      results(1)%energy%scf = -10.0_dp
      results(2)%energy%scf = -12.0_dp
      results(3)%energy%scf = -11.0_dp
      results(1:3)%has_energy = .true.

      ! Dimers: [1,2], [1,3], [2,3]
      polymers(4, 1:2) = [1, 2]
      polymers(5, 1:2) = [1, 3]
      polymers(6, 1:2) = [2, 3]
      results(4)%energy%scf = -22.5_dp  ! E(1) + E(2) = -22, so delta = -0.5
      results(5)%energy%scf = -21.3_dp  ! E(1) + E(3) = -21, so delta = -0.3
      results(6)%energy%scf = -23.4_dp  ! E(2) + E(3) = -23, so delta = -0.4
      results(4:6)%has_energy = .true.

      ! Trimer: [1,2,3]
      polymers(7, 1:3) = [1, 2, 3]
      results(7)%energy%scf = -33.5_dp
      results(7)%has_energy = .true.

      call compute_mbe(polymers, fragment_count, max_level, results, mbe_result)

      ! 1-body: -10 + -12 + -11 = -33
      ! 2-body deltas: -0.5 + -0.3 + -0.4 = -1.2
      ! 3-body delta: E(1,2,3) - E(1) - E(2) - E(3) - delta(1,2) - delta(1,3) - delta(2,3)
      !             = -33.5 - (-10) - (-12) - (-11) - (-0.5) - (-0.3) - (-0.4)
      !             = -33.5 + 10 + 12 + 11 + 0.5 + 0.3 + 0.4 = 0.7
      ! Total = -33 + (-1.2) + 0.7 = -33.5
      call check(error, mbe_result%total_energy, -33.5_dp, thr=1.0e-10_dp, &
                 message="MBE with fragments in sorted order")
      call mbe_result%destroy()
      if (allocated(error)) return

      if (allocated(results)) deallocate (results)
      if (allocated(polymers)) deallocate (polymers)

   end subroutine test_mbe_sorted_order

   subroutine test_mbe_reverse_order(error)
      !! Test MBE with fragments in REVERSE size order (trimer first, monomers last)
      !! This tests that internal sorting works correctly
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), allocatable :: results(:)
      type(mbe_result_t) :: mbe_result
      integer, allocatable :: polymers(:, :)
      integer(int64) :: fragment_count
      integer :: max_level

      fragment_count = 7  ! 3 monomers + 3 dimers + 1 trimer
      max_level = 3
      allocate (polymers(fragment_count, max_level))
      allocate (results(fragment_count))
      polymers = 0

      ! REVERSE ORDER: Trimer first
      polymers(1, 1:3) = [1, 2, 3]
      results(1)%energy%scf = -33.5_dp
      results(1)%has_energy = .true.

      ! Dimers: [1,2], [1,3], [2,3]
      polymers(2, 1:2) = [1, 2]
      polymers(3, 1:2) = [1, 3]
      polymers(4, 1:2) = [2, 3]
      results(2)%energy%scf = -22.5_dp
      results(3)%energy%scf = -21.3_dp
      results(4)%energy%scf = -23.4_dp
      results(2:4)%has_energy = .true.

      ! Monomers last: [1], [2], [3]
      polymers(5, 1) = 1
      polymers(6, 1) = 2
      polymers(7, 1) = 3
      results(5)%energy%scf = -10.0_dp
      results(6)%energy%scf = -12.0_dp
      results(7)%energy%scf = -11.0_dp
      results(5:7)%has_energy = .true.

      call compute_mbe(polymers, fragment_count, max_level, results, mbe_result)

      ! Should get same answer as sorted order test
      call check(error, mbe_result%total_energy, -33.5_dp, thr=1.0e-10_dp, &
                 message="MBE with fragments in reverse order (tests internal sorting)")
      call mbe_result%destroy()
      if (allocated(error)) return

      if (allocated(results)) deallocate (results)
      if (allocated(polymers)) deallocate (polymers)

   end subroutine test_mbe_reverse_order

   subroutine test_mbe_random_order(error)
      !! Test MBE with fragments in mixed/random order
      !! Order: [dimer], [monomer], [trimer], [monomer], [dimer], [dimer], [monomer]
      type(error_type), allocatable, intent(out) :: error
      type(calculation_result_t), allocatable :: results(:)
      type(mbe_result_t) :: mbe_result
      integer, allocatable :: polymers(:, :)
      integer(int64) :: fragment_count
      integer :: max_level

      fragment_count = 7
      max_level = 3
      allocate (polymers(fragment_count, max_level))
      allocate (results(fragment_count))
      polymers = 0

      ! Mixed order
      ! Fragment 1: dimer [1,2]
      polymers(1, 1:2) = [1, 2]
      results(1)%energy%scf = -22.5_dp
      results(1)%has_energy = .true.

      ! Fragment 2: monomer [1]
      polymers(2, 1) = 1
      results(2)%energy%scf = -10.0_dp
      results(2)%has_energy = .true.

      ! Fragment 3: trimer [1,2,3]
      polymers(3, 1:3) = [1, 2, 3]
      results(3)%energy%scf = -33.5_dp
      results(3)%has_energy = .true.

      ! Fragment 4: monomer [2]
      polymers(4, 1) = 2
      results(4)%energy%scf = -12.0_dp
      results(4)%has_energy = .true.

      ! Fragment 5: dimer [1,3]
      polymers(5, 1:2) = [1, 3]
      results(5)%energy%scf = -21.3_dp
      results(5)%has_energy = .true.

      ! Fragment 6: dimer [2,3]
      polymers(6, 1:2) = [2, 3]
      results(6)%energy%scf = -23.4_dp
      results(6)%has_energy = .true.

      ! Fragment 7: monomer [3]
      polymers(7, 1) = 3
      results(7)%energy%scf = -11.0_dp
      results(7)%has_energy = .true.

      call compute_mbe(polymers, fragment_count, max_level, results, mbe_result)

      ! Should get same answer regardless of input order
      call check(error, mbe_result%total_energy, -33.5_dp, thr=1.0e-10_dp, &
                 message="MBE with fragments in random/mixed order")
      call mbe_result%destroy()
      if (allocated(error)) return

      if (allocated(results)) deallocate (results)
      if (allocated(polymers)) deallocate (polymers)

   end subroutine test_mbe_random_order

   subroutine test_unconverged_carry_their_monomers(error)
      !! The failed fragments come back with what they were built from
      !!
      !! An identifier on its own cannot be re-run: a dimer is only
      !! reconstructible if you know which two monomers it was. So the pair has
      !! to survive together, and in the same order as the fragments they came
      !! from -- a follow-up job built from a list that is right in aggregate
      !! and shuffled in detail would re-run the wrong dimers and look fine.
      !!
      !! Five fragments here: three monomers and two dimers, of which the second
      !! monomer and the second dimer failed. The third monomer is left
      !! `SCF_UNKNOWN` **deliberately and in the same array as real failures**:
      !! a case built only from unknowns never reaches the selection at all,
      !! because there is nothing to collect and the routine returns early. It
      !! takes a mixture to find out whether "did not converge" and "did not
      !! say" are being told apart.
      type(error_type), allocatable, intent(out) :: error

      integer :: status(5), polymers(5, 2)
      integer(int64), allocatable :: ids(:)
      integer, allocatable :: monomers(:, :)

      polymers = 0
      polymers(1, 1) = 1
      polymers(2, 1) = 2
      polymers(3, 1) = 3
      polymers(4, :) = [1, 2]
      polymers(5, :) = [2, 3]

      status = SCF_CONVERGED
      status(2) = SCF_NOT_CONVERGED
      status(3) = SCF_UNKNOWN
      status(5) = SCF_NOT_CONVERGED

      call collect_unconverged(status, polymers, 5_int64, ids, monomers)

      call check(error, size(ids), 2, &
                 "two fragments failed; the silent one is not a third")
      if (allocated(error)) return
      call check(error, int(ids(1)), 2, "the first failure is fragment 2")
      if (allocated(error)) return
      call check(error, int(ids(2)), 5, "the second failure is fragment 5")
      if (allocated(error)) return

      ! The monomer keeps its single member and its padding; the dimer keeps
      ! both, in order.
      call check(error, monomers(1, 1), 2, "fragment 2 is monomer 2")
      if (allocated(error)) return
      call check(error, monomers(1, 2), 0, "and nothing else")
      if (allocated(error)) return
      call check(error, monomers(2, 1), 2, "fragment 5 is the dimer of 2")
      if (allocated(error)) return
      call check(error, monomers(2, 2), 3, "and 3")
   end subroutine test_unconverged_carry_their_monomers

   subroutine test_unconverged_ignores_silent_methods(error)
      !! A method that never reports convergence has not failed everywhere
      !!
      !! `SCF_UNKNOWN` is what a method that does not report leaves on every
      !! fragment. Listing those would fill the follow-up job with the entire
      !! calculation, in exactly the runs this exists to help with, so only
      !! genuine failures are collected, and nothing is allocated at all --
      !! because "nothing failed" and "nobody said" are different claims and
      !! the caller distinguishes them by allocation. An empty list is what a
      !! reporting method returns when everything converged, so returning one
      !! here too would make the two indistinguishable and the JSON section
      !! would assert success on a method that never said anything.
      type(error_type), allocatable, intent(out) :: error

      integer :: status(3), polymers(3, 1)
      integer(int64), allocatable :: ids(:)
      integer, allocatable :: monomers(:, :)

      polymers(:, 1) = [1, 2, 3]
      status = SCF_UNKNOWN

      call collect_unconverged(status, polymers, 3_int64, ids, monomers)
      call check(error,.not. allocated(ids), &
                 "a silent method has not failed everywhere, and has not said so either")
      if (allocated(error)) return
      call check(error,.not. allocated(monomers), &
                 "the composition list goes the same way as the identifiers")
   end subroutine test_unconverged_ignores_silent_methods

   subroutine test_converged_run_says_so(error)
      !! A method that reported and lost nothing returns an empty list, not nothing
      !!
      !! The other half of the silent-method case, and the reason that one
      !! cannot simply return an empty list too. Here the list is allocated and
      !! zero-length, which the writer turns into `"count": 0` -- an assertion
      !! that every fragment converged. Absence would say only that nobody
      !! asked. Both are pinned because a change to either collapses the pair,
      !! and the collapse is invisible until somebody reads the JSON.
      type(error_type), allocatable, intent(out) :: error

      integer :: status(3), polymers(3, 1)
      integer(int64), allocatable :: ids(:)
      integer, allocatable :: monomers(:, :)

      polymers(:, 1) = [1, 2, 3]
      status = SCF_CONVERGED

      call collect_unconverged(status, polymers, 3_int64, ids, monomers)
      call check(error, allocated(ids), "a reporting method returns a list")
      if (allocated(error)) return
      call check(error, size(ids), 0, "and nothing is on it")
      if (allocated(error)) return
      call check(error, allocated(monomers), "the composition list comes back too")
   end subroutine test_converged_run_says_so

   subroutine test_failures_name_their_culprit(error)
      !! The monomer that keeps turning up, and what the failures cost
      !!
      !! One bad monomer drags down every fragment it belongs to, so failures
      !! cluster on their cause. Here monomer 2 is in all three failed
      !! fragments and monomers 3 and 4 in one each -- which is one problem
      !! wearing three disguises, and the ordering is what says so. A list in
      !! fragment order would put monomer 2 first by accident, so the fixture
      !! deliberately does not: the failures are ordered such that monomer 5
      !! is seen first and must still come last.
      type(error_type), allocatable, intent(out) :: error

      integer(int64) :: ids(3)
      integer :: monomers(3, 2)
      real(dp) :: delta_energies(8)
      real(dp), allocatable :: deltas(:)
      integer, allocatable :: culprits(:)
      integer(int64), allocatable :: counts(:)
      integer :: i

      ids = [3_int64, 5_int64, 8_int64]
      monomers(1, :) = [5, 2]     ! monomer 5 appears first...
      monomers(2, :) = [2, 3]
      monomers(3, :) = [2, 4]

      do i = 1, 8
         delta_energies(i) = -0.001_dp*real(i, dp)
      end do

      call score_unconverged(ids, monomers, delta_energies, deltas, culprits, counts)

      ! The contributions are gathered by fragment index, not by position.
      call check(error, deltas(1), -0.003_dp, thr=1.0e-12_dp, &
                 more="the first failure's contribution came from the wrong fragment")
      if (allocated(error)) return
      call check(error, deltas(3), -0.008_dp, thr=1.0e-12_dp, &
                 more="the last failure's contribution came from the wrong fragment")
      if (allocated(error)) return

      call check(error, size(culprits), 4, "four distinct monomers are involved")
      if (allocated(error)) return
      call check(error, culprits(1), 2, &
                 "the monomer in every failure should be named first")
      if (allocated(error)) return
      call check(error, int(counts(1)), 3, "and it is in three of them")
      if (allocated(error)) return
      ! ...but must not be ranked first merely for having been seen first.
      call check(error, int(counts(size(counts))), 1, &
                 "the monomers in one failure each should rank last")
   end subroutine test_failures_name_their_culprit

end module test_mqc_mbe

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mbe, only: collect_mqc_mbe_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_mbe", collect_mqc_mbe_tests) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if

end program tester
