!! The Mulliken arithmetic, on inputs small enough to check by hand
!!
!! This is the module both integrals backends now call, so a mistake here is a
!! mistake on CPU and GPU at once, and the two would agree with each other
!! while disagreeing with everyone else. That is the case worth testing: the
!! backend-specific parts -- how each arrives at the AO-to-atom map -- are
!! checked where they live, and what is checked here is the part they share.
!!
!! The matrices are written out rather than computed, so every expected number
!! below is arithmetic a reader can redo.
module test_mqc_population_analysis
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mqc_error, only: error_t
   use mqc_population_analysis, only: ao_owner_from_counts, mulliken_atomic_charges, &
                                      mulliken_atomic_spin_populations
   implicit none
   private

   public :: collect_population_analysis

   real(dp), parameter :: THR = 1.0e-12_dp

contains

   subroutine collect_population_analysis(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("ao owner map expands a count per atom", test_owner_map), &
                  new_unittest("ao owner map tolerates an atom with no functions", test_owner_bare), &
                  new_unittest("mulliken charges on an orthonormal basis", test_orthonormal), &
                  new_unittest("mulliken charges sum to the molecular charge", test_sum_rule), &
                  new_unittest("mulliken splits an overlap evenly", test_even_split), &
                  new_unittest("spin populations sum to the number of unpaired electrons", test_spin), &
                  new_unittest("a mis-sized density is refused", test_size_mismatch) &
                  ]
   end subroutine collect_population_analysis

   subroutine test_owner_map(error)
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: owner(:)

      call ao_owner_from_counts([2, 1, 3], owner)
      call check(error, size(owner), 6)
      if (allocated(error)) return
      call check(error, all(owner == [1, 1, 2, 3, 3, 3]), "owner map is not grouped by atom")
   end subroutine test_owner_map

   subroutine test_owner_bare(error)
      !! A bare charge -- a point with no basis on it -- contributes no columns
      !!
      !! It must still occupy an atom index, which is the reason
      !! `mulliken_atomic_spin_populations` takes `natm` rather than deriving
      !! it from the map: a trailing empty atom is invisible in `owner`.
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: owner(:)

      call ao_owner_from_counts([2, 0, 1], owner)
      call check(error, all(owner == [1, 1, 3]), "an empty atom shifted the map")
   end subroutine test_owner_bare

   subroutine test_orthonormal(error)
      !! With S = I the population is the diagonal of D, so the answer is
      !! nuclear charge minus the electrons already sitting on each atom.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: density(3, 3), overlap(3, 3)
      real(dp), allocatable :: q(:)
      integer, allocatable :: owner(:)
      type(error_t) :: err
      integer :: i

      overlap = 0.0_dp
      do i = 1, 3
         overlap(i, i) = 1.0_dp
      end do
      density = 0.0_dp
      density(1, 1) = 2.0_dp     !! two electrons on atom 1's first function
      density(2, 2) = 1.0_dp     !! one on its second
      density(3, 3) = 4.0_dp     !! four on atom 2

      call ao_owner_from_counts([2, 1], owner)
      call mulliken_atomic_charges(owner, [4.0_dp, 5.0_dp], density, overlap, q, err)
      call check(error, err%has_error(), .false.)
      if (allocated(error)) return
      call check(error, q(1), 1.0_dp, thr=THR)      !! 4 - (2 + 1)
      if (allocated(error)) return
      call check(error, q(2), 1.0_dp, thr=THR)      !! 5 - 4
   end subroutine test_orthonormal

   subroutine test_sum_rule(error)
      !! sum_A q_A = sum_A Z_A - tr(D S), whatever the overlap is
      !!
      !! The identity Mulliken cannot violate however arbitrary the partition,
      !! and the one that catches a function counted on the wrong atom only if
      !! it is *also* counted twice -- which is why the split is checked
      !! separately below.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: density(3, 3), overlap(3, 3)
      real(dp), allocatable :: q(:)
      integer, allocatable :: owner(:)
      type(error_t) :: err
      real(dp) :: trace
      integer :: i, j

      overlap = reshape([1.0_dp, 0.31_dp, 0.12_dp, &
                         0.31_dp, 1.0_dp, 0.24_dp, &
                         0.12_dp, 0.24_dp, 1.0_dp], [3, 3])
      density = reshape([1.7_dp, 0.42_dp, -0.13_dp, &
                         0.42_dp, 0.9_dp, 0.28_dp, &
                         -0.13_dp, 0.28_dp, 1.4_dp], [3, 3])

      call ao_owner_from_counts([2, 1], owner)
      call mulliken_atomic_charges(owner, [3.0_dp, 2.0_dp], density, overlap, q, err)
      call check(error, err%has_error(), .false.)
      if (allocated(error)) return

      trace = 0.0_dp
      do i = 1, 3
         do j = 1, 3
            trace = trace + density(i, j)*overlap(j, i)
         end do
      end do
      call check(error, sum(q), 5.0_dp - trace, thr=THR)
   end subroutine test_sum_rule

   subroutine test_even_split(error)
      !! The halving that *is* Mulliken
      !!
      !! One function per atom, a density whose only off-diagonal is the bond,
      !! and a symmetric overlap: the shared population must land half on each
      !! atom. An implementation that took `(D S)` the other way round, or that
      !! summed a row instead of the diagonal, gives the same total and the
      !! wrong split -- so the total is not enough to catch it.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: density(2, 2), overlap(2, 2)
      real(dp), allocatable :: q(:)
      integer, allocatable :: owner(:)
      type(error_t) :: err

      overlap = reshape([1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2, 2])
      density = reshape([1.0_dp, 0.4_dp, 0.4_dp, 1.0_dp], [2, 2])

      call ao_owner_from_counts([1, 1], owner)
      call mulliken_atomic_charges(owner, [1.0_dp, 1.0_dp], density, overlap, q, err)
      call check(error, err%has_error(), .false.)
      if (allocated(error)) return
      ! identical atoms, identical environment: no charge separation at all
      call check(error, q(1), q(2), thr=THR)
      if (allocated(error)) return
      ! and each carries 1 - (1*1 + 0.4*0.5) = -0.2
      call check(error, q(1), -0.2_dp, thr=THR)
   end subroutine test_even_split

   subroutine test_spin(error)
      !! Populations sum to tr((D_a - D_b) S), the unpaired electron count
      !!
      !! No nuclear charge enters: a nucleus carries no spin, so this sums to
      !! n_alpha - n_beta rather than to the molecular charge.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: spin_density(2, 2), overlap(2, 2)
      real(dp), allocatable :: p(:)
      integer, allocatable :: owner(:)
      type(error_t) :: err

      overlap = reshape([1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2, 2])
      ! one unpaired electron, entirely on atom 1
      spin_density = 0.0_dp
      spin_density(1, 1) = 1.0_dp

      call ao_owner_from_counts([1, 1], owner)
      call mulliken_atomic_spin_populations(owner, 2, spin_density, overlap, p, err)
      call check(error, err%has_error(), .false.)
      if (allocated(error)) return
      call check(error, sum(p), 1.0_dp, thr=THR)
      if (allocated(error)) return
      call check(error, p(1), 1.0_dp, thr=THR)
      if (allocated(error)) return
      call check(error, p(2), 0.0_dp, thr=THR)
   end subroutine test_spin

   subroutine test_size_mismatch(error)
      !! A density that is not the size of the basis is an error, not a crash
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: density(2, 2), overlap(2, 2)
      real(dp), allocatable :: q(:)
      integer, allocatable :: owner(:)
      type(error_t) :: err

      density = 0.0_dp
      overlap = 0.0_dp
      call ao_owner_from_counts([2, 1], owner)   !! three functions, two supplied
      call mulliken_atomic_charges(owner, [1.0_dp, 1.0_dp], density, overlap, q, err)
      call check(error, err%has_error(), .true.)
   end subroutine test_size_mismatch

end module test_mqc_population_analysis

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_population_analysis, only: collect_population_analysis
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_population_analysis", collect_population_analysis)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
