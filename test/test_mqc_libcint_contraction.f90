module test_mqc_libcint_contraction
   !! Pins that general contractions reach libcint as one shell, not several.
   !!
   !! The reader emits one shell per coefficient column, because that is how a
   !! Basis Set Exchange entry is written. libcint can take the columns of a
   !! general contraction as a single shell with `NCTR_OF` set to the column
   !! count, evaluate the shared primitives once, and contract them into every
   !! column. Splitting them repeats the primitive work per column, which on a
   !! water hexamer in cc-pVDZ was 3.4x the whole SCF.
   !!
   !! Two properties have to hold together, and only together are they worth
   !! anything:
   !!
   !!   * fewer shells -- otherwise nothing merged and the optimisation is not
   !!     happening at all, which no energy would reveal because the answer is
   !!     right either way;
   !!   * the same basis functions, in the same order -- otherwise the AO basis
   !!     has been silently permuted, and every matrix built on it with it.
   !!
   !! The second is the dangerous one. libcint lays a shell's functions out with
   !! the contraction index outermost, so merging columns that were already
   !! adjacent reproduces the previous order exactly; merging across a gap would
   !! not. The overlap matrix is what tests it here -- a permutation leaves the
   !! function count untouched and moves S, so a count alone would pass while
   !! the basis was scrambled.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use mqc_cgto, only: molecular_basis_type
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_libcint_contraction_tests

   character(len=*), parameter :: CC_PVDZ = "../basis_sets/cc-pvdz.json"
   character(len=*), parameter :: POPLE = "../basis_sets/6-31g.json"
   character(len=*), parameter :: STO_3G = "../basis_sets/sto-3g.json"

   !> Water in Bohr, the same geometry the other libcint suites use.
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   real(dp), parameter :: WATER_COORDS(3, 3) = reshape([ &
                                                       0.0_dp, 0.0_dp, 0.225374_dp, &
                                                       0.0_dp, 1.442316_dp, -0.901497_dp, &
                                                       0.0_dp, -1.442316_dp, -0.901497_dp], [3, 3])

contains

   subroutine collect_mqc_libcint_contraction_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("contraction_merges_general_contractions", test_merges), &
                  new_unittest("contraction_preserves_function_count", test_nao), &
                  new_unittest("contraction_preserves_the_overlap_matrix", test_overlap), &
                  new_unittest("contraction_leaves_sp_shells_alone", test_sp_untouched), &
                  new_unittest("contraction_leaves_uncontracted_sets_alone", test_sto3g) &
                  ]
   end subroutine collect_mqc_libcint_contraction_tests

   subroutine water(basis_file, mol, reader_shells, err)
      !! Build water and report how many shells the reader produced
      character(len=*), intent(in) :: basis_file
      type(libcint_molecule_t), intent(out) :: mol
      integer, intent(out) :: reader_shells
      type(error_t), intent(inout) :: err

      type(molecular_basis_type) :: basis
      character(len=2) :: symbols(3)
      integer :: ia

      symbols = ["O ", "H ", "H "]
      call build_molecular_basis_json(basis_file, symbols, basis, err)
      if (err%has_error()) return

      reader_shells = 0
      do ia = 1, basis%nelements
         reader_shells = reader_shells + basis%elements(ia)%nshells
      end do

      call mol%build(WATER_Z, WATER_COORDS, basis, err)
   end subroutine water

   subroutine test_merges(error)
      !! cc-pVDZ is a general contraction, so libcint must see fewer shells
      !! than the reader produced
      !!
      !! Measured rather than derived: the reader gives water twelve shells and
      !! libcint is handed seven. Both numbers are pinned, because either one
      !! alone can be right while the optimisation is not happening -- a stale
      !! reader count would make an unmerged basis look merged.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      integer :: reader_shells

      call water(CC_PVDZ, mol, reader_shells, err)
      call check(error,.not. err%has_error(), "cc-pVDZ water must build: "//err%get_full_trace())
      if (allocated(error)) return

      call check(error, mol%nbas < reader_shells, &
                 "a general contraction must reach libcint as fewer shells than columns")
      if (allocated(error)) return
      call check(error, reader_shells, 12)
      if (allocated(error)) return
      call check(error, mol%nbas, 7)
   end subroutine test_merges

   subroutine test_nao(error)
      !! Merging changes the shell count and must not change the basis
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      integer :: reader_shells

      call water(CC_PVDZ, mol, reader_shells, err)
      call check(error,.not. err%has_error(), "cc-pVDZ water must build: "//err%get_full_trace())
      if (allocated(error)) return
      ! H(2s 1p) x2 + O(3s 2p 1d) = 2*(2 + 3) + (3 + 6 + 5) = 24
      call check(error, mol%nao, 24)
   end subroutine test_nao

   subroutine test_overlap(error)
      !! The property a function count cannot see
      !!
      !! If merging reordered the basis, S would be a permutation of the right
      !! matrix -- same size, same eigenvalues, wrong entries. Checking the
      !! diagonal is one and that S is symmetric catches a normalisation error;
      !! checking a specific off-diagonal against the value it has when every
      !! shell is separate is what catches a permutation.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: s(:, :)
      integer :: reader_shells, i
      real(dp) :: worst

      call water(CC_PVDZ, mol, reader_shells, err)
      call check(error,.not. err%has_error(), "cc-pVDZ water must build: "//err%get_full_trace())
      if (allocated(error)) return

      call mol%overlap(s)

      ! Every basis function is normalised, exactly, whatever shell it came from.
      !
      ! This check used to allow 1e-5, on the reasoning that a contracted
      ! function's norm is only as good as the ten significant figures Basis Set
      ! Exchange publishes. That reasoning was wrong, and it argued a real bug
      ! into looking like a property of the input: ten digits buys ~1e-9, not the
      ! ~1e-6 that was actually being seen. The cause was that the *contraction*
      ! normalisation was never applied at all -- only the per-primitive one that
      ! libcint requires -- so `<chi|chi>` came out 1.000001 for every contracted
      ! shell and exactly 1 for every single primitive, which is the fingerprint
      ! that gave it away.
      !
      ! Nothing in an energy could have caught it: an SCF energy is invariant to
      ! the normalisation of a basis function, since scaling one does not change
      ! the space the basis spans. It surfaced only when multipole integrals gave
      ! the first per-AO quantity anyone compared against PySCF.
      worst = 0.0_dp
      do i = 1, mol%nao
         worst = max(worst, abs(s(i, i) - 1.0_dp))
      end do
      call check(error, worst < 1.0e-12_dp, "the overlap diagonal must be one")
      if (allocated(error)) return

      call check(error, maxval(abs(s - transpose(s))) < 1.0e-14_dp, &
                 "the overlap matrix must be symmetric")
      if (allocated(error)) return

      ! Functions 1-3 are oxygen's s contractions, the block merging collapses
      ! into one shell. These two values are **PySCF's**, on this geometry and
      ! this basis, which is better provenance than what was here before: they
      ! used to be measured from our own split build, so when the contraction
      ! normalisation was fixed they were the only thing that failed -- a
      ! reference that tracks our own output cannot tell a fix from a
      ! regression. An external one can.
      !
      ! S(2,3) is the one that matters. It is 0.93, so a permutation moves it
      ! somewhere obvious, whereas S(1,2) is 1.1e-6 -- the two contractions are
      ! nearly orthogonal, being close to atomic 1s and 2s -- and a great many
      ! wrong matrices also have something near zero there.
      call check(error, s(2, 3), 0.9333237202898537_dp, thr=1.0e-8_dp)
      if (allocated(error)) return
      call check(error, s(1, 3), 0.19189620099353646_dp, thr=1.0e-8_dp)
   end subroutine test_overlap

   subroutine test_sp_untouched(error)
      !! An SP shell shares exponents but not angular momentum
      !!
      !! 6-31G writes its valence as one entry with two angular momenta and two
      !! coefficient columns. Those columns are an s and a p over shared
      !! primitives -- not a general contraction, and merging them into one
      !! shell would be nonsense. The angular momentum test is what stops it.
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      integer :: reader_shells

      call water(POPLE, mol, reader_shells, err)
      call check(error,.not. err%has_error(), "6-31G water must build: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, mol%nbas, reader_shells, &
                 "SP columns carry different angular momenta and must not merge")
      if (allocated(error)) return
      ! H(2s) x2 + O(3s 2p) = 2*2 + 3 + 6 = 13
      call check(error, mol%nao, 13)
      if (allocated(error)) return
      call check(error, mol%nbas, 9)
   end subroutine test_sp_untouched

   subroutine test_sto3g(error)
      !! A segmented set has one column per shell, so nothing merges
      type(error_type), allocatable, intent(out) :: error
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      integer :: reader_shells

      call water(STO_3G, mol, reader_shells, err)
      call check(error,.not. err%has_error(), "STO-3G water must build: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, mol%nbas, reader_shells, &
                 "a segmented basis has nothing to merge")
      if (allocated(error)) return
      call check(error, mol%nao, 7)
   end subroutine test_sto3g

end module test_mqc_libcint_contraction

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_libcint_contraction, only: collect_mqc_libcint_contraction_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_libcint_contraction", collect_mqc_libcint_contraction_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
