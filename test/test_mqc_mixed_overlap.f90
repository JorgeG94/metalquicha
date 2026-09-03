!! The overlap between two different bases on the same atoms
module test_mqc_mixed_overlap
   !! `< chi_mu | phi_k >` between the orbital basis an SCF ran in and a second
   !! basis on the same nuclei. The quasi-atomic construction is built entirely
   !! out of singular value decompositions of this matrix, so it is the one
   !! genuinely new integral the analysis needs.
   !!
   !! The load-bearing test is the first one: asking for the overlap of a basis
   !! with *itself* must reproduce the ordinary overlap matrix. That is a real
   !! check rather than a tautology because the two answers come from different
   !! code -- one walks a single shell table, the other walks a concatenated one
   !! with the second basis's `env` offsets shifted -- and the shift is exactly
   !! the thing that would be wrong.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_cgto, only: molecular_basis_type
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use mqc_aambs, only: aambs_file, aambs_dimensions, aambs_dimensions_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule, &
                                mixed_basis_overlap
   implicit none
   private

   public :: collect_mqc_mixed_overlap_tests

   integer, parameter :: N_DIM = 3
   !! Water, the geometry the exchange-correlation kernel checks use, in Bohr.
   real(dp), parameter :: WATER(N_DIM, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.0_dp, 1.4308_dp, 1.1078_dp], [N_DIM, 3])
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]

contains

   subroutine collect_mqc_mixed_overlap_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("self_overlap_reproduces_S", test_self_overlap), &
                  new_unittest("aambs_shells_are_normalized", test_aambs_normalization), &
                  new_unittest("shape_matches_the_counting", test_shape), &
                  new_unittest("angular_form_mismatch_is_refused", test_refusal) &
                  ]
   end subroutine collect_mqc_mixed_overlap_tests

   subroutine build_aambs(mol, error, ok)
      !! Water in the accurate atomic minimal basis
      type(czt_molecule_t), intent(out) :: mol
      type(error_t), intent(inout) :: error
      logical, intent(out) :: ok

      type(molecular_basis_type) :: basis
      character(len=:), allocatable :: path

      ok = .false.
      call aambs_file(path, error)
      if (error%has_error()) return
      call build_molecular_basis_json(path, WATER_SYM, basis, error)
      if (error%has_error()) return
      call mol%build(WATER_Z, WATER, basis, error)
      call basis%destroy()
      ok = .not. error%has_error()
   end subroutine build_aambs

   subroutine test_self_overlap(error)
      !! The overlap of a basis with itself is the overlap matrix
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: mixed(:, :), s(:, :)

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "cc-pvdz", mol, err)
      call check(error,.not. err%has_error(), "water in cc-pVDZ should build")
      if (allocated(error)) return

      call mol%overlap(s)
      call mixed_basis_overlap(mol, mol, mixed, err)
      call check(error,.not. err%has_error(), "a basis with itself should be allowed")
      if (allocated(error)) return

      call check(error, size(mixed, 1), size(s, 1), "row count should match")
      if (allocated(error)) return
      call check(error, size(mixed, 2), size(s, 2), "column count should match")
      if (allocated(error)) return

      ! Machine precision, not a tolerance: these are the same integrals from
      ! the same library, so any difference at all is a bookkeeping error in the
      ! concatenated shell table rather than numerical noise.
      call check(error, maxval(abs(mixed - s)) < 1.0e-13_dp, &
                 "the mixed overlap of a basis with itself must reproduce S exactly")
      if (allocated(error)) return
      call mol%destroy()
   end subroutine test_self_overlap

   subroutine test_aambs_normalization(error)
      !! Every minimal-basis function is normalized
      !!
      !! The extraction verified this against the raw contraction coefficients;
      !! this verifies it survives the whole path into libcint -- the reader,
      !! the primitive normalization libcint wants, and the spherical transform.
      !! A diagonal that is not 1 means the free-atom orbitals are being scaled
      !! somewhere between the JSON and the integral.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      real(dp), allocatable :: s(:, :)
      logical :: ok
      integer :: i
      real(dp) :: worst

      call build_aambs(mol, err, ok)
      call check(error, ok, "water in the atomic minimal basis should build")
      if (allocated(error)) return

      call mol%overlap(s)
      worst = 0.0_dp
      do i = 1, size(s, 1)
         worst = max(worst, abs(s(i, i) - 1.0_dp))
      end do
      call check(error, worst < 1.0e-6_dp, &
                 "the minimal-basis functions should come out normalized")
      if (allocated(error)) return
      call mol%destroy()
   end subroutine test_aambs_normalization

   subroutine test_shape(error)
      !! The overlap is (n_ao, n_mbs), with n_mbs what the counting predicts
      !!
      !! Two independent routes to the same number: the minimal-basis dimension
      !! summed from the per-element tables, and the actual size of the basis
      !! libcint built from the same file. They are only equal if the shell data
      !! and the orbital counts describe the same basis, which is the thing most
      !! likely to drift if either is edited.
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: orb, aambs
      type(aambs_dimensions_t) :: dims
      type(error_t) :: err
      real(dp), allocatable :: mixed(:, :)
      logical :: ok

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "cc-pvdz", orb, err)
      call check(error,.not. err%has_error(), "water in cc-pVDZ should build")
      if (allocated(error)) return

      call build_aambs(aambs, err, ok)
      call check(error, ok, "water in the atomic minimal basis should build")
      if (allocated(error)) return

      call aambs_dimensions(WATER_Z, 10, dims, err)
      call check(error,.not. err%has_error(), "the counting should succeed")
      if (allocated(error)) return

      ! Water: oxygen contributes 1 core + 4 valence, each hydrogen 1. Seven.
      call check(error, dims%n_mbs, 7, "water has a seven-orbital minimal basis")
      if (allocated(error)) return
      call check(error, aambs%nao, dims%n_mbs, &
                 "the built minimal basis must have as many functions as the "// &
                 "per-element counts say")
      if (allocated(error)) return

      call mixed_basis_overlap(orb, aambs, mixed, err)
      call check(error,.not. err%has_error(), "the mixed overlap should build")
      if (allocated(error)) return
      call check(error, size(mixed, 1), orb%nao, "rows are the orbital basis")
      if (allocated(error)) return
      call check(error, size(mixed, 2), dims%n_mbs, "columns are the minimal basis")
      if (allocated(error)) return

      ! A projection onto a minimal basis that lands nowhere would leave this
      ! zero, and every downstream singular value with it.
      call check(error, maxval(abs(mixed)) > 0.1_dp, &
                 "the two bases should actually overlap")
      if (allocated(error)) return

      call orb%destroy()
      call aambs%destroy()
   end subroutine test_shape

   subroutine test_refusal(error)
      !! Cartesian against spherical is refused rather than silently integrated
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: spherical, cartesian
      type(error_t) :: err
      real(dp), allocatable :: mixed(:, :)

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "cc-pvdz", spherical, err)
      call check(error,.not. err%has_error(), "the spherical build should succeed")
      if (allocated(error)) return
      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "cc-pvdz", cartesian, err, &
                              force_cartesian=.true.)
      call check(error,.not. err%has_error(), "the Cartesian build should succeed")
      if (allocated(error)) return

      call mixed_basis_overlap(spherical, cartesian, mixed, err)
      call check(error, err%has_error(), &
                 "one basis spherical and the other Cartesian should be refused: "// &
                 "libcint takes a single convention for the whole call")

      call spherical%destroy()
      call cartesian%destroy()
   end subroutine test_refusal

end module test_mqc_mixed_overlap

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mixed_overlap, only: collect_mqc_mixed_overlap_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_mixed_overlap", collect_mqc_mixed_overlap_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
