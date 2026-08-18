!! The Hamiltonian diagonal over an occupation-restricted active space
module test_mqc_ormas_ci
   !! Two things are being checked here at once, and they are worth separating.
   !!
   !! The *physics* is not new. A restricted space contains fewer determinants
   !! than a complete one, not different ones, so every diagonal element is the
   !! same Slater rule it always was. `matches_slater_rules` recomputes each
   !! one from the two occupation strings by the textbook expression and
   !! compares.
   !!
   !! The *layout* is new, and is the reason this test exists at all.
   !! `ormas_diagonal` fills its array by walking the space in storage order
   !! with a running counter, never once evaluating an address; the test reaches
   !! the same elements through `determinant_address`. Agreement at every
   !! determinant means the counter and the formula describe the same layout,
   !! which is a far sharper statement than any energy could make and is
   !! available before a single excitation has been written. If the two ever
   !! disagree, the offsets are wrong and everything built on them would be
   !! wrong in a way that shows up as a plausible number much later.
   !!
   !! `agrees_with_the_cas_diagonal` closes the loop from the other side: a
   !! partition that restricts nothing must reproduce, element for element, the
   !! diagonal that `mqc_ci` was validated against PySCF for.
   !!
   !! The traces come from `tools/ormas_reference/ormas_reference.py`, which
   !! builds the same model Hamiltonian and the same partitions in Python and
   !! is itself checked against `pyscf.fci`. It shares no code, no ordering and
   !! no language with any of this, so agreement is not two readings of one
   !! mistake -- and unlike the element-by-element comparison, a trace notices
   !! a determinant that was visited twice.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp, int64
   use mqc_error, only: error_t
   use mqc_determinants, only: link_table_t, build_link_table
   use mqc_ci, only: ci_diagonal
   use mqc_ormas_space, only: ormas_space_t, build_ormas_space, ormas_strings, &
                              determinant_address
   use mqc_ormas_ci, only: ormas_diagonal
   implicit none
   private

   public :: collect_mqc_ormas_ci_tests

   integer, parameter :: NCHOL = 3
   real(dp), parameter :: TOL = 1.0e-12_dp

contains

   subroutine collect_mqc_ormas_ci_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("matches_slater_rules", test_slater), &
                  new_unittest("agrees_with_the_cas_diagonal", test_against_cas), &
                  new_unittest("rejects_mismatched_integrals", test_refusals) &
                  ]
   end subroutine collect_mqc_ormas_ci_tests

   subroutine model_integrals(norb, h1e, eri)
      !! A Hamiltonian with the right symmetries and no physics
      !!
      !! `h_pq = -1/(p+q)` and a Cholesky-factorised two-electron tensor
      !! `(pq|rs) = sum_L B_pqL B_rsL` with `B_pqL = 1/(p+q+L)`, on zero-based
      !! indices. Factorising rather than filling entries at random is what
      !! makes `(pq|rs)` carry the eightfold permutational symmetry a real
      !! integral has, which the diagonal expression quietly assumes.
      integer, intent(in) :: norb
      real(dp), allocatable, intent(out) :: h1e(:, :)
      real(dp), allocatable, intent(out) :: eri(:, :, :, :)

      real(dp) :: b(norb, norb, NCHOL)
      integer :: p, q, r, s, l

      allocate (h1e(norb, norb), eri(norb, norb, norb, norb))
      do q = 1, norb
         do p = 1, norb
            h1e(p, q) = -1.0_dp/real((p - 1) + (q - 1) + 2, dp)
         end do
      end do
      do l = 1, NCHOL
         do q = 1, norb
            do p = 1, norb
               b(p, q, l) = 1.0_dp/real((p - 1) + (q - 1) + (l - 1) + 3, dp)
            end do
         end do
      end do

      eri = 0.0_dp
      do s = 1, norb
         do r = 1, norb
            do q = 1, norb
               do p = 1, norb
                  do l = 1, NCHOL
                     eri(p, q, r, s) = eri(p, q, r, s) + b(p, q, l)*b(r, s, l)
                  end do
               end do
            end do
         end do
      end do
   end subroutine model_integrals

   pure function slater_diagonal(alpha, beta, norb, h1e, eri) result(element)
      !! `<D|H|D>` written out directly from the two strings
      !!
      !! Deliberately the naive form -- every pair, no hoisting, no
      !! precomputed vectors -- so that it shares nothing with the
      !! implementation but the integrals themselves.
      integer(int64), intent(in) :: alpha, beta
      integer, intent(in) :: norb
      real(dp), intent(in) :: h1e(:, :), eri(:, :, :, :)
      real(dp) :: element

      integer :: p, q

      element = 0.0_dp

      do p = 1, norb
         if (btest(alpha, p - 1)) element = element + h1e(p, p)
         if (btest(beta, p - 1)) element = element + h1e(p, p)
      end do

      do p = 1, norb
         do q = p + 1, norb
            if (btest(alpha, p - 1) .and. btest(alpha, q - 1)) &
               element = element + eri(p, p, q, q) - eri(p, q, q, p)
            if (btest(beta, p - 1) .and. btest(beta, q - 1)) &
               element = element + eri(p, p, q, q) - eri(p, q, q, p)
         end do
      end do

      do p = 1, norb
         if (.not. btest(alpha, p - 1)) cycle
         do q = 1, norb
            if (btest(beta, q - 1)) element = element + eri(p, p, q, q)
         end do
      end do
   end function slater_diagonal

   subroutine check_one_partition(first_orbital, norb, na, nb, min_e, max_e, trace, error)
      !! Every determinant of one partition, by Slater rules and by address
      !!
      !! `trace` is the sum of the diagonal, computed independently by
      !! `tools/ormas_reference/ormas_reference.py`. Being a sum it does not
      !! care what order the determinants were enumerated in, which is what
      !! lets a number from another language check this one at all -- and it
      !! catches a determinant counted twice or missed, which comparing
      !! elements one at a time cannot.
      integer, intent(in) :: first_orbital(:), min_e(:), max_e(:)
      integer, intent(in) :: norb, na, nb
      real(dp), intent(in) :: trace
      type(error_type), allocatable, intent(inout) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err
      integer(int64), allocatable :: alpha(:), beta(:)
      real(dp), allocatable :: h1e(:, :), eri(:, :, :, :), diagonal(:)
      logical, allocatable :: seen(:)
      integer :: ia, ib, ga, gb
      integer(int64) :: at
      real(dp) :: expected

      call model_integrals(norb, h1e, eri)
      call build_ormas_space(first_orbital, norb, na, nb, min_e, max_e, space, err)
      call ormas_strings(space, alpha, beta, err)
      call ormas_diagonal(space, alpha, beta, h1e, eri, diagonal, err)
      call check(error,.not. err%has_error(), "building the diagonal failed")
      if (allocated(error)) return

      call check(error, size(diagonal, kind=int64) == space%n_determinants, &
                 "the diagonal is the wrong length")
      if (allocated(error)) return

      allocate (seen(space%n_determinants))
      seen = .false.

      do ia = 1, size(alpha)
         do ib = 1, size(beta)
            ga = space%alpha_string_class(ia)
            gb = space%beta_string_class(ib)
            if (.not. space%compatible(gb, ga)) cycle

            at = determinant_address(space, ia, ib)
            expected = slater_diagonal(alpha(ia), beta(ib), norb, h1e, eri)
            call check(error, diagonal(at), expected, thr=TOL, &
                       message="a diagonal element disagrees with the Slater rules")
            if (allocated(error)) return
            seen(at) = .true.
         end do
      end do

      call check(error, all(seen), "the walk skipped a determinant the address "// &
                 "formula reaches")
      if (allocated(error)) return

      call check(error, sum(diagonal), trace, thr=1.0e-9_dp, &
                 message="the trace disagrees with the python reference")
      if (allocated(error)) return

      deallocate (seen, diagonal, h1e, eri, alpha, beta)
      call space%destroy()
   end subroutine check_one_partition

   subroutine test_slater(error)
      !! Four partitions, every determinant of each
      type(error_type), allocatable, intent(out) :: error

      call check_one_partition([1, 3], 4, 2, 2, [2, 0], [4, 2], &
                               -20.010664176292_dp, error)
      if (allocated(error)) return

      call check_one_partition([1, 3, 5], 6, 3, 3, [2, 0, 0], [4, 4, 2], &
                               -173.027624139105_dp, error)
      if (allocated(error)) return

      ! Subspaces of unequal width.
      call check_one_partition([1, 4], 7, 3, 3, [4, 0], [6, 2], &
                               -155.090056159540_dp, error)
      if (allocated(error)) return

      ! Unequal spins, so the alpha and beta class lists differ and the layout
      ! is not symmetric.
      call check_one_partition([1, 3, 5], 6, 3, 1, [1, 0, 0], [4, 4, 2], &
                               -67.145428591886_dp, error)
      if (allocated(error)) return

      ! Windows the tightening pass rewrites.
      call check_one_partition([1, 4], 6, 2, 2, [0, 3], [6, 6], &
                               -28.034091825652_dp, error)
      if (allocated(error)) return
   end subroutine test_slater

   subroutine test_against_cas(error)
      !! A partition that restricts nothing reproduces `mqc_ci`'s diagonal
      !!
      !! Both spellings of an unrestricted space are checked: one subspace, and
      !! two subspaces with open windows -- the second has three occupation
      !! classes per spin and a full compatibility grid, so it exercises the
      !! whole layout and must still come out as the plain rectangle.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: NORB = 4, NA = 2, NB = 2
      type(ormas_space_t) :: space
      type(error_t) :: err
      type(link_table_t) :: alpha_table, beta_table
      integer(int64), allocatable :: alpha(:), beta(:)
      real(dp), allocatable :: h1e(:, :), eri(:, :, :, :)
      real(dp), allocatable :: diagonal(:), reference(:, :)
      integer :: ia, ib
      integer(int64) :: at

      call model_integrals(NORB, h1e, eri)
      call build_link_table(NORB, NA, alpha_table, err)
      call build_link_table(NORB, NB, beta_table, err)
      call ci_diagonal(h1e, eri, alpha_table, beta_table, reference, err)
      call check(error,.not. err%has_error(), "the reference diagonal failed")
      if (allocated(error)) return

      call build_ormas_space([1], NORB, NA, NB, [NA + NB], [NA + NB], space, err)
      call ormas_strings(space, alpha, beta, err)
      call ormas_diagonal(space, alpha, beta, h1e, eri, diagonal, err)
      call check(error,.not. err%has_error(), "the restricted diagonal failed")
      if (allocated(error)) return

      do ib = 1, size(beta)
         do ia = 1, size(alpha)
            at = determinant_address(space, ia, ib)
            call check(error, diagonal(at), reference(ia, ib), thr=TOL, &
                       message="one open subspace differs from the CAS diagonal")
            if (allocated(error)) return
         end do
      end do
      deallocate (diagonal, alpha, beta)
      call space%destroy()

      call build_ormas_space([1, 3], NORB, NA, NB, [0, 0], [4, 4], space, err)
      call ormas_strings(space, alpha, beta, err)
      call ormas_diagonal(space, alpha, beta, h1e, eri, diagonal, err)
      call check(error,.not. err%has_error(), "the restricted diagonal failed")
      if (allocated(error)) return
      call check(error, space%n_alpha_classes, 3, "the partition should have classes")
      if (allocated(error)) return

      ! The strings are grouped by class here and so are not in the order
      ! `mqc_ci` uses, which is exactly why the comparison goes through the
      ! address rather than through a shared index.
      do ib = 1, size(beta)
         do ia = 1, size(alpha)
            at = determinant_address(space, ia, ib)
            call check(error, diagonal(at), &
                       slater_diagonal(alpha(ia), beta(ib), NORB, h1e, eri), thr=TOL, &
                       message="open windows differ from the Slater rules")
            if (allocated(error)) return
         end do
      end do

      call check(error, abs(sum(diagonal) - sum(reference)) < TOL, &
                 "the two spellings of an unrestricted space have different traces")
      if (allocated(error)) return

      call space%destroy()
      call alpha_table%destroy()
      call beta_table%destroy()
   end subroutine test_against_cas

   subroutine test_refusals(error)
      !! Integrals that do not span the partition are refused
      type(error_type), allocatable, intent(out) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err
      integer(int64), allocatable :: alpha(:), beta(:)
      real(dp), allocatable :: h1e(:, :), eri(:, :, :, :), diagonal(:)

      call model_integrals(5, h1e, eri)
      call build_ormas_space([1, 3], 4, 2, 2, [2, 0], [4, 2], space, err)
      call ormas_strings(space, alpha, beta, err)
      call ormas_diagonal(space, alpha, beta, h1e, eri, diagonal, err)
      call check(error, err%has_error(), "integrals of the wrong size were accepted")
      if (allocated(error)) return

      call space%destroy()
   end subroutine test_refusals

end module test_mqc_ormas_ci

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_ormas_ci, only: collect_mqc_ormas_ci_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_ormas_ci", collect_mqc_ormas_ci_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
