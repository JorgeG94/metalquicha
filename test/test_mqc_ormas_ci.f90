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
   use mqc_ormas_space, only: ormas_space_t, build_ormas_space, ormas_strings, &
                              determinant_address, ormas_closure_t, build_ormas_closure
   use mqc_ci, only: ci_diagonal, absorb_one_electron
   use mqc_rdm, only: rdm_energy, active_space_rdms
   use mqc_ormas_ci, only: ormas_diagonal, full_space_index, ormas_lowest, &
                           ormas_sigma, ormas_sigma_direct, ormas_solve, &
                           ormas_density_matrices
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
                  new_unittest("energies_match_the_python_reference", test_energies), &
                  new_unittest("direct_sigma_matches_the_projected_one", test_direct), &
                  new_unittest("iterative_solver_finds_the_same_states", test_iterative), &
                  new_unittest("densities_rebuild_the_energy", test_densities), &
                  new_unittest("rejects_mismatched_integrals", test_refusals) &
                  ]
   end subroutine collect_mqc_ormas_ci_tests

   subroutine model_integrals(norb, h1e, eri, dominant)
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
      logical, intent(in), optional :: dominant
         !! Separate the orbitals in energy, so the determinants separate too.
         !! The plain model puts the whole spectrum inside a few times 1e-4,
         !! which is fine for checking a Hamiltonian and useless for checking a
         !! solver that preconditions with the diagonal -- there the diagonal
         !! says almost nothing about where a determinant's energy sits, and
         !! Davidson stalls for reasons that have nothing to do with the code
         !! under test.

      real(dp) :: b(norb, norb, NCHOL)
      integer :: p, q, r, s, l

      allocate (h1e(norb, norb), eri(norb, norb, norb, norb))
      do q = 1, norb
         do p = 1, norb
            h1e(p, q) = -1.0_dp/real((p - 1) + (q - 1) + 2, dp)
         end do
      end do
      if (present(dominant)) then
         if (dominant) then
            do p = 1, norb
               h1e(p, p) = -2.0_dp*real(norb - p + 1, dp)
            end do
         end if
      end if
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

   subroutine energies_of(first_orbital, norb, na, nb, min_e, max_e, expected, error)
      !! The three lowest roots of one partition, against transcribed values
      integer, intent(in) :: first_orbital(:), min_e(:), max_e(:)
      integer, intent(in) :: norb, na, nb
      real(dp), intent(in) :: expected(:)
      type(error_type), allocatable, intent(inout) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err
      type(link_table_t) :: alpha_table, beta_table
      integer(int64), allocatable :: alpha(:), beta(:)
      integer, allocatable :: in_alpha(:), in_beta(:)
      real(dp), allocatable :: h1e(:, :), eri(:, :, :, :), folded(:, :), energies(:)
      integer :: root

      call model_integrals(norb, h1e, eri)
      call absorb_one_electron(h1e, eri, na + nb, folded, err)
      call build_link_table(norb, na, alpha_table, err)
      call build_link_table(norb, nb, beta_table, err)

      call build_ormas_space(first_orbital, norb, na, nb, min_e, max_e, space, err)
      call ormas_strings(space, alpha, beta, err)
      call full_space_index(space, alpha, beta, in_alpha, in_beta)
      call ormas_lowest(space, folded, alpha_table, beta_table, in_alpha, in_beta, &
                        size(expected), energies, err)
      call check(error,.not. err%has_error(), "solving the restricted space failed")
      if (allocated(error)) return

      do root = 1, size(expected)
         call check(error, energies(root), expected(root), thr=1.0e-10_dp, &
                    message="a root disagrees with the python reference")
         if (allocated(error)) return
      end do

      call space%destroy()
      call alpha_table%destroy()
      call beta_table%destroy()
   end subroutine energies_of

   subroutine test_energies(error)
      !! Every partition, against energies computed outside this language
      !!
      !! The numbers come from `tools/ormas_reference/ormas_reference.py`, which
      !! enumerates the determinants its own way, builds the Hamiltonian by the
      !! Slater-Condon rules and diagonalises it, and which reproduces
      !! `pyscf.fci` exactly when the windows exclude nothing. Nothing about the
      !! occupation classes, the compatibility grid or the addressing exists on
      !! that side, so agreement here is a statement about the space and the
      !! physics rather than about a shared convention.
      !!
      !! Three roots rather than one on purpose. A ground state can come out
      !! right from a Hamiltonian that is wrong on everything above it -- the
      !! lowest eigenvalue is the least sensitive thing the matrix has to say.
      type(error_type), allocatable, intent(out) :: error

      call energies_of([1, 3], 4, 2, 2, [2, 0], [4, 2], &
                       [-1.058450665249_dp, -1.049170312395_dp, -1.046492881715_dp], error)
      if (allocated(error)) return

      call energies_of([1, 3, 5], 6, 3, 3, [2, 0, 0], [4, 4, 2], &
                       [-1.186905707931_dp, -1.186602623927_dp, -1.186573312464_dp], error)
      if (allocated(error)) return

      call energies_of([1, 4], 7, 3, 3, [4, 0], [6, 2], &
                       [-1.229988119132_dp, -1.214487093893_dp, -1.212223692704_dp], error)
      if (allocated(error)) return

      call energies_of([1, 3, 5], 6, 3, 1, [1, 0, 0], [4, 4, 2], &
                       [-1.185387147112_dp, -1.184328705476_dp, -1.175651771473_dp], error)
      if (allocated(error)) return

      call energies_of([1, 4], 6, 2, 2, [0, 3], [6, 6], &
                       [-1.121243956258_dp, -1.121240844637_dp, -1.121137600250_dp], error)
      if (allocated(error)) return
   end subroutine test_energies

   subroutine direct_matches(first_orbital, norb, na, nb, min_e, max_e, error)
      !! One partition: the two sigma builds, on several vectors
      integer, intent(in) :: first_orbital(:), min_e(:), max_e(:)
      integer, intent(in) :: norb, na, nb
      type(error_type), allocatable, intent(inout) :: error

      type(ormas_space_t) :: space
      type(ormas_closure_t) :: closure
      type(error_t) :: err
      type(link_table_t) :: alpha_table, beta_table
      integer(int64), allocatable :: alpha(:), beta(:)
      integer, allocatable :: in_alpha(:), in_beta(:), from_alpha(:), from_beta(:)
      real(dp), allocatable :: h1e(:, :), eri(:, :, :, :), folded(:, :)
      real(dp), allocatable :: ci(:), projected(:), direct(:)
      integer :: ndet, d, trial

      call model_integrals(norb, h1e, eri)
      call absorb_one_electron(h1e, eri, na + nb, folded, err)
      call build_link_table(norb, na, alpha_table, err)
      call build_link_table(norb, nb, beta_table, err)

      call build_ormas_space(first_orbital, norb, na, nb, min_e, max_e, space, err)
      call ormas_strings(space, alpha, beta, err)
      call full_space_index(space, alpha, beta, in_alpha, in_beta, from_alpha, from_beta)
      call build_ormas_closure(space, closure, err)
      call check(error,.not. err%has_error(), "building the closure failed")
      if (allocated(error)) return

      ! Wider than the space, because an excitation leaves it; narrower than the
      ! rectangle, because that is the saving the whole exercise is for.
      call check(error, closure%n_determinants >= space%n_determinants, &
                 "the closure does not contain the space")
      if (allocated(error)) return
      call check(error, closure%n_determinants <= &
                 int(size(alpha), int64)*int(size(beta), int64), &
                 "the closure exceeds the complete space")
      if (allocated(error)) return

      ndet = int(space%n_determinants)
      allocate (ci(ndet), projected(ndet), direct(ndet))

      do trial = 1, 3
         ! Deterministic and spread over several orders of magnitude, so a term
         ! dropped anywhere shows rather than hiding under a larger one.
         do d = 1, ndet
            ci(d) = sin(real(d*trial, dp))/real(d + trial, dp)
         end do

         call ormas_sigma(space, folded, alpha_table, beta_table, in_alpha, in_beta, &
                          ci, projected, err)
         call ormas_sigma_direct(space, closure, folded, alpha_table, beta_table, &
                                 in_alpha, in_beta, from_alpha, from_beta, ci, direct, err)
         call check(error,.not. err%has_error(), "a sigma build failed")
         if (allocated(error)) return

         do d = 1, ndet
            call check(error, direct(d), projected(d), thr=1.0e-11_dp, &
                       message="the direct sigma differs from the projected one")
            if (allocated(error)) return
         end do
      end do

      ! And symmetric, which is a different claim from matching the projected
      ! build and fails differently: an asymmetric sigma gives a Davidson that
      ! stalls at a residual it cannot reduce rather than one that reports
      ! anything wrong.
      block
         real(dp), allocatable :: matrix(:, :)
         real(dp) :: worst
         integer :: r, c
         allocate (matrix(ndet, ndet))
         do d = 1, ndet
            ci = 0.0_dp
            ci(d) = 1.0_dp
            call ormas_sigma_direct(space, closure, folded, alpha_table, beta_table, &
                                    in_alpha, in_beta, from_alpha, from_beta, ci, direct, err)
            matrix(:, d) = direct
         end do
         worst = 0.0_dp
         do c = 1, ndet
            do r = 1, ndet
               worst = max(worst, abs(matrix(r, c) - matrix(c, r)))
            end do
         end do
         call check(error, worst < 1.0e-11_dp, "the direct sigma is not symmetric")
         deallocate (matrix)
      end block
      if (allocated(error)) return

      deallocate (ci, projected, direct)
      call space%destroy()
      call closure%destroy()
      call alpha_table%destroy()
      call beta_table%destroy()
   end subroutine direct_matches

   subroutine test_direct(error)
      !! The sigma that stays inside the closure, against the one that does not
      !!
      !! The projected build is already known good -- its energies reproduce an
      !! independent reference that itself reproduces PySCF. So this asks a
      !! single sharp question: does restricting the intermediate to the space
      !! plus one excitation lose anything? If the closure were too small the
      !! two would differ, and differ in exactly the terms that are hardest to
      !! reason about -- a double excitation whose intermediate step is outside
      !! the space, which no count and no symmetry would miss.
      !!
      !! Compared vector by vector rather than through an energy. An eigenvalue
      !! is a scalar summary and a wrong sigma can still produce a plausible
      !! one; two sigma vectors agreeing at every determinant, on several
      !! unrelated inputs, cannot be a coincidence.
      type(error_type), allocatable, intent(out) :: error

      call direct_matches([1, 3], 4, 2, 2, [2, 0], [4, 2], error)
      if (allocated(error)) return

      call direct_matches([1, 3, 5], 6, 3, 3, [2, 0, 0], [4, 4, 2], error)
      if (allocated(error)) return

      call direct_matches([1, 4], 7, 3, 3, [4, 0], [6, 2], error)
      if (allocated(error)) return

      call direct_matches([1, 3, 5], 6, 3, 1, [1, 0, 0], [4, 4, 2], error)
      if (allocated(error)) return

      call direct_matches([1, 4], 6, 2, 2, [0, 3], [6, 6], error)
      if (allocated(error)) return

      ! And a partition that restricts nothing, where the closure is the whole
      ! rectangle and the two builds must agree trivially.
      call direct_matches([1, 3], 4, 2, 2, [0, 0], [4, 4], error)
      if (allocated(error)) return
   end subroutine test_direct

   subroutine iterative_of(first_orbital, norb, na, nb, min_e, max_e, expected, error)
      !! Davidson over one partition, against the transcribed reference
      integer, intent(in) :: first_orbital(:), min_e(:), max_e(:)
      integer, intent(in) :: norb, na, nb
      real(dp), intent(in) :: expected(:)
      type(error_type), allocatable, intent(inout) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err
      real(dp), allocatable :: h1e(:, :), eri(:, :, :, :)
      real(dp), allocatable :: energies(:), vectors(:, :)
      integer :: root
      real(dp) :: norm

      call model_integrals(norb, h1e, eri, dominant=.true.)
      call build_ormas_space(first_orbital, norb, na, nb, min_e, max_e, space, err)
      call ormas_solve(space, h1e, eri, size(expected), energies, vectors, err)
      call check(error,.not. err%has_error(), "the iterative solve failed")
      if (allocated(error)) return

      do root = 1, size(expected)
         call check(error, energies(root), expected(root), thr=1.0e-9_dp, &
                    message="an iterative root disagrees with the reference")
         if (allocated(error)) return

         norm = sqrt(sum(vectors(:, root)**2))
         call check(error, norm, 1.0_dp, thr=1.0e-10_dp, &
                    message="an eigenvector came back unnormalised")
         if (allocated(error)) return
      end do

      call space%destroy()
   end subroutine iterative_of

   subroutine test_iterative(error)
      !! Davidson finds what the dense solve found
      !!
      !! Same reference values as `energies_match_the_python_reference`, reached
      !! without ever forming a matrix. Comparing against the same numbers on
      !! purpose: the dense route and this one share the space, the addressing
      !! and the integrals but not the sigma -- one goes through the complete
      !! space and the other through the closure -- so agreement is a statement
      !! about the two sigma builds, and about the solver finding the states
      !! rather than merely something low.
      !!
      !! Three roots each, because a preconditioner that is subtly wrong still
      !! finds a ground state and stalls or lands elsewhere above it.
      !!
      !! Against a separated one-electron diagonal, not the one the other tests
      !! use. With the plain model every state of these spaces lies within a few
      !! times 1e-4 of every other -- two of them nine digits apart -- and a
      !! diagonal preconditioner has nothing to tell the determinants apart
      !! with. Davidson then stalls at a residual it cannot reduce while the
      !! energies are already right to eight digits, which is a statement about
      !! the model and not about the solver. Separating the orbitals in energy
      !! puts the problem in the regime any real CI is in, where the gaps here
      !! run from 0.08 to 1.9.
      type(error_type), allocatable, intent(out) :: error

      call iterative_of([1, 3], 4, 2, 2, [2, 0], [4, 2], &
                        [-27.420851433945_dp, -25.535119017717_dp, -25.422516051785_dp], error)
      if (allocated(error)) return

      call iterative_of([1, 3, 5], 6, 3, 3, [2, 0, 0], [4, 4, 2], &
                        [-59.024917557510_dp, -57.103222181291_dp, -57.026760563957_dp], error)
      if (allocated(error)) return

      call iterative_of([1, 4], 7, 3, 3, [4, 0], [6, 2], &
                        [-71.029902648387_dp, -69.101862579285_dp, -69.026073765100_dp], error)
      if (allocated(error)) return

      call iterative_of([1, 3, 5], 6, 3, 1, [1, 0, 0], [4, 4, 2], &
                        [-41.547254882025_dp, -39.740658079124_dp, -39.547023100484_dp], error)
      if (allocated(error)) return

      call iterative_of([1, 4], 6, 2, 2, [0, 3], [6, 6], &
                        [-27.880210126288_dp, -27.766902778851_dp, -25.947416001524_dp], error)
      if (allocated(error)) return
   end subroutine test_iterative

   subroutine densities_of(first_orbital, norb, na, nb, min_e, max_e, error)
      !! One partition: solve, take densities, rebuild the energy from them
      integer, intent(in) :: first_orbital(:), min_e(:), max_e(:)
      integer, intent(in) :: norb, na, nb
      type(error_type), allocatable, intent(inout) :: error

      type(ormas_space_t) :: space
      type(error_t) :: err
      real(dp), allocatable :: h1e(:, :), eri(:, :, :, :)
      real(dp), allocatable :: energies(:), vectors(:, :)
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :)
      real(dp) :: electrons
      integer :: p

      call model_integrals(norb, h1e, eri, dominant=.true.)
      call build_ormas_space(first_orbital, norb, na, nb, min_e, max_e, space, err)
      call ormas_solve(space, h1e, eri, 1, energies, vectors, err)
      call ormas_density_matrices(space, vectors(:, 1), dm1, dm2, err)
      call check(error,.not. err%has_error(), "building the densities failed")
      if (allocated(error)) return

      ! The electron count, which a wrong intermediate would usually still get
      ! right -- worth asserting, and worth not stopping at.
      electrons = 0.0_dp
      do p = 1, norb
         electrons = electrons + dm1(p, p)
      end do
      call check(error, electrons, real(na + nb, dp), thr=1.0e-10_dp, &
                 message="the one-particle density does not hold the electrons")
      if (allocated(error)) return

      ! The energy, which it would not. Reached by arithmetic that shares
      ! nothing with the eigensolver: a contraction of two density matrices
      ! against the integrals, against a Davidson eigenvalue.
      call check(error, rdm_energy(h1e, eri, dm1, dm2), energies(1), thr=1.0e-9_dp, &
                 message="the densities do not rebuild the CI energy")
      if (allocated(error)) return

      call space%destroy()
   end subroutine densities_of

   subroutine test_densities(error)
      !! The density matrices, by the one check that catches everything
      !!
      !! Contracting them back against the integrals has to give the energy the
      !! CI already found. The two numbers are reached by completely different
      !! arithmetic, so a transposed index, a factor of two, or an intermediate
      !! summed over the space instead of its closure will not survive it --
      !! while every trace identity would still hold in all three cases.
      !!
      !! That last one is the reason this test is here rather than left to a
      !! trace. The two-particle matrix is an inner product of two
      !! intermediates, both of which have weight outside the space; restricting
      !! the sum to the space loses those products quietly and leaves densities
      !! that look entirely healthy.
      type(error_type), allocatable, intent(out) :: error

      call densities_of([1, 3], 4, 2, 2, [2, 0], [4, 2], error)
      if (allocated(error)) return

      call densities_of([1, 3, 5], 6, 3, 3, [2, 0, 0], [4, 4, 2], error)
      if (allocated(error)) return

      call densities_of([1, 4], 7, 3, 3, [4, 0], [6, 2], error)
      if (allocated(error)) return

      call densities_of([1, 3, 5], 6, 3, 1, [1, 0, 0], [4, 4, 2], error)
      if (allocated(error)) return

      call densities_of([1, 4], 6, 2, 2, [0, 3], [6, 6], error)
      if (allocated(error)) return

      ! And against the complete-space densities where the two must agree.
      call densities_against_cas(error)
   end subroutine test_densities

   subroutine densities_against_cas(error)
      !! With windows that restrict nothing, `mqc_rdm`'s answer, element by
      !! element
      type(error_type), allocatable, intent(inout) :: error

      integer, parameter :: NORB = 4, NA = 2, NB = 2
      type(ormas_space_t) :: space
      type(error_t) :: err
      type(link_table_t) :: alpha_table, beta_table
      real(dp), allocatable :: h1e(:, :), eri(:, :, :, :)
      real(dp), allocatable :: energies(:), vectors(:, :)
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :)
      real(dp), allocatable :: ref1(:, :), ref2(:, :, :, :), square(:, :)
      integer :: p, q, r, s

      call model_integrals(NORB, h1e, eri, dominant=.true.)
      call build_ormas_space([1], NORB, NA, NB, [NA + NB], [NA + NB], space, err)
      call ormas_solve(space, h1e, eri, 1, energies, vectors, err)
      call ormas_density_matrices(space, vectors(:, 1), dm1, dm2, err)
      call check(error,.not. err%has_error(), "the restricted densities failed")
      if (allocated(error)) return

      call build_link_table(NORB, NA, alpha_table, err)
      call build_link_table(NORB, NB, beta_table, err)
      allocate (square(alpha_table%n_strings, beta_table%n_strings))
      square = reshape(vectors(:, 1), [alpha_table%n_strings, beta_table%n_strings])
      call active_space_rdms(square, alpha_table, beta_table, ref1, ref2, err)
      call check(error,.not. err%has_error(), "the complete-space densities failed")
      if (allocated(error)) return

      do q = 1, NORB
         do p = 1, NORB
            call check(error, dm1(p, q), ref1(p, q), thr=1.0e-10_dp, &
                       message="the one-particle densities differ")
            if (allocated(error)) return
         end do
      end do
      do s = 1, NORB
         do r = 1, NORB
            do q = 1, NORB
               do p = 1, NORB
                  call check(error, dm2(p, q, r, s), ref2(p, q, r, s), thr=1.0e-10_dp, &
                             message="the two-particle densities differ")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do

      deallocate (square)
      call space%destroy()
      call alpha_table%destroy()
      call beta_table%destroy()
   end subroutine densities_against_cas

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
