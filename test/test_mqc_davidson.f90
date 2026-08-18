!! The iterative eigensolver, against a dense diagonalisation
module test_mqc_davidson
   !! Every reference here is the dense solve of the same Hamiltonian, not
   !! PySCF. That separation is the point of having built `ci_hamiltonian`: the
   !! Hamiltonian itself is already validated against PySCF in
   !! `test_mqc_ci.f90`, so what is left to establish is that the iterative
   !! solver finds the eigenvalues that are actually there. Those are two
   !! different failures and they deserve two different tests -- a solver that
   !! converges neatly onto the wrong number would pass a comparison against
   !! PySCF only if the Hamiltonian were also wrong in the same way.
   !!
   !! **The model here is not the model in `test_mqc_ci.f90`, and the difference
   !! matters.** That one is built from reciprocals and comes out with its
   !! eigenvalues below every diagonal element and 3e-6 apart from each other.
   !! For testing a sigma build that is fine -- the contraction does not care
   !! whether the numbers are physical. For testing a *preconditioned* solver it
   !! is close to worthless: the preconditioner divides by `theta - H_DD`, which
   !! there is nearly the same number for every determinant, so Davidson
   !! degenerates into steepest descent and stalls at 1e-4 on a near-degenerate
   !! ground state. Measured, before this model replaced it.
   !!
   !! A real CI Hamiltonian is strongly diagonally dominant -- a determinant far
   !! from the reference is far from it in energy -- which is the property the
   !! method is built on. So the model below puts orbital energies on the
   !! diagonal of `h1e` and weak coupling off it, and the two-electron part is
   !! scaled to be a perturbation rather than the whole Hamiltonian. Testing a
   !! solver against a system its central approximation does not describe
   !! measures nothing.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t
   use mqc_determinants, only: link_table_t, build_link_table
   use mqc_ci, only: absorb_one_electron, ci_hamiltonian, ci_diagonal
   use mqc_davidson, only: davidson_lowest, davidson_result_t
   implicit none
   private

   public :: collect_mqc_davidson_tests

   integer, parameter :: NORB = 6
   integer, parameter :: NCHOL = 3

contains

   subroutine collect_mqc_davidson_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ground_state_matches_dense", test_ground_state), &
                  new_unittest("several_roots_match_dense", test_several_roots), &
                  new_unittest("eigenvectors_are_eigenvectors", test_vectors), &
                  new_unittest("subspace_collapse", test_collapse), &
                  new_unittest("guess_is_used", test_guess), &
                  new_unittest("refusals", test_refusals) &
                  ]
   end subroutine collect_mqc_davidson_tests

   subroutine model(n_alpha, n_beta, folded, diagonal, alpha, beta, dense, err, ok)
      !! A model Hamiltonian, in both forms
      integer, intent(in) :: n_alpha, n_beta
      real(dp), allocatable, intent(out) :: folded(:, :), diagonal(:, :), dense(:, :)
      type(link_table_t), intent(out) :: alpha, beta
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      real(dp) :: h1e(NORB, NORB), eri(NORB, NORB, NORB, NORB), b(NORB, NORB, NCHOL)
      integer :: p, q, r, s, l

      ok = .false.
      do q = 1, NORB
         do p = 1, NORB
            if (p == q) then
               h1e(p, q) = -real(NORB - p + 1, dp)
            else
               h1e(p, q) = -0.1_dp/real(abs(p - q) + 1, dp)
            end if
         end do
      end do
      do l = 1, NCHOL
         do q = 1, NORB
            do p = 1, NORB
               b(p, q, l) = 0.45_dp/real((p - 1) + (q - 1) + (l - 1) + 3, dp)
            end do
         end do
      end do
      eri = 0.0_dp
      do s = 1, NORB
         do r = 1, NORB
            do q = 1, NORB
               do p = 1, NORB
                  do l = 1, NCHOL
                     eri(p, q, r, s) = eri(p, q, r, s) + b(p, q, l)*b(r, s, l)
                  end do
               end do
            end do
         end do
      end do

      call absorb_one_electron(h1e, eri, n_alpha + n_beta, folded, err)
      call build_link_table(NORB, n_alpha, alpha, err)
      call build_link_table(NORB, n_beta, beta, err)
      if (err%has_error()) return
      call ci_diagonal(h1e, eri, alpha, beta, diagonal, err)
      call ci_hamiltonian(folded, alpha, beta, dense, err)
      ok = .not. err%has_error()
   end subroutine model

   subroutine dense_spectrum(dense, values)
      real(dp), intent(in) :: dense(:, :)
      real(dp), allocatable, intent(out) :: values(:)

      real(dp), allocatable :: work(:, :)
      integer :: info

      allocate (work(size(dense, 1), size(dense, 2)), values(size(dense, 1)))
      work = dense
      call pic_syev(work, values, jobz="N", uplo="U", info=info)
      deallocate (work)
   end subroutine dense_spectrum

   subroutine test_ground_state(error)
      !! The lowest eigenvalue, to the residual tolerance asked for
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(davidson_result_t) :: result
      type(link_table_t) :: alpha, beta
      real(dp), allocatable :: folded(:, :), diagonal(:, :), dense(:, :), values(:)
      logical :: ok

      call model(3, 3, folded, diagonal, alpha, beta, dense, err, ok)
      call check(error, ok, "the model should build")
      if (allocated(error)) return
      call dense_spectrum(dense, values)

      call davidson_lowest(folded, diagonal, alpha, beta, 1, result, err)
      call check(error,.not. err%has_error(), "the solver should run")
      if (allocated(error)) return
      call check(error, result%converged, "it should converge")
      if (allocated(error)) return
      call check(error, result%values(1), values(1), &
                 "the lowest eigenvalue should be the lowest eigenvalue", &
                 thr=1.0e-11_dp)
      if (allocated(error)) return

      ! It should take far fewer matrix-vector products than the space has
      ! dimensions, or there was no point building it. 400 determinants here.
      call check(error, result%sigma_products < 100, &
                 "the whole point is convergence in far fewer sigma products "// &
                 "than the determinant space has dimensions")
      if (allocated(error)) return
      call alpha%destroy()
      call beta%destroy()
   end subroutine test_ground_state

   subroutine test_several_roots(error)
      !! Four roots at once, all matching the dense spectrum
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(davidson_result_t) :: result
      type(link_table_t) :: alpha, beta
      real(dp), allocatable :: folded(:, :), diagonal(:, :), dense(:, :), values(:)
      logical :: ok
      integer :: i

      call model(3, 3, folded, diagonal, alpha, beta, dense, err, ok)
      call check(error, ok, "the model should build")
      if (allocated(error)) return
      call dense_spectrum(dense, values)

      call davidson_lowest(folded, diagonal, alpha, beta, 4, result, err)
      call check(error,.not. err%has_error(), "the solver should run")
      if (allocated(error)) return
      call check(error, result%converged, "four roots should converge")
      if (allocated(error)) return

      ! Ascending, and each one the right eigenvalue. Excited roots are the
      ! harder case: the ground state can be found by almost any descent, while
      ! a wrong subspace orthogonalisation shows up first in root four.
      do i = 1, 4
         call check(error, result%values(i), values(i), &
                    "root "//char(48 + i)//" against the dense spectrum", &
                    thr=1.0e-11_dp)
         if (allocated(error)) return
         if (i > 1) then
            call check(error, result%values(i) >= result%values(i - 1), &
                       "roots should come back ascending")
            if (allocated(error)) return
         end if
      end do
      call alpha%destroy()
      call beta%destroy()
   end subroutine test_several_roots

   subroutine test_vectors(error)
      !! The vectors are normalised, orthogonal, and satisfy H c = E c
      !!
      !! The eigenvalue can be right while the vector is not -- a Rayleigh
      !! quotient is stationary, so a vector with a small error gives an energy
      !! with a much smaller one. Every density matrix downstream is built from
      !! the vector, not the energy, so the vector is what has to be checked.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(davidson_result_t) :: result
      type(link_table_t) :: alpha, beta
      real(dp), allocatable :: folded(:, :), diagonal(:, :), dense(:, :)
      real(dp), allocatable :: flat(:), image(:)
      logical :: ok
      integer :: na, nb, ndet, i, j

      call model(3, 3, folded, diagonal, alpha, beta, dense, err, ok)
      call check(error, ok, "the model should build")
      if (allocated(error)) return
      na = alpha%n_strings
      nb = beta%n_strings
      ndet = na*nb

      call davidson_lowest(folded, diagonal, alpha, beta, 3, result, err)
      call check(error, result%converged .and. .not. err%has_error(), &
                 "three roots should converge")
      if (allocated(error)) return

      do i = 1, 3
         flat = reshape(result%vectors(:, :, i), [ndet])
         call check(error, abs(dot_product(flat, flat) - 1.0_dp) < 1.0e-10_dp, &
                    "each eigenvector should be normalised")
         if (allocated(error)) return

         image = matmul(dense, flat)
         call check(error, maxval(abs(image - result%values(i)*flat)) < 1.0e-9_dp, &
                    "H c should equal E c, element by element")
         if (allocated(error)) return

         do j = 1, i - 1
            call check(error, abs(dot_product(reshape(result%vectors(:, :, j), [ndet]), &
                                              flat)) < 1.0e-9_dp, &
                       "eigenvectors of different roots should be orthogonal")
            if (allocated(error)) return
         end do
      end do
      call alpha%destroy()
      call beta%destroy()
   end subroutine test_vectors

   subroutine test_collapse(error)
      !! A subspace too small to grow into still converges
      !!
      !! Capped at the smallest subspace that can hold the roots plus a little,
      !! so the solver has to collapse onto its Ritz vectors and restart
      !! repeatedly. The answer must not depend on that: collapsing is a
      !! memory strategy, not an approximation. It is also the path least
      !! likely to be taken by a small test and most likely to be wrong,
      !! because the sigma vectors have to be carried through the collapse
      !! rather than recomputed.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(davidson_result_t) :: roomy, cramped
      type(link_table_t) :: alpha, beta
      real(dp), allocatable :: folded(:, :), diagonal(:, :), dense(:, :)
      logical :: ok
      integer :: i

      call model(3, 3, folded, diagonal, alpha, beta, dense, err, ok)
      call check(error, ok, "the model should build")
      if (allocated(error)) return

      call davidson_lowest(folded, diagonal, alpha, beta, 2, roomy, err)
      call davidson_lowest(folded, diagonal, alpha, beta, 2, cramped, err, &
                           max_subspace=5)
      call check(error,.not. err%has_error(), "both should run")
      if (allocated(error)) return
      call check(error, cramped%converged, "the cramped solve should still converge")
      if (allocated(error)) return

      do i = 1, 2
         call check(error, cramped%values(i), roomy%values(i), &
                    "collapsing the subspace must not change the answer", &
                    thr=1.0e-10_dp)
         if (allocated(error)) return
      end do

      ! It should have cost more, or the cap did nothing and the test is empty.
      call check(error, cramped%sigma_products > roomy%sigma_products, &
                 "a cramped subspace should need more matrix-vector products, "// &
                 "otherwise the collapse path was never taken")
      call alpha%destroy()
      call beta%destroy()
   end subroutine test_collapse

   subroutine test_guess(error)
      !! A supplied starting vector reaches the same answer
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(davidson_result_t) :: plain, guided
      type(link_table_t) :: alpha, beta
      real(dp), allocatable :: folded(:, :), diagonal(:, :), dense(:, :)
      real(dp), allocatable :: guess(:, :, :)
      logical :: ok

      call model(3, 3, folded, diagonal, alpha, beta, dense, err, ok)
      call check(error, ok, "the model should build")
      if (allocated(error)) return

      call davidson_lowest(folded, diagonal, alpha, beta, 1, plain, err)
      call check(error, plain%converged, "the default start should converge")
      if (allocated(error)) return

      ! Start from the converged answer. It should be recognised as converged
      ! almost at once, which is what a macro-iteration of CASSCF depends on --
      ! the CI is re-solved every time the orbitals move, from an answer that is
      ! nearly right.
      allocate (guess(alpha%n_strings, beta%n_strings, 1))
      guess(:, :, 1) = plain%vectors(:, :, 1)
      call davidson_lowest(folded, diagonal, alpha, beta, 1, guided, err, guess=guess)
      call check(error,.not. err%has_error(), "the guided solve should run")
      if (allocated(error)) return
      call check(error, guided%converged, "and converge")
      if (allocated(error)) return
      call check(error, guided%values(1), plain%values(1), &
                 "to the same eigenvalue", thr=1.0e-11_dp)
      if (allocated(error)) return
      call check(error, guided%sigma_products < plain%sigma_products, &
                 "starting from the answer should cost less than starting from a "// &
                 "unit vector")
      call alpha%destroy()
      call beta%destroy()
   end subroutine test_guess

   subroutine test_refusals(error)
      !! Requests that do not describe an eigenproblem
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(davidson_result_t) :: result
      type(link_table_t) :: alpha, beta
      real(dp), allocatable :: folded(:, :), diagonal(:, :), dense(:, :)
      logical :: ok

      call model(3, 3, folded, diagonal, alpha, beta, dense, err, ok)
      call check(error, ok, "the model should build")
      if (allocated(error)) return

      call davidson_lowest(folded, diagonal, alpha, beta, 0, result, err)
      call check(error, err%has_error(), "zero roots should be refused")
      if (allocated(error)) return
      call err%clear()

      call davidson_lowest(folded, diagonal, alpha, beta, 10, result, err, &
                           max_subspace=4)
      call check(error, err%has_error(), &
                 "a subspace too small to hold the roots asked for should be refused")
      call err%clear()
      call alpha%destroy()
      call beta%destroy()
   end subroutine test_refusals

end module test_mqc_davidson

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_davidson, only: collect_mqc_davidson_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_davidson", collect_mqc_davidson_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
