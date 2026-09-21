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
   use mqc_davidson, only: davidson_lowest, davidson_result_t, davidson_flat, &
                           sigma_operator_t
   implicit none
   private

   public :: collect_mqc_davidson_tests

   integer, parameter :: NORB = 6
   integer, parameter :: NCHOL = 3
   integer, parameter :: NDENSE = 40
      !! Dimension of the dense model the block-product cases run on. Small
      !! enough for a `matmul` per product, large enough that the subspace
      !! never reaches it and the solver has to iterate.
   integer, parameter :: MAX_BLOCK_CALLS = 256
      !! How many block widths the counting operator has room to record. A
      !! solve needing more than this has already failed the test it is in.

   type, extends(sigma_operator_t) :: dense_operator_t
      !! A dense symmetric matrix as an operator, keeping count of what it is
      !! asked for
      !!
      !! It inherits the default `apply_many`, so it is also the case that
      !! proves the default still works.
      real(dp), allocatable :: matrix(:, :)
      integer :: apply_calls = 0
         !! Single-vector products performed, whichever route reached them
      integer :: block_calls = 0
         !! `apply_many` calls, which only an override is in a position to count
      integer :: width(MAX_BLOCK_CALLS) = 0
         !! Columns handed to each `apply_many` call, in call order
   contains
      procedure :: apply => dense_apply
   end type dense_operator_t

   type, extends(dense_operator_t) :: blocked_operator_t
      !! The same matrix, with the block product overridden rather than
      !! inherited
   contains
      procedure :: apply_many => blocked_apply_many
   end type blocked_operator_t

contains

   subroutine dense_apply(this, vector, image, error)
      !! One matrix-vector product, and one more on the counter
      class(dense_operator_t), intent(inout) :: this
      real(dp), intent(in) :: vector(:)
      real(dp), intent(out) :: image(:)
      type(error_t), intent(inout) :: error

      if (error%has_error()) return
      this%apply_calls = this%apply_calls + 1
      image = matmul(this%matrix, vector)
   end subroutine dense_apply

   subroutine blocked_apply_many(this, vectors, images, error)
      !! The whole block in one call, recording how wide it was
      !!
      !! The arithmetic inside is deliberately the single-vector arithmetic,
      !! column by column. What this override exists to test is that the
      !! solver hands its new vectors over in one call; keeping the floating
      !! point identical on both sides is what lets the eigenpairs of the
      !! overriding and the inheriting operator be compared for equality
      !! rather than to a threshold.
      class(blocked_operator_t), intent(inout) :: this
      real(dp), intent(in) :: vectors(:, :)
      real(dp), intent(out) :: images(:, :)
      type(error_t), intent(inout) :: error

      integer :: i

      if (error%has_error()) return
      this%block_calls = this%block_calls + 1
      if (this%block_calls <= MAX_BLOCK_CALLS) then
         this%width(this%block_calls) = size(vectors, 2)
      end if
      do i = 1, size(vectors, 2)
         call dense_apply(this, vectors(:, i), images(:, i), error)
         if (error%has_error()) return
      end do
   end subroutine blocked_apply_many

   subroutine dense_model(matrix, diagonal)
      !! A diagonally dominant symmetric matrix, and the diagonal to
      !! precondition with
      real(dp), allocatable, intent(out) :: matrix(:, :), diagonal(:)

      integer :: i, j

      allocate (matrix(NDENSE, NDENSE), diagonal(NDENSE))
      do j = 1, NDENSE
         do i = 1, NDENSE
            if (i == j) then
               matrix(i, j) = real(i, dp)
            else
               matrix(i, j) = 0.15_dp/real(abs(i - j) + 1, dp)
            end if
         end do
      end do
      do i = 1, NDENSE
         diagonal(i) = matrix(i, i)
      end do
   end subroutine dense_model

   subroutine collect_mqc_davidson_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ground_state_matches_dense", test_ground_state), &
                  new_unittest("several_roots_match_dense", test_several_roots), &
                  new_unittest("eigenvectors_are_eigenvectors", test_vectors), &
                  new_unittest("subspace_collapse", test_collapse), &
                  new_unittest("guess_is_used", test_guess), &
                  new_unittest("dependent_guess_is_replaced", test_dependent_guess), &
                  new_unittest("guess_scale_is_not_dependence", test_guess_scale), &
                  new_unittest("one_block_call_per_iteration", test_block_calls), &
                  new_unittest("block_matches_loop", test_block_matches_loop), &
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

   subroutine test_dependent_guess(error)
      !! A guess whose second vector repeats the first still solves
      !!
      !! Orthogonalising the repeat against the original leaves nothing, and a
      !! subspace started on a near-null column produces Ritz vectors that mean
      !! nothing -- so the column has to be replaced rather than kept. A CASSCF
      !! macro-iteration hands the CI the previous cycle's vectors, and two of
      !! those can collapse onto each other when roots cross.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(dense_operator_t) :: plain, guided
      real(dp), allocatable :: diagonal(:), guess_diagonal(:), guess(:, :)
      real(dp), allocatable :: plain_values(:), plain_vectors(:, :), plain_residuals(:)
      real(dp), allocatable :: values(:), vectors(:, :), residuals(:)
      integer :: iterations, products, i
      logical :: converged

      call dense_model(plain%matrix, diagonal)
      call dense_model(guided%matrix, guess_diagonal)

      call davidson_flat(plain, diagonal, 2, plain_values, plain_vectors, &
                         plain_residuals, iterations, products, converged, err)
      call check(error, converged .and. .not. err%has_error(), &
                 "the unguided solve should converge")
      if (allocated(error)) return

      ! Both columns the same vector, and not one of the unit vectors the
      ! fallback picks, so the replacement has to survive orthogonalisation.
      allocate (guess(NDENSE, 2))
      guess(:, 1) = 1.0_dp/sqrt(real(NDENSE, dp))
      guess(:, 2) = guess(:, 1)
      call davidson_flat(guided, guess_diagonal, 2, values, vectors, residuals, &
                         iterations, products, converged, err, guess=guess)
      call check(error,.not. err%has_error(), "the guided solve should run")
      if (allocated(error)) return
      call check(error, converged, "and converge despite the repeated column")
      if (allocated(error)) return
      do i = 1, 2
         call check(error, values(i), plain_values(i), &
                    "to the eigenvalues the unguided solve found", thr=1.0e-10_dp)
         if (allocated(error)) return
      end do
   end subroutine test_dependent_guess

   subroutine test_guess_scale(error)
      !! The same guess at two magnitudes gives the same solve
      !!
      !! Whether a starting vector is worth keeping is a question about its
      !! direction, so the linear-dependence test asks what *fraction* of a
      !! column survives projection rather than how long what survives is.
      !! Scaling the whole guess must therefore change nothing. Tested with a
      !! factor small enough that an absolute test would throw the guess away
      !! -- `2**-40` is nine orders below `LINEAR_DEPENDENCE` -- and a power
      !! of two, so that every operation on the scaled guess is the unscaled
      !! one with the exponent moved and the comparison below can be equality
      !! rather than a threshold.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(dense_operator_t) :: plain, scaled
      real(dp), allocatable :: plain_diagonal(:), scaled_diagonal(:), guess(:, :)
      real(dp), allocatable :: plain_values(:), plain_vectors(:, :), plain_residuals(:)
      real(dp), allocatable :: values(:), vectors(:, :), residuals(:)
      integer :: plain_iterations, plain_products, iterations, products, i
      logical :: plain_converged, converged
      real(dp), parameter :: SMALL_SCALE = 2.0_dp**(-40)

      call dense_model(plain%matrix, plain_diagonal)
      call dense_model(scaled%matrix, scaled_diagonal)

      ! Two independent columns, neither of them a unit vector the fallback
      ! could reach for, and not orthogonal to each other either -- so the
      ! second one has to survive the projection on its own merits.
      allocate (guess(NDENSE, 2))
      guess = 0.0_dp
      guess(:, 1) = 1.0_dp/sqrt(real(NDENSE, dp))
      guess(1, 2) = 1.0_dp

      call davidson_flat(plain, plain_diagonal, 2, plain_values, plain_vectors, &
                         plain_residuals, plain_iterations, plain_products, &
                         plain_converged, err, guess=guess)
      call davidson_flat(scaled, scaled_diagonal, 2, values, vectors, residuals, &
                         iterations, products, converged, err, &
                         guess=SMALL_SCALE*guess)
      call check(error,.not. err%has_error(), "both solves should run")
      if (allocated(error)) return
      call check(error, plain_converged .and. converged, "and both converge")
      if (allocated(error)) return

      call check(error, iterations, plain_iterations, &
                 "a scaled guess should take the same number of iterations")
      if (allocated(error)) return
      call check(error, products, plain_products, &
                 "and the same number of sigma products")
      if (allocated(error)) return
      do i = 1, 2
         call check(error, values(i) == plain_values(i), &
                    "root "//char(48 + i)//" should be bit-identical")
         if (allocated(error)) return
      end do
      call check(error, all(vectors == plain_vectors), &
                 "and so should the eigenvectors")
   end subroutine test_guess_scale

   subroutine test_block_calls(error)
      !! One block product per iteration, and nothing applied outside a block
      !!
      !! This is the property a batched operator is built on: if the solver
      !! reached `apply` directly anywhere, the vectors it applied there would
      !! not be in any block, and an operator whose saving comes from sharing
      !! work across a block would silently lose it.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(blocked_operator_t) :: operator
      real(dp), allocatable :: diagonal(:), values(:), vectors(:, :), residuals(:)
      integer :: iterations, products, i
      logical :: converged

      call dense_model(operator%matrix, diagonal)
      call davidson_flat(operator, diagonal, 3, values, vectors, residuals, &
                         iterations, products, converged, err)
      call check(error,.not. err%has_error(), "the solver should run")
      if (allocated(error)) return
      call check(error, converged, "three roots should converge")
      if (allocated(error)) return

      ! One call before the loop for the starting vectors, then one per
      ! iteration that expands the subspace. The last iteration finds every
      ! root converged and leaves without expanding, so the two counts meet.
      call check(error, operator%block_calls, iterations, &
                 "the operator should be asked once per Davidson iteration")
      if (allocated(error)) return
      call check(error, operator%width(1), 3, &
                 "the first block is the starting vectors, one per root")
      if (allocated(error)) return
      call check(error, maxval(operator%width(1:operator%block_calls)), 3, &
                 "and an expansion block is as wide as the roots still running")
      if (allocated(error)) return
      do i = 1, operator%block_calls
         call check(error, operator%width(i) >= 1, "no block should be empty")
         if (allocated(error)) return
      end do
      call check(error, sum(operator%width(1:operator%block_calls)), products, &
                 "every sigma product should have arrived inside a block")
      if (allocated(error)) return
      call check(error, operator%apply_calls, products, &
                 "and the blocks should hold exactly those vectors")
   end subroutine test_block_calls

   subroutine test_block_matches_loop(error)
      !! Overriding the block product changes nothing about the answer
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(dense_operator_t) :: looped
      type(blocked_operator_t) :: blocked
      real(dp), allocatable :: loop_diagonal(:), block_diagonal(:)
      real(dp), allocatable :: loop_values(:), loop_vectors(:, :), loop_residuals(:)
      real(dp), allocatable :: block_values(:), block_vectors(:, :), block_residuals(:)
      integer :: loop_iterations, loop_products, block_iterations, block_products
      logical :: loop_converged, block_converged
      integer :: i

      call dense_model(looped%matrix, loop_diagonal)
      call dense_model(blocked%matrix, block_diagonal)

      call davidson_flat(looped, loop_diagonal, 3, loop_values, loop_vectors, &
                         loop_residuals, loop_iterations, loop_products, &
                         loop_converged, err)
      call davidson_flat(blocked, block_diagonal, 3, block_values, block_vectors, &
                         block_residuals, block_iterations, block_products, &
                         block_converged, err)
      call check(error,.not. err%has_error(), "both solves should run")
      if (allocated(error)) return
      call check(error, loop_converged .and. block_converged, "and both converge")
      if (allocated(error)) return

      ! `looped` inherits `apply_many`, so nothing it owns can ever increment
      ! `block_calls` and asserting it is zero asserts nothing. What the
      ! inherited route is on the hook for is one single-vector product per
      ! sigma product counted -- which is also what makes the comparison
      ! below one between two different routes rather than one route twice.
      call check(error, looped%apply_calls, loop_products, &
                 "the inherited default should apply one vector per sigma product")
      if (allocated(error)) return
      call check(error, blocked%block_calls > 0, "and the override is")
      if (allocated(error)) return
      call check(error, loop_iterations, block_iterations, &
                 "the two routes should take the same number of iterations")
      if (allocated(error)) return
      call check(error, loop_products, block_products, &
                 "and the same number of sigma products")
      if (allocated(error)) return

      ! Equality, not a threshold. The inherited loop and this override do the
      ! same arithmetic in the same order, so anything other than bit-identical
      ! output means the solver took a different path through one of them.
      do i = 1, 3
         call check(error, loop_values(i) == block_values(i), &
                    "root "//char(48 + i)//" should be bit-identical")
         if (allocated(error)) return
      end do
      call check(error, all(loop_vectors == block_vectors), &
                 "and so should the eigenvectors")
   end subroutine test_block_matches_loop

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
