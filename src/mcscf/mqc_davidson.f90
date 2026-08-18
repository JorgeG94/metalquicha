!! The lowest eigenpairs of the CI Hamiltonian, without the Hamiltonian
module mqc_davidson
   !! Davidson's method: build a small subspace, diagonalise the Hamiltonian
   !! inside it, and grow it in the direction the residual points. The only
   !! thing it asks of the operator is a matrix-vector product, which is exactly
   !! what `sigma_vector` provides and all that is affordable -- the matrix
   !! whose lowest eigenvalue is wanted is 32 GB at CAS(10,10) and unthinkable
   !! past that, while the subspace here is a few dozen vectors.
   !!
   !! **The preconditioner is the whole method.** Without it this is Lanczos
   !! with extra steps and converges at the rate the spectrum dictates. Dividing
   !! the residual by `theta - H_DD` before adding it approximates the
   !! correction that would solve the problem exactly if the Hamiltonian were
   !! diagonal, and for a CI Hamiltonian -- strongly diagonally dominant, since
   !! a determinant far from the reference is far in energy too -- that is a
   !! good enough approximation to converge in tens of iterations rather than
   !! thousands.
   !!
   !! A wrong preconditioner does not give a wrong answer. It gives the right
   !! answer slowly, or no answer at all, and neither symptom points back here;
   !! that is why `ci_diagonal` is checked against the assembled matrix
   !! separately rather than trusted because the solver seems to work.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_determinants, only: link_table_t
   use mqc_ci, only: sigma_vector
   implicit none
   private

   public :: davidson_lowest
   public :: davidson_result_t

   real(dp), parameter :: DEFAULT_TOLERANCE = 1.0e-10_dp
      !! On the residual norm, not the energy. The energy error is second order
      !! in the vector error near an eigenvector, so an energy-based test stops
      !! about half way through the digits the vector has -- and the vector is
      !! what the density matrices are built from.
   integer, parameter :: DEFAULT_MAX_ITERATIONS = 200
   real(dp), parameter :: PRECONDITION_FLOOR = 1.0e-8_dp
      !! A determinant whose diagonal sits on the current eigenvalue estimate
      !! would divide by nothing. It happens: the reference determinant's own
      !! diagonal approaches the ground state energy as the expansion converges.
   real(dp), parameter :: LINEAR_DEPENDENCE = 1.0e-8_dp
      !! How much of a new direction has to survive projection out of the
      !! subspace for it to be worth adding, as a fraction of its own length.
      !!
      !! **The fraction is the point.** Tested against the raw correction
      !! instead, this threshold scales with the residual, so it stops being a
      !! statement about linear dependence and becomes a floor on the residual:
      !! near convergence every correction looks tiny and gets dropped, the
      !! subspace stops growing, and the solve halts a long way short of the
      !! tolerance it was asked for. Normalising first makes the test scale-free
      !! and is the difference between converging to 1e-10 and stalling at 1e-6.

   type :: davidson_result_t
      real(dp), allocatable :: values(:)       !! (n_roots), ascending
      real(dp), allocatable :: vectors(:, :, :)
         !! (n_alpha, n_beta, n_roots), each normalised
      real(dp), allocatable :: residuals(:)    !! (n_roots), final norms
      integer :: iterations = 0
      integer :: sigma_products = 0            !! What the cost is actually made of
      logical :: converged = .false.
   end type davidson_result_t

contains

   subroutine davidson_lowest(folded, diagonal, alpha, beta, n_roots, result, error, &
                              tolerance, max_iterations, max_subspace, guess)
      !! The `n_roots` lowest eigenpairs of the CI Hamiltonian
      real(dp), intent(in) :: folded(:, :)      !! From `absorb_one_electron`
      real(dp), intent(in) :: diagonal(:, :)    !! From `ci_diagonal`
      type(link_table_t), intent(in) :: alpha, beta
      integer, intent(in) :: n_roots
      type(davidson_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      integer, intent(in), optional :: max_subspace
         !! Vectors kept before the subspace is collapsed onto the current Ritz
         !! vectors. Defaults to `max(2*n_roots + 8, 16)`, bounded by the size of
         !! the determinant space.
      real(dp), intent(in), optional :: guess(:, :, :)
         !! (n_alpha, n_beta, n_roots) starting vectors. Absent takes unit
         !! vectors on the lowest diagonal elements, which for a CI is the
         !! reference determinant and its nearest neighbours in energy.

      real(dp), allocatable :: basis(:, :), sigma(:, :), small(:, :), small_values(:)
      real(dp), allocatable :: ritz(:, :), hritz(:, :), residual(:), correction(:)
      real(dp), allocatable :: work(:, :), sigma_work(:, :)
      real(dp) :: tol, norm, denominator, overlap
      integer :: na, nb, ndet, nsub, nmax, iterations, iroot, i, j, info, added
      integer :: iteration, pass
      logical, allocatable :: root_converged(:)

      if (error%has_error()) return
      na = alpha%n_strings
      nb = beta%n_strings
      ndet = na*nb

      tol = DEFAULT_TOLERANCE
      if (present(tolerance)) tol = tolerance
      iterations = DEFAULT_MAX_ITERATIONS
      if (present(max_iterations)) iterations = max_iterations
      nmax = max(2*n_roots + 8, 16)
      if (present(max_subspace)) nmax = max_subspace
      nmax = min(nmax, ndet)

      if (n_roots < 1 .or. n_roots > ndet) then
         call error%set(ERROR_VALIDATION, "asked for "//to_char(n_roots)// &
                        " roots of a Hamiltonian with "//to_char(ndet)// &
                        " determinants.")
         return
      end if
      if (nmax < n_roots) then
         call error%set(ERROR_VALIDATION, "the subspace is capped at "// &
                        to_char(nmax)//" vectors, which cannot hold "// &
                        to_char(n_roots)//" roots.")
         return
      end if

      allocate (basis(ndet, nmax), sigma(ndet, nmax))
      allocate (ritz(ndet, n_roots), hritz(ndet, n_roots))
      allocate (residual(ndet), correction(ndet))
      allocate (work(na, nb), sigma_work(na, nb))
      allocate (root_converged(n_roots))
      allocate (result%values(n_roots), result%residuals(n_roots))
      allocate (result%vectors(na, nb, n_roots))

      call initial_basis(diagonal, n_roots, na, nb, ndet, basis, guess)
      nsub = n_roots

      ! Sigma for the starting vectors.
      do i = 1, nsub
         work = reshape(basis(:, i), [na, nb])
         call sigma_vector(folded, work, alpha, beta, sigma_work, error)
         if (error%has_error()) return
         sigma(:, i) = reshape(sigma_work, [ndet])
         result%sigma_products = result%sigma_products + 1
      end do

      do iteration = 1, iterations
         result%iterations = iteration

         ! The Hamiltonian projected into the subspace, and its eigenpairs.
         allocate (small(nsub, nsub), small_values(nsub))
         call pic_gemm(basis(:, 1:nsub), sigma(:, 1:nsub), small, transa="T")
         ! Symmetric to rounding; the average is what a symmetric solver would
         ! read anyway and keeps a drifting subspace from biasing a root.
         small = 0.5_dp*(small + transpose(small))
         call pic_syev(small, small_values, jobz="V", uplo="U", info=info)
         if (info /= 0) then
            call error%set(ERROR_VALIDATION, "the Davidson subspace matrix could "// &
                           "not be diagonalized (info = "//to_char(info)//")")
            return
         end if

         ! Ritz vectors and their images, both from the subspace, so no new
         ! sigma products are needed.
         call pic_gemm(basis(:, 1:nsub), small(:, 1:n_roots), ritz)
         call pic_gemm(sigma(:, 1:nsub), small(:, 1:n_roots), hritz)
         result%values = small_values(1:n_roots)

         do iroot = 1, n_roots
            residual = hritz(:, iroot) - small_values(iroot)*ritz(:, iroot)
            result%residuals(iroot) = sqrt(dot_product(residual, residual))
            root_converged(iroot) = result%residuals(iroot) < tol
         end do

         if (all(root_converged)) then
            result%converged = .true.
            do iroot = 1, n_roots
               result%vectors(:, :, iroot) = reshape(ritz(:, iroot), [na, nb])
            end do
            deallocate (small, small_values)
            exit
         end if

         ! Collapse before growing, so the expansion vectors below are added to
         ! a fresh subspace rather than discarded by it.
         if (nsub + count(.not. root_converged) > nmax) then
            basis(:, 1:n_roots) = ritz
            sigma(:, 1:n_roots) = hritz
            nsub = n_roots
         end if

         added = 0
         do iroot = 1, n_roots
            if (root_converged(iroot)) cycle
            if (nsub + added >= nmax) exit

            residual = hritz(:, iroot) - small_values(iroot)*ritz(:, iroot)
            do i = 1, ndet
               denominator = small_values(iroot) - diagonal(mod(i - 1, na) + 1, &
                                                            (i - 1)/na + 1)
               if (abs(denominator) < PRECONDITION_FLOOR) then
                  denominator = sign(PRECONDITION_FLOOR, denominator)
               end if
               correction(i) = residual(i)/denominator
            end do

            ! Normalise before projecting, so the linear-dependence test below
            ! measures direction rather than magnitude.
            norm = sqrt(dot_product(correction, correction))
            if (norm < tiny(1.0_dp)) cycle
            correction = correction/norm

            ! Twice, because one pass of Gram-Schmidt loses orthogonality when
            ! the new direction is nearly in the subspace already -- which near
            ! convergence it always is.
            do pass = 1, 2
               do j = 1, nsub + added
                  overlap = dot_product(basis(:, j), correction)
                  correction = correction - overlap*basis(:, j)
               end do
            end do
            norm = sqrt(dot_product(correction, correction))
            if (norm < LINEAR_DEPENDENCE) cycle

            added = added + 1
            basis(:, nsub + added) = correction/norm
            work = reshape(basis(:, nsub + added), [na, nb])
            call sigma_vector(folded, work, alpha, beta, sigma_work, error)
            if (error%has_error()) return
            sigma(:, nsub + added) = reshape(sigma_work, [ndet])
            result%sigma_products = result%sigma_products + 1
         end do

         deallocate (small, small_values)

         ! Nothing left to expand with: every direction the residuals point in
         ! is already spanned. The answer is as converged as this subspace can
         ! make it, which for a full space is exactly converged.
         if (added == 0) then
            do iroot = 1, n_roots
               result%vectors(:, :, iroot) = reshape(ritz(:, iroot), [na, nb])
            end do
            result%converged = all(root_converged)
            exit
         end if
         nsub = nsub + added
      end do

      if (.not. result%converged) then
         do iroot = 1, n_roots
            result%vectors(:, :, iroot) = reshape(ritz(:, iroot), [na, nb])
         end do
      end if

      deallocate (basis, sigma, ritz, hritz, residual, correction, work, sigma_work)
      deallocate (root_converged)
   end subroutine davidson_lowest

   subroutine initial_basis(diagonal, n_roots, na, nb, ndet, basis, guess)
      !! Starting vectors: the supplied ones, or the lowest determinants
      real(dp), intent(in) :: diagonal(:, :)
      integer, intent(in) :: n_roots, na, nb, ndet
      real(dp), intent(inout) :: basis(:, :)
      real(dp), intent(in), optional :: guess(:, :, :)

      real(dp), allocatable :: flat(:)
      logical, allocatable :: taken(:)
      real(dp) :: best, norm, overlap
      integer :: iroot, i, pick, j

      if (present(guess)) then
         do iroot = 1, n_roots
            basis(:, iroot) = reshape(guess(:, :, iroot), [ndet])
            do j = 1, iroot - 1
               overlap = dot_product(basis(:, j), basis(:, iroot))
               basis(:, iroot) = basis(:, iroot) - overlap*basis(:, j)
            end do
            norm = sqrt(dot_product(basis(:, iroot), basis(:, iroot)))
            if (norm > LINEAR_DEPENDENCE) basis(:, iroot) = basis(:, iroot)/norm
         end do
         return
      end if

      allocate (flat(ndet), taken(ndet))
      flat = reshape(diagonal, [ndet])
      taken = .false.
      basis(:, 1:n_roots) = 0.0_dp
      do iroot = 1, n_roots
         pick = 0
         best = huge(1.0_dp)
         do i = 1, ndet
            if (taken(i)) cycle
            if (flat(i) < best) then
               best = flat(i)
               pick = i
            end if
         end do
         taken(pick) = .true.
         basis(pick, iroot) = 1.0_dp
      end do
      deallocate (flat, taken)
   end subroutine initial_basis

end module mqc_davidson
