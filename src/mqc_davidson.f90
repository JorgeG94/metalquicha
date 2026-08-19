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
   public :: davidson_flat
   public :: davidson_result_t
   public :: sigma_operator_t

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

   type, abstract :: sigma_operator_t
      !! Whatever can multiply a vector by the Hamiltonian
      !!
      !! Davidson asks the operator for one thing and the determinant space for
      !! another: a matrix-vector product, and a diagonal to precondition with.
      !! Neither needs the vector to be shaped, and a restricted active space
      !! cannot shape it -- its determinants are not a rectangle. So the method
      !! is written against a flat vector and this type, and the two spaces
      !! differ only in what they put behind it.
   contains
      procedure(apply_sigma), deferred :: apply
   end type sigma_operator_t

   abstract interface
      subroutine apply_sigma(this, vector, image, error)
         import :: sigma_operator_t, dp, error_t
         implicit none
         class(sigma_operator_t), intent(inout) :: this
         real(dp), intent(in) :: vector(:)
         real(dp), intent(out) :: image(:)
         type(error_t), intent(inout) :: error
      end subroutine apply_sigma
   end interface

   type, extends(sigma_operator_t) :: cas_operator_t
      !! The complete-active-space Hamiltonian, as an operator
      real(dp), pointer :: folded(:, :) => null()
      type(link_table_t), pointer :: alpha => null(), beta => null()
      integer :: na = 0, nb = 0
   contains
      procedure :: apply => cas_apply
   end type cas_operator_t

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

   subroutine cas_apply(this, vector, image, error)
      !! `sigma = H c` for a complete active space, shaping and unshaping
      class(cas_operator_t), intent(inout) :: this
      real(dp), intent(in) :: vector(:)
      real(dp), intent(out) :: image(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: work(:, :), sigma_work(:, :)

      allocate (work(this%na, this%nb), sigma_work(this%na, this%nb))
      work = reshape(vector, [this%na, this%nb])
      call sigma_vector(this%folded, work, this%alpha, this%beta, sigma_work, error)
      if (.not. error%has_error()) image = reshape(sigma_work, [this%na*this%nb])
      deallocate (work, sigma_work)
   end subroutine cas_apply

   subroutine davidson_lowest(folded, diagonal, alpha, beta, n_roots, result, error, &
                              tolerance, max_iterations, max_subspace, guess)
      !! The `n_roots` lowest eigenpairs of the CI Hamiltonian
      !!
      !! The complete-space spelling of `davidson_flat`: the same method, with
      !! the vector shaped back into the alpha-by-beta rectangle its
      !! determinants form.
      real(dp), intent(in), target :: folded(:, :)      !! From `absorb_one_electron`
      real(dp), intent(in) :: diagonal(:, :)            !! From `ci_diagonal`
      type(link_table_t), intent(in), target :: alpha, beta
      integer, intent(in) :: n_roots
      type(davidson_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      integer, intent(in), optional :: max_subspace
      real(dp), intent(in), optional :: guess(:, :, :)
         !! (n_alpha, n_beta, n_roots) starting vectors

      type(cas_operator_t) :: operator
      real(dp), allocatable :: vectors(:, :), flat_guess(:, :)
      integer :: na, nb, ndet, iroot

      if (error%has_error()) return
      na = alpha%n_strings
      nb = beta%n_strings
      ndet = na*nb

      operator%folded => folded
      operator%alpha => alpha
      operator%beta => beta
      operator%na = na
      operator%nb = nb

      if (present(guess)) then
         allocate (flat_guess(ndet, size(guess, 3)))
         do iroot = 1, size(guess, 3)
            flat_guess(:, iroot) = reshape(guess(:, :, iroot), [ndet])
         end do
         call davidson_flat(operator, reshape(diagonal, [ndet]), n_roots, &
                            result%values, vectors, result%residuals, &
                            result%iterations, result%sigma_products, &
                            result%converged, error, tolerance, max_iterations, &
                            max_subspace, flat_guess)
      else
         call davidson_flat(operator, reshape(diagonal, [ndet]), n_roots, &
                            result%values, vectors, result%residuals, &
                            result%iterations, result%sigma_products, &
                            result%converged, error, tolerance, max_iterations, &
                            max_subspace)
      end if
      if (error%has_error()) return

      allocate (result%vectors(na, nb, n_roots))
      do iroot = 1, n_roots
         result%vectors(:, :, iroot) = reshape(vectors(:, iroot), [na, nb])
      end do
   end subroutine davidson_lowest

   subroutine davidson_flat(operator, diagonal, n_roots, values, vectors, residuals, &
                            iterations_taken, sigma_products, converged, error, &
                            tolerance, max_iterations, max_subspace, guess)
      !! The `n_roots` lowest eigenpairs of anything that can multiply a vector
      !!
      !! The method itself, over a flat vector. Everything that knows what a
      !! determinant is lives in `operator` and in `diagonal`.
      class(sigma_operator_t), intent(inout) :: operator
      real(dp), intent(in) :: diagonal(:)       !! `<D|H|D>` for every determinant
      integer, intent(in) :: n_roots
      real(dp), allocatable, intent(out) :: values(:), vectors(:, :), residuals(:)
      integer, intent(out) :: iterations_taken, sigma_products
      logical, intent(out) :: converged
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      integer, intent(in), optional :: max_subspace
         !! Vectors kept before the subspace is collapsed onto the current Ritz
         !! vectors. Defaults to `max(2*n_roots + 8, 16)`, bounded by the size of
         !! the determinant space.
      real(dp), intent(in), optional :: guess(:, :)
         !! (n_determinants, n_roots) starting vectors. Absent takes unit vectors
         !! on the lowest diagonal elements, which for a CI is the reference
         !! determinant and its nearest neighbours in energy.

      real(dp), allocatable :: basis(:, :), sigma(:, :), small(:, :), small_values(:)
      real(dp), allocatable :: ritz(:, :), hritz(:, :), residual(:), correction(:)
      real(dp) :: tol, norm, denominator, overlap
      integer :: ndet, nsub, nmax, iterations, iroot, i, j, info, added
      integer :: iteration, pass
      logical, allocatable :: root_converged(:)

      if (error%has_error()) return
      ndet = size(diagonal)
      iterations_taken = 0
      sigma_products = 0
      converged = .false.

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
      allocate (root_converged(n_roots))
      allocate (values(n_roots), residuals(n_roots))
      allocate (vectors(ndet, n_roots))

      call initial_basis(diagonal, n_roots, ndet, basis, guess)
      nsub = n_roots

      ! Sigma for the starting vectors.
      do i = 1, nsub
         call operator%apply(basis(:, i), sigma(:, i), error)
         if (error%has_error()) return
         sigma_products = sigma_products + 1
      end do

      do iteration = 1, iterations
         iterations_taken = iteration

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
         values = small_values(1:n_roots)

         do iroot = 1, n_roots
            residual = hritz(:, iroot) - small_values(iroot)*ritz(:, iroot)
            residuals(iroot) = sqrt(dot_product(residual, residual))
            root_converged(iroot) = residuals(iroot) < tol
         end do

         if (all(root_converged)) then
            converged = .true.
            vectors = ritz
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
               denominator = small_values(iroot) - diagonal(i)
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
            call operator%apply(basis(:, nsub + added), sigma(:, nsub + added), error)
            if (error%has_error()) return
            sigma_products = sigma_products + 1
         end do

         deallocate (small, small_values)

         ! Nothing left to expand with: every direction the residuals point in
         ! is already spanned. The answer is as converged as this subspace can
         ! make it, which for a full space is exactly converged.
         if (added == 0) then
            vectors = ritz
            converged = all(root_converged)
            exit
         end if
         nsub = nsub + added
      end do

      if (.not. converged) vectors = ritz

      deallocate (basis, sigma, ritz, hritz, residual, correction)
      deallocate (root_converged)
   end subroutine davidson_flat

   subroutine initial_basis(diagonal, n_roots, ndet, basis, guess)
      !! Starting vectors: the supplied ones, or the lowest determinants
      real(dp), intent(in) :: diagonal(:)
      integer, intent(in) :: n_roots, ndet
      real(dp), intent(inout) :: basis(:, :)
      real(dp), intent(in), optional :: guess(:, :)

      logical, allocatable :: taken(:)
      real(dp) :: best, norm, overlap
      integer :: iroot, i, pick, j

      if (present(guess)) then
         do iroot = 1, n_roots
            basis(:, iroot) = guess(:, iroot)
            do j = 1, iroot - 1
               overlap = dot_product(basis(:, j), basis(:, iroot))
               basis(:, iroot) = basis(:, iroot) - overlap*basis(:, j)
            end do
            norm = sqrt(dot_product(basis(:, iroot), basis(:, iroot)))
            if (norm > LINEAR_DEPENDENCE) basis(:, iroot) = basis(:, iroot)/norm
         end do
         return
      end if

      allocate (taken(ndet))
      taken = .false.
      basis(:, 1:n_roots) = 0.0_dp
      do iroot = 1, n_roots
         pick = 0
         best = huge(1.0_dp)
         do i = 1, ndet
            if (taken(i)) cycle
            if (diagonal(i) < best) then
               best = diagonal(i)
               pick = i
            end if
         end do
         taken(pick) = .true.
         basis(pick, iroot) = 1.0_dp
      end do
      deallocate (taken)
   end subroutine initial_basis

end module mqc_davidson
