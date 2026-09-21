!! The lowest eigenpairs of the CI Hamiltonian, without the Hamiltonian
module mqc_davidson
   !! Davidson's method: build a small subspace, diagonalise the Hamiltonian
   !! inside it, and grow it in the direction the residual points. The only
   !! thing it asks of the operator is a matrix-vector product, which is what
   !! `sigma_vector` provides.
   !!
   !! **The preconditioner is the whole method.** Without it this is Lanczos
   !! with extra steps. Dividing the residual by `theta - H_DD` before adding it
   !! is the correction that would be exact if the Hamiltonian were diagonal,
   !! and a CI Hamiltonian is diagonally dominant enough for that to converge in
   !! tens of iterations rather than thousands. A wrong preconditioner does not
   !! give a wrong answer, only a slow one or none at all, so `ci_diagonal` is
   !! checked separately rather than trusted because the solver seems to work.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_determinants, only: link_table_t
   use mqc_ci, only: sigma_vector
   use, intrinsic :: iso_fortran_env, only: output_unit, int64
   use pic_logger, only: logger => global_logger
   use mqc_convergence_report, only: convergence_header, convergence_footer
   implicit none
   private

   public :: davidson_lowest
   public :: davidson_flat
   public :: davidson_result_t
   public :: sigma_operator_t

   real(dp), parameter :: DEFAULT_TOLERANCE = 1.0e-10_dp
      !! On the residual norm, not the energy: near an eigenvector the energy
      !! error is second order in the vector error, so an energy test stops
      !! about half way through the digits the vector has, and the vector is
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
      !! **The fraction is the point.** Tested against the raw correction the
      !! threshold scales with the residual, so near convergence every
      !! correction looks tiny, the subspace stops growing, and the solve stalls
      !! short of its tolerance. Normalising first makes the test scale-free.
      !! `initial_basis` reaches the same property from the other side, by
      !! measuring what survives against the column's own length, so that a
      !! guess it means to keep is handed on with its arithmetic untouched.

   type, abstract :: sigma_operator_t
      !! Whatever can multiply a vector by the Hamiltonian
      !!
      !! Davidson needs a matrix-vector product and a diagonal to precondition
      !! with, neither of which needs the vector shaped -- and a restricted
      !! active space cannot shape it, its determinants not being a rectangle.
      !! Hence a flat vector and this type.
      !!
      !! The solver never applies a single vector on its own: it collects the
      !! vectors it wants and hands them over as a block through `apply_many`,
      !! whose default is a loop over `apply`. An operator whose cost is
      !! dominated by something a block can share -- a Fock build over several
      !! densities, a quadrature grid walked once -- overrides `apply_many` and
      !! gets that sharing without the solver knowing.
   contains
      procedure(apply_sigma), deferred :: apply
      procedure :: apply_many => apply_many_default
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
      integer :: sigma_products = 0
         !! Columns handed to the operator, summed over the block calls
         !!
         !! That is the cost only for an operator that inherits the default
         !! `apply_many`, which is a loop over single-vector products. An
         !! operator that overrides it pays once per call for whatever the
         !! block shares, so `iterations` is the count to read there: the
         !! solver makes one block call before the loop and one more per
         !! iteration that expands the subspace.
      logical :: converged = .false.
   end type davidson_result_t

contains

   subroutine apply_many_default(this, vectors, images, error)
      !! `sigma = H c` for a block of vectors, one column at a time
      !!
      !! What every operator gets unless it says otherwise, and the whole of
      !! the contract: column `i` of `images` is `apply` of column `i` of
      !! `vectors`, in that order. An override is free to compute the block
      !! however it likes, but it may not reorder the columns -- the solver
      !! pairs them by position.
      class(sigma_operator_t), intent(inout) :: this
      real(dp), intent(in) :: vectors(:, :)   !! (n_determinants, n_vectors)
      real(dp), intent(out) :: images(:, :)   !! (n_determinants, n_vectors)
      type(error_t), intent(inout) :: error

      integer :: i

      do i = 1, size(vectors, 2)
         call this%apply(vectors(:, i), images(:, i), error)
         if (error%has_error()) return
      end do
   end subroutine apply_many_default

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
                              tolerance, max_iterations, max_subspace, guess, verbose, energy_offset)
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
      logical, intent(in), optional :: verbose
         !! Print a line per iteration, off by default
      real(dp), intent(in), optional :: energy_offset
         !! Added to the eigenvalue before printing, and to nothing else. The
         !! Davidson solves the active-space problem, so its own eigenvalue is
         !! the active energy alone; the caller knows the inactive-plus-nuclear
         !! constant this adds back.

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
                            max_subspace, flat_guess, verbose, energy_offset)
      else
         call davidson_flat(operator, reshape(diagonal, [ndet]), n_roots, &
                            result%values, vectors, result%residuals, &
                            result%iterations, result%sigma_products, &
                            result%converged, error, tolerance, max_iterations, &
                            max_subspace, verbose=verbose, &
                            energy_offset=energy_offset)
      end if
      if (error%has_error()) return

      allocate (result%vectors(na, nb, n_roots))
      do iroot = 1, n_roots
         result%vectors(:, :, iroot) = reshape(vectors(:, iroot), [na, nb])
      end do
   end subroutine davidson_lowest

   subroutine davidson_flat(operator, diagonal, n_roots, values, vectors, residuals, &
                            iterations_taken, sigma_products, converged, error, &
                            tolerance, max_iterations, max_subspace, guess, verbose, &
                            energy_offset, label, value_label)
      !! The `n_roots` lowest eigenpairs of anything that can multiply a vector
      !!
      !! The method itself, over a flat vector. Everything that knows what a
      !! determinant is lives in `operator` and in `diagonal`.
      class(sigma_operator_t), intent(inout) :: operator
      real(dp), intent(in) :: diagonal(:)       !! `<D|H|D>` for every determinant
      integer, intent(in) :: n_roots
      real(dp), allocatable, intent(out) :: values(:), vectors(:, :), residuals(:)
      integer, intent(out) :: iterations_taken
      integer, intent(out) :: sigma_products
         !! Columns handed to the operator, not calls made to it. The two
         !! differ for an operator that overrides `apply_many`; see
         !! `davidson_result_t`.
      logical, intent(out) :: converged
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      integer, intent(in), optional :: max_subspace
         !! Vectors kept before the subspace is collapsed onto the current Ritz
         !! vectors. Defaults to `max(2*n_roots + 8, 16)`, bounded by the size of
         !! the determinant space.
      real(dp), intent(in), optional :: guess(:, :)
         !! (n_determinants, n_start) starting vectors, `n_start >= n_roots`.
         !! Absent takes `n_roots` unit vectors on the lowest diagonal
         !! elements, which for a CI is the reference determinant and its
         !! nearest neighbours in energy.
         !!
         !! More columns than roots widens the **starting subspace** and does
         !! not ask for more roots: `n_roots` eigenpairs still come back. A
         !! root the `n_roots` lowest diagonal elements cannot reach -- one
         !! the off-diagonal pushes below others whose diagonal sits lower --
         !! is otherwise missed silently, the solver converging cleanly onto
         !! whatever its subspace does span. Columns past `max_subspace` are
         !! dropped.
      logical, intent(in), optional :: verbose
         !! Print a line per iteration, off by default
      real(dp), intent(in), optional :: energy_offset
         !! Added to the eigenvalue before printing, and to nothing else. The
         !! Davidson solves the active-space problem, so its own eigenvalue is
         !! the active energy alone; the caller knows the inactive-plus-nuclear
         !! constant this adds back.
      character(len=*), intent(in), optional :: label
         !! What the iteration table is called. Defaults to the CI this solver
         !! was written for; the electronic Hessian is the other operator that
         !! reaches it, and a verbose run of that has no business claiming to
         !! be a CI.
      character(len=*), intent(in), optional :: value_label
         !! What the eigenvalue column is called, `energy` by default. The same
         !! reason: the lowest eigenvalue of an orbital-rotation Hessian is a
         !! curvature, not an energy.

      real(dp), allocatable :: basis(:, :), sigma(:, :), small(:, :), small_values(:)
      character(len=128) :: line, header, name_of_value
      integer(int64) :: tick, last, rate
      real(dp) :: shift
      logical :: loud
      real(dp), allocatable :: ritz(:, :), hritz(:, :), residual(:), correction(:)
      real(dp) :: tol, norm, denominator, overlap
      integer :: ndet, nsub, nmax, nstart, iterations, iroot, i, j, info, added
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

      ! A guess wider than the roots asked for is a wider starting subspace,
      ! bounded by what the subspace can hold.
      nstart = n_roots
      if (present(guess)) nstart = min(max(size(guess, 2), n_roots), nmax)

      allocate (basis(ndet, nmax), sigma(ndet, nmax))
      allocate (ritz(ndet, n_roots), hritz(ndet, n_roots))
      allocate (residual(ndet), correction(ndet))
      allocate (root_converged(n_roots))
      allocate (values(n_roots), residuals(n_roots))
      allocate (vectors(ndet, n_roots))

      loud = .false.
      if (present(verbose)) loud = verbose
      shift = 0.0_dp
      if (present(energy_offset)) shift = energy_offset
      name_of_value = "energy"
      if (present(value_label)) name_of_value = value_label
      ! Built rather than a literal, so the column heading lines up with the
      ! `f23.12` the rows are written in whatever the caller named it.
      write (header, "(a,a23,a)") "    iter", trim(name_of_value), &
         "    residual   subspace     sigma       time"
      if (present(label)) then
         call convergence_header(loud, label, trim(header), 74)
      else
         call convergence_header(loud, "CI iterations", trim(header), 74)
      end if

      call system_clock(last, rate)
      call initial_basis(diagonal, nstart, ndet, basis, error, guess)
      if (error%has_error()) return
      nsub = nstart

      ! Sigma for the starting vectors, as one block.
      call operator%apply_many(basis(:, 1:nsub), sigma(:, 1:nsub), error)
      if (error%has_error()) return
      sigma_products = sigma_products + nsub

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

         ! Time per iteration rather than a total, as the SCF prints it.
         if (loud) then
            call system_clock(tick)
            write (line, "(i8,f23.12,es12.3,i11,i10,f10.2,a)") &
               iteration, values(1) + shift, maxval(residuals), nsub, sigma_products, &
               real(tick - last, dp)/real(rate, dp), " s"
            call logger%info(trim(line))
            ! Redirected output is block buffered, so without this the rows sit
            ! in a 4 kB buffer and arrive in lumps.
            flush (output_unit)
            last = tick
         end if

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

         ! The whole expansion in one call. Each correction was orthogonalised
         ! against the ones before it as it was built, so the block is already
         ! independent and nothing is lost by deferring the products to here.
         call operator%apply_many(basis(:, nsub + 1:nsub + added), &
                                  sigma(:, nsub + 1:nsub + added), error)
         if (error%has_error()) return
         sigma_products = sigma_products + added

         nsub = nsub + added
      end do

      call convergence_footer(loud, converged, iterations_taken, "iterations", 74)
      if (.not. converged) vectors = ritz

      deallocate (basis, sigma, ritz, hritz, residual, correction)
      deallocate (root_converged)
   end subroutine davidson_flat

   subroutine initial_basis(diagonal, n_start, ndet, basis, error, guess)
      !! Starting vectors: the supplied ones, or the lowest determinants
      !!
      !! The linear-dependence test here is on the *fraction* of a guess
      !! column that survives projection, for the reason `LINEAR_DEPENDENCE`
      !! gives: a test against an absolute length is a test of how long the
      !! caller happened to make its vectors. Comparing the surviving length
      !! against the column's own length is scale-free without rescaling the
      !! column, which would round every well-conditioned guess for the sake
      !! of the degenerate ones.
      !!
      !! Without a guess, `lowest_free` cannot come back empty: this branch
      !! takes `n_start` determinants out of `ndet`, `davidson_flat` refuses
      !! `n_roots > ndet` before calling here, and the widening it does after
      !! that is bounded by a subspace already capped at `ndet`.
      real(dp), intent(in) :: diagonal(:)
      integer, intent(in) :: n_start
         !! Columns of `basis` to fill. More than the roots asked for is a
         !! wider starting subspace; a guess narrower than this is topped up
         !! from the lowest free diagonal elements, as an absent one is.
      integer, intent(in) :: ndet
      real(dp), intent(inout) :: basis(:, :)
      type(error_t), intent(inout) :: error
         !! Set when a guess column cannot be replaced by any determinant
      real(dp), intent(in), optional :: guess(:, :)

      logical, allocatable :: taken(:)
      real(dp) :: norm, length
      integer :: iroot, pick

      allocate (taken(ndet))
      taken = .false.

      if (present(guess)) then
         do iroot = 1, n_start
            if (iroot <= size(guess, 2)) then
               basis(:, iroot) = guess(:, iroot)
            else
               ! Past the end of the guess. A zero column has zero length, so
               ! the fallback below takes it as dependent and fills it from the
               ! determinants, which is what an absent guess would have done.
               basis(:, iroot) = 0.0_dp
            end if
            length = sqrt(dot_product(basis(:, iroot), basis(:, iroot)))
            norm = project_out_earlier(basis, iroot)
            ! A supplied vector the earlier ones already span leaves nothing
            ! behind to normalise, and a near-null column makes every Ritz
            ! vector built on it meaningless. Fall back to the lowest
            ! determinant still free, which is where an absent guess would
            ! have started, and keep falling back until something survives.
            do while (norm <= LINEAR_DEPENDENCE*length)
               pick = lowest_free(diagonal, taken)
               if (pick == 0) then
                  ! Every determinant has been spent and the column is still
                  ! null. Keeping it would seed the subspace with noise and
                  ! report the resulting Ritz values as eigenvalues, so the
                  ! solve stops instead.
                  call error%set(ERROR_VALIDATION, "starting vector "//to_char(iroot)// &
                                 " is spanned by the ones before it, and all "// &
                                 to_char(ndet)//" determinants have already been "// &
                                 "used to replace one.")
                  deallocate (taken)
                  return
               end if
               taken(pick) = .true.
               basis(:, iroot) = 0.0_dp
               basis(pick, iroot) = 1.0_dp
               length = 1.0_dp
               norm = project_out_earlier(basis, iroot)
            end do
            basis(:, iroot) = basis(:, iroot)/norm
         end do
         deallocate (taken)
         return
      end if

      basis(:, 1:n_start) = 0.0_dp
      do iroot = 1, n_start
         ! Never zero here; see the note on the routine.
         pick = lowest_free(diagonal, taken)
         taken(pick) = .true.
         basis(pick, iroot) = 1.0_dp
      end do
      deallocate (taken)
   end subroutine initial_basis

   function lowest_free(diagonal, taken) result(pick)
      !! The determinant of lowest diagonal energy not already spoken for
      real(dp), intent(in) :: diagonal(:)
      logical, intent(in) :: taken(:)
      integer :: pick
         !! Zero when every determinant has been taken

      real(dp) :: best
      integer :: i

      pick = 0
      best = huge(1.0_dp)
      do i = 1, size(diagonal)
         if (taken(i)) cycle
         if (diagonal(i) < best) then
            best = diagonal(i)
            pick = i
         end if
      end do
   end function lowest_free

   function project_out_earlier(basis, iroot) result(norm)
      !! Project columns 1 to `iroot - 1` out of column `iroot`, in place
      !!
      !! Twice, for the same reason the expansion loop does it twice: one
      !! pass of Gram-Schmidt loses orthogonality when the column is nearly
      !! in the subspace already, and a column nearly in the subspace is the
      !! entire input domain of the fallback this feeds. A single pass leaves
      !! a residue of the earlier columns behind, which is then divided by a
      !! small norm and blown up into the basis.
      real(dp), intent(inout) :: basis(:, :)
      integer, intent(in) :: iroot
      real(dp) :: norm
         !! What is left of the column's length, before it is normalised

      real(dp) :: overlap
      integer :: j, pass

      do pass = 1, 2
         do j = 1, iroot - 1
            overlap = dot_product(basis(:, j), basis(:, iroot))
            basis(:, iroot) = basis(:, iroot) - overlap*basis(:, j)
         end do
      end do
      norm = sqrt(dot_product(basis(:, iroot), basis(:, iroot)))
   end function project_out_earlier

end module mqc_davidson
