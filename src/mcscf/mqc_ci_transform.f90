!! Carrying a CI vector from one orbital basis to another
module mqc_ci_transform
   !! A complete active space is invariant under rotation of its active
   !! orbitals: the same wave function, expanded in a different set of
   !! determinants. This module performs that re-expansion, so that a CI which
   !! has been solved in one basis need not be solved again in the other.
   !!
   !! The identity is the one every textbook derives and nobody writes down.
   !! With orbitals related by an orthogonal `T`,
   !!
   !!     |q> = sum_p |p> T(p,q)
   !!
   !! a string of creation operators factorizes into a sum over strings of the
   !! other basis with a **minor determinant** as its coefficient:
   !!
   !!     M(I,J) = det( T[ occ(I), occ(J) ] )
   !!
   !! -- rows the orbitals occupied in `I`, columns those occupied in `J`, both
   !! in ascending order. Alpha and beta factorize independently, so the whole
   !! transformation of a rectangular CI vector is two matrix products:
   !!
   !!     C' = M_alpha^T C M_beta
   !!
   !! **`M` is exactly the size of the CI vector, and the transformation costs
   !! one dimension more.** For a CAS(14,14) the vector is 3432 x 3432 and the
   !! products are 3432^3, which is dense BLAS and a few seconds. Compare that
   !! with re-solving: a Davidson in a basis with no dominant determinant starts
   !! from a guess with negligible overlap and preconditions on a diagonal that
   !! no longer approximates the spectrum, and may take hundreds of iterations
   !! or none that converge.
   !!
   !! **`M` inherits `T`'s orthogonality**, which is the check worth making and
   !! `rotation_is_orthogonal` makes it. If `T` is orthogonal then so is `M`,
   !! for every electron count at once, and a violation means the rotation
   !! handed in was not what it claimed to be rather than that this arithmetic
   !! is wrong.
   !!
   !! The asymptotically better route is Malmqvist's: factorize `T` into
   !! elementary one-orbital operations and apply each straight to the vector,
   !! `O(n_str^2 n_orb^2)` and no `M` stored. It is considerably more code, and
   !! the active spaces a full CI can be solved in at all are small enough that
   !! the matrix form is not the constraint.
   use, intrinsic :: iso_fortran_env, only: int64
   use pic_types, only: dp
   use pic_io, only: to_char
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_determinants, only: generate_strings, n_strings, string_address
   use mqc_ormas_space, only: ormas_space_t, ormas_strings, determinant_address
   implicit none
   private

   public :: string_rotation_matrix
   public :: transform_ci_vector
   public :: scatter_restricted_ci
   public :: rotation_is_orthogonal
   public :: orthogonality_defect

   real(dp), parameter :: TOL_ORTHOGONAL = 1.0e-10_dp
      !! How far `T^T T` may sit from the identity, per element. The
      !! transformations this is fed are built from a symmetric
      !! orthogonalization followed by Jacobi rotations, so they are orthogonal
      !! to rounding and nothing else is expected to reach here.
   integer(int64), parameter :: MAX_RECTANGLE = 2_int64**28
      !! Determinants, about 2.1 GB of coefficients. The ceiling on
      !! `scatter_restricted_ci`, and it is a real one: a restricted space is
      !! chosen precisely because the complete space is unaffordable, and this
      !! routine spends the complete space anyway. Ethane's full valence
      !! CAS(14,14) is 11,778,624 determinants -- ninety-four megabytes, and
      !! twenty times inside this -- so the molecules the analysis can reach at
      !! all are reachable. Benzene's is 2.4e16 and is not, by this route or any
      !! other.
   real(dp), parameter :: TOL_NORM = 1.0e-8_dp
      !! How far the transformed vector's norm may drift from the original's.
      !! An orthogonal transformation preserves it exactly, so this tolerates
      !! the accumulated rounding of two matrix products and no more.

contains

   pure function rotation_is_orthogonal(rotation) result(orthogonal)
      !! Whether `T^T T` is the identity to `TOL_ORTHOGONAL`
      !!
      !! Cheap -- `n_orb^3` on a matrix whose dimension is the active space --
      !! and the precondition for everything below, since a non-orthogonal
      !! rotation does not relate two orthonormal determinant bases and the
      !! minors would not be a transformation of anything.
      real(dp), intent(in) :: rotation(:, :)
      logical :: orthogonal

      orthogonal = orthogonality_defect(rotation) <= TOL_ORTHOGONAL
   end function rotation_is_orthogonal

   pure function orthogonality_defect(rotation) result(defect)
      !! The largest element of `T^T T - I`, or `huge` if `T` is not square
      !!
      !! Reported rather than only tested because the size of the miss is the
      !! diagnosis. Two orthonormal bases of the *same* space give rounding;
      !! two bases of spaces that merely have the same dimension give something
      !! far larger, and a caller deciding whether to trust an offered wave
      !! function wants to be able to say which it is looking at.
      real(dp), intent(in) :: rotation(:, :)
      real(dp) :: defect

      real(dp) :: element
      integer :: n, i, j, k

      n = size(rotation, 1)
      if (size(rotation, 2) /= n) then
         defect = huge(1.0_dp)
         return
      end if
      defect = 0.0_dp
      do j = 1, n
         do i = 1, n
            element = 0.0_dp
            do k = 1, n
               element = element + rotation(k, i)*rotation(k, j)
            end do
            if (i == j) element = element - 1.0_dp
            defect = max(defect, abs(element))
         end do
      end do
   end function orthogonality_defect

   subroutine string_rotation_matrix(rotation, n_electrons, matrix, error)
      !! The determinant-basis image of an orbital rotation, for one spin
      !!
      !! `matrix(I,J)` is the amplitude of target-basis string `J` in
      !! source-basis string `I`, so a source-basis vector is carried over by
      !! `matrix^T`. Strings are indexed as `generate_strings` orders them,
      !! which is ascending bit pattern, and the orbitals within a string are
      !! taken in ascending order -- the same convention the creation operators
      !! are written in throughout `mqc_determinants`, which is what makes the
      !! minor's sign the determinant's sign with no correction.
      real(dp), intent(in) :: rotation(:, :)
         !! (n_orbitals, n_orbitals) orthogonal. Column `q` expands target
         !! orbital `q` in the source orbitals, which is the direction
         !! `|q> = sum_p |p> T(p,q)` is written in.
      integer, intent(in) :: n_electrons
         !! Of one spin. The two spins of a CI vector need one call each unless
         !! their counts agree, in which case one matrix serves both.
      real(dp), allocatable, intent(out) :: matrix(:, :)
         !! (n_strings, n_strings)
      type(error_t), intent(inout) :: error

      integer(int64), allocatable :: strings(:)
      integer, allocatable :: occupancy(:, :)
      real(dp), allocatable :: rows(:, :), work(:, :)
      integer :: n_orb, n_str, i, j, k, p, filled

      if (error%has_error()) return

      n_orb = size(rotation, 1)
      if (size(rotation, 2) /= n_orb) then
         call error%set(ERROR_VALIDATION, "the orbital rotation is "// &
                        to_char(n_orb)//" by "//to_char(size(rotation, 2))// &
                        ", and a rotation between two bases of the same space "// &
                        "is square.")
         return
      end if
      if (.not. rotation_is_orthogonal(rotation)) then
         call error%set(ERROR_VALIDATION, "the orbital rotation is not orthogonal, so "// &
                        "it does not relate two orthonormal sets of determinants and "// &
                        "the transformation would not preserve the wave function.")
         return
      end if

      call generate_strings(n_orb, n_electrons, strings, error)
      if (error%has_error()) return
      n_str = size(strings)
      allocate (matrix(n_str, n_str))

      if (n_electrons == 0) then
         ! One string, holding nothing, and the empty determinant is the same
         ! object in either basis. The minor of an empty submatrix is one.
         matrix = 1.0_dp
         deallocate (strings)
         return
      end if

      ! Which orbitals each string occupies, ascending, once rather than once
      ! per pair: the pair count is the square of the string count, so scanning
      ! bits inside the double loop would repeat every scan `n_str` times.
      allocate (occupancy(n_electrons, n_str))
      do i = 1, n_str
         filled = 0
         do p = 1, n_orb
            if (btest(strings(i), p - 1)) then
               filled = filled + 1
               occupancy(filled, i) = p
            end if
         end do
      end do

      ! `rows` is hoisted out of the inner loop: for a fixed source string the
      ! rows of the submatrix are the same for every target string, and only
      ! the columns gathered from them change.
      !$omp parallel default(shared) private(i, j, k, rows, work)
      allocate (rows(n_electrons, n_orb), work(n_electrons, n_electrons))
      !$omp do schedule(static)
      do i = 1, n_str
         do k = 1, n_electrons
            rows(k, :) = rotation(occupancy(k, i), :)
         end do
         do j = 1, n_str
            do k = 1, n_electrons
               work(:, k) = rows(:, occupancy(k, j))
            end do
            matrix(i, j) = determinant(work)
         end do
      end do
      !$omp end do
      deallocate (rows, work)
      !$omp end parallel

      deallocate (strings, occupancy)
   end subroutine string_rotation_matrix

   pure function determinant(a) result(value)
      !! Gaussian elimination with partial pivoting, on a copy
      !!
      !! Small -- the dimension is the electron count of one spin -- so LAPACK
      !! would spend more in its own call overhead than this spends computing,
      !! and it is called once per pair of strings.
      real(dp), intent(in) :: a(:, :)
      real(dp) :: value

      real(dp) :: matrix(size(a, 1), size(a, 1))
      real(dp) :: row(size(a, 1))
      real(dp) :: factor
      integer :: n, i, j, pivot

      n = size(a, 1)
      matrix = a
      value = 1.0_dp
      do i = 1, n
         pivot = i
         do j = i + 1, n
            if (abs(matrix(j, i)) > abs(matrix(pivot, i))) pivot = j
         end do
         if (matrix(pivot, i) == 0.0_dp) then
            ! A singular submatrix is not an error: it means these two strings
            ! do not overlap, which most pairs do not.
            value = 0.0_dp
            return
         end if
         if (pivot /= i) then
            row = matrix(i, :)
            matrix(i, :) = matrix(pivot, :)
            matrix(pivot, :) = row
            value = -value
         end if
         value = value*matrix(i, i)
         do j = i + 1, n
            factor = matrix(j, i)/matrix(i, i)
            matrix(j, i + 1:n) = matrix(j, i + 1:n) - factor*matrix(i, i + 1:n)
         end do
      end do
   end function determinant

   subroutine scatter_restricted_ci(space, flat, ci, error)
      !! Lay a restricted-space CI vector out over the complete space
      !!
      !! An occupation-restricted space (ORMAS) keeps a chosen subset of the
      !! determinants, so its coefficients arrive as a flat list rather than an
      !! alpha-by-beta rectangle. This writes them into the rectangle, zero
      !! everywhere the restriction excluded, which is the same wave function
      !! written in the complete basis.
      !!
      !! **This exists so that a restricted wave function can be rotated at
      !! all.** A restricted space is not closed under rotation of its active
      !! orbitals -- mixing two active orbitals stops being redundant the moment
      !! the space is incomplete -- so a rotation carries such a wave function
      !! out of its own space and there is nowhere in that space to put the
      !! answer. The complete space is where the answer lives, for the same
      !! reason that a subspace rotated inside a larger space lands on a
      !! different subspace of it.
      !!
      !! GAMESS solves the same problem the other way, by *forbidding* the
      !! rotation to leave the space: `LOCAL_PPASVD` builds the quasi-atomic
      !! transformation block-diagonal over the subspaces and
      !! `LOCAL_BLOCK_ORMAS` asserts it six times over. The expansion done here
      !! is named there as `MALMQ_CIDRIV` option 2, "expands ORMAS CI vector
      !! from NSPACE>1 to NSPACE=1 form", and declined -- "THIS OPTION IS NOT
      !! YET IMPLEMENTED AND MIGHT NOT BE". That is a reasonable decision for
      !! Malmqvist's algorithm, which acts on the vector in its own space and
      !! would lose its whole advantage by expanding; it is not one for a
      !! transformation that is already dense matrix products over the complete
      !! rectangle.
      !!
      !! The cost is honest and is the reason for `MAX_RECTANGLE`: the
      !! restriction buys a cheap *solve* and this spends the complete space to
      !! store the result. For ethane that is 16,500 determinants converged and
      !! 11.7 million stored, which is the right trade; past a point it is not
      !! available at all.
      type(ormas_space_t), intent(in) :: space
      real(dp), intent(in) :: flat(:)          !! (space%n_determinants)
      real(dp), allocatable, intent(out) :: ci(:, :)
      type(error_t), intent(inout) :: error

      integer(int64), allocatable :: alpha(:), beta(:)
      integer(int64) :: n_full_alpha, n_full_beta, address
      integer :: n_orb, ia, ib, fa, fb, ga, gb
      real(dp) :: before, after

      if (error%has_error()) return

      n_orb = space%n_active_orbitals
      n_full_alpha = n_strings(n_orb, space%n_alpha)
      n_full_beta = n_strings(n_orb, space%n_beta)
      if (n_full_alpha*n_full_beta > MAX_RECTANGLE) then
         call error%set(ERROR_VALIDATION, "the restricted space holds "// &
                        to_char(space%n_determinants)//" determinants, but the "// &
                        "complete space it would have to be written out over holds "// &
                        to_char(n_full_alpha*n_full_beta)//", which is past what this "// &
                        "can allocate. The restriction makes the wave function "// &
                        "affordable to solve for and not to rotate.")
         return
      end if
      if (int(size(flat), int64) /= space%n_determinants) then
         call error%set(ERROR_VALIDATION, "the coefficient list holds "// &
                        to_char(size(flat))//" values and the restricted space has "// &
                        to_char(space%n_determinants)//" determinants.")
         return
      end if

      call ormas_strings(space, alpha, beta, error)
      if (error%has_error()) return

      allocate (ci(n_full_alpha, n_full_beta))
      ci = 0.0_dp
      do ia = 1, size(alpha)
         ga = space%alpha_string_class(ia)
         fa = string_address(n_orb, space%n_alpha, alpha(ia))
         do ib = 1, size(beta)
            gb = space%beta_string_class(ib)
            ! `determinant_address` returns a plausible number for an
            ! incompatible pair rather than a signal, so the gate is here.
            if (.not. space%compatible(gb, ga)) cycle
            fb = string_address(n_orb, space%n_beta, beta(ib))
            address = determinant_address(space, ia, ib)
            ci(fa, fb) = flat(address)
         end do
      end do

      ! Every coefficient placed exactly once, nothing dropped and nothing
      ! written twice. A collision in the addressing would show up here and
      ! nowhere else downstream.
      before = sum(flat**2)
      after = sum(ci**2)
      if (abs(after - before) > TOL_NORM*max(1.0_dp, before)) then
         call error%set(ERROR_VALIDATION, "the restricted coefficients did not "// &
                        "survive being written out over the complete space, so two "// &
                        "determinants were given the same address or one was missed.")
         return
      end if

      deallocate (alpha, beta)
   end subroutine scatter_restricted_ci

   subroutine transform_ci_vector(m_alpha, m_beta, ci, transformed, error)
      !! Re-expand a CI vector in the target basis: `C' = M_alpha^T C M_beta`
      !!
      !! The amplitudes change and the wave function does not. What comes back
      !! is the same state, so its energy against the target-basis Hamiltonian
      !! must equal the energy it had in the basis it was solved in -- which is
      !! the check a caller should make, because it tests the whole chain
      !! rather than this arithmetic alone.
      real(dp), intent(in) :: m_alpha(:, :)   !! From `string_rotation_matrix`, alpha
      real(dp), intent(in) :: m_beta(:, :)    !! The same for beta
      real(dp), intent(in) :: ci(:, :)        !! (n_alpha_strings, n_beta_strings)
      real(dp), allocatable, intent(out) :: transformed(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: work(:, :)
      real(dp) :: before, after
      integer :: na, nb

      if (error%has_error()) return
      na = size(ci, 1)
      nb = size(ci, 2)
      if (size(m_alpha, 1) /= na .or. size(m_alpha, 2) /= na .or. &
          size(m_beta, 1) /= nb .or. size(m_beta, 2) /= nb) then
         call error%set(ERROR_VALIDATION, "the string transformations are "// &
                        to_char(size(m_alpha, 1))//" and "//to_char(size(m_beta, 1))// &
                        " but the coefficient array is "//to_char(na)//" by "// &
                        to_char(nb)//", so they were built for a different active space.")
         return
      end if

      allocate (work(na, nb), transformed(na, nb))
      call pic_gemm(m_alpha, ci, work, transa="T")
      call pic_gemm(work, m_beta, transformed)

      ! Orthogonal in, orthogonal out. This does not prove the convention is
      ! right -- a sign error survives it -- but it does catch a transformation
      ! built for the wrong space, which is the mistake with no other symptom.
      before = sum(ci**2)
      after = sum(transformed**2)
      if (abs(after - before) > TOL_NORM*max(1.0_dp, before)) then
         call error%set(ERROR_VALIDATION, "the transformed CI vector has a different "// &
                        "norm from the one it came from, which an orthogonal "// &
                        "transformation cannot do.")
         return
      end if

      deallocate (work)
   end subroutine transform_ci_vector

end module mqc_ci_transform
