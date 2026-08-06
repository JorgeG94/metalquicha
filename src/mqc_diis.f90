!! DIIS subspace acceleration for SCF
module mqc_diis
   !! Holds the DIIS history and produces the extrapolated Fock matrix.
   !!
   !! Split out of the cuEST SCF driver for two reasons. It is pure linear
   !! algebra over flat vectors with no integrals backend behind it, so it builds
   !! and is unit-tested without a GPU -- which the SCF driver is not, since that
   !! needs libcuest and an sm_80 card. And it gives the device port a host
   !! reference to diff against: the failure mode of a device DIIS is a silently
   !! stale history producing plausible-but-wrong coefficients, which only shows
   !! up against a known-good implementation.
   !!
   !! Vectors are stored flat, one column per subspace entry, which serves the
   !! restricted case (one n_ao*n_ao Fock) and the unrestricted case (both spins
   !! stacked into one vector) through the same code, and is the layout a device
   !! reduction wants: contiguous and coalesced.
   use pic_types, only: dp
   implicit none
   private

   public :: diis_state_t
   public :: diis_slot_of_age   !! Ring slot holding the age-th oldest entry
   public :: diis_coefficients  !! Extrapolation weights from a cached overlap matrix

   !> The two procedures above are the parts of DIIS that have nothing to do
   !  with where the vectors live: which slot an entry occupies, and what
   !  weights a given overlap matrix implies. A device implementation shares
   !  them rather than reimplementing them, so a host-versus-device divergence
   !  can only have come from the vectors themselves -- which is the whole
   !  point of keeping a host reference.

   real(dp), parameter :: PIVOT_FLOOR = 1.0e-14_dp
      !! Below this a pivot is treated as singular and extrapolation is skipped

   type :: diis_state_t
      !! One SCF's worth of DIIS history
      integer :: max_vectors = 0   !! Subspace size
      integer :: n_fock = 0        !! Elements in one Fock vector
      integer :: n_error = 0       !! Elements in one error vector
      integer :: n_stored = 0      !! Entries currently held (<= max_vectors)
      integer :: newest = 0        !! Slot holding the most recent entry

      real(dp), allocatable :: fock_history(:, :)   !! (n_fock, max_vectors)
      real(dp), allocatable :: error_history(:, :)  !! (n_error, max_vectors)

      real(dp), allocatable :: overlap(:, :)
         !! (max_vectors, max_vectors) cached <e_i|e_j> in slot coordinates.
         !! Kept across iterations: pushing one vector changes exactly one row
         !! and column, so rebuilding the whole matrix each iteration costs
         !! O(n_stored^2 * n_error) where O(n_stored * n_error) suffices.
   contains
      procedure :: init => diis_init
      procedure :: push => diis_push
      procedure :: extrapolate => diis_extrapolate
      procedure :: destroy => diis_destroy
      procedure :: count => diis_count
   end type diis_state_t

contains

   subroutine diis_init(this, max_vectors, n_fock, n_error)
      !! Allocate the history for a subspace of max_vectors entries
      class(diis_state_t), intent(inout) :: this
      integer, intent(in) :: max_vectors  !! Subspace size
      integer, intent(in) :: n_fock       !! Elements in one Fock vector
      integer, intent(in) :: n_error      !! Elements in one error vector

      call this%destroy()

      this%max_vectors = max(1, max_vectors)
      this%n_fock = n_fock
      this%n_error = n_error
      this%n_stored = 0
      this%newest = 0

      allocate (this%fock_history(n_fock, this%max_vectors))
      allocate (this%error_history(n_error, this%max_vectors))
      allocate (this%overlap(this%max_vectors, this%max_vectors))
      this%overlap = 0.0_dp
   end subroutine diis_init

   pure function diis_count(this) result(n)
      !! Number of entries currently in the subspace
      class(diis_state_t), intent(in) :: this
      integer :: n
      n = this%n_stored
   end function diis_count

   pure function diis_slot_of_age(newest, n_stored, max_vectors, age) result(slot)
      !! Slot holding the age-th oldest entry, age = 1 .. n_stored
      !!
      !! The history is a ring: pushing past the end overwrites the oldest entry
      !! in place, so nothing is ever shifted. The previous implementation moved
      !! the whole history down by one every iteration once full, which is
      !! O(max_vectors * n_fock) of pure copying per SCF step.
      integer, intent(in) :: newest, n_stored, max_vectors
      integer, intent(in) :: age
      integer :: slot
      slot = modulo(newest - n_stored + age - 1, max_vectors) + 1
   end function diis_slot_of_age

   pure function slot_of_age(this, age) result(slot)
      !! `diis_slot_of_age` for this state's ring
      class(diis_state_t), intent(in) :: this
      integer, intent(in) :: age
      integer :: slot
      slot = diis_slot_of_age(this%newest, this%n_stored, this%max_vectors, age)
   end function slot_of_age

   subroutine diis_push(this, fock, error_vector)
      !! Add one (Fock, error) pair, evicting the oldest entry when full
      class(diis_state_t), intent(inout) :: this
      real(dp), intent(in) :: fock(this%n_fock)
      real(dp), intent(in) :: error_vector(this%n_error)

      integer :: slot, age, other

      this%newest = modulo(this%newest, this%max_vectors) + 1
      if (this%n_stored < this%max_vectors) this%n_stored = this%n_stored + 1
      slot = this%newest

      this%fock_history(:, slot) = fock
      this%error_history(:, slot) = error_vector

      ! Only the new entry's row and column of the overlap matrix have changed.
      do age = 1, this%n_stored
         other = slot_of_age(this, age)
         this%overlap(slot, other) = sum(this%error_history(:, slot)* &
                                         this%error_history(:, other))
         this%overlap(other, slot) = this%overlap(slot, other)
      end do
   end subroutine diis_push

   subroutine diis_extrapolate(this, fock, ok)
      !! Replace fock with the DIIS combination of the stored Fock matrices
      !!
      !! fock is left untouched when the subspace is too small or the linear
      !! system is singular, so a caller can use the result unconditionally.
      class(diis_state_t), intent(in) :: this
      real(dp), intent(inout) :: fock(this%n_fock)
      logical, intent(out) :: ok

      real(dp), allocatable :: coefficients(:)
      integer :: i

      call diis_coefficients(this%overlap, this%newest, this%n_stored, &
                             this%max_vectors, coefficients, ok)
      if (ok) then
         fock = 0.0_dp
         do i = 1, this%n_stored
            fock = fock + coefficients(i)*this%fock_history(:, slot_of_age(this, i))
         end do
      end if

      if (allocated(coefficients)) deallocate (coefficients)
   end subroutine diis_extrapolate

   subroutine diis_coefficients(overlap, newest, n_stored, max_vectors, coefficients, ok)
      !! Extrapolation weights, oldest-first, from a cached error overlap matrix
      !!
      !! `coefficients(i)` weights the age-i-oldest entry, so a caller walks its
      !! history through `diis_slot_of_age` in the same order.
      !!
      !! The B matrix is assembled oldest-first rather than in slot order. The
      !! DIIS solution is invariant to permutation of the subspace, but the
      !! elimination is not bit-invariant, and age order is what the pre-ring
      !! implementation solved -- so results stay reproducible across that
      !! change, and across the host/device split.
      real(dp), intent(in) :: overlap(:, :)   !! (max_vectors, max_vectors), slot coordinates
      integer, intent(in) :: newest, n_stored, max_vectors
      real(dp), allocatable, intent(out) :: coefficients(:)
         !! (n_stored+1), oldest first. The trailing entry is the Lagrange
         !! multiplier enforcing sum(c) = 1 and is not a weight -- read
         !! 1..n_stored only.
      logical, intent(out) :: ok              !! False when the subspace is too small or singular

      real(dp), allocatable :: b_matrix(:, :)
      integer :: n, i, j

      ok = .false.
      if (n_stored < 2) return

      n = n_stored
      allocate (b_matrix(n + 1, n + 1))
      b_matrix = -1.0_dp
      b_matrix(n + 1, n + 1) = 0.0_dp
      do j = 1, n
         do i = 1, n
            b_matrix(i, j) = overlap(diis_slot_of_age(newest, n_stored, max_vectors, i), &
                                     diis_slot_of_age(newest, n_stored, max_vectors, j))
         end do
      end do

      call solve_diis(b_matrix, coefficients, ok)
      deallocate (b_matrix)
   end subroutine diis_coefficients

   subroutine diis_destroy(this)
      !! Release the history
      class(diis_state_t), intent(inout) :: this
      if (allocated(this%fock_history)) deallocate (this%fock_history)
      if (allocated(this%error_history)) deallocate (this%error_history)
      if (allocated(this%overlap)) deallocate (this%overlap)
      this%max_vectors = 0
      this%n_fock = 0
      this%n_error = 0
      this%n_stored = 0
      this%newest = 0
   end subroutine diis_destroy

   subroutine solve_diis(b_matrix, coefficients, ok)
      !! Solve the small DIIS linear system by Gaussian elimination
      !!
      !! The system is at most (diis_size+1) square, so a dedicated LAPACK
      !! solve would cost more in dependencies than it saves in time.
      real(dp), intent(in) :: b_matrix(:, :)
      real(dp), allocatable, intent(out) :: coefficients(:)
      logical, intent(out) :: ok

      real(dp), allocatable :: augmented(:, :)
      integer :: n, i, j, pivot_row
      real(dp) :: pivot, factor

      n = size(b_matrix, 1)
      allocate (augmented(n, n + 1), coefficients(n))
      augmented(:, 1:n) = b_matrix
      augmented(:, n + 1) = 0.0_dp
      augmented(n, n + 1) = -1.0_dp
      ok = .true.

      do i = 1, n
         pivot_row = i
         do j = i + 1, n
            if (abs(augmented(j, i)) > abs(augmented(pivot_row, i))) pivot_row = j
         end do
         if (pivot_row /= i) then
            augmented([i, pivot_row], :) = augmented([pivot_row, i], :)
         end if

         pivot = augmented(i, i)
         if (abs(pivot) < PIVOT_FLOOR) then
            ok = .false.
            return
         end if

         do j = i + 1, n
            factor = augmented(j, i)/pivot
            augmented(j, i:n + 1) = augmented(j, i:n + 1) - factor*augmented(i, i:n + 1)
         end do
      end do

      do i = n, 1, -1
         coefficients(i) = (augmented(i, n + 1) - &
                            sum(augmented(i, i + 1:n)*coefficients(i + 1:n)))/augmented(i, i)
      end do
   end subroutine solve_diis

end module mqc_diis
