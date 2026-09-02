!! DIIS subspace acceleration for SCF
module mqc_diis
   !! Holds the DIIS history and produces the extrapolated Fock matrix.
   !!
   !! Pure linear algebra over flat vectors, with no integrals backend behind
   !! it. Vectors are stored flat, one column per subspace entry, which serves
   !! the restricted case (one n_ao*n_ao Fock) and the unrestricted case (both
   !! spins stacked into one vector) through the same code.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_ediis, only: ediis_coefficients, adiis_coefficients
   implicit none
   private

   public :: diis_state_t
   public :: ACCEL_DIIS, ACCEL_EDIIS, ACCEL_ADIIS
   public :: accelerator_name
   public :: parse_accelerator_name
   public :: diis_slot_of_age   !! Ring slot holding the age-th oldest entry
   public :: diis_coefficients  !! Extrapolation weights from a cached overlap matrix

   ! The two procedures above are the parts of DIIS that have nothing to do
   ! with where the vectors live, so a device implementation shares them rather
   ! than reimplementing them.

   real(dp), parameter :: PIVOT_FLOOR = 1.0e-14_dp
      !! Below this a pivot is treated as singular and extrapolation is skipped

   real(dp), parameter, public :: ACCEL_SWITCH = 1.0e-2_dp
      !! Error-vector norm below which an energy-based scheme hands over to
      !! DIIS. Applied by the SCF driver, not by this module.
      !!
      !! **A deliberate deviation from PySCF**, which has no handover: setting
      !! `mf.DIIS` to `ADIIS` or `EDIIS` there runs that class from the first
      !! iteration to the last. Neither energy scheme is usable alone, and the
      !! reason is structural: the convex constraint that keeps them inside the
      !! hull of the densities seen so far is what stops them reaching a
      !! solution lying outside it.
      !!
      !! The handover is a *performance* choice; nothing rests on it for
      !! correctness.

   integer, parameter :: ACCEL_DIIS = 1   !! Pulay, on the error vectors
   integer, parameter :: ACCEL_EDIIS = 2  !! Kudin/Scuseria/Cances, on the energy
   integer, parameter :: ACCEL_ADIIS = 3  !! Hu/Yang, augmented Roothaan-Hall

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

      logical :: energy_based = .false.
         !! Whether the history below is carried. Off by default and allocated
         !! only on request: a density is the size of a Fock matrix, so this
         !! doubles the history's footprint.

      real(dp), allocatable :: density_history(:, :)  !! (n_fock, max_vectors)
      real(dp), allocatable :: energy_history(:)      !! (max_vectors)

      real(dp), allocatable :: df(:, :)
         !! (max_vectors, max_vectors) cached <D_i | F_j> in slot coordinates,
         !! which is everything EDIIS and ADIIS need from the history.
         !!
         !! **Not symmetric.** `<D_i|F_j>` and `<D_j|F_i>` differ off the
         !! diagonal and both energy models depend on the difference, so a new
         !! entry updates a full row *and* a full column rather than one
         !! triangle -- unlike `overlap` above, which is mirrored.
   contains
      procedure :: init => diis_init
      procedure :: push => diis_push
      procedure :: extrapolate => diis_extrapolate
      procedure :: extrapolate_with => diis_extrapolate_with
      procedure :: destroy => diis_destroy
      procedure :: count => diis_count
   end type diis_state_t

contains

   subroutine diis_init(this, max_vectors, n_fock, n_error, energy_based)
      !! Allocate the history for a subspace of max_vectors entries
      class(diis_state_t), intent(inout) :: this
      integer, intent(in) :: max_vectors  !! Subspace size
      integer, intent(in) :: n_fock       !! Elements in one Fock vector
      integer, intent(in) :: n_error      !! Elements in one error vector
      logical, intent(in), optional :: energy_based
         !! Carry densities and energies too, for EDIIS and ADIIS. Absent is
         !! false, so a plain DIIS caller allocates exactly what it did before.

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

      this%energy_based = .false.
      if (present(energy_based)) this%energy_based = energy_based
      if (this%energy_based) then
         allocate (this%density_history(n_fock, this%max_vectors))
         allocate (this%energy_history(this%max_vectors))
         allocate (this%df(this%max_vectors, this%max_vectors))
         this%energy_history = 0.0_dp
         this%df = 0.0_dp
      end if
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
      !! in place, so nothing is ever shifted.
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

   subroutine diis_push(this, fock, error_vector, density, energy, error)
      !! Add one entry, evicting the oldest when full
      class(diis_state_t), intent(inout) :: this
      real(dp), intent(in) :: fock(this%n_fock)
      real(dp), intent(in) :: error_vector(this%n_error)
      real(dp), intent(in), optional :: density(this%n_fock)
         !! Required when the state was initialised `energy_based`; ignored
         !! otherwise. Same layout as `fock` -- for an unrestricted SCF both
         !! spins stacked, so `<D|F>` comes out as the sum over spins, which is
         !! the trace both energy models want.
      real(dp), intent(in), optional :: energy
         !! The total energy of this entry, for EDIIS. ADIIS does not use it.
      type(error_t), intent(inout), optional :: error
         !! Set when an `energy_based` state is pushed without a `density`.

      integer :: slot, age, other

      ! Rejected before anything is written, not half-applied. The history is a
      ! ring that overwrites in place, so a push that skipped the density would
      ! leave the *evicted* entry's density in the slot -- and the `df` update
      ! below immediately contracts that slot against every stored Fock, so the
      ! stale density reaches both EDIIS and ADIIS as if it belonged there.
      if (this%energy_based .and. .not. present(density)) then
         if (present(error)) then
            call error%set(ERROR_VALIDATION, "DIIS: an energy-based subspace needs a "// &
                           "density with every entry, and this push has none")
         end if
         return
      end if

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

      if (this%energy_based) then
         this%density_history(:, slot) = density
         if (present(energy)) this%energy_history(slot) = energy

         ! A full row *and* a full column, because `df` is not symmetric: the
         ! new density against every stored Fock, and every stored density
         ! against the new Fock. Mirroring one onto the other, as `overlap`
         ! above correctly does, would quietly symmetrise both energy models.
         do age = 1, this%n_stored
            other = slot_of_age(this, age)
            this%df(slot, other) = sum(this%density_history(:, slot)* &
                                       this%fock_history(:, other))
            this%df(other, slot) = sum(this%density_history(:, other)* &
                                       this%fock_history(:, slot))
         end do
      end if
   end subroutine diis_push

   subroutine diis_extrapolate_with(this, scheme, fock, ok)
      !! Replace fock with the combination the named scheme asks for
      !!
      !! DIIS extrapolates on the error vectors and can leave the hull of the
      !! densities it was built from; EDIIS and ADIIS interpolate inside it.
      !! That is the whole difference: the first is fast near the solution and
      !! unreliable far from it, the second pair the other way round.
      !!
      !! One call runs one scheme. Choosing which, and handing over to DIIS
      !! below `ACCEL_SWITCH`, is the caller's job -- see `mqc_libcint_rhf`.
      class(diis_state_t), intent(in) :: this
      integer, intent(in) :: scheme          !! One of the ACCEL_* parameters
      real(dp), intent(inout) :: fock(this%n_fock)
      logical, intent(out) :: ok             !! False leaves fock untouched

      real(dp), allocatable :: coefficients(:), df_age(:, :), e_age(:)
      integer :: i, j, n

      ok = .false.
      if (scheme == ACCEL_DIIS) then
         call diis_extrapolate(this, fock, ok)
         return
      end if
      if (.not. this%energy_based) return

      ! Both energy models work oldest-first, matching what `diis_coefficients`
      ! returns, so the caches are transposed out of slot coordinates once here
      ! rather than indexed through `slot_of_age` inside two minimisations.
      n = this%n_stored
      allocate (df_age(n, n), e_age(n))
      do j = 1, n
         do i = 1, n
            df_age(i, j) = this%df(slot_of_age(this, i), slot_of_age(this, j))
         end do
         e_age(j) = this%energy_history(slot_of_age(this, j))
      end do

      select case (scheme)
      case (ACCEL_EDIIS)
         call ediis_coefficients(e_age, df_age, n, coefficients, ok)
      case (ACCEL_ADIIS)
         call adiis_coefficients(df_age, n, coefficients, ok)
      case default
         ! ACCEL_DIIS returned above, so reaching here is an unknown scheme.
         ! Leaving `ok` false means the caller keeps its own Fock matrix, which
         ! is the one safe answer -- inventing coefficients for a scheme nobody
         ! named would be a silently different SCF.
         ok = .false.
      end select

      if (ok) then
         fock = 0.0_dp
         do i = 1, n
            fock = fock + coefficients(i)*this%fock_history(:, slot_of_age(this, i))
         end do
      end if

      if (allocated(coefficients)) deallocate (coefficients)
      deallocate (df_age, e_age)
   end subroutine diis_extrapolate_with

   subroutine parse_accelerator_name(name, scheme, ok)
      !! `keywords.scf.accelerator` to one of the ACCEL_* parameters
      !!
      !! Unknown names are refused rather than defaulted: `scheme` still comes
      !! back `ACCEL_DIIS`, but `ok` is false so the caller can say so.
      character(len=*), intent(in) :: name
      integer, intent(out) :: scheme
      logical, intent(out) :: ok

      character(len=len(name)) :: lower
      integer :: i, c

      lower = name
      do i = 1, len_trim(lower)
         c = iachar(lower(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) lower(i:i) = achar(c + 32)
      end do

      ok = .true.
      select case (trim(adjustl(lower)))
      case ("", "diis")
         scheme = ACCEL_DIIS
      case ("ediis")
         scheme = ACCEL_EDIIS
      case ("adiis")
         scheme = ACCEL_ADIIS
      case default
         scheme = ACCEL_DIIS
         ok = .false.
      end select
   end subroutine parse_accelerator_name

   pure function accelerator_name(scheme) result(name)
      !! For the log line that says which one ran
      integer, intent(in) :: scheme
      character(len=5) :: name
      select case (scheme)
      case (ACCEL_EDIIS)
         name = "EDIIS"
      case (ACCEL_ADIIS)
         name = "ADIIS"
      case default
         name = "DIIS "
      end select
   end function accelerator_name

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
      !! DIIS solution is invariant to permutation of the subspace but the
      !! elimination is not bit-invariant, so fixing the order keeps results
      !! reproducible across the host/device split.
      real(dp), intent(in) :: overlap(:, :)   !! (max_vectors, max_vectors), slot coordinates
      integer, intent(in) :: newest, n_stored, max_vectors
      real(dp), allocatable, intent(out) :: coefficients(:)
         !! (n_stored+1), oldest first. The trailing entry is the Lagrange
         !! multiplier enforcing sum(c) = 1 and is not a weight -- read
         !! 1..n_stored only.
      logical, intent(out) :: ok              !! False when the subspace is too small or singular

      real(dp), allocatable :: b_matrix(:, :)
      real(dp) :: scale
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

      ! Scale the error block to O(1) against the -1 constraint border. Its
      ! entries are ~|e|^2, around 1e-16 by the time the SCF should stop, and
      ! eliminating that unscaled makes a saturated subspace's weights noise --
      ! without a visible blow-up: the density drifts away from the fixed point
      ! with the energy flat. The weights are invariant under the scaling, which
      ! rescales the Lagrange multiplier alone.
      scale = maxval(abs(b_matrix(1:n, 1:n)))
      if (scale > 0.0_dp) b_matrix(1:n, 1:n) = b_matrix(1:n, 1:n)/scale

      call solve_diis(b_matrix, coefficients, ok)
      deallocate (b_matrix)
   end subroutine diis_coefficients

   subroutine diis_destroy(this)
      !! Release the history
      class(diis_state_t), intent(inout) :: this
      if (allocated(this%fock_history)) deallocate (this%fock_history)
      if (allocated(this%error_history)) deallocate (this%error_history)
      if (allocated(this%overlap)) deallocate (this%overlap)
      if (allocated(this%density_history)) deallocate (this%density_history)
      if (allocated(this%energy_history)) deallocate (this%energy_history)
      if (allocated(this%df)) deallocate (this%df)
      this%energy_based = .false.
      this%max_vectors = 0
      this%n_fock = 0
      this%n_error = 0
      this%n_stored = 0
      this%newest = 0
   end subroutine diis_destroy

   subroutine solve_diis(b_matrix, coefficients, ok)
      !! Solve the small DIIS linear system by Gaussian elimination
      !!
      !! The system is at most (diis_size+1) square.
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
