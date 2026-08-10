!! DIIS subspace acceleration over device-resident vectors
module mqc_diis_device
   !! The device twin of `mqc_diis`: same ring, same overlap cache, same
   !! coefficients, but the Fock and error histories live on the GPU and are
   !! never copied to the host.
   !!
   !! What is shared with the host version, deliberately. `diis_slot_of_age`
   !! and `diis_coefficients` come straight out of `mqc_diis` rather than being
   !! reimplemented here, so the ring indexing and the small linear solve are
   !! literally the same code. That narrows a host-versus-device divergence to
   !! the vector arithmetic, which is the only part that genuinely differs --
   !! and the characteristic DIIS failure (a stale history producing plausible
   !! coefficients, no crash, no error) is otherwise very hard to localise.
   !!
   !! What is not shared: the overlap matrix stays on the host in both. It is
   !! max_vectors square -- 64 doubles at the usual subspace size -- and the
   !! solve over it is a Gaussian elimination that would need cuSOLVER to move,
   !! for no measurable gain.
   !!
   !! Two cuBLAS habits this module exists to enforce:
   !!
   !!   * The overlap row is **one `cublasDgemv`** of the new error vector
   !!     against the whole history, not `n_stored` calls to `cublasDdot`.
   !!     `Ddot` returns through a host pointer and therefore blocks; the naive
   !!     loop costs one synchronise per stored vector, every iteration.
   !!   * History slots are offsets into one allocation, not separate buffers.
   !!     A `cudaMalloc` per slot per fragment would dominate a fragmented run.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_int64_t
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_diis, only: diis_slot_of_age, diis_coefficients
   use mqc_cuest_context, only: cuest_context_t
   use mqc_cuest_runtime, only: cublas_status_check, copy_to_host, device_offset
   use cublas, only: cublasDaxpy, cublasDcopy, cublasDscal, cublasDgemv, CUBLAS_OP_T
   implicit none
   private

   public :: diis_device_t

   type :: diis_device_t
      !! One SCF's worth of DIIS history, held on the device
      integer :: max_vectors = 0   !! Subspace size
      integer :: n_fock = 0        !! Elements in one Fock vector
      integer :: n_error = 0       !! Elements in one error vector
      integer :: n_stored = 0      !! Entries currently held (<= max_vectors)
      integer :: newest = 0        !! Slot holding the most recent entry

      type(c_ptr) :: cublas = c_null_ptr  !! Borrowed handle, not owned here

      ! Borrowed from the rank's pools, exactly as `cuest_system_t` borrows its
      ! scratch: `destroy` drops the references rather than freeing, so the next
      ! fragment reuses the allocation.
      type(c_ptr) :: d_fock_history = c_null_ptr   !! (n_fock, max_vectors)
      type(c_ptr) :: d_error_history = c_null_ptr  !! (n_error, max_vectors)
      type(c_ptr) :: d_row = c_null_ptr            !! (max_vectors), gemv output

      real(dp), allocatable :: overlap(:, :)
         !! (max_vectors, max_vectors) cached <e_i|e_j> in slot coordinates,
         !! on the HOST. Kept across iterations for the same reason as in
         !! `mqc_diis`: a push changes one row and column, not the whole matrix.
      real(dp), allocatable :: row_host(:)
         !! Staging for the one row the gemv produces, max_vectors doubles
   contains
      procedure :: init => diis_device_init
      procedure :: push => diis_device_push
      procedure :: extrapolate => diis_device_extrapolate
      procedure :: destroy => diis_device_destroy
      procedure :: count => diis_device_count
   end type diis_device_t

contains

   subroutine diis_device_init(this, context, max_vectors, n_fock, n_error, error)
      !! Size the history for a subspace of max_vectors entries
      class(diis_device_t), intent(inout) :: this
      type(cuest_context_t), intent(inout) :: context  !! Per-rank handle and pools
      integer, intent(in) :: max_vectors  !! Subspace size
      integer, intent(in) :: n_fock       !! Elements in one Fock vector
      integer, intent(in) :: n_error      !! Elements in one error vector
      type(error_t), intent(inout) :: error

      call this%destroy()

      this%max_vectors = max(1, max_vectors)
      this%n_fock = n_fock
      this%n_error = n_error
      this%n_stored = 0
      this%newest = 0
      this%cublas = context%cublas_handle

      call context%scratch_diis_fock%ensure(int(n_fock, c_int64_t)* &
                                            int(this%max_vectors, c_int64_t), &
                                            "DIIS Fock history", error)
      call context%scratch_diis_error%ensure(int(n_error, c_int64_t)* &
                                             int(this%max_vectors, c_int64_t), &
                                             "DIIS error history", error)
      call context%scratch_diis_row%ensure(int(this%max_vectors, c_int64_t), &
                                           "DIIS overlap row", error)
      if (error%has_error()) then
         call error%add_context("mqc_diis_device:init")
         return
      end if

      this%d_fock_history = context%scratch_diis_fock%ptr
      this%d_error_history = context%scratch_diis_error%ptr
      this%d_row = context%scratch_diis_row%ptr

      allocate (this%overlap(this%max_vectors, this%max_vectors))
      allocate (this%row_host(this%max_vectors))
      this%overlap = 0.0_dp
      this%row_host = 0.0_dp
   end subroutine diis_device_init

   pure function diis_device_count(this) result(n)
      !! Number of entries currently in the subspace
      class(diis_device_t), intent(in) :: this
      integer :: n
      n = this%n_stored
   end function diis_device_count

   subroutine diis_device_push(this, d_fock, d_error, error)
      !! Add one (Fock, error) pair, evicting the oldest entry when full
      !!
      !! Both arguments are device addresses and stay on the device: the copies
      !! into the ring are device-to-device. The only thing that crosses back is
      !! the new row of the overlap matrix, `n_stored` doubles.
      class(diis_device_t), intent(inout) :: this
      type(c_ptr), intent(in) :: d_fock    !! (n_fock) on device
      type(c_ptr), intent(in) :: d_error   !! (n_error) on device
      type(error_t), intent(inout) :: error

      integer :: slot, other

      if (error%has_error()) return

      this%newest = modulo(this%newest, this%max_vectors) + 1
      if (this%n_stored < this%max_vectors) this%n_stored = this%n_stored + 1
      slot = this%newest

      call cublas_status_check(cublasDcopy(this%cublas, int(this%n_fock, c_int), &
                                           d_fock, 1, &
                                           device_offset(this%d_fock_history, &
                                                         int(slot - 1, c_int64_t)* &
                                                         int(this%n_fock, c_int64_t)), 1), &
                               "cublasDcopy(Fock -> DIIS history)", error)
      call cublas_status_check(cublasDcopy(this%cublas, int(this%n_error, c_int), &
                                           d_error, 1, &
                                           device_offset(this%d_error_history, &
                                                         int(slot - 1, c_int64_t)* &
                                                         int(this%n_error, c_int64_t)), 1), &
                               "cublasDcopy(error -> DIIS history)", error)
      if (error%has_error()) return

      ! One gemv gives every <e_new|e_i> at once. Columns 1..n_stored are
      ! exactly the occupied slots: the ring fills 1,2,3,... in order, so until
      ! it wraps the used slots are the leading ones, and once it wraps every
      ! slot is used. Nothing past n_stored is read, so unused columns never
      ! need clearing.
      call cublas_status_check(cublasDgemv(this%cublas, CUBLAS_OP_T, &
                                           int(this%n_error, c_int), &
                                           int(this%n_stored, c_int), &
                                           1.0_dp, this%d_error_history, &
                                           int(this%n_error, c_int), &
                                           d_error, 1, 0.0_dp, this%d_row, 1), &
                               "cublasDgemv(DIIS overlap row)", error)
      if (error%has_error()) return

      ! A synchronous D2H on the default stream is ordered after the gemv, so
      ! no explicit synchronise is needed -- and adding one here would put a
      ! device-wide stall in the middle of every SCF iteration.
      call copy_to_host(this%row_host(1:this%n_stored), this%d_row, &
                        "DIIS overlap row", error)
      if (error%has_error()) return

      ! Only the new entry's row and column have changed.
      do other = 1, this%n_stored
         this%overlap(slot, other) = this%row_host(other)
         this%overlap(other, slot) = this%row_host(other)
      end do
   end subroutine diis_device_push

   subroutine diis_device_extrapolate(this, d_fock, ok, error)
      !! Replace the device Fock with the DIIS combination of the history
      !!
      !! `d_fock` is left untouched when the subspace is too small or the linear
      !! system is singular, so a caller can use the result unconditionally.
      class(diis_device_t), intent(inout) :: this
      type(c_ptr), intent(in) :: d_fock   !! (n_fock) on device, updated in place
      logical, intent(out) :: ok
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: coefficients(:)
      integer :: i, slot

      ok = .false.
      if (error%has_error()) return

      call diis_coefficients(this%overlap, this%newest, this%n_stored, &
                             this%max_vectors, coefficients, ok)
      if (.not. ok) return

      ! Accumulated oldest-first, matching `mqc_diis`, so the two agree to the
      ! last bit rather than merely to within rounding.
      call cublas_status_check(cublasDscal(this%cublas, int(this%n_fock, c_int), &
                                           0.0_dp, d_fock, 1), &
                               "cublasDscal(clear Fock)", error)
      do i = 1, this%n_stored
         slot = diis_slot_of_age(this%newest, this%n_stored, this%max_vectors, i)
         call cublas_status_check(cublasDaxpy(this%cublas, int(this%n_fock, c_int), &
                                              coefficients(i), &
                                              device_offset(this%d_fock_history, &
                                                            int(slot - 1, c_int64_t)* &
                                                            int(this%n_fock, c_int64_t)), 1, &
                                              d_fock, 1), &
                                  "cublasDaxpy(DIIS extrapolation)", error)
      end do

      if (allocated(coefficients)) deallocate (coefficients)
      if (error%has_error()) ok = .false.
   end subroutine diis_device_extrapolate

   subroutine diis_device_destroy(this)
      !! Drop the borrowed history and release the host-side cache
      class(diis_device_t), intent(inout) :: this

      ! Borrowed from the context's pools -- drop the references, do not free.
      this%d_fock_history = c_null_ptr
      this%d_error_history = c_null_ptr
      this%d_row = c_null_ptr
      this%cublas = c_null_ptr

      if (allocated(this%overlap)) deallocate (this%overlap)
      if (allocated(this%row_host)) deallocate (this%row_host)
      this%max_vectors = 0
      this%n_fock = 0
      this%n_error = 0
      this%n_stored = 0
      this%newest = 0
   end subroutine diis_device_destroy

end module mqc_diis_device
