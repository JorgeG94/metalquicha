!! Stand-in for the HDF5 checkpoint when the build has no HDF5
module mqc_hdf5_checkpoint
   !! Same name, same type, same procedures as the real backend, and every one
   !! of them declines.
   !!
   !! A stub rather than a preprocessor conditional at each call site, matching
   !! `mqc_cuest_bridge_stub` next door: the code that uses a checkpoint reads
   !! the same either way, and there is exactly one place where "this build
   !! cannot do that" is decided. `hdf5_checkpoint_available` is how a caller
   !! asks in advance, so the refusal can name the build option instead of
   !! failing somewhere less useful.
   use pic_types, only: dp, int64
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: hdf5_checkpoint_t
   public :: hdf5_checkpoint_available

   type :: hdf5_checkpoint_t
      logical :: active = .false.
      integer(int64) :: n_loaded = 0
   contains
      procedure :: open => h5ck_open
      procedure :: record => h5ck_record
      procedure :: lookup => h5ck_lookup
      procedure :: close => h5ck_close
   end type hdf5_checkpoint_t

contains

   pure function hdf5_checkpoint_available() result(available)
      !! Whether this build can write HDF5 checkpoints
      logical :: available
      available = .false.
   end function hdf5_checkpoint_available

   subroutine h5ck_open(this, path, fingerprint, max_level, error)
      class(hdf5_checkpoint_t), intent(inout) :: this
      character(len=*), intent(in) :: path
      character(len=*), intent(in) :: fingerprint
      integer, intent(in) :: max_level
      type(error_t), intent(inout) :: error

      this%active = .false.
      call error%set(ERROR_VALIDATION, "this build has no HDF5, so "//trim(path)// &
                     " cannot be written; configure with -DMQC_ENABLE_HDF5=ON. "// &
                     "Energy-only runs need no HDF5 and use the text format.")
      if (len_trim(fingerprint) == 0 .or. max_level < 0) return
   end subroutine h5ck_open

   subroutine h5ck_record(this, term, energy, scf_status, n_atoms, gradient, hessian)
      class(hdf5_checkpoint_t), intent(inout) :: this
      integer, intent(in) :: term(:)
      real(dp), intent(in) :: energy
      integer, intent(in) :: scf_status
      integer, intent(in) :: n_atoms
      real(dp), intent(in), optional :: gradient(:, :)
      real(dp), intent(in), optional :: hessian(:, :)

      if (this%active .or. size(term) < 0 .or. energy /= energy .or. &
          scf_status < 0 .or. n_atoms < 0) return
      if (present(gradient) .or. present(hessian)) return
   end subroutine h5ck_record

   subroutine h5ck_lookup(this, term, found, energy, scf_status, n_atoms, gradient, hessian)
      class(hdf5_checkpoint_t), intent(in) :: this
      integer, intent(in) :: term(:)
      logical, intent(out) :: found
      real(dp), intent(out) :: energy
      integer, intent(out) :: scf_status
      integer, intent(out) :: n_atoms
      real(dp), allocatable, intent(out), optional :: gradient(:, :)
      real(dp), allocatable, intent(out), optional :: hessian(:, :)

      found = .false.
      energy = 0.0_dp
      scf_status = 0
      n_atoms = 0
      if (this%active .and. size(term) > 0) return
      if (present(gradient) .or. present(hessian)) return
   end subroutine h5ck_lookup

   subroutine h5ck_close(this)
      class(hdf5_checkpoint_t), intent(inout) :: this
      this%active = .false.
   end subroutine h5ck_close

end module mqc_hdf5_checkpoint
