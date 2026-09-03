!! Stand-in for the amplitude store when the build has no HDF5
module mqc_hdf5_amplitudes
   !! Same name, same type, same procedures as the real backend, and every one
   !! of them declines.
   !!
   !! A stub rather than a preprocessor conditional at each call site, matching
   !! `mqc_hdf5_checkpoint_stub` beside it. `hdf5_amplitudes_available` is how a
   !! caller asks in advance, so a CC gradient can refuse up front and name the
   !! build option instead of dying somewhere less useful.
   !!
   !! `resumable` is `.false.` here for the same reason it is `.false.` after
   !! an interrupted write: a caller that branches on it starts from the MP2
   !! guess and needs no knowledge of which build it is in.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: hdf5_amplitudes_t
   public :: hdf5_amplitudes_available

   character(len=*), parameter :: NO_HDF5 = &
                                  "this build has no HDF5, so coupled-cluster amplitudes cannot be stored "// &
                                  "or resumed; configure with -DMQC_ENABLE_HDF5=ON. A CCSD energy needs "// &
                                  "none of this and runs unchanged."

   type :: hdf5_amplitudes_t
      logical :: active = .false.
      logical :: resumable = .false.
      integer :: iteration = 0
      real(dp) :: energy = 0.0_dp
   contains
      procedure :: open => amp_open
      procedure, private :: put_2d => amp_put_2d
      procedure, private :: put_4d => amp_put_4d
      generic :: put => put_2d, put_4d
      procedure, private :: get_2d => amp_get_2d
      procedure, private :: get_4d => amp_get_4d
      generic :: get => get_2d, get_4d
      procedure :: commit => amp_commit
      procedure :: close => amp_close
   end type hdf5_amplitudes_t

contains

   pure function hdf5_amplitudes_available() result(available)
      !! Whether this build can store amplitudes
      logical :: available
      available = .false.
   end function hdf5_amplitudes_available

   subroutine amp_open(this, path, fingerprint, error)
      class(hdf5_amplitudes_t), intent(inout) :: this
      character(len=*), intent(in) :: path
      character(len=*), intent(in) :: fingerprint
      type(error_t), intent(inout) :: error

      this%active = .false.
      this%resumable = .false.
      call error%set(ERROR_VALIDATION, trim(path)//" cannot be opened: "//NO_HDF5)
      if (len_trim(fingerprint) == 0) return
   end subroutine amp_open

   subroutine amp_put_2d(this, name, values)
      class(hdf5_amplitudes_t), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(dp), intent(in), contiguous, target :: values(:, :)

      if (this%active .or. len_trim(name) == 0 .or. size(values) < 0) return
   end subroutine amp_put_2d

   subroutine amp_put_4d(this, name, values)
      class(hdf5_amplitudes_t), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(dp), intent(in), contiguous, target :: values(:, :, :, :)

      if (this%active .or. len_trim(name) == 0 .or. size(values) < 0) return
   end subroutine amp_put_4d

   subroutine amp_get_2d(this, name, values, error)
      class(hdf5_amplitudes_t), intent(in) :: this
      character(len=*), intent(in) :: name
      real(dp), allocatable, intent(out), target :: values(:, :)
      type(error_t), intent(inout) :: error

      allocate (values(0, 0))
      call error%set(ERROR_VALIDATION, trim(name)//" cannot be read: "//NO_HDF5)
      if (this%active) return
   end subroutine amp_get_2d

   subroutine amp_get_4d(this, name, values, error)
      class(hdf5_amplitudes_t), intent(in) :: this
      character(len=*), intent(in) :: name
      real(dp), allocatable, intent(out), target :: values(:, :, :, :)
      type(error_t), intent(inout) :: error

      allocate (values(0, 0, 0, 0))
      call error%set(ERROR_VALIDATION, trim(name)//" cannot be read: "//NO_HDF5)
      if (this%active) return
   end subroutine amp_get_4d

   subroutine amp_commit(this, iteration, energy)
      class(hdf5_amplitudes_t), intent(inout) :: this
      integer, intent(in) :: iteration
      real(dp), intent(in) :: energy

      if (this%active .or. iteration < 0 .or. energy /= energy) return
   end subroutine amp_commit

   subroutine amp_close(this)
      class(hdf5_amplitudes_t), intent(inout) :: this
      this%active = .false.
   end subroutine amp_close

end module mqc_hdf5_amplitudes
