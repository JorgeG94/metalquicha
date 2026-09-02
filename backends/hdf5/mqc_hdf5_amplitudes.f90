!! HDF5 store for coupled-cluster amplitudes, so a killed run can resume
module mqc_hdf5_amplitudes
   !! A small number of large fixed-shape arrays, overwritten in place: `t1`,
   !! `t2`, `l1` and `l2` for one calculation, rewritten every macro-iteration.
   !! `mqc_hdf5_checkpoint` next door is append-only and keyed by a monomer
   !! tuple, which is why the two are separate modules.
   !!
   !! **The T amplitudes are stored, not only Lambda.** The Lambda equations
   !! are built from `t1` and `t2` at every iteration, so a file holding `l1`
   !! and `l2` alone cannot restart anything.
   !!
   !! **A kill must not be able to produce a plausible file.** `complete` is an
   !! attribute that starts at zero; every array is written, the file is
   !! flushed, and only then does `complete` become one. A reader that finds
   !! zero treats the file as holding nothing, so a kill loses a whole
   !! macro-iteration rather than part of an array.
   !!
   !! **The fingerprint is a refusal, not a hint.** Amplitudes belong to a
   !! geometry, a basis, a frozen-core choice and an orbital count. A mismatch
   !! stops the run rather than silently starting fresh, because a set from
   !! another calculation converges to a stationary point of the wrong
   !! equations.
   use, intrinsic :: iso_c_binding, only: c_loc, c_char, c_size_t
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_IO, ERROR_VALIDATION
   use mqc_hdf5_bindings, only: hid_t, hsize_t, herr_t, h5_start, c_string, &
                                H5F_ACC_TRUNC, H5F_ACC_RDWR, H5P_DEFAULT, H5S_ALL, &
                                H5F_SCOPE_LOCAL, H5T_NATIVE_DOUBLE, H5T_NATIVE_INT, &
                                H5T_STD_I32LE, H5T_C_S1, &
                                H5Fcreate, H5Fopen, H5Fclose, H5Fflush, &
                                H5Screate, H5Screate_simple, H5Sclose, &
                                H5Sget_simple_extent_dims, H5Sget_simple_extent_ndims, &
                                H5Dcreate2, H5Dopen2, H5Dclose, H5Dwrite, H5Dread, &
                                H5Dget_space, H5Lexists, &
                                H5Acreate2, H5Aopen, H5Awrite, H5Aread, H5Aclose, H5Aexists, &
                                H5Tcopy, H5Tset_size, H5Tclose
   implicit none
   private

   public :: hdf5_amplitudes_t
   public :: hdf5_amplitudes_available

   integer, parameter :: FP_LEN = 17   !! 16 hex digits and a terminator

   type :: hdf5_amplitudes_t
      !! An open amplitude file, and what a resumed run needs to know about it
      logical :: active = .false.
      logical :: resumable = .false.
         !! A complete set was found and the fingerprint matched. False on a
         !! new file and false on one whose last write was interrupted; a
         !! caller reads this to decide between resuming and starting from the
         !! MP2 guess, and never has to ask why.
      integer :: iteration = 0
      real(dp) :: energy = 0.0_dp
      integer(hid_t) :: file = -1
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
      available = .true.
   end function hdf5_amplitudes_available

   subroutine amp_open(this, path, fingerprint, error)
      !! Open `path`, creating it if absent and validating it if not
      class(hdf5_amplitudes_t), intent(inout) :: this
      character(len=*), intent(in) :: path
      character(len=*), intent(in) :: fingerprint
      type(error_t), intent(inout) :: error

      logical :: exists
      character(len=FP_LEN) :: stored
      integer :: complete

      this%active = .false.
      this%resumable = .false.
      this%iteration = 0
      this%energy = 0.0_dp

      if (.not. h5_start()) then
         call error%set(ERROR_IO, "the HDF5 library would not initialise, so "// &
                        trim(path)//" cannot be used for amplitudes")
         return
      end if

      inquire (file=path, exist=exists)
      if (exists) then
         this%file = H5Fopen(c_string(path), H5F_ACC_RDWR, H5P_DEFAULT)
         if (this%file < 0) then
            call error%set(ERROR_IO, trim(path)//" exists but is not a readable HDF5 file")
            return
         end if

         call read_string_attr(this%file, "fingerprint", stored)
         if (trim(stored) /= trim(fingerprint)) then
            call error%set(ERROR_VALIDATION, "amplitude file "//trim(path)//" was written "// &
                           "by a different calculation (fingerprint "//trim(stored)// &
                           ", this run is "//trim(fingerprint)//"). Its amplitudes solve "// &
                           "another geometry, basis or frozen-core choice; resuming from "// &
                           "them would converge to a stationary point of the wrong "// &
                           "equations. Delete it, or point the run at another file.")
            return
         end if

         ! Zero means the last run was killed while writing. Whatever is in
         ! the datasets may be a complete iteration or half of one, and there
         ! is no way to tell from here -- so it is worth nothing either way.
         call read_int_attr(this%file, "complete", complete)
         if (complete == 1) then
            this%resumable = .true.
            call read_int_attr(this%file, "iteration", this%iteration)
            call read_double_attr(this%file, "energy", this%energy)
         end if
      else
         this%file = H5Fcreate(c_string(path), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)
         if (this%file < 0) then
            call error%set(ERROR_IO, "could not create amplitude file "//trim(path))
            return
         end if
         call write_string_attr(this%file, "fingerprint", fingerprint)
         call write_int_attr(this%file, "complete", 0)
      end if

      this%active = .true.
   end subroutine amp_open

   subroutine amp_put_2d(this, name, values)
      !! Store a rank-2 amplitude block, overwriting any previous one
      class(hdf5_amplitudes_t), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(dp), intent(in), contiguous, target :: values(:, :)

      integer(hsize_t) :: dims(2)

      if (.not. this%active) return
      dims = int(shape(values), hsize_t)
      call put_block(this, name, 2, dims, c_loc(values))
   end subroutine amp_put_2d

   subroutine amp_put_4d(this, name, values)
      !! Store a rank-4 amplitude block, overwriting any previous one
      class(hdf5_amplitudes_t), intent(inout) :: this
      character(len=*), intent(in) :: name
      real(dp), intent(in), contiguous, target :: values(:, :, :, :)

      integer(hsize_t) :: dims(4)

      if (.not. this%active) return
      dims = int(shape(values), hsize_t)
      call put_block(this, name, 4, dims, c_loc(values))
   end subroutine amp_put_4d

   subroutine put_block(this, name, rank, dims, buffer)
      !! The shared body: create on first write, reopen and overwrite after
      !!
      !! `complete` drops to zero here rather than in `commit`: from the first
      !! byte of the first array until the flush in `commit` lands, the file
      !! advertises itself as unusable.
      use, intrinsic :: iso_c_binding, only: c_ptr
      type(hdf5_amplitudes_t), intent(inout) :: this
      character(len=*), intent(in) :: name
      integer, intent(in) :: rank
      integer(hsize_t), intent(in) :: dims(:)
      type(c_ptr), intent(in) :: buffer

      integer(hid_t) :: space, dset
      integer(hsize_t) :: file_dims(size(dims))
      integer(herr_t) :: ignored

      call write_int_attr(this%file, "complete", 0)

      ! The extents are reversed on the way in, and this is not cosmetic.
      ! HDF5 declares a dataset in C order -- first extent slowest-varying --
      ! while Fortran's leftmost index is the fastest. Written straight
      ! through, `t1(3,5)` becomes a 3x5 dataset whose rows interleave the two
      ! indices: the round trip through this module is exact and every other
      ! reader gets nonsense. Reversed, the dataset is 5x3 and
      ! `t1_h5[a, i] == t1(i, a)` in h5py, h5dump and anything else. HDF5's
      ! own Fortran API reverses here too.
      file_dims = dims(rank:1:-1)

      if (H5Lexists(this%file, c_string(name), H5P_DEFAULT) > 0) then
         dset = H5Dopen2(this%file, c_string(name), H5P_DEFAULT)
      else
         space = H5Screate_simple(rank, file_dims, file_dims)
         dset = H5Dcreate2(this%file, c_string(name), H5T_NATIVE_DOUBLE, space, &
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)
         ignored = H5Sclose(space)
      end if
      if (dset < 0) return

      ignored = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer)
      ignored = H5Dclose(dset)
   end subroutine put_block

   subroutine amp_get_2d(this, name, values, error)
      !! Read a rank-2 block back, allocated to the shape the file records
      class(hdf5_amplitudes_t), intent(in) :: this
      character(len=*), intent(in) :: name
      real(dp), allocatable, intent(out), target :: values(:, :)
      type(error_t), intent(inout) :: error

      integer(hsize_t) :: dims(2)
      integer(hid_t) :: dset

      call open_block(this, name, 2, dims, dset, error)
      if (error%has_error() .or. dset < 0) return
      allocate (values(dims(1), dims(2)))
      call read_block(dset, c_loc(values))
   end subroutine amp_get_2d

   subroutine amp_get_4d(this, name, values, error)
      !! Read a rank-4 block back, allocated to the shape the file records
      class(hdf5_amplitudes_t), intent(in) :: this
      character(len=*), intent(in) :: name
      real(dp), allocatable, intent(out), target :: values(:, :, :, :)
      type(error_t), intent(inout) :: error

      integer(hsize_t) :: dims(4)
      integer(hid_t) :: dset

      call open_block(this, name, 4, dims, dset, error)
      if (error%has_error() .or. dset < 0) return
      allocate (values(dims(1), dims(2), dims(3), dims(4)))
      call read_block(dset, c_loc(values))
   end subroutine amp_get_4d

   subroutine open_block(this, name, want_rank, dims, dset, error)
      !! Locate a block and report its shape, refusing one of the wrong rank
      !!
      !! Blocks are written by name, so a caller that asks for `l1` and is
      !! handed `l2` would otherwise get an allocation of the right total size
      !! and the wrong meaning.
      class(hdf5_amplitudes_t), intent(in) :: this
      character(len=*), intent(in) :: name
      integer, intent(in) :: want_rank
      integer(hsize_t), intent(out) :: dims(:)
      integer(hid_t), intent(out) :: dset
      type(error_t), intent(inout) :: error

      integer(hid_t) :: space
      integer(hsize_t) :: file_dims(size(dims)), maxdims(size(dims))
      integer :: rank
      integer(herr_t) :: ignored

      dset = -1
      dims = 0

      if (.not. this%active) then
         call error%set(ERROR_IO, "no amplitude file is open, so "//trim(name)// &
                        " cannot be read")
         return
      end if

      if (H5Lexists(this%file, c_string(name), H5P_DEFAULT) <= 0) then
         call error%set(ERROR_VALIDATION, "the amplitude file holds no block named "// &
                        trim(name))
         return
      end if

      dset = H5Dopen2(this%file, c_string(name), H5P_DEFAULT)
      if (dset < 0) then
         call error%set(ERROR_IO, "could not open the amplitude block "//trim(name))
         return
      end if

      space = H5Dget_space(dset)
      rank = int(H5Sget_simple_extent_ndims(space))
      if (rank /= want_rank) then
         ignored = H5Sclose(space)
         ignored = H5Dclose(dset)
         dset = -1
         call error%set(ERROR_VALIDATION, "the amplitude block "//trim(name)//" is stored "// &
                        "with a rank this reader did not ask for; the file was written by "// &
                        "a different version of this program.")
         return
      end if
      ! Reversed back, to undo the C-order declaration made on write.
      ignored = H5Sget_simple_extent_dims(space, file_dims, maxdims)
      dims = file_dims(rank:1:-1)
      ignored = H5Sclose(space)
   end subroutine open_block

   subroutine read_block(dset, buffer)
      use, intrinsic :: iso_c_binding, only: c_ptr
      integer(hid_t), intent(in) :: dset
      type(c_ptr), intent(in) :: buffer

      integer(herr_t) :: ignored

      ignored = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer)
      ignored = H5Dclose(dset)
   end subroutine read_block

   subroutine amp_commit(this, iteration, energy)
      !! Declare everything written since the last commit to be a usable set
      !!
      !! The flush has to land before `complete` is raised: a `complete` that
      !! reaches the disk ahead of the arrays it describes is exactly the file
      !! this design exists to prevent.
      class(hdf5_amplitudes_t), intent(inout) :: this
      integer, intent(in) :: iteration
      real(dp), intent(in) :: energy

      integer(herr_t) :: ignored

      if (.not. this%active) return

      ignored = H5Fflush(this%file, H5F_SCOPE_LOCAL)
      call write_int_attr(this%file, "iteration", iteration)
      call write_double_attr(this%file, "energy", energy)
      call write_int_attr(this%file, "complete", 1)
      ignored = H5Fflush(this%file, H5F_SCOPE_LOCAL)

      this%iteration = iteration
      this%energy = energy
      this%resumable = .true.
   end subroutine amp_commit

   subroutine amp_close(this)
      class(hdf5_amplitudes_t), intent(inout) :: this
      integer(herr_t) :: ignored

      if (this%file >= 0) ignored = H5Fclose(this%file)
      this%file = -1
      this%active = .false.
   end subroutine amp_close

   ! -- attributes ----------------------------------------------------------
   ! Scalars live as attributes rather than one-element datasets: HDF5 keeps
   ! them in the object header, so reading `complete` on open costs no seek
   ! into the data, which is the whole point of checking it first.

   subroutine write_int_attr(file, name, value)
      integer(hid_t), intent(in) :: file
      character(len=*), intent(in) :: name
      integer, intent(in) :: value

      integer(hid_t) :: space, attr
      integer, target :: buffer(1)
      integer(herr_t) :: ignored

      buffer(1) = value
      if (H5Aexists(file, c_string(name)) > 0) then
         attr = H5Aopen(file, c_string(name), H5P_DEFAULT)
      else
         space = H5Screate(0)
         attr = H5Acreate2(file, c_string(name), H5T_STD_I32LE, space, &
                           H5P_DEFAULT, H5P_DEFAULT)
         ignored = H5Sclose(space)
      end if
      ignored = H5Awrite(attr, H5T_NATIVE_INT, c_loc(buffer))
      ignored = H5Aclose(attr)
   end subroutine write_int_attr

   subroutine read_int_attr(file, name, value)
      integer(hid_t), intent(in) :: file
      character(len=*), intent(in) :: name
      integer, intent(out) :: value

      integer(hid_t) :: attr
      integer, target :: buffer(1)
      integer(herr_t) :: ignored

      value = 0
      if (H5Aexists(file, c_string(name)) <= 0) return
      attr = H5Aopen(file, c_string(name), H5P_DEFAULT)
      buffer(1) = 0
      ignored = H5Aread(attr, H5T_NATIVE_INT, c_loc(buffer))
      ignored = H5Aclose(attr)
      value = buffer(1)
   end subroutine read_int_attr

   subroutine write_double_attr(file, name, value)
      integer(hid_t), intent(in) :: file
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: value

      integer(hid_t) :: space, attr
      real(dp), target :: buffer(1)
      integer(herr_t) :: ignored

      buffer(1) = value
      if (H5Aexists(file, c_string(name)) > 0) then
         attr = H5Aopen(file, c_string(name), H5P_DEFAULT)
      else
         space = H5Screate(0)
         attr = H5Acreate2(file, c_string(name), H5T_NATIVE_DOUBLE, space, &
                           H5P_DEFAULT, H5P_DEFAULT)
         ignored = H5Sclose(space)
      end if
      ignored = H5Awrite(attr, H5T_NATIVE_DOUBLE, c_loc(buffer))
      ignored = H5Aclose(attr)
   end subroutine write_double_attr

   subroutine read_double_attr(file, name, value)
      integer(hid_t), intent(in) :: file
      character(len=*), intent(in) :: name
      real(dp), intent(out) :: value

      integer(hid_t) :: attr
      real(dp), target :: buffer(1)
      integer(herr_t) :: ignored

      value = 0.0_dp
      if (H5Aexists(file, c_string(name)) <= 0) return
      attr = H5Aopen(file, c_string(name), H5P_DEFAULT)
      buffer(1) = 0.0_dp
      ignored = H5Aread(attr, H5T_NATIVE_DOUBLE, c_loc(buffer))
      ignored = H5Aclose(attr)
      value = buffer(1)
   end subroutine read_double_attr

   subroutine write_string_attr(file, name, text)
      ! TODO(mqc): creates unconditionally, where `write_int_attr` and
      ! `write_double_attr` reopen an existing attribute. Writing the same
      ! string attribute twice fails rather than overwriting.
      integer(hid_t), intent(in) :: file
      character(len=*), intent(in) :: name, text

      integer(hid_t) :: space, atype, attr
      character(kind=c_char), target :: buffer(FP_LEN)
      integer(herr_t) :: ignored

      buffer = c_string(text)
      atype = H5Tcopy(H5T_C_S1)
      ignored = H5Tset_size(atype, int(FP_LEN, c_size_t))
      space = H5Screate(0)
      attr = H5Acreate2(file, c_string(name), atype, space, H5P_DEFAULT, H5P_DEFAULT)
      ignored = H5Awrite(attr, atype, c_loc(buffer))
      ignored = H5Aclose(attr)
      ignored = H5Sclose(space)
      ignored = H5Tclose(atype)
   end subroutine write_string_attr

   subroutine read_string_attr(file, name, text)
      integer(hid_t), intent(in) :: file
      character(len=*), intent(in) :: name
      character(len=*), intent(out) :: text

      integer(hid_t) :: atype, attr
      character(kind=c_char), target :: buffer(FP_LEN)
      integer(herr_t) :: ignored
      integer :: i

      text = ""
      if (H5Aexists(file, c_string(name)) <= 0) return
      atype = H5Tcopy(H5T_C_S1)
      ignored = H5Tset_size(atype, int(FP_LEN, c_size_t))
      attr = H5Aopen(file, c_string(name), H5P_DEFAULT)
      ignored = H5Aread(attr, atype, c_loc(buffer))
      ignored = H5Aclose(attr)
      ignored = H5Tclose(atype)
      do i = 1, min(FP_LEN - 1, len(text))
         if (iachar(buffer(i)) == 0) exit
         text(i:i) = buffer(i)
      end do
   end subroutine read_string_attr

end module mqc_hdf5_amplitudes
