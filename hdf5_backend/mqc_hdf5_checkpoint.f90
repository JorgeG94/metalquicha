!! HDF5 backend for checkpoints that carry derivatives
module mqc_hdf5_checkpoint
   !! The same append-and-resume contract as the text checkpoint, in a format
   !! that can hold a gradient and a Hessian per fragment.
   !!
   !! **Why not text.** For MBE(3) over a thousand monomers the energies alone
   !! are 10 GB written as decimal and 1.3 GB as doubles. Gradients are 75 GB
   !! against 24, and Hessians 1350 GB against 432. At that point text stops
   !! being a slow format and starts being a refusal.
   !!
   !! **A few big datasets, never one per fragment.** HDF5 keeps metadata per
   !! object; a hundred million datasets would cost more than the numbers in
   !! them and take hours to open. Everything here is a handful of extendible
   !! one-dimensional datasets plus an offset index, so a fragment's gradient
   !! is a slice rather than an object.
   !!
   !! **`n_valid` is what makes a kill survivable, and it is the whole trick.**
   !! HDF5 is less forgiving than a text file: a process killed mid-write can
   !! leave metadata the library then refuses to open *at all*, losing every
   !! record rather than the last one. So records are appended, the file is
   !! flushed, and only then is `n_valid` raised. A reader trusts `n_valid`
   !! records and ignores whatever lies past it. The window that a kill can
   !! lose is therefore bounded by the flush interval, and what it loses is
   !! always a suffix -- never the file.
   use, intrinsic :: iso_c_binding, only: c_loc, c_null_ptr, c_size_t, c_char, c_int
   use pic_types, only: dp, int64
   use mqc_error, only: error_t, ERROR_IO, ERROR_VALIDATION
   use mqc_hdf5_bindings, only: hid_t, hsize_t, herr_t, h5_start, h5_version_text, c_string, &
                                H5F_ACC_TRUNC, H5F_ACC_RDWR, H5P_DEFAULT, H5S_UNLIMITED, &
                                H5S_SELECT_SET, H5F_SCOPE_LOCAL, H5T_NATIVE_DOUBLE, &
                                H5T_NATIVE_INT, H5T_NATIVE_LLONG, H5T_STD_I32LE, &
                                H5T_STD_I64LE, H5T_C_S1, H5P_CLS_DATASET_CREATE_ID, &
                                H5Fcreate, H5Fopen, H5Fclose, H5Fflush, &
                                H5Screate, H5Screate_simple, H5Sclose, H5Sselect_hyperslab, &
                                H5Dcreate2, H5Dopen2, H5Dclose, H5Dwrite, H5Dread, &
                                H5Dset_extent, H5Dget_space, H5Pcreate, H5Pset_chunk, H5Pclose, &
                                H5Acreate2, H5Aopen, H5Awrite, H5Aread, H5Aclose, H5Aexists, &
                                H5Tcopy, H5Tset_size, H5Tclose
   implicit none
   private

   public :: hdf5_checkpoint_t
   public :: hdf5_checkpoint_available

   integer, parameter :: FP_LEN = 17          !! 16 hex digits and a terminator
   integer(hsize_t), parameter :: CHUNK = 1024
      !! Dataset chunk, a storage decision: big enough that HDF5's per-chunk
      !! metadata stays small against the data.

   integer(int64), parameter :: COMMIT_EVERY = 1024_int64
   real(dp), parameter :: COMMIT_SECONDS = 30.0_dp
      !! How much a kill may cost, which is a different question from how the
      !! data is laid out and must not share its answer. Whichever comes
      !! first: a thousand records, or half a minute.
      !!
      !! Count alone is wrong at both ends. Fragments that take milliseconds
      !! make a commit-per-record cost more than the science -- at a hundred
      !! million records a millisecond flush each is a day of flushing. And
      !! fragments that take seconds, which is what DFT on a trimer costs,
      !! make a thousand records several hours, so a count-only rule would
      !! quietly reintroduce exactly the loss this file exists to prevent.
      !! The clock covers the expensive case, the count covers the cheap one.

   type :: hdf5_checkpoint_t
      !! One open HDF5 checkpoint
      logical :: active = .false.
      integer(hid_t) :: file = -1
      integer(hid_t) :: d_terms = -1, d_energy = -1, d_status = -1, d_natoms = -1
      integer(hid_t) :: d_goff = -1, d_grad = -1, d_hoff = -1, d_hess = -1
      integer :: max_level = 0

      integer(int64) :: n_written = 0     !! Records appended this session
      integer(int64) :: n_committed = 0   !! Records a reader would be told about
      integer(int64) :: last_tick = 0     !! system_clock at the last commit
      integer(int64) :: n_grad = 0        !! Doubles in /gradients
      integer(int64) :: n_hess = 0        !! Doubles in /hessians

      ! What was loaded, sorted by term for bisection
      integer(int64) :: n_loaded = 0
      integer, allocatable :: terms(:, :)
      real(dp), allocatable :: energies(:)
      integer, allocatable :: scf_status(:)
      integer, allocatable :: natoms(:)
      integer(int64), allocatable :: gstart(:), gcount(:)
      integer(int64), allocatable :: hstart(:), hcount(:)
         !! Per record, not cumulative. On disk the offsets are cumulative,
         !! which is compact; in memory the records get sorted by term, and a
         !! cumulative array cannot be permuted -- record i's slice would end
         !! up described by its neighbour's bounds. Converted once on load.
      real(dp), allocatable :: grad(:), hess(:)
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
      available = .true.
   end function hdf5_checkpoint_available

   subroutine h5ck_open(this, path, fingerprint, max_level, error)
      !! Load an existing file if there is one, then open for appending
      class(hdf5_checkpoint_t), intent(inout) :: this
      character(len=*), intent(in) :: path
      character(len=*), intent(in) :: fingerprint
      integer, intent(in) :: max_level
      type(error_t), intent(inout) :: error

      logical :: exists

      this%max_level = max_level
      if (.not. h5_start()) then
         call error%set(ERROR_IO, "HDF5 would not start, or is older than 1.10 "// &
                        "(found "//h5_version_text()//"); its ids would be the wrong width")
         return
      end if

      inquire (file=path, exist=exists)
      if (exists) then
         call load(this, path, fingerprint, error)
         if (error%has_error()) return
      else
         call create(this, path, fingerprint, error)
         if (error%has_error()) return
      end if
      this%active = .true.
   end subroutine h5ck_open

   subroutine create(this, path, fingerprint, error)
      !! A new file with empty extendible datasets
      class(hdf5_checkpoint_t), intent(inout) :: this
      character(len=*), intent(in) :: path, fingerprint
      type(error_t), intent(inout) :: error

      this%file = H5Fcreate(c_string(path), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)
      if (this%file < 0) then
         call error%set(ERROR_IO, "could not create checkpoint "//trim(path))
         return
      end if

      call write_string_attr(this%file, "fingerprint", fingerprint)
      call write_count(this%file, 0_int64)
      call system_clock(this%last_tick)

      this%d_terms = make_2d("terms", this%file, int(this%max_level, hsize_t), H5T_STD_I32LE)
      this%d_energy = make_1d("energies", this%file, H5T_NATIVE_DOUBLE)
      this%d_status = make_1d("scf_status", this%file, H5T_STD_I32LE)
      this%d_natoms = make_1d("n_atoms", this%file, H5T_STD_I32LE)
      this%d_goff = make_1d("grad_offset", this%file, H5T_STD_I64LE)
      this%d_grad = make_1d("gradients", this%file, H5T_NATIVE_DOUBLE)
      this%d_hoff = make_1d("hess_offset", this%file, H5T_STD_I64LE)
      this%d_hess = make_1d("hessians", this%file, H5T_NATIVE_DOUBLE)

      if (min(this%d_terms, this%d_energy, this%d_status, this%d_natoms, &
              this%d_goff, this%d_grad, this%d_hoff, this%d_hess) < 0) then
         call error%set(ERROR_IO, "could not create the checkpoint datasets")
      end if
   end subroutine create

   subroutine load(this, path, fingerprint, error)
      !! Read back a previous run, trusting only the first `n_valid` records
      class(hdf5_checkpoint_t), intent(inout) :: this
      character(len=*), intent(in) :: path, fingerprint
      type(error_t), intent(inout) :: error

      character(len=FP_LEN) :: stored
      integer(int64) :: n

      this%file = H5Fopen(c_string(path), H5F_ACC_RDWR, H5P_DEFAULT)
      if (this%file < 0) then
         call error%set(ERROR_IO, trim(path)//" exists but is not a readable HDF5 file")
         return
      end if

      call read_string_attr(this%file, "fingerprint", stored)
      if (trim(stored) /= trim(fingerprint)) then
         call error%set(ERROR_VALIDATION, "checkpoint "//trim(path)//" was written by a "// &
                        "different calculation (fingerprint "//trim(stored)//", this run is "// &
                        trim(fingerprint)//"). Its energies belong to another geometry, "// &
                        "basis or method; reusing them would give a converged total for "// &
                        "neither. Point system.checkpoint at another file, or delete this one.")
         return
      end if

      ! Anything past n_valid was being written when the last run stopped. It
      ! may be complete, it may be half a chunk of zeros; it is not trusted
      ! either way, which is what keeps a kill from costing more than a suffix.
      call read_count(this%file, n)

      this%d_terms = H5Dopen2(this%file, c_string("terms"), H5P_DEFAULT)
      this%d_energy = H5Dopen2(this%file, c_string("energies"), H5P_DEFAULT)
      this%d_status = H5Dopen2(this%file, c_string("scf_status"), H5P_DEFAULT)
      this%d_natoms = H5Dopen2(this%file, c_string("n_atoms"), H5P_DEFAULT)
      this%d_goff = H5Dopen2(this%file, c_string("grad_offset"), H5P_DEFAULT)
      this%d_grad = H5Dopen2(this%file, c_string("gradients"), H5P_DEFAULT)
      this%d_hoff = H5Dopen2(this%file, c_string("hess_offset"), H5P_DEFAULT)
      this%d_hess = H5Dopen2(this%file, c_string("hessians"), H5P_DEFAULT)
      if (min(this%d_terms, this%d_energy, this%d_natoms) < 0) then
         call error%set(ERROR_VALIDATION, trim(path)//" is not an mqc checkpoint")
         return
      end if

      this%n_loaded = n
      this%n_written = n
      this%n_committed = n
      call system_clock(this%last_tick)
      if (n <= 0) return

      allocate (this%terms(n, this%max_level))
      allocate (this%energies(n), this%scf_status(n), this%natoms(n))
      allocate (this%gstart(n), this%gcount(n), this%hstart(n), this%hcount(n))
      call read_2d_int(this%d_terms, n, int(this%max_level, hsize_t), this%terms)
      call read_1d_double(this%d_energy, n, this%energies)
      call read_1d_int(this%d_status, n, this%scf_status)
      call read_1d_int(this%d_natoms, n, this%natoms)
      block
         integer(int64), allocatable :: goff(:), hoff(:)
         integer(int64) :: i
         allocate (goff(n + 1), hoff(n + 1))
         call read_1d_long(this%d_goff, n + 1, goff)
         call read_1d_long(this%d_hoff, n + 1, hoff)
         do i = 1, n
            this%gstart(i) = goff(i)
            this%gcount(i) = goff(i + 1) - goff(i)
            this%hstart(i) = hoff(i)
            this%hcount(i) = hoff(i + 1) - hoff(i)
         end do
         this%n_grad = goff(n + 1)
         this%n_hess = hoff(n + 1)
      end block
      if (this%n_grad > 0) then
         allocate (this%grad(this%n_grad))
         call read_1d_double(this%d_grad, this%n_grad, this%grad)
      end if
      if (this%n_hess > 0) then
         allocate (this%hess(this%n_hess))
         call read_1d_double(this%d_hess, this%n_hess, this%hess)
      end if

      call sort_loaded(this)
   end subroutine load

   subroutine h5ck_record(this, term, energy, scf_status, n_atoms, gradient, hessian)
      !! Append one finished fragment
      !!
      !! The order is append, flush, then raise the count -- never the other
      !! way round. A count raised before the data is on disk would promise a
      !! record that is not there, which is exactly the lie the whole scheme
      !! exists to prevent.
      class(hdf5_checkpoint_t), intent(inout) :: this
      integer, intent(in) :: term(:)
      real(dp), intent(in) :: energy
      integer, intent(in) :: scf_status
      integer, intent(in) :: n_atoms
      real(dp), intent(in), optional :: gradient(:, :)
      real(dp), intent(in), optional :: hessian(:, :)

      integer(int64) :: at, ng, nh
      integer :: row(this%max_level)
      real(dp) :: scalar(1)
      integer :: ivalue(1)
      integer(int64) :: lvalue(1)

      if (.not. this%active) return
      at = this%n_written
      row = 0
      row(1:min(size(term), this%max_level)) = term(1:min(size(term), this%max_level))

      call append_2d_int(this%d_terms, at, int(this%max_level, hsize_t), row)
      scalar(1) = energy
      call append_1d_double(this%d_energy, at, 1_int64, scalar)
      ivalue(1) = scf_status
      call append_1d_int(this%d_status, at, 1_int64, ivalue)
      ivalue(1) = n_atoms
      call append_1d_int(this%d_natoms, at, 1_int64, ivalue)

      ! Offsets are written as [start, end] so that record i owns
      ! goff(i)..goff(i+1)-1 whatever order the reader arrives in.
      if (at == 0_int64) then
         lvalue(1) = 0_int64
         call append_1d_long(this%d_goff, 0_int64, 1_int64, lvalue)
         call append_1d_long(this%d_hoff, 0_int64, 1_int64, lvalue)
      end if

      ng = 0
      if (present(gradient)) then
         ng = size(gradient, kind=int64)
         call append_1d_double(this%d_grad, this%n_grad, ng, reshape(gradient, [ng]))
      end if
      this%n_grad = this%n_grad + ng
      lvalue(1) = this%n_grad
      call append_1d_long(this%d_goff, at + 1_int64, 1_int64, lvalue)

      nh = 0
      if (present(hessian)) then
         nh = size(hessian, kind=int64)
         call append_1d_double(this%d_hess, this%n_hess, nh, reshape(hessian, [nh]))
      end if
      this%n_hess = this%n_hess + nh
      lvalue(1) = this%n_hess
      call append_1d_long(this%d_hoff, at + 1_int64, 1_int64, lvalue)

      this%n_written = at + 1_int64

      if (due(this)) call commit(this)
   end subroutine h5ck_record

   function due(this) result(now)
      !! Whether enough records, or enough time, have passed to commit
      class(hdf5_checkpoint_t), intent(in) :: this
      logical :: now

      integer(int64) :: tick, rate

      now = (this%n_written - this%n_committed >= COMMIT_EVERY)
      if (now) return

      call system_clock(tick, rate)
      if (rate <= 0) return
      now = (real(tick - this%last_tick, dp)/real(rate, dp) >= COMMIT_SECONDS)
   end function due

   subroutine commit(this)
      !! Get the data down, then say how much of it is real
      !!
      !! This order is the entire guarantee. A count raised before the flush
      !! would promise records that are not on disk, and a reader that
      !! believed it would splice uninitialised numbers into an expansion.
      class(hdf5_checkpoint_t), intent(inout) :: this

      integer(int64) :: rate

      if (H5Fflush(this%file, H5F_SCOPE_LOCAL) < 0) return
      call write_count(this%file, this%n_written)
      this%n_committed = this%n_written
      call system_clock(this%last_tick, rate)
   end subroutine commit

   subroutine h5ck_lookup(this, term, found, energy, scf_status, n_atoms, gradient, hessian)
      !! Is this term already done, and with what
      class(hdf5_checkpoint_t), intent(in) :: this
      integer, intent(in) :: term(:)
      logical, intent(out) :: found
      real(dp), intent(out) :: energy
      integer, intent(out) :: scf_status
      integer, intent(out) :: n_atoms
      real(dp), allocatable, intent(out), optional :: gradient(:, :)
      real(dp), allocatable, intent(out), optional :: hessian(:, :)

      integer(int64) :: lo, hi, mid, a, b, n3
      integer :: order

      found = .false.
      energy = 0.0_dp
      scf_status = 0
      n_atoms = 0
      if (this%n_loaded <= 0) return

      lo = 1
      hi = this%n_loaded
      do while (lo <= hi)
         mid = (lo + hi)/2
         order = compare(this%terms(mid, :), term, this%max_level)
         if (order == 0) then
            found = .true.
            energy = this%energies(mid)
            scf_status = this%scf_status(mid)
            n_atoms = this%natoms(mid)
            if (present(gradient) .and. this%gcount(mid) > 0) then
               a = this%gstart(mid) + 1
               b = a + this%gcount(mid) - 1
               gradient = reshape(this%grad(a:b), [3, n_atoms])
            end if
            if (present(hessian) .and. this%hcount(mid) > 0) then
               a = this%hstart(mid) + 1
               b = a + this%hcount(mid) - 1
               n3 = 3*n_atoms
               hessian = reshape(this%hess(a:b), [n3, n3])
            end if
            return
         end if
         if (order < 0) then
            lo = mid + 1
         else
            hi = mid - 1
         end if
      end do
   end subroutine h5ck_lookup

   subroutine h5ck_close(this)
      !! Commit what has not been committed, then close everything
      class(hdf5_checkpoint_t), intent(inout) :: this

      integer(herr_t) :: ignored

      if (.not. this%active) return
      call commit(this)
      ignored = H5Dclose(this%d_terms)
      ignored = H5Dclose(this%d_energy)
      ignored = H5Dclose(this%d_status)
      ignored = H5Dclose(this%d_natoms)
      ignored = H5Dclose(this%d_goff)
      ignored = H5Dclose(this%d_grad)
      ignored = H5Dclose(this%d_hoff)
      ignored = H5Dclose(this%d_hess)
      ignored = H5Fclose(this%file)
      this%active = .false.
      this%file = -1
   end subroutine h5ck_close

   ! ==========================================================================
   !  Dataset plumbing
   ! ==========================================================================

   function make_1d(name, file, type_id) result(dset)
      character(len=*), intent(in) :: name
      integer(hid_t), intent(in) :: file, type_id
      integer(hid_t) :: dset

      integer(hid_t) :: space, dcpl
      integer(hsize_t) :: dims(1), maxdims(1), chunk_dims(1)
      integer(herr_t) :: ignored

      dims = [0_hsize_t]
      maxdims = [H5S_UNLIMITED]
      space = H5Screate_simple(1, dims, maxdims)
      chunk_dims = [CHUNK]
      dcpl = H5Pcreate(H5P_CLS_DATASET_CREATE_ID)
      ignored = H5Pset_chunk(dcpl, 1, chunk_dims)
      dset = H5Dcreate2(file, c_string(name), type_id, space, H5P_DEFAULT, dcpl, H5P_DEFAULT)
      ignored = H5Pclose(dcpl)
      ignored = H5Sclose(space)
   end function make_1d

   function make_2d(name, file, width, type_id) result(dset)
      character(len=*), intent(in) :: name
      integer(hid_t), intent(in) :: file, type_id
      integer(hsize_t), intent(in) :: width
      integer(hid_t) :: dset

      integer(hid_t) :: space, dcpl
      integer(hsize_t) :: dims(2), maxdims(2), chunk_dims(2)
      integer(herr_t) :: ignored

      ! Row-major on the C side, so a record is contiguous and one append
      ! writes one row.
      dims = [0_hsize_t, width]
      maxdims = [H5S_UNLIMITED, width]
      space = H5Screate_simple(2, dims, maxdims)
      chunk_dims = [CHUNK, width]
      dcpl = H5Pcreate(H5P_CLS_DATASET_CREATE_ID)
      ignored = H5Pset_chunk(dcpl, 2, chunk_dims)
      dset = H5Dcreate2(file, c_string(name), type_id, space, H5P_DEFAULT, dcpl, H5P_DEFAULT)
      ignored = H5Pclose(dcpl)
      ignored = H5Sclose(space)
   end function make_2d

   subroutine append_1d_double(dset, at, n, buffer)
      integer(hid_t), intent(in) :: dset
      integer(int64), intent(in) :: at, n
      real(dp), intent(in), target :: buffer(:)

      integer(hid_t) :: fspace, mspace
      integer(hsize_t) :: dims(1), start(1), count(1)
      integer(herr_t) :: ignored

      if (n <= 0) return
      dims = [int(at + n, hsize_t)]
      ignored = H5Dset_extent(dset, dims)
      fspace = H5Dget_space(dset)
      start = [int(at, hsize_t)]
      count = [int(n, hsize_t)]
      ignored = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, c_null_ptr, count, c_null_ptr)
      mspace = H5Screate_simple(1, count, count)
      ignored = H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, c_loc(buffer))
      ignored = H5Sclose(mspace)
      ignored = H5Sclose(fspace)
   end subroutine append_1d_double

   subroutine append_1d_int(dset, at, n, buffer)
      integer(hid_t), intent(in) :: dset
      integer(int64), intent(in) :: at, n
      integer, intent(in), target :: buffer(:)

      integer(hid_t) :: fspace, mspace
      integer(hsize_t) :: dims(1), start(1), count(1)
      integer(herr_t) :: ignored

      if (n <= 0) return
      dims = [int(at + n, hsize_t)]
      ignored = H5Dset_extent(dset, dims)
      fspace = H5Dget_space(dset)
      start = [int(at, hsize_t)]
      count = [int(n, hsize_t)]
      ignored = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, c_null_ptr, count, c_null_ptr)
      mspace = H5Screate_simple(1, count, count)
      ignored = H5Dwrite(dset, H5T_NATIVE_INT, mspace, fspace, H5P_DEFAULT, c_loc(buffer))
      ignored = H5Sclose(mspace)
      ignored = H5Sclose(fspace)
   end subroutine append_1d_int

   subroutine append_1d_long(dset, at, n, buffer)
      integer(hid_t), intent(in) :: dset
      integer(int64), intent(in) :: at, n
      integer(int64), intent(in), target :: buffer(:)

      integer(hid_t) :: fspace, mspace
      integer(hsize_t) :: dims(1), start(1), count(1)
      integer(herr_t) :: ignored

      if (n <= 0) return
      dims = [int(at + n, hsize_t)]
      ignored = H5Dset_extent(dset, dims)
      fspace = H5Dget_space(dset)
      start = [int(at, hsize_t)]
      count = [int(n, hsize_t)]
      ignored = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, c_null_ptr, count, c_null_ptr)
      mspace = H5Screate_simple(1, count, count)
      ignored = H5Dwrite(dset, H5T_NATIVE_LLONG, mspace, fspace, H5P_DEFAULT, c_loc(buffer))
      ignored = H5Sclose(mspace)
      ignored = H5Sclose(fspace)
   end subroutine append_1d_long

   subroutine append_2d_int(dset, at, width, row)
      integer(hid_t), intent(in) :: dset
      integer(int64), intent(in) :: at
      integer(hsize_t), intent(in) :: width
      integer, intent(in), target :: row(:)

      integer(hid_t) :: fspace, mspace
      integer(hsize_t) :: dims(2), start(2), count(2)
      integer(herr_t) :: ignored

      dims = [int(at + 1, hsize_t), width]
      ignored = H5Dset_extent(dset, dims)
      fspace = H5Dget_space(dset)
      start = [int(at, hsize_t), 0_hsize_t]
      count = [1_hsize_t, width]
      ignored = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, c_null_ptr, count, c_null_ptr)
      mspace = H5Screate_simple(2, count, count)
      ignored = H5Dwrite(dset, H5T_NATIVE_INT, mspace, fspace, H5P_DEFAULT, c_loc(row))
      ignored = H5Sclose(mspace)
      ignored = H5Sclose(fspace)
   end subroutine append_2d_int

   subroutine read_1d_double(dset, n, buffer)
      integer(hid_t), intent(in) :: dset
      integer(int64), intent(in) :: n
      real(dp), intent(out), target :: buffer(:)

      integer(hid_t) :: fspace, mspace
      integer(hsize_t) :: start(1), count(1)
      integer(herr_t) :: ignored

      if (n <= 0) return
      fspace = H5Dget_space(dset)
      start = [0_hsize_t]
      count = [int(n, hsize_t)]
      ignored = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, c_null_ptr, count, c_null_ptr)
      mspace = H5Screate_simple(1, count, count)
      ignored = H5Dread(dset, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, c_loc(buffer))
      ignored = H5Sclose(mspace)
      ignored = H5Sclose(fspace)
   end subroutine read_1d_double

   subroutine read_1d_int(dset, n, buffer)
      integer(hid_t), intent(in) :: dset
      integer(int64), intent(in) :: n
      integer, intent(out), target :: buffer(:)

      integer(hid_t) :: fspace, mspace
      integer(hsize_t) :: start(1), count(1)
      integer(herr_t) :: ignored

      if (n <= 0) return
      fspace = H5Dget_space(dset)
      start = [0_hsize_t]
      count = [int(n, hsize_t)]
      ignored = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, c_null_ptr, count, c_null_ptr)
      mspace = H5Screate_simple(1, count, count)
      ignored = H5Dread(dset, H5T_NATIVE_INT, mspace, fspace, H5P_DEFAULT, c_loc(buffer))
      ignored = H5Sclose(mspace)
      ignored = H5Sclose(fspace)
   end subroutine read_1d_int

   subroutine read_1d_long(dset, n, buffer)
      integer(hid_t), intent(in) :: dset
      integer(int64), intent(in) :: n
      integer(int64), intent(out), target :: buffer(:)

      integer(hid_t) :: fspace, mspace
      integer(hsize_t) :: start(1), count(1)
      integer(herr_t) :: ignored

      buffer = 0_int64
      if (n <= 0) return
      fspace = H5Dget_space(dset)
      start = [0_hsize_t]
      count = [int(n, hsize_t)]
      ignored = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, c_null_ptr, count, c_null_ptr)
      mspace = H5Screate_simple(1, count, count)
      ignored = H5Dread(dset, H5T_NATIVE_LLONG, mspace, fspace, H5P_DEFAULT, c_loc(buffer))
      ignored = H5Sclose(mspace)
      ignored = H5Sclose(fspace)
   end subroutine read_1d_long

   subroutine read_2d_int(dset, n, width, buffer)
      integer(hid_t), intent(in) :: dset
      integer(int64), intent(in) :: n
      integer(hsize_t), intent(in) :: width
      integer, intent(out) :: buffer(:, :)

      integer, allocatable, target :: flat(:)
      integer(hid_t) :: fspace, mspace
      integer(hsize_t) :: start(2), count(2)
      integer(herr_t) :: ignored
      integer(int64) :: i
      integer :: j

      if (n <= 0) return
      allocate (flat(n*width))
      fspace = H5Dget_space(dset)
      start = [0_hsize_t, 0_hsize_t]
      count = [int(n, hsize_t), width]
      ignored = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, c_null_ptr, count, c_null_ptr)
      mspace = H5Screate_simple(2, count, count)
      ignored = H5Dread(dset, H5T_NATIVE_INT, mspace, fspace, H5P_DEFAULT, c_loc(flat))
      ignored = H5Sclose(mspace)
      ignored = H5Sclose(fspace)

      ! HDF5 hands back rows contiguously; Fortran wants columns.
      do i = 1, n
         do j = 1, int(width)
            buffer(i, j) = flat((i - 1)*width + j)
         end do
      end do
      deallocate (flat)
   end subroutine read_2d_int

   subroutine write_count(file, n)
      !! Raise `n_valid`, the number of records a reader may trust
      integer(hid_t), intent(in) :: file
      integer(int64), intent(in) :: n

      integer(hid_t) :: space, attr
      integer(int64), target :: value(1)
      integer(herr_t) :: ignored

      value(1) = n
      if (H5Aexists(file, c_string("n_valid")) > 0) then
         attr = H5Aopen(file, c_string("n_valid"), H5P_DEFAULT)
      else
         space = H5Screate(0)
         attr = H5Acreate2(file, c_string("n_valid"), H5T_STD_I64LE, space, &
                           H5P_DEFAULT, H5P_DEFAULT)
         ignored = H5Sclose(space)
      end if
      ignored = H5Awrite(attr, H5T_NATIVE_LLONG, c_loc(value))
      ignored = H5Aclose(attr)
   end subroutine write_count

   subroutine read_count(file, n)
      integer(hid_t), intent(in) :: file
      integer(int64), intent(out) :: n

      integer(hid_t) :: attr
      integer(int64), target :: value(1)
      integer(herr_t) :: ignored

      n = 0
      if (H5Aexists(file, c_string("n_valid")) <= 0) return
      attr = H5Aopen(file, c_string("n_valid"), H5P_DEFAULT)
      value(1) = 0
      ignored = H5Aread(attr, H5T_NATIVE_LLONG, c_loc(value))
      ignored = H5Aclose(attr)
      n = value(1)
   end subroutine read_count

   subroutine write_string_attr(file, name, text)
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

   pure function compare(a, b, n) result(order)
      integer, intent(in) :: a(:), b(:)
      integer, intent(in) :: n
      integer :: order

      integer :: i

      order = 0
      do i = 1, n
         if (a(i) < b(i)) then
            order = -1
            return
         end if
         if (a(i) > b(i)) then
            order = 1
            return
         end if
      end do
   end function compare

   subroutine sort_loaded(this)
      !! Sort every loaded array together, by term
      !!
      !! An index permutation rather than moving the payload: the gradients
      !! and Hessians are addressed through the offsets, so those stay put and
      !! only the small per-record arrays are reordered.
      class(hdf5_checkpoint_t), intent(inout) :: this

      integer(int64) :: i, j
      integer :: key_term(this%max_level)
      integer :: key_status, key_natoms
      real(dp) :: key_energy
      integer(int64) :: key_goff(2), key_hoff(2)

      do i = 2, this%n_loaded
         key_term = this%terms(i, :)
         key_energy = this%energies(i)
         key_status = this%scf_status(i)
         key_natoms = this%natoms(i)
         key_goff = [this%gstart(i), this%gcount(i)]
         key_hoff = [this%hstart(i), this%hcount(i)]
         j = i - 1
         do while (j >= 1)
            if (compare(this%terms(j, :), key_term, this%max_level) <= 0) exit
            this%terms(j + 1, :) = this%terms(j, :)
            this%energies(j + 1) = this%energies(j)
            this%scf_status(j + 1) = this%scf_status(j)
            this%natoms(j + 1) = this%natoms(j)
            this%gstart(j + 1) = this%gstart(j)
            this%gcount(j + 1) = this%gcount(j)
            this%hstart(j + 1) = this%hstart(j)
            this%hcount(j + 1) = this%hcount(j)
            j = j - 1
         end do
         this%terms(j + 1, :) = key_term
         this%energies(j + 1) = key_energy
         this%scf_status(j + 1) = key_status
         this%natoms(j + 1) = key_natoms
         this%gstart(j + 1) = key_goff(1)
         this%gcount(j + 1) = key_goff(2)
         this%hstart(j + 1) = key_hoff(1)
         this%hcount(j + 1) = key_hoff(2)
      end do
   end subroutine sort_loaded

end module mqc_hdf5_checkpoint
