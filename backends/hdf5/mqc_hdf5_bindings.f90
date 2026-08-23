!! Hand-written HDF5 C bindings
module mqc_hdf5_bindings
   !! The slice of the HDF5 C API this program needs, bound directly.
   !!
   !! **The C API rather than HDF5's Fortran interface, deliberately.**
   !! `hdf5.mod` is compiled by, and only readable by, the compiler that built
   !! it -- a gfortran-built module cannot be used from ifx or nvfortran, and
   !! this project claims all of them. Binding the C entry points sidesteps
   !! that entirely: the ABI is stable and compiler-independent, and the
   !! machine that has libhdf5 but no matching module is the normal case
   !! rather than the broken one. The cuBLAS and cuSOLVER bindings in this
   !! tree were written the same way for the same reason.
   !!
   !! **The type constants are variables, not parameters.** `H5T_NATIVE_DOUBLE`
   !! is a C macro over a global that the library fills in during
   !! initialisation, so reading it before `H5open` yields zero -- and a zero
   !! type id is not an error, it is a different type. The Fortran interface
   !! and the C macros both hide an implicit `H5open`; binding the raw globals
   !! does not, so `h5_start` exists and must be called first.
   !!
   !! Sizes taken from the installed headers rather than assumed:
   !! `hid_t` is int64_t, `herr_t` is int, `hsize_t` is unsigned long long.
   !! Getting one of those wrong is not a compile error, it is silent memory
   !! corruption, which is why they are written down here with their source.
   use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_char, c_ptr, c_size_t, c_null_char
   implicit none
   private

   public :: h5_version, h5_version_text
   public :: hid_t, herr_t, hsize_t, H5F_ACC_RDONLY, H5F_ACC_RDWR, H5F_ACC_TRUNC, &
             H5F_ACC_EXCL, H5P_DEFAULT, H5S_ALL, H5S_UNLIMITED, H5D_CHUNKED, H5S_SELECT_SET, &
             H5F_SCOPE_LOCAL, H5T_NATIVE_DOUBLE, H5T_NATIVE_INT, H5T_NATIVE_LLONG, &
             H5T_STD_I32LE, H5T_STD_I64LE, H5T_C_S1, H5P_CLS_DATASET_CREATE_ID, h5_start, &
             c_string, H5Aclose, H5Acreate2, H5Aexists, H5Aopen, H5Aread, H5Awrite, &
             H5Dclose, H5Dcreate2, H5Dget_space, H5Dopen2, H5Dread, H5Dset_extent, H5Dwrite, &
             H5Fclose, H5Fcreate, H5Fflush, H5Fopen, H5Lexists, H5Pclose, H5Pcreate, &
             H5Pset_chunk, H5Sclose, H5Screate, H5Screate_simple, H5Sselect_hyperslab, &
             H5Sget_simple_extent_dims, H5Sget_simple_extent_ndims, &
             H5Tclose, H5Tcopy, H5Tset_size, H5open

   !> Assumed-size dummies below are allowed deliberately and individually.
   !  `dimension(*)` is the correct interop declaration for a C parameter that
   !  is a plain pointer: the assumed-shape alternative the rule recommends is
   !  passed as a descriptor under `bind(C)`, which is not what HDF5 expects
   !  and would corrupt every call. The rule is right about Fortran and wrong
   !  about this boundary.

   integer, parameter :: hid_t = c_int64_t     !! H5Ipublic.h: typedef int64_t hid_t
   integer, parameter :: herr_t = c_int        !! H5public.h:  typedef int herr_t
   integer, parameter :: hsize_t = c_int64_t   !! H5public.h:  unsigned long long

   !> Flags, from H5Fpublic.h. Unsigned in C; these values all fit in a
   !  positive default integer, so the kind mismatch cannot bite.
   integer(c_int), parameter :: H5F_ACC_RDONLY = 0
   integer(c_int), parameter :: H5F_ACC_RDWR = 1
   integer(c_int), parameter :: H5F_ACC_TRUNC = 2
   integer(c_int), parameter :: H5F_ACC_EXCL = 4

   integer(hid_t), parameter :: H5P_DEFAULT = 0_hid_t   !! H5Ppublic.h
   integer(hid_t), parameter :: H5S_ALL = 0_hid_t       !! H5Spublic.h
   integer(hsize_t), parameter :: H5S_UNLIMITED = -1_hsize_t
      !! HSIZE_UNDEF is ULLONG_MAX; as a signed 64-bit value that is -1, and
      !! the bit pattern is what crosses the boundary.

   integer(c_int), parameter :: H5D_CHUNKED = 2        !! H5D_layout_t
   integer(c_int), parameter :: H5S_SELECT_SET = 0     !! H5S_seloper_t
   integer(c_int), parameter :: H5F_SCOPE_LOCAL = 0    !! H5F_scope_t

   !> Datatype ids, filled in by the library at start-up. See the note above:
   !  these are zero until `h5_start` has run.
   integer(hid_t), bind(C, name="H5T_NATIVE_DOUBLE_g") :: H5T_NATIVE_DOUBLE
   integer(hid_t), bind(C, name="H5T_NATIVE_INT_g") :: H5T_NATIVE_INT
   integer(hid_t), bind(C, name="H5T_NATIVE_LLONG_g") :: H5T_NATIVE_LLONG
   integer(hid_t), bind(C, name="H5T_STD_I32LE_g") :: H5T_STD_I32LE
   integer(hid_t), bind(C, name="H5T_STD_I64LE_g") :: H5T_STD_I64LE
   integer(hid_t), bind(C, name="H5T_C_S1_g") :: H5T_C_S1

   !> Property-list classes are globals too, for the same reason and with the
   !  same trap. `H5P_DATASET_CREATE` is a macro over this in C.
   integer(hid_t), bind(C, name="H5P_CLS_DATASET_CREATE_ID_g") :: H5P_CLS_DATASET_CREATE_ID

   interface

      ! -- library ----------------------------------------------------------
      function H5open() result(status) bind(C, name="H5open")
         import :: herr_t
         implicit none
         integer(herr_t) :: status
      end function H5open

      function H5get_libversion(major, minor, release) result(status) &
         bind(C, name="H5get_libversion")
         import :: herr_t, c_int
         implicit none
         integer(c_int), intent(out) :: major, minor, release
         integer(herr_t) :: status
      end function H5get_libversion

      ! -- files ------------------------------------------------------------
      function H5Fcreate(name, flags, fcpl, fapl) result(file) bind(C, name="H5Fcreate")
         import :: hid_t, c_int, c_char
         implicit none
         ! allow(assumed-size)
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), value :: flags
         integer(hid_t), value :: fcpl, fapl
         integer(hid_t) :: file
      end function H5Fcreate

      function H5Fopen(name, flags, fapl) result(file) bind(C, name="H5Fopen")
         import :: hid_t, c_int, c_char
         implicit none
         ! allow(assumed-size)
         character(kind=c_char), intent(in) :: name(*)
         integer(c_int), value :: flags
         integer(hid_t), value :: fapl
         integer(hid_t) :: file
      end function H5Fopen

      function H5Fclose(file) result(status) bind(C, name="H5Fclose")
         import :: hid_t, herr_t
         implicit none
         integer(hid_t), value :: file
         integer(herr_t) :: status
      end function H5Fclose

      function H5Fflush(object, scope) result(status) bind(C, name="H5Fflush")
         !! Pushes the library's cache to the file. This is what bounds how
         !! much a killed process loses -- without it, HDF5 can leave metadata
         !! the library will later refuse to open at all, losing not the last
         !! record but every one of them.
         import :: hid_t, herr_t, c_int
         implicit none
         integer(hid_t), value :: object
         integer(c_int), value :: scope
         integer(herr_t) :: status
      end function H5Fflush

      ! -- dataspaces -------------------------------------------------------
      function H5Screate_simple(rank, dims, maxdims) result(space) &
         bind(C, name="H5Screate_simple")
         import :: hid_t, hsize_t, c_int
         implicit none
         integer(c_int), value :: rank
         ! allow(assumed-size)
         integer(hsize_t), intent(in) :: dims(*)
         ! allow(assumed-size)
         integer(hsize_t), intent(in) :: maxdims(*)
         integer(hid_t) :: space
      end function H5Screate_simple

      function H5Screate(class) result(space) bind(C, name="H5Screate")
         import :: hid_t, c_int
         implicit none
         integer(c_int), value :: class
         integer(hid_t) :: space
      end function H5Screate

      function H5Sclose(space) result(status) bind(C, name="H5Sclose")
         import :: hid_t, herr_t
         implicit none
         integer(hid_t), value :: space
         integer(herr_t) :: status
      end function H5Sclose

      !> Shape of a stored dataset, so a reader can allocate before it reads.
      !  `maxdims` is a real array rather than an optional C null: nothing here
      !  needs the unlimited bounds, and a dummy array costs one stack slot
      !  against a c_ptr dance at every call site.
      function H5Sget_simple_extent_ndims(space) result(rank) &
         bind(C, name="H5Sget_simple_extent_ndims")
         import :: hid_t, c_int
         implicit none
         integer(hid_t), value :: space
         integer(c_int) :: rank
      end function H5Sget_simple_extent_ndims

      function H5Sget_simple_extent_dims(space, dims, maxdims) result(rank) &
         bind(C, name="H5Sget_simple_extent_dims")
         import :: hid_t, hsize_t, c_int
         implicit none
         integer(hid_t), value :: space
         ! allow(assumed-size)
         integer(hsize_t), intent(out) :: dims(*)
         ! allow(assumed-size)
         integer(hsize_t), intent(out) :: maxdims(*)
         integer(c_int) :: rank
      end function H5Sget_simple_extent_dims

      function H5Sselect_hyperslab(space, op, start, stride, count, block) result(status) &
         bind(C, name="H5Sselect_hyperslab")
         import :: hid_t, herr_t, hsize_t, c_int, c_ptr
         implicit none
         integer(hid_t), value :: space
         integer(c_int), value :: op
         ! allow(assumed-size)
         integer(hsize_t), intent(in) :: start(*)
         type(c_ptr), value :: stride       !! C_NULL_PTR for contiguous
         ! allow(assumed-size)
         integer(hsize_t), intent(in) :: count(*)
         type(c_ptr), value :: block        !! C_NULL_PTR for single elements
         integer(herr_t) :: status
      end function H5Sselect_hyperslab

      ! -- datasets ---------------------------------------------------------
      function H5Dcreate2(loc, name, type_id, space, lcpl, dcpl, dapl) result(dset) &
         bind(C, name="H5Dcreate2")
         import :: hid_t, c_char
         implicit none
         integer(hid_t), value :: loc
         ! allow(assumed-size)
         character(kind=c_char), intent(in) :: name(*)
         integer(hid_t), value :: type_id, space, lcpl, dcpl, dapl
         integer(hid_t) :: dset
      end function H5Dcreate2

      function H5Dopen2(loc, name, dapl) result(dset) bind(C, name="H5Dopen2")
         import :: hid_t, c_char
         implicit none
         integer(hid_t), value :: loc
         ! allow(assumed-size)
         character(kind=c_char), intent(in) :: name(*)
         integer(hid_t), value :: dapl
         integer(hid_t) :: dset
      end function H5Dopen2

      function H5Dclose(dset) result(status) bind(C, name="H5Dclose")
         import :: hid_t, herr_t
         implicit none
         integer(hid_t), value :: dset
         integer(herr_t) :: status
      end function H5Dclose

      function H5Dget_space(dset) result(space) bind(C, name="H5Dget_space")
         import :: hid_t
         implicit none
         integer(hid_t), value :: dset
         integer(hid_t) :: space
      end function H5Dget_space

      function H5Dset_extent(dset, size) result(status) bind(C, name="H5Dset_extent")
         import :: hid_t, herr_t, hsize_t
         implicit none
         integer(hid_t), value :: dset
         ! allow(assumed-size)
         integer(hsize_t), intent(in) :: size(*)
         integer(herr_t) :: status
      end function H5Dset_extent

      function H5Dwrite(dset, mem_type, mem_space, file_space, dxpl, buf) result(status) &
         bind(C, name="H5Dwrite")
         import :: hid_t, herr_t, c_ptr
         implicit none
         integer(hid_t), value :: dset, mem_type, mem_space, file_space, dxpl
         type(c_ptr), value :: buf
         integer(herr_t) :: status
      end function H5Dwrite

      function H5Dread(dset, mem_type, mem_space, file_space, dxpl, buf) result(status) &
         bind(C, name="H5Dread")
         import :: hid_t, herr_t, c_ptr
         implicit none
         integer(hid_t), value :: dset, mem_type, mem_space, file_space, dxpl
         type(c_ptr), value :: buf
         integer(herr_t) :: status
      end function H5Dread

      ! -- property lists ---------------------------------------------------
      function H5Pcreate(class) result(plist) bind(C, name="H5Pcreate")
         import :: hid_t
         implicit none
         integer(hid_t), value :: class
         integer(hid_t) :: plist
      end function H5Pcreate

      function H5Pset_chunk(plist, ndims, dim) result(status) bind(C, name="H5Pset_chunk")
         import :: hid_t, herr_t, hsize_t, c_int
         implicit none
         integer(hid_t), value :: plist
         integer(c_int), value :: ndims
         ! allow(assumed-size)
         integer(hsize_t), intent(in) :: dim(*)
         integer(herr_t) :: status
      end function H5Pset_chunk

      function H5Pclose(plist) result(status) bind(C, name="H5Pclose")
         import :: hid_t, herr_t
         implicit none
         integer(hid_t), value :: plist
         integer(herr_t) :: status
      end function H5Pclose

      ! -- attributes -------------------------------------------------------
      function H5Acreate2(loc, name, type_id, space, acpl, aapl) result(attr) &
         bind(C, name="H5Acreate2")
         import :: hid_t, c_char
         implicit none
         integer(hid_t), value :: loc
         ! allow(assumed-size)
         character(kind=c_char), intent(in) :: name(*)
         integer(hid_t), value :: type_id, space, acpl, aapl
         integer(hid_t) :: attr
      end function H5Acreate2

      function H5Aopen(loc, name, aapl) result(attr) bind(C, name="H5Aopen")
         import :: hid_t, c_char
         implicit none
         integer(hid_t), value :: loc
         ! allow(assumed-size)
         character(kind=c_char), intent(in) :: name(*)
         integer(hid_t), value :: aapl
         integer(hid_t) :: attr
      end function H5Aopen

      function H5Awrite(attr, type_id, buf) result(status) bind(C, name="H5Awrite")
         import :: hid_t, herr_t, c_ptr
         implicit none
         integer(hid_t), value :: attr, type_id
         type(c_ptr), value :: buf
         integer(herr_t) :: status
      end function H5Awrite

      function H5Aread(attr, type_id, buf) result(status) bind(C, name="H5Aread")
         import :: hid_t, herr_t, c_ptr
         implicit none
         integer(hid_t), value :: attr, type_id
         type(c_ptr), value :: buf
         integer(herr_t) :: status
      end function H5Aread

      function H5Aclose(attr) result(status) bind(C, name="H5Aclose")
         import :: hid_t, herr_t
         implicit none
         integer(hid_t), value :: attr
         integer(herr_t) :: status
      end function H5Aclose

      function H5Aexists(loc, name) result(present) bind(C, name="H5Aexists")
         import :: hid_t, herr_t, c_char
         implicit none
         integer(hid_t), value :: loc
         ! allow(assumed-size)
         character(kind=c_char), intent(in) :: name(*)
         integer(herr_t) :: present   !! >0 yes, 0 no, <0 error
      end function H5Aexists

      ! -- datatypes --------------------------------------------------------
      function H5Tcopy(type_id) result(copy) bind(C, name="H5Tcopy")
         import :: hid_t
         implicit none
         integer(hid_t), value :: type_id
         integer(hid_t) :: copy
      end function H5Tcopy

      function H5Tset_size(type_id, size) result(status) bind(C, name="H5Tset_size")
         import :: hid_t, herr_t, c_size_t
         implicit none
         integer(hid_t), value :: type_id
         integer(c_size_t), value :: size
         integer(herr_t) :: status
      end function H5Tset_size

      function H5Tclose(type_id) result(status) bind(C, name="H5Tclose")
         import :: hid_t, herr_t
         implicit none
         integer(hid_t), value :: type_id
         integer(herr_t) :: status
      end function H5Tclose

      ! -- links ------------------------------------------------------------
      function H5Lexists(loc, name, lapl) result(present) bind(C, name="H5Lexists")
         import :: hid_t, herr_t, c_char
         implicit none
         integer(hid_t), value :: loc
         ! allow(assumed-size)
         character(kind=c_char), intent(in) :: name(*)
         integer(hid_t), value :: lapl
         integer(herr_t) :: present
      end function H5Lexists

   end interface

contains

   subroutine h5_version(major, minor, release)
      !! Which HDF5 is actually loaded
      integer, intent(out) :: major, minor, release

      integer(c_int) :: c_major, c_minor, c_release

      major = 0
      minor = 0
      release = 0
      if (H5get_libversion(c_major, c_minor, c_release) < 0) return
      major = int(c_major)
      minor = int(c_minor)
      release = int(c_release)
   end subroutine h5_version

   function h5_version_text() result(text)
      !! The loaded version, for a log line or an error message
      character(len=:), allocatable :: text
      character(len=32) :: buffer
      integer :: major, minor, release

      call h5_version(major, minor, release)
      write (buffer, "(i0,'.',i0,'.',i0)") major, minor, release
      text = trim(buffer)
   end function h5_version_text

   function h5_start() result(ok)
      !! Initialise the library, so the datatype globals are real ids
      !!
      !! Must run before any of `H5T_NATIVE_*` is read. The C macros and the
      !! Fortran interface each smuggle this in; a raw binding does not, and
      !! the failure -- a zero type id -- is accepted by the library as a
      !! different type rather than rejected as an error.
      !! It also refuses HDF5 older than 1.10, and that check is not
      !! decoration. `hid_t` was `int` until 1.10 and has been `int64_t`
      !! since; the declarations above assume the latter. Built or run
      !! against 1.8, every id would be assembled from the wrong eight bytes
      !! and the failure would be a wrong file rather than a bad status. A
      !! newer library is fine -- the C ABI has been stable since 1.10 for
      !! everything used here -- so this is a floor, not a pin.
      logical :: ok

      integer :: major, minor, release

      ok = (H5open() >= 0)
      if (.not. ok) return

      call h5_version(major, minor, release)
      ok = (major > 1) .or. (major == 1 .and. minor >= 10)
   end function h5_start

   pure function c_string(text) result(buffer)
      !! A null-terminated copy, for the name arguments above
      character(len=*), intent(in) :: text
      character(kind=c_char) :: buffer(len_trim(text) + 1)

      integer :: i

      do i = 1, len_trim(text)
         buffer(i) = text(i:i)
      end do
      buffer(len_trim(text) + 1) = c_null_char
   end function c_string

end module mqc_hdf5_bindings
