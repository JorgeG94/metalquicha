!! Manual check that the hand-written HDF5 bindings actually work
!!
!!     cmake -B build -DMQC_ENABLE_HDF5=ON && ./build/check_hdf5
!!     h5dump check_hdf5.h5      # readable by the real tool, not just by us
!!
!! Bindings are the one place where being wrong is silent: an integer kind
!! that does not match the C type is not a compile error, it is a value the
!! library reads out of the wrong bytes. So this does not check return codes
!! and stop -- it writes an extendible dataset, grows it, writes again, reads
!! everything back and compares. If `hid_t` or `hsize_t` were wrong, the
!! numbers would not survive that.
!!
!! It also exercises the two things the checkpoint depends on and nothing else
!! would catch: growing a dataset after the file has been flushed, and a
!! string attribute round-tripping.
program check_hdf5
   use, intrinsic :: iso_c_binding, only: c_loc, c_null_ptr, c_size_t, c_char
   use pic_types, only: dp
   use mqc_hdf5_bindings, only: hid_t, hsize_t, h5_start, c_string, h5_version_text, &
                                H5F_ACC_TRUNC, H5F_ACC_RDONLY, H5P_DEFAULT, H5S_ALL, &
                                H5S_UNLIMITED, H5S_SELECT_SET, H5F_SCOPE_LOCAL, &
                                H5T_NATIVE_DOUBLE, H5T_C_S1, H5P_CLS_DATASET_CREATE_ID, &
                                H5Fcreate, H5Fopen, H5Fclose, H5Fflush, &
                                H5Screate, H5Screate_simple, H5Sclose, H5Sselect_hyperslab, &
                                H5Dcreate2, H5Dopen2, H5Dclose, H5Dwrite, H5Dread, &
                                H5Dset_extent, H5Dget_space, &
                                H5Pcreate, H5Pset_chunk, H5Pclose, &
                                H5Acreate2, H5Aopen, H5Awrite, H5Aread, H5Aclose, &
                                H5Tcopy, H5Tset_size, H5Tclose
   implicit none

   character(len=*), parameter :: PATH = "check_hdf5.h5"
   integer, parameter :: CHUNK = 4
   integer, parameter :: N_FIRST = 6      !! values written before the flush
   integer, parameter :: N_SECOND = 4     !! values written after it
   integer, parameter :: N_TOTAL = N_FIRST + N_SECOND
   integer, parameter :: FP_LEN = 17      !! 16 hex digits and a terminator

   integer(hid_t) :: file, space, dset, dcpl, mem_space, attr, atype, ascalar
   integer(hsize_t) :: dims(1), maxdims(1), chunk_dims(1), start(1), count(1)
   real(dp), target :: values(N_FIRST), readback(N_TOTAL)
   real(dp), target :: more(N_SECOND)
   character(kind=c_char), target :: fingerprint(FP_LEN)
   character(kind=c_char), target :: fp_read(FP_LEN)
   integer :: i, failures
   logical :: ok

   failures = 0

   ! Without this the datatype globals below are zero, which the library
   ! accepts as a different type rather than rejecting.
   ok = h5_start()
   call expect(ok, "H5open, and HDF5 at least 1.10")
   write (*, "(A)") "[hdf5] linked against "//h5_version_text()// &
      " (hid_t is int64_t from 1.10 on; these bindings assume that)"
   if (.not. ok) error stop 1

   ! ---- create, with an unlimited first dimension ---------------------------
   file = H5Fcreate(c_string(PATH), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)
   call expect(file >= 0, "H5Fcreate")

   dims = [0_hsize_t]
   maxdims = [H5S_UNLIMITED]
   space = H5Screate_simple(1, dims, maxdims)
   call expect(space >= 0, "H5Screate_simple")

   ! Chunking is what makes a dataset extendible at all, and the chunk size is
   ! also the granularity a killed writer loses.
   chunk_dims = [int(CHUNK, hsize_t)]
   dcpl = H5Pcreate(H5P_CLS_DATASET_CREATE_ID)
   call expect(dcpl >= 0, "H5Pcreate")
   call expect(H5Pset_chunk(dcpl, 1, chunk_dims) >= 0, "H5Pset_chunk")

   dset = H5Dcreate2(file, c_string("energies"), H5T_NATIVE_DOUBLE, space, &
                     H5P_DEFAULT, dcpl, H5P_DEFAULT)
   call expect(dset >= 0, "H5Dcreate2")

   ! ---- first append: six values --------------------------------------------
   do i = 1, N_FIRST
      values(i) = -1.0_dp*i - 0.25_dp
   end do
   call append(dset, values, 0_hsize_t, int(N_FIRST, hsize_t))

   ! Flushed mid-stream, then grown again. This is exactly what the checkpoint
   ! does and the case where HDF5 is least forgiving.
   call expect(H5Fflush(file, H5F_SCOPE_LOCAL) >= 0, "H5Fflush")

   ! ---- second append, after the flush --------------------------------------
   do i = 1, N_SECOND
      more(i) = -100.0_dp*i
   end do
   call append(dset, more, int(N_FIRST, hsize_t), int(N_SECOND, hsize_t))

   ! ---- a string attribute, as the fingerprint will be --------------------
   fingerprint = c_string("f5f16b57fa5252ab")
   atype = H5Tcopy(H5T_C_S1)
   call expect(H5Tset_size(atype, int(FP_LEN, c_size_t)) >= 0, "H5Tset_size")
   ascalar = H5Screate(0)   ! H5S_SCALAR
   attr = H5Acreate2(file, c_string("fingerprint"), atype, ascalar, &
                     H5P_DEFAULT, H5P_DEFAULT)
   call expect(attr >= 0, "H5Acreate2")
   call expect(H5Awrite(attr, atype, c_loc(fingerprint)) >= 0, "H5Awrite")
   call expect(H5Aclose(attr) >= 0, "H5Aclose")

   call expect(H5Dclose(dset) >= 0, "H5Dclose")
   call expect(H5Sclose(space) >= 0, "H5Sclose")
   call expect(H5Pclose(dcpl) >= 0, "H5Pclose")
   call expect(H5Fclose(file) >= 0, "H5Fclose")

   ! ---- reopen and compare --------------------------------------------------
   file = H5Fopen(c_string(PATH), H5F_ACC_RDONLY, H5P_DEFAULT)
   call expect(file >= 0, "H5Fopen")

   dset = H5Dopen2(file, c_string("energies"), H5P_DEFAULT)
   call expect(dset >= 0, "H5Dopen2")
   readback = 0.0_dp
   call expect(H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &
                       c_loc(readback)) >= 0, "H5Dread")

   do i = 1, N_FIRST
      call expect(abs(readback(i) - (-1.0_dp*i - 0.25_dp)) < 1.0e-15_dp, &
                  "value before the flush survived")
   end do
   do i = 1, N_SECOND
      call expect(abs(readback(N_FIRST + i) - (-100.0_dp*i)) < 1.0e-15_dp, &
                  "value written after the flush survived")
   end do

   attr = H5Aopen(file, c_string("fingerprint"), H5P_DEFAULT)
   call expect(attr >= 0, "H5Aopen")
   call expect(H5Aread(attr, atype, c_loc(fp_read)) >= 0, "H5Aread")
   call expect(all(fp_read(1:FP_LEN - 1) == fingerprint(1:FP_LEN - 1)), "fingerprint round-tripped")

   call expect(H5Aclose(attr) >= 0, "H5Aclose (read)")
   call expect(H5Tclose(atype) >= 0, "H5Tclose")
   call expect(H5Sclose(ascalar) >= 0, "H5Sclose (scalar)")
   call expect(H5Dclose(dset) >= 0, "H5Dclose (read)")
   call expect(H5Fclose(file) >= 0, "H5Fclose (read)")

   write (*, "(A)") ""
   if (failures == 0) then
      write (*, "(A)") "[hdf5] all ok -- bindings agree with the C ABI"
      write (*, "(A)") "       check "//PATH//" with h5dump to confirm the file is"
      write (*, "(A)") "       readable by tools that did not write it."
   else
      write (*, "(A,I0,A)") "[hdf5] ", failures, " FAILURE(S)"
      error stop 1
   end if

contains

   subroutine append(dataset, buffer, offset, n)
      !! Grow the dataset and write `n` values at `offset`
      integer(hid_t), intent(in) :: dataset
      real(dp), intent(in), target :: buffer(:)
      integer(hsize_t), intent(in) :: offset, n

      integer(hid_t) :: fspace, mspace
      integer(hsize_t) :: new_dims(1), sel_start(1), sel_count(1)

      new_dims = [offset + n]
      call expect(H5Dset_extent(dataset, new_dims) >= 0, "H5Dset_extent")

      fspace = H5Dget_space(dataset)
      sel_start = [offset]
      sel_count = [n]
      call expect(H5Sselect_hyperslab(fspace, H5S_SELECT_SET, sel_start, c_null_ptr, &
                                      sel_count, c_null_ptr) >= 0, "H5Sselect_hyperslab")
      mspace = H5Screate_simple(1, sel_count, sel_count)
      call expect(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, &
                           c_loc(buffer)) >= 0, "H5Dwrite")
      call expect(H5Sclose(mspace) >= 0, "H5Sclose (mem)")
      call expect(H5Sclose(fspace) >= 0, "H5Sclose (file)")
   end subroutine append

   subroutine expect(condition, what)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: what

      if (.not. condition) then
         write (*, "(A)") "  FAILED: "//what
         failures = failures + 1
      end if
   end subroutine expect

end program check_hdf5
