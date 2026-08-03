!! Low-level scaffolding shared by all cuEST-backed code in metalquicha
module mqc_cuest_runtime
   !! Wraps the raw cuEST and CUDA runtime calls in metalquicha conventions:
   !! statuses become `error_t` rather than `error stop`, and the three
   !! recurring chores (workspace allocation, device buffers, host transfers)
   !! get one implementation each.
   !!
   !! cuEST object creation follows a fixed **query -> allocate -> create**
   !! dance: the `...WorkspaceQuery` variant fills two workspace descriptors,
   !! scratch is allocated to those sizes, then the real call runs. The
   !! *persistent* workspace must outlive the object it created; the
   !! *temporary* one can be freed immediately afterwards.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_int64_t, &
                                                                             c_size_t, c_intptr_t, c_double, c_loc, &
                                                                             c_associated
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use cuest, only: cuestWorkspace_t, cuestWorkspaceDescriptor_t, CUEST_STATUS_SUCCESS
   use cuest_helpers, only: cuest_status_name
   use cuda_runtime, only: cudaMalloc, cudaFree, cudaMemcpy, cudaSuccess, &
                           cudaMemcpyHostToDevice, cudaMemcpyDeviceToHost, &
                           cudaDeviceSynchronize
   use cuda_helpers, only: cuda_error_name, cuda_error_string
   implicit none
   private

   public :: cuest_status_check   !! Turn a cuEST status into an error_t
   public :: cuda_status_check    !! Turn a CUDA status into an error_t
   public :: device_sync          !! cudaDeviceSynchronize with error_t
   public :: workspace_alloc, workspace_free  !! cuestWorkspace_t lifecycle
   public :: device_alloc, device_free        !! Device buffers counted in doubles
   public :: copy_to_device, copy_to_host     !! Host <-> device transfers

   integer(c_size_t), parameter :: BYTES_PER_DOUBLE = 8_c_size_t

   !> Plain C malloc/free for the host half of a cuEST workspace. The buffer
   !  must outlive the scoping unit that created it and is freed through a raw
   !  address, so a Fortran allocatable is not a substitute.
   interface
      type(c_ptr) function c_malloc(nbytes) bind(C, name="malloc")
         import :: c_ptr, c_size_t
         integer(c_size_t), value :: nbytes
      end function c_malloc

      subroutine c_free(ptr) bind(C, name="free")
         import :: c_ptr
         type(c_ptr), value :: ptr
      end subroutine c_free
   end interface

contains

   subroutine cuest_status_check(status, context, error)
      !! Record a failed cuEST status as an error, decoding the status name
      integer(c_int), intent(in) :: status      !! Status returned by a cuEST call
      character(len=*), intent(in) :: context   !! What was being attempted
      type(error_t), intent(inout) :: error

      if (status == CUEST_STATUS_SUCCESS) return
      if (error%has_error()) return  ! keep the first failure, not the last

      call error%set(ERROR_VALIDATION, "cuEST call failed: "//trim(context)// &
                     " -> "//trim(cuest_status_name(status)))
   end subroutine cuest_status_check

   subroutine cuda_status_check(status, context, error)
      !! Record a failed CUDA runtime status as an error
      integer(c_int), intent(in) :: status      !! Status returned by a CUDA call
      character(len=*), intent(in) :: context   !! What was being attempted
      type(error_t), intent(inout) :: error

      if (status == cudaSuccess) return
      if (error%has_error()) return

      call error%set(ERROR_VALIDATION, "CUDA call failed: "//trim(context)// &
                     " -> "//trim(cuda_error_name(status))//": "// &
                     trim(cuda_error_string(status)))
   end subroutine cuda_status_check

   subroutine device_sync(context, error)
      !! Block until all queued GPU work has completed
      character(len=*), intent(in) :: context
      type(error_t), intent(inout) :: error

      call cuda_status_check(cudaDeviceSynchronize(), &
                             "cudaDeviceSynchronize ("//trim(context)//")", error)
   end subroutine device_sync

   subroutine workspace_alloc(workspace, descriptor, error)
      !! Allocate host and device scratch sized by a workspace descriptor
      !!
      !! The `hostBuffer`/`deviceBuffer` fields are `uintptr_t` in C rather
      !! than pointers, so addresses are stored with TRANSFER.
      type(cuestWorkspace_t), intent(out) :: workspace
      type(cuestWorkspaceDescriptor_t), intent(in) :: descriptor
      type(error_t), intent(inout) :: error

      type(c_ptr) :: host_ptr, device_ptr

      workspace%hostBufferSizeInBytes = descriptor%hostBufferSizeInBytes
      workspace%deviceBufferSizeInBytes = descriptor%deviceBufferSizeInBytes
      workspace%hostBuffer = 0_c_intptr_t
      workspace%deviceBuffer = 0_c_intptr_t

      if (descriptor%hostBufferSizeInBytes > 0) then
         host_ptr = c_malloc(descriptor%hostBufferSizeInBytes)
         if (.not. c_associated(host_ptr)) then
            call error%set(ERROR_VALIDATION, "Failed to allocate cuEST host workspace")
            return
         end if
         workspace%hostBuffer = transfer(host_ptr, workspace%hostBuffer)
      end if

      if (descriptor%deviceBufferSizeInBytes > 0) then
         call cuda_status_check(cudaMalloc(device_ptr, descriptor%deviceBufferSizeInBytes), &
                                "cudaMalloc(cuEST device workspace)", error)
         if (error%has_error()) return
         workspace%deviceBuffer = transfer(device_ptr, workspace%deviceBuffer)
      end if
   end subroutine workspace_alloc

   subroutine workspace_free(workspace)
      !! Release the host and device scratch held by a workspace
      type(cuestWorkspace_t), intent(inout) :: workspace
      integer(c_int) :: status

      if (workspace%hostBuffer /= 0_c_intptr_t) then
         call c_free(transfer(workspace%hostBuffer, c_null_ptr))
      end if
      if (workspace%deviceBuffer /= 0_c_intptr_t) then
         status = cudaFree(transfer(workspace%deviceBuffer, c_null_ptr))
      end if

      workspace%hostBuffer = 0_c_intptr_t
      workspace%deviceBuffer = 0_c_intptr_t
      workspace%hostBufferSizeInBytes = 0_c_size_t
      workspace%deviceBufferSizeInBytes = 0_c_size_t
   end subroutine workspace_free

   subroutine device_alloc(device_ptr, n_doubles, context, error)
      !! Allocate `n_doubles` doubles of device memory
      type(c_ptr), intent(out) :: device_ptr
      integer(c_int64_t), intent(in) :: n_doubles
      character(len=*), intent(in) :: context
      type(error_t), intent(inout) :: error

      device_ptr = c_null_ptr
      if (n_doubles <= 0) return

      call cuda_status_check(cudaMalloc(device_ptr, &
                                        int(n_doubles, c_size_t)*BYTES_PER_DOUBLE), &
                             "cudaMalloc ("//trim(context)//")", error)
   end subroutine device_alloc

   subroutine device_free(device_ptr)
      !! Release a device buffer, tolerating a null pointer
      type(c_ptr), intent(inout) :: device_ptr
      integer(c_int) :: status

      if (c_associated(device_ptr)) status = cudaFree(device_ptr)
      device_ptr = c_null_ptr
   end subroutine device_free

   subroutine copy_to_device(device_ptr, host_array, context, error)
      !! Copy a host array of doubles onto the device
      type(c_ptr), intent(in) :: device_ptr
      real(dp), intent(in), target :: host_array(:)
      character(len=*), intent(in) :: context
      type(error_t), intent(inout) :: error

      if (size(host_array) == 0) return
      call cuda_status_check(cudaMemcpy(device_ptr, c_loc(host_array), &
                                        int(size(host_array), c_size_t)*BYTES_PER_DOUBLE, &
                                        cudaMemcpyHostToDevice), &
                             "cudaMemcpy H2D ("//trim(context)//")", error)
   end subroutine copy_to_device

   subroutine copy_to_host(host_array, device_ptr, context, error)
      !! Copy a device buffer of doubles back into a host array
      real(dp), intent(out), target :: host_array(:)
      type(c_ptr), intent(in) :: device_ptr
      character(len=*), intent(in) :: context
      type(error_t), intent(inout) :: error

      if (size(host_array) == 0) return
      call cuda_status_check(cudaMemcpy(c_loc(host_array), device_ptr, &
                                        int(size(host_array), c_size_t)*BYTES_PER_DOUBLE, &
                                        cudaMemcpyDeviceToHost), &
                             "cudaMemcpy D2H ("//trim(context)//")", error)
   end subroutine copy_to_host

end module mqc_cuest_runtime
