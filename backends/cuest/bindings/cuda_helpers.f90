! ============================================================================
!  cuda_helpers.f90 -- optional convenience layer over the generated
!  `cuda_runtime` module.
!
!  Hand-written (the generated modules are pure, faithful interfaces and stay
!  that way). Everything here is a thin wrapper that forwards to the raw call
!  and returns the cudaError_t it produced; nothing here hides an error.
!
!  What it provides:
!    * cuda_error_string / cuda_error_name -- C string -> Fortran character(:)
!    * cuda_check                          -- abort with a decoded message
!    * cuda_malloc / cuda_free             -- element-counted device allocation
!    * cuda_memcpy_*                       -- typed Fortran-array copies, so
!                                             callers need no C_LOC/C_SIZEOF
!    * cuda_device_name / cuda_capability  -- readable device properties
! ============================================================================
module cuda_helpers
   use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_int, c_int32_t, c_int64_t, &
                                                                             c_size_t, c_float, c_double, &
                                                                             c_loc, c_f_pointer, c_associated, &
                                                                             c_null_ptr, c_null_char
   use cuda_runtime, only: cudaGetErrorString, cudaGetErrorName, cudaSuccess, &
                           cudaMalloc, cudaFree, cudaMemcpy, &
                           cudaMemcpyHostToDevice, cudaMemcpyDeviceToHost, &
                           cudaDeviceProp, cudaGetDeviceCount, cudaGetDeviceProperties
   use pic_logger, only: logger => global_logger
   use mqc_program_limits, only: MAX_LINE_LENGTH
   implicit none
   private

   public :: cuda_error_string, cuda_error_name, cuda_check
   public :: cuda_malloc, cuda_free
   public :: cuda_memcpy_to_device, cuda_memcpy_to_host
   public :: cuda_device_name, cuda_capability, cuda_report_devices

   !> Copy a Fortran host array to a device buffer.
   interface cuda_memcpy_to_device
      module procedure memcpy_h2d_r64_1, memcpy_h2d_r64_2
      module procedure memcpy_h2d_r32_1, memcpy_h2d_i32_1
   end interface

   !> Copy a device buffer back into a Fortran host array.
   interface cuda_memcpy_to_host
      module procedure memcpy_d2h_r64_1, memcpy_d2h_r64_2
      module procedure memcpy_d2h_r32_1, memcpy_d2h_i32_1
   end interface

contains

   ! ---- error reporting --------------------------------------------------

   !> Copy a NUL-terminated C string into a Fortran allocatable string.
   function c_string(p) result(s)
      type(c_ptr), intent(in) :: p
      character(len=:), allocatable :: s
      character(kind=c_char), pointer :: buf(:)
      integer :: n
      if (.not. c_associated(p)) then
         s = "<null>"
         return
      end if
      call c_f_pointer(p, buf, [huge(0_c_int)])
      n = 0
      do while (buf(n + 1) /= c_null_char)
         n = n + 1
         if (n > 4096) exit
      end do
      allocate (character(len=n) :: s)
      block
         integer :: i
         do i = 1, n
            s(i:i) = buf(i)
         end do
      end block
   end function c_string

   function cuda_error_string(code) result(s)
      integer(c_int), intent(in) :: code
      character(len=:), allocatable :: s
      s = c_string(cudaGetErrorString(code))
   end function cuda_error_string

   function cuda_error_name(code) result(s)
      integer(c_int), intent(in) :: code
      character(len=:), allocatable :: s
      s = c_string(cudaGetErrorName(code))
   end function cuda_error_name

   !> Abort with a decoded message if a CUDA call did not return cudaSuccess.
   subroutine cuda_check(code, what)
      integer(c_int), intent(in) :: code
      character(*), intent(in) :: what
      character(len=MAX_LINE_LENGTH) :: line

      if (code /= cudaSuccess) then
         write (line, "(A,A,A,I0,A,A,A,A)") "CUDA FAILED: ", what, " (code ", code, ", ", &
            cuda_error_name(code), "): ", cuda_error_string(code)
         call logger%info(trim(line))
         error stop 1
      end if
   end subroutine cuda_check

   ! ---- allocation -------------------------------------------------------

   !> Allocate `nbytes` of device memory.
   integer(c_int) function cuda_malloc(dptr, nbytes) result(ist)
      type(c_ptr), intent(out) :: dptr
      integer(c_int64_t), intent(in)  :: nbytes
      ist = cudaMalloc(dptr, int(nbytes, c_size_t))
   end function cuda_malloc

   integer(c_int) function cuda_free(dptr) result(ist)
      type(c_ptr), intent(inout) :: dptr
      ist = cudaSuccess
      if (c_associated(dptr)) ist = cudaFree(dptr)
      dptr = c_null_ptr
   end function cuda_free

   ! ---- typed host <-> device copies -------------------------------------
   ! The arrays are declared TARGET here so C_LOC is legal regardless of how
   ! the caller declared them.

   integer(c_int) function memcpy_h2d_r64_1(dst, src) result(ist)
      type(c_ptr), intent(in) :: dst
      real(c_double), intent(in), target :: src(:)
      ist = cudaMemcpy(dst, c_loc(src), &
                       int(size(src), c_size_t)*8_c_size_t, &
                       cudaMemcpyHostToDevice)
   end function memcpy_h2d_r64_1

   integer(c_int) function memcpy_h2d_r64_2(dst, src) result(ist)
      type(c_ptr), intent(in) :: dst
      real(c_double), intent(in), target :: src(:, :)
      ist = cudaMemcpy(dst, c_loc(src), &
                       int(size(src), c_size_t)*8_c_size_t, &
                       cudaMemcpyHostToDevice)
   end function memcpy_h2d_r64_2

   integer(c_int) function memcpy_h2d_r32_1(dst, src) result(ist)
      type(c_ptr), intent(in) :: dst
      real(c_float), intent(in), target :: src(:)
      ist = cudaMemcpy(dst, c_loc(src), &
                       int(size(src), c_size_t)*4_c_size_t, &
                       cudaMemcpyHostToDevice)
   end function memcpy_h2d_r32_1

   integer(c_int) function memcpy_h2d_i32_1(dst, src) result(ist)
      type(c_ptr), intent(in) :: dst
      integer(c_int32_t), intent(in), target :: src(:)
      ist = cudaMemcpy(dst, c_loc(src), &
                       int(size(src), c_size_t)*4_c_size_t, &
                       cudaMemcpyHostToDevice)
   end function memcpy_h2d_i32_1

   integer(c_int) function memcpy_d2h_r64_1(dst, src) result(ist)
      real(c_double), intent(out), target :: dst(:)
      type(c_ptr), intent(in) :: src
      ist = cudaMemcpy(c_loc(dst), src, &
                       int(size(dst), c_size_t)*8_c_size_t, &
                       cudaMemcpyDeviceToHost)
   end function memcpy_d2h_r64_1

   integer(c_int) function memcpy_d2h_r64_2(dst, src) result(ist)
      real(c_double), intent(out), target :: dst(:, :)
      type(c_ptr), intent(in) :: src
      ist = cudaMemcpy(c_loc(dst), src, &
                       int(size(dst), c_size_t)*8_c_size_t, &
                       cudaMemcpyDeviceToHost)
   end function memcpy_d2h_r64_2

   integer(c_int) function memcpy_d2h_r32_1(dst, src) result(ist)
      real(c_float), intent(out), target :: dst(:)
      type(c_ptr), intent(in) :: src
      ist = cudaMemcpy(c_loc(dst), src, &
                       int(size(dst), c_size_t)*4_c_size_t, &
                       cudaMemcpyDeviceToHost)
   end function memcpy_d2h_r32_1

   integer(c_int) function memcpy_d2h_i32_1(dst, src) result(ist)
      integer(c_int32_t), intent(out), target :: dst(:)
      type(c_ptr), intent(in) :: src
      ist = cudaMemcpy(c_loc(dst), src, &
                       int(size(dst), c_size_t)*4_c_size_t, &
                       cudaMemcpyDeviceToHost)
   end function memcpy_d2h_i32_1

   ! ---- device properties ------------------------------------------------

   !> The device name as a Fortran string (cudaDeviceProp%name is char(256)).
   function cuda_device_name(prop) result(s)
      type(cudaDeviceProp), intent(in) :: prop
      character(len=:), allocatable :: s
      integer :: i, n
      n = 0
      do i = 1, size(prop%name)
         if (prop%name(i) == c_null_char) exit
         n = i
      end do
      allocate (character(len=n) :: s)
      do i = 1, n
         s(i:i) = prop%name(i)
      end do
   end function cuda_device_name

   !> Compute capability as the familiar integer, e.g. 70 for Volta, 90 for
   !  Hopper -- the form CMAKE_CUDA_ARCHITECTURES and cuEST both use.
   integer function cuda_capability(prop) result(sm)
      type(cudaDeviceProp), intent(in) :: prop
      sm = 10*prop%major + prop%minor
   end function cuda_capability

   !> Print a one-line summary of every visible device.
   subroutine cuda_report_devices()
      integer(c_int) :: ist, ndev, i
      type(cudaDeviceProp) :: prop
      character(len=MAX_LINE_LENGTH) :: line

      ist = cudaGetDeviceCount(ndev)
      call cuda_check(ist, "cudaGetDeviceCount")
      write (line, "(A,I0,A)") "CUDA devices visible: ", ndev
      call logger%info(trim(line))
      do i = 0, ndev - 1
         ist = cudaGetDeviceProperties(prop, i)
         call cuda_check(ist, "cudaGetDeviceProperties")
         write (line, "(2X,A,I0,A,A,A,I0,A,F6.1,A)") "[", i, "] ", cuda_device_name(prop), "  sm_", &
            cuda_capability(prop), "  ", real(prop%totalGlobalMem)/1024.0_c_double**3, " GiB"
         call logger%info(trim(line))
      end do
   end subroutine cuda_report_devices

end module cuda_helpers
