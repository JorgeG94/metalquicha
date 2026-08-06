! ============================================================================
!  cublas.f90 -- Fortran 2008 iso_c_binding interface to NVIDIA cuBLAS
!
!  Hand-written, covering only the calls the SCF needs. cuBLAS has hundreds of
!  entry points; binding the six used here keeps the surface reviewable and the
!  file short enough to extend by hand when something else is needed.
!
!  Why the C API rather than CUDA Fortran
!  -------------------------------------
!  `cudafor` would tie every file that touches the device to nvfortran. That is
!  a hard constraint here, not a preference: toml-f -- reached through tblite --
!  does not compile with nvfortran (a backslash literal it rejects, with no flag
!  to disable), so an nvfortran-only device path would force a choice between
!  the xTB backend and the GPU backend in any single build. Binding the C ABI
!  through iso_c_binding compiles under gfortran, ifx and nvfortran alike, and
!  the device memory it operates on already comes from the CUDA runtime
!  bindings in cuda_runtime.f90.
!
!  Conventions
!  -----------
!  * Every function returns cublasStatus_t, exposed as INTEGER(c_int).
!  * The handle is C  cublasHandle_t  (an opaque pointer):
!      - created  -> TYPE(c_ptr), INTENT(OUT)   (C cublasHandle_t*)
!      - passed in-> TYPE(c_ptr), VALUE
!  * Vector and matrix arguments are DEVICE pointers, TYPE(c_ptr), VALUE.
!    Passing a host address silently produces wrong results or a launch
!    failure rather than a clean status, exactly as with cuEST.
!  * The _v2 names are the current ABI. The unsuffixed names are legacy macros
!    in the C header and are not exported symbols, so bind to _v2 explicitly.
!  * cuBLAS is COLUMN-major, which is Fortran's own layout -- no transpose
!    dance, unlike the row-major convention cuEST uses (see the note at the top
!    of mqc_cuest_integrals.f90).
!  * Scalars (alpha, beta) and scalar RESULTS (dot, nrm2) default to HOST
!    pointers under CUBLAS_POINTER_MODE_HOST. They are bound here as ordinary
!    Fortran references, which is what that mode expects. Switching the handle
!    to CUBLAS_POINTER_MODE_DEVICE would make these device addresses and break
!    every call below.
!
!  Build:  gfortran -c cublas.f90     (produces cublas.mod)
!  Link :  ... -lcublas
! ============================================================================
module cublas
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double
   implicit none
   public

   ! ---- cublasStatus_t --------------------------------------------------
   integer(c_int), parameter :: CUBLAS_STATUS_SUCCESS = 0
   integer(c_int), parameter :: CUBLAS_STATUS_NOT_INITIALIZED = 1
   integer(c_int), parameter :: CUBLAS_STATUS_ALLOC_FAILED = 3
   integer(c_int), parameter :: CUBLAS_STATUS_INVALID_VALUE = 7
   integer(c_int), parameter :: CUBLAS_STATUS_ARCH_MISMATCH = 8
   integer(c_int), parameter :: CUBLAS_STATUS_MAPPING_ERROR = 11
   integer(c_int), parameter :: CUBLAS_STATUS_EXECUTION_FAILED = 13
   integer(c_int), parameter :: CUBLAS_STATUS_INTERNAL_ERROR = 14
   integer(c_int), parameter :: CUBLAS_STATUS_NOT_SUPPORTED = 15
   integer(c_int), parameter :: CUBLAS_STATUS_LICENSE_ERROR = 16

   ! ---- cublasOperation_t -----------------------------------------------
   integer(c_int), parameter :: CUBLAS_OP_N = 0  !! no transpose
   integer(c_int), parameter :: CUBLAS_OP_T = 1  !! transpose
   integer(c_int), parameter :: CUBLAS_OP_C = 2  !! conjugate transpose

   ! ---- cublasPointerMode_t ---------------------------------------------
   integer(c_int), parameter :: CUBLAS_POINTER_MODE_HOST = 0
   integer(c_int), parameter :: CUBLAS_POINTER_MODE_DEVICE = 1

   interface

      ! ---- handle lifetime ----------------------------------------------

      integer(c_int) function cublasCreate(handle) &
         bind(C, name="cublasCreate_v2")
         import :: c_ptr, c_int
         type(c_ptr), intent(out) :: handle  !! cublasHandle_t*
      end function cublasCreate

      integer(c_int) function cublasDestroy(handle) &
         bind(C, name="cublasDestroy_v2")
         import :: c_ptr, c_int
         type(c_ptr), value :: handle
      end function cublasDestroy

      integer(c_int) function cublasSetStream(handle, streamId) &
         bind(C, name="cublasSetStream_v2")
         !! Put cuBLAS on the same stream as the rest of the work, so its
         !! calls order against them without a device-wide synchronize.
         import :: c_ptr, c_int
         type(c_ptr), value :: handle
         type(c_ptr), value :: streamId  !! cudaStream_t; null means the default stream
      end function cublasSetStream

      integer(c_int) function cublasSetPointerMode(handle, mode) &
         bind(C, name="cublasSetPointerMode_v2")
         !! Whether scalar arguments and results live on the host or the device.
         !! The bindings below assume HOST, which is the default.
         import :: c_ptr, c_int
         type(c_ptr), value :: handle
         integer(c_int), value :: mode
      end function cublasSetPointerMode

      ! ---- level 1 -------------------------------------------------------

      integer(c_int) function cublasDaxpy(handle, n, alpha, x, incx, y, incy) &
         bind(C, name="cublasDaxpy_v2")
         !! y := alpha*x + y
         !!
         !! The Fock extrapolation and the h + J - K + Vxc assembly are both
         !! this, over flattened matrices.
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         integer(c_int), value :: n
         real(c_double), intent(in) :: alpha   !! host scalar under POINTER_MODE_HOST
         type(c_ptr), value :: x               !! device
         integer(c_int), value :: incx
         type(c_ptr), value :: y               !! device, updated in place
         integer(c_int), value :: incy
      end function cublasDaxpy

      integer(c_int) function cublasDcopy(handle, n, x, incx, y, incy) &
         bind(C, name="cublasDcopy_v2")
         !! y := x, device to device
         !!
         !! Pushing a Fock or error matrix into the DIIS history, without a
         !! round trip through the host.
         import :: c_ptr, c_int
         type(c_ptr), value :: handle
         integer(c_int), value :: n
         type(c_ptr), value :: x               !! device
         integer(c_int), value :: incx
         type(c_ptr), value :: y               !! device
         integer(c_int), value :: incy
      end function cublasDcopy

      integer(c_int) function cublasDscal(handle, n, alpha, x, incx) &
         bind(C, name="cublasDscal_v2")
         !! x := alpha*x
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         integer(c_int), value :: n
         real(c_double), intent(in) :: alpha   !! host scalar
         type(c_ptr), value :: x               !! device
         integer(c_int), value :: incx
      end function cublasDscal

      integer(c_int) function cublasDdot(handle, n, x, incx, y, incy, result) &
         bind(C, name="cublasDdot_v2")
         !! result := x . y
         !!
         !! One DIIS overlap entry. Note that the result crossing to the host
         !! forces a synchronize: for a subspace of any size prefer cublasDgemv
         !! of the new error vector against the whole history, one call instead
         !! of n_stored round trips.
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         integer(c_int), value :: n
         type(c_ptr), value :: x               !! device
         integer(c_int), value :: incx
         type(c_ptr), value :: y               !! device
         integer(c_int), value :: incy
         real(c_double), intent(out) :: result !! host scalar under POINTER_MODE_HOST
      end function cublasDdot

      integer(c_int) function cublasDnrm2(handle, n, x, incx, result) &
         bind(C, name="cublasDnrm2_v2")
         !! result := ||x||_2 -- the DIIS convergence measure
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         integer(c_int), value :: n
         type(c_ptr), value :: x               !! device
         integer(c_int), value :: incx
         real(c_double), intent(out) :: result !! host scalar
      end function cublasDnrm2

      ! ---- level 2 -------------------------------------------------------

      integer(c_int) function cublasDgemv(handle, trans, m, n, alpha, A, lda, &
                                          x, incx, beta, y, incy) &
         bind(C, name="cublasDgemv_v2")
         !! y := alpha*op(A)*x + beta*y
         !!
         !! The batched form of the overlap update: A is the error history as
         !! (n_error, n_stored), x the new error vector, and y the whole new
         !! row of the overlap matrix in one call.
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         integer(c_int), value :: trans
         integer(c_int), value :: m, n
         real(c_double), intent(in) :: alpha
         type(c_ptr), value :: A               !! device
         integer(c_int), value :: lda
         type(c_ptr), value :: x               !! device
         integer(c_int), value :: incx
         real(c_double), intent(in) :: beta
         type(c_ptr), value :: y               !! device
         integer(c_int), value :: incy
      end function cublasDgemv

      ! ---- level 3 -------------------------------------------------------

      integer(c_int) function cublasDgemm(handle, transa, transb, m, n, k, &
                                          alpha, A, lda, B, ldb, beta, C, ldc) &
         bind(C, name="cublasDgemm_v2")
         !! C := alpha*op(A)*op(B) + beta*C
         !!
         !! The commutator FDS - SDF and its projection into the orthogonal
         !! basis: four gemms for the commutator, two for the transform.
         !! Column-major, so the leading dimensions are the Fortran ones and no
         !! transposition is needed to match Fortran storage.
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: handle
         integer(c_int), value :: transa, transb
         integer(c_int), value :: m, n, k
         real(c_double), intent(in) :: alpha
         type(c_ptr), value :: A               !! device
         integer(c_int), value :: lda
         type(c_ptr), value :: B               !! device
         integer(c_int), value :: ldb
         real(c_double), intent(in) :: beta
         type(c_ptr), value :: C               !! device
         integer(c_int), value :: ldc
      end function cublasDgemm

   end interface

contains

   pure function cublas_status_name(status) result(name)
      !! Decode a cuBLAS status for an error message
      integer(c_int), intent(in) :: status
      character(len=:), allocatable :: name

      select case (status)
      case (CUBLAS_STATUS_SUCCESS); name = "CUBLAS_STATUS_SUCCESS"
      case (CUBLAS_STATUS_NOT_INITIALIZED); name = "CUBLAS_STATUS_NOT_INITIALIZED"
      case (CUBLAS_STATUS_ALLOC_FAILED); name = "CUBLAS_STATUS_ALLOC_FAILED"
      case (CUBLAS_STATUS_INVALID_VALUE); name = "CUBLAS_STATUS_INVALID_VALUE"
      case (CUBLAS_STATUS_ARCH_MISMATCH); name = "CUBLAS_STATUS_ARCH_MISMATCH"
      case (CUBLAS_STATUS_MAPPING_ERROR); name = "CUBLAS_STATUS_MAPPING_ERROR"
      case (CUBLAS_STATUS_EXECUTION_FAILED); name = "CUBLAS_STATUS_EXECUTION_FAILED"
      case (CUBLAS_STATUS_INTERNAL_ERROR); name = "CUBLAS_STATUS_INTERNAL_ERROR"
      case (CUBLAS_STATUS_NOT_SUPPORTED); name = "CUBLAS_STATUS_NOT_SUPPORTED"
      case (CUBLAS_STATUS_LICENSE_ERROR); name = "CUBLAS_STATUS_LICENSE_ERROR"
      case default; name = "CUBLAS_STATUS_UNKNOWN"
      end select
   end function cublas_status_name

end module cublas
