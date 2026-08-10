! ============================================================================
!  cusolver.f90 -- Fortran 2008 iso_c_binding interface to NVIDIA cuSOLVER
!
!  Hand-written, and deliberately tiny: the SCF needs exactly one dense
!  eigensolver, so this binds `cusolverDnDsyevd` and the four calls around it
!  and nothing else. cuSOLVER's surface is enormous; binding it wholesale
!  would be unreviewable and none of the rest would ever be called.
!
!  Why Dsyevd and not Dsygvd
!  -------------------------
!  The SCF solves FC = SCe, which looks like the generalized problem
!  `cusolverDnDsygvd` exists for. It is not solved that way here. The overlap
!  is canonically orthogonalized once, up front, into X = U s^(-1/2) with the
!  near-null modes dropped -- which is what keeps a basis with near-linear
!  dependence from producing garbage, and which a generalized solver would not
!  do for us. After that transform the problem is the ordinary symmetric one,
!  F' = X^T F X, and Dsyevd is the right and cheaper call.
!
!  Conventions -- the same ones as cublas.f90, which see for the reasoning
!  -----------------------------------------------------------------------
!  * Every function returns cusolverStatus_t, exposed as INTEGER(c_int).
!  * The handle is C  cusolverDnHandle_t  (an opaque pointer):
!      - created  -> TYPE(c_ptr), INTENT(OUT)   (C cusolverDnHandle_t*)
!      - passed in-> TYPE(c_ptr), VALUE
!  * Matrix and vector arguments are DEVICE pointers, TYPE(c_ptr), VALUE.
!  * cuSOLVER is COLUMN-major, Fortran's own layout, so leading dimensions are
!    the Fortran ones.
!  * `lwork` comes back to the HOST; `info` is written to the DEVICE. That
!    asymmetry is not a mistake in the transcription -- it is how the C API is,
!    and mixing the two up gives a segfault rather than a status.
!  * Eigenvalues come back ascending, matching LAPACK's dsyev and therefore
!    `pic_syev`, so the orbital ordering the SCF assumes is unchanged.
!  * Eigenvectors are determined only up to sign, and cuSOLVER need not pick
!    the same signs as LAPACK. Nothing downstream cares: the density
!    D = 2 C_occ C_occ^T, the exchange K[C_occ] and the energy-weighted density
!    are all quadratic in C.
!
!  Build:  gfortran -c cusolver.f90   (produces cusolver.mod)
!  Link :  ... -lcusolver
! ============================================================================
module cusolver
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int
   implicit none
   private

   public :: cusolver_status_name
   public :: cusolverDnCreate, cusolverDnDestroy, cusolverDnSetStream
   public :: cusolverDnDsyevd_bufferSize, cusolverDnDsyevd
   public :: CUSOLVER_STATUS_SUCCESS, CUSOLVER_STATUS_NOT_INITIALIZED
   public :: CUSOLVER_STATUS_ALLOC_FAILED, CUSOLVER_STATUS_INVALID_VALUE
   public :: CUSOLVER_STATUS_ARCH_MISMATCH, CUSOLVER_STATUS_MAPPING_ERROR
   public :: CUSOLVER_STATUS_EXECUTION_FAILED, CUSOLVER_STATUS_INTERNAL_ERROR
   public :: CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED, CUSOLVER_STATUS_NOT_SUPPORTED
   public :: CUSOLVER_STATUS_ZERO_PIVOT, CUSOLVER_STATUS_INVALID_LICENSE
   public :: CUSOLVER_STATUS_INVALID_WORKSPACE
   public :: CUSOLVER_EIG_MODE_NOVECTOR, CUSOLVER_EIG_MODE_VECTOR

   ! ---- cusolverStatus_t -------------------------------------------------
   integer(c_int), parameter :: CUSOLVER_STATUS_SUCCESS = 0
   integer(c_int), parameter :: CUSOLVER_STATUS_NOT_INITIALIZED = 1
   integer(c_int), parameter :: CUSOLVER_STATUS_ALLOC_FAILED = 2
   integer(c_int), parameter :: CUSOLVER_STATUS_INVALID_VALUE = 3
   integer(c_int), parameter :: CUSOLVER_STATUS_ARCH_MISMATCH = 4
   integer(c_int), parameter :: CUSOLVER_STATUS_MAPPING_ERROR = 5
   integer(c_int), parameter :: CUSOLVER_STATUS_EXECUTION_FAILED = 6
   integer(c_int), parameter :: CUSOLVER_STATUS_INTERNAL_ERROR = 7
   integer(c_int), parameter :: CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED = 8
   integer(c_int), parameter :: CUSOLVER_STATUS_NOT_SUPPORTED = 9
   integer(c_int), parameter :: CUSOLVER_STATUS_ZERO_PIVOT = 10
   integer(c_int), parameter :: CUSOLVER_STATUS_INVALID_LICENSE = 11
   integer(c_int), parameter :: CUSOLVER_STATUS_INVALID_WORKSPACE = 31

   ! ---- cusolverEigMode_t ------------------------------------------------
   integer(c_int), parameter :: CUSOLVER_EIG_MODE_NOVECTOR = 0
   integer(c_int), parameter :: CUSOLVER_EIG_MODE_VECTOR = 1

   ! `uplo` is cublasFillMode_t, declared in cublas.f90 -- CUBLAS_FILL_MODE_*.

   interface

      ! ---- handle lifetime ----------------------------------------------

      integer(c_int) function cusolverDnCreate(handle) result(stat) &
         bind(C, name="cusolverDnCreate")
         import :: c_ptr, c_int
         implicit none
         type(c_ptr), intent(out) :: handle  !! cusolverDnHandle_t*
      end function cusolverDnCreate

      integer(c_int) function cusolverDnDestroy(handle) result(stat) &
         bind(C, name="cusolverDnDestroy")
         import :: c_ptr, c_int
         implicit none
         type(c_ptr), value :: handle
      end function cusolverDnDestroy

      integer(c_int) function cusolverDnSetStream(handle, streamId) result(stat) &
         bind(C, name="cusolverDnSetStream")
         !! Put cuSOLVER on the same stream as everything else, so its work
         !! orders against the cuBLAS around it without a device-wide sync.
         import :: c_ptr, c_int
         implicit none
         type(c_ptr), value :: handle
         type(c_ptr), value :: streamId  !! cudaStream_t; null means the default stream
      end function cusolverDnSetStream

      ! ---- symmetric eigensolver, divide and conquer ---------------------

      integer(c_int) function cusolverDnDsyevd_bufferSize(handle, jobz, uplo, n, &
                                                          A, lda, W, lwork) result(stat) &
         bind(C, name="cusolverDnDsyevd_bufferSize")
         !! How many doubles of device workspace the solve below will want
         !!
         !! Depends on `n` but not on the contents of A or W, so one query
         !! before the SCF loop covers every iteration of it.
         import :: c_ptr, c_int
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: jobz     !! CUSOLVER_EIG_MODE_*
         integer(c_int), value :: uplo     !! CUBLAS_FILL_MODE_*
         integer(c_int), value :: n
         type(c_ptr), value :: A           !! device, (lda, n)
         integer(c_int), value :: lda
         type(c_ptr), value :: W           !! device, (n); read for sizing only
         integer(c_int), intent(out) :: lwork  !! HOST, in doubles
      end function cusolverDnDsyevd_bufferSize

      integer(c_int) function cusolverDnDsyevd(handle, jobz, uplo, n, A, lda, W, &
                                               work, lwork, info) result(stat) &
         bind(C, name="cusolverDnDsyevd")
         !! Eigenvalues and eigenvectors of a real symmetric matrix
         !!
         !! A is OVERWRITTEN with the eigenvectors, one per column, so a caller
         !! that still needs the matrix must hand over a copy.
         import :: c_ptr, c_int
         implicit none
         type(c_ptr), value :: handle
         integer(c_int), value :: jobz     !! CUSOLVER_EIG_MODE_*
         integer(c_int), value :: uplo     !! CUBLAS_FILL_MODE_*
         integer(c_int), value :: n
         type(c_ptr), value :: A           !! device, (lda, n), overwritten
         integer(c_int), value :: lda
         type(c_ptr), value :: W           !! device, (n), eigenvalues ascending
         type(c_ptr), value :: work        !! device, lwork doubles
         integer(c_int), value :: lwork
         type(c_ptr), value :: info
            !! DEVICE int. 0 on success; <0 means argument -info was bad, >0
            !! means that many off-diagonal elements failed to converge.
      end function cusolverDnDsyevd

   end interface

contains

   pure function cusolver_status_name(status) result(name)
      !! Decode a cuSOLVER status for an error message
      integer(c_int), intent(in) :: status
      character(len=:), allocatable :: name

      select case (status)
      case (CUSOLVER_STATUS_SUCCESS)
         name = "CUSOLVER_STATUS_SUCCESS"
      case (CUSOLVER_STATUS_NOT_INITIALIZED)
         name = "CUSOLVER_STATUS_NOT_INITIALIZED"
      case (CUSOLVER_STATUS_ALLOC_FAILED)
         name = "CUSOLVER_STATUS_ALLOC_FAILED"
      case (CUSOLVER_STATUS_INVALID_VALUE)
         name = "CUSOLVER_STATUS_INVALID_VALUE"
      case (CUSOLVER_STATUS_ARCH_MISMATCH)
         name = "CUSOLVER_STATUS_ARCH_MISMATCH"
      case (CUSOLVER_STATUS_MAPPING_ERROR)
         name = "CUSOLVER_STATUS_MAPPING_ERROR"
      case (CUSOLVER_STATUS_EXECUTION_FAILED)
         name = "CUSOLVER_STATUS_EXECUTION_FAILED"
      case (CUSOLVER_STATUS_INTERNAL_ERROR)
         name = "CUSOLVER_STATUS_INTERNAL_ERROR"
      case (CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED)
         name = "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED"
      case (CUSOLVER_STATUS_NOT_SUPPORTED)
         name = "CUSOLVER_STATUS_NOT_SUPPORTED"
      case (CUSOLVER_STATUS_ZERO_PIVOT)
         name = "CUSOLVER_STATUS_ZERO_PIVOT"
      case (CUSOLVER_STATUS_INVALID_LICENSE)
         name = "CUSOLVER_STATUS_INVALID_LICENSE"
      case (CUSOLVER_STATUS_INVALID_WORKSPACE)
         name = "CUSOLVER_STATUS_INVALID_WORKSPACE"
      case default
         name = "CUSOLVER_STATUS_UNKNOWN"
      end select
   end function cusolver_status_name

end module cusolver
