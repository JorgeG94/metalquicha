module test_mqc_czt_gemm_threads
   !! The threaded linear algebra against the BLAS it is threaded over.
   !!
   !! Every routine in `mqc_czt_gemm_threads` recomputes something a single BLAS
   !! or LAPACK call already computes, on more cores. So the test is the
   !! identity: the same product to rounding, the same solution from the LU.
   !! Sizes straddle the panel width so that a partial last panel, a single
   !! panel and several full panels are all exercised.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_getrf, pic_getrs
   use mqc_czt_gemm_threads, only: gemm_over_columns, gemm_over_inner, getrf_threaded
   implicit none
   private
   public :: collect_mqc_czt_gemm_threads_tests

   real(dp), parameter :: TOL = 1.0e-11_dp

contains

   subroutine collect_mqc_czt_gemm_threads_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("column_split_gemm_equals_the_blas", test_columns), &
                  new_unittest("inner_split_gemm_equals_the_blas", test_inner), &
                  new_unittest("threaded_lu_solves_like_getrf", test_lu), &
                  new_unittest("threaded_lu_reports_a_singular_matrix", test_lu_singular) &
                  ]
   end subroutine collect_mqc_czt_gemm_threads_tests

   subroutine fill(a, seed)
      !! Deterministic pseudo-random entries in [-1, 1)
      real(dp), intent(out) :: a(:, :)
      integer, intent(in) :: seed
      integer :: i, j
      real(dp) :: x

      x = real(seed, dp)*0.618033988749895_dp
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            x = mod(x*9821.0_dp + 0.211327_dp, 1.0_dp)
            a(i, j) = 2.0_dp*x - 1.0_dp
         end do
      end do
   end subroutine fill

   subroutine test_columns(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: a(:, :), b(:, :), c(:, :), ref(:, :)

      allocate (a(37, 53), b(53, 301), c(37, 301), ref(37, 301))
      call fill(a, 1)
      call fill(b, 2)
      call pic_gemm(a, b, ref)
      call gemm_over_columns(a, b, c)
      call check(error, maxval(abs(c - ref)) < TOL, "C = A B by columns must equal the BLAS")
      if (allocated(error)) return

      deallocate (a, b, c, ref)
      allocate (a(53, 37), b(53, 301), c(37, 301), ref(37, 301))
      call fill(a, 3)
      call fill(b, 4)
      call pic_gemm(a, b, ref, transa="T")
      call gemm_over_columns(a, b, c, transa="T")
      call check(error, maxval(abs(c - ref)) < TOL, "C = A^T B by columns must equal the BLAS")
      if (allocated(error)) return

      deallocate (a, b, c, ref)
      allocate (a(37, 53), b(301, 53), c(37, 301), ref(37, 301))
      call fill(a, 5)
      call fill(b, 6)
      call pic_gemm(a, b, ref, transb="T")
      call gemm_over_columns(a, b, c, transb="T")
      call check(error, maxval(abs(c - ref)) < TOL, "C = A B^T by columns must equal the BLAS")
   end subroutine test_columns

   subroutine test_inner(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp), allocatable :: a(:, :), b(:, :), c(:, :), ref(:, :)

      allocate (a(1201, 17), b(1201, 23), c(17, 23), ref(17, 23))
      call fill(a, 7)
      call fill(b, 8)
      call pic_gemm(a, b, ref, transa="T")
      call gemm_over_inner(a, b, c)
      call check(error, maxval(abs(c - ref)) < 1.0e-9_dp, &
                 "C = A^T B over the inner index must equal the BLAS")
   end subroutine test_inner

   subroutine test_lu(error)
      !! Three sizes: below one panel, exactly two, and two and a fraction
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: SIZES(3) = [100, 512, 700]
      real(dp), allocatable :: a(:, :), lu_ref(:, :), lu_thr(:, :), x_ref(:, :), x_thr(:, :)
      integer, allocatable :: ipiv_ref(:), ipiv_thr(:)
      integer :: n, i, s, info_ref, info_thr
      character(len=64) :: what

      do s = 1, size(SIZES)
         n = SIZES(s)
         allocate (a(n, n), lu_ref(n, n), lu_thr(n, n), x_ref(n, 3), x_thr(n, 3))
         allocate (ipiv_ref(n), ipiv_thr(n))
         call fill(a, 10 + s)
         ! Diagonally weighted enough to be well conditioned, not enough to make
         ! the pivoting trivial.
         do i = 1, n
            a(i, i) = a(i, i) + 2.0_dp
         end do
         call fill(x_ref, 20 + s)
         x_thr = x_ref
         lu_ref = a
         lu_thr = a
         call pic_getrf(lu_ref, ipiv_ref, info_ref)
         call getrf_threaded(n, lu_thr, ipiv_thr, info_thr, threads=4)
         write (what, "(A,I0)") "n = ", n
         call check(error, info_ref == 0 .and. info_thr == 0, &
                    "both factorizations must succeed at "//trim(what))
         if (allocated(error)) return
         call pic_getrs(lu_ref, ipiv_ref, x_ref)
         call pic_getrs(lu_thr, ipiv_thr, x_thr)
         call check(error, maxval(abs(x_thr - x_ref)) < TOL*maxval(abs(x_ref)), &
                    "the threaded LU must solve the system like getrf at "//trim(what))
         if (allocated(error)) return
         deallocate (a, lu_ref, lu_thr, x_ref, x_thr, ipiv_ref, ipiv_thr)
      end do
   end subroutine test_lu

   subroutine test_lu_singular(error)
      !! A zero column has to come back as a zero pivot, not as a factor
      !!
      !! Exactly zero, as LAPACK's own `info` requires: a *repeated* column
      !! leaves a pivot of 1e-13 that neither `getrf` reports.
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: N = 300
      real(dp), allocatable :: a(:, :)
      integer, allocatable :: ipiv(:)
      integer :: info

      allocate (a(N, N), ipiv(N))
      call fill(a, 31)
      a(:, 270) = 0.0_dp
      call getrf_threaded(N, a, ipiv, info, threads=4)
      call check(error, info > 0, "a singular matrix must be reported through info")
   end subroutine test_lu_singular

end module test_mqc_czt_gemm_threads

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_gemm_threads, only: collect_mqc_czt_gemm_threads_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_gemm_threads", collect_mqc_czt_gemm_threads_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
