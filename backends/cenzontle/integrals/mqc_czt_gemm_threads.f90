module mqc_czt_gemm_threads
   !! Matrix products threaded from outside the BLAS call
   !!
   !! The BLAS this program links is sequential by design: `tools/run.sh` pins it
   !! to one thread so that fragment ranks do not compete for cores. A product
   !! large enough to matter therefore has to be split here, or it runs on one
   !! core while the rest of the node waits -- the `(A-B)(A+B)` product of the
   !! response solver cost twelve seconds of a twenty-second stage that way.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use omp_lib, only: omp_get_max_threads
   implicit none
   private

   public :: gemm_over_columns

contains

   subroutine gemm_over_columns(a, b, c, transa)
      !! C = op(A) B, with C's columns split across threads
      !!
      !! Thread t computes `C(:, s0:s1)` from all of `op(A)` and `B(:, s0:s1)`.
      !! Every operand slice is a column range, so it is contiguous and nothing
      !! is copied, and the ranges are disjoint, so there is no reduction and
      !! no ordering question. `transa` is "N" (the default) or "T", passed
      !! straight to the BLAS. `C` is overwritten.
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), intent(inout) :: c(:, :)
      character(len=1), intent(in), optional :: transa

      integer :: n, nchunk, width, chunk, s0, s1

      n = size(c, 2)
      nchunk = omp_get_max_threads()
      if (nchunk < 1) nchunk = 1
      if (nchunk > n) nchunk = n
      width = (n + nchunk - 1)/nchunk

      !$omp parallel do default(none) &
      !$omp    shared(a, b, c, n, nchunk, width, transa) &
      !$omp    private(chunk, s0, s1) schedule(static)
      do chunk = 1, nchunk
         s0 = (chunk - 1)*width + 1
         s1 = min(chunk*width, n)
         if (s0 > s1) cycle
         call pic_gemm(a, b(:, s0:s1), c(:, s0:s1), transa=transa)
      end do
      !$omp end parallel do
   end subroutine gemm_over_columns

end module mqc_czt_gemm_threads
