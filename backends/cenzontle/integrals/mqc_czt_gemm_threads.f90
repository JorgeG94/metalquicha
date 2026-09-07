module mqc_czt_gemm_threads
   !! Matrix products threaded from outside the BLAS call
   !!
   !! The BLAS this program links is sequential by design: `tools/run.sh` pins it
   !! to one thread so that fragment ranks do not compete for cores. A product
   !! large enough to matter therefore has to be split here, or it runs on one
   !! core while the rest of the node waits -- the `(A-B)(A+B)` product of the
   !! response solver cost twelve seconds of a twenty-second stage that way.
   use pic_types, only: dp, int64
   use pic_blas_interfaces, only: pic_gemm
   use omp_lib, only: omp_get_max_threads
   implicit none
   private

   public :: gemm_over_columns
   public :: gemm_over_inner

   integer(int64), parameter :: PARTIAL_LIMIT = 2_int64*1024_int64**3
      !! Bytes all threads' partial products of `gemm_over_inner` may hold together

contains

   subroutine gemm_over_columns(a, b, c, transa, transb)
      !! C = op(A) op(B), with C's columns split across threads
      !!
      !! Thread t computes `C(:, s0:s1)` from all of `op(A)` and the slice of
      !! `op(B)` that feeds those columns. With `B` untransposed that slice is
      !! `B(:, s0:s1)`, contiguous, and nothing is copied; with `transb` "T"
      !! it is the row range `B(s0:s1, :)`, which the BLAS call copies once. The
      !! ranges are disjoint, so there is no reduction and no ordering question.
      !! `transa` and `transb` are "N" (the default) or "T", passed straight to
      !! the BLAS. `C` is overwritten.
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), intent(inout) :: c(:, :)
      character(len=1), intent(in), optional :: transa
      character(len=1), intent(in), optional :: transb

      integer :: n, nchunk, width, chunk, s0, s1
      logical :: b_rows

      n = size(c, 2)
      nchunk = omp_get_max_threads()
      if (nchunk < 1) nchunk = 1
      if (nchunk > n) nchunk = n
      width = (n + nchunk - 1)/nchunk
      b_rows = .false.
      if (present(transb)) b_rows = transb == "T" .or. transb == "t"

      !$omp parallel do default(none) &
      !$omp    shared(a, b, c, n, nchunk, width, transa, transb, b_rows) &
      !$omp    private(chunk, s0, s1) schedule(static)
      do chunk = 1, nchunk
         s0 = (chunk - 1)*width + 1
         s1 = min(chunk*width, n)
         if (s0 > s1) cycle
         if (b_rows) then
            call pic_gemm(a, b(s0:s1, :), c(:, s0:s1), transa=transa, transb=transb)
         else
            call pic_gemm(a, b(:, s0:s1), c(:, s0:s1), transa=transa)
         end if
      end do
      !$omp end parallel do
   end subroutine gemm_over_columns

   subroutine gemm_over_inner(a, b, c)
      !! C = A^T B, with the summed index split across threads
      !!
      !! For the products whose result is small and whose shared index is long:
      !! `a` is `(k, m)`, `b` is `(k, n)`, and `c(m, n)` sums over `k`. Splitting
      !! `c`'s columns for a shape like this makes every thread read all of `a`
      !! -- 1.1 GB, a hundred and twelve times over, for the virtual block of the
      !! MP2 correction density -- so the row range of both operands is split
      !! instead, each thread forms its own `m x n` partial product from its
      !! range, and the partials are summed once at the end. A row range of a
      !! column-major array is strided, so the BLAS call packs each chunk; that
      !! reads every element of either operand once, which is the point. The
      !! thread count is capped so the partials together stay under
      !! `PARTIAL_LIMIT`. `C` is overwritten.
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), intent(inout) :: c(:, :)

      real(dp), allocatable :: part(:, :)
      integer :: k, m, n, nchunk, width, chunk, s0, s1
      integer(int64) :: per_partial

      k = size(a, 1)
      m = size(a, 2)
      n = size(b, 2)
      c = 0.0_dp
      if (k < 1 .or. m < 1 .or. n < 1) return
      per_partial = int(m, int64)*int(n, int64)*8_int64
      nchunk = int(min(int(omp_get_max_threads(), int64), &
                       max(1_int64, PARTIAL_LIMIT/max(per_partial, 1_int64))))
      if (nchunk > k) nchunk = k
      width = (k + nchunk - 1)/nchunk

      !$omp parallel default(none) num_threads(nchunk) &
      !$omp    shared(a, b, c, k, m, n, nchunk, width) private(chunk, s0, s1, part)
      allocate (part(m, n))
      part = 0.0_dp
      !$omp do schedule(static)
      do chunk = 1, nchunk
         s0 = (chunk - 1)*width + 1
         s1 = min(chunk*width, k)
         if (s0 > s1) cycle
         call pic_gemm(a(s0:s1, :), b(s0:s1, :), part, transa="T", beta=1.0_dp)
      end do
      !$omp end do
      !$omp critical
      c = c + part
      !$omp end critical
      deallocate (part)
      !$omp end parallel
   end subroutine gemm_over_inner

end module mqc_czt_gemm_threads
