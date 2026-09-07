module mqc_czt_gemm_threads
   !! Dense linear algebra threaded from outside the BLAS call
   !!
   !! The BLAS this program links is sequential by design: `tools/run.sh` pins it
   !! to one thread so that fragment ranks do not compete for cores. A product
   !! or a factorization large enough to matter therefore has to be split here,
   !! or it runs on one core while the rest of the node waits -- the `(A-B)(A+B)`
   !! product of the response solver cost twelve seconds of a twenty-second
   !! stage that way, and its thirteen LU factorizations ran on thirteen of
   !! 128 cores for 72 s at 17710 pairs.
   use pic_types, only: dp, int64, default_int
   use pic_blas_interfaces, only: pic_gemm, pic_dgemm_x, pic_trsm
   use pic_lapack_interfaces, only: pic_getrf
   use omp_lib, only: omp_get_max_threads
   implicit none
   private

   public :: gemm_over_columns
   public :: gemm_over_inner
   public :: getrf_threaded

   integer, parameter :: LU_PANEL = 256
      !! Columns factorized per panel of `getrf_threaded`. The panel is the
      !! serial part, `n^2 nb` flops in all, so smaller is faster there; the
      !! trailing update's GEMMs have `nb` as their inner dimension, so smaller
      !! is slower there. 256 keeps the panel under two seconds at 17710 while
      !! the update runs the BLAS at its full rate.

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

   subroutine getrf_threaded(n, a, ipiv, info, threads)
      !! LU with partial pivoting, the trailing updates split across threads
      !!
      !! What `getrf` returns -- `L` unit lower below the diagonal, `U` on and
      !! above it, `ipiv` the row interchanges, all in LAPACK's own layout, so
      !! `getrs` reads the result unchanged -- built as a right-looking blocked
      !! factorization over the sequential BLAS. Each panel of `LU_PANEL`
      !! columns is factorized by `getrf` on one thread; its interchanges are
      !! applied to the rest of the rows, to the left as LAPACK does; the `U`
      !! block above the trailing matrix is triangular-solved; and the trailing
      !! matrix is updated by column blocks on `threads` threads. The update is
      !! `2/3 n^3` flops and the panel `n^2 nb`, so the serial part is the
      !! `nb/n` fraction. `info` is `getrf`'s: zero, or the first zero pivot.
      !! `a` is explicit-shape so its elements can be handed to the BLAS with a
      !! leading dimension and nothing is copied for the update; the panel and
      !! the triangular solve go through sections, which the compiler packs,
      !! and those are `n nb` each.
      integer, intent(in) :: n
      real(dp), intent(inout) :: a(n, n)
      integer, intent(out) :: ipiv(n)
      integer, intent(out) :: info
      integer, intent(in), optional :: threads
         !! The team size for the threaded parts; absent is the OpenMP maximum.
         !! Passed rather than read because the caller may itself sit in a
         !! parallel region, with this as the inner level.

      integer :: nthr, k, kb, i, j, p, m22, n22, nchunk, width, chunk, j0, j1, pinfo
      integer, allocatable :: piv(:)
      real(dp) :: swap
      integer(default_int) :: bm, bn, bk, lda

      info = 0
      nthr = omp_get_max_threads()
      if (present(threads)) nthr = max(1, threads)
      lda = int(n, default_int)
      allocate (piv(LU_PANEL))

      do k = 1, n, LU_PANEL
         kb = min(LU_PANEL, n - k + 1)

         ! The panel, on one thread. Its pivots are local to the panel's rows.
         call pic_getrf(a(k:n, k:k + kb - 1), piv(1:kb), pinfo)
         if (pinfo > 0 .and. info == 0) info = pinfo + k - 1
         do i = 1, kb
            ipiv(k + i - 1) = piv(i) + k - 1
         end do

         ! The interchanges, applied to every column outside the panel -- the
         ! `L` columns already done on the left as well as the trailing ones,
         ! which is what makes the result LAPACK's.
         !$omp parallel do default(none) num_threads(nthr) schedule(static) &
         !$omp    shared(a, ipiv, n, k, kb) private(i, j, p, swap)
         do j = 1, n
            if (j >= k .and. j < k + kb) cycle
            do i = k, k + kb - 1
               p = ipiv(i)
               if (p /= i) then
                  swap = a(i, j)
                  a(i, j) = a(p, j)
                  a(p, j) = swap
               end if
            end do
         end do
         !$omp end parallel do

         if (k + kb > n) exit
         m22 = n - k - kb + 1
         n22 = m22
         nchunk = min(nthr, n22)
         width = (n22 + nchunk - 1)/nchunk
         bm = int(m22, default_int)
         bk = int(kb, default_int)

         ! `U12 = L11^-1 A12` and `A22 -= L21 U12`, both by column blocks of the
         ! trailing matrix; each block's solve feeds only its own update, so the
         ! two share one region and one barrier.
         !$omp parallel do default(none) num_threads(nthr) schedule(static) &
         !$omp    shared(a, n, k, kb, m22, n22, nchunk, width, lda, bm, bk) &
         !$omp    private(chunk, j0, j1, bn)
         do chunk = 1, nchunk
            j0 = k + kb + (chunk - 1)*width
            j1 = min(k + kb + chunk*width - 1, n)
            if (j0 > j1) cycle
            call pic_trsm(a(k:k + kb - 1, k:k + kb - 1), a(k:k + kb - 1, j0:j1), &
                          side="L", uplo="L", transa="N", diag="U")
            bn = int(j1 - j0 + 1, default_int)
            call pic_dgemm_x("N", "N", bm, bn, bk, -1.0_dp, a(k + kb, k), lda, &
                            a(k, j0), lda, 1.0_dp, a(k + kb, j0), lda)
         end do
         !$omp end parallel do
      end do
      deallocate (piv)
   end subroutine getrf_threaded

end module mqc_czt_gemm_threads
