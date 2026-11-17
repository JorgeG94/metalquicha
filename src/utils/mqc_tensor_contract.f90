!! Pairwise tensor contraction named by an index string
module mqc_tensor_contract
   !! `call contract("ijac,cb->ijab", a, b, c, scale)` adds `scale` times the
   !! contraction of `a` with `b` into `c`.
   !!
   !! The string is the one an einsum would take, restricted to what a single
   !! matrix product can do: every label of an operand appears either in the
   !! other operand, and is summed, or in the result, and is kept. A label in
   !! all three, a label traced within one operand and a repeated label are
   !! refused. An outer product, with nothing summed, is allowed.
   !!
   !! Each call permutes both operands into matrices, multiplies them and
   !! permutes the product into the result's order, so it copies its operands
   !! and holds one array the size of the result. It is for writing a set of
   !! equations down as they are published; a contraction that dominates a run
   !! wants its own loop.
   !!
   !! Operands and result are rank two or rank four, in any combination.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use omp_lib, only: omp_get_max_threads
   implicit none
   private

   public :: contract
   public :: transposed

   integer, parameter :: MAX_RANK = 4
      !! Highest rank of an operand or a result
   real(dp), parameter :: THREAD_ABOVE = 2.0e6_dp
      !! A product is split over threads when it takes more multiplications
      !! than this, and a permutation when it moves more elements

   interface contract
      module procedure contract_44_4, contract_44_2, contract_42_4, contract_24_4
      module procedure contract_42_2, contract_24_2, contract_22_2, contract_22_4
   end interface contract

contains

   subroutine contract_44_4(spec, a, b, c, scale)
      !! Rank four with rank four into rank four
      character(len=*), intent(in) :: spec
      real(dp), intent(in) :: a(:, :, :, :), b(:, :, :, :)
      real(dp), intent(inout) :: c(:, :, :, :)
      real(dp), intent(in) :: scale

      real(dp), allocatable :: r(:)

      call contract_flat(spec, reshape(a, [size(a)]), shape(a), &
                         reshape(b, [size(b)]), shape(b), shape(c), r)
      c = c + scale*reshape(r, shape(c))
   end subroutine contract_44_4

   subroutine contract_44_2(spec, a, b, c, scale)
      !! Rank four with rank four into rank two
      character(len=*), intent(in) :: spec
      real(dp), intent(in) :: a(:, :, :, :), b(:, :, :, :)
      real(dp), intent(inout) :: c(:, :)
      real(dp), intent(in) :: scale

      real(dp), allocatable :: r(:)

      call contract_flat(spec, reshape(a, [size(a)]), shape(a), &
                         reshape(b, [size(b)]), shape(b), shape(c), r)
      c = c + scale*reshape(r, shape(c))
   end subroutine contract_44_2

   subroutine contract_42_4(spec, a, b, c, scale)
      !! Rank four with rank two into rank four
      character(len=*), intent(in) :: spec
      real(dp), intent(in) :: a(:, :, :, :), b(:, :)
      real(dp), intent(inout) :: c(:, :, :, :)
      real(dp), intent(in) :: scale

      real(dp), allocatable :: r(:)

      call contract_flat(spec, reshape(a, [size(a)]), shape(a), &
                         reshape(b, [size(b)]), shape(b), shape(c), r)
      c = c + scale*reshape(r, shape(c))
   end subroutine contract_42_4

   subroutine contract_24_4(spec, a, b, c, scale)
      !! Rank two with rank four into rank four
      character(len=*), intent(in) :: spec
      real(dp), intent(in) :: a(:, :), b(:, :, :, :)
      real(dp), intent(inout) :: c(:, :, :, :)
      real(dp), intent(in) :: scale

      real(dp), allocatable :: r(:)

      call contract_flat(spec, reshape(a, [size(a)]), shape(a), &
                         reshape(b, [size(b)]), shape(b), shape(c), r)
      c = c + scale*reshape(r, shape(c))
   end subroutine contract_24_4

   subroutine contract_42_2(spec, a, b, c, scale)
      !! Rank four with rank two into rank two
      character(len=*), intent(in) :: spec
      real(dp), intent(in) :: a(:, :, :, :), b(:, :)
      real(dp), intent(inout) :: c(:, :)
      real(dp), intent(in) :: scale

      real(dp), allocatable :: r(:)

      call contract_flat(spec, reshape(a, [size(a)]), shape(a), &
                         reshape(b, [size(b)]), shape(b), shape(c), r)
      c = c + scale*reshape(r, shape(c))
   end subroutine contract_42_2

   subroutine contract_24_2(spec, a, b, c, scale)
      !! Rank two with rank four into rank two
      character(len=*), intent(in) :: spec
      real(dp), intent(in) :: a(:, :), b(:, :, :, :)
      real(dp), intent(inout) :: c(:, :)
      real(dp), intent(in) :: scale

      real(dp), allocatable :: r(:)

      call contract_flat(spec, reshape(a, [size(a)]), shape(a), &
                         reshape(b, [size(b)]), shape(b), shape(c), r)
      c = c + scale*reshape(r, shape(c))
   end subroutine contract_24_2

   subroutine contract_22_2(spec, a, b, c, scale)
      !! Rank two with rank two into rank two
      character(len=*), intent(in) :: spec
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), intent(inout) :: c(:, :)
      real(dp), intent(in) :: scale

      real(dp), allocatable :: r(:)

      call contract_flat(spec, reshape(a, [size(a)]), shape(a), &
                         reshape(b, [size(b)]), shape(b), shape(c), r)
      c = c + scale*reshape(r, shape(c))
   end subroutine contract_22_2

   subroutine contract_22_4(spec, a, b, c, scale)
      !! Rank two with rank two into rank four, an outer product
      character(len=*), intent(in) :: spec
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), intent(inout) :: c(:, :, :, :)
      real(dp), intent(in) :: scale

      real(dp), allocatable :: r(:)

      call contract_flat(spec, reshape(a, [size(a)]), shape(a), &
                         reshape(b, [size(b)]), shape(b), shape(c), r)
      c = c + scale*reshape(r, shape(c))
   end subroutine contract_22_4

   function transposed(x, order) result(y)
      !! `x` with its indices reordered: index `k` of the result is index
      !! `order(k)` of `x`
      real(dp), intent(in) :: x(:, :, :, :)
      integer, intent(in) :: order(4)   !! A permutation of 1 to 4
      real(dp), allocatable :: y(:, :, :, :)

      real(dp), allocatable :: flat(:)
      integer :: extents(4)

      extents = shape(x)
      call permute(reshape(x, [size(x)]), extents, order, flat)
      y = reshape(flat, [extents(order(1)), extents(order(2)), extents(order(3)), &
                         extents(order(4))])
   end function transposed

   subroutine contract_flat(spec, a, shape_a, b, shape_b, shape_c, r)
      !! The contraction over flattened operands, in the result's index order
      character(len=*), intent(in) :: spec
      real(dp), intent(in) :: a(:), b(:)
      integer, intent(in) :: shape_a(:), shape_b(:), shape_c(:)
      real(dp), allocatable, intent(out) :: r(:)   !! (product of `shape_c`)

      character(len=MAX_RANK) :: la, lb, lc
      integer :: na, nb, nc
      integer :: perm_a(MAX_RANK), perm_b(MAX_RANK), perm_r(MAX_RANK)
      integer :: shape_r(2*MAX_RANK)
      character(len=2*MAX_RANK) :: lr
      integer :: n_free_a, n_free_b, n_sum
      integer :: rows, cols, inner
      integer :: k, at, where_b
      real(dp), allocatable :: am(:), bm(:), rm(:)

      call split_spec(spec, la, na, lb, nb, lc, nc)
      if (na /= size(shape_a) .or. nb /= size(shape_b) .or. nc /= size(shape_c)) then
         error stop "contract: '"//spec//"' does not match the ranks it was given"
      end if
      call refuse_repeats(spec, la(1:na))
      call refuse_repeats(spec, lb(1:nb))
      call refuse_repeats(spec, lc(1:nc))

      ! The first operand goes to (kept, summed) and the second to (summed,
      ! kept), with the summed labels in the first operand's order in both.
      n_free_a = 0
      n_sum = 0
      rows = 1
      inner = 1
      lr = ""
      do k = 1, na
         if (index(lc(1:nc), la(k:k)) > 0) then
            if (index(lb(1:nb), la(k:k)) > 0) then
               error stop "contract: '"//spec//"' has a label in both operands and the result"
            end if
            n_free_a = n_free_a + 1
            perm_a(n_free_a) = k
            rows = rows*shape_a(k)
            lr(n_free_a:n_free_a) = la(k:k)
            shape_r(n_free_a) = shape_a(k)
         end if
      end do
      at = n_free_a
      do k = 1, na
         if (index(lc(1:nc), la(k:k)) > 0) cycle
         where_b = index(lb(1:nb), la(k:k))
         if (where_b == 0) then
            error stop "contract: '"//spec//"' has a label that is neither summed nor kept"
         end if
         if (shape_a(k) /= shape_b(where_b)) then
            error stop "contract: '"//spec//"' sums over two extents that differ"
         end if
         n_sum = n_sum + 1
         at = at + 1
         perm_a(at) = k
         perm_b(n_sum) = where_b
         inner = inner*shape_a(k)
      end do
      n_free_b = 0
      cols = 1
      do k = 1, nb
         if (index(la(1:na), lb(k:k)) > 0) cycle
         if (index(lc(1:nc), lb(k:k)) == 0) then
            error stop "contract: '"//spec//"' has a label that is neither summed nor kept"
         end if
         n_free_b = n_free_b + 1
         perm_b(n_sum + n_free_b) = k
         cols = cols*shape_b(k)
         lr(n_free_a + n_free_b:n_free_a + n_free_b) = lb(k:k)
         shape_r(n_free_a + n_free_b) = shape_b(k)
      end do
      if (n_free_a + n_free_b /= nc) then
         error stop "contract: '"//spec//"' asks for a label neither operand has"
      end if

      ! Where each index of the result sits in the product, and that the two
      ! agree on its extent.
      do k = 1, nc
         perm_r(k) = index(lr(1:nc), lc(k:k))
         if (shape_r(perm_r(k)) /= shape_c(k)) then
            error stop "contract: '"//spec//"' gives a result of another shape"
         end if
      end do

      call permute(a, shape_a, perm_a(1:na), am)
      call permute(b, shape_b, perm_b(1:nb), bm)
      allocate (rm(rows*cols))
      call multiply(am, bm, rows, inner, cols, rm)
      call permute(rm, shape_r(1:nc), perm_r(1:nc), r)
   end subroutine contract_flat

   subroutine multiply(am, bm, rows, inner, cols, rm)
      !! `rm(rows, cols) = am(rows, inner) bm(inner, cols)` over flat storage
      !!
      !! The BLAS this is built against is sequential, so a product worth
      !! threading is split by hand into one block of columns for each thread.
      integer, intent(in) :: rows, inner, cols
      real(dp), intent(in) :: am(rows, inner), bm(inner, cols)
      real(dp), intent(out) :: rm(rows, cols)

      integer :: n_blocks, block, first, last

      if (rows == 0 .or. cols == 0 .or. inner == 0) then
         rm = 0.0_dp
         return
      end if
      n_blocks = min(omp_get_max_threads(), cols)
      if (n_blocks < 2 .or. real(rows, dp)*real(inner, dp)*real(cols, dp) < THREAD_ABOVE) then
         call pic_gemm(am, bm, rm)
         return
      end if
      !$omp parallel do default(none) shared(am, bm, rm, n_blocks, cols) &
      !$omp    private(block, first, last) schedule(static)
      do block = 1, n_blocks
         first = (block - 1)*cols/n_blocks + 1
         last = block*cols/n_blocks
         if (last >= first) call pic_gemm(am, bm(:, first:last), rm(:, first:last))
      end do
      !$omp end parallel do
   end subroutine multiply

   subroutine split_spec(spec, la, na, lb, nb, lc, nc)
      !! The three label strings of `"a,b->c"`
      character(len=*), intent(in) :: spec
      character(len=MAX_RANK), intent(out) :: la, lb, lc
      integer, intent(out) :: na, nb, nc

      integer :: comma, arrow

      comma = index(spec, ",")
      arrow = index(spec, "->")
      if (comma < 2 .or. arrow < comma + 2 .or. arrow + 2 > len_trim(spec)) then
         error stop "contract: '"//spec//"' is not of the form 'ab,bc->ac'"
      end if
      na = comma - 1
      nb = arrow - comma - 1
      nc = len_trim(spec) - arrow - 1
      if (na > MAX_RANK .or. nb > MAX_RANK .or. nc > MAX_RANK) then
         error stop "contract: '"//spec//"' has an operand of too high a rank"
      end if
      la = spec(1:comma - 1)
      lb = spec(comma + 1:arrow - 1)
      lc = spec(arrow + 2:len_trim(spec))
   end subroutine split_spec

   subroutine refuse_repeats(spec, labels)
      !! Stop on a label that appears twice within one operand
      character(len=*), intent(in) :: spec, labels

      integer :: k

      do k = 2, len(labels)
         if (index(labels(1:k - 1), labels(k:k)) > 0) then
            error stop "contract: '"//spec//"' repeats a label within one operand"
         end if
      end do
   end subroutine refuse_repeats

   subroutine permute(src, extents, perm, dst)
      !! `dst`, whose k-th index is index `perm(k)` of `src`, both flattened
      real(dp), intent(in) :: src(:)
      integer, intent(in) :: extents(:)   !! Shape of `src`
      integer, intent(in) :: perm(:)
      real(dp), allocatable, intent(out) :: dst(:)

      integer :: n, k, at, first, step, offset, slab, outer, n_outer
      integer :: stride(MAX_RANK), jump(MAX_RANK), top(MAX_RANK), counter(MAX_RANK)

      n = size(extents)
      allocate (dst(size(src)))
      if (size(src) == 0) return

      ! `jump(k)` is how far `src` moves when the result's k-th index advances.
      stride(1) = 1
      do k = 2, n
         stride(k) = stride(k - 1)*extents(k - 1)
      end do
      do k = 1, n
         jump(k) = stride(perm(k))
         top(k) = extents(perm(k))
      end do
      first = top(1)
      step = jump(1)

      ! One value of the result's last index is one contiguous slab of it, and
      ! the slabs are independent.
      n_outer = top(n)
      slab = size(src)/n_outer
      !$omp parallel do default(none) if (real(size(src), dp) > THREAD_ABOVE) &
      !$omp    shared(src, dst, n, jump, top, first, step, slab, n_outer) &
      !$omp    private(outer, counter, at, offset, k) schedule(static)
      do outer = 0, n_outer - 1
         counter(1:n) = 0
         at = outer*slab
         if (n == 1) then
            dst(outer + 1) = src(1 + outer*jump(1))
            cycle
         end if
         do
            offset = 1 + outer*jump(n)
            do k = 2, n - 1
               offset = offset + counter(k)*jump(k)
            end do
            do k = 1, first
               dst(at + k) = src(offset + (k - 1)*step)
            end do
            at = at + first
            ! Advance every index but the first, which the loop above ran, and
            ! the last, which the slab fixes.
            k = 2
            do
               if (k > n - 1) exit
               counter(k) = counter(k) + 1
               if (counter(k) < top(k)) exit
               counter(k) = 0
               k = k + 1
            end do
            if (k > n - 1) exit
         end do
      end do
      !$omp end parallel do
   end subroutine permute

end module mqc_tensor_contract
