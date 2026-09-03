!! Small exact-arithmetic helpers with no domain of their own
module mqc_math_utils
   !! Pure combinatorial arithmetic, shared rather than reimplemented: a
   !! fragment count and a determinant count want the same function.
   use pic_types, only: int32, int64
   implicit none
   private

   public :: binomial

   interface binomial
      !! `C(n, k)`, exactly, in 64-bit integers, for either integer kind.
      !!
      !! Two specifics because `pic`'s `default_int` is not a fixed kind: it is
      !! `int64` in a USE_INT8 build and `c_int` otherwise, so one specific
      !! would compile in one configuration and not the other.
      module procedure binomial_i32
      module procedure binomial_i64
   end interface binomial

contains

   pure function binomial_i32(n, k) result(c)
      !! `C(n, k)` from default-kind arguments
      integer(int32), intent(in) :: n, k
      integer(int64) :: c

      c = binomial_i64(int(n, int64), int(k, int64))
   end function binomial_i32

   pure function binomial_i64(n, k) result(c)
      !! `C(n, k)`, exactly, in 64-bit integers
      !!
      !! One factor at a time, multiplying before dividing, so no intermediate
      !! exceeds the answer. Out-of-range arguments give zero.
      integer(int64), intent(in) :: n, k
      integer(int64) :: c

      integer(int64) :: i, kk

      c = 0_int64
      if (k < 0_int64 .or. k > n .or. n < 0_int64) return
      kk = min(k, n - k)
      c = 1_int64
      do i = 1_int64, kk
         c = c*(n - kk + i)/i
      end do
   end function binomial_i64

end module mqc_math_utils
