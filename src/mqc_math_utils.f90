!! Small exact-arithmetic helpers with no domain of their own
module mqc_math_utils
   !! Pure combinatorial arithmetic, shared rather than reimplemented.
   !!
   !! `binomial` existed three times before this module did: once in
   !! `mqc_combinatorics` for counting fragments and twice in `src/mcscf` for
   !! sizing CI vectors. The two in `src/mcscf` were identical to the character;
   !! the one in `mqc_combinatorics` was subtly weaker, and that is the reason
   !! this lives here rather than in either of them. Neither domain owns it,
   !! and a fragment count and a determinant count want exactly the same
   !! function.
   use pic_types, only: int32, int64
   implicit none
   private

   public :: binomial

   interface binomial
      !! `C(n, k)`, exactly, in 64-bit integers, for either integer kind.
      !!
      !! Two specifics rather than one because `pic`'s `default_int` is not a
      !! fixed kind -- it is `int64` in a USE_INT8 build and `c_int` otherwise --
      !! so a single specific would compile in one configuration and not the
      !! other.
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
      !! One factor at a time, multiplying before dividing, so the running value
      !! is always the binomial coefficient of a smaller problem and never
      !! exceeds the answer -- which the naive factorial form does by an enormous
      !! margin, overflowing at `n = 21` for a result that fits easily.
      !!
      !! `kk = min(k, n - k)` is not only a saving: it is what keeps the
      !! intermediate products small. Out-of-range arguments give zero rather
      !! than a wrong number, which the count of an empty set should be.
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
