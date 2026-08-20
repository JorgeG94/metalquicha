!! The status codes every C entry point returns
module mqc_capi_status
   !! One declaration of the C API's return values, for all of it.
   !!
   !! These are the contract with anything that is not Fortran -- the ctypes
   !! wrapper in `python/`, and whatever else links `libmqc.so`. That makes them
   !! the worst constants in the program to hold four copies of: nothing inside
   !! this program compares a code from one module against a code from another,
   !! so the copies could drift apart for a long time without a test noticing,
   !! and the symptom would land in somebody else's language.
   !!
   !! `MQC_BAD_HANDLE` is separate from `MQC_FAIL` because a caller can act on
   !! it: a bad handle is a bug in the caller's bookkeeping, where a failure is
   !! usually a bad calculation.
   use, intrinsic :: iso_c_binding, only: c_int
   implicit none
   private

   public :: MQC_OK, MQC_FAIL, MQC_BAD_HANDLE

   integer(c_int), parameter :: MQC_OK = 0          !! The call did what was asked
   integer(c_int), parameter :: MQC_FAIL = 1        !! It did not; see the last message
   integer(c_int), parameter :: MQC_BAD_HANDLE = 2  !! The handle names nothing

end module mqc_capi_status
