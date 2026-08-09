!! C entry points for building and filtering a fragment term list
module mqc_capi_fraglist
   !! The fragment list, reachable from C and so from Python.
   !!
   !! Three verbs: create a list, read it out, put a different one back. A
   !! screen -- by distance, by a previous run's energies, by whatever a caller
   !! invents -- is read-then-put-back, and needs nothing here to know about
   !! it. That is the whole point: the criterion lives with the caller and the
   !! combinatorics stay in Fortran, where they already are.
   !!
   !! **Handles are opaque and owned by the caller.** `mqc_fraglist_new`
   !! allocates, `mqc_fraglist_free` releases, and nothing here keeps a
   !! registry -- a leaked handle leaks its term list, which for a
   !! twenty-million-term run is not a small leak.
   !!
   !! **Reading out is two calls, deliberately.** Ask for the count, allocate,
   !! then fill. Returning an allocated buffer across the boundary would put
   !! the free on the wrong side of it, and a term list is large enough that
   !! getting that wrong matters.
   !!
   !! The array handed back is packed **one term per row**: term `i` occupies
   !! elements `i*max_level .. i*max_level + max_level - 1`, zero-padded. That
   !! is a C-contiguous `(n_terms, max_level)` array, and the same shape the
   !! per-fragment CSV writes its `m1..m<max_level>` columns in, so a caller
   !! screening against a previous run's CSV can index the two alike. Note the
   !! Fortran store is `(n_terms, max_level)`, so this is a transpose of its
   !! memory order, not a reinterpretation of it.
   !!
   !! Every entry point returns 0 for success and non-zero for failure rather
   !! than reporting through `error_t`: a C caller cannot read one, and a
   !! Python wrapper turns a status into an exception anyway. The message is
   !! retrievable through `mqc_fraglist_last_error` while it is still the most
   !! recent thing that happened.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_int64_t, c_double, &
                                          c_char, c_null_char, c_f_pointer, c_loc, c_associated
   use pic_types, only: dp, default_int, int64
   use mqc_fraglist, only: fraglist_t
   use mqc_error, only: error_t
   implicit none
   private

   public :: mqc_fraglist_new, mqc_fraglist_free
   public :: mqc_fraglist_generate, mqc_fraglist_count, mqc_fraglist_max_level
   public :: mqc_fraglist_get, mqc_fraglist_set, mqc_fraglist_close_subsets
   public :: mqc_fraglist_last_error

   integer(c_int), parameter :: MQC_OK = 0
   integer(c_int), parameter :: MQC_ERROR = 1
   integer(c_int), parameter :: MQC_BAD_HANDLE = 2

   integer, parameter :: MESSAGE_LEN = 512
   character(len=MESSAGE_LEN), save :: last_message = ""
      !! Most recent failure, for a caller that wants to say why. Process-wide
      !! and overwritten by the next failure -- read it immediately or not at
      !! all. Python only ever runs on one rank here, so there is no thread or
      !! rank to race with.

contains

   function mqc_fraglist_new() result(handle) bind(C, name="mqc_fraglist_new")
      !! Allocate an empty term list and return its handle
      type(c_ptr) :: handle

      type(fraglist_t), pointer :: list

      allocate (list)
      handle = c_loc(list)
   end function mqc_fraglist_new

   subroutine mqc_fraglist_free(handle) bind(C, name="mqc_fraglist_free")
      !! Release a term list. Safe on a null handle.
      type(c_ptr), value :: handle

      type(fraglist_t), pointer :: list

      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, list)
      call list%destroy()
      deallocate (list)
   end subroutine mqc_fraglist_free

   function mqc_fraglist_generate(handle, n_monomers, max_level) result(status) &
      bind(C, name="mqc_fraglist_generate")
      !! Fill the list with every combination from pairs up to max_level
      type(c_ptr), value :: handle
      integer(c_int), value :: n_monomers
      integer(c_int), value :: max_level
      integer(c_int) :: status

      type(fraglist_t), pointer :: list
      type(error_t) :: error

      status = MQC_BAD_HANDLE
      if (.not. c_associated(handle)) then
         last_message = "null fragment list handle"
         return
      end if
      call c_f_pointer(handle, list)

      call list%create(int(n_monomers, default_int), int(max_level, default_int), error)
      status = report(error)
   end function mqc_fraglist_generate

   function mqc_fraglist_count(handle) result(n) bind(C, name="mqc_fraglist_count")
      !! Terms currently in the list, or -1 for a bad handle
      type(c_ptr), value :: handle
      integer(c_int64_t) :: n

      type(fraglist_t), pointer :: list

      n = -1_c_int64_t
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, list)
      n = int(list%n_terms, c_int64_t)
   end function mqc_fraglist_count

   function mqc_fraglist_max_level(handle) result(level) &
      bind(C, name="mqc_fraglist_max_level")
      !! Row count of the term array, or -1 for a bad handle
      type(c_ptr), value :: handle
      integer(c_int) :: level

      type(fraglist_t), pointer :: list

      level = -1_c_int
      if (.not. c_associated(handle)) return
      call c_f_pointer(handle, list)
      level = int(list%max_level, c_int)
   end function mqc_fraglist_max_level

   function mqc_fraglist_get(handle, terms, n_terms, max_level) result(status) &
      bind(C, name="mqc_fraglist_get")
      !! Copy the terms into a caller's buffer
      !!
      !! `n_terms` and `max_level` are what the caller allocated for, and are
      !! checked rather than trusted: a buffer sized from a stale count is the
      !! obvious way to use this wrongly, and the result would be a silent
      !! overrun rather than a short read.
      type(c_ptr), value :: handle
      ! Sizes are declared before the array they shape: a dimension expression
      ! can only name symbols already typed. Shaping from the caller's own
      ! arguments rather than leaving the array assumed-size makes a mismatch a
      ! bounds error under a checking build instead of a silent overrun -- and
      ! assumed-size arrays are on this project's forbidden list for that
      ! reason.
      integer(c_int64_t), value :: n_terms
      integer(c_int), value :: max_level
      integer(c_int), intent(inout) :: terms(n_terms*int(max_level, c_int64_t))
         !! One term per row, row-major
      integer(c_int) :: status

      type(fraglist_t), pointer :: list
      integer(int64) :: iterm, base
      integer :: irow

      status = MQC_BAD_HANDLE
      if (.not. c_associated(handle)) then
         last_message = "null fragment list handle"
         return
      end if
      call c_f_pointer(handle, list)

      if (n_terms < list%n_terms .or. int(max_level, default_int) < list%max_level) then
         last_message = "mqc_fraglist_get: the supplied buffer is smaller than the list"
         status = MQC_ERROR
         return
      end if

      do iterm = 1, list%n_terms
         base = (iterm - 1)*int(max_level, int64)
         do irow = 1, list%max_level
            terms(base + irow) = int(list%terms(iterm, irow), c_int)
         end do
         do irow = list%max_level + 1, max_level
            terms(base + irow) = 0_c_int
         end do
      end do

      status = MQC_OK
   end function mqc_fraglist_get

   function mqc_fraglist_set(handle, terms, n_terms, max_level) result(status) &
      bind(C, name="mqc_fraglist_set")
      !! Replace the list with the caller's terms
      !!
      !! This is where a screened list comes back in, and where a restart hands
      !! over the terms a dead run never reached.
      type(c_ptr), value :: handle
      integer(c_int64_t), value :: n_terms
      integer(c_int), value :: max_level
      integer(c_int), intent(in) :: terms(n_terms*int(max_level, c_int64_t))
         !! One term per row, row-major
      integer(c_int) :: status

      type(fraglist_t), pointer :: list
      integer(default_int), allocatable :: buffer(:, :)
      type(error_t) :: error
      integer(int64) :: iterm, base
      integer :: irow

      status = MQC_BAD_HANDLE
      if (.not. c_associated(handle)) then
         last_message = "null fragment list handle"
         return
      end if
      call c_f_pointer(handle, list)

      if (n_terms < 0 .or. max_level < 1) then
         last_message = "mqc_fraglist_set: term count and level must be positive"
         status = MQC_ERROR
         return
      end if

      allocate (buffer(max(n_terms, 1_c_int64_t), max_level))
      buffer = 0
      do iterm = 1, n_terms
         base = (iterm - 1)*int(max_level, int64)
         do irow = 1, max_level
            buffer(iterm, irow) = int(terms(base + irow), default_int)
         end do
      end do

      call list%replace(buffer, int(n_terms, int64), int(max_level, default_int), error)
      deallocate (buffer)
      status = report(error)
   end function mqc_fraglist_set

   function mqc_fraglist_close_subsets(handle) result(status) &
      bind(C, name="mqc_fraglist_close_subsets")
      !! Add every proper subset the list is missing
      !!
      !! The other half of screening, and the half a caller does not think to
      !! ask for. An n-body term's delta is its energy less the delta of each
      !! proper subset, so a screen that keeps a trimer and drops one of its
      !! dimers has not approximated the expansion -- it has made one that
      !! cannot be assembled. Every route that supplies a list checks for this
      !! and refuses; this is what makes the list acceptable rather than
      !! merely diagnosing it.
      !!
      !! It only ever grows the list, so applying it twice does nothing the
      !! first did not, and applying it to a full list changes nothing at all.
      type(c_ptr), value :: handle
      integer(c_int) :: status

      type(fraglist_t), pointer :: list
      type(error_t) :: error

      status = MQC_BAD_HANDLE
      if (.not. c_associated(handle)) then
         last_message = "null fragment list handle"
         return
      end if
      call c_f_pointer(handle, list)

      call list%close_subsets(error)
      status = report(error)
   end function mqc_fraglist_close_subsets

   subroutine mqc_fraglist_last_error(buffer, buffer_len) &
      bind(C, name="mqc_fraglist_last_error")
      !! Copy the most recent failure message out as a C string
      integer(c_int), value :: buffer_len
      character(kind=c_char), intent(inout) :: buffer(buffer_len)

      integer :: n, i

      if (buffer_len <= 0) return
      n = min(len_trim(last_message), int(buffer_len) - 1)
      do i = 1, n
         buffer(i) = last_message(i:i)
      end do
      buffer(n + 1) = c_null_char
   end subroutine mqc_fraglist_last_error

   function report(error) result(status)
      !! Turn an error_t into a status, keeping the message for the caller
      type(error_t), intent(in) :: error
      integer(c_int) :: status

      if (error%has_error()) then
         last_message = error%get_message()
         status = MQC_ERROR
      else
         last_message = ""
         status = MQC_OK
      end if
   end function report

end module mqc_capi_fraglist
