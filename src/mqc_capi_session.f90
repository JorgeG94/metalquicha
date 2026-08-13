!! C entry points for the MPI session
module mqc_capi_session
   !! Starting and stopping MPI from outside, and the asymmetry that makes a
   !! driven run possible.
   !!
   !! **`mqc_session_begin` does not return on every rank.** On rank 0 it
   !! returns and the caller carries on; on every other rank it enters a loop,
   !! waits for rank 0, and never comes back -- it leaves only to finalize and
   !! stop the process. Under `mpirun -np N python script.py` that means N
   !! interpreters start, one runs the script, and the rest disappear into
   !! Fortran having executed nothing of consequence.
   !!
   !! A caller does not need to arrange any of that, but does need to know it,
   !! because two reasonable-looking habits break:
   !!
   !!   * **Anything before `begin` runs N times.** Imports, argument parsing,
   !!     reading a structure file -- all of it, on every rank, on a shared
   !!     filesystem. Call `begin` first.
   !!   * **Anything after `end` runs only on rank 0.** Which is usually what
   !!     is wanted, but it is not obvious from reading the script.
   !!
   !! There is one session per process and it is held here rather than handed
   !! out as a handle. MPI is process-wide, two of them would be a bug rather
   !! than a feature, and a handle would imply otherwise.
   use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
   use pic_types, only: int32
   use mqc_session, only: mqc_session_t
   use mqc_error, only: error_t
   implicit none
   private

   public :: mqc_session_begin, mqc_session_end
   public :: mqc_session_rank, mqc_session_size, mqc_session_is_root
   public :: mqc_session_active, mqc_session_last_error
   public :: current_session

   integer(c_int), parameter :: MQC_OK = 0
   integer(c_int), parameter :: MQC_FAIL = 1

   integer, parameter :: MESSAGE_LEN = 512
   character(len=MESSAGE_LEN), save :: last_message = ""

   type(mqc_session_t), save, target :: the_session
      !! The one session. Not exposed as a handle: MPI is process-wide, so a
      !! second one would be a bug that a handle-shaped interface would invite.

contains

   function current_session() result(session)
      !! The live session, for Fortran code that needs its communicators
      !!
      !! Not a C entry point -- `mqc_run` and anything else on this side of the
      !! boundary needs the resources the session holds, and this is how they
      !! reach them without another copy of the state.
      type(mqc_session_t), pointer :: session
      session => the_session
   end function current_session

   function mqc_session_begin() result(status) bind(C, name="mqc_session_begin")
      !! Start MPI. Returns on rank 0 only; workers never come back.
      integer(c_int) :: status

      type(error_t) :: error

      call the_session%begin(error)
      if (error%has_error()) then
         last_message = error%get_message()
         status = MQC_FAIL
         return
      end if
      status = MQC_OK
   end function mqc_session_begin

   function mqc_session_end() result(status) bind(C, name="mqc_session_end")
      !! Release the workers and shut MPI down
      !!
      !! Must be called. The workers are blocked waiting for exactly this, and
      !! a caller that exits without it leaves every other rank in the job
      !! hanging on a message that never arrives -- which on a batch system
      !! means burning the wall clock rather than failing.
      integer(c_int) :: status

      type(error_t) :: error

      call the_session%end_session(error)
      if (error%has_error()) then
         last_message = error%get_message()
         status = MQC_FAIL
         return
      end if
      status = MQC_OK
   end function mqc_session_end

   function mqc_session_rank() result(rank) bind(C, name="mqc_session_rank")
      !! This rank in the world communicator, or -1 without a session
      integer(c_int) :: rank
      rank = int(the_session%rank, c_int)
   end function mqc_session_rank

   function mqc_session_size() result(n) bind(C, name="mqc_session_size")
      !! Ranks in the job, or 0 without a session
      integer(c_int) :: n
      n = int(the_session%n_ranks, c_int)
   end function mqc_session_size

   function mqc_session_is_root() result(root) bind(C, name="mqc_session_is_root")
      !! Whether this is the rank the caller drives from
      !!
      !! Always true where a caller can observe it, since no other rank returns
      !! from `begin`. It exists so a script can assert that rather than assume
      !! it, and so a Fortran caller -- which can reach the workers' path -- can
      !! ask.
      integer(c_int) :: root

      root = 0_c_int
      if (the_session%is_root()) root = 1_c_int
   end function mqc_session_is_root

   function mqc_session_active() result(active) bind(C, name="mqc_session_active")
      !! Whether MPI is up
      integer(c_int) :: active

      active = 0_c_int
      if (the_session%active) active = 1_c_int
   end function mqc_session_active

   subroutine mqc_session_last_error(buffer_len, buffer) &
      bind(C, name="mqc_session_last_error")
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
   end subroutine mqc_session_last_error

end module mqc_capi_session
