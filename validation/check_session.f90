!! Manual check that an mqc session splits rank 0 from the workers
!!
!! There is no MPI test harness in ctest, so this is a standalone driver:
!!
!!     mpirun -np 4 ./build/check_session
!!
!! What it should print, in any order for the worker lines:
!!
!!     [rank 0] begin returned, 4 ranks
!!     [rank 0] issuing RUN
!!     [rank 0] issuing STOP and finalizing
!!     ok
!!
!! and nothing at all from ranks 1..3 -- their `begin` never returns, so they
!! never reach a print. If a worker line appears after "begin returned", the
!! asymmetric return is broken and every rank is running the driver script,
!! which is exactly the failure this arrangement exists to prevent.
!!
!! Run it at -np 1 too: with no workers there is nobody to broadcast to, and
!! the session must still start and stop cleanly rather than deadlocking on
!! its own collective.
!!
!! Expect one "IEEE_INVALID_FLAG IEEE_DIVIDE_BY_ZERO" note per WORKER on exit.
!! It is gfortran reporting raised floating-point flags at termination, and it
!! appears on the workers' path and not on rank 0's or on `mqc`'s. Nothing here
!! does arithmetic, so the flags come from below -- worth tracking down before
!! this is put in front of users, since a warning on every rank of a clean run
!! reads as a fault. Not suppressed, because clearing the flags here would also
!! clear genuine ones once a worker does real numerics.
program check_session
   use pic_types, only: int32
   use mqc_session, only: mqc_session_t, MQC_CMD_RUN
   use mqc_error, only: error_t
   implicit none

   type(mqc_session_t) :: session
   type(error_t) :: error

   call session%begin(error)
   if (error%has_error()) then
      write (*, "(A)") "FAILED to begin: "//error%get_message()
      stop 1
   end if

   ! Only rank 0 can be here. A worker is blocked inside begin and will leave
   ! it only to finalize and stop.
   if (.not. session%is_root()) then
      write (*, "(A)") "FAILED: a worker returned from begin"
      stop 1
   end if

   write (*, "(A,I0,A)") "[rank 0] begin returned, ", session%n_ranks, " ranks"

   write (*, "(A)") "[rank 0] issuing RUN"
   call session%command(MQC_CMD_RUN, error)
   if (error%has_error()) then
      write (*, "(A)") "FAILED to command: "//error%get_message()
      stop 1
   end if

   write (*, "(A)") "[rank 0] issuing STOP and finalizing"
   call session%end_session(error)
   if (error%has_error()) then
      write (*, "(A)") "FAILED to end: "//error%get_message()
      stop 1
   end if

   write (*, "(A)") "ok"
end program check_session
