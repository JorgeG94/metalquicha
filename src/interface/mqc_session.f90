!! MPI session for an externally driven run
module mqc_session
   !! Owns the MPI lifecycle when metalquicha is driven from outside rather
   !! than by `main`.
   !!
   !! **`begin` does not return on every rank, and that is the point.** Under
   !! `mpirun -np N python script.py` every rank starts an interpreter, but
   !! only rank 0 should execute the caller's script -- the others have no
   !! business deciding which fragments to compute, and running the setup N
   !! times would at best waste it and at worst have N ranks disagree. So on
   !! rank 0 `begin` returns and the caller carries on; on every other rank it
   !! enters a loop, waits for rank 0 to say what to do, and never comes back.
   !! A worker leaves that loop only to finalize and stop the process.
   !!
   !! The consequence to keep in mind: **anything the workers need must be
   !! sent to them.** They never read the input file, never see the caller's
   !! geometry, and never run a line of the driving language. Code that used to
   !! rely on every rank parsing the same file independently -- which is what
   !! `main` does -- does not work here.
   !!
   !! The command channel is one broadcast integer. That is deliberately the
   !! smallest thing that works: a worker blocked in `bcast` costs nothing, and
   !! the alternative of polling with `iprobe` would burn a core per rank for
   !! the whole setup phase.
   use pic_types, only: int32
   use pic_mpi_lib, only: comm_world, bcast, pic_mpi_init, pic_mpi_finalize
   use mqc_resources, only: resources_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: mqc_session_t
   public :: MQC_CMD_STOP, MQC_CMD_RUN

   integer(int32), parameter :: MQC_CMD_STOP = 0
      !! Leave the worker loop, finalize, and stop
   integer(int32), parameter :: MQC_CMD_RUN = 1
      !! Take part in a calculation rank 0 is about to drive

   type :: mqc_session_t
      !! A live MPI session
      logical :: active = .false.
      integer(int32) :: rank = -1
      integer(int32) :: n_ranks = 0
      type(resources_t) :: resources
         !! World and node communicators, in the shape `run_calculation` wants
   contains
      procedure :: begin => session_begin
      procedure :: command => session_command
      procedure :: end_session => session_end
      procedure :: is_root => session_is_root
   end type mqc_session_t

contains

   subroutine session_begin(this, error)
      !! Start MPI. Returns on rank 0; on any other rank, never returns.
      !!
      !! A worker spends the rest of the process inside `worker_loop` and exits
      !! from there. Callers on rank 0 get control back with a live session.
      class(mqc_session_t), intent(inout) :: this
      type(error_t), intent(inout) :: error

      if (this%active) then
         call error%set(ERROR_VALIDATION, "mqc_session: already started")
         return
      end if

      call pic_mpi_init()
      this%resources%mpi_comms%world_comm = comm_world()
      this%resources%mpi_comms%node_comm = this%resources%mpi_comms%world_comm%split()
      this%rank = this%resources%mpi_comms%world_comm%rank()
      this%n_ranks = this%resources%mpi_comms%world_comm%size()
      this%active = .true.

      if (this%rank /= 0) call worker_loop(this)
   end subroutine session_begin

   subroutine worker_loop(this)
      !! Wait for rank 0's instructions until told to stop, then stop
      !!
      !! This does not return to its caller. A worker that returned would carry
      !! on executing the driving script -- the very thing the session exists
      !! to prevent -- so the process ends here instead.
      type(mqc_session_t), intent(inout) :: this

      integer(int32) :: command

      do
         command = MQC_CMD_STOP
         call bcast(this%resources%mpi_comms%world_comm, command, 1_int32, 0_int32)

         select case (command)
         case (MQC_CMD_STOP)
            exit
         case (MQC_CMD_RUN)
            call worker_run(this)
         case default
            ! An unrecognised command can only mean rank 0 and the workers are
            ! running different code. Carrying on would deadlock at the next
            ! collective; stopping at least fails where the cause is.
            exit
         end select
      end do

      call this%resources%mpi_comms%node_comm%finalize()
      call this%resources%mpi_comms%world_comm%finalize()
      call pic_mpi_finalize()
      stop
   end subroutine worker_loop

   subroutine worker_run(this)
      !! A worker's half of a calculation
      !!
      !! Empty until the term list and geometry are broadcast: a worker cannot
      !! take part in a calculation it has not been told the shape of. It is
      !! separate from the loop already so that adding that is a change to one
      !! routine rather than to the dispatch.
      type(mqc_session_t), intent(inout) :: this

      ! Referenced so the argument is not merely decorative while this is a
      ! placeholder; the session is what the real body will work from.
      if (.not. this%active) return
   end subroutine worker_run

   subroutine session_command(this, command, error)
      !! Tell the workers what to do next. Rank 0 only.
      class(mqc_session_t), intent(inout) :: this
      integer(int32), intent(in) :: command
      type(error_t), intent(inout) :: error

      integer(int32) :: buffer

      if (.not. this%active) then
         call error%set(ERROR_VALIDATION, "mqc_session: no session to command")
         return
      end if
      if (this%rank /= 0) then
         ! Only rank 0 ever reaches here in a correct program; a worker issuing
         ! a command would leave every other rank waiting on a broadcast that
         ! never comes.
         call error%set(ERROR_VALIDATION, "mqc_session: only rank 0 may issue commands")
         return
      end if

      buffer = command
      call bcast(this%resources%mpi_comms%world_comm, buffer, 1_int32, 0_int32)
   end subroutine session_command

   subroutine session_end(this, error)
      !! Release the workers and shut MPI down
      class(mqc_session_t), intent(inout) :: this
      type(error_t), intent(inout) :: error

      if (.not. this%active) return

      ! The workers are blocked in a broadcast; this is what wakes them for the
      ! last time. Skipping it hangs every one of them.
      if (this%rank == 0 .and. this%n_ranks > 1) then
         call this%command(MQC_CMD_STOP, error)
         if (error%has_error()) return
      end if

      call this%resources%mpi_comms%node_comm%finalize()
      call this%resources%mpi_comms%world_comm%finalize()
      call pic_mpi_finalize()
      this%active = .false.
      this%rank = -1
      this%n_ranks = 0
   end subroutine session_end

   pure function session_is_root(this) result(root)
      !! Whether this rank is the one the caller drives from
      class(mqc_session_t), intent(in) :: this
      logical :: root
      root = (this%rank == 0)
   end function session_is_root

end module mqc_session
