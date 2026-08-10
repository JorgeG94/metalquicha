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
   use, intrinsic :: ieee_exceptions, only: ieee_set_flag, ieee_all
   use pic_types, only: int32, int64
   use pic_mpi_lib, only: comm_world, bcast, pic_mpi_init, pic_mpi_finalize
   use mqc_resources, only: resources_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_bcast, only: bcast_string
   use mqc_bcast_system, only: bcast_system_geometry
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_config_types, only: mqc_config_t
   use mqc_json_config_reader, only: read_json_config_text
   use mqc_config_adapter, only: driver_config_t, config_to_driver, get_logger_level
   use mqc_driver, only: run_calculation
   use mqc_result_types, only: calculation_result_t
   use mqc_io_helpers, only: set_output_json_filename
   use mqc_validate, only: validate_system, validate_terms
   use pic_logger, only: logger => global_logger
   use mqc_fraglist, only: fraglist_t
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
      procedure :: run => session_run
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

      ! MPI start-up raises IEEE_INVALID and IEEE_DIVIDE_BY_ZERO on every rank
      ! whenever there is more than one of them -- measured, not guessed: the
      ! flags are already set on return from pic_mpi_init and are untouched at
      ! -np 1. They are the MPI implementation's arithmetic, not ours.
      !
      ! Clearing them here rather than before the workers exit is what keeps
      ! them useful. gfortran reports raised flags at STOP but not at END
      ! PROGRAM, so the noise surfaced only on the workers, and suppressing it
      ! at that end would have discarded genuine flags from whatever a worker
      ! computes. Cleared at the start, a flag set later belongs to us.
      call ieee_set_flag(ieee_all, .false.)

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
      !! Receives what rank 0 sends, in the order rank 0 sends it, and then
      !! calls the same `run_calculation` rank 0 does. The two sides are the
      !! same three steps written twice rather than shared, because they are
      !! mirror images and a shared routine would be mostly `if (is_root)`.
      !! Keep them in step: a broadcast added to one and not the other hangs
      !! the job rather than failing it.
      type(mqc_session_t), intent(inout) :: this

      character(len=:), allocatable :: settings
      type(system_geometry_t) :: sys_geom
      type(mqc_config_t) :: config
      type(driver_config_t) :: driver
      type(error_t) :: error
      integer :: dummy_terms(1, 1)

      if (.not. this%active) return

      call bcast_string(this%resources%mpi_comms%world_comm, settings, 0_int32)
      call bcast_system_geometry(this%resources%mpi_comms%world_comm, sys_geom, 0_int32)

      ! Rank 0 parsed this before it committed to the run, so a worker reaching
      ! here has a document already known to be good. Parsing it again rather
      ! than receiving a filled struct is the whole reason the settings travel
      ! as text: there is no field list to keep in step, so there is no field
      ! to forget, and a worker cannot end up running a different method than
      ! was asked for.
      call read_json_config_text(settings, config, error)
      if (error%has_error()) return
      ! Same count as rank 0, from the geometry that just arrived. See
      ! `session_run` for why leaving this out is not a crash.
      call config_to_driver(config, driver, n_fragments=sys_geom%n_monomers)
      call apply_log_level(config)

      ! No term list: the coordinator hands out fragment definitions as tasks
      ! and workers build them from sys_geom, so the list never leaves rank 0.
      dummy_terms = 0

      ! No result and no output. Both are rank-0-only inside `run_calculation`
      ! and neither is collective, so this asks for nothing rather than asking
      ! for something that would be discarded.
      call run_calculation(this%resources, driver, sys_geom, sys_geom%bonds, &
                           supplied_terms=dummy_terms, n_supplied_terms=0_int64, &
                           write_output=.false.)
   end subroutine worker_run

   subroutine session_run(this, settings, sys_geom, terms, n_terms, label, &
                          write_output, result, error)
      !! Drive a calculation across the job. Rank 0 only.
      !!
      !! The order here is load-bearing. Everything that can fail on rank 0 --
      !! parsing the settings, rejecting a bad document -- happens *before*
      !! `MQC_CMD_RUN` goes out. Once that command is sent the workers are
      !! committed to a set of broadcasts, and a rank 0 that discovered a
      !! problem afterwards could only return an error to its caller while N-1
      !! ranks sat waiting for a geometry that would never arrive.
      class(mqc_session_t), intent(inout) :: this
      character(len=*), intent(in) :: settings
         !! A JSON document of everything but the molecules
      type(system_geometry_t), intent(inout) :: sys_geom
         !! Read, not modified; `intent(inout)` because the broadcast it is
         !! passed to writes on the receiving ranks.
      integer, intent(in) :: terms(:, :)
         !! (n_terms, max_level), 1-based monomer indices, zero-padded
      integer(int64), intent(in) :: n_terms
         !! Rows of `terms` to use; 0 to let the driver generate its own list
      character(len=*), intent(in) :: label
         !! Names the output files. Empty leaves the default.
      logical, intent(in) :: write_output
      type(calculation_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: payload
      type(mqc_config_t) :: config
      type(driver_config_t) :: driver

      if (.not. this%active) then
         call error%set(ERROR_VALIDATION, "mqc_session: no session to run in")
         return
      end if
      if (this%rank /= 0) then
         call error%set(ERROR_VALIDATION, "mqc_session: only rank 0 may drive a calculation")
         return
      end if

      call read_json_config_text(settings, config, error)
      if (error%has_error()) return
      ! The fragment count comes from the system, not the document. Without it
      ! a settings-only config reports zero fragments and the driver computes
      ! the whole thing unfragmented -- converged, plausible, and not what was
      ! asked for. Both sides must pass the same count, and they do because
      ! both take it from the geometry they hold.
      call config_to_driver(config, driver, n_fragments=sys_geom%n_monomers)
      ! `main` sets the verbosity and a session never runs `main`, so without
      ! this the level in the settings is read, validated, and ignored -- and
      ! the default is loud enough that a script running a few thousand
      ! fragments gets a few thousand lines it did not ask for.
      call apply_log_level(config)

      ! Checked here as well as in the driver, and the difference is what
      ! happens on failure. The driver aborts the job, which is right for
      ! `main` -- there is no one to tell -- and wrong for a driven run, where
      ! it would take down the caller's interpreter and everything it had
      ! built. Checking before the command goes out means a bad system or a
      ! bad list comes back as an error the caller can act on, and the workers
      ! never learn anything was attempted. The driver's copy stays as the
      ! backstop for every other entry path.
      call validate_system(sys_geom,.not. config%unchecked_input, error, &
                           check_bonds=allocated(sys_geom%bonds))
      if (error%has_error()) return
      if (n_terms > 0) then
         call check_supplied_terms(terms, n_terms, sys_geom, &
                                   .not. config%unchecked_input, error)
         if (error%has_error()) return
      end if

      call this%command(MQC_CMD_RUN, error)
      if (error%has_error()) return

      payload = settings
      call bcast_string(this%resources%mpi_comms%world_comm, payload, 0_int32)
      call bcast_system_geometry(this%resources%mpi_comms%world_comm, sys_geom, 0_int32)

      if (len_trim(label) > 0) call set_output_json_filename(trim(label))

      if (n_terms > 0) then
         call run_calculation(this%resources, driver, sys_geom, sys_geom%bonds, &
                              result_out=result, supplied_terms=terms, &
                              n_supplied_terms=n_terms, write_output=write_output)
      else
         call run_calculation(this%resources, driver, sys_geom, sys_geom%bonds, &
                              result_out=result, write_output=write_output)
      end if
   end subroutine session_run

   subroutine apply_log_level(config)
      !! Set the logger verbosity the settings asked for
      type(mqc_config_t), intent(in) :: config

      if (allocated(config%log_level)) then
         call logger%configure(get_logger_level(config%log_level))
      end if
   end subroutine apply_log_level

   subroutine check_supplied_terms(terms, n_terms, sys_geom, strict, error)
      !! Run a supplied term list past the validator without committing to it
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: terms(:, :)
      integer(int64), intent(in) :: n_terms
      logical, intent(in) :: strict
      type(error_t), intent(inout) :: error

      type(fraglist_t) :: list

      call list%replace(terms, n_terms, size(terms, 2), error)
      if (.not. error%has_error()) then
         call validate_terms(list, sys_geom, strict, error)
      end if
      call list%destroy()
   end subroutine check_supplied_terms

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
