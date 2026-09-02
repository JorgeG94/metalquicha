!! MPI session for an externally driven run
module mqc_session
   !! Owns the MPI lifecycle when metalquicha is driven from outside rather
   !! than by `main`.
   !!
   !! **`begin` does not return on every rank.** Under `mpirun -np N python
   !! script.py` every rank starts an interpreter, but only rank 0 executes the
   !! caller's script: on rank 0 `begin` returns, and on every other rank it
   !! enters a loop, waits for rank 0 to say what to do, and never comes back.
   !! A worker leaves that loop only to finalize and stop the process.
   !!
   !! The consequence to keep in mind: **anything the workers need must be
   !! sent to them.** They never read the input file, never see the caller's
   !! geometry, and never run a line of the driving language, so nothing here
   !! can rely on every rank parsing the same file the way `main` does.
   !!
   !! The command channel is one broadcast integer, which a worker blocks in
   !! for free.
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
   use mqc_calc_types, only: CALC_TYPE_OPTIMIZE, CALC_TYPE_CONFORMERS
   use mqc_method_factory, only: method_backend_built
   use mqc_method_types, only: method_type_to_string
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
      ! whenever there is more than one of them: the flags are already set on
      ! return from `pic_mpi_init` and are untouched at -np 1, so they are the
      ! MPI implementation's arithmetic and not ours. Cleared here rather than
      ! before the workers exit, so that a flag set later belongs to us.
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
      !! Does not return to its caller: a worker that returned would carry on
      !! executing the driving script, so the process ends here instead.
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
      !! calls the same `run_calculation` rank 0 does. **Keep this and
      !! `session_run` in step**: a broadcast added to one and not the other
      !! hangs the job rather than failing it.
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
      ! here has a document already known to be good. The settings travel as
      ! text and are parsed again rather than sent as a filled struct, so there
      ! is no field list to keep in step.
      call read_json_config_text(settings, config, error)
      if (error%has_error()) return
      ! Same count as rank 0, from the geometry that just arrived. See
      ! `session_run` for why leaving this out is not a crash.
      call config_to_driver(config, driver, n_fragments=sys_geom%n_monomers)
      call apply_log_level(config)

      ! No term list: the coordinator hands out fragment definitions as tasks
      ! and workers build them from sys_geom, so the list never leaves rank 0.
      dummy_terms = 0

      ! No result and no output: both are rank-0-only inside `run_calculation`
      ! and neither is collective.
      call run_calculation(this%resources, driver, sys_geom, sys_geom%bonds, &
                           supplied_terms=dummy_terms, n_supplied_terms=0_int64, &
                           write_output=.false.)
   end subroutine worker_run

   subroutine session_run(this, settings, sys_geom, terms, n_terms, label, &
                          write_output, result, error)
      !! Drive a calculation across the job. Rank 0 only.
      !!
      !! **The order here is load-bearing.** Everything that can fail on rank 0
      !! happens *before* `MQC_CMD_RUN` goes out, because once that command is
      !! sent the workers are committed to a set of broadcasts and a failure
      !! afterwards leaves N-1 ranks waiting for a geometry that never arrives.
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
      ! The fragment count comes from the system, not the document: a
      ! settings-only config reports zero fragments and the driver would
      ! compute the whole thing unfragmented. Both sides pass the same count
      ! because both take it from the geometry they hold.
      call config_to_driver(config, driver, n_fragments=sys_geom%n_monomers)
      ! Effective fragment potentials ride with the geometry, the settings
      ! document having no molecules block for them to live in. Copied after
      ! `config_to_driver`, which fills that field from a deck's molecules and
      ! leaves it unallocated when, as here, there are none.
      if (allocated(sys_geom%fragment_potentials)) then
         driver%fragment_potentials = sys_geom%fragment_potentials
      end if
      ! `main` sets the verbosity and a session never runs `main`, so without
      ! this the level in the settings is read, validated, and ignored.
      call apply_log_level(config)

      ! Everything from here to the command below is checked in the driver too,
      ! and the difference is what happens on failure: the driver aborts the
      ! job, which is right for `main` and would take down a caller's
      ! interpreter here. Refused before the command goes out, so the workers
      ! never learn anything was attempted.
      !
      ! A method whose backend this build does not carry reaches a factory with
      ! no branch for it and stops the process. Asked here rather than in the
      ! reader, so that parsing a deck means the same thing on every build.
      if (.not. method_backend_built(driver%method_config%method_type)) then
         call error%set(ERROR_VALIDATION, "method '"// &
                        method_type_to_string(driver%method_config%method_type)// &
                        "' needs the tblite library, and this build does not have "// &
                        "it. Reconfigure with -DMQC_ENABLE_TBLITE=ON.")
         return
      end if

      if (driver%calc_type == CALC_TYPE_CONFORMERS) then
         call error%set(ERROR_VALIDATION, "driver 'conformers' is not available "// &
                        "through a session or the C API. The sampling drives "// &
                        "run_calculation rather than being driven by it, so it "// &
                        "works from an input deck for a single molecule only. "// &
                        "Run the deck through the mqc executable.")
         return
      end if

      if (driver%calc_type == CALC_TYPE_OPTIMIZE) then
         call error%set(ERROR_VALIDATION, "driver 'optimize' is not available "// &
                        "through a session or the C API. The optimizer drives "// &
                        "run_calculation rather than being driven by it, so it "// &
                        "works from an input deck for a single molecule only. "// &
                        "Ask for 'energy' or 'gradient' here, or run the deck "// &
                        "through the mqc executable.")
         return
      end if

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
         ! A worker issuing a command would leave every other rank waiting on a
         ! broadcast that never comes.
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
