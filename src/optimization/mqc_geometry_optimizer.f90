!! Drive a geometry optimization over any method this program can gradient
module mqc_geometry_optimizer
   !! An optimization is a loop *over* calculations, so it sits above
   !! `run_calculation` rather than inside it. That is also why the dispatch is
   !! in `main` and not in `mqc_driver`: the driver would otherwise have to use
   !! this module while this module uses the driver, which Fortran will not
   !! compile.
   !!
   !! **What makes this backend-agnostic.** The evaluator below asks
   !! `run_calculation` for a gradient and takes back a full-system one. Which
   !! engine produced it -- tblite, libcint on the CPU, cuEST on a GPU -- and
   !! whether it came from one calculation or from an MBE over twenty thousand
   !! fragments are both settled inside that call. Nothing here, and nothing in
   !! the optimizer engine, has a branch on the method.
   !!
   !! **MPI shape.** Rank 0 runs the optimizer and decides every step; the other
   !! ranks sit in `worker_loop` waiting to be told a geometry. Steps are never
   !! recomputed per rank -- two ranks that disagreed in the last bits of a
   !! coordinate would build different fragments and the job would hang rather
   !! than fail. This mirrors `mqc_session`, which splits rank 0 from its
   !! workers the same way and for the same reason.
   use pic_types, only: dp, int32, int64
   use pic_mpi_lib, only: comm_t, bcast
   use pic_logger, only: logger => global_logger, warning_level, verbose_level
   use pic_io, only: to_char
   use mqc_optimizer_types, only: optimizer_settings_t, &
                                  OPT_COORDS_CARTESIAN, OPT_COORDS_UNKNOWN, OPT_ALGO_UNKNOWN, &
                                  coordinates_to_string, algorithm_to_string
   use mqc_physical_fragment, only: system_geometry_t, to_angstrom
   use mqc_config_adapter, only: driver_config_t
   use mqc_config_types, only: bond_t
   use mqc_resources, only: resources_t
   use mqc_result_types, only: calculation_result_t
   use mqc_calc_types, only: CALC_TYPE_GRADIENT, CALC_TYPE_ENERGY
   use mqc_error, only: error_t, ERROR_GENERIC, ERROR_VALIDATION
   use mqc_elements, only: element_number_to_symbol
   use mqc_dlfind_bridge, only: dlfind_available, dlfind_optimize
   use mqc_optimization_output, only: optimization_record_t, write_optimization_json, &
                                      write_trajectory_xyz
   use mqc_frag_utils, only: generate_mbe_term_list
   implicit none
   private

   public :: run_geometry_optimization

   ! What rank 0 tells the workers to do next. One broadcast integer, for the
   ! reason mqc_session gives: a worker blocked in bcast costs nothing.
   integer(int32), parameter :: OPT_CMD_STOP = 0
   integer(int32), parameter :: OPT_CMD_EVALUATE = 1
   integer(int32), parameter :: OPT_CMD_FINAL = 2

   ! The optimization in progress, reachable from the evaluator.
   !
   ! Module state, and not because it is convenient: the engine calls the
   ! evaluator through a plain procedure interface with no argument to carry a
   ! context, so there is nowhere else for the geometry and the settings to
   ! live. One optimization at a time, which is what `active` guards -- a
   ! second one entered while the first is running would silently take over the
   ! first one's system.
   logical, save :: active = .false.
   type(resources_t), save :: ctx_resources
   type(driver_config_t), save :: ctx_config
   type(system_geometry_t), save :: ctx_sys_geom
   type(bond_t), allocatable, save :: ctx_bonds(:)
   logical, save :: ctx_has_bonds = .false.
   integer, save :: ctx_n_evaluations = 0
   logical, save :: ctx_failed = .false.
   logical, save :: ctx_probing = .false.
      !! True during the capability probe, which is an evaluation the user did
      !! not ask for and should not see numbered among their steps.
   ! The path taken, grown a step at a time. Held until the run ends so the
   ! whole thing can be written as one document rather than a file appended to
   ! a hundred times, which is also what makes `trajectory: false` worth having
   ! for a large system.
   real(dp), allocatable, save :: ctx_trajectory(:, :, :)
   real(dp), allocatable, save :: ctx_trajectory_energies(:)
   integer, save :: ctx_n_steps = 0
   ! The term list, chosen once and kept. See `freeze_term_list`.
   integer, allocatable, save :: ctx_terms(:, :)
   integer(int64), save :: ctx_n_terms = 0
   logical, save :: ctx_terms_frozen = .false.
   real(dp), save :: ctx_last_gradient_max = huge(1.0_dp)
      !! Largest gradient component of the last successful step, which is what
      !! says whether this converged. DL-FIND reports geometries but not a
      !! converged flag, and running out of steps looks from the outside
      !! exactly like finishing -- see `report_result`.

contains

   subroutine run_geometry_optimization(resources, config, sys_geom, bonds, error)
      !! Optimize `sys_geom` in place. Every rank calls this.
      !!
      !! On return rank 0 holds the optimized coordinates in `sys_geom`; the
      !! other ranks hold the last geometry they were asked to compute, which
      !! is the same thing for every step but the bookkeeping of the last one.
      !! Rank 0 is the only rank whose copy should be written out.
      type(resources_t), intent(in) :: resources
      type(driver_config_t), intent(in) :: config
      type(system_geometry_t), intent(inout) :: sys_geom
      type(bond_t), intent(in), optional :: bonds(:)
      type(error_t), intent(inout) :: error

      integer :: n_atoms
      integer :: rank
      integer(int32) :: refused
      logical :: converged
      integer :: probe_status
      real(dp) :: final_energy
      real(dp), allocatable :: coords(:, :)
      integer, allocatable :: atomic_numbers(:)
      integer, allocatable :: residues(:)

      rank = resources%mpi_comms%world_comm%rank()
      n_atoms = sys_geom%total_atoms

      if (active) then
         call error%set(ERROR_GENERIC, &
                        "run_geometry_optimization: an optimization is already running")
         return
      end if

      ! The engine has to exist before anything is broadcast: a refusal
      ! discovered after the workers were told to expect a geometry would
      ! leave them waiting for one that never comes. Same ordering rule as
      ! mqc_session%run.
      if (rank == 0) then
         if (.not. dlfind_available()) then
            call error%set(ERROR_VALIDATION, &
                           'driver "Optimize" needs the DL-FIND backend; this build has none. '// &
                           "Configure with -DMQC_ENABLE_DLFIND=ON.")
         else if (config%optimization%coordinates == OPT_COORDS_UNKNOWN) then
            ! A backstop for the callers that convert a config without passing
            ! an error to `config_to_driver`. Reaching here means the spelling
            ! was refused there and nobody was listening; optimizing in
            ! whatever the default happened to be would be worse.
            call error%set(ERROR_VALIDATION, &
                           "Unknown keywords.optimization.coordinates. Use cartesian, hdlc or dlc.")
         else if (config%optimization%algorithm == OPT_ALGO_UNKNOWN) then
            call error%set(ERROR_VALIDATION, &
                           "Unknown keywords.optimization.algorithm. Use lbfgs, cg, sd or prfo.")
         else if (n_atoms < 2) then
            call error%set(ERROR_VALIDATION, &
                           "Nothing to optimize: a single atom has no geometry.")
         end if
      end if

      ! Every rank has to learn the answer, not just the one that worked it
      ! out. A rank 0 that refused and returned on its own would leave the
      ! workers heading into `worker_loop` to block on a broadcast that is
      ! never coming -- the job would hang instead of printing the message
      ! that was already sitting in `error`.
      refused = 0_int32
      if (rank == 0 .and. error%has_error()) refused = 1_int32
      call bcast(resources%mpi_comms%world_comm, refused, 1_int32, 0_int32)
      if (refused /= 0_int32) return

      call begin_context(resources, config, sys_geom, bonds)

      if (rank == 0) then
         call report_settings(config, n_atoms)
         call freeze_term_list(sys_geom, config)

         allocate (coords(3, n_atoms), source=sys_geom%coordinates)
         allocate (atomic_numbers(n_atoms), source=sys_geom%element_numbers)
         allocate (residues(n_atoms))
         call build_residues(sys_geom, residues)

         ! One gradient before handing over, to find out whether this method
         ! can produce one at all.
         !
         ! Not every method here can: the CPU ab initio backend computes
         ! energies and refuses gradients, and there is no way to ask in
         ! advance -- the refusal is raised where the gradient would have been
         ! built. Without this probe that refusal reaches DL-FIND as a failed
         ! first step, and DL-FIND responds by calling dlf_error, which ends
         ! the process with a Fortran backtrace. The message was correct and
         ! the presentation was a crash.
         !
         ! It costs one evaluation on a run that works. That is a few percent
         ! of a short optimization and nothing on a long one, against turning
         ! an unsupported method from a backtrace into one line.
         call probe_gradient(n_atoms, coords, probe_status)
         if (probe_status /= 0) then
            call error%set(ERROR_VALIDATION, &
                           "This method cannot produce a gradient, so it cannot be "// &
                           "optimized. The message above says why.")
            call stop_workers(resources%mpi_comms%world_comm)
            call end_context()
            return
         end if

         call dlfind_optimize(config%optimization, n_atoms, atomic_numbers, residues, &
                              coords, evaluate_energy_gradient, record_step, &
                              final_energy, error)

         ! Stop the workers whether or not that succeeded. A rank 0 that
         ! returned an error without sending this would leave N-1 ranks
         ! blocked in bcast for the rest of the job.
         call stop_workers(resources%mpi_comms%world_comm)

         if (.not. error%has_error()) then
            sys_geom%coordinates = coords
            converged = ctx_last_gradient_max <= config%optimization%gradient_tolerance
            call write_final_single_point(sys_geom)
            call report_result(final_energy, converged, config%optimization%gradient_tolerance)
            ! Written either way. A run that ran out of steps left the geometry
            ! closer to the minimum than it started, and that is exactly what
            ! someone wants to restart from.
            call write_optimized_xyz(sys_geom, final_energy, converged, error)
            if (.not. error%has_error()) then
               call write_optimization_record(sys_geom, config, final_energy, converged, error)
            end if
            if (.not. converged .and. .not. error%has_error()) then
               call error%set(ERROR_GENERIC, &
                              "Geometry optimization did not converge in "// &
                              to_char(config%optimization%max_steps)//" steps. "// &
                              "Largest gradient component "// &
                              to_char(ctx_last_gradient_max)//" Hartree/Bohr against a "// &
                              "target of "//to_char(config%optimization%gradient_tolerance)// &
                              ". The last geometry was still written.")
            end if
         end if
      else
         call worker_loop(resources%mpi_comms%world_comm)
      end if

      call end_context()

   end subroutine run_geometry_optimization

   subroutine evaluate_energy_gradient(n_atoms, coords, energy, gradient, status)
      !! One energy and gradient, at the geometry the engine asked for
      !!
      !! Rank 0 only -- this is the engine's callback. It brings the workers to
      !! the same geometry and then every rank enters `run_calculation`
      !! together, which is what the fragmented path requires and what the
      !! unfragmented path tolerates (the workers no-op inside it).
      use mqc_driver, only: run_calculation
      integer, intent(in) :: n_atoms
      real(dp), intent(in) :: coords(3, n_atoms)
      real(dp), intent(out) :: energy
      real(dp), intent(out) :: gradient(3, n_atoms)
      integer, intent(out) :: status

      type(calculation_result_t) :: result
      type(driver_config_t) :: gradient_config
      integer :: saved_level

      energy = 0.0_dp
      gradient = 0.0_dp
      status = 0

      ctx_sys_geom%coordinates = coords
      call send_geometry(ctx_resources%mpi_comms%world_comm, ctx_sys_geom)

      ! The deck said "Optimize"; each step of one is a gradient. Copied and
      ! overridden rather than mutated, so the caller's config still says what
      ! the user asked for when this is over.
      gradient_config = ctx_config
      gradient_config%calc_type = CALC_TYPE_GRADIENT

      ! No files per step. A hundred-step optimization would otherwise write a
      ! hundred output documents over the top of each other, and the fragment
      ! breakdown of an intermediate geometry is not something anyone wants.
      !
      ! And no single-point report per step either. `run_calculation` ends by
      ! logging its energy, dipole, HOMO-LUMO gap and gradient norm, which is
      ! the right thing to say about a calculation somebody asked for and the
      ! wrong thing to say a hundred times about the steps of an optimization.
      ! Lowered to warnings so an unconverged fragment still speaks up, and
      ! left alone entirely when the user asked for verbose -- at that point
      ! they want the per-step detail.
      saved_level = logger%log_level
      if (saved_level < verbose_level) call logger%configure(level=warning_level)
      call run_step(gradient_config, result)
      call logger%configure(level=saved_level)

      if (.not. ctx_probing) ctx_n_evaluations = ctx_n_evaluations + 1

      ! A fragment that failed leaves the driver going -- it reports every bad
      ! term at the end rather than the first -- so the result carries the
      ! failure and the gradient array carries whatever had accumulated. Handing
      ! that to an optimizer is the worst shape a failure can take: it is a
      ! perfectly plausible small gradient and the step taken on it looks fine.
      if (result%has_error) then
         call logger%error(step_label()//" failed: "//result%error%get_message())
         ctx_failed = .true.
         status = 1
         return
      end if

      if (.not. result%has_gradient) then
         call logger%error(step_label()//" returned no gradient. The method named "// &
                           "in this deck cannot produce one.")
         ctx_failed = .true.
         status = 1
         return
      end if

      energy = result%energy%total()
      gradient = result%gradient
      ctx_last_gradient_max = maxval(abs(gradient))

      if (.not. ctx_probing) then
         call logger%info("  step "//to_char(ctx_n_evaluations)// &
                          "  E = "//to_char(energy)// &
                          "  |g|max = "//to_char(maxval(abs(gradient))))
      end if

   end subroutine evaluate_energy_gradient

   subroutine worker_loop(comm)
      !! A worker's half of an optimization: compute what rank 0 asks for
      !!
      !! Keep this in step with `evaluate_energy_gradient` above. The two are
      !! mirror images written twice rather than shared, as `mqc_session` does
      !! it -- a broadcast added to one and not the other hangs the job instead
      !! of failing it.
      use mqc_driver, only: run_calculation
      type(comm_t), intent(in) :: comm

      type(calculation_result_t) :: result
      type(driver_config_t) :: gradient_config, step_config
      integer(int32) :: command
      integer :: saved_level

      gradient_config = ctx_config
      gradient_config%calc_type = CALC_TYPE_GRADIENT

      do
         command = OPT_CMD_STOP
         call bcast(comm, command, 1_int32, 0_int32)
         if (command == OPT_CMD_STOP) exit

         call receive_geometry(comm, ctx_sys_geom)

         ! The last round is an energy rather than a gradient, and both sides
         ! have to agree on which: the calculation type decides what each
         ! fragment computes, so a worker still doing gradients while rank 0
         ! collects energies produces a result assembled from two different
         ! calculations.
         step_config = gradient_config
         if (command == OPT_CMD_FINAL) step_config%calc_type = CALC_TYPE_ENERGY

         saved_level = logger%log_level
         if (saved_level < verbose_level) call logger%configure(level=warning_level)
         call run_step(step_config, result)
         call logger%configure(level=saved_level)
      end do

   end subroutine worker_loop

   subroutine send_final_geometry(comm, sys_geom)
      !! Move a geometry to the workers, with the command already sent
      !!
      !! Separate from `send_geometry` only in that the caller has already
      !! chosen and broadcast the command -- EVALUATE or FINAL. Splitting the
      !! two keeps the "one broadcast integer then the coordinates" order in
      !! one place, which is the order `worker_loop` reads them in.
      type(comm_t), intent(in) :: comm
      type(system_geometry_t), intent(in) :: sys_geom

      real(dp), allocatable :: buffer(:)
      integer(int32) :: n

      ! Flat, because pic-mpi broadcasts rank-1 arrays and the coordinates are
      ! rank-2. reshape rather than a pointer remap: this is 3N doubles once
      ! per step, against an SCF, and the copy is not measurable.
      n = int(3*sys_geom%total_atoms, int32)
      allocate (buffer(n))
      buffer = reshape(sys_geom%coordinates, [3*sys_geom%total_atoms])
      call bcast(comm, buffer, n, 0_int32)
      deallocate (buffer)

   end subroutine send_final_geometry

   subroutine send_geometry(comm, sys_geom)
      !! Rank 0's half of an ordinary optimization step
      type(comm_t), intent(in) :: comm
      type(system_geometry_t), intent(in) :: sys_geom

      integer(int32) :: command

      command = OPT_CMD_EVALUATE
      call bcast(comm, command, 1_int32, 0_int32)
      call send_final_geometry(comm, sys_geom)

   end subroutine send_geometry

   subroutine receive_geometry(comm, sys_geom)
      !! A worker's half of moving a geometry from rank 0
      type(comm_t), intent(in) :: comm
      type(system_geometry_t), intent(inout) :: sys_geom

      real(dp), allocatable :: buffer(:)
      integer(int32) :: n

      n = int(3*sys_geom%total_atoms, int32)
      allocate (buffer(n))
      call bcast(comm, buffer, n, 0_int32)
      sys_geom%coordinates = reshape(buffer, [3, sys_geom%total_atoms])
      deallocate (buffer)

   end subroutine receive_geometry

   subroutine stop_workers(comm)
      !! Release the workers from `worker_loop`
      type(comm_t), intent(in) :: comm

      integer(int32) :: command

      command = OPT_CMD_STOP
      call bcast(comm, command, 1_int32, 0_int32)

   end subroutine stop_workers

   subroutine begin_context(resources, config, sys_geom, bonds)
      !! Put this optimization where the evaluator can reach it
      type(resources_t), intent(in) :: resources
      type(driver_config_t), intent(in) :: config
      type(system_geometry_t), intent(in) :: sys_geom
      type(bond_t), intent(in), optional :: bonds(:)

      ctx_resources = resources
      ctx_config = config
      ctx_sys_geom = sys_geom
      ctx_has_bonds = present(bonds)
      if (ctx_has_bonds) then
         if (allocated(ctx_bonds)) deallocate (ctx_bonds)
         allocate (ctx_bonds, source=bonds)
      end if
      ctx_n_evaluations = 0
      ctx_failed = .false.
      ctx_last_gradient_max = huge(1.0_dp)
      ctx_n_steps = 0
      if (allocated(ctx_trajectory)) deallocate (ctx_trajectory)
      if (allocated(ctx_trajectory_energies)) deallocate (ctx_trajectory_energies)
      if (allocated(ctx_terms)) deallocate (ctx_terms)
      ctx_n_terms = 0
      ctx_terms_frozen = .false.
      active = .true.

   end subroutine begin_context

   subroutine end_context()
      !! Let go of the system, so the next optimization starts clean
      if (allocated(ctx_bonds)) deallocate (ctx_bonds)
      if (allocated(ctx_trajectory)) deallocate (ctx_trajectory)
      if (allocated(ctx_trajectory_energies)) deallocate (ctx_trajectory_energies)
      if (allocated(ctx_terms)) deallocate (ctx_terms)
      ctx_terms_frozen = .false.
      call ctx_sys_geom%destroy()
      ctx_has_bonds = .false.
      active = .false.

   end subroutine end_context

   subroutine write_final_single_point(sys_geom)
      !! One energy at the optimized geometry, written out properly
      !!
      !! Every step runs with `write_output=.false.`, so without this an
      !! optimization produces no output document at all -- the energy exists
      !! only in the log, and anything downstream that reads
      !! `output_<name>.json` (a driving script, run_validation.py) finds
      !! nothing. One more energy against a hundred gradients is not a cost
      !! worth avoiding, and it also gets the dipole and the HOMO-LUMO gap
      !! reported at the geometry they belong to rather than the input one.
      use mqc_driver, only: run_calculation
      type(system_geometry_t), intent(in) :: sys_geom

      type(calculation_result_t) :: result
      type(driver_config_t) :: energy_config
      integer(int32) :: command

      command = OPT_CMD_FINAL
      call bcast(ctx_resources%mpi_comms%world_comm, command, 1_int32, 0_int32)

      ctx_sys_geom%coordinates = sys_geom%coordinates
      call send_final_geometry(ctx_resources%mpi_comms%world_comm, ctx_sys_geom)

      energy_config = ctx_config
      energy_config%calc_type = CALC_TYPE_ENERGY

      call run_step(energy_config, result, write_output=.true.)

   end subroutine write_final_single_point

   subroutine record_step(n_atoms, coords, energy)
      !! Keep one accepted geometry, growing the trajectory to fit
      !!
      !! Doubling rather than one reallocation per step: an optimization of a
      !! few hundred steps would otherwise copy the whole path every time it
      !! took another one.
      integer, intent(in) :: n_atoms
      real(dp), intent(in) :: coords(3, n_atoms)
      real(dp), intent(in) :: energy

      real(dp), allocatable :: bigger(:, :, :)
      real(dp), allocatable :: bigger_energies(:)
      integer :: capacity, new_capacity

      ctx_n_steps = ctx_n_steps + 1
      if (.not. ctx_config%optimization%trajectory) return

      if (.not. allocated(ctx_trajectory)) then
         allocate (ctx_trajectory(3, n_atoms, 16))
         allocate (ctx_trajectory_energies(16))
      end if

      capacity = size(ctx_trajectory, 3)
      if (ctx_n_steps > capacity) then
         new_capacity = 2*capacity
         allocate (bigger(3, n_atoms, new_capacity))
         allocate (bigger_energies(new_capacity))
         bigger(:, :, 1:capacity) = ctx_trajectory
         bigger_energies(1:capacity) = ctx_trajectory_energies
         call move_alloc(bigger, ctx_trajectory)
         call move_alloc(bigger_energies, ctx_trajectory_energies)
      end if

      ctx_trajectory(:, :, ctx_n_steps) = coords
      ctx_trajectory_energies(ctx_n_steps) = energy

   end subroutine record_step

   subroutine write_optimization_record(sys_geom, config, final_energy, converged, error)
      !! Write what the optimization did, for a reader that is not a person
      !!
      !! Both files, or neither: a trajectory without the record that says
      !! whether it converged is a path to an unknown place.
      use mqc_io_helpers, only: get_output_json_filename
      type(system_geometry_t), intent(in) :: sys_geom
      type(driver_config_t), intent(in) :: config
      real(dp), intent(in) :: final_energy
      logical, intent(in) :: converged
      type(error_t), intent(inout) :: error

      type(optimization_record_t) :: record

      record%converged = converged
      record%n_steps = ctx_n_steps
      record%n_evaluations = ctx_n_evaluations
      record%final_energy = final_energy
      record%final_gradient_max = ctx_last_gradient_max
      record%gradient_tolerance = config%optimization%gradient_tolerance
      record%coordinates = coordinates_to_string(config%optimization%coordinates)
      record%algorithm = algorithm_to_string(config%optimization%algorithm)
      record%element_numbers = sys_geom%element_numbers
      record%final_coordinates = sys_geom%coordinates

      ! Trimmed to the steps actually taken. The array is grown by doubling, so
      ! its last frames are whatever the previous allocation held -- writing
      ! those out would put geometries in the trajectory that the optimization
      ! never visited.
      if (allocated(ctx_trajectory) .and. ctx_n_steps > 0) then
         record%trajectory = ctx_trajectory(:, :, 1:ctx_n_steps)
         record%trajectory_energies = ctx_trajectory_energies(1:ctx_n_steps)
      end if

      call write_optimization_json(optimization_path("_optimization.json"), record, error)
      if (error%has_error()) return
      call write_trajectory_xyz(optimization_path("_trajectory.xyz"), record, error)

   end subroutine write_optimization_record

   function optimization_path(suffix) result(path)
      !! A file beside the output document, named for it
      use mqc_io_helpers, only: get_output_json_filename
      character(len=*), intent(in) :: suffix
      character(len=:), allocatable :: path

      character(len=:), allocatable :: json_name
      integer :: dot

      json_name = get_output_json_filename()
      dot = index(json_name, ".", back=.true.)
      if (dot > 1) then
         path = json_name(1:dot - 1)//suffix
      else
         path = trim(json_name)//suffix
      end if

   end function optimization_path

   function step_label() result(text)
      !! How to name the evaluation that just failed
      !!
      !! The probe is not a step and numbering it as one reads as "step 0",
      !! which invites the question of what step zero is.
      character(len=:), allocatable :: text

      if (ctx_probing) then
         text = "The initial gradient check"
      else
         text = "Optimization step "//to_char(ctx_n_evaluations)
      end if
   end function step_label

   subroutine probe_gradient(n_atoms, coords, status)
      !! Ask for one gradient, and report only whether it worked
      !!
      !! Collective in the same way an ordinary step is: rank 0 calls this and
      !! the workers service it from `worker_loop`, which they are already
      !! sitting in.
      integer, intent(in) :: n_atoms
      real(dp), intent(in) :: coords(3, n_atoms)
      integer, intent(out) :: status

      real(dp) :: energy
      real(dp), allocatable :: gradient(:, :)

      allocate (gradient(3, n_atoms))
      ctx_probing = .true.
      call evaluate_energy_gradient(n_atoms, coords, energy, gradient, status)
      ctx_probing = .false.
      deallocate (gradient)

   end subroutine probe_gradient

   subroutine run_step(step_config, result, write_output)
      !! One calculation at the current geometry, on the frozen term list
      !!
      !! Every rank calls this, and every rank must reach it with the same
      !! configuration -- the four combinations of bonds and frozen terms are
      !! written out here once rather than at each of the three call sites,
      !! which is where they had already started to drift.
      use mqc_driver, only: run_calculation
      type(driver_config_t), intent(in) :: step_config
      type(calculation_result_t), intent(out) :: result
      logical, intent(in), optional :: write_output

      logical :: wants_output

      wants_output = .false.
      if (present(write_output)) wants_output = write_output

      if (ctx_terms_frozen) then
         if (ctx_has_bonds) then
            call run_calculation(ctx_resources, step_config, ctx_sys_geom, ctx_bonds, &
                                 result, write_output=wants_output, &
                                 supplied_terms=ctx_terms, n_supplied_terms=ctx_n_terms)
         else
            call run_calculation(ctx_resources, step_config, ctx_sys_geom, &
                                 result_out=result, write_output=wants_output, &
                                 supplied_terms=ctx_terms, n_supplied_terms=ctx_n_terms)
         end if
      else
         if (ctx_has_bonds) then
            call run_calculation(ctx_resources, step_config, ctx_sys_geom, ctx_bonds, &
                                 result, write_output=wants_output)
         else
            call run_calculation(ctx_resources, step_config, ctx_sys_geom, &
                                 result_out=result, write_output=wants_output)
         end if
      end if

   end subroutine run_step

   subroutine build_residues(sys_geom, residues)
      !! Group the atoms the way an internal-coordinate scheme needs them
      !!
      !! A residue is a group of atoms that an optimizer may describe with
      !! bonds, angles and torsions among themselves. This program already
      !! knows that grouping -- it is the monomers -- so HDLC gets the
      !! chemically right answer for free rather than having to perceive it.
      !!
      !! Everything in one residue when there is only one monomer, which is
      !! what an ordinary molecule wants.
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(out) :: residues(:)

      integer :: imon, k, atom

      residues = 1
      if (sys_geom%n_monomers <= 1) return

      if (allocated(sys_geom%fragment_atoms) .and. allocated(sys_geom%fragment_sizes)) then
         do imon = 1, sys_geom%n_monomers
            do k = 1, sys_geom%fragment_sizes(imon)
               ! fragment_atoms is 0-based, as the deck writes it
               atom = sys_geom%fragment_atoms(k, imon) + 1
               if (atom >= 1 .and. atom <= size(residues)) residues(atom) = imon
            end do
         end do
      else if (sys_geom%atoms_per_monomer > 0) then
         do atom = 1, size(residues)
            residues(atom) = (atom - 1)/sys_geom%atoms_per_monomer + 1
         end do
      end if

   end subroutine build_residues

   subroutine report_settings(config, n_atoms)
      !! Say what is about to be optimized, and how
      type(driver_config_t), intent(in) :: config
      integer, intent(in) :: n_atoms

      call logger%info(" ")
      call logger%info("============================================")
      call logger%info("Geometry optimization")
      call logger%info("  Atoms: "//to_char(n_atoms))
      call logger%info("  Coordinates: "//coordinates_to_string(config%optimization%coordinates))
      call logger%info("  Algorithm: "//algorithm_to_string(config%optimization%algorithm))
      call logger%info("  Max steps: "//to_char(config%optimization%max_steps))
      call logger%info("  Gradient tolerance: "// &
                       to_char(config%optimization%gradient_tolerance)//" Hartree/Bohr")
      if (config%nlevel > 0) then
         call logger%info("  Energy and gradient from MBE("//to_char(config%nlevel)//")")
      end if
      call logger%info("============================================")

   end subroutine report_settings

   subroutine freeze_term_list(sys_geom, config)
      !! Choose the MBE term list once, here, and keep it for the whole run
      !!
      !! **The problem.** A fragmented run generates its term list from the
      !! geometry it was handed, screening out the n-mers whose monomers are
      !! further apart than the cutoff. Regenerate that at every step of an
      !! optimization and the list changes as the molecules move: a dimer
      !! crossing the cutoff enters or leaves the sum, and the total energy
      !! jumps by that dimer's whole interaction. The optimizer has no way to
      !! know that the function changed shape rather than the geometry being
      !! bad -- it reads the jump as real, rejects the step, shrinks the trust
      !! radius, and either stalls or converges to a point that is not a
      !! stationary point of anything. Nothing fails; the numbers all look
      !! plausible.
      !!
      !! **The fix.** Generate the list once at the starting geometry and hand
      !! the same one to every step through `supplied_terms`, which
      !! `run_calculation` takes as given and does not re-screen. The energy is
      !! then a fixed sum of a fixed set of subsystems -- a smooth function of
      !! the coordinates -- and its gradient is that function's actual
      !! derivative, which is what makes the optimization well posed at all.
      !!
      !! **What is given up.** A term screened out at the start stays out even
      !! if the molecules later approach. That is the right trade: it is a
      !! definite approximation, fixed and reportable, rather than a surface
      !! that changes underneath the optimizer. A run whose fragments move far
      !! is better served by looser cutoffs, or by `freeze_terms: false` and
      !! the knowledge that convergence is then not guaranteed.
      !!
      !! Only for plain MBE. GMBE derives its PIE terms from overlapping
      !! primaries rather than taking a supplied list, so there is nothing to
      !! freeze through this route and the honest thing is to say so.
      type(system_geometry_t), intent(in) :: sys_geom
      type(driver_config_t), intent(in) :: config

      integer, allocatable :: polymers(:, :)
      integer(int64) :: n_generated

      ctx_terms_frozen = .false.
      ctx_n_terms = 0

      if (config%nlevel <= 0) return

      if (config%allow_overlapping_fragments) then
         if (allocated(config%fragment_cutoffs)) then
            call logger%warning(" ")
            call logger%warning("GMBE term lists cannot be frozen, and screening is active.")
            call logger%warning("  The intersection terms are rederived at every step, so a")
            call logger%warning("  fragment crossing a cutoff changes the energy expression")
            call logger%warning("  mid-optimization. Remove the cutoffs if this stalls.")
            call logger%warning(" ")
         end if
         return
      end if

      if (.not. config%optimization%freeze_terms) then
         if (allocated(config%fragment_cutoffs)) then
            call logger%warning(" ")
            call logger%warning("freeze_terms is off and distance screening is active.")
            call logger%warning("  The term list will be rebuilt at every geometry, which puts")
            call logger%warning("  steps in the energy the optimizer reads as real.")
            call logger%warning(" ")
         end if
         return
      end if

      ! The same call the driver makes, so the frozen list is the list this
      ! run would have used anyway -- not a reimplementation that agrees today.
      call generate_mbe_term_list(sys_geom, config, config%nlevel, polymers, n_generated)

      if (n_generated <= 0) return

      allocate (ctx_terms(n_generated, config%nlevel))
      ctx_terms = polymers(1:n_generated, 1:config%nlevel)
      ctx_n_terms = n_generated
      ctx_terms_frozen = .true.

      call logger%info("  Term list frozen at the starting geometry: "// &
                       to_char(ctx_n_terms)//" terms")
      if (allocated(config%fragment_cutoffs)) then
         call logger%info("    (distance screening applied once, not per step)")
      end if

   end subroutine freeze_term_list

   subroutine report_result(final_energy, converged, tolerance)
      !! Say where it ended up, and whether that is a minimum
      !!
      !! The convergence test is this program's own, on the largest gradient
      !! component of the last step, because DL-FIND hands back geometries and
      !! not a verdict. Without it a run that exhausted `max_steps` printed
      !! "Optimization finished" and a final energy, which is indistinguishable
      !! from a converged one to anything reading the log or the exit code.
      real(dp), intent(in) :: final_energy
      logical, intent(in) :: converged
      real(dp), intent(in) :: tolerance

      call logger%info(" ")
      call logger%info("============================================")
      if (converged) then
         call logger%info("Optimization converged")
      else
         call logger%info("Optimization STOPPED WITHOUT CONVERGING")
      end if
      call logger%info("  Steps accepted: "//to_char(ctx_n_steps))
      call logger%info("  Energy evaluations: "//to_char(ctx_n_evaluations))
      call logger%info("  Final energy: "//to_char(final_energy)//" Hartree")
      call logger%info("  Largest gradient component: "// &
                       to_char(ctx_last_gradient_max)//" Hartree/Bohr")
      call logger%info("  Target: "//to_char(tolerance)//" Hartree/Bohr")
      call logger%info("============================================")
      call logger%info(" ")

   end subroutine report_result

   subroutine write_optimized_xyz(sys_geom, final_energy, converged, error)
      !! Write the optimized geometry beside the input, in Angstrom
      !!
      !! Angstrom and .xyz rather than the JSON document this program usually
      !! writes, because the thing a person does next with an optimized
      !! structure is open it in a viewer or feed it to another program, and
      !! neither reads Bohr out of a JSON field.
      use mqc_io_helpers, only: get_output_json_filename
      type(system_geometry_t), intent(in) :: sys_geom
      real(dp), intent(in) :: final_energy
      logical, intent(in) :: converged
      type(error_t), intent(inout) :: error

      integer :: unit, stat, i
      character(len=:), allocatable :: status_text
      character(len=:), allocatable :: path

      path = optimized_xyz_path()

      open (newunit=unit, file=path, status="replace", action="write", iostat=stat)
      if (stat /= 0) then
         call error%set(ERROR_GENERIC, "Could not write optimized geometry to "//path)
         return
      end if

      ! The comment line says which it is, because an .xyz outlives the log it
      ! came from and a geometry that never converged looks like one that did.
      if (converged) then
         status_text = "converged"
      else
         status_text = "NOT CONVERGED"
      end if

      write (unit, "(i0)") sys_geom%total_atoms
      write (unit, "(a,f20.12,a)") "metalquicha "//status_text//", E = ", final_energy, " Hartree"
      do i = 1, sys_geom%total_atoms
         write (unit, "(a4,3f20.12)") element_number_to_symbol(sys_geom%element_numbers(i)), &
            to_angstrom(sys_geom%coordinates(1, i)), &
            to_angstrom(sys_geom%coordinates(2, i)), &
            to_angstrom(sys_geom%coordinates(3, i))
      end do
      close (unit)

      call logger%info("Optimized geometry written to "//path)

   end subroutine write_optimized_xyz

   function optimized_xyz_path() result(path)
      !! Where the optimized geometry goes: beside the output document
      use mqc_io_helpers, only: get_output_json_filename
      character(len=:), allocatable :: path

      character(len=:), allocatable :: json_name
      integer :: dot

      json_name = get_output_json_filename()
      dot = index(json_name, ".", back=.true.)
      if (dot > 1) then
         path = json_name(1:dot - 1)//"_optimized.xyz"
      else
         path = trim(json_name)//"_optimized.xyz"
      end if

   end function optimized_xyz_path

end module mqc_geometry_optimizer
