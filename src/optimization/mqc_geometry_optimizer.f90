!! Drive a geometry optimization over any method this program can gradient
module mqc_geometry_optimizer
   !! An optimization is a loop *over* calculations, so it sits above
   !! `run_calculation`, and the dispatch is in `main` rather than in
   !! `mqc_driver` to keep the two out of a circular dependency.
   !!
   !! **MPI shape.** Rank 0 runs the optimizer and decides every step; the other
   !! ranks sit in `worker_loop` waiting to be told a geometry. Steps are never
   !! recomputed per rank -- two ranks that disagreed in the last bits of a
   !! coordinate would build different fragments and the job would hang rather
   !! than fail.
   use pic_types, only: dp, int32, int64
   use pic_mpi_lib, only: comm_t, bcast
   use pic_logger, only: logger => global_logger, warning_level, verbose_level
   use pic_io, only: to_char
   use pic_timer, only: timer_type
   use mqc_convergence_report, only: convergence_header
   use mqc_optimizer_types, only: optimizer_settings_t, &
                                  OPT_COORDS_CARTESIAN, OPT_COORDS_UNKNOWN, OPT_ALGO_UNKNOWN, &
                                  coordinates_to_string, algorithm_to_string, &
                                  algorithm_needs_hessian
   use mqc_physical_fragment, only: system_geometry_t, to_angstrom, to_bohr
   use mqc_bond_perception, only: connected_components
   use mqc_config_adapter, only: driver_config_t
   use mqc_config_types, only: bond_t
   use mqc_resources, only: resources_t
   use mqc_result_types, only: calculation_result_t
   use mqc_calc_types, only: CALC_TYPE_GRADIENT, CALC_TYPE_ENERGY, CALC_TYPE_HESSIAN
   use mqc_method_types, only: method_type_to_string, METHOD_TYPE_GFN1, METHOD_TYPE_GFN2
   use mqc_error, only: error_t, ERROR_GENERIC, ERROR_VALIDATION
   use mqc_elements, only: element_number_to_symbol
   use mqc_physical_constants, only: HARTREE_TO_KCALMOL
   use mqc_dlfind_bridge, only: dlfind_available, dlfind_optimize, &
                                dlfind_connected_minima
   use mqc_optimization_output, only: optimization_record_t, write_optimization_json, &
                                      write_trajectory_xyz
   use mqc_frag_utils, only: generate_mbe_term_list
   implicit none

   real(dp), parameter :: ZERO_FREQ_CM = 10.0_dp
      !! Below this a mode is translation or rotation rather than a vibration,
      !! and the same number the analysis prints its own summary against.
   real(dp), parameter :: NOISE_FREQ_CM = 0.1_dp
      !! Under this a negative frequency is arithmetic, not physics: projection
      !! zeroes the translation and rotation rows of the mass-weighted Hessian,
      !! and their eigenvalues come back on either side of zero at the level of
      !! rounding. Only negatives above this are worth mentioning as possibly
      !! being a very flat mode.
   integer, parameter :: MAX_REACTION_MODE_ATOMS = 6
      !! How many atoms of the imaginary mode to name.
   real(dp), parameter :: REACTION_MODE_FLOOR = 0.02_dp
      !! Stop listing once an atom carries under 2% of the motion.
   integer, parameter :: LINE_LEN = 256
   private

   public :: run_geometry_optimization

   ! What rank 0 tells the workers to do next, as one broadcast integer.
   integer(int32), parameter :: OPT_CMD_STOP = 0
   integer(int32), parameter :: OPT_CMD_EVALUATE = 1
   integer(int32), parameter :: OPT_CMD_FINAL = 2
   integer(int32), parameter :: OPT_CMD_HESSIAN = 3

   ! The optimization in progress, reachable from the evaluator. Module state
   ! because the engine calls the evaluator through a plain procedure interface
   ! with no argument to carry a context. One optimization at a time, which is
   ! what `active` guards.
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
   ! The path taken, grown a step at a time and held until the run ends so the
   ! whole thing can be written as one document.
   real(dp), allocatable, save :: ctx_trajectory(:, :, :)
   real(dp), allocatable, save :: ctx_trajectory_energies(:)
   integer, save :: ctx_n_steps = 0
   ! The term list, chosen once and kept. See `freeze_term_list`.
   integer, allocatable, save :: ctx_terms(:, :)
   integer(int64), save :: ctx_n_terms = 0
   logical, save :: ctx_terms_frozen = .false.
   type(timer_type), save :: ctx_wall
      !! Wall clock over the whole engine run, so each step can report how its
      !! time split between the quantum-chemistry evaluation and the optimizer.
   real(dp), save :: ctx_qc_accum = 0.0_dp
      !! Time in `run_step` since the last accepted step was printed,
      !! accumulated across the line-search evaluations that do not print, so
      !! the row's `qc_time` covers all of them and `opt_time` is the
      !! remainder.
   real(dp), save :: ctx_last_print_wall = 0.0_dp
      !! `ctx_wall` reading at the last printed step, the other end of the
      !! interval whose optimizer time is `wall_gap - qc_time`.
   real(dp), save :: ctx_last_gradient_max = huge(1.0_dp)
      !! Largest gradient component of the last successful step, Hartree/Bohr.
      !! This is what says whether the run converged: DL-FIND reports
      !! geometries but not a verdict.

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
      real(dp), allocatable :: endpoint_geom(:, :)
      integer, allocatable :: atomic_numbers(:)
      integer, allocatable :: residues(:)

      rank = resources%mpi_comms%world_comm%rank()
      n_atoms = sys_geom%total_atoms

      if (active) then
         call error%set(ERROR_GENERIC, &
                        "run_geometry_optimization: an optimization is already running")
         return
      end if

      ! Every refusal is decided before anything is broadcast: one discovered
      ! after the workers were told to expect a geometry leaves them waiting
      ! for one that never comes.
      if (rank == 0) then
         if (.not. dlfind_available()) then
            call error%set(ERROR_VALIDATION, &
                           'driver "Optimize" needs the DL-FIND backend; this build has none. '// &
                           "Configure with -DMQC_ENABLE_DLFIND=ON.")
         else if (config%optimization%coordinates == OPT_COORDS_UNKNOWN) then
            ! A backstop for the callers that convert a config without passing
            ! an error to `config_to_driver`.
            call error%set(ERROR_VALIDATION, &
                           "Unknown keywords.optimization.coordinates. Use cartesian, "// &
                           "hdlc, hdlc-tc, dlc or dlc-tc.")
         else if (config%optimization%algorithm == OPT_ALGO_UNKNOWN) then
            call error%set(ERROR_VALIDATION, &
                           "Unknown keywords.optimization.algorithm. Use lbfgs, cg, cg-auto, "// &
                           "sd, prfo, nr or damped.")
         else if (n_atoms < 2) then
            call error%set(ERROR_VALIDATION, &
                           "Nothing to optimize: a single atom has no geometry.")
         else if (config%optimization%hess_end .and. .not. is_restricted_hf(config)) then
            ! Refused before the steps are paid for rather than after.
            call error%set(ERROR_VALIDATION, &
                           "keywords.optimization.hess_end is restricted to restricted "// &
                           "Hartree-Fock for now. Drop the keyword, or run the "// &
                           'optimization and a separate "Hessian" deck on the geometry '// &
                           "it writes.")
         end if
      end if

      ! Every rank has to learn the answer, not just the one that worked it
      ! out: a rank 0 that refused and returned on its own would leave the
      ! workers blocked in `worker_loop` and the job would hang rather than
      ! print the message already sitting in `error`.
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
         ! HDLC's residues decide the coordinate system, and a wrong partition
         ! does not fail -- it converges early, slightly above the minimum. One
         ! residue for a system that visibly comes apart is the sign of that.
         call logger%verbose("  optimizer residues: "//to_char(maxval(residues))// &
                             " for "//to_char(n_atoms)//" atoms")

         ! One gradient before handing over, to find out whether this method
         ! can produce one at all. There is no way to ask in advance -- the
         ! refusal is raised where the gradient would have been built -- and
         ! without the probe it reaches DL-FIND as a failed first step, which
         ! DL-FIND answers with dlf_error and a Fortran backtrace. It costs one
         ! evaluation on a run that works.
         call probe_gradient(n_atoms, coords, probe_status)
         if (probe_status /= 0) then
            call error%set(ERROR_VALIDATION, &
                           "This method cannot produce a gradient, so it cannot be "// &
                           "optimized. The message above says why.")
            call stop_workers(resources%mpi_comms%world_comm)
            call end_context()
            return
         end if

         ! The step table's frame, shared with the SCF and CCSD tables. Rows are
         ! written from the engine's callback as each step is accepted; the
         ! verdict and the totals come from `report_result` below.
         ctx_qc_accum = 0.0_dp
         ctx_last_print_wall = 0.0_dp
         call ctx_wall%start()
         call convergence_header(.true., "optimization steps", &
                                 "    step                 energy        |g|max     qc_time    opt_time", 69)
         ! `endpoint_geom` stays unallocated for an ordinary optimization, and
         ! an unallocated allocatable is absent at an optional dummy, so one
         ! pair of call sites serves both the band and the single structure.
         ! The analytic Hessian below is offered only where there is one; any
         ! other method has DL-FIND build its own from `6N` gradients.
         call load_endpoint(config%optimization%endpoint, atomic_numbers, endpoint_geom, error)
         ! TODO(mqc): this return skips `stop_workers` and `end_context`, so a
         ! bad `endpoint` file hangs every other rank in `worker_loop` and
         ! leaves `active` set for the rest of the process.
         if (error%has_error()) return

         if (algorithm_needs_hessian(config%optimization%algorithm) .and. &
             is_restricted_hf(config)) then
            call dlfind_optimize(config%optimization, n_atoms, atomic_numbers, residues, &
                                 coords, evaluate_energy_gradient, record_step, &
                                 final_energy, error, hessian=evaluate_hessian, &
                                 endpoint=endpoint_geom)
         else
            call dlfind_optimize(config%optimization, n_atoms, atomic_numbers, residues, &
                                 coords, evaluate_energy_gradient, record_step, &
                                 final_energy, error, endpoint=endpoint_geom)
         end if
         call logger%info("  "//repeat("-", 69))

         ! Stop the workers whether or not that succeeded: a rank 0 that
         ! returned an error without sending this would leave N-1 ranks blocked
         ! in bcast for the rest of the job.
         call stop_workers(resources%mpi_comms%world_comm)

         if (.not. error%has_error()) then
            sys_geom%coordinates = coords
            converged = ctx_last_gradient_max <= config%optimization%gradient_tolerance
            call write_final_single_point(sys_geom)
            call report_result(final_energy, converged, config%optimization%gradient_tolerance)
            ! Written either way: a run that ran out of steps still left the
            ! geometry closer to the minimum than it started.
            call write_optimized_xyz(sys_geom, final_energy, converged, error)
            if (.not. error%has_error()) then
               call write_optimization_record(sys_geom, config, final_energy, converged, error)
            end if

            ! Only on a converged geometry: the curvature at a point the
            ! optimizer was still moving away from describes that point and not
            ! the minimum.
            if (config%optimization%connect .and. converged .and. .not. error%has_error()) then
               call report_connected_minima(sys_geom)
            end if
            if (config%optimization%hess_end .and. converged .and. .not. error%has_error()) then
               call run_final_hessian(sys_geom, config%optimization%target)
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
      !! together, which the fragmented path requires and the unfragmented one
      !! tolerates.
      ! TODO(mqc): `run_calculation` is not referenced here, nor in
      ! `worker_loop` or `write_final_single_point`; all three reach it through
      ! `run_step`. Three imports that only make the dependency look wider than
      ! it is.
      use mqc_driver, only: run_calculation
      integer, intent(in) :: n_atoms
      real(dp), intent(in) :: coords(3, n_atoms)
      real(dp), intent(out) :: energy
      real(dp), intent(out) :: gradient(3, n_atoms)
      integer, intent(out) :: status

      type(calculation_result_t) :: result
      type(driver_config_t) :: gradient_config
      integer :: saved_level
      character(len=100) :: line
      type(timer_type) :: step_clock
      real(dp) :: step_time, now_wall, qc_time, opt_time

      energy = 0.0_dp
      gradient = 0.0_dp
      status = 0

      ctx_sys_geom%coordinates = coords
      call send_geometry(ctx_resources%mpi_comms%world_comm, ctx_sys_geom)

      ! The deck said "Optimize"; each step of one is a gradient. Copied rather
      ! than mutated, so the caller's config still says what the user asked for.
      gradient_config = ctx_config
      gradient_config%calc_type = CALC_TYPE_GRADIENT

      ! No files and no single-point report per step: `run_calculation` ends by
      ! logging its energy, dipole, HOMO-LUMO gap and gradient norm, which is
      ! the wrong thing to say a hundred times over. Lowered to warnings so an
      ! unconverged fragment still speaks up, and left alone entirely when the
      ! user asked for verbose.
      saved_level = logger%log_level
      if (saved_level < verbose_level) call logger%configure(level=warning_level)
      call step_clock%start()
      call run_step(gradient_config, result)
      step_time = step_clock%get_elapsed_time()
      call logger%configure(level=saved_level)

      ! Every evaluation is quantum chemistry, printed or not, so a printed step
      ! reports its own cost plus that of the line-search probes since the last.
      ctx_qc_accum = ctx_qc_accum + step_time

      if (.not. ctx_probing) ctx_n_evaluations = ctx_n_evaluations + 1

      ! A fragment that failed leaves the driver going -- it reports every bad
      ! term at the end rather than the first -- so the result carries the
      ! failure and the gradient array carries whatever had accumulated, which
      ! reads as a perfectly plausible small gradient.
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

      ! One row of the step table, the same shape as the SCF and CCSD tables.
      ! Probe evaluations do not print: they are not steps.
      if (.not. ctx_probing) then
         now_wall = ctx_wall%get_elapsed_time()
         qc_time = ctx_qc_accum
         ! Whatever the interval was not spent doing quantum chemistry, DL-FIND
         ! spent on the step. Floored at zero against timer granularity.
         opt_time = max(0.0_dp, (now_wall - ctx_last_print_wall) - qc_time)
         write (line, "(i8,f23.12,es14.3,f10.2,a,f10.2,a)") ctx_n_evaluations, energy, &
            maxval(abs(gradient)), qc_time, " s", opt_time, " s"
         call logger%info(trim(line))
         ctx_last_print_wall = now_wall
         ctx_qc_accum = 0.0_dp
      end if

   end subroutine evaluate_energy_gradient

   subroutine evaluate_hessian(n_atoms, coords, hessian, status)
      !! Second derivatives at a geometry, for an engine that climbs
      !!
      !! The mirror of `evaluate_energy_gradient` for the algorithms that need
      !! curvature rather than slope -- P-RFO and Newton-Raphson. Same MPI
      !! shape: the command goes out, the geometry follows, and every rank runs
      !! the same calculation type, which for a finite-difference Hessian is
      !! what distributes the displacements across them.
      !!
      !! A failure here is reported and not raised: `status` non-zero sends
      !! DL-FIND to its own finite differences.
      integer, intent(in) :: n_atoms
      real(dp), intent(in) :: coords(3, n_atoms)
      real(dp), intent(out) :: hessian(3*n_atoms, 3*n_atoms)
      integer, intent(out) :: status

      type(calculation_result_t) :: result
      type(driver_config_t) :: hess_config
      integer(int32) :: command
      integer :: saved_level

      hessian = 0.0_dp
      status = 1

      command = OPT_CMD_HESSIAN
      call bcast(ctx_resources%mpi_comms%world_comm, command, 1_int32, 0_int32)

      ctx_sys_geom%coordinates = coords
      call send_final_geometry(ctx_resources%mpi_comms%world_comm, ctx_sys_geom)

      hess_config = ctx_config
      hess_config%calc_type = CALC_TYPE_HESSIAN

      saved_level = logger%log_level
      if (saved_level < verbose_level) call logger%configure(level=warning_level)
      call run_step(hess_config, result)
      call logger%configure(level=saved_level)

      if (result%has_error .or. .not. result%has_hessian) return
      if (size(result%hessian, 1) /= 3*n_atoms) return

      hessian = result%hessian
      status = 0

   end subroutine evaluate_hessian

   subroutine worker_loop(comm)
      !! A worker's half of an optimization: compute what rank 0 asks for
      !!
      !! **Keep this in step with `evaluate_energy_gradient` above.** The two
      !! are mirror images written twice rather than shared: a broadcast added
      !! to one and not the other hangs the job instead of failing it.
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
         if (command == OPT_CMD_HESSIAN) step_config%calc_type = CALC_TYPE_HESSIAN

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
      !! chosen and broadcast the command. The order -- one integer, then the
      !! coordinates -- is the order `worker_loop` reads them in.
      type(comm_t), intent(in) :: comm
      type(system_geometry_t), intent(in) :: sys_geom

      real(dp), allocatable :: buffer(:)
      integer(int32) :: n

      ! Flat, because pic-mpi broadcasts rank-1 arrays and the coordinates are
      ! rank-2.
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
      !! optimization produces no `output_<name>.json` at all. It also puts the
      !! dipole and the HOMO-LUMO gap at the optimized geometry rather than the
      !! input one.
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

   pure function is_restricted_hf(config) result(restricted)
      !! Whether this deck is the one case `hess_end` is allowed for
      !!
      !! Restricted Hartree-Fock and nothing else, for now -- a restriction on
      !! what has been checked rather than on what would run.
      use mqc_method_types, only: METHOD_TYPE_HF
      type(driver_config_t), intent(in) :: config
      logical :: restricted

      restricted = config%method_config%method_type == METHOD_TYPE_HF &
                   .and. .not. config%method_config%scf%unrestricted
   end function is_restricted_hf

   subroutine run_final_hessian(sys_geom, want)
      !! A Hessian at the converged geometry, and the verdict it settles
      !!
      !! "Converged" names the condition a minimiser met -- a vanishing
      !! gradient -- and not the thing it found, which a saddle satisfies just
      !! as exactly. The second derivatives are what tell the two apart.
      !!
      !! The frequency table, thermochemistry and warnings are printed by the
      !! Hessian workflow itself; what is added here is the one line saying
      !! which stationary point this is. `write_output` is false, so the
      !! `output_<name>.json` the final energy wrote is left alone.
      use mqc_vibrational_analysis, only: compute_vibrational_analysis
      use mqc_optimizer_types, only: OPT_TARGET_SADDLE
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: want
         !! `OPT_TARGET_MINIMUM` or `OPT_TARGET_SADDLE`. The curvature test is
         !! the same either way; which count is the pass is not.

      type(calculation_result_t) :: result
      type(driver_config_t) :: hess_config
      integer(int32) :: command
      real(dp), allocatable :: frequencies(:), reduced_masses(:), force_constants(:)
      real(dp), allocatable :: cart_disp(:, :)
      integer :: k, n_imag
      real(dp) :: worst, softest

      call logger%info(" ")
      call logger%info("  Hessian at the converged geometry (keywords.optimization.hess_end)")

      command = OPT_CMD_HESSIAN
      call bcast(ctx_resources%mpi_comms%world_comm, command, 1_int32, 0_int32)

      ctx_sys_geom%coordinates = sys_geom%coordinates
      call send_final_geometry(ctx_resources%mpi_comms%world_comm, ctx_sys_geom)

      hess_config = ctx_config
      hess_config%calc_type = CALC_TYPE_HESSIAN

      call run_step(hess_config, result, write_output=.false.)

      if (result%has_error .or. .not. result%has_hessian) then
         ! Not fatal: the optimization converged and its geometry is written.
         ! What failed is the check on it.
         call logger%warning("  the final Hessian could not be computed; the optimized "// &
                             "geometry stands, unverified")
         return
      end if

      call compute_vibrational_analysis(result%hessian, sys_geom%element_numbers, &
                                        frequencies, reduced_masses, force_constants, &
                                        cart_disp, coordinates=sys_geom%coordinates, &
                                        project_trans_rot=.true.)
      if (.not. allocated(frequencies)) then
         call logger%warning("  no frequencies came back from the final Hessian")
         return
      end if

      ! Negative is how an imaginary frequency is carried here -- the square
      ! root of a negative eigenvalue is stored with its sign rather than as a
      ! complex number, which is why this is a comparison and not an `aimag`.
      n_imag = 0
      worst = 0.0_dp
      softest = 0.0_dp
      do k = 1, size(frequencies)
         if (frequencies(k) < -ZERO_FREQ_CM) then
            n_imag = n_imag + 1
            worst = min(worst, frequencies(k))
         else if (frequencies(k) < -NOISE_FREQ_CM) then
            ! Negative, but under the floor that separates a vibration from a
            ! projected translation. Kept because a very flat saddle looks
            ! exactly like this.
            softest = min(softest, frequencies(k))
         end if
      end do

      ! The same three cases, read against what was asked for: a first-order
      ! saddle is the failure when a minimum was wanted and the success when it
      ! was not.
      if (want == OPT_TARGET_SADDLE) then
         if (n_imag == 1) then
            call logger%info("  one imaginary frequency, "//to_char(abs(worst))// &
                             " cm-1: this is a first-order saddle point")
            call report_reaction_mode(frequencies, cart_disp, sys_geom)
         else if (n_imag == 0) then
            call logger%warning("  no imaginary frequencies: this is a minimum, not a "// &
                                "saddle. The search fell off the ridge -- start closer "// &
                                "to the barrier, or use a path method to find a guess.")
            if (softest < 0.0_dp) then
               call logger%warning("     but the lowest mode is "//to_char(abs(softest))// &
                                   " cm-1 imaginary, under the floor that separates a "// &
                                   "vibration from a projected translation. A very flat "// &
                                   "saddle looks exactly like this, so read the frequency "// &
                                   "table before taking the verdict above.")
            end if
         else
            call logger%warning("  "//to_char(n_imag)//" imaginary frequencies, the largest "// &
                                to_char(abs(worst))//" cm-1: this is a stationary point of "// &
                                "order "//to_char(n_imag)//", not a transition state")
         end if
      else
         if (n_imag == 0) then
            call logger%info("  no imaginary frequencies: this is a minimum")
         else if (n_imag == 1) then
            call logger%warning("  one imaginary frequency, "//to_char(abs(worst))// &
                                " cm-1: this is a first-order saddle point, not a minimum")
         else
            call logger%warning("  "//to_char(n_imag)//" imaginary frequencies, the largest "// &
                                to_char(abs(worst))//" cm-1: this is not a minimum")
         end if
      end if
      call logger%info(" ")

   end subroutine run_final_hessian

   subroutine report_connected_minima(sys_geom)
      !! The two minima the saddle falls to, and the barrier from each side
      !!
      !! One imaginary frequency proves a first-order saddle, but not that the
      !! saddle joins the two structures anybody had in mind. Displacing along
      !! the imaginary mode and relaxing in each direction is what settles it.
      !!
      !! Not an intrinsic reaction coordinate: the path downhill is a
      !! minimisation path rather than the steepest-descent one. The endpoints
      !! are the same.
      type(system_geometry_t), intent(in) :: sys_geom

      real(dp), allocatable :: side_a(:, :), side_b(:, :)
      real(dp) :: e_a, e_b, saddle_energy
      logical :: found
      character(len=LINE_LEN) :: line

      call dlfind_connected_minima(side_a, e_a, side_b, e_b, found, saddle_energy)
      if (.not. found) then
         call logger%warning("  the downhill runs did not reach two minima, so what this "// &
                             "saddle connects is unknown")
         return
      end if

      call logger%info("  the saddle relaxes to two minima:")
      write (line, "(a,f20.12,a)") "     forward   ", e_a, " Hartree"
      call logger%info(trim(line))
      write (line, "(a,f20.12,a)") "     backward  ", e_b, " Hartree"
      call logger%info(trim(line))
      write (line, "(a,f14.6,a,f14.6,a)") "     barrier from each side: ", &
         (saddle_energy - e_a)*HARTREE_TO_KCALMOL, " and ", &
         (saddle_energy - e_b)*HARTREE_TO_KCALMOL, " kcal/mol"
      call logger%info(trim(line))
      call write_side_xyz("connected_forward.xyz", sys_geom, side_a, e_a)
      call write_side_xyz("connected_backward.xyz", sys_geom, side_b, e_b)
      call logger%info("  written to connected_forward.xyz and connected_backward.xyz")
      call logger%info(" ")
   end subroutine report_connected_minima

   subroutine write_side_xyz(filename, sys_geom, coords, energy)
      !! One connected minimum, in Angstroms, for whatever opens it next
      character(len=*), intent(in) :: filename
      type(system_geometry_t), intent(in) :: sys_geom
      real(dp), intent(in) :: coords(:, :)
      real(dp), intent(in) :: energy

      integer :: unit, stat, i

      open (newunit=unit, file=filename, status="replace", action="write", iostat=stat)
      if (stat /= 0) return
      write (unit, "(i0)") size(sys_geom%element_numbers)
      write (unit, "(a,f20.12,a)") "metalquicha connected minimum, E = ", energy, " Hartree"
      do i = 1, size(sys_geom%element_numbers)
         write (unit, "(a4,3f20.12)") element_number_to_symbol(sys_geom%element_numbers(i)), &
            to_angstrom(coords(1, i)), to_angstrom(coords(2, i)), to_angstrom(coords(3, i))
      end do
      close (unit)
   end subroutine write_side_xyz

   subroutine load_endpoint(path, atomic_numbers, endpoint, error)
      !! Read the second structure of a chain-of-states run
      !!
      !! Unallocated on return when no endpoint was asked for, which is how the
      !! caller passes "absent" to an optional dummy without branching.
      !!
      !! Atom count and atom *order* are both checked. DL-FIND interpolates
      !! image `i` between coordinate `i` of each structure, so a product
      !! written with two atoms swapped describes a reaction in which those
      !! atoms trade places, and the band converges happily on it.
      use mqc_xyz_reader, only: read_xyz_file
      use mqc_geometry, only: geometry_type
      use mqc_elements, only: element_symbol_to_number
      character(len=:), allocatable, intent(in) :: path
      integer, intent(in) :: atomic_numbers(:)
      real(dp), allocatable, intent(out) :: endpoint(:, :)
      type(error_t), intent(inout) :: error

      type(geometry_type) :: geom
      type(error_t) :: read_error
      integer :: i, z

      if (.not. allocated(path)) return

      call read_xyz_file(path, geom, read_error)
      if (read_error%has_error()) then
         call error%set(ERROR_VALIDATION, "keywords.optimization.endpoint could not be "// &
                        "read: "//read_error%get_message())
         return
      end if

      if (geom%natoms /= size(atomic_numbers)) then
         call error%set(ERROR_VALIDATION, "keywords.optimization.endpoint has "// &
                        trim(to_char(geom%natoms))//" atoms and the starting geometry has "// &
                        trim(to_char(size(atomic_numbers)))// &
                        ". A reaction path runs between two structures of the same system.")
         return
      end if

      do i = 1, geom%natoms
         z = element_symbol_to_number(trim(geom%elements(i)))
         if (z /= atomic_numbers(i)) then
            call error%set(ERROR_VALIDATION, "keywords.optimization.endpoint disagrees with "// &
                           "the starting geometry at atom "//trim(to_char(i))//": "// &
                           trim(element_number_to_symbol(atomic_numbers(i)))//" there, "// &
                           trim(geom%elements(i))//" here. The images are interpolated atom "// &
                           "by atom, so the two files have to list the same atoms in the "// &
                           "same order.")
            return
         end if
      end do

      allocate (endpoint(3, geom%natoms))
      do i = 1, geom%natoms
         endpoint(:, i) = to_bohr(geom%coords(:, i))
      end do
   end subroutine load_endpoint

   subroutine report_reaction_mode(frequencies, cart_disp, sys_geom)
      !! Name the atoms the imaginary mode actually moves
      !!
      !! The same molecule has more than one first-order saddle and P-RFO finds
      !! whichever is nearest the guess; which atoms move along the mode is
      !! what separates them.
      !!
      !! Atoms are ranked by the norm of their own three displacement
      !! components, at most `MAX_REACTION_MODE_ATOMS` of them, stopping below
      !! `REACTION_MODE_FLOOR` of the total motion.
      real(dp), intent(in) :: frequencies(:)
      real(dp), intent(in) :: cart_disp(:, :)
         !! `(3N, mode)`, Cartesian, one column per mode
      type(system_geometry_t), intent(in) :: sys_geom

      real(dp), allocatable :: amplitude(:)
      integer, allocatable :: order(:)
      integer :: n_atoms, k, imode, i, j, swap
      real(dp) :: total
      character(len=LINE_LEN) :: line

      imode = 0
      do k = 1, size(frequencies)
         if (frequencies(k) < -ZERO_FREQ_CM) then
            imode = k
            exit
         end if
      end do
      if (imode == 0) return
      if (size(cart_disp, 2) < imode) return

      n_atoms = size(sys_geom%element_numbers)
      if (size(cart_disp, 1) < 3*n_atoms) return

      allocate (amplitude(n_atoms), order(n_atoms))
      do i = 1, n_atoms
         amplitude(i) = norm2(cart_disp(3*i - 2:3*i, imode))
         order(i) = i
      end do
      total = sum(amplitude)
      if (total <= 0.0_dp) return

      ! A selection sort: this runs once per optimization.
      do i = 1, n_atoms - 1
         do j = i + 1, n_atoms
            if (amplitude(order(j)) > amplitude(order(i))) then
               swap = order(i)
               order(i) = order(j)
               order(j) = swap
            end if
         end do
      end do

      ! Atoms are numbered from one, as everywhere else in this log, and from
      ! zero in a deck.
      call logger%info("  the imaginary mode moves, most to least (atoms numbered from 1):")
      do k = 1, min(n_atoms, MAX_REACTION_MODE_ATOMS)
         i = order(k)
         if (amplitude(i)/total < REACTION_MODE_FLOOR) exit
         write (line, "(a,i0,a,a,a,f5.1,a)") "     atom ", i, " ", &
            trim(element_number_to_symbol(sys_geom%element_numbers(i))), &
            "  ", 100.0_dp*amplitude(i)/total, "% of the motion"
         call logger%info(trim(line))
      end do
   end subroutine report_reaction_mode

   subroutine record_step(n_atoms, coords, energy)
      !! Keep one accepted geometry, growing the trajectory to fit
      !!
      !! The store doubles when it fills, rather than growing by one step.
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
      ! TODO(mqc): `get_output_json_filename` is not referenced here, nor in
      ! `write_optimized_xyz`; both go through `optimization_path` and
      ! `optimized_xyz_path`, which import it themselves.
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

      ! Trimmed to the steps actually taken: the array is grown by doubling, so
      ! its last frames are whatever the previous allocation held.
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
      !! The probe is not a step, so it is named rather than numbered.
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
      !! configuration.
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
      !! A residue is a group of atoms an optimizer may describe with bonds,
      !! angles and torsions among themselves, numbered from 1 in `residues`.
      !!
      !! A fragmented calculation uses its monomers. An unfragmented one
      !! perceives its own connectivity, so a cluster comes back as one residue
      !! per molecule and HDLC does not try to build coordinates across the gap
      !! between them; a single molecule comes back as one residue.
      !!
      !! This changes no energy, only which internal coordinates the optimizer
      !! works in.
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(out) :: residues(:)

      integer, allocatable :: component(:)
      integer :: imon, k, atom, n_components

      residues = 1
      if (sys_geom%n_monomers <= 1) then
         if (size(residues) < 2) return
         call connected_components(sys_geom, component, n_components)
         ! One component is an ordinary molecule and already has its answer.
         if (n_components > 1) residues = component(1:size(residues))
         return
      end if

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
      ! The method and basis up front: an optimization suppresses the per-step
      ! single-point banner, so otherwise the theory appears only in the final
      ! single point, too late to notice a deck that asked for the wrong one.
      call logger%info("  Method: "//trim(method_type_to_string(config%method_config%method_type)))
      if (config%method_config%method_type /= METHOD_TYPE_GFN1 .and. &
          config%method_config%method_type /= METHOD_TYPE_GFN2) then
         call logger%info("  Basis set: "//trim(config%method_config%basis_set))
      end if
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
      !! Screening a fragmented run's term list at every geometry lets a dimer
      !! cross the cutoff mid-optimization, and the energy then jumps by that
      !! dimer's whole interaction. The optimizer reads the jump as real:
      !! nothing fails, and it stalls or converges to a point that is not a
      !! stationary point of anything.
      !!
      !! Generating the list once and passing it to every step through
      !! `supplied_terms`, which `run_calculation` does not re-screen, makes
      !! the energy a smooth function of the coordinates again. The cost is
      !! that a term screened out at the start stays out even if the molecules
      !! later approach.
      !!
      !! Plain MBE only. GMBE derives its PIE terms from overlapping primaries
      !! rather than from a supplied list, so there is nothing to freeze; a
      !! warning is logged instead.
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

      ! The same call the driver makes, so the frozen list is the list this run
      ! would have used anyway.
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
      !! not a verdict.
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
      !! Angstrom and .xyz, whatever the run converged to; the comment line
      !! says which of the two it was.
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

      ! An .xyz outlives the log it came from, and a geometry that never
      ! converged looks like one that did.
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
