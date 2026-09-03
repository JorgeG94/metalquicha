!! Conformer sampling, driven through CREST with mqc supplying the gradients
module mqc_crest_driver
   !! CREST keeps its own algorithms, defaults and bookkeeping; what changes is
   !! only where the energies come from.
   !!
   !! **Two calculation levels.** Level one is xTB in process, sampling; level
   !! two is mqc, applied to the ensemble the search leaves behind. An iMTD-GC
   !! run asks for tens of thousands of gradients and almost none of them need
   !! to be good.
   !!
   !! **CREST is configured by its own parser, not field by field here.** An
   !! argument vector is built in memory -- no command line, no configuration
   !! file -- and handed to `parseflags`, exactly as `crest_main` does with the
   !! real argv.
   use, intrinsic :: iso_fortran_env, only: output_unit
   use omp_lib, only: omp_get_max_threads, omp_set_num_threads
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_mpi_lib, only: abort_comm
   use mqc_error, only: error_t, ERROR_VALIDATION, ERROR_IO
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_config_adapter, only: driver_config_t
   use mqc_resources, only: resources_t
   use mqc_calculation_interface, only: compute_energy_and_forces
   use mqc_calc_types, only: CALC_TYPE_ENERGY
   use mqc_physical_constants, only: BOHR_TO_ANGSTROM
   use mqc_elements, only: element_number_to_symbol
   use crest_data, only: systemdata, timer, refine
   use crest_calculator, only: jobtype, calculation_settings
   implicit none
   private

   public :: run_conformer_search
   public :: crest_available

   interface
      subroutine parseflags(env, arg, nra)
         import :: systemdata
         implicit none
         type(systemdata), intent(inout) :: env
         integer, intent(in) :: nra
         ! allow(missing-intent)
         character(len=*) :: arg(nra)
         !! No intent, matching the definition in confparse.f90; intent is part
         !! of a procedure's characteristics, so declaring one here would make
         !! this interface disagree with the procedure it describes.
      end subroutine parseflags

      subroutine crest_search_imtdgc(env, tim)
         import :: systemdata, timer
         implicit none
         type(systemdata), intent(inout) :: env
         type(timer), intent(inout) :: tim
      end subroutine crest_search_imtdgc

      subroutine new_ompautoset(env, modus, maxjobs, parallel_jobs, cores_per_job)
         import :: systemdata
         implicit none
         type(systemdata), intent(inout) :: env
         character(len=*), intent(in) :: modus
         integer, intent(in) :: maxjobs
         integer, intent(out) :: parallel_jobs
         integer, intent(out) :: cores_per_job
      end subroutine new_ompautoset
   end interface

   ! What the refinement level needs, waiting for a callback that is handed a
   ! geometry and nothing else. Module state because that callback takes no
   ! context of its own, so only one search can be in flight -- true anyway,
   ! since CREST runs on a single rank.
   !
   ! One context, and CREST reaches the callback from concurrent OpenMP tasks,
   ! so `crest_engrad` serialises itself rather than duplicating this per
   ! thread. See the comment there for why duplicating it would not be enough.
   type(system_geometry_t), save :: ctx_geom
   type(driver_config_t), save :: ctx_config
   type(resources_t), save :: ctx_resources
   logical, save :: ctx_ready = .false.
   integer, save :: ctx_calls = 0

contains

   pure function crest_available() result(available)
      !! Whether this build can sample conformers
      logical :: available
      available = .true.
   end function crest_available

   subroutine run_conformer_search(resources, config, sys_geom, error)
      !! Sample conformers, refining the survivors with this program
      type(resources_t), intent(in) :: resources
      type(driver_config_t), intent(in) :: config
      type(system_geometry_t), intent(in) :: sys_geom
      type(error_t), intent(inout) :: error

      type(systemdata) :: env
      type(timer) :: tim
      type(calculation_settings) :: refine_level
      character(len=:), allocatable :: argv(:)
      character(len=:), allocatable :: threads_text
      character(len=*), parameter :: START_FILE = "crest_input.xyz"
      character(len=*), parameter :: THREAD_FLAG = "-T"
      integer :: unit, i, io
      integer :: n_threads, parallel_jobs, cores_per_job

      ! The thread count the launcher asked for, read before anything below
      ! changes it. CREST defaults to one thread and consults OMP_NUM_THREADS
      ! only when `-T` was absent *and* nothing set the count manually, so an
      ! argument vector without it samples serially however many threads the
      ! job was given.
      n_threads = max(1, omp_get_max_threads())

      ! **Single rank, refused rather than deadlocked.** CREST samples on one
      ! rank with OpenMP threads underneath, while `compute_energy_and_forces`
      ! expects every rank to call it. The two cannot both be satisfied, and
      ! the failure if this is not caught is a hang rather than an error.
      !
      ! Checked here rather than at build time so one binary serves both and
      ! `mpirun -n 1` still works.
      if (resources%mpi_comms%world_comm%size() > 1) then
         call error%set(ERROR_VALIDATION, "conformer sampling runs on a single rank. "// &
                        "CREST samples on one rank with OpenMP underneath, and the "// &
                        "gradient interface it calls expects every rank to take part, "// &
                        "so the two cannot both be satisfied. Re-run without mpirun, "// &
                        "or with -n 1.")
         return
      end if

      if (sys_geom%total_atoms < 1) then
         call error%set(ERROR_VALIDATION, "conformer sampling needs a molecule")
         return
      end if

      ! The starting structure, written once. This is not the file-based
      ! gradient exchange the callback exists to avoid: CREST reads it at
      ! setup, the way any program reads its input geometry, and every
      ! gradient after this crosses in memory.
      open (newunit=unit, file=START_FILE, status="replace", action="write", iostat=io)
      if (io /= 0) then
         call error%set(ERROR_IO, "could not write "//START_FILE)
         return
      end if
      write (unit, "(i0)") sys_geom%total_atoms
      write (unit, "(a)") "written by mqc for conformer sampling"
      ! Symbols, not atomic numbers. CREST's reader takes an element symbol
      ! here, and handed a bare integer it parses nothing and leaves `mol`
      ! unusable -- which surfaces as a segfault later, inside `axis`, rather
      ! than as a read error.
      do i = 1, sys_geom%total_atoms
         write (unit, "(a2,3f20.12)") element_number_to_symbol(sys_geom%element_numbers(i)), &
            sys_geom%coordinates(:, i)*BOHR_TO_ANGSTROM
      end do
      close (unit)

      call install_context(sys_geom, config, resources)

      call tim%init(20)
      ! `-T n`, which is how CREST is told a thread count on a real command
      ! line; going through its parser rather than assigning `env%threads` also
      ! sets the two flags that stop it overriding the value later.
      threads_text = int_text(n_threads)
      allocate (character(len=max(len(START_FILE), len(threads_text))) :: argv(3))
      argv(1) = START_FILE
      argv(2) = THREAD_FLAG
      argv(3) = threads_text
      call parseflags(env, argv, 3)

      ! The OpenMP setup `crest_main` does after its own parse, which nothing
      ! else performs: every sampling stage divides `env%threads`, but the
      ! process thread count and OMP_NUM_THREADS are established here.
      !
      ! `env%omp_allow_nested` is deliberately left at CREST's default, which
      ! is on, and it is worth saying why turning it off is not the safety
      ! measure the name suggests. CREST's parallel loops gate the *per-thread*
      ! `ompmklset(Tn)` on that flag -- the call that gives each worker the
      ! thread count it is meant to compute with, one apiece whenever there are
      ! at least as many structures as threads. Switching it off does not make
      ! the inner work serial; it leaves every worker at the outer count and
      ! takes away the only thing that was pinning the engines down. The
      ! parallelism that matters here is over structures either way.
      call new_ompautoset(env, "max", 0, parallel_jobs, cores_per_job)

      if (env%calc%ncalculations < 1) then
         call error%set(ERROR_VALIDATION, "CREST set up no calculation level to sample with")
         return
      end if

      ! Level two: this program, on what survives. Same arrangement CREST uses
      ! for two xTB levels in `legacy_wrappers.f90`.
      refine_level%id = jobtype%callback
      refine_level%external_engrad => crest_engrad
      refine_level%refine_lvl = refine%singlepoint
      call env%calc%add(refine_level)
      if (allocated(env%refine_queue)) deallocate (env%refine_queue)
      call env%addrefine(refine%singlepoint)

      call logger%info("")
      call logger%info("  conformer sampling: CREST samples, this program refines")
      call logger%info("  sampling threads: "//int_text(parallel_jobs)// &
                       ", "//int_text(cores_per_job)//" per calculation")

      call crest_search_imtdgc(env, tim)

      if (env%iostatus_meta /= 0) then
         call error%set(ERROR_VALIDATION, "the conformer search did not complete; "// &
                        "CREST reported status "//int_text(env%iostatus_meta)// &
                        ". Its own output above says more than this can.")
         return
      end if

      ! A search that never reached the refinement level looks like a success
      ! from every other angle, and the ensemble energies would silently be
      ! the sampling method's rather than this program's.
      if (ctx_calls < 1) then
         call error%set(ERROR_VALIDATION, "the search finished without ever refining. "// &
                        "The ensemble energies belong to the sampling method, not to "// &
                        "the requested one.")
         return
      end if

      call logger%info("  refinement calculations: "//int_text(ctx_calls))
      flush (output_unit)
   end subroutine run_conformer_search

   subroutine install_context(sys_geom, config, resources)
      type(system_geometry_t), intent(in) :: sys_geom
      type(driver_config_t), intent(in) :: config
      type(resources_t), intent(in) :: resources

      ctx_geom = sys_geom
      ctx_config = config
      ! **An energy, not a gradient.** The level this installs is registered
      ! as `refine%singlepoint`, and CREST uses it to re-rank an ensemble: it
      ! reads the energy and never looks at the gradient. A level registered as
      ! `refine%geoopt` does need gradients, and would make this conditional on
      ! the stage.
      ctx_config%calc_type = CALC_TYPE_ENERGY
      ctx_resources = resources
      ctx_ready = .true.
      ctx_calls = 0
   end subroutine install_context

   subroutine crest_engrad(nat, at, xyz, energy, grad, iostatus)
      !! CREST's `engrad_callback`, answered by this program
      !!
      !! Coordinates arrive in bohr and the gradient is returned in Eh/bohr,
      !! which is what every in-process calculator on CREST's side hands back.
      !! `system_geometry_t` also holds bohr, so nothing is converted.
      integer, intent(in) :: nat
      integer, intent(in) :: at(nat)
      real(dp), intent(in) :: xyz(3, nat)
      real(dp), intent(out) :: energy
      real(dp), intent(out) :: grad(3, nat)
      integer, intent(out) :: iostatus

      integer :: saved_threads

      iostatus = 0
      energy = 0.0_dp
      grad = 0.0_dp

      if (.not. ctx_ready) then
         iostatus = 1
         return
      end if
      ! A sampling run may not change what the molecule is. If it ever does,
      ! the installed context is the wrong one rather than a resizable one.
      if (nat /= ctx_geom%total_atoms) then
         iostatus = 2
         return
      end if
      if (any(at /= ctx_geom%element_numbers)) then
         iostatus = 3
         return
      end if

      ! **One refinement at a time.** CREST re-ranks an ensemble from concurrent
      ! OpenMP tasks in `crest_sploop`, and every one of them arrives here.
      ! What this calls is `run_calculation`, the whole program's driver: it
      ! writes the global logger, clamps the process thread count when the
      ! method is xTB, and keeps module state of its own. None of that is
      ! re-entrant, and the context above is single besides, so two calls in
      ! flight would interleave two geometries into one calculation. Giving
      ! each thread its own context would fix the second of those and none of
      ! the first, so the lock is the honest answer.
      !
      ! The threads are not wasted by this: the tens of thousands of gradients
      ! a search spends its time on are CREST's own sampling level, which runs
      ! them concurrently and never reaches this pointer. Only the few hundred
      ! survivors are refined here.
      !$omp critical(mqc_crest_engrad)
      ctx_calls = ctx_calls + 1

      ! `run_calculation` sets the thread count to one for xTB and, by its own
      ! comment, does not put it back. Left that way the calling thread would
      ! carry the clamp into whatever CREST does next.
      saved_threads = omp_get_max_threads()

      ctx_geom%coordinates = xyz
      ! `write_output` off: the driver writes results.json on every call
      ! otherwise, so a search asking for a few thousand energies produced a
      ! few thousand writes of a file each one immediately overwrote.
      !
      ! No `gradient` argument, so none is computed. `grad` stays zero, which
      ! is what CREST does with it at this refinement stage anyway.
      call compute_energy_and_forces(ctx_geom, ctx_config, ctx_resources, &
                                     energy, write_output=.false.)

      call omp_set_num_threads(saved_threads)
      !$omp end critical(mqc_crest_engrad)
   end subroutine crest_engrad

   pure function int_text(value) result(text)
      integer, intent(in) :: value
      character(len=:), allocatable :: text
      character(len=12) :: buffer

      write (buffer, "(i0)") value
      text = trim(buffer)
   end function int_text

end module mqc_crest_driver
