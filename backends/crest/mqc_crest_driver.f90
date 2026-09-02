!! Conformer sampling, driven through CREST with mqc supplying the gradients
module mqc_crest_driver
   !! CREST keeps its own algorithms, defaults and bookkeeping; what changes is
   !! only where the energies come from.
   !!
   !! **Two calculation levels, and the split is the point.** An iMTD-GC run
   !! asks for tens of thousands of gradients and almost none of them need to
   !! be good, so level one is xTB in process and level two is mqc, applied to
   !! the ensemble the search leaves behind. On water with an HF/STO-3G
   !! refinement that measured 1432 mqc calls against roughly 56000 gradients.
   !! Sending every gradient to mqc would work and nobody would wait for it.
   !!
   !! **CREST is configured by its own parser, not field by field here.**
   !! `systemdata` carries on the order of a hundred fields and the algorithms
   !! read more of them than their entry points name; filling it by hand would
   !! be a slow way of reimplementing `confparse` and getting the defaults
   !! subtly wrong. An argument vector is built in memory instead -- no command
   !! line, no configuration file -- and handed to `parseflags`, exactly as
   !! `crest_main` does with the real argv.
   use, intrinsic :: iso_fortran_env, only: output_unit
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
         !! No intent, matching the definition in confparse.f90; intent is part
         !! of a procedure's characteristics, so declaring one here would make
         !! this interface disagree with the procedure it describes.
         ! allow(missing-intent)
         character(len=*) :: arg(nra)
      end subroutine parseflags

      subroutine crest_search_imtdgc(env, tim)
         import :: systemdata, timer
         implicit none
         type(systemdata), intent(inout) :: env
         type(timer), intent(inout) :: tim
      end subroutine crest_search_imtdgc
   end interface

   !! What the refinement level needs, waiting for a callback that is handed a
   !! geometry and nothing else.
   !!
   !! Module state for the same reason `mqc_geometry_optimizer` keeps its
   !! `ctx_config` and `ctx_resources` that way. It also means only one search
   !! can be in flight, which is true anyway: CREST runs on a single rank.
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
      character(len=*), parameter :: START_FILE = "crest_input.xyz"
      integer :: unit, i, io

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
      allocate (character(len=len(START_FILE)) :: argv(1))
      argv(1) = START_FILE
      call parseflags(env, argv, 1)

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
      ! as `refine%singlepoint`, and CREST uses it to re-rank an ensemble --
      ! it reads the energy and never looks at the gradient. Computing one
      ! anyway is invisible at HF/STO-3G, where the arithmetic rounds to zero
      ! against the per-call overhead, and is most of the cost at anything
      ! worth refining with: an RI-MP2 gradient is several times its energy.
      !
      ! If a level is ever registered as `refine%geoopt`, which does need
      ! gradients, this has to become conditional on the stage.
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

      !$omp atomic
      ctx_calls = ctx_calls + 1

      ctx_geom%coordinates = xyz
      ! `write_output` off: the driver writes results.json on every call
      ! otherwise, so a search asking for a few thousand energies produced a
      ! few thousand writes of a file each one immediately overwrote.
      !
      ! No `gradient` argument, so none is computed. `grad` stays zero, which
      ! is what CREST does with it at this refinement stage anyway.
      call compute_energy_and_forces(ctx_geom, ctx_config, ctx_resources, &
                                     energy, write_output=.false.)
   end subroutine crest_engrad

   pure function int_text(value) result(text)
      integer, intent(in) :: value
      character(len=:), allocatable :: text
      character(len=12) :: buffer

      write (buffer, "(i0)") value
      text = trim(buffer)
   end function int_text

end module mqc_crest_driver
