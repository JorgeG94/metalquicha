!! DL-FIND, reached through libdlfind's C entry point
module mqc_dlfind_bridge
   !! Turns an `optimizer_settings_t` and an evaluator into a call to
   !! `api_dl_find`, and turns DL-FIND's seven callbacks back into calls to
   !! that evaluator.
   !!
   !! **Why the state is at module scope.** DL-FIND's callbacks are plain
   !! function pointers with no user-data argument -- there is nowhere to hang
   !! a context, so the system being optimized has to be somewhere the
   !! callbacks can see without being handed it. One optimization at a time
   !! follows from that, and `mqc_geometry_optimizer` is what enforces it.
   !!
   !! **The interface is C, deliberately.** `api_dl_find` is bind(c), so
   !! nothing here needs libdlfind's .mod files and the two projects may be
   !! built by different compilers. The interface block below is this
   !! program's own declaration of it; there is no `use` of anything from
   !! libdlfind anywhere in this file.
   !!
   !! **A DL-FIND failure ends the process.** libdlfind's `dlf_error` calls
   !! back here and then executes `error stop` regardless of what this does, so
   !! `cb_error` gets one chance to say what was being optimized before the job
   !! dies. There is no way to turn that into a returned error without patching
   !! libdlfind.
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_funptr, c_funloc
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_optimizer_types, only: optimizer_settings_t, energy_gradient_i, step_callback_i, &
                                  OPT_COORDS_CARTESIAN, OPT_COORDS_HDLC, OPT_COORDS_DLC, &
                                  OPT_COORDS_HDLC_TC, OPT_COORDS_DLC_TC, &
                                  OPT_ALGO_SD, OPT_ALGO_CG, OPT_ALGO_LBFGS, OPT_ALGO_PRFO
   use mqc_error, only: error_t, ERROR_GENERIC, ERROR_VALIDATION
   implicit none
   private

   public :: dlfind_available
   public :: dlfind_optimize

   ! DL-FIND's own numbering, which this file is the only place to know.
   ! `icoord`: mod 10 selects the coordinate system, and the tens and hundreds
   ! digits select multi-image methods this bridge does not offer.
   integer(c_int), parameter :: DLF_COORDS_CARTESIAN = 0
   integer(c_int), parameter :: DLF_COORDS_HDLC = 1
   integer(c_int), parameter :: DLF_COORDS_DLC = 3
   integer(c_int), parameter :: DLF_COORDS_HDLC_TC = 2
   integer(c_int), parameter :: DLF_COORDS_DLC_TC = 4

   integer(c_int), parameter :: DLF_OPT_SD = 0
   integer(c_int), parameter :: DLF_OPT_CG = 2
   integer(c_int), parameter :: DLF_OPT_LBFGS = 3
   integer(c_int), parameter :: DLF_OPT_PRFO = 10

   ! The optimization in progress. See the note above on why this is here.
   procedure(energy_gradient_i), pointer, save :: evaluator => null()
   procedure(step_callback_i), pointer, save :: on_step => null()
   type(optimizer_settings_t), save :: settings
   integer, save :: n_atoms = 0
   integer, allocatable, save :: atomic_numbers(:)
   integer, allocatable, save :: atom_residues(:)
   real(dp), allocatable, save :: initial_coords(:, :)
   real(dp), allocatable, save :: latest_coords(:, :)
   real(dp), save :: latest_energy = 0.0_dp
   logical, save :: have_result = .false.
   logical, save :: evaluation_failed = .false.

   interface
      subroutine api_dl_find(nvarin, nvarin2, nspec, master, &
                             c_dlf_error, c_dlf_get_gradient, c_dlf_get_hessian, &
                             c_dlf_get_multistate_gradients, c_dlf_get_params, &
                             c_dlf_put_coords, c_dlf_update) bind(C, name="api_dl_find")
         !! libdlfind's entry point, declared here rather than used from there
         import :: c_int, c_funptr
         implicit none
         integer(c_int), value :: nvarin   !! 3*nat
         integer(c_int), value :: nvarin2  !! Size of the second coordinate array
         integer(c_int), value :: nspec    !! Size of the integer specification array
         integer(c_int), value :: master   !! 1 if this task leads a parallel run
         type(c_funptr), value :: c_dlf_error
         type(c_funptr), value :: c_dlf_get_gradient
         type(c_funptr), value :: c_dlf_get_hessian
         type(c_funptr), value :: c_dlf_get_multistate_gradients
         type(c_funptr), value :: c_dlf_get_params
         type(c_funptr), value :: c_dlf_put_coords
         type(c_funptr), value :: c_dlf_update
      end subroutine api_dl_find
   end interface

contains

   pure function dlfind_available() result(available)
      !! Whether this build can optimize a geometry
      logical :: available
      available = .true.
   end function dlfind_available

   subroutine dlfind_optimize(opt_settings, natoms, znuc, residues, coords, &
                              energy_gradient, step_taken, final_energy, error)
      !! Minimize `coords` in place, calling `energy_gradient` for each geometry
      !!
      !! Coordinates are Bohr in and Bohr out. On success `coords` holds the
      !! optimized geometry and `final_energy` its energy; on failure both are
      !! left as they were and `error` says why.
      type(optimizer_settings_t), intent(in) :: opt_settings
      integer, intent(in) :: natoms
      integer, intent(in) :: znuc(natoms)  !! Atomic numbers, for internal coordinates
      integer, intent(in) :: residues(natoms)
         !! Which residue each atom belongs to, 1-based. HDLC builds internal
         !! coordinates within a residue and keeps Cartesians between them, so
         !! this is the setting that decides whether HDLC helps at all.
      real(dp), intent(inout) :: coords(3, natoms)
      procedure(energy_gradient_i) :: energy_gradient
      procedure(step_callback_i) :: step_taken
         !! Called once per accepted geometry, for the trajectory
      real(dp), intent(out) :: final_energy
      type(error_t), intent(inout) :: error

      integer(c_int) :: nvar, nspec

      final_energy = 0.0_dp

      if (natoms < 1) then
         call error%set(ERROR_VALIDATION, "dlfind_optimize: nothing to optimize")
         return
      end if

      ! P-RFO climbs to a saddle rather than descending to a minimum, and it
      ! needs a Hessian this bridge does not supply. Refused here, naming the
      ! setting, rather than left to converge on something nobody asked for.
      if (opt_settings%algorithm == OPT_ALGO_PRFO) then
         call error%set(ERROR_VALIDATION, &
                        'keywords.optimization.algorithm "prfo" searches for a transition '// &
                        "state and is not wired up yet. Use lbfgs, cg or sd.")
         return
      end if

      settings = opt_settings
      n_atoms = natoms
      evaluator => energy_gradient
      on_step => step_taken
      have_result = .false.
      evaluation_failed = .false.

      if (allocated(atomic_numbers)) deallocate (atomic_numbers)
      if (allocated(atom_residues)) deallocate (atom_residues)
      if (allocated(initial_coords)) deallocate (initial_coords)
      if (allocated(latest_coords)) deallocate (latest_coords)
      allocate (atomic_numbers(natoms), source=znuc)
      allocate (atom_residues(natoms), source=residues)
      allocate (initial_coords(3, natoms), source=coords)
      allocate (latest_coords(3, natoms), source=coords)

      nvar = int(3*natoms, c_int)

      ! nspec is fixed by DL-FIND as nat + nz + 5*ncons + 2*nconn + nat, and it
      ! checks: a mismatch is refused as an interface error rather than read
      ! past. With nz = nat and no constraints that is 3*nat, laid out as
      ! [frozen(nat), znuc(nat), micspec(nat)]. Passing znuc rather than
      ! leaving nz at zero is what lets HDLC and DLC perceive connectivity --
      ! without it those coordinate systems have no elements to work from.
      nspec = int(3*natoms, c_int)

      ! nvarin2 = 0: the second coordinate array carries NEB images, a reaction
      ! path or per-atom masses and weights, none of which a minimization uses.
      ! DL-FIND still hands the callback a length of max(nvarin2,1).
      call api_dl_find(nvar, 0_c_int, nspec, 1_c_int, &
                       c_funloc(cb_error), &
                       c_funloc(cb_get_gradient), &
                       c_funloc(cb_get_hessian), &
                       c_funloc(cb_get_multistate_gradients), &
                       c_funloc(cb_get_params), &
                       c_funloc(cb_put_coords), &
                       c_funloc(cb_update))

      if (evaluation_failed) then
         call error%set(ERROR_GENERIC, &
                        "Geometry optimization stopped: an energy and gradient step failed. "// &
                        "The message above says which.")
      else if (.not. have_result) then
         ! DL-FIND reports geometries through dlf_put_coords, and only when
         ! printf >= 1 -- which cb_get_params sets. Reaching here means it
         ! returned without ever offering one, so there is nothing to hand
         ! back and overwriting the caller's geometry with the starting one
         ! dressed as a result would be worse than saying so.
         call error%set(ERROR_GENERIC, &
                        "DL-FIND returned no geometry. Nothing was optimized.")
      else
         coords = latest_coords
         final_energy = latest_energy
      end if

      evaluator => null()
      on_step => null()
      if (allocated(atom_residues)) deallocate (atom_residues)
      if (allocated(atomic_numbers)) deallocate (atomic_numbers)
      if (allocated(initial_coords)) deallocate (initial_coords)
      if (allocated(latest_coords)) deallocate (latest_coords)
      n_atoms = 0

   end subroutine dlfind_optimize

   subroutine cb_get_params(nvar, nvar2, nspec, coords, coords2, spec, ierr, tolerance, &
                            printl, maxcycle, maxene, tatoms, icoord, iopt, iline, &
                            maxstep, scalestep, lbfgs_mem, nimage, nebk, dump, restart, &
                            nz, ncons, nconn, update, maxupd, delta, soft, inithessian, &
                            carthessian, tsrel, maxrot, tolrot, nframe, nmass, nweight, &
                            timestep, fric0, fricfac, fricp, imultistate, state_i, &
                            state_j, pf_c1, pf_c2, gp_c3, gp_c4, ln_t1, ln_t2, printf, &
                            tolerance_e, distort, massweight, minstep, maxdump, task, &
                            temperature, po_pop_size, po_radius, po_contraction, &
                            po_tolerance_r, po_tolerance_g, po_distribution, po_maxcycle, &
                            po_init_pop_size, po_reset, po_mutation_rate, po_death_rate, &
                            po_scalefac, po_nsave, ntasks, tdlf_farm, n_po_scaling, &
                            neb_climb_test, neb_freeze_test, nzero, coupled_states, &
                            qtsflag, imicroiter, maxmicrocycle, micro_esp_fit) bind(c)
      !! Hand DL-FIND the starting geometry and the settings we care about
      !!
      !! Every argument arrives already holding DL-FIND's own default --
      !! `dlf_default_init` runs before this is called -- so this overwrites
      !! only what the deck actually asked for and leaves the other eighty
      !! alone. That is why there is no second table of defaults here to drift
      !! out of step with DL-FIND's.
      integer(c_int), intent(in), value :: nvar, nvar2, nspec
      real(c_double), intent(inout) :: coords(nvar), coords2(nvar2)
      integer(c_int), intent(inout) :: spec(nspec)
      integer(c_int), intent(out) :: ierr
      real(c_double), intent(inout) :: tolerance
      integer(c_int), intent(inout) :: printl, maxcycle, maxene, tatoms, icoord, iopt, iline
      real(c_double), intent(inout) :: maxstep, scalestep
      integer(c_int), intent(inout) :: lbfgs_mem, nimage
      real(c_double), intent(inout) :: nebk
      integer(c_int), intent(inout) :: dump, restart, nz, ncons, nconn, update, maxupd
      real(c_double), intent(inout) :: delta, soft
      integer(c_int), intent(inout) :: inithessian, carthessian, tsrel, maxrot
      real(c_double), intent(inout) :: tolrot
      integer(c_int), intent(inout) :: nframe, nmass, nweight
      real(c_double), intent(inout) :: timestep, fric0, fricfac, fricp
      integer(c_int), intent(inout) :: imultistate, state_i, state_j
      real(c_double), intent(inout) :: pf_c1, pf_c2, gp_c3, gp_c4, ln_t1, ln_t2
      integer(c_int), intent(inout) :: printf
      real(c_double), intent(inout) :: tolerance_e, distort
      integer(c_int), intent(inout) :: massweight
      real(c_double), intent(inout) :: minstep
      integer(c_int), intent(inout) :: maxdump, task
      real(c_double), intent(inout) :: temperature
      integer(c_int), intent(inout) :: po_pop_size
      real(c_double), intent(inout) :: po_radius, po_contraction, po_tolerance_r, po_tolerance_g
      integer(c_int), intent(inout) :: po_distribution, po_maxcycle, po_init_pop_size, po_reset
      real(c_double), intent(inout) :: po_mutation_rate, po_death_rate, po_scalefac
      integer(c_int), intent(inout) :: po_nsave, ntasks, tdlf_farm, n_po_scaling
      real(c_double), intent(inout) :: neb_climb_test, neb_freeze_test
      integer(c_int), intent(inout) :: nzero, coupled_states, qtsflag
      integer(c_int), intent(inout) :: imicroiter, maxmicrocycle, micro_esp_fit

      integer :: i

      ierr = 0_c_int

      coords(1:nvar) = reshape(initial_coords, [3*n_atoms])

      ! [residue(nat), znuc(nat), micspec(nat)]. The first block is not a
      ! frozen/free flag, which is what it looks like and what this passed at
      ! first: a negative entry freezes an atom, but a positive one is the
      ! *residue* the atom belongs to. All zeros means no atom is in any
      ! residue, which Cartesians ignore and HDLC cannot survive -- it builds
      ! internal coordinates per residue, found none, and segfaulted.
      !
      ! The residues are this program's monomers, which is the mapping that
      ! makes HDLC worth having: internals inside each molecule, Cartesians
      ! between them. Pure DLC is the other extreme, one residue over
      ! everything, and it fails outright on a cluster ("cyclic failure at
      ! residue 1") because the molecules are not bonded to each other and
      ! there is no connected internal-coordinate system to build.
      !
      ! micspec must be 0 or 1 and DL-FIND validates it, so 1 (this atom takes
      ! part) is the only safe fill whether or not microiterations run.
      do i = 1, n_atoms
         spec(i) = int(atom_residues(i), c_int)
      end do
      do i = 1, n_atoms
         spec(n_atoms + i) = int(atomic_numbers(i), c_int)
      end do
      spec(2*n_atoms + 1:3*n_atoms) = 1_c_int

      nz = int(n_atoms, c_int)
      ncons = 0_c_int
      nconn = 0_c_int
      tatoms = 1_c_int  !! The variables really are atoms; HDLC needs to know

      ! Not optional, whatever the comment above about defaults says.
      ! `dlf_default_init` sets most of these arguments and then dlf_read_in
      ! sets nz and nweight itself -- but nframe, nmass and n_po_scaling are
      ! plain uninitialized locals in dlf_read_in, passed here intent(inout),
      ! and ntasks is initialized to -1 on purpose ("why??? no real need for
      ! this not to be 1" is the comment in DL-FIND's own source) and then
      ! checked for being positive. Leaving any of them alone means reading an
      ! uninitialized integer: nframe and nmass are validated against
      ! nvarin2 = nframe*nat*3 + nweight + nmass, which is 0 here, so a garbage
      ! value is a hard failure with a confusing message, and ntasks at -1 is
      ! the "Number of task farms must be positive" refusal.
      !
      ! One task farm: this program does its own MPI, and DL-FIND is linked
      ! against its serial stubs. Splitting the job again underneath would be
      ! two schedulers dividing the same ranks.
      ntasks = 1_c_int
      nframe = 0_c_int
      nmass = 0_c_int
      nweight = 0_c_int

      icoord = dlfind_coordinates(settings%coordinates)
      iopt = dlfind_algorithm(settings%algorithm)

      tolerance = real(settings%gradient_tolerance, c_double)
      maxcycle = int(settings%max_steps, c_int)

      ! Negative in the settings means "DL-FIND's own default", so these are
      ! the only ones that touch the argument at all.
      if (settings%energy_tolerance > 0.0_dp) tolerance_e = real(settings%energy_tolerance, c_double)
      if (settings%max_step > 0.0_dp) maxstep = real(settings%max_step, c_double)
      if (settings%lbfgs_memory > 0) lbfgs_mem = int(settings%lbfgs_memory, c_int)
      printl = dlfind_print_level()

      ! printf >= 1 is what makes DL-FIND report geometries through
      ! dlf_put_coords. Without it the optimization runs to convergence and
      ! hands back nothing, which is the one failure that looks like success.
      printf = 1_c_int

      ! No restart files. DL-FIND writes them into the working directory under
      ! fixed names, so two optimizations sharing a directory would read each
      ! other's -- and this program checkpoints through mqc_checkpoint anyway.
      dump = 0_c_int
      restart = 0_c_int

   end subroutine cb_get_params

   subroutine cb_get_gradient(nvar, coords, energy, gradient, iimage, kiter, status) bind(c)
      !! DL-FIND asking for the energy and gradient at a geometry
      integer(c_int), intent(in), value :: nvar
      real(c_double), intent(in) :: coords(nvar)
      real(c_double), intent(out) :: energy
      real(c_double), intent(out) :: gradient(nvar)
      integer(c_int), intent(in), value :: iimage  !! Image index; single-image here
      integer(c_int), intent(in), value :: kiter   !! Microiteration counter, unused
      integer(c_int), intent(out) :: status

      real(dp) :: step_energy
      real(dp) :: step_gradient(3, n_atoms)
      integer :: step_status

      call evaluator(n_atoms, reshape(real(coords, dp), [3, n_atoms]), &
                     step_energy, step_gradient, step_status)

      if (step_status /= 0) then
         ! Recorded as well as reported: DL-FIND may stop on a non-zero status
         ! or may try a shorter step, and either way the run is no longer
         ! optimizing what the deck asked for.
         evaluation_failed = .true.
         energy = 0.0_c_double
         gradient = 0.0_c_double
         status = 1_c_int
         return
      end if

      energy = real(step_energy, c_double)
      gradient = reshape(real(step_gradient, c_double), [nvar])
      status = 0_c_int

   end subroutine cb_get_gradient

   subroutine cb_put_coords(nvar, switch, energy, coords, iam) bind(c)
      !! DL-FIND reporting a geometry back
      !!
      !! `switch` is 1 for a geometry along the way and 3 for a final one. The
      !! last of either is the answer: DL-FIND tests convergence after
      !! evaluating, so the last geometry it reports is the converged one.
      integer(c_int), intent(in), value :: nvar
      integer(c_int), intent(in), value :: switch
      real(c_double), intent(in), value :: energy
      real(c_double), intent(in) :: coords(nvar)
      integer(c_int), intent(in), value :: iam  !! Task id; only the first reports

      if (iam /= 0_c_int) return
      if (switch /= 1_c_int .and. switch /= 3_c_int) return
      if (nvar /= 3*n_atoms) return

      latest_coords = reshape(real(coords, dp), [3, n_atoms])
      latest_energy = real(energy, dp)
      have_result = .true.

      ! Reported here rather than from the gradient callback on purpose: this
      ! fires for a geometry the optimizer accepted, while the gradient
      ! callback also fires for the trial points of a line search. A trajectory
      ! built from the latter doubles back on itself.
      call on_step(n_atoms, latest_coords, latest_energy)

   end subroutine cb_put_coords

   subroutine cb_get_hessian(nvar, coords, hessian, status) bind(c)
      !! Declining to supply an analytic Hessian
      !!
      !! A non-zero status is how DL-FIND is told to build its own by finite
      !! difference or by update, which is what the minimizers here want. This
      !! program can compute a Hessian, but doing it per optimization step
      !! costs 6N gradients for a curvature that L-BFGS approximates for free.
      integer(c_int), intent(in), value :: nvar
      real(c_double), intent(in) :: coords(nvar)
      real(c_double), intent(out) :: hessian(nvar, nvar)
      integer(c_int), intent(out) :: status

      hessian = 0.0_c_double
      status = 1_c_int

   end subroutine cb_get_hessian

   subroutine cb_get_multistate_gradients(nvar, coords, energy, gradient, coupling, &
                                          needcoupling, iimage, status) bind(c)
      !! Conical intersection search, which this bridge does not offer
      integer(c_int), intent(in), value :: nvar
      real(c_double), intent(in) :: coords(nvar)
      real(c_double), intent(out) :: energy(2)
      real(c_double), intent(out) :: gradient(nvar, 2)
      real(c_double), intent(out) :: coupling(nvar)
      integer(c_int), intent(in), value :: needcoupling
      integer(c_int), intent(in), value :: iimage
      integer(c_int), intent(out) :: status

      energy = 0.0_c_double
      gradient = 0.0_c_double
      coupling = 0.0_c_double
      status = 1_c_int

   end subroutine cb_get_multistate_gradients

   subroutine cb_error() bind(c)
      !! DL-FIND's last word before libdlfind executes `error stop`
      !!
      !! There is no returning from this: api.f90 stops the process as soon as
      !! this comes back. Say what was being optimized while there is still
      !! somewhere to say it.
      call logger%error("DL-FIND failed and is about to stop the process.")
      call logger%error("  Atoms: "//to_char(n_atoms))
      call logger%error("  Energy evaluations completed: see the step log above.")

   end subroutine cb_error

   subroutine cb_update() bind(c)
      !! Called between steps. Nothing here needs to happen then.
   end subroutine cb_update

   function dlfind_print_level() result(printl)
      !! How loud DL-FIND should be, from how loud this program was asked to be
      !!
      !! DL-FIND writes with `write(stdout,...)` throughout and cannot be
      !! routed into pic's logger without editing every one of its thirty
      !! source files, so the next best thing is to make its own verbosity dial
      !! follow ours. `keywords.optimization.print_level` still overrides this
      !! outright, which is the escape hatch for debugging an optimization.
      !!
      !! The default is quiet at info level rather than DL-FIND's own 2,
      !! because this program already prints a line per step with the energy
      !! and the largest gradient component. DL-FIND's convergence table on top
      !! of that is the same information a second time, six lines per step.
      use pic_logger, only: debug_level, verbose_level, info_level
      integer(c_int) :: printl

      if (settings%print_level >= 0) then
         printl = int(settings%print_level, c_int)
      else if (logger%log_level >= debug_level) then
         printl = 4_c_int
      else if (logger%log_level >= verbose_level) then
         printl = 2_c_int
      else if (logger%log_level >= info_level) then
         printl = 0_c_int
      else
         printl = 0_c_int
      end if

   end function dlfind_print_level

   pure function dlfind_coordinates(coordinates) result(icoord)
      !! This program's coordinate-system numbering, in DL-FIND's
      integer, intent(in) :: coordinates
      integer(c_int) :: icoord

      select case (coordinates)
      case (OPT_COORDS_HDLC)
         icoord = DLF_COORDS_HDLC
      case (OPT_COORDS_DLC)
         icoord = DLF_COORDS_DLC
      case (OPT_COORDS_HDLC_TC)
         icoord = DLF_COORDS_HDLC_TC
      case (OPT_COORDS_DLC_TC)
         icoord = DLF_COORDS_DLC_TC
      case default
         icoord = DLF_COORDS_CARTESIAN
      end select
   end function dlfind_coordinates

   pure function dlfind_algorithm(algorithm) result(iopt)
      !! This program's algorithm numbering, in DL-FIND's
      integer, intent(in) :: algorithm
      integer(c_int) :: iopt

      select case (algorithm)
      case (OPT_ALGO_SD)
         iopt = DLF_OPT_SD
      case (OPT_ALGO_CG)
         iopt = DLF_OPT_CG
      case (OPT_ALGO_PRFO)
         iopt = DLF_OPT_PRFO
      case default
         iopt = DLF_OPT_LBFGS
      end select
   end function dlfind_algorithm

end module mqc_dlfind_bridge
