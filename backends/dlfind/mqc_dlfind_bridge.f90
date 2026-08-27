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
                                  OPT_ALGO_SD, OPT_ALGO_CG, OPT_ALGO_LBFGS, OPT_ALGO_PRFO, &
                                  OPT_ALGO_CG_AUTO, OPT_ALGO_NR, OPT_ALGO_DAMPED, &
                                  OPT_HESSIAN_UPDATE_ENGINE, hessian_i, &
                                  algorithm_needs_hessian, constraint_atom_count, &
                                  NEB_ENDS_FREE, NEB_ENDS_PERPENDICULAR, &
                                  SADDLE_METHOD_DIMER
   use mqc_error, only: error_t, ERROR_GENERIC, ERROR_VALIDATION
   implicit none
   private

   public :: dlfind_available
   public :: dlfind_optimize
   public :: dlfind_connected_minima

   ! DL-FIND's own numbering, which this file is the only place to know.
   ! `icoord`: mod 10 selects the coordinate system, and the tens and hundreds
   ! digits select multi-image methods this bridge does not offer.
   integer(c_int), parameter :: DLF_COORDS_CARTESIAN = 0
   integer(c_int), parameter :: DLF_COORDS_HDLC = 1
   integer(c_int), parameter :: DLF_COORDS_DLC = 3
   integer(c_int), parameter :: DLF_COORDS_HDLC_TC = 2
   integer(c_int), parameter :: DLF_COORDS_DLC_TC = 4

   ! Chain-of-states runs are selected through `icoord`, not through `iopt`:
   ! DL-FIND reads the unit place as the coordinate system and the hundreds and
   ! tens as what the band is. 10X is NEB with free endpoints, 11X endpoints
   ! moving only perpendicular to their tangent, 12X frozen endpoints. 190 is
   ! not NEB at all -- it is the quantum transition state search -- which is why
   ! the offsets stop at 120 and are added to a unit place rather than used
   ! whole. The implementation is improved-tangent NEB with a climbing image.
   ! The dimer occupies the 200 band, with the tens place choosing how the pair
   ! is rotated: 20X leaves translation and rotation to the optimizer named by
   ! `iopt`, 21X and 22X do the rotation by a line search inside DL-FIND's own
   ! dimer module at the cost of extra energy evaluations per rotation. 20X is
   ! taken here because it is the one that respects the algorithm the deck
   ! chose.
   integer(c_int), parameter :: DLF_DIMER = 200

   ! DL-FIND's task numbers. 1011 is "TS search, then find both downhill
   ! structures": it runs the saddle search as configured, stores the geometry
   ! and the imaginary mode, then displaces along that mode in each direction
   ! and minimises, which is the practical answer to what the saddle connects.
   ! DL-FIND has no IRC -- there is no `dlf_irc.f90` -- so this is the closest
   ! thing to one it offers, and it answers the question an IRC is usually run
   ! to answer.
   integer(c_int), parameter :: DLF_TASK_TS_AND_DOWNHILL = 1011

   ! Which structure `cb_put_coords` is being handed. 1 is an ordinary step, 2
   ! the transition mode, 3 the transition structure, and 4 and 5 the two minima
   ! reached by relaxing downhill from it.
   integer(c_int), parameter :: DLF_PUT_STEP = 1
   integer(c_int), parameter :: DLF_PUT_TS = 3
   integer(c_int), parameter :: DLF_PUT_DOWNHILL_A = 4
   integer(c_int), parameter :: DLF_PUT_DOWNHILL_B = 5

   integer(c_int), parameter :: DLF_NEB_ENDS_FREE = 100
   integer(c_int), parameter :: DLF_NEB_ENDS_PERPENDICULAR = 110
   integer(c_int), parameter :: DLF_NEB_ENDS_FROZEN = 120

   integer(c_int), parameter :: DLF_OPT_SD = 0
   integer(c_int), parameter :: DLF_OPT_CG_AUTO = 1
   integer(c_int), parameter :: DLF_OPT_CG = 2
      !! Two spellings of Polak-Ribiere: 1 restarts on the Powell-Beale test,
      !! 2 every ten steps. DL-FIND's own source calls 2 "better
      !! implementation", which is why it is what `cg` has always meant here.
   integer(c_int), parameter :: DLF_OPT_LBFGS = 3
   integer(c_int), parameter :: DLF_OPT_PRFO = 10
   integer(c_int), parameter :: DLF_OPT_NR = 20
   integer(c_int), parameter :: DLF_OPT_DAMPED = 30

   ! `inithessian`: where the first Hessian comes from. 0 asks the external
   ! program -- this bridge -- and DL-FIND falls back to two-point finite
   ! differences by itself if the callback declines, so 0 is safe to set even
   ! when nothing can supply one.
   integer(c_int), parameter :: DLF_INITHESS_EXTERNAL = 0
   integer(c_int), parameter :: DLF_INITHESS_TWO_POINT = 2

   integer(c_int), parameter :: DLF_SPEC_FROZEN = -1
      !! A `spec` entry of -1 removes the atom's three coordinates from the
      !! optimization entirely. DL-FIND also has -2/-3/-4 for freezing single
      !! Cartesian components and -23/-24/-34 for pairs; only the whole-atom
      !! case is offered here, because the others are a coordinate-frame choice
      !! disguised as a chemistry one -- they mean nothing unless the caller
      !! knows which way the input axes point.

   ! The optimization in progress. See the note above on why this is here.
   procedure(energy_gradient_i), pointer, save :: evaluator => null()
   procedure(hessian_i), pointer, save :: hessian_evaluator => null()
      !! Null when the caller has none to offer, which is not an error: the
      !! Hessian callback then declines and DL-FIND builds its own.
   procedure(step_callback_i), pointer, save :: on_step => null()
   type(optimizer_settings_t), save :: settings
   integer, save :: n_atoms = 0
   integer, allocatable, save :: atomic_numbers(:)
   integer, allocatable, save :: atom_residues(:)
   real(dp), allocatable, save :: initial_coords(:, :)
   real(dp), allocatable, save :: downhill_a(:, :), downhill_b(:, :)
   real(dp), save :: energy_a = 0.0_dp, energy_b = 0.0_dp
   real(dp), save :: saddle_energy = 0.0_dp
   logical, save :: have_saddle = .false.
      !! The transition structure's own energy. A `connect` run does not stop
      !! there -- it carries on downhill twice -- so the energy the optimizer
      !! finishes with belongs to the last minimum, not to the saddle.
   logical, save :: have_downhill_a = .false., have_downhill_b = .false.
      !! The two minima a `connect` run falls to, kept for the caller to report.
   real(dp), allocatable, save :: endpoint_coords(:, :)
      !! The second structure of a chain-of-states run, held here for the same
      !! reason the geometry is: `cb_get_params` is a C callback and takes no
      !! context of its own.
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
                              energy_gradient, step_taken, final_energy, error, &
                              hessian, endpoint)
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
      procedure(hessian_i), optional :: hessian
      real(dp), intent(in), optional :: endpoint(:, :)
         !! The product geometry, Bohr, same atom order as `coords`. Present
         !! turns this into a chain-of-states run.
         !! Second derivatives, for the algorithms that climb rather than
         !! descend. Absent is not a refusal: DL-FIND falls back to two-point
         !! finite differences of the gradients it is already asking for, which
         !! costs `6N` of them per Hessian and gets the same answer.

      integer(c_int) :: nvar, nspec, nvar2_in
      integer :: k, n_cons_atoms, ncons_here

      final_energy = 0.0_dp

      if (natoms < 1) then
         call error%set(ERROR_VALIDATION, "dlfind_optimize: nothing to optimize")
         return
      end if

      ! Frozen atoms are written into `spec` as negative entries, and pure
      ! internal coordinates cannot express one: DL-FIND checks for it and
      ! calls `dlf_fail`, which ends the process. Refused here so it is a
      ! message rather than a Fortran backtrace.
      if (allocated(opt_settings%frozen_atoms)) then
         if (size(opt_settings%frozen_atoms) > 0 .and. &
             (opt_settings%coordinates == OPT_COORDS_DLC .or. &
              opt_settings%coordinates == OPT_COORDS_DLC_TC)) then
            call error%set(ERROR_VALIDATION, &
                           "keywords.optimization.frozen_atoms needs a coordinate system "// &
                           "that can hold an atom still. Pure internals cannot -- use "// &
                           '"hdlc" or "cartesian".')
            return
         end if
         if (any(opt_settings%frozen_atoms < 1) .or. &
             any(opt_settings%frozen_atoms > natoms)) then
            call error%set(ERROR_VALIDATION, &
                           "keywords.optimization.frozen_atoms names an atom outside the "// &
                           "system.")
            return
         end if
      end if

      ! A constraint is a primitive internal coordinate held fixed, so there
      ! has to be an internal coordinate system holding it. In Cartesians
      ! DL-FIND has nowhere to put one and ignores it, which is the failure
      ! that looks like success -- the run converges to the unconstrained
      ! answer and says nothing.
      if (allocated(opt_settings%constraints)) then
         if (size(opt_settings%constraints) > 0 .and. &
             opt_settings%coordinates == OPT_COORDS_CARTESIAN) then
            call error%set(ERROR_VALIDATION, &
                           "keywords.optimization.constraints needs an internal coordinate "// &
                           'system to constrain. Use "hdlc"; in cartesian the constraint '// &
                           "would be silently ignored.")
            return
         end if
         do k = 1, size(opt_settings%constraints)
            n_cons_atoms = constraint_atom_count(opt_settings%constraints(k)%kind)
            if (n_cons_atoms == 0) then
               call error%set(ERROR_VALIDATION, &
                              "keywords.optimization.constraints has an entry whose type "// &
                              "is not one this program knows.")
               return
            end if
            if (any(opt_settings%constraints(k)%atoms(1:n_cons_atoms) < 1) .or. &
                any(opt_settings%constraints(k)%atoms(1:n_cons_atoms) > natoms)) then
               call error%set(ERROR_VALIDATION, &
                              "keywords.optimization.constraints names an atom outside "// &
                              "the system.")
               return
            end if
         end do
      end if

      settings = opt_settings
      n_atoms = natoms
      evaluator => energy_gradient
      on_step => step_taken

      ! Only for the algorithms that hold one. Pointing it at a live evaluator
      ! for L-BFGS would be harmless -- nothing would call it -- but it would
      ! also be a lie about what the run is going to do.
      nullify (hessian_evaluator)
      if (present(hessian) .and. algorithm_needs_hessian(opt_settings%algorithm)) then
         hessian_evaluator => hessian
      end if
      have_result = .false.
      evaluation_failed = .false.

      if (allocated(atomic_numbers)) deallocate (atomic_numbers)
      if (allocated(atom_residues)) deallocate (atom_residues)
      if (allocated(initial_coords)) deallocate (initial_coords)
      if (allocated(latest_coords)) deallocate (latest_coords)
      allocate (atomic_numbers(natoms), source=znuc)
      allocate (atom_residues(natoms), source=residues)
      allocate (initial_coords(3, natoms), source=coords)
      ! Held for `cb_get_params`, and cleared when absent: this module's state is
      ! `save`, so a second optimization in the same process would otherwise
      ! inherit the previous run's path and quietly become a chain-of-states run.
      have_downhill_a = .false.
      have_downhill_b = .false.
      have_saddle = .false.
      if (allocated(downhill_a)) deallocate (downhill_a)
      if (allocated(downhill_b)) deallocate (downhill_b)
      if (allocated(endpoint_coords)) deallocate (endpoint_coords)
      if (present(endpoint)) allocate (endpoint_coords(3, natoms), source=endpoint)
      allocate (latest_coords(3, natoms), source=coords)

      nvar = int(3*natoms, c_int)

      ! nspec is fixed by DL-FIND as nat + nz + 5*ncons + 2*nconn + nat, and it
      ! checks: a mismatch is refused as an interface error rather than read
      ! past. With nz = nat and no constraints that is 3*nat, laid out as
      ! [frozen(nat), znuc(nat), micspec(nat)]. Passing znuc rather than
      ! leaving nz at zero is what lets HDLC and DLC perceive connectivity --
      ! without it those coordinate systems have no elements to work from.
      ! [residue(nat), znuc(nat), icons(5*ncons), iconn(2*nconn), micspec(nat)].
      ! The constraint block sits between the charges and the microiterative
      ! flags, so adding constraints moves micspec -- which is why the offset
      ! below is computed rather than written as 2*nat.
      !
      ! `nconn` stays zero. DL-FIND reads `iconn` from `nat+nz+ncons` while it
      ! reads `micspec` from `nat+nz+5*ncons+2*nconn`; the two disagree about
      ! how long the constraint block is, so a deck using both would have its
      ! connections read out of the middle of the constraints. Constraints are
      ! the half worth having, so they are the half that is offered.
      ncons_here = 0
      if (allocated(opt_settings%constraints)) ncons_here = size(opt_settings%constraints)
      nspec = int(3*natoms + 5*ncons_here, c_int)

      ! The second coordinate array carries NEB images, a reaction endpoint or
      ! per-atom masses and weights. A minimization uses none of them and
      ! passes zero, and DL-FIND still hands the callback a length of
      ! max(nvarin2,1). A chain-of-states run passes one frame: the product.
      ! DL-FIND interpolates the images between it and the geometry it was
      ! given, so a path is two structures and a count rather than a file of
      ! guesses.
      nvar2_in = 0_c_int
      if (allocated(endpoint_coords)) nvar2_in = int(3*natoms, c_int)
      if (opt_settings%saddle_method == SADDLE_METHOD_DIMER .and. &
          .not. allocated(endpoint_coords)) nvar2_in = 0_c_int

      call api_dl_find(nvar, nvar2_in, nspec, 1_c_int, &
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

      integer :: i, k, ic, ispec, n_constraints

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

      ! A frozen atom's residue is overwritten rather than annotated: the entry
      ! is one number and DL-FIND reads a negative one as "held still". So a
      ! frozen atom belongs to no residue, which is DL-FIND's model and not a
      ! choice made here -- it has no coordinates for a residue to contain.
      if (allocated(settings%frozen_atoms)) then
         do i = 1, size(settings%frozen_atoms)
            spec(settings%frozen_atoms(i)) = DLF_SPEC_FROZEN
         end do
      end if

      do i = 1, n_atoms
         spec(n_atoms + i) = int(atomic_numbers(i), c_int)
      end do

      ! [type, atom1..atom4], 1-based, unused atom slots zero. DL-FIND reshapes
      ! this block column-major into icons(5, ncons), so the five entries of one
      ! constraint are contiguous.
      ispec = 2*n_atoms
      n_constraints = 0
      if (allocated(settings%constraints)) n_constraints = size(settings%constraints)
      if (allocated(settings%constraints)) then
         do k = 1, size(settings%constraints)
            spec(ispec + 1) = int(settings%constraints(k)%kind, c_int)
            do ic = 1, 4
               spec(ispec + 1 + ic) = int(settings%constraints(k)%atoms(ic), c_int)
            end do
            ispec = ispec + 5
         end do
      end if

      spec(ispec + 1:ispec + n_atoms) = 1_c_int

      nz = int(n_atoms, c_int)
      ncons = int(n_constraints, c_int)
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
      if (settings%connect) then
         task = DLF_TASK_TS_AND_DOWNHILL
         if (settings%connect_distort > 0.0_dp) then
            distort = real(settings%connect_distort, c_double)
         end if
      end if

      nframe = 0_c_int
      nmass = 0_c_int
      nweight = 0_c_int

      ! A band only when a second structure was actually given. The deck asking
      ! for images or a spring constant without an endpoint is not a path, and
      ! silently running one from a single geometry would optimize something
      ! nobody described.
      if (settings%saddle_method == SADDLE_METHOD_DIMER) then
         ! The dimer needs a direction to start from, and randomises one when
         ! it is given nothing -- which on any real system means starting the
         ! search along an arbitrary direction in 3N space. An endpoint, when
         ! there is one, points the pair along the reaction instead.
         !
         ! Handed over as the geometry, not as the displacement to it. DL-FIND
         ! reads `coords2` as an absolute structure and takes the difference
         ! itself unless `tsrelative` says otherwise, so passing a difference
         ! while leaving that flag alone has it subtract the coordinates a
         ! second time. The resulting axis points somewhere arbitrary, the
         ! first translation lands on a geometry with no SCF, and the run dies
         ! as "energy evaluation failed" several cycles later.
         icoord = dlfind_coordinates(settings%coordinates, DLF_DIMER)
         if (settings%dimer_separation > 0.0_dp) delta = real(settings%dimer_separation, c_double)
         if (settings%dimer_max_rotations > 0) maxrot = int(settings%dimer_max_rotations, c_int)
         if (settings%dimer_rotation_tolerance > 0.0_dp) then
            tolrot = real(settings%dimer_rotation_tolerance, c_double)
         end if
         if (allocated(endpoint_coords)) then
            nframe = 1_c_int
            coords2(1:3*n_atoms) = reshape(endpoint_coords, [3*n_atoms])
         end if
      else if (allocated(endpoint_coords)) then
         icoord = dlfind_coordinates(settings%coordinates, dlfind_neb_band(settings%neb_ends))
         if (settings%n_images > 0) nimage = int(settings%n_images, c_int)
         if (settings%neb_spring >= 0.0_dp) nebk = real(settings%neb_spring, c_double)
         nframe = 1_c_int
         coords2(1:3*n_atoms) = reshape(endpoint_coords, [3*n_atoms])
      else
         icoord = dlfind_coordinates(settings%coordinates)
      end if
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
      ! A `connect` run needs 4. DL-FIND hands back the transition structure and
      ! the two minima it relaxed to only from inside `printf >= 4` blocks --
      ! the same blocks that write TS.xyz and minimum_+.xyz -- so at the default
      ! the downhill runs happen, report themselves as finished, and are never
      ! passed to the caller. It costs a handful of files in the working
      ! directory, which is a fair price for the answer being reachable.
      if (settings%connect) printf = 4_c_int

      ! Curvature settings, which only P-RFO and Newton-Raphson read. Asking
      ! for the external Hessian is safe whether or not one can be produced:
      ! the callback declines and DL-FIND says "External Hessian not
      ! available, using two point FD" and carries on.
      if (algorithm_needs_hessian(settings%algorithm)) then
         if (associated(hessian_evaluator)) then
            inithessian = DLF_INITHESS_EXTERNAL
         else
            inithessian = DLF_INITHESS_TWO_POINT
         end if
         if (settings%hessian_update /= OPT_HESSIAN_UPDATE_ENGINE) then
            update = int(settings%hessian_update, c_int)
         end if
      end if

      ! Damped dynamics only. Negative in the settings means the engine's own,
      ! as everywhere else here.
      if (settings%algorithm == OPT_ALGO_DAMPED) then
         if (settings%timestep > 0.0_dp) timestep = real(settings%timestep, c_double)
         if (settings%friction >= 0.0_dp) fric0 = real(settings%friction, c_double)
         if (settings%friction_factor > 0.0_dp) fricfac = real(settings%friction_factor, c_double)
         if (settings%friction_rising >= 0.0_dp) fricp = real(settings%friction_rising, c_double)
      end if

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
      ! The two downhill minima arrive here and nowhere else, so they are kept
      ! rather than filtered out with everything that is not a step.
      if (nvar == 3*n_atoms) then
         if (switch == DLF_PUT_DOWNHILL_A) then
            if (.not. allocated(downhill_a)) allocate (downhill_a(3, n_atoms))
            downhill_a = reshape(coords, [3, n_atoms])
            energy_a = energy
            have_downhill_a = .true.
            return
         end if
         if (switch == DLF_PUT_DOWNHILL_B) then
            if (.not. allocated(downhill_b)) allocate (downhill_b(3, n_atoms))
            downhill_b = reshape(coords, [3, n_atoms])
            energy_b = energy
            have_downhill_b = .true.
            return
         end if
      end if
      if (switch == DLF_PUT_TS .and. nvar == 3*n_atoms) then
         saddle_energy = energy
         have_saddle = .true.
      end if
      if (switch /= DLF_PUT_STEP .and. switch /= DLF_PUT_TS) return
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
      !! DL-FIND asking for second derivatives at a geometry
      !!
      !! **Declining is a supported answer.** A non-zero status tells DL-FIND
      !! to build its own by two-point finite difference of the gradients it is
      !! already requesting, and it says so and carries on. So this returns 1
      !! whenever there is no evaluator rather than failing the run -- a
      !! minimiser never asks, and an algorithm that does asks for something it
      !! can get another way.
      !!
      !! Cartesian on the way out. DL-FIND calls `dlf_coords_hessian_xtoi` on
      !! whatever it gets, so handing it a matrix already in internals would be
      !! transformed a second time.
      integer(c_int), intent(in), value :: nvar
      real(c_double), intent(in) :: coords(nvar)
      real(c_double), intent(out) :: hessian(nvar, nvar)
      integer(c_int), intent(out) :: status

      real(dp), allocatable :: geometry(:, :), h(:, :)
      integer :: local_status

      hessian = 0.0_c_double
      status = 1_c_int

      if (.not. associated(hessian_evaluator)) return
      if (nvar /= 3*n_atoms) return

      allocate (geometry(3, n_atoms), h(3*n_atoms, 3*n_atoms))
      geometry = reshape(real(coords, dp), [3, n_atoms])

      call hessian_evaluator(n_atoms, geometry, h, local_status)
      if (local_status == 0) then
         hessian = real(h, c_double)
         status = 0_c_int
      else
         ! Left at 1. The evaluation that failed is not this optimization's to
         ! recover from, and DL-FIND's own fallback is a correct Hessian by a
         ! slower route rather than a degraded one.
         call logger%verbose("dlfind: no analytic Hessian at this geometry, "// &
                             "DL-FIND will build one by finite difference")
      end if

      deallocate (geometry, h)

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

   pure function dlfind_coordinates(coordinates, band) result(icoord)
      !! This program's coordinate-system numbering, in DL-FIND's
      !!
      !! `band` is added on top: zero for an ordinary single-structure run, and
      !! one of the `DLF_NEB_ENDS_*` offsets for a chain of states. The unit
      !! place stays the coordinate system either way, which is exactly how
      !! DL-FIND decodes it.
      integer, intent(in) :: coordinates
      integer(c_int), intent(in), optional :: band
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
      if (present(band)) icoord = icoord + band
   end function dlfind_coordinates

   subroutine dlfind_connected_minima(a, energy_first, b, energy_second, found, e_saddle)
      !! The two minima a `connect` run relaxed to, if it ran and reached them
      real(dp), allocatable, intent(out) :: a(:, :), b(:, :)
      real(dp), intent(out) :: energy_first, energy_second
      logical, intent(out) :: found
      real(dp), intent(out) :: e_saddle
         !! The saddle's own energy, which is not the one the run finishes with.

      e_saddle = saddle_energy
      found = have_downhill_a .and. have_downhill_b .and. have_saddle
      energy_first = energy_a
      energy_second = energy_b
      if (.not. found) return
      a = downhill_a
      b = downhill_b
   end subroutine dlfind_connected_minima

   pure function dlfind_neb_band(ends) result(band)
      !! This program's endpoint treatment, as DL-FIND's `icoord` offset
      integer, intent(in) :: ends
      integer(c_int) :: band

      select case (ends)
      case (NEB_ENDS_FREE)
         band = DLF_NEB_ENDS_FREE
      case (NEB_ENDS_PERPENDICULAR)
         band = DLF_NEB_ENDS_PERPENDICULAR
      case default
         band = DLF_NEB_ENDS_FROZEN
      end select
   end function dlfind_neb_band

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
      case (OPT_ALGO_CG_AUTO)
         iopt = DLF_OPT_CG_AUTO
      case (OPT_ALGO_NR)
         iopt = DLF_OPT_NR
      case (OPT_ALGO_DAMPED)
         iopt = DLF_OPT_DAMPED
      case default
         iopt = DLF_OPT_LBFGS
      end select
   end function dlfind_algorithm

end module mqc_dlfind_bridge
