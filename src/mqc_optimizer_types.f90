!! What a geometry optimization needs to know, independent of the optimizer
module mqc_optimizer_types
   !! The vocabulary shared by the optimization driver and whichever engine
   !! takes the steps. It knows about neither: no MPI, no DL-FIND, no quantum
   !! chemistry. That is what lets `mqc_geometry_optimizer` and
   !! `mqc_dlfind_bridge` reach each other without a circular dependency -- the
   !! driver hands the engine a procedure matching `energy_gradient_i`, and the
   !! engine calls it without ever learning what is on the other end.
   !!
   !! The same arrangement is what makes the optimizer backend-agnostic. The
   !! evaluator behind that interface is `run_calculation` with a gradient
   !! request, so xTB through tblite, Hartree-Fock through libcint and HF/DFT
   !! through cuEST all arrive here as the same three arrays, fragmented or
   !! not, and no engine has to be taught about any of them.
   use pic_types, only: dp
   use mqc_calculation_defaults, only: DEFAULT_OPT_MAX_STEPS, &
                                       DEFAULT_OPT_GRADIENT_TOLERANCE, &
                                       DEFAULT_OPT_ENERGY_TOLERANCE, &
                                       DEFAULT_OPT_MAX_STEP, &
                                       DEFAULT_OPT_LBFGS_MEMORY, &
                                       DEFAULT_OPT_PRINT_LEVEL
   implicit none
   private

   public :: optimizer_settings_t
   public :: opt_constraint_t
   public :: energy_gradient_i
   public :: hessian_i
   public :: step_callback_i
   public :: OPT_COORDS_UNKNOWN, OPT_COORDS_CARTESIAN, OPT_COORDS_HDLC, OPT_COORDS_DLC
   public :: OPT_COORDS_HDLC_TC, OPT_COORDS_DLC_TC
   public :: OPT_ALGO_UNKNOWN, OPT_ALGO_SD, OPT_ALGO_CG, OPT_ALGO_LBFGS, OPT_ALGO_PRFO
   public :: OPT_ALGO_CG_AUTO, OPT_ALGO_NR, OPT_ALGO_DAMPED
   public :: OPT_HESSIAN_UPDATE_ENGINE, OPT_HESSIAN_UPDATE_NONE
   public :: OPT_HESSIAN_UPDATE_POWELL, OPT_HESSIAN_UPDATE_BOFILL
   public :: OPT_CONSTRAIN_BOND, OPT_CONSTRAIN_ANGLE, OPT_CONSTRAIN_TORSION
   public :: OPT_CONSTRAIN_CARTESIAN, OPT_CONSTRAIN_BOND_DIFFERENCE
   public :: coordinates_from_string, coordinates_to_string
   public :: algorithm_from_string, algorithm_to_string
   public :: hessian_update_from_string, hessian_update_to_string
   public :: constraint_from_string, constraint_atom_count
   public :: OPT_TARGET_MINIMUM, OPT_TARGET_SADDLE
   public :: target_from_string, target_to_string, algorithm_finds_saddle
   public :: NEB_ENDS_FROZEN, NEB_ENDS_PERPENDICULAR, NEB_ENDS_FREE, MIN_NEB_IMAGES
   public :: neb_ends_from_string, neb_ends_to_string
   public :: SADDLE_METHOD_PRFO, SADDLE_METHOD_DIMER
   public :: saddle_method_from_string, saddle_method_to_string
   public :: algorithm_needs_hessian

   ! This program's own numbering rather than DL-FIND's, even though DL-FIND is
   ! the only engine today and two of these would have mapped to themselves.
   ! The translation lives in the bridge so that a deck saying "hdlc" keeps
   ! meaning HDLC if a second engine numbers them differently, and so that an
   ! engine which cannot do internal coordinates at all refuses in its own
   ! terms rather than silently optimizing in Cartesians.
   integer, parameter :: OPT_COORDS_UNKNOWN = 0
   integer, parameter :: OPT_COORDS_CARTESIAN = 1
   integer, parameter :: OPT_COORDS_HDLC = 2  !! Hybrid delocalised internal coordinates
   integer, parameter :: OPT_COORDS_DLC = 3   !! Delocalised internal coordinates
   integer, parameter :: OPT_COORDS_HDLC_TC = 4
      !! HDLC over a total connection scheme
      !!
      !! The primitive internals of a residue are built from every atom pair in
      !! it rather than from the perceived bonds. That is what makes it useful on
      !! a system whose connectivity is not what a covalent-radius rule says: a
      !! weakly bound complex the perception splits, an ion pair it merges, a
      !! transition state where a bond is half formed. The cost is quadratic in
      !! residue size where a bond list is linear, so it is for small residues.
   integer, parameter :: OPT_COORDS_DLC_TC = 5
      !! DLC over a total connection scheme
      !!
      !! One residue, totally connected. Unlike plain `dlc` this does not fail on
      !! a cluster -- the total connection supplies the coordinates spanning the
      !! molecules that no bond provides -- but it pays that quadratic cost over
      !! the whole system rather than within a residue.

   integer, parameter :: OPT_ALGO_UNKNOWN = 0
   integer, parameter :: OPT_ALGO_SD = 1      !! Steepest descent
   integer, parameter :: OPT_ALGO_CG = 2      !! Conjugate gradient
   integer, parameter :: OPT_ALGO_LBFGS = 3   !! Limited-memory BFGS
   integer, parameter :: OPT_ALGO_PRFO = 4    !! Partitioned rational function, for a saddle
   integer, parameter :: OPT_ALGO_CG_AUTO = 5
      !! Conjugate gradient restarting on the Powell-Beale test rather than
      !! every ten steps. The other spelling of the same method: `cg` restarts
      !! on a fixed schedule, this one when the search directions stop being
      !! usefully conjugate.
   integer, parameter :: OPT_ALGO_NR = 6
      !! Newton-Raphson. Needs a Hessian, and descends to whatever stationary
      !! point it starts nearest -- which is a minimum only if it began in that
      !! basin.
   integer, parameter :: OPT_ALGO_DAMPED = 7
      !! Damped molecular dynamics. Not a minimiser in the line-search sense:
      !! it integrates motion and bleeds energy out through a friction term, so
      !! it crosses small barriers a downhill method would stop at.

   ! What the optimization is looking for. Not the same question as which
   ! algorithm runs: P-RFO can be asked for either, and a Newton-Raphson run
   ! settles on whatever stationary point it started nearest. Naming the target
   ! separately is what lets the end-of-run Hessian say whether the thing found
   ! is the thing wanted, rather than always reporting against "minimum".
   integer, parameter :: OPT_TARGET_MINIMUM = 0
   integer, parameter :: OPT_TARGET_SADDLE = 1
      !! A first-order saddle: one imaginary frequency, and exactly one. Zero
      !! means the search fell off the ridge into a basin; two or more means it
      !! is a higher-order stationary point and not a transition state.

   ! How a chain-of-states run treats the two structures it was given. The
   ! endpoints are usually already optimized minima, so freezing them is the
   ! default: relaxing a minimum again costs images and moves nothing.
   ! How a saddle is hunted, once `target` says one is wanted. Not the same
   ! axis as `algorithm`: the dimer is a way of *searching*, and the algorithm
   ! named alongside it is what translates and rotates the pair.
   integer, parameter :: SADDLE_METHOD_PRFO = 0
      !! Follow a Hessian eigenvector uphill. Needs curvature, and pays for it.
   integer, parameter :: SADDLE_METHOD_DIMER = 1
      !! Two images a short distance apart, rotated to find the softest
      !! direction and translated with the force along it inverted. Never forms
      !! a Hessian at all, which is the point: on a system where the analytic
      !! Hessian is the bottleneck, this is the method that runs.

   integer, parameter :: MIN_NEB_IMAGES = 3
      !! DL-FIND's own floor. Two images are the two endpoints, so a band of
      !! two has nothing between them to relax.
   integer, parameter :: NEB_ENDS_FROZEN = 0
   integer, parameter :: NEB_ENDS_PERPENDICULAR = 1
      !! Endpoints move only perpendicular to their own tangent, which lets a
      !! slightly-off endpoint settle onto the path without sliding along it.
   integer, parameter :: NEB_ENDS_FREE = 2
      !! Fully relaxed endpoints. For a pair of structures that are not yet
      !! minima, at the cost of the path being free to shorten itself.

   ! How the engine carries curvature between steps once it has a Hessian.
   ! Only the algorithms that hold one pay attention to this.
   integer, parameter :: OPT_HESSIAN_UPDATE_ENGINE = -1  !! Leave DL-FIND's own choice
   integer, parameter :: OPT_HESSIAN_UPDATE_NONE = 0
   integer, parameter :: OPT_HESSIAN_UPDATE_POWELL = 1
   integer, parameter :: OPT_HESSIAN_UPDATE_BOFILL = 2
      !! The usual choice for a saddle point: Powell's update keeps the matrix
      !! symmetric but not positive definite, which is what a transition state
      !! needs, and Bofill mixes it with BFGS by how much the step looks like
      !! either case.

   ! DL-FIND's constraint numbering, which the bridge writes into `spec`.
   integer, parameter :: OPT_CONSTRAIN_BOND = 1
   integer, parameter :: OPT_CONSTRAIN_ANGLE = 2
   integer, parameter :: OPT_CONSTRAIN_TORSION = 3
   integer, parameter :: OPT_CONSTRAIN_CARTESIAN = 4
   integer, parameter :: OPT_CONSTRAIN_BOND_DIFFERENCE = 5

   type :: opt_constraint_t
      !! One geometric quantity held fixed during the optimization
      !!
      !! `atoms` is 1-based here and 0-based in the deck, the same as every
      !! other atom index this program reads. Only the first
      !! `constraint_atom_count(kind)` entries mean anything.
      integer :: kind = 0
      integer :: atoms(4) = 0
   end type opt_constraint_t

   type :: optimizer_settings_t
      !! One optimization's settings, as the deck asked for them
      !!
      !! A negative real or integer means "leave the engine's own default
      !! alone", which is DL-FIND's own convention internally and is why it is
      !! used here rather than transcribing every default twice. The two
      !! settings a user actually reads back -- the step cap and the gradient
      !! threshold -- are given real defaults anyway, because a run that stops
      !! should be able to say what it stopped on.
      integer :: max_steps = DEFAULT_OPT_MAX_STEPS
         !! Give up after this many steps without converging
      real(dp) :: gradient_tolerance = DEFAULT_OPT_GRADIENT_TOLERANCE
         !! Convergence on the largest gradient component, Hartree/Bohr
      real(dp) :: energy_tolerance = DEFAULT_OPT_ENERGY_TOLERANCE
         !! Convergence on the energy change, Hartree. Negative: engine default
      real(dp) :: max_step = DEFAULT_OPT_MAX_STEP
         !! Longest step the engine may take, Bohr. Negative: engine default
      integer :: coordinates = OPT_COORDS_CARTESIAN
      integer :: algorithm = OPT_ALGO_LBFGS
      integer :: lbfgs_memory = DEFAULT_OPT_LBFGS_MEMORY
         !! Steps of curvature history. Negative: engine default
      integer :: print_level = DEFAULT_OPT_PRINT_LEVEL
         !! How much the engine itself prints. Negative: engine default
      logical :: freeze_terms = .true.
         !! Choose the MBE term list once, at the starting geometry, and keep
         !! it for every step. See `mqc_geometry_optimizer` for why an
         !! optimization on a re-screened list is optimizing a moving target.
      integer :: target = OPT_TARGET_MINIMUM
         !! Whether this run is looking for a minimum or a first-order saddle.
         !! Independent of `algorithm` on purpose -- see `OPT_TARGET_SADDLE`.
      character(len=:), allocatable :: endpoint
         !! Path to a second geometry: the product, for a chain-of-states run.
         !! Unallocated is the ordinary single-structure optimization, and is
         !! what makes NEB opt-in rather than inferred.
      integer :: n_images = 0
         !! Images along the path, endpoints included. Zero means the engine's
         !! own default. DL-FIND refuses fewer than three, since two images are
         !! just the endpoints and there is no path between them to relax.
      real(dp) :: neb_spring = -1.0_dp
         !! Spring constant holding neighbouring images apart. Negative means
         !! the engine default.
      integer :: neb_ends = NEB_ENDS_FROZEN
      integer :: saddle_method = SADDLE_METHOD_PRFO
      real(dp) :: dimer_separation = -1.0_dp
         !! Distance between the two dimer images, Bohr. Negative is the engine
         !! default. Small enough to sample curvature, large enough that the
         !! gradient difference is not noise.
      integer :: dimer_max_rotations = -1
         !! Rotation steps per translation step. Negative is the engine default.
      real(dp) :: dimer_rotation_tolerance = -1.0_dp
         !! Stop rotating once the angle falls below this. Negative is the
         !! engine default.
      integer :: zero_modes = -1
         !! How many Hessian eigenvalues to treat as zero rather than as
         !! vibrations. Negative is the engine default. Six for a non-linear
         !! molecule in Cartesian coordinates -- three translations and three
         !! rotations -- and five for a linear one. Internal coordinates have
         !! none, which is why this is off by default.
      real(dp) :: soft_mode_threshold = -1.0_dp
         !! Eigenvalues smaller than this in magnitude are soft: P-RFO steps
         !! along them are damped rather than trusted. Negative is the engine
         !! default.
      logical :: connect = .false.
         !! After the saddle is found, follow its imaginary mode downhill in
         !! both directions and report the two minima it falls to. One
         !! imaginary frequency proves a first-order saddle; it does not prove
         !! the saddle joins the two structures anybody had in mind, and this
         !! is the check that settles it.
      real(dp) :: connect_distort = -1.0_dp
         !! How far to displace along the imaginary mode before relaxing.
         !! Negative is the engine default. Too small and both sides fall back
         !! to the saddle; too large and they skip past the adjacent minima.
      integer :: hessian_update = OPT_HESSIAN_UPDATE_ENGINE
         !! How curvature is carried between steps once a Hessian exists.
         !! Ignored by the algorithms that never build one.
      real(dp) :: timestep = -1.0_dp
         !! Damped dynamics only. Negative: engine default
      real(dp) :: friction = -1.0_dp
         !! Damped dynamics only: the starting friction. Negative: engine default
      real(dp) :: friction_factor = -1.0_dp
         !! Damped dynamics only: what the friction is multiplied by while the
         !! energy is falling, so the run coasts further the better it is
         !! going. Negative: engine default
      real(dp) :: friction_rising = -1.0_dp
         !! Damped dynamics only: the friction used on a step where the energy
         !! went up. Negative: engine default
      integer, allocatable :: frozen_atoms(:)
         !! Atoms held at their input position, 1-based. Frozen atoms are
         !! removed from the optimization's variables rather than restrained,
         !! so a frozen pair contributes no coordinates at all.
      type(opt_constraint_t), allocatable :: constraints(:)
         !! Internal coordinates held fixed. Needs an internal coordinate
         !! system to hold them in -- see `mqc_dlfind_bridge`.
      logical :: hess_end = .false.
         !! Compute a Hessian at the converged geometry and say whether any
         !! frequency came back imaginary.
         !!
         !! A minimiser converges on the gradient alone, and a vanishing
         !! gradient is a stationary point rather than a minimum -- a saddle
         !! satisfies it exactly as well. Nothing in the optimization can tell
         !! the two apart, so a run that ends "converged" has not established
         !! what it converged *to*. This is what establishes it.
         !!
         !! Off by default because it is not free: one Hessian on top of the
         !! optimization, analytic for restricted Hartree-Fock and central
         !! differences of gradients otherwise.
      logical :: trajectory = .true.
         !! Record every accepted geometry, for analysing the path afterwards.
         !! Worth turning off for a large system optimized over many steps:
         !! the geometries are held until the run ends so they can be written
         !! as one document, which is 3N doubles per step of memory.
   end type optimizer_settings_t

   abstract interface
      subroutine hessian_i(n_atoms, coords, hessian, status)
         !! Second derivatives at a geometry, Hartree/Bohr^2
         !!
         !! Cartesian, `(3N, 3N)`, in the system's own atom order. Cartesian
         !! regardless of the coordinate system the optimization runs in: the
         !! engine transforms it to whatever internals it is using, and cannot
         !! do that if it is handed something already transformed.
         !!
         !! `status` is 0 for success. Anything else means no Hessian was
         !! produced, which is not a failure -- an engine that asked for one is
         !! expected to fall back to finite differences of gradients, which is
         !! slower and correct. That is why this interface has no `error_t`:
         !! there is nothing here a caller could not recover from.
         import :: dp
         implicit none
         integer, intent(in) :: n_atoms
         real(dp), intent(in) :: coords(3, n_atoms)
         real(dp), intent(out) :: hessian(3*n_atoms, 3*n_atoms)
         integer, intent(out) :: status
      end subroutine hessian_i
   end interface

   abstract interface
      subroutine energy_gradient_i(n_atoms, coords, energy, gradient, status)
         !! Energy and gradient at a geometry, in atomic units
         !!
         !! Coordinates are Bohr and the gradient is Hartree/Bohr, both in the
         !! system's own atom order. Fragmentation and hydrogen caps have
         !! already been folded back in by the time an engine sees this, so an
         !! engine never learns whether the number it was handed came from one
         !! calculation or from twenty thousand.
         !!
         !! `status` is 0 for success. An engine must stop on anything else
         !! rather than step on a gradient that was never computed -- a failed
         !! fragment leaves the array at whatever it held, which reads as a
         !! perfectly plausible small gradient.
         import :: dp
         implicit none
         integer, intent(in) :: n_atoms
         real(dp), intent(in) :: coords(3, n_atoms)
         real(dp), intent(out) :: energy
         real(dp), intent(out) :: gradient(3, n_atoms)
         integer, intent(out) :: status
      end subroutine energy_gradient_i
   end interface

   abstract interface
      subroutine step_callback_i(n_atoms, coords, energy)
         !! One accepted geometry, as the optimization passes through it
         !!
         !! Called by the engine for each step it accepts -- not for every
         !! energy evaluation, which includes the trial points of a line
         !! search and would make a trajectory that goes backwards. What the
         !! driver does with it (writes a frame, records it for the output
         !! document) is not the engine's business.
         import :: dp
         implicit none
         integer, intent(in) :: n_atoms
         real(dp), intent(in) :: coords(3, n_atoms)
         real(dp), intent(in) :: energy
      end subroutine step_callback_i
   end interface

contains

   pure function coordinates_from_string(text) result(coordinates)
      !! Parse the coordinate system a deck asked for
      !!
      !! Returns OPT_COORDS_UNKNOWN for anything unrecognised, which the reader
      !! turns into a refusal naming the spellings that do work.
      character(len=*), intent(in) :: text
      integer :: coordinates

      select case (lowercase(text))
      case ("cartesian", "cart", "xyz")
         coordinates = OPT_COORDS_CARTESIAN
      case ("hdlc")
         coordinates = OPT_COORDS_HDLC
      case ("dlc", "internal")
         coordinates = OPT_COORDS_DLC
      case ("hdlc-tc", "hdlc_tc", "hdlctc")
         coordinates = OPT_COORDS_HDLC_TC
      case ("dlc-tc", "dlc_tc", "dlctc")
         coordinates = OPT_COORDS_DLC_TC
      case default
         coordinates = OPT_COORDS_UNKNOWN
      end select
   end function coordinates_from_string

   pure function coordinates_to_string(coordinates) result(text)
      !! Name a coordinate system, for the log and the output document
      integer, intent(in) :: coordinates
      character(len=:), allocatable :: text

      select case (coordinates)
      case (OPT_COORDS_CARTESIAN)
         text = "cartesian"
      case (OPT_COORDS_HDLC)
         text = "hdlc"
      case (OPT_COORDS_DLC)
         text = "dlc"
      case (OPT_COORDS_HDLC_TC)
         text = "hdlc-tc"
      case (OPT_COORDS_DLC_TC)
         text = "dlc-tc"
      case default
         text = "unknown"
      end select
   end function coordinates_to_string

   pure function algorithm_from_string(text) result(algorithm)
      !! Parse the stepping algorithm a deck asked for
      character(len=*), intent(in) :: text
      integer :: algorithm

      select case (lowercase(text))
      case ("lbfgs", "l-bfgs")
         algorithm = OPT_ALGO_LBFGS
      case ("cg", "conjugate-gradient")
         algorithm = OPT_ALGO_CG
      case ("sd", "steepest-descent")
         algorithm = OPT_ALGO_SD
      case ("prfo", "p-rfo")
         algorithm = OPT_ALGO_PRFO
      case ("cg-auto", "conjugate-gradient-auto")
         algorithm = OPT_ALGO_CG_AUTO
      case ("nr", "newton-raphson")
         algorithm = OPT_ALGO_NR
      case ("damped", "damped-dynamics")
         algorithm = OPT_ALGO_DAMPED
      case default
         algorithm = OPT_ALGO_UNKNOWN
      end select
   end function algorithm_from_string

   pure function algorithm_to_string(algorithm) result(text)
      !! Name a stepping algorithm, for the log and the output document
      integer, intent(in) :: algorithm
      character(len=:), allocatable :: text

      select case (algorithm)
      case (OPT_ALGO_SD)
         text = "steepest-descent"
      case (OPT_ALGO_CG)
         text = "conjugate-gradient"
      case (OPT_ALGO_LBFGS)
         text = "lbfgs"
      case (OPT_ALGO_PRFO)
         text = "prfo"
      case (OPT_ALGO_CG_AUTO)
         text = "cg-auto"
      case (OPT_ALGO_NR)
         text = "newton-raphson"
      case (OPT_ALGO_DAMPED)
         text = "damped-dynamics"
      case default
         text = "unknown"
      end select
   end function algorithm_to_string

   pure function algorithm_needs_hessian(algorithm) result(needs)
      !! Whether this algorithm holds a Hessian at all
      !!
      !! What it decides is whether the driver bothers to offer one. An
      !! algorithm that never builds one would be handed a callback it never
      !! calls, and the Hessian settings beside it -- the update scheme above
      !! all -- would read as though they were doing something.
      integer, intent(in) :: algorithm
      logical :: needs

      needs = algorithm == OPT_ALGO_PRFO .or. algorithm == OPT_ALGO_NR
   end function algorithm_needs_hessian

   pure function saddle_method_from_string(text) result(method)
      !! Parse how a deck wants the saddle hunted
      character(len=*), intent(in) :: text
      integer :: method

      select case (lowercase(text))
      case ("prfo", "p-rfo", "hessian")
         method = SADDLE_METHOD_PRFO
      case ("dimer")
         method = SADDLE_METHOD_DIMER
      case default
         method = -2
      end select
   end function saddle_method_from_string

   pure function saddle_method_to_string(method) result(text)
      !! Name a saddle-search method, for the log
      integer, intent(in) :: method
      character(len=:), allocatable :: text

      select case (method)
      case (SADDLE_METHOD_PRFO)
         text = "p-rfo"
      case (SADDLE_METHOD_DIMER)
         text = "dimer"
      case default
         text = "unknown"
      end select
   end function saddle_method_to_string

   pure function neb_ends_from_string(text) result(ends)
      !! Parse how a deck wants the path endpoints treated
      character(len=*), intent(in) :: text
      integer :: ends

      select case (lowercase(text))
      case ("frozen", "fixed")
         ends = NEB_ENDS_FROZEN
      case ("perpendicular", "perp")
         ends = NEB_ENDS_PERPENDICULAR
      case ("free", "relaxed")
         ends = NEB_ENDS_FREE
      case default
         ends = -2
      end select
   end function neb_ends_from_string

   pure function neb_ends_to_string(ends) result(text)
      !! Name an endpoint treatment, for the log
      integer, intent(in) :: ends
      character(len=:), allocatable :: text

      select case (ends)
      case (NEB_ENDS_FROZEN)
         text = "frozen"
      case (NEB_ENDS_PERPENDICULAR)
         text = "perpendicular"
      case (NEB_ENDS_FREE)
         text = "free"
      case default
         text = "unknown"
      end select
   end function neb_ends_to_string

   pure function algorithm_finds_saddle(algorithm) result(can)
      !! Whether this algorithm can converge on a first-order saddle
      !!
      !! P-RFO by construction -- it maximises along one eigenvector and
      !! minimises along the others. Newton-Raphson because it seeks a
      !! vanishing gradient and does not care about the curvature it lands in,
      !! so from a guess near the ridge it will settle there. Everything else
      !! here goes downhill and will find a minimum however the deck spells the
      !! target.
      integer, intent(in) :: algorithm
      logical :: can

      can = algorithm == OPT_ALGO_PRFO .or. algorithm == OPT_ALGO_NR
   end function algorithm_finds_saddle

   pure function target_from_string(text) result(want)
      !! Parse what the deck says it is looking for
      !!
      !! The result is `want` and not `target`: the component these two read
      !! and write keeps the keyword's own spelling, but a local of that name
      !! shadows the Fortran attribute, which is legal and reads badly.
      character(len=*), intent(in) :: text
      integer :: want

      select case (lowercase(text))
      case ("minimum", "min")
         want = OPT_TARGET_MINIMUM
      case ("saddle", "ts", "transition_state", "transition-state")
         want = OPT_TARGET_SADDLE
      case default
         want = -2   ! not a spelling anybody recognises
      end select
   end function target_from_string

   pure function target_to_string(want) result(text)
      !! Name a target, for the log
      integer, intent(in) :: want
      character(len=:), allocatable :: text

      select case (want)
      case (OPT_TARGET_MINIMUM)
         text = "minimum"
      case (OPT_TARGET_SADDLE)
         text = "saddle"
      case default
         text = "unknown"
      end select
   end function target_to_string

   pure function hessian_update_from_string(text) result(update)
      !! Parse the Hessian update scheme a deck asked for
      character(len=*), intent(in) :: text
      integer :: update

      select case (lowercase(text))
      case ("none")
         update = OPT_HESSIAN_UPDATE_NONE
      case ("powell")
         update = OPT_HESSIAN_UPDATE_POWELL
      case ("bofill")
         update = OPT_HESSIAN_UPDATE_BOFILL
      case ("auto", "engine", "default")
         update = OPT_HESSIAN_UPDATE_ENGINE
      case default
         update = -2   ! not a spelling anybody recognises
      end select
   end function hessian_update_from_string

   pure function hessian_update_to_string(update) result(text)
      !! Name a Hessian update scheme, for the log
      integer, intent(in) :: update
      character(len=:), allocatable :: text

      select case (update)
      case (OPT_HESSIAN_UPDATE_NONE)
         text = "none"
      case (OPT_HESSIAN_UPDATE_POWELL)
         text = "powell"
      case (OPT_HESSIAN_UPDATE_BOFILL)
         text = "bofill"
      case (OPT_HESSIAN_UPDATE_ENGINE)
         text = "engine default"
      case default
         text = "unknown"
      end select
   end function hessian_update_to_string

   pure function constraint_from_string(text) result(kind)
      !! Parse a constrained coordinate's type
      character(len=*), intent(in) :: text
      integer :: kind

      select case (lowercase(text))
      case ("bond", "distance")
         kind = OPT_CONSTRAIN_BOND
      case ("angle")
         kind = OPT_CONSTRAIN_ANGLE
      case ("torsion", "dihedral")
         kind = OPT_CONSTRAIN_TORSION
      case ("cartesian", "position")
         kind = OPT_CONSTRAIN_CARTESIAN
      case ("bond-difference", "bond_difference")
         kind = OPT_CONSTRAIN_BOND_DIFFERENCE
      case default
         kind = 0
      end select
   end function constraint_from_string

   pure function constraint_atom_count(kind) result(n)
      !! How many atoms a constraint of this type names
      !!
      !! Zero for a type nobody recognises, which is how a deck naming one is
      !! caught: the count is what the reader validates the atom list against.
      integer, intent(in) :: kind
      integer :: n

      select case (kind)
      case (OPT_CONSTRAIN_BOND)
         n = 2
      case (OPT_CONSTRAIN_ANGLE)
         n = 3
      case (OPT_CONSTRAIN_TORSION)
         n = 4
      case (OPT_CONSTRAIN_CARTESIAN)
         n = 1
      case (OPT_CONSTRAIN_BOND_DIFFERENCE)
         n = 3
      case default
         n = 0
      end select
   end function constraint_atom_count

   pure function lowercase(text) result(lowered)
      !! Case-insensitive comparison, the way mqc_calc_types does it
      character(len=*), intent(in) :: text
      character(len=len_trim(adjustl(text))) :: lowered

      integer :: i

      lowered = trim(adjustl(text))
      do i = 1, len(lowered)
         if (lowered(i:i) >= "A" .and. lowered(i:i) <= "Z") then
            lowered(i:i) = achar(iachar(lowered(i:i)) + 32)
         end if
      end do
   end function lowercase

end module mqc_optimizer_types
