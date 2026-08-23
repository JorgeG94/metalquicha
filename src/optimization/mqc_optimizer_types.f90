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
   public :: energy_gradient_i
   public :: step_callback_i
   public :: OPT_COORDS_UNKNOWN, OPT_COORDS_CARTESIAN, OPT_COORDS_HDLC, OPT_COORDS_DLC
   public :: OPT_COORDS_HDLC_TC, OPT_COORDS_DLC_TC
   public :: OPT_ALGO_UNKNOWN, OPT_ALGO_SD, OPT_ALGO_CG, OPT_ALGO_LBFGS, OPT_ALGO_PRFO
   public :: coordinates_from_string, coordinates_to_string
   public :: algorithm_from_string, algorithm_to_string

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
      case default
         text = "unknown"
      end select
   end function algorithm_to_string

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
