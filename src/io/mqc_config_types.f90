!! Configuration types an input file parses into
module mqc_config_types
   !! `mqc_config_t` and the pieces it is made of: one molecule, one input
   !! fragment, and the bonds between them.
   !!
   !! These live apart from any reader on purpose. Nearly everything that
   !! touches a config only wants the *type* -- the driver, the adapter, the
   !! MBE engine and the tests all import `mqc_config_t` and never parse
   !! anything -- so tying the type to one file format would make every
   !! consumer depend on that format's parser. It also means a second reader
   !! is a new module rather than a change to an existing one.
   !!
   !! Defaults are declared here, once, and a reader that finds no value for a
   !! key simply leaves the component alone. That is the whole mechanism: there
   !! is no separate table of defaults to drift out of step.
   use pic_types, only: dp, int32
   use mqc_method_types, only: METHOD_TYPE_GFN2
   use mqc_calc_types, only: CALC_TYPE_ENERGY
   use mqc_geometry, only: geometry_type
   use mqc_physical_fragment, only: bond_t
   use mqc_calculation_defaults, only: DEFAULT_DISPLACEMENT, DEFAULT_TEMPERATURE, &
                                       DEFAULT_PRESSURE, DEFAULT_AIMD_DT, &
                                       DEFAULT_AIMD_NSTEPS, DEFAULT_AIMD_TEMPERATURE, &
                                       DEFAULT_AIMD_OUTPUT_FREQ, DEFAULT_SCF_CONV, &
                                       DEFAULT_SCF_DENSITY_CONV, DEFAULT_VDW_SCALE, &
                                       DEFAULT_DYNAMIC_TOL, DEFAULT_DYNAMIC_MAXITER, &
                                       DEFAULT_RESPONSE_BATCH, &
                                       EFP_RESPONSE_AUTO, &
                                       DEFAULT_OPT_MAX_STEPS, DEFAULT_OPT_GRADIENT_TOLERANCE, &
                                       DEFAULT_OPT_ENERGY_TOLERANCE, DEFAULT_OPT_MAX_STEP, &
                                       DEFAULT_OPT_LBFGS_MEMORY, DEFAULT_OPT_PRINT_LEVEL, &
                                       DEFAULT_CPCM_NANG, DEFAULT_CPCM_RSCALE, &
                                       DEFAULT_PCM_NANG, DEFAULT_PCM_RSCALE, &
                                       DEFAULT_PCM_ZETA, &
                                       DEFAULT_FRAG_LEVEL, DEFAULT_MAX_INTERSECTION
   implicit none
   private

   public :: mqc_config_t       !! Everything one input file specifies
   public :: molecule_t         !! One molecule of a multi-molecule input
   public :: input_fragment_t   !! One fragment definition, as written in the input
   public :: bond_t             !! Re-exported so consumers need only this module
   public :: guess_step_t       !! One rung of a basis-set-projection ladder

   type :: input_fragment_t
      !! Input fragment definition with charge, multiplicity, and atom indices
      !! This is the parsed representation from the input file, not the computational fragment
      integer :: charge = 0
      integer :: multiplicity = 1
      integer, allocatable :: indices(:)  !! Atom indices in this fragment
      !> The effective fragment potential describing this fragment, if it has one.
      !> A fragment carrying one is an EFP fragment and a fragment without one stays
      !> quantum, which is how a mixed calculation is spelled -- there is no second way
      !> to declare a fragment. The path is resolved relative to the deck, like `xyz`.
      !>
      !> **A path, not a format.** GAMESS's `.efp` text is what exists today, but the
      !> potential is moving to JSON, so nothing here may require an extension or infer
      !> a format from one. Whoever loads it dispatches on what the file actually is;
      !> the deck only ever names it.
      character(len=:), allocatable :: potential
   contains
      procedure :: destroy => input_fragment_destroy
   end type input_fragment_t

   type :: molecule_t
      !! Single molecule definition with structure, geometry, fragments, and connectivity
      character(len=:), allocatable :: name  !! Optional molecule name

      ! Structure information
      integer :: charge = 0
      integer :: multiplicity = 1

      ! Geometry
      type(geometry_type) :: geometry

      ! Fragments
      integer :: nfrag = 0
      type(input_fragment_t), allocatable :: fragments(:)
      !> Every fragment is the same species, so one potential describes all of them.
      !> Set by `uniform_system` in the deck, which is what lets `fragment_potentials`
      !> be a single name rather than ten thousand copies of it. It is an assertion the
      !> reader checks rather than takes on trust: the fragments must agree in size, so
      !> a deck claiming uniformity without it is refused rather than quietly handed
      !> the wrong potential.
      logical :: uniform_system = .false.

      ! Connectivity
      integer :: nbonds = 0
      integer :: nbroken = 0
      type(bond_t), allocatable :: bonds(:)

   contains
      procedure :: destroy => molecule_destroy
   end type molecule_t

   type :: guess_step_t
      !! One rung of a basis-set-projection ladder
      !!
      !! Each rung converges an SCF in its own basis and hands the result to the
      !! next as a starting density. The convergence settings are per rung
      !! because the rungs are not doing the same job: an early one only has to
      !! land in the right basin and can stop loosely, where the last one before
      !! the target may as well be tight since it is cheap relative to what
      !! follows.
      character(len=:), allocatable :: basis
      integer :: maxiter = 50
      real(dp) :: tolerance = 1.0e-5_dp
   end type guess_step_t

   type :: mqc_config_t
      !! Complete configuration from .mqc file

      ! Schema information
      character(len=:), allocatable :: schema_name
      character(len=:), allocatable :: schema_version
      integer :: index_base = 0  !! 0-based or 1-based indexing
      character(len=:), allocatable :: units  !! angstrom or bohr

      ! Properties asked for alongside the energy. Not settings for *how* to
      ! compute -- those are `keywords` -- but requests for something extra to
      ! be reported once the wave function exists.
      character(len=:), allocatable :: bonding_analysis
         !! From `properties.bonding_analysis.type`. Absent or "none" means no
         !! analysis; "gms_quao" means the quasi-atomic bonding picture.
      character(len=:), allocatable :: fukui_population
         !! From `properties.fukui.population`. Allocated means the deck asked
         !! for Fukui indices; "chelpg" or "mulliken" chooses how the density
         !! difference is condensed onto atoms. Required inside the object
         !! rather than defaulted, so that asking for the analysis and choosing
         !! how it is condensed are one decision -- the choice changes the
         !! numbers enough that it should not be made by omission.
      character(len=:), allocatable :: charges_scheme
         !! From `properties.charges.scheme`. Allocated means the deck asked for
         !! atomic partial charges; "mulliken" or "chelpg" chooses how the
         !! converged density is partitioned. Unlike `fukui_population` this
         !! carries a default, because here the scheme is the whole request
         !! rather than one half of it -- see `charges_keys` in the schema.
      logical :: bonding_energy = .false.
         !! From `properties.bonding_analysis.energy_decomposition`. Off by
         !! default because the two-electron term needs the dense `n_ao^4`
         !! integral array, and the bonding tables on their own need only
         !! one-electron integrals. Asking for the decomposition is asking for
         !! that memory.
      logical :: bonding_no_sharing = .false.
         !! From `properties.bonding_analysis.no_sharing`. Asks for the
         !! no-sharing wave function, which needs a full valence CI over the
         !! quasi-atomic orbitals and is therefore far too expensive to do by
         !! default.
      logical :: bonding_restrict_localization = .false.
         !! From `properties.bonding_analysis.restrict_localization`. Confine
         !! the quasi-atomic localization to the wave function's
         !! occupation-restricted subspaces, so that no rotation mixes two of
         !! them and the wave function stays invariant under the
         !! transformation. Off by default, and that is a measurement rather
         !! than a preference -- see `mqc_docs/source/bonding_analysis.rst`.
      character(len=:), allocatable :: bonding_no_sharing_ci
         !! From `properties.bonding_analysis.no_sharing_ci`. How the CI
         !! expansion over quasi-atomic orbitals is obtained: `"transform"`
         !! solves in the molecular orbital basis and carries the vector across
         !! with the orbital transformation, `"resolve"` runs a second Davidson
         !! in the quasi-atomic basis. The two give the same wave function and
         !! this chooses how it is computed, not which one it is.
      real(dp) :: bonding_threshold = 1.0_dp
         !! From `properties.bonding_analysis.energy_threshold`, in kcal/mol.
         !! Orbital pairs whose kinetic bond order is weaker than this are
         !! counted and then not printed. Formyl chloride has 27 pairs above
         !! GAMESS's own cutoff and eight above one kcal/mol, and the other
         !! nineteen are tenths of a kcal/mol that bury the four bonds.

      ! Model information
      integer(int32) :: method = METHOD_TYPE_GFN2
      character(len=:), allocatable :: basis
      character(len=:), allocatable :: ecp
         !! `model.ecp`. Unallocated means no potential, which is not the
         !! same as an empty one: a deck that names a set gets it applied to
         !! whichever atoms the file covers, and a deck that names none gets
         !! an all-electron calculation.
      character(len=:), allocatable :: aux_basis
      logical :: cartesian = .false.
         !! Read the basis in Cartesian form whatever its file declares.
         !!
         !! Not a preference: for a basis with a shell above p the two forms
         !! span different spaces and give different energies, so this changes
         !! the model rather than its representation. It exists because BSE is
         !! inconsistent about the Pople sets -- 6-31G* is marked Cartesian and
         !! cc-pVDZ spherical -- and because GAMESS assumes Cartesian for the
         !! Pople sets throughout, so reproducing a GAMESS number needs a way
         !! to say so. Defaults false, which honours the file.
      character(len=:), allocatable :: functional
      character(len=:), allocatable :: scf_guess
         !! Initial guess name from `keywords.scf.guess`.
         !!
         !! Superseded by `keywords.guess.type`. Still read so existing decks
         !! keep working, and setting both is an error rather than a precedence
         !! rule: two keys naming one thing is confusing enough without the
         !! answer depending on which one the reader happened to reach first.
      character(len=:), allocatable :: guess_type
         !! Initial guess name from `keywords.guess.type`, the current spelling.
      type(guess_step_t), allocatable :: guess_steps(:)
         !! The projection ladder, from `keywords.guess.subscf.steps`.
         !!
         !! One entry per preliminary SCF, in order. The target basis is
         !! `model.basis` and is not repeated here, so three SCFs -- STO-3G, then
         !! 6-31G, then cc-pVTZ -- is two steps and a model basis.
      logical :: scf_unrestricted = .false.
      logical :: allow_crap_scf = .false.
      logical :: scf_density_fitting = .false.
      ! keywords.correlation -- post-Hartree-Fock, kept apart from the SCF ones
      logical :: corr_freeze_core = .true.
      integer :: corr_n_frozen_core = -1   !! -1 counts the core from the elements
      logical :: corr_density_fitting = .false.  !! RI for the correlation step
      logical :: corr_scs = .false.        !! Spin-component scaling
      real(dp) :: corr_scs_ss = 1.0_dp/3.0_dp
      real(dp) :: corr_scs_os = 1.2_dp
      ! keywords.cc -- what only an iterative correlation method needs
      integer :: cc_maxiter = 100
      real(dp) :: cc_tolerance = 1.0e-8_dp
      logical :: cc_triples_set = .false.
         !! Whether a deck named `triples` explicitly. Without this there is no
         !! way to tell "triples: false" from "the key was absent", and the
         !! method name has to lose to an explicit keyword while winning over an
         !! absent one.
      logical :: cc_triples = .false.
      logical :: cc_diis = .true.
      integer :: cc_diis_size = 8
      logical :: cc_spin_adapted = .true.
         !! Which coupled-cluster formulation runs: spin-adapted (spatial
         !! orbitals) when true, spin orbitals when false.
         !!
         !! Defaults to spin-adapted because for the closed-shell reference
         !! that is the only one either formulation supports, it is strictly
         !! better on both axes -- about sixteen times less memory and several
         !! times faster. The spin-orbital path stays reachable, and is what the
         !! two are cross-checked against: they must agree to machine precision,
         !! so a deck can settle any doubt about a number by running it twice.
      ! keywords.mcscf -- the active space, for CASSCF and CASCI
      character(len=:), allocatable :: mcscf_avas_orbitals(:)
         !! Atomic orbital labels from `keywords.mcscf.avas.orbitals`, e.g.
         !! "N 2p". Unallocated means the active space was given by counts.
      logical :: mcscf_full_valence = .false.
         !! `keywords.mcscf.full_valence` -- take the whole valence shell as the
         !! active space, sized by counting the free-atom minimal basis rather
         !! than by any threshold.
      real(dp) :: mcscf_avas_threshold = 0.2_dp
      integer, allocatable :: mcscf_ormas_subspaces(:)
         !! Active orbital each subspace starts at, from
         !! `keywords.mcscf.ormas.subspaces`. Absent leaves a complete active
         !! space, which is the one subspace covering everything.
      integer, allocatable :: mcscf_ormas_min_electrons(:)
      integer, allocatable :: mcscf_ormas_max_electrons(:)
         !! Electrons allowed in each subspace, both spins together
      integer :: mcscf_n_active_electrons = 0
      integer :: mcscf_n_active_orbitals = 0
      integer :: mcscf_n_inactive_orbitals = -1
         !! -1 derives it from the electron count: every electron the active
         !! space does not hold is in a doubly occupied orbital.
      logical :: mcscf_optimize_orbitals = .true.
         !! CASSCF when true, CASCI when false. "casci" and "casscf" parse to
         !! the same method type, so the spelling is the only place the
         !! distinction survives; the reader defaults this from it and lets
         !! `keywords.mcscf.optimize_orbitals` override.
      integer :: mcscf_max_macro_iter = 100
      real(dp) :: mcscf_orbital_convergence = 1.0e-6_dp
      ! keywords.dft -- the quadrature, not the functional
      real(dp) :: dft_screening_tolerance = 1.0e-12_dp
         !! From `keywords.dft.screening_tolerance`. The AO value below which a
         !! shell is dropped from a grid block. Zero or negative evaluates the
         !! whole basis, which is how a run measures what the screen is worth.
      integer :: dft_block_size = -1
         !! From `keywords.dft.block_size`. Grid points per block; -1 keeps the
         !! backend default. Sets both how tightly the screen bites and how many
         !! blocks there are to spread over threads.
      integer :: dft_grid_level = 3
      integer :: dft_nlc_grid_level = -1
         !! `keywords.dft.nlc_grid_level`, the quadrature VV10's double integral
         !! uses. Separate from `grid_level` because the non-local term costs
         !! the *product* of two grid sizes while everything else costs one, so
         !! the level that is right for exchange is an order of magnitude too
         !! expensive here. Negative means the backend's own default.
      integer :: dft_radial_points = -1    !! -1 leaves grid_level in charge
      integer :: dft_angular_points = -1
         !! Fit J and K against the auxiliary basis, on the CPU backend
         !! Let a non-converged SCF into the expansion instead of stopping.
         !! Off by default: the energy of an SCF that ran out of iterations is
         !! of the right magnitude and nothing downstream can tell, so silence
         !! is the one thing it must not be. On, every offender is named at the
         !! end of the run and marked in the fragment table.
         !! Force UHF/UKS even when the shell is closed
         !! XC functional name, only meaningful when method = dft

      ! keywords.pcm -- the polarizable continuum on the ab initio backends.
      !
      ! Deliberately not the xTB solvation keys below. Those configure tblite's
      ! own CPCM, which builds its own cavity with its own defaults
      ! (cpcm_rscale is 1.0 there, against the 1.2 conventional for a van der
      ! Waals cavity), and routing one implementation's settings into another's
      ! would make two different models look like one keyword.
      ! Which integral backend to run on: "auto" (default), "cuest"/"gpu", or
      ! "libcint"/"cpu". A request that the build or the method cannot honour is
      ! refused, not substituted.
      character(len=16) :: backend = "auto"
      ! system.gpu -- the same choice said the way people ask for it. It resolves
      ! into `backend` above rather than being carried alongside it, so there is
      ! one answer downstream however the deck spelled the question; naming both
      ! and disagreeing is refused rather than given a precedence rule.
      logical :: gpu = .false.
      logical :: gpu_set = .false.
         !! Whether the deck named `system.gpu`. "Absent" and "false" have to be
         !! told apart: absent leaves `backend` alone, false pins it to the CPU.

      logical :: pcm_enabled = .false.
      character(len=16) :: pcm_method = "cpcm"
         !! Continuum model: "cpcm" or "iefpcm". See `pcm_config_t`.
      real(dp) :: pcm_dielectric = -1.0_dp
         !! Solvent dielectric constant. Required: there is no solvent-name
         !! table on this path, because inventing one that disagreed with
         !! tblite's would make the same deck mean two things.
      integer :: pcm_angular_points = DEFAULT_PCM_NANG
      real(dp) :: pcm_radii_scale = DEFAULT_PCM_RSCALE
      real(dp) :: pcm_zeta = DEFAULT_PCM_ZETA
      real(dp) :: pcm_tolerance = 1.0e-8_dp
      integer :: pcm_max_iter = 100

      ! XTB solvation settings
      character(len=:), allocatable :: solvent  !! Solvent name (e.g., "water", "ethanol") or empty for gas phase
      character(len=:), allocatable :: solvation_model  !! Solvation model: "alpb" (default), "gbsa", or "cpcm"
      logical :: use_cds = .true.               !! Include non-polar CDS terms in solvation (not for CPCM)
      logical :: use_shift = .true.             !! Include solution state shift in solvation (not for CPCM)
      ! CPCM-specific settings
      real(dp) :: dielectric = -1.0_dp              !! Direct dielectric constant (-1 = use solvent lookup)
      integer :: cpcm_nang = DEFAULT_CPCM_NANG      !! Number of angular grid points for CPCM cavity
      real(dp) :: cpcm_rscale = DEFAULT_CPCM_RSCALE  !! Radii scaling factor for CPCM cavity

      ! Driver information
      integer(int32) :: calc_type = CALC_TYPE_ENERGY

      ! Multiple molecules support
      integer :: nmol = 0  !! Number of molecules (0 = single molecule mode for backward compatibility)
      type(molecule_t), allocatable :: molecules(:)  !! Array of molecules (if nmol > 0)

      ! Single molecule fields (backward compatibility - used if nmol == 0)
      ! Structure information
      integer :: charge = 0
      integer :: multiplicity = 1

      ! Geometry
      type(geometry_type) :: geometry

      ! Fragments
      integer :: nfrag = 0
      type(input_fragment_t), allocatable :: fragments(:)
      !> As on `molecule_t`: every fragment the same species, which is what lets
      !> `fragment_potentials` be one name rather than one per fragment.
      logical :: uniform_system = .false.

      ! Connectivity
      integer :: nbonds = 0
      integer :: nbroken = 0
      type(bond_t), allocatable :: bonds(:)

      ! SCF settings
      integer :: scf_maxiter = 300              !! Using 300 (parser-specific, different from DEFAULT_SCF_MAXITER)
      real(dp) :: scf_tolerance = DEFAULT_SCF_CONV
      logical :: scf_tolerance_set = .false.
         !! Whether the deck named it. Needed because "the default" and "the
         !! user asked for the default" are different requests to a caller that
         !! has a stricter default of its own -- MAKEFP converges to 1e-10/1e-8
         !! because the multipoles are taken from that density, and must keep
         !! doing so unless a deck says otherwise.
      real(dp) :: scf_density_tolerance = DEFAULT_SCF_DENSITY_CONV
      real(dp) :: scf_linear_dependence = 0.0_dp
         !! `keywords.scf.linear_dependence_threshold`. Overlap eigenvalues at
         !! or below this are dropped as linearly dependent. Raise it to shed
         !! more of a diffuse basis, lower it to keep modes the default would
         !! discard.
         !!
         !! Zero means "not set", and the orthogonaliser then uses
         !! `LINEAR_DEPENDENCE_TOL`. A sentinel rather than spelling 1e-7
         !! again here: this module would otherwise carry a second copy of a
         !! constant that already exists in `mqc_scf_common`, which is the
         !! arrangement that let the two backends drift apart the last time
         !! -- see the note on `LINEAR_DEPENDENCE_TOL` itself. Zero is not a
         !! meaningful cutoff on its own, so nothing is lost by spending it.
      logical :: scf_density_tolerance_set = .false.
      real(dp) :: scf_level_shift = 0.0_dp
      logical :: scf_use_diis = .true.
         !! `keywords.scf.diis`. Off is a diagnostic rather than a setting: an
         !! SCF without DIIS is how you find out whether DIIS is the thing
         !! hiding a problem, not a way to run one.
      integer :: scf_diis_size = 8
         !! `keywords.scf.diis_size`, the subspace the extrapolation is drawn
         !! from. Eight is the usual default and was, until this keyword, the
         !! only value obtainable -- the field existed all the way down to
         !! `scf_config_t` and nothing read it from a deck. Widening it to
         !! 12-20 is the first thing to try on an open-shell SCF that converges
         !! monotonically but slowly, where a level shift would only make it
         !! slower.
         !! Hartree added to the virtual orbitals before each diagonalisation.
         !! Zero is off, which is the default: a shift costs iterations near the
         !! solution and only earns them far from it.

      ! EFP settings, read from `keywords.efp` and used only by the MAKEFP driver
      !
      ! Apart from `scf` on purpose. The SCF a potential runs is configured by
      ! `keywords.scf` and reads those keys already; adding a second spelling of a
      ! tolerance here would let one deck set the same thing twice and then have to
      ! decide which copy wins. These four are the stages after the SCF, which is
      ! where the wall clock of a MAKEFP run actually goes.
      real(dp) :: efp_dynamic_tolerance = DEFAULT_DYNAMIC_TOL
      integer :: efp_dynamic_maxiter = DEFAULT_DYNAMIC_MAXITER
      logical :: efp_allow_crap_response = .false.
      integer :: efp_response_batch = DEFAULT_RESPONSE_BATCH
      integer :: efp_response = EFP_RESPONSE_AUTO
         !! Held as a code rather than as the spelling, the same as `calc_type`:
         !! the reader is the only place that has to know the three words, and it
         !! refuses a fourth there rather than somewhere an SCF has already run.
      real(dp) :: efp_vdw_scale = DEFAULT_VDW_SCALE

      ! Hessian settings
      real(dp) :: hessian_displacement = DEFAULT_DISPLACEMENT  !! Finite difference displacement (Bohr)
      real(dp) :: hessian_temperature = DEFAULT_TEMPERATURE    !! Temperature for thermochemistry (K)
      real(dp) :: hessian_pressure = DEFAULT_PRESSURE          !! Pressure for thermochemistry (atm)

      ! AIMD settings
      real(dp) :: aimd_dt = DEFAULT_AIMD_DT                          !! Timestep (femtoseconds)
      integer :: aimd_nsteps = DEFAULT_AIMD_NSTEPS                   !! Number of MD steps (0 = no AIMD)
      real(dp) :: aimd_initial_temperature = DEFAULT_AIMD_TEMPERATURE  !! Initial temperature for velocity init (K)
      integer :: aimd_output_frequency = DEFAULT_AIMD_OUTPUT_FREQ    !! Write output every N steps

      ! Geometry optimization settings
      !
      ! The two spelled ones stay strings here and are parsed in the adapter,
      ! where a misspelling can be refused by name. Parsing them in the reader
      ! instead would put the vocabulary in two places.
      integer :: opt_max_steps = DEFAULT_OPT_MAX_STEPS
      real(dp) :: opt_gradient_tolerance = DEFAULT_OPT_GRADIENT_TOLERANCE  !! Hartree/Bohr
      real(dp) :: opt_energy_tolerance = DEFAULT_OPT_ENERGY_TOLERANCE      !! Hartree; <0 = engine default
      real(dp) :: opt_max_step = DEFAULT_OPT_MAX_STEP                      !! Bohr; <0 = engine default
      integer :: opt_lbfgs_memory = DEFAULT_OPT_LBFGS_MEMORY               !! <0 = engine default
      integer :: opt_print_level = DEFAULT_OPT_PRINT_LEVEL                 !! <0 = engine default
      character(len=:), allocatable :: opt_coordinates  !! cartesian, hdlc, dlc
      character(len=:), allocatable :: opt_algorithm    !! lbfgs, cg, sd, prfo
      logical :: opt_trajectory = .true.               !! Record the path taken
      logical :: opt_freeze_terms = .true.             !! Fix the MBE term list for the run
      logical :: opt_hess_end = .false.                !! Hessian at the converged geometry
      character(len=:), allocatable :: opt_hessian_update  !! none, powell, bofill, auto
      character(len=:), allocatable :: opt_target  !! minimum or saddle
      character(len=:), allocatable :: opt_endpoint  !! xyz of the product, for NEB
      character(len=:), allocatable :: opt_neb_ends  !! frozen, perpendicular, free
      integer :: opt_n_images = 0                      !! Images along the path
      real(dp) :: opt_neb_spring = -1.0_dp             !! <0 = engine default
      character(len=:), allocatable :: opt_saddle_method  !! prfo or dimer
      real(dp) :: opt_dimer_separation = -1.0_dp       !! <0 = engine default
      integer :: opt_dimer_max_rotations = -1          !! <0 = engine default
      real(dp) :: opt_dimer_rot_tol = -1.0_dp          !! <0 = engine default
      integer :: opt_zero_modes = -1                   !! <0 = engine default
      real(dp) :: opt_soft_modes = -1.0_dp             !! <0 = engine default
      logical :: opt_connect = .false.                 !! Relax downhill from the saddle
      real(dp) :: opt_connect_distort = -1.0_dp        !! <0 = engine default
      real(dp) :: opt_timestep = -1.0_dp               !! Damped dynamics step
      real(dp) :: opt_friction = -1.0_dp               !! Damped dynamics start friction
      real(dp) :: opt_friction_factor = -1.0_dp        !! Friction decay while descending
      real(dp) :: opt_friction_rising = -1.0_dp        !! Friction on an uphill step
      integer, allocatable :: opt_frozen_atoms(:)      !! 0-based in the deck, 1-based here
      integer, allocatable :: opt_constraint_kinds(:)  !! One per constraint
      integer, allocatable :: opt_constraint_atoms(:, :)
         !! (4, n_constraints), 1-based, unused slots zero

      ! Fragmentation settings
      character(len=:), allocatable :: frag_method  !! MBE, etc.
      integer :: frag_level = DEFAULT_FRAG_LEVEL
      logical :: allow_overlapping_fragments = .false.
      character(len=:), allocatable :: expansion_kind
         !! Which expansion runs: "mbe" (default), "fmo" or "ee-mbe".
         !! GMBE is selected by `allow_overlapping_fragments` instead, which
         !! predates this and stays as it was.
      character(len=:), allocatable :: counterpoise
         !! How basis-set superposition error is handled: "none" (default) or
         !! "vmfc".
         !!
         !! Under "none" each subfragment is solved in its own basis, so a pair
         !! borrows functions its monomers did not have and looks more bound
         !! than it is. That error lands in the pair term and survives
         !! truncation. "vmfc" solves every subfragment in its parent's basis
         !! instead -- valiron-mayer function counterpoise -- so the borrowing
         !! appears on both sides of each difference and cancels.
         !!
         !! It is not free: a monomer's energy stops being one number reusable
         !! across every pair it belongs to and becomes one number per pair, so
         !! a level-2 expansion goes from N + C(N,2) subcalculations to
         !! 3*C(N,2).

      character(len=:), allocatable :: fmo_far_field
         !! What a distant fragment contributes: mulliken, chelpg or ignore
      real(dp) :: fmo_resppc = 2.0_dp
         !! Separation past which a fragment becomes point charges; negative
         !! turns the approximation off and makes every neighbour exact
      integer :: fmo_max_outer = 50
      real(dp) :: fmo_tolerance = 1.0e-7_dp
         !! Outer (monomer) SCF convergence, on the monomer energy sum
      integer :: fmo_scf_max_iter = 100
      real(dp) :: fmo_scf_energy_tol = 1.0e-9_dp
      real(dp) :: fmo_scf_density_tol = 1.0e-7_dp
         !! The inner per-fragment SCF: iteration cap and the energy/density
         !! convergence each monomer and n-mer SCF is held to. Independent of the
         !! outer loop above, and of a top-level `keywords.scf`, so a fragment run
         !! can be converged more loosely than a whole-system one would be.
      integer :: max_intersection_level = DEFAULT_MAX_INTERSECTION  !! Maximum k-way intersection depth for GMBE
      character(len=:), allocatable :: bond_breaking
         !! How a fragment represents a covalent bond the partition cut.
         !!
         !! `"none"` refuses a partition that cuts one, which is the only honest
         !! answer for a method with no way to represent the dangling valence.
         !! `"caps"` closes it with a hydrogen. `"afo"` is the adaptive frozen
         !! orbital, and is not built yet.
         !!
         !! Separate from `embedding` on purpose: what closes a cut bond and what
         !! field a fragment sits in are independent choices, and only some of
         !! their pairings are sound -- a cap puts an electron where an exact
         !! embedding already supplies the neighbour's density, so the two
         !! double-count and the combination is refused rather than approximated.
      real(dp) :: cap_scale = 1.0_dp
         !! Where the cap goes along the cut bond: `R_H = R_kept + s (R_gone - R_kept)`.
         !!
         !! `1.0` puts it exactly where the removed atom was, which is what this
         !! program has always done, and is why redistributing the cap's gradient
         !! onto that atom is the whole chain rule rather than an approximation.
         !! A physical C-H wants roughly `0.71` of a C-C, and the gradient then
         !! splits `s` and `1 - s` between the two atoms -- still exact, one term
         !! longer. The default keeps existing numbers reproducible; it is not a
         !! recommendation.
      character(len=:), allocatable :: embedding
      character(len=:), allocatable :: cutoff_method
      character(len=:), allocatable :: distance_metric
      real(dp), allocatable :: fragment_cutoffs(:)  !! Distance cutoffs indexed by n-mer level (2=dimer, 3=trimer, etc.)
      integer :: global_groups = 0
      integer :: nodes_per_group = 0

      ! Logger settings (kept for compatibility)
      character(len=:), allocatable :: log_level

      ! Output control
      logical :: skip_json_output = .false.  !! Skip JSON output for large calculations
      character(len=:), allocatable :: checkpoint_file
         !! Where to append each fragment energy as it finishes, and where to
         !! resume from. Absent means neither. An existing file whose
         !! fingerprint disagrees with this run is refused, not overwritten.
      logical :: unchecked_input = .false.
         !! Downgrade semantic validation to warnings.
         !!
         !! Off by default and meant to stay that way. The checks it disables
         !! are for inputs that produce a converged, plausible, wrong answer --
         !! monomer charges that do not sum, atoms in no fragment, a term list
         !! the expansion cannot be assembled from. Some are conventions rather
         !! than laws and a user may be breaking one deliberately, which is why
         !! the escape exists; it is deliberately not a comfortable thing to
         !! reach for.
      character(len=:), allocatable :: fragment_breakdown
         !! Where the per-fragment MBE table goes: "csv", "json" or "none"

   contains
      procedure :: destroy => config_destroy
   end type mqc_config_t

contains

   subroutine molecule_destroy(this)
      !! Clean up allocated memory in molecule_t
      class(molecule_t), intent(inout) :: this
      integer :: i

      if (allocated(this%name)) deallocate (this%name)

      call this%geometry%destroy()

      if (allocated(this%fragments)) then
         do i = 1, size(this%fragments)
            call this%fragments(i)%destroy()
         end do
         deallocate (this%fragments)
      end if

      if (allocated(this%bonds)) deallocate (this%bonds)

   end subroutine molecule_destroy

   subroutine config_destroy(this)
      !! Clean up allocated memory in mqc_config_t
      class(mqc_config_t), intent(inout) :: this
      integer :: i

      if (allocated(this%schema_name)) deallocate (this%schema_name)
      if (allocated(this%schema_version)) deallocate (this%schema_version)
      if (allocated(this%units)) deallocate (this%units)
      if (allocated(this%basis)) deallocate (this%basis)
      if (allocated(this%aux_basis)) deallocate (this%aux_basis)
      if (allocated(this%functional)) deallocate (this%functional)
      if (allocated(this%charges_scheme)) deallocate (this%charges_scheme)
      if (allocated(this%log_level)) deallocate (this%log_level)
      if (allocated(this%fragment_breakdown)) deallocate (this%fragment_breakdown)
      if (allocated(this%frag_method)) deallocate (this%frag_method)
      if (allocated(this%embedding)) deallocate (this%embedding)
      if (allocated(this%bond_breaking)) deallocate (this%bond_breaking)
      if (allocated(this%cutoff_method)) deallocate (this%cutoff_method)
      if (allocated(this%distance_metric)) deallocate (this%distance_metric)
      if (allocated(this%fragment_cutoffs)) deallocate (this%fragment_cutoffs)

      call this%geometry%destroy()

      if (allocated(this%fragments)) then
         do i = 1, size(this%fragments)
            call this%fragments(i)%destroy()
         end do
         deallocate (this%fragments)
      end if

      if (allocated(this%bonds)) deallocate (this%bonds)

      ! Clean up molecules array (multi-molecule mode)
      if (allocated(this%molecules)) then
         do i = 1, size(this%molecules)
            call this%molecules(i)%destroy()
         end do
         deallocate (this%molecules)
      end if

   end subroutine config_destroy

   subroutine input_fragment_destroy(this)
      !! Clean up allocated memory in input_fragment_t
      class(input_fragment_t), intent(inout) :: this
      if (allocated(this%indices)) deallocate (this%indices)
   end subroutine input_fragment_destroy
end module mqc_config_types
