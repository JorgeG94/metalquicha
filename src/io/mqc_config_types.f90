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
      character(len=:), allocatable :: aux_basis
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
      integer :: dft_grid_level = 3
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
      if (allocated(this%log_level)) deallocate (this%log_level)
      if (allocated(this%fragment_breakdown)) deallocate (this%fragment_breakdown)
      if (allocated(this%frag_method)) deallocate (this%frag_method)
      if (allocated(this%embedding)) deallocate (this%embedding)
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
