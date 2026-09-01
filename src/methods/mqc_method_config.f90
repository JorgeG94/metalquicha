!! Unified configuration type for quantum chemistry method creation
module mqc_method_config
   !! Provides configuration types for all quantum chemistry methods.
   !! Uses composition pattern: method_config_t contains nested config types
   !! for each method family. The factory reads from the appropriate nested type.
   use pic_types, only: int32, dp
   use mqc_program_limits, only: MAX_ORBITAL_LABEL_LEN
   use mqc_scf_types, only: guess_step_t, scf_numerics_t, deltascf_options_t
   use mqc_method_types, only: METHOD_TYPE_UNKNOWN
   use mqc_calculation_defaults, only: DEFAULT_DISPLACEMENT, DEFAULT_VDW_SCALE, DEFAULT_DYNAMIC_TOL, &
                                       DEFAULT_DYNAMIC_MAXITER, EFP_RESPONSE_AUTO, &
                                       DEFAULT_RESPONSE_BATCH
   implicit none
   private

   public :: method_config_t
   public :: scf_config_t, xtb_config_t, dft_config_t, mcscf_config_t
   public :: scf_options_t
   public :: scf_numerics_t, deltascf_options_t  !! Re-exported from mqc_config_types
   public :: correlation_config_t, cc_config_t, f12_config_t
   public :: efp_config_t
   public :: pcm_config_t
   public :: properties_config_t

   !============================================================================
   ! SCF Configuration (shared by HF and DFT)
   !============================================================================
   type :: scf_config_t
      !! Shared SCF settings for HF and DFT methods
      logical :: cartesian = .false.
         !! `model.cartesian`; see `mqc_config_t`. Carried here rather than
         !! beside the basis name because it reaches the backend by the same
         !! route the rest of the SCF settings do.
      integer :: max_iter = 100
         !! Maximum SCF iterations
      real(dp) :: energy_convergence = 1.0e-8_dp
         !! Energy convergence threshold (Hartree)
      real(dp) :: density_convergence = 1.0e-6_dp
         !! Density matrix convergence threshold
      real(dp) :: gradient_convergence = 0.0_dp
         !! `keywords.scf.gradient_tolerance` -- pyscf's `conv_tol_grad`.
         !!
         !! Zero means "not set", and the SCF then derives `sqrt(tolerance)`,
         !! which is pyscf's own default. Set it when what you want from the
         !! run is the *density* rather than the energy: the energy's error
         !! goes as the commutator squared, the density's only as the
         !! commutator, so the derived value bounds the second about three
         !! orders more loosely than the first.
      real(dp) :: linear_dependence = 0.0_dp
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
      logical :: max_iter_set = .false.
         !! Whether the deck named `keywords.scf.maxiter`.
         !!
         !! For the same reason the two convergence flags below exist: a path
         !! whose own default differs from the shared one needs to tell "the
         !! deck asked for 100" from "the deck said nothing". MakeFP is that
         !! path -- it runs 200 iterations by default, and cutting a silent
         !! deck back to 100 would be a regression dressed as a fix.
      logical :: energy_convergence_set = .false.
      logical :: density_convergence_set = .false.
         !! Whether the deck named these, as opposed to inheriting them. A
         !! caller with a stricter default of its own -- MAKEFP -- keeps that
         !! default when they are false and honours the deck when they are true.
      real(dp) :: level_shift = 0.0_dp
         !! Hartree added to the virtual block of the Fock matrix before each
         !! diagonalisation, to widen the gap the next density is built through.
         !! Zero is off. Tapered to zero as the SCF converges, so that the
         !! orbitals and orbital energies reported at exit belong to the
         !! unshifted operator -- see `mqc_libcint_rhf`, where that matters more
         !! than it looks.
      logical :: use_diis = .true.
         !! Use DIIS acceleration
      integer :: diis_size = 8
         !! Number of Fock matrices for DIIS
      logical :: incremental_fock = .true.
         !! Build each iteration's Fock matrix from the density *change* since
         !! the last full build, rather than from the density.
         !!
         !! On by default: it is exact to the convergence threshold and several
         !! times cheaper late in an SCF. False forces a full build every
         !! iteration, which is what to reach for when an SCF misbehaves. The
         !! incremental path accumulates a correction, so it is the first thing
         !! to rule out when a run stalls or wanders -- and rebuilding the binary
         !! to test that is not a debugging step anyone should have to take. It
         !! is also the honest setting for timing a Fock build, since an
         !! incremental one gets cheaper with every iteration.
      character(len=32) :: accelerator = "diis"
      character(len=32) :: convergence_metric = "standard"
         !! `keywords.scf.accelerator`: 'diis' (the default), 'adiis' or
         !! 'ediis'. The energy-based pair runs only while the error is large
         !! and hands over to DIIS, so naming one asks for a different opening,
         !! not a different endgame.
      character(len=32) :: guess = "auto"
         !! Initial guess: 'core', 'gwh', 'sac', 'sad', 'basis_set_projection',
         !! or 'auto'
         !!
         !! 'auto' means the backend picks, because the best starting point
         !! is a property of the backend rather than of the request: the CPU
         !! path resolves it to 'sad', and cuEST to 'gwh', each having
         !! measured its own. An explicit spelling always wins over both.
      type(guess_step_t), allocatable :: guess_steps(:)
         !! The basis ladder for 'basis_set_projection', one entry per
         !! preliminary SCF in order. The target basis is the model's and is not
         !! repeated here.
      logical :: allow_crap_scf = .false.
         !! Accept a non-converged SCF rather than stopping. See mqc_config_types.
      logical :: unrestricted = .false.
         !! Force an unrestricted (UHF/UKS) treatment even for a closed shell.
         !! Needed for broken-symmetry singlets, and the cleanest check that
         !! the unrestricted code reduces to the restricted result.
      character(len=32) :: aux_basis_set = "def2-universal-jkfit"
      logical :: aux_basis_named = .false.
         !! Whether a deck actually asked for the set above, as opposed to
         !! inheriting the default.
         !!
         !! The default is not optional -- cuEST is always density fitted and an
         !! unset auxiliary basis fails there at run time -- which means
         !! `len_trim(aux_basis_set) > 0` is true whether or not anyone asked,
         !! and it cannot be used to decide anything. A double hybrid's
         !! perturbative term did use it that way, and so fitted itself with a
         !! JKFIT set nobody named: 3e-9 on water/STO-3G, and against exactly the
         !! kind of auxiliary `correlation_aux_basis` warns about fitting an
         !! (ia|jb) block with.
         !! Auxiliary (JKFIT) basis for the density-fitted J and K.
         !! Required, not optional, for the cuEST backend: cuEST exposes no
         !! conventional four-index ERI path, so J/K are always fitted.
      logical :: density_fitting = .false.
         !! Fit J and K against `aux_basis_set` on the CPU backend.
         !!
         !! Off by default, and only the CPU backend reads it: cuEST has no
         !! four-index path, so it fits whether or not this is set. libcint has
         !! both, so which one runs has to be asked for rather than inferred
         !! from an auxiliary basis being present -- that name carries a
         !! default, so inferring would mean every calculation silently fitted.
   end type scf_config_t

   !============================================================================
   ! EFP Configuration (MAKEFP)
   !============================================================================
   type :: efp_config_t
      !! What `keywords.efp` carries into a MAKEFP run
      !!
      !! Its own group rather than more keys under `scf` because none of it is an
      !! SCF setting: the SCF a potential runs is already configured by
      !! `keywords.scf.tolerance` and `keywords.scf.density_tolerance`, which reach
      !! it and must keep meaning only that. What is here belongs to the stages
      !! after the SCF -- the response solve and the screening fit -- and those are
      !! where a MAKEFP run spends its time.
      real(dp) :: dynamic_tolerance = DEFAULT_DYNAMIC_TOL
         !! Residual target for the frequency-dependent response solve.
         !!
         !! The one setting in this group that changes the numbers written rather
         !! than the route taken to them. The dynamic polarizabilities, and every
         !! dispersion energy computed from them downstream, are converged to
         !! whatever this says and no further.
      integer :: dynamic_maxiter = DEFAULT_DYNAMIC_MAXITER
         !! Iterations that solve gets per system before it declines to converge.
      integer :: response_batch = DEFAULT_RESPONSE_BATCH
         !! Densities per integral pass in the matrix-free response. A tuning
         !! knob, not physics -- the answer does not change, only how long it
         !! takes. See `DEFAULT_RESPONSE_BATCH`.
      logical :: allow_crap_response = .false.
         !! Take whatever the response solve reached when it runs out of
         !! iterations, instead of refusing.
         !!
         !! Named after `allow_crap_scf` and meant as literally: the potential
         !! that comes out is wrong, and how wrong is not bounded by anything.
         !! It exists because profiling wants a run that reaches every stage
         !! without converging any of them -- `dynamic_maxiter: 1` with this on
         !! is one pass through the whole pipeline -- and because a solve that
         !! stalls on the last few systems is worth seeing the rest of rather
         !! than losing entirely.
      integer :: response = EFP_RESPONSE_AUTO
         !! Build the response operator, never build it, or decide on size.
      real(dp) :: vdw_scale = DEFAULT_VDW_SCALE
         !! Innermost layer of the charge-penetration screening grid, as a
         !! fraction of a van der Waals radius. GAMESS's `VDWSCL`.
   end type efp_config_t

   !============================================================================
   ! XTB Configuration (GFN1, GFN2)
   !============================================================================
   type :: xtb_config_t
      !! Configuration for semi-empirical xTB methods
      real(dp) :: accuracy = 0.01_dp
         !! Numerical accuracy parameter
      real(dp) :: electronic_temp = 300.0_dp
         !! Electronic temperature in Kelvin (Fermi smearing)

      ! Solvation
      character(len=32) :: solvent = ""
         !! Solvent name: "water", "ethanol", etc. Empty for gas phase
      character(len=16) :: solvation_model = ""
         !! Solvation model: "alpb", "gbsa", "cpcm"
      logical :: use_cds = .true.
         !! Include non-polar CDS terms
      logical :: use_shift = .true.
         !! Include solution state shift
      real(dp) :: dielectric = -1.0_dp
         !! Dielectric constant (-1 = use solvent table)
      integer :: cpcm_nang = 110
         !! Angular grid points for CPCM
      real(dp) :: cpcm_rscale = 1.0_dp
         !! Radii scaling for CPCM
   contains
      procedure :: has_solvation => xtb_has_solvation
      procedure :: configure => xtb_configure
      procedure :: get_solvation_info => xtb_get_solvation_info
   end type xtb_config_t

   !============================================================================
   ! DFT Configuration (uses scf_config_t for SCF settings)
   !============================================================================
   type :: dft_config_t
      !! Configuration for Kohn-Sham DFT method
      !! Note: SCF settings (convergence, DIIS) come from scf_config_t
      character(len=32) :: functional = "b3lyp"
         !! XC functional: "lda", "pbe", "b3lyp", "m06-2x", etc.

      ! Integration grid
      character(len=16) :: grid_type = "medium"
         !! Grid quality: "coarse", "medium", "fine", "ultrafine"
      integer :: grid_level = 3
      integer :: nlc_grid_level = -1
         !! `keywords.dft.nlc_grid_level`, the quadrature VV10's double integral
         !! uses. Separate from `grid_level` because the non-local term costs
         !! the *product* of two grid sizes while everything else costs one, so
         !! the level that is right for exchange is an order of magnitude too
         !! expensive here. Negative means the backend's own default.
      real(dp) :: screening_tolerance = 1.0e-12_dp
      integer :: block_size = -1
         !! Standard grid tables, 0 to 9. Three is the usual production default.
      integer :: radial_points = 75
         !! Radial grid points per atom
      integer :: angular_points = 302
         !! Angular grid points (Lebedev)

      ! Density fitting
      logical :: use_density_fitting = .false.
         !! Use RI-J approximation
      character(len=32) :: aux_basis_set = ""
         !! Auxiliary basis for density fitting

      ! Dispersion correction
      logical :: use_dispersion = .false.
         !! Add empirical dispersion
      character(len=8) :: dispersion_type = "d3bj"
         !! Dispersion type: "d3", "d3bj", "d4"
   end type dft_config_t

   !============================================================================
   ! MCSCF/CASSCF Configuration
   !============================================================================
   type :: pcm_config_t
      !! A polarizable continuum: the cavity, the solvent, and the charge solve
      !!
      !! Backend-neutral in shape. The cuEST path hands the cavity and the charge
      !! solve to the library; the CPU path builds both itself in
      !! `mqc_libcint_pcm`. tblite's CPCM is configured through `xtb_config_t`
      !! and builds its own cavity; the two are separate models and share no
      !! settings.
      logical :: enabled = .false.
      character(len=16) :: method = "cpcm"
         !! Which continuum model the surface charges solve: "cpcm"
         !! (conductor-like, scale factor (eps-1)/eps) or "iefpcm" (the integral
         !! equation formalism, scale (eps-1)/(eps+1) with the D-matrix terms).
         !! The two are different models with different energies -- water differs
         !! by ~1% of the solvation energy -- so the choice is named, not
         !! defaulted per backend. Only the CPU path reads it; cuEST's solver is
         !! fixed and refuses "iefpcm" rather than substituting.
      real(dp) :: dielectric = -1.0_dp
         !! Solvent dielectric. No solvent-name table on this path, on purpose:
         !! see `mqc_config_types`.
      integer :: angular_points = 110
         !! Lebedev points per atom on the cavity surface.
      real(dp) :: radii_scale = 1.2_dp
         !! Scaling from van der Waals radii, per `mqc_pcm_radii`.
      real(dp) :: zeta = 2.0_dp
         !! Gaussian switching prefactor for the smooth surface. Unverified
         !! against cuEST's convention -- see `DEFAULT_PCM_ZETA`.
      real(dp) :: tolerance = 1.0e-8_dp
         !! Residual the surface-charge solve must reach.
      integer :: max_iter = 100
         !! Iterations that solve is allowed.
   end type pcm_config_t

   type :: properties_config_t
      !! Analyses to run once the wave function exists
      character(len=32) :: bonding_analysis = "none"
      character(len=:), allocatable :: fukui_population
      character(len=:), allocatable :: fukui_guess
         !! "neutral" (default) or "independent". See `mqc_config_t`.
      type(deltascf_options_t) :: fukui_scf
         !! How the ion SCFs converge. Already resolved against
         !! `keywords.scf` by the time it reaches here -- the reader seeds it
         !! and the deck overwrites what it names -- so a consumer reads these
         !! fields directly and never tests a sentinel.
      character(len=:), allocatable :: charges_scheme
         !! Allocated when `properties.charges` asked for partial charges;
         !! "mulliken" or "chelpg". Unallocated means no charges, which is why
         !! this is allocatable rather than a "none" sentinel like
         !! `bonding_analysis` -- there is no scheme that means "do not".
      logical :: bonding_energy = .false.
      logical :: bonding_no_sharing = .false.
      logical :: bonding_restrict_localization = .false.
      character(len=32) :: bonding_no_sharing_ci = "transform"
      real(dp) :: bonding_threshold = 1.0_dp
   end type properties_config_t

   type, extends(scf_numerics_t) :: scf_options_t
      !! What every self-consistent-field method carries, defined once
      !!
      !! **This exists because the copies drifted and one of them cost a bug.**
      !! `hf_options_t` and `dft_options_t` each held these fields by hand,
      !! and each had to be updated in three places -- the type, its
      !! `configure_*` in the factory, and the `settings%... = options%...`
      !! block in the method. Miss the third and a keyword is read from the
      !! deck, passes the schema, reaches the config and is silently dropped:
      !! that is what happened to `allow_crap_scf` on the Kohn-Sham path,
      !! which was honoured for Hartree-Fock and ignored for DFT. The two
      !! copies had also drifted in spelling -- `conv_tol` against
      !! `energy_tol`, `density_fitting` against `use_density_fitting` -- for
      !! fields meant to mean the same thing.
      !!
      !! **The convergence knobs live one level down, in `scf_numerics_t`.**
      !! `max_iter`, the tolerances, the level shift, DIIS and the accelerator
      !! are inherited from there rather than declared here, so that a second
      !! SCF in the same run -- the deltaSCF ions behind a Fukui analysis --
      !! can carry them without also carrying a basis, a `pcm` and a
      !! `properties` block it has no use for. Reading `options%max_iter` is
      !! unchanged by that; the field is inherited, not moved away.
      !!
      !! Method options **extend** this rather than holding it as a component,
      !! so `options%max_iter` keeps working and a method cannot forget a
      !! field it does not own. Anything specific to one reference -- a
      !! functional, a grid, coupled-cluster settings -- belongs in the
      !! extending type, not here.
      !!
      !! **Adding a field here means adding it in two more places**:
      !! `configure_scf` in `mqc_method_factory`, which fills it from the deck,
      !! and `apply_scf_settings` in `mqc_cuest_iface`, which hands it to the
      !! backend. `test/test_mqc_scf_options.f90` asserts the second of those
      !! for every field, so add it there too -- a field nothing asserts is a
      !! field that can go missing again, and the way it goes missing is
      !! silently.
      character(len=32) :: basis_set = "sto-3g"
         !! Orbital basis set name
      character(len=32) :: ecp_set = ""
         !! Effective core potential set, empty for an all-electron run
      character(len=32) :: aux_basis_set = "def2-universal-jkfit"
         !! Auxiliary (JKFIT) basis for the density-fitted J and K
      logical :: aux_basis_named = .false.
         !! Whether the deck asked for `aux_basis_set`. See `scf_config_t`.
      logical :: cartesian = .false.
         !! `model.cartesian`; see `mqc_config_t`.
      logical :: spherical = .true.
         !! Use spherical (true) or Cartesian (false) basis
      logical :: unrestricted = .false.
         !! Force UHF/UKS even for a closed shell
      logical :: density_fitting = .false.
         !! Fit J and K rather than computing exact integrals
      logical :: freeze_core = .false.
         !! Exclude core orbitals from a post-SCF correlation treatment
      integer :: n_frozen_core = -1
         !! Core orbitals to freeze; -1 counts them from the elements
      character(len=16) :: backend = "auto"
         !! Which backend evaluates the integrals
      integer :: device_rank = 0
         !! Node-local MPI rank, for spreading ranks across a node's GPUs
      logical :: verbose = .false.
         !! Print SCF iterations
      real(dp) :: hessian_displacement = DEFAULT_DISPLACEMENT
         !! Step for a finite-difference Hessian, in Bohr.
         !!
         !! Here because the method is what runs the displacements, and the
         !! deck's value had no route to it: `keywords.hessian.displacement`
         !! reached `driver_config%hessian` and stopped, so every
         !! finite-difference Hessian ran at the default while the log printed
         !! that default back as though it had been asked for.
      type(pcm_config_t) :: pcm
         !! Continuum solvation. A property of the reference rather than of
         !! the functional, which is why it sits here: Hartree-Fock, DFT, MP2
         !! and coupled cluster all reach the backend through an extending
         !! type, and leaving it off one of them meant `keywords.pcm` was read,
         !! validated and then silently dropped for that path.
      type(properties_config_t) :: properties
         !! Population analysis and other post-SCF properties
   end type scf_options_t

   type :: mcscf_config_t
      !! Configuration for MCSCF/CASSCF method

      ! Active space definition
      character(len=MAX_ORBITAL_LABEL_LEN), allocatable :: avas_orbitals(:)
         !! Atomic orbital labels the active space should be built from, e.g.
         !! "N 2p". Unallocated means the space was given by counts instead.
      logical :: full_valence = .false.
         !! Take the whole valence shell as the active space
      real(dp) :: avas_threshold = 0.2_dp
      integer :: n_active_electrons = 0
         !! Number of active electrons
      integer :: n_active_orbitals = 0
         !! Number of active orbitals
      integer :: n_inactive_orbitals = -1
         !! Inactive orbitals (-1 = auto from nelec)
      logical :: optimize_orbitals = .true.
         !! Move the orbitals as well as the CI coefficients.
         !!
         !! True is CASSCF, false is CASCI on whatever the reference SCF
         !! produced. Both spellings parse to `METHOD_TYPE_MCSCF`, so by the
         !! time a method is built the distinction no longer exists in the type
         !! -- it is carried here instead, defaulted from the method name by the
         !! reader and overridable by `keywords.mcscf.optimize_orbitals`. Same
         !! arrangement as coupled cluster's triples.

      ! State averaging
      integer :: n_states = 1
         !! Number of states for SA-CASSCF
      real(dp), allocatable :: state_weights(:)
         !! State weights (must sum to 1)

      ! Convergence
      integer :: max_macro_iter = 100
         !! Maximum orbital optimization iterations
      integer :: max_micro_iter = 50
         !! Maximum CI iterations per macro
      real(dp) :: orbital_convergence = 1.0e-6_dp
         !! Orbital gradient threshold
      real(dp) :: ci_convergence = 1.0e-8_dp
         !! CI energy threshold

      ! Perturbative corrections
      logical :: use_pt2 = .false.
         !! Apply CASPT2/NEVPT2 after CASSCF
      character(len=16) :: pt2_type = "nevpt2"
         !! PT2 flavor: "caspt2", "nevpt2"
      real(dp) :: ipea_shift = 0.25_dp
         !! IPEA shift for CASPT2
      real(dp) :: imaginary_shift = 0.0_dp
         !! Imaginary shift for intruder states
      integer, allocatable :: ormas_subspaces(:)
         !! Active orbital each subspace starts at, from
         !! `keywords.mcscf.ormas`. Unallocated is a complete active space.
      integer, allocatable :: ormas_min_electrons(:)
      integer, allocatable :: ormas_max_electrons(:)
         !! Electrons allowed in each subspace, both spins together
   end type mcscf_config_t

   !============================================================================
   ! Correlation Configuration (shared by MP2, CC, etc.)
   !============================================================================
   type :: correlation_config_t
      !! Shared settings for all post-HF correlation methods
      real(dp) :: energy_convergence = 1.0e-8_dp
         !! Correlation energy convergence threshold

      ! Frozen core
      integer :: n_frozen_core = -1
         !! Number of frozen core orbitals (-1 = auto from elements)
      logical :: freeze_core = .true.
         !! Whether to freeze core orbitals

      ! Density fitting for correlation
      logical :: use_df = .true.
         !! Use density fitting (RI) for correlation integrals

      ! Local correlation
      logical :: use_local = .false.
         !! Use local correlation approximation
      character(len=16) :: local_type = "dlpno"
         !! Local correlation type: "pno", "dlpno", "lmp2", "lno"
      real(dp) :: pno_threshold = 1.0e-7_dp
         !! PNO occupation threshold for truncation

      ! Spin-component scaling (for MP2)
      logical :: use_scs = .false.
         !! Use spin-component scaled MP2
      real(dp) :: scs_ss = 1.0_dp/3.0_dp
         !! Same-spin scaling factor (default: 1/3 for SCS-MP2)
      real(dp) :: scs_os = 1.2_dp
         !! Opposite-spin scaling factor (default: 6/5 for SCS-MP2)
   end type correlation_config_t

   !============================================================================
   ! Coupled Cluster Configuration
   !============================================================================
   type :: cc_config_t
      !! Coupled-cluster specific settings (CCSD, CCSD(T), CC2, CC3, etc.)
      integer :: max_iter = 100
         !! Maximum CC iterations
      real(dp) :: amplitude_convergence = 1.0e-7_dp
         !! T-amplitude convergence threshold

      ! Excitation level
      logical :: include_triples = .false.
         !! Include (T) triples correction
      logical :: perturbative_triples = .true.
         !! Use perturbative (T) vs full CCSDT

      ! DIIS for CC
      logical :: use_diis = .true.
         !! Use DIIS for amplitude equations
      integer :: diis_size = 8
         !! DIIS subspace size
      logical :: spin_adapted = .true.
         !! Spatial-orbital (spin-adapted) formulation rather than spin orbitals.
         !! See mqc_config_types for why this is the default.

      ! EOM-CC for excited states
      integer :: n_roots = 0
         !! Number of EOM-CC roots (0 = ground state only)
      character(len=8) :: eom_type = "ee"
         !! EOM type: "ee" (excitation), "ip" (ionization), "ea" (attachment)
   end type cc_config_t

   !============================================================================
   ! F12 Explicitly Correlated Configuration
   !============================================================================
   type :: f12_config_t
      !! Settings for explicitly correlated F12 methods (MP2-F12, CCSD-F12, etc.)
      real(dp) :: geminal_exponent = 1.0_dp
         !! Slater-type geminal exponent (beta)
      character(len=8) :: ansatz = "3c"
         !! F12 ansatz: "3c", "3c(fix)", "2b", "2a"

      ! Auxiliary basis sets for F12
      character(len=32) :: cabs_basis = ""
         !! Complementary auxiliary basis (CABS) for RI
      character(len=32) :: optri_basis = ""
         !! Optional RI basis for F12 intermediates

      ! Approximations
      logical :: use_exponent_fit = .false.
         !! Fit geminal exponent to basis set
      logical :: scale_triples = .true.
         !! Apply F12 scaling to (T) correction
   end type f12_config_t

   !============================================================================
   ! Main Configuration Type (Composition)
   !============================================================================
   type :: method_config_t
      !! Master configuration containing all method-specific configs
      !!
      !! Usage:
      !!   config%method_type = METHOD_TYPE_DFT
      !!   config%basis_set = 'cc-pvdz'
      !!   config%dft%functional = 'pbe0'
      !!   config%dft%use_dispersion = .true.

      !----- Common settings (all ab-initio methods) -----
      integer(int32) :: method_type = METHOD_TYPE_UNKNOWN
         !! Method type constant
      logical :: verbose = .false.
         !! Enable verbose output
      real(dp) :: hessian_displacement = DEFAULT_DISPLACEMENT
         !! Step for a finite-difference Hessian, in Bohr
      character(len=32) :: basis_set = "sto-3g"
         !! Basis set name (HF, DFT, MCSCF)
      character(len=32) :: ecp_set = ""
         !! Effective core potential set, empty for an all-electron run
      logical :: use_spherical = .true.
         !! Spherical vs Cartesian basis functions
      character(len=16) :: backend = "auto"
         !! Integral backend request; see `parse_backend_name`.
      integer :: device_rank = 0
         !! Node-local MPI rank, used to spread ranks across the GPUs on a
         !! node. Zero is correct for a serial run; a fragmented run must set
         !! the real value or every rank binds to device 0.

      !----- Shared configurations -----
      type(scf_config_t) :: scf
      type(pcm_config_t) :: pcm
      type(properties_config_t) :: properties
         !! Shared SCF settings (used by HF and DFT)
      type(correlation_config_t) :: corr
         !! Shared correlation settings (used by MP2, CC, etc.)

      !----- Method-specific configurations -----
      type(xtb_config_t) :: xtb
         !! XTB settings (GFN1, GFN2)
      type(dft_config_t) :: dft
         !! DFT-specific settings (functional, grid, dispersion)
      type(mcscf_config_t) :: mcscf
         !! MCSCF/CASSCF settings
      type(cc_config_t) :: cc
         !! Coupled-cluster specific settings (CCSD, CCSD(T), etc.)
      type(f12_config_t) :: f12
         !! F12 explicitly correlated settings
      type(efp_config_t) :: efp
         !! MAKEFP settings: the response solve and the screening grid

   contains
      procedure :: reset => config_reset
      procedure :: log_settings => config_log_settings
   end type method_config_t

contains

   pure function xtb_has_solvation(this) result(has_solvation)
      !! Check if solvation is configured for XTB
      class(xtb_config_t), intent(in) :: this
      logical :: has_solvation
      has_solvation = len_trim(this%solvent) > 0 .or. this%dielectric > 0.0_dp
   end function xtb_has_solvation

   subroutine xtb_configure(this, use_cds, use_shift, dielectric, cpcm_nang, cpcm_rscale, solvent, solvation_model)
      !! Configure XTB solvation settings from driver configuration
      !!
      !! Sets up all XTB-specific parameters and applies default solvation model
      !! logic (alpb if solvent or dielectric is specified but no model given).
      class(xtb_config_t), intent(inout) :: this
      logical, intent(in) :: use_cds             !! Include CDS non-polar terms
      logical, intent(in) :: use_shift           !! Include solution state shift
      real(dp), intent(in) :: dielectric         !! Direct dielectric constant (-1 = use solvent lookup)
      integer, intent(in) :: cpcm_nang           !! Angular grid points for CPCM
      real(dp), intent(in) :: cpcm_rscale        !! Radii scaling for CPCM
      character(len=*), intent(in), optional :: solvent          !! Solvent name
      character(len=*), intent(in), optional :: solvation_model  !! Solvation model name

      this%use_cds = use_cds
      this%use_shift = use_shift
      this%dielectric = dielectric
      this%cpcm_nang = cpcm_nang
      this%cpcm_rscale = cpcm_rscale

      if (present(solvent)) this%solvent = solvent
      if (present(solvation_model)) then
         this%solvation_model = solvation_model
      else if (len_trim(this%solvent) > 0 .or. dielectric > 0.0_dp) then
         this%solvation_model = "alpb"  ! Default solvation model
      end if
   end subroutine xtb_configure

   subroutine xtb_get_solvation_info(this, info_lines, n_lines)
      !! Get solvation configuration info for logging
      !!
      !! Returns array of info strings describing the solvation setup.
      !! If no solvation is configured, n_lines = 0.
      use pic_io, only: to_char
      class(xtb_config_t), intent(in) :: this
      character(len=128), intent(out) :: info_lines(4)  !! Up to 4 lines of info
      integer, intent(out) :: n_lines                    !! Number of lines populated

      n_lines = 0
      info_lines = ""

      if (.not. this%has_solvation()) return

      n_lines = 1
      if (trim(this%solvation_model) == "cpcm") then
         if (this%dielectric > 0.0_dp) then
            info_lines(1) = "XTB solvation enabled: cpcm with dielectric = "//to_char(this%dielectric)
         else
            info_lines(1) = "XTB solvation enabled: cpcm with "//trim(this%solvent)
         end if
         n_lines = 3
         info_lines(2) = "  CPCM grid points (nang): "//to_char(this%cpcm_nang)
         info_lines(3) = "  CPCM radii scale: "//to_char(this%cpcm_rscale)
      else
         info_lines(1) = "XTB solvation enabled: "//trim(this%solvation_model)//" with "//trim(this%solvent)
         if (this%use_cds) then
            n_lines = n_lines + 1
            info_lines(n_lines) = "  CDS (non-polar) terms: enabled"
         end if
         if (this%use_shift) then
            n_lines = n_lines + 1
            info_lines(n_lines) = "  Solution state shift: enabled"
         end if
      end if
   end subroutine xtb_get_solvation_info

   subroutine config_reset(this)
      !! Reset all configuration values to defaults
      class(method_config_t), intent(inout) :: this

      ! Common settings
      this%method_type = METHOD_TYPE_UNKNOWN
      this%verbose = .false.
      this%hessian_displacement = DEFAULT_DISPLACEMENT
      this%basis_set = "sto-3g"
      this%use_spherical = .true.

      ! XTB defaults
      this%xtb%accuracy = 0.01_dp
      this%xtb%electronic_temp = 300.0_dp
      this%xtb%solvent = ""
      this%xtb%solvation_model = ""
      this%xtb%use_cds = .true.
      this%xtb%use_shift = .true.
      this%xtb%dielectric = -1.0_dp
      this%xtb%cpcm_nang = 110
      this%xtb%cpcm_rscale = 1.0_dp

      ! SCF defaults (shared by HF and DFT)
      this%scf%max_iter = 100
      this%scf%energy_convergence = 1.0e-8_dp
      this%scf%density_convergence = 1.0e-6_dp
      this%scf%use_diis = .true.
      this%scf%diis_size = 8

      ! DFT-specific defaults
      this%dft%functional = "b3lyp"
      this%dft%grid_type = "medium"
      this%dft%radial_points = 75
      this%dft%angular_points = 302
      this%dft%use_density_fitting = .false.
      this%dft%aux_basis_set = ""
      this%scf%aux_basis_named = .false.
      this%dft%use_dispersion = .false.
      this%dft%dispersion_type = "d3bj"

      ! MCSCF defaults
      if (allocated(this%mcscf%ormas_subspaces)) deallocate (this%mcscf%ormas_subspaces)
      if (allocated(this%mcscf%ormas_min_electrons)) then
         deallocate (this%mcscf%ormas_min_electrons)
      end if
      if (allocated(this%mcscf%ormas_max_electrons)) then
         deallocate (this%mcscf%ormas_max_electrons)
      end if
      this%mcscf%full_valence = .false.
      this%mcscf%n_active_electrons = 0
      this%mcscf%n_active_orbitals = 0
      this%mcscf%n_inactive_orbitals = -1
      this%mcscf%optimize_orbitals = .true.
      this%mcscf%n_states = 1
      if (allocated(this%mcscf%state_weights)) deallocate (this%mcscf%state_weights)
      this%mcscf%max_macro_iter = 100
      this%mcscf%max_micro_iter = 50
      this%mcscf%orbital_convergence = 1.0e-6_dp
      this%mcscf%ci_convergence = 1.0e-8_dp
      this%mcscf%use_pt2 = .false.
      this%mcscf%pt2_type = "nevpt2"
      this%mcscf%ipea_shift = 0.25_dp
      this%mcscf%imaginary_shift = 0.0_dp

      ! Correlation defaults (shared by MP2, CC, etc.)
      this%corr%energy_convergence = 1.0e-8_dp
      this%corr%n_frozen_core = -1
      this%corr%freeze_core = .true.
      this%corr%use_df = .true.
      this%corr%use_local = .false.
      this%corr%local_type = "dlpno"
      this%corr%pno_threshold = 1.0e-7_dp
      this%corr%use_scs = .false.
      this%corr%scs_ss = 1.0_dp/3.0_dp
      this%corr%scs_os = 1.2_dp

      ! Coupled-cluster defaults
      this%cc%max_iter = 100
      this%cc%amplitude_convergence = 1.0e-7_dp
      this%cc%include_triples = .false.
      this%cc%perturbative_triples = .true.
      this%cc%use_diis = .true.
      this%cc%diis_size = 8
      this%cc%n_roots = 0
      this%cc%eom_type = "ee"

      ! F12 defaults
      this%f12%geminal_exponent = 1.0_dp
      this%f12%ansatz = "3c"
      this%f12%cabs_basis = ""
      this%f12%optri_basis = ""
      this%f12%use_exponent_fit = .false.
      this%f12%scale_triples = .true.
   end subroutine config_reset

   subroutine config_log_settings(this)
      !! Log method-specific settings based on method type
      !!
      !! Dispatches to appropriate logging based on the configured method.
      !! This allows the driver to log settings without knowing method details.
      use mqc_method_types, only: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2, METHOD_TYPE_HF, &
                                  METHOD_TYPE_DFT, METHOD_TYPE_MCSCF, method_type_to_string
      use pic_logger, only: logger => global_logger
      class(method_config_t), intent(in) :: this

      character(len=128) :: info_lines(4)
      integer :: n_lines, i

      select case (this%method_type)
      case (METHOD_TYPE_GFN1, METHOD_TYPE_GFN2)
         ! Semi-empirical: no basis to report, but solvation matters
         call this%xtb%get_solvation_info(info_lines, n_lines)
         do i = 1, n_lines
            call logger%info(trim(info_lines(i)))
         end do

      case (METHOD_TYPE_HF, METHOD_TYPE_DFT, METHOD_TYPE_MCSCF)
         ! Which basis a run actually used is the first thing anyone needs when
         ! comparing two energies, so report it rather than leaving it implicit.
         call logger%info("Method settings:")
         call logger%info("  Method:          "//trim(method_type_to_string(this%method_type)))
         if (this%method_type == METHOD_TYPE_DFT) then
            call logger%info("  Functional:      "//trim(this%dft%functional))
         end if
         call logger%info("  Basis set:       "//trim(this%basis_set))
         ! Only when it is actually fitting something. The aux basis has a
         ! non-empty default, so printing it unconditionally told every exact run
         ! it had an auxiliary basis it was not using -- which read as "this is
         ! density-fitted" when it was not.
         if (this%scf%density_fitting) then
            call logger%info("  Density fitting: "//trim(this%scf%aux_basis_set))
         end if
         ! The angular form is not reported here any more, and deliberately.
         ! `use_spherical` is the flag cuEST builds its AO shells with, not a
         ! statement about the basis: the basis file decides, through its
         ! `function_type` entries, and 6-31G* is Cartesian whatever this says.
         ! Printing this line beside a Cartesian basis claimed "spherical
         ! (pure)" for a run that was not, which is the reporting half of the
         ! bug it sat next to. The SCF prints the form it actually used, along
         ! with the function count that distinguishes the two.
         if (this%method_type == METHOD_TYPE_DFT) then
            block
               character(len=80) :: grid_line
               ! Reported as the level unless a deck overrode it with explicit
               ! counts, because the level is what the CPU path actually uses --
               ! printing "75 radial x 302" for a level-3 run was true of neither
               ! backend and is the reporting half of a bug rather than a bug's
               ! symptom. The counts still appear when they are the thing in force.
               if (this%dft%radial_points > 0 .and. this%dft%angular_points > 0 &
                   .and. this%dft%grid_level < 0) then
                  write (grid_line, "(a,i0,a,i0)") "  XC grid:         ", &
                     this%dft%radial_points, " radial x ", this%dft%angular_points
               else
                  write (grid_line, "(a,i0)") "  XC grid:         level ", &
                     this%dft%grid_level
               end if
               call logger%info(trim(grid_line))
            end block
         end if

      case default
         return
      end select

   end subroutine config_log_settings

end module mqc_method_config
