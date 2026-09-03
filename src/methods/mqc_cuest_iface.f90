!! Backend-independent description of a cuEST SCF calculation
module mqc_cuest_iface
   !! Holds `cuest_scf_settings_t`, the settings an HF or DFT method hands to
   !! the cuEST backend.
   ! Here in `src/` rather than with the backend because the backend cannot be
   ! part of the fpm build: fpm cannot link cuEST, and its module-dependency
   ! scanner is preprocessor-blind, so it would demand the backend modules even
   ! from behind an `#ifdef`.
   use pic_types, only: dp
   use mqc_scf_types, only: guess_step_t, deltascf_options_t
   use mqc_method_config, only: pcm_config_t, mcscf_config_t, scf_options_t, properties_config_t
   implicit none
   private

   public :: cuest_scf_settings_t
   public :: apply_scf_settings
   public :: apply_properties_settings
   public :: BACKEND_AUTO, BACKEND_CUEST, BACKEND_CZT
   public :: parse_backend_name
   public :: method_runs_on_cuest

   integer, parameter :: BACKEND_AUTO = 0
   integer, parameter :: BACKEND_CUEST = 1
   integer, parameter :: BACKEND_CZT = 2
      !! Which integral backend a deck asked for. `auto` is the default: cuEST
      !! when the build has it, the CPU path otherwise. The other two are
      !! requests, and a request that cannot be honoured is refused rather than
      !! quietly substituted.

   type :: cuest_scf_settings_t
      !! Method-independent description of one cuEST SCF calculation
      character(len=32) :: basis_set = "sto-3g"
         !! Orbital basis set name
      logical :: cartesian = .false.
         !! Read the basis in Cartesian form whatever its file declares; see
         !! `mqc_config_t`. Only the cenzontle path acts on it.
      character(len=32) :: ecp_set = ""
         !! `model.ecp`, the effective core potential set, empty for none. A
         !! basis does not imply a potential, and asking for one on a light
         !! element is not an error -- the reader returns no channels and the
         !! term is zero.
      character(len=32) :: aux_basis_set = "def2-universal-jkfit"
         !! Auxiliary (JKFIT) basis. Required by cuEST, which fits J and K
         !! always.
      logical :: density_fitting = .false.
      ! `density_fitting` and `corr_density_fitting` are separate flags on
      ! purpose: they have to be able to disagree.
      logical :: run_mp2 = .false.
      logical :: aux_basis_named = .false.
         !! Whether `aux_basis_set` was asked for or merely defaulted. The
         !! default exists because cuEST needs one, so its presence says nothing.
      logical :: freeze_core = .true.
         !! Matches `correlation_config_t`; see `scf_options_t`.
      integer :: n_frozen_core = -1     !! -1 means count it from the elements
      logical :: corr_density_fitting = .false.  !! RI for the correlation step
      real(dp) :: scs_ss = 1.0_dp       !! Spin-component scaling, one for plain MP2
      real(dp) :: scs_os = 1.0_dp
      ! Coupled cluster. CPU backend only; cuEST refuses it rather than
      ! silently handing back Hartree-Fock.
      logical :: run_cc = .false.
      logical :: cc_triples = .false.   !! Run the perturbative (T) after CCSD
      integer :: cc_max_iter = 100
      real(dp) :: cc_tolerance = 1.0e-8_dp
      integer :: cc_diis_size = 8
      logical :: cc_spin_adapted = .true.
         !! Run coupled cluster in spatial orbitals rather than spin orbitals.
         !! CPU backend only; cuEST has no coupled cluster at all.
      ! Kohn-Sham grid. `grid_level` picks per-element radial and angular counts
      ! from the standard tables; radial_points/angular_points override it for
      ! every atom.
      integer :: grid_level = 3
      integer :: nlc_grid_level = -1
         !! VV10's quadrature level; negative means the backend default.
      real(dp) :: screening_tolerance = 1.0e-12_dp
      integer :: block_size = -1
      integer :: backend = BACKEND_AUTO
         !! One of BACKEND_*. Resolved from the deck's `backend` key.

      character(len=32) :: bonding_analysis = "none"
         !! Post-SCF bonding analysis from `properties.bonding_analysis`, run
         !! once the orbitals are converged. "none" for none.
      character(len=:), allocatable :: fukui_population
      character(len=:), allocatable :: fukui_guess
         !! "neutral" (default) or "independent"; see `mqc_config_t`.
      type(deltascf_options_t) :: fukui_scf
         !! How the ion SCFs converge, already resolved against
         !! `keywords.scf`; see `scf_numerics_t`.
      character(len=:), allocatable :: charges_scheme
         !! Atomic partial charges from `properties.charges`, unallocated when
         !! none were asked for.
      logical :: bonding_energy = .false.
      logical :: bonding_no_sharing = .false.
      logical :: bonding_restrict_localization = .false.
      character(len=32) :: bonding_no_sharing_ci = "transform"
      real(dp) :: bonding_threshold = 1.0_dp
         !! kcal/mol; pairs weaker than this are counted rather than printed.

      type(pcm_config_t) :: pcm
         !! The polarizable continuum, when one was asked for. Read by the CPU
         !! backend only.

      type(mcscf_config_t) :: mcscf
         !! The active space, when a multiconfigurational method is being run.
         !! Read by the CPU backend only; cuEST has no CI.
      character(len=32) :: functional = ""
         !! Exchange-correlation functional; empty means Hartree-Fock
      logical :: spherical = .true.
         !! Pure (spherical) vs Cartesian angular functions
      logical :: verbose = .false.
         !! Print the SCF iteration table
      character(len=32) :: accelerator = "diis"
         !! `keywords.scf.accelerator`: 'diis' (the default), 'adiis' or
         !! 'ediis'. The energy-based pair runs only while the error is large
         !! and hands over to DIIS, so naming one chooses a different opening,
         !! not a different endgame.
      character(len=32) :: convergence_metric = "standard"
         !! See `mqc_scf_convergence`.
      character(len=32) :: guess = "auto"
         !! Initial guess: 'core', 'gwh', 'sac', 'sad', 'basis_set_projection',
         !! or 'auto'. 'auto' lets the backend pick -- the CPU path resolves it
         !! to 'sad' and cuEST to 'gwh'; an explicit spelling wins over both.
      type(guess_step_t), allocatable :: guess_steps(:)
         !! The basis ladder for 'basis_set_projection', one entry per
         !! preliminary SCF in order. The target basis is `basis_set` and is not
         !! repeated, so STO-3G then 6-31G then cc-pVTZ is two steps.
         !!
         !! Only the cenzontle backend acts on it; the cuEST path refuses the
         !! guess rather than silently running a different one.
      integer :: device_rank = 0
         !! Node-local MPI rank; decides which GPU this rank binds to
      logical :: unrestricted = .false.
         !! Force UHF/UKS even for a closed shell

      logical :: allow_crap_scf = .false.
         !! Keep a non-converged SCF instead of failing the fragment. Same
         !! meaning and same default as the xTB path.
      integer :: max_iter = 100
      real(dp) :: energy_tol = 1.0e-8_dp
      real(dp) :: density_tol = 1.0e-6_dp
      real(dp) :: grad_tol = 0.0_dp
         !! Commutator threshold; zero derives `sqrt(energy_tol)`.
      real(dp) :: level_shift = 0.0_dp
         !! Hartree added to the virtual block before each diagonalisation, and
         !! tapered to zero before the SCF exits. See `scf_config_t`.
      real(dp) :: linear_dependence = 0.0_dp
         !! Overlap eigenvalues at or below this are dropped as linearly
         !! dependent. Zero means the orthogonaliser's own cutoff. See
         !! `scf_config_t`.
      logical :: use_diis = .true.
      integer :: diis_size = 8
      logical :: incremental_fock = .true.
         !! Build from the density change rather than the density. See
         !! `scf_config_t`; false forces a full build every iteration.

      integer :: radial_points = 75    !! XC grid radial points per atom
      integer :: angular_points = 302  !! XC grid Lebedev order
   end type cuest_scf_settings_t

contains

   subroutine apply_properties_settings(settings, properties)
      !! Hand a backend the post-SCF properties the deck asked for
      !!
      !! `fukui_population`, `fukui_guess` and `charges_scheme` are allocatable
      !! and stay guarded: an unset one must not overwrite what a backend
      !! already has.
      type(cuest_scf_settings_t), intent(inout) :: settings
      type(properties_config_t), intent(in) :: properties

      settings%bonding_analysis = properties%bonding_analysis
      settings%bonding_threshold = properties%bonding_threshold
      settings%bonding_energy = properties%bonding_energy
      settings%bonding_no_sharing = properties%bonding_no_sharing
      settings%bonding_no_sharing_ci = properties%bonding_no_sharing_ci
      settings%bonding_restrict_localization = properties%bonding_restrict_localization
      if (allocated(properties%fukui_population)) then
         settings%fukui_population = properties%fukui_population
      end if
      if (allocated(properties%fukui_guess)) then
         settings%fukui_guess = properties%fukui_guess
      end if
      settings%fukui_scf = properties%fukui_scf
      if (allocated(properties%charges_scheme)) then
         settings%charges_scheme = properties%charges_scheme
      end if
   end subroutine apply_properties_settings

   subroutine apply_scf_settings(settings, options)
      !! Hand a backend everything the SCF half of a method decided
      !!
      !! **`backend` is not copied here.** Parsing its name can fail, so the
      !! caller calls `parse_backend_name` itself and reports the error.
      type(cuest_scf_settings_t), intent(inout) :: settings
      class(scf_options_t), intent(in) :: options

      settings%basis_set = options%basis_set
      settings%ecp_set = options%ecp_set
      settings%aux_basis_set = options%aux_basis_set
      settings%aux_basis_named = options%aux_basis_named
      settings%cartesian = options%cartesian
      settings%spherical = options%spherical
      settings%density_fitting = options%density_fitting
      settings%freeze_core = options%freeze_core
      settings%n_frozen_core = options%n_frozen_core
      settings%unrestricted = options%unrestricted
      settings%guess = options%guess
      if (allocated(options%guess_steps)) settings%guess_steps = options%guess_steps
      settings%max_iter = options%max_iter
      settings%allow_crap_scf = options%allow_crap_scf
      settings%energy_tol = options%energy_tol
      settings%density_tol = options%density_tol
      settings%grad_tol = options%grad_tol
      settings%level_shift = options%level_shift
      settings%linear_dependence = options%linear_dependence
      settings%use_diis = options%use_diis
      settings%diis_size = options%diis_size
      settings%incremental_fock = options%incremental_fock
      settings%accelerator = options%accelerator
      settings%convergence_metric = options%convergence_metric
      settings%verbose = options%verbose
      settings%device_rank = options%device_rank
      settings%pcm = options%pcm
   end subroutine apply_scf_settings

   subroutine parse_backend_name(name, kind, error)
      !! A deck's backend name to one of BACKEND_*
      !!
      !! Accepts 'auto', 'cuest' (or 'gpu') and 'libcint' (or 'cpu'), in any
      !! case. An unknown name is refused rather than treated as 'auto'.
      use mqc_error, only: error_t, ERROR_VALIDATION
      character(len=*), intent(in) :: name
      integer, intent(out) :: kind
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: lower
      integer :: i

      lower = trim(adjustl(name))
      do i = 1, len(lower)
         if (lower(i:i) >= "A" .and. lower(i:i) <= "Z") then
            lower(i:i) = achar(iachar(lower(i:i)) + 32)
         end if
      end do

      kind = BACKEND_AUTO
      select case (lower)
      case ("", "auto")
         kind = BACKEND_AUTO
      case ("cuest", "gpu")
         kind = BACKEND_CUEST
      case ("libcint", "cpu")
         kind = BACKEND_CZT
      case default
         call error%set(ERROR_VALIDATION, "unknown backend '"//lower//"'; expected "// &
                        "'auto', 'cuest' (or 'gpu'), or 'libcint' (or 'cpu')")
      end select
   end subroutine parse_backend_name

   pure function method_runs_on_cuest(method_type) result(offloadable)
      !! Whether cuEST has an implementation of this method
      !!
      !! An allow-list: Hartree-Fock and Kohn-Sham are the two, so a method
      !! added later is refused on the GPU until someone writes it there.
      !! Everything beyond the SCF -- MP2, coupled cluster, MCSCF -- and the
      !! semi-empirical methods are CPU-only.
      use mqc_method_types, only: METHOD_TYPE_HF, METHOD_TYPE_DFT
      integer, intent(in) :: method_type
      logical :: offloadable

      offloadable = (method_type == METHOD_TYPE_HF .or. &
                     method_type == METHOD_TYPE_DFT)
   end function method_runs_on_cuest

end module mqc_cuest_iface
