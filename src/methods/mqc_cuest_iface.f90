!! Backend-independent description of a cuEST SCF calculation
module mqc_cuest_iface
   !! Holds `cuest_scf_settings_t`, the settings an HF or DFT method hands to
   !! the cuEST backend.
   !!
   !! It lives here, in `src/`, rather than with the backend for one reason:
   !! the backend cannot be part of the fpm build (fpm cannot link the cuEST
   !! library, and its module-dependency scanner is preprocessor-blind, so it
   !! would demand the backend modules even from behind an `#ifdef`). Keeping
   !! the *type* here and reaching the backend through `mqc_cuest_bridge` --
   !! which has a stub form fpm compiles and a real form CMake compiles -- lets
   !! the method files carry no preprocessor conditionals at all.
   use pic_types, only: dp
   use mqc_config_types, only: guess_step_t, deltascf_options_t
   use mqc_method_config, only: pcm_config_t, mcscf_config_t, scf_options_t, properties_config_t
   implicit none
   private

   public :: cuest_scf_settings_t
   public :: apply_scf_settings
   public :: apply_properties_settings
   public :: BACKEND_AUTO, BACKEND_CUEST, BACKEND_LIBCINT
   public :: parse_backend_name
   public :: method_runs_on_cuest

   !> Which integral backend a deck asked for.
   !>
   !> `auto` is the historical behaviour and the default: cuEST when the build
   !> has it, the CPU path otherwise. The other two are requests, and a request
   !> that cannot be honoured is refused rather than quietly substituted -- which
   !> is the whole reason for naming one. A deck that says `cuest` and silently
   !> got the CPU path would report timings and a provenance that were not true.
   integer, parameter :: BACKEND_AUTO = 0
   integer, parameter :: BACKEND_CUEST = 1
   integer, parameter :: BACKEND_LIBCINT = 2

   type :: cuest_scf_settings_t
      !! Method-independent description of one cuEST SCF calculation
      character(len=32) :: basis_set = "sto-3g"
         !! Orbital basis set name
      logical :: cartesian = .false.
         !! Read the basis in Cartesian form whatever its file declares; see
         !! `mqc_config_t`. Only the libcint path acts on it.
      character(len=32) :: ecp_set = ""
         !! `model.ecp`, the effective core potential set, empty for none.
         !!
         !! Named separately from `basis_set` because the two are separate
         !! files and a basis does not imply a potential: def2-SVP is used with
         !! def2-ECP above krypton and with nothing below it, and asking for
         !! the potential on a light element is not an error -- the reader
         !! returns no channels and the term is zero.
      character(len=32) :: aux_basis_set = "def2-universal-jkfit"
      logical :: density_fitting = .false.
      ! Post-Hartree-Fock. Kept beside the SCF settings rather than inside
      ! them because they are not SCF settings: a density-fitted reference and
      ! a conventional correlation treatment is a combination someone will ask
      ! for, and the two flags have to be able to disagree.
      logical :: run_mp2 = .false.
      logical :: aux_basis_named = .false.
         !! Whether `aux_basis_set` was asked for or merely defaulted. The
         !! default exists because cuEST needs one, so its presence says nothing.
      logical :: freeze_core = .false.
      integer :: n_frozen_core = -1     !! -1 means count it from the elements
      logical :: corr_density_fitting = .false.  !! RI for the correlation step
      real(dp) :: scs_ss = 1.0_dp       !! Spin-component scaling, one for plain MP2
      real(dp) :: scs_os = 1.0_dp
      ! Coupled cluster. CPU backend only -- cuEST has no CC -- and refused
      ! rather than ignored there, so a deck cannot ask the GPU for CCSD and be
      ! handed Hartree-Fock.
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
      ! every atom, which is what a convergence study wants.
      integer :: grid_level = 3
      integer :: nlc_grid_level = -1
         !! VV10's quadrature level; negative means the backend default.
      real(dp) :: screening_tolerance = 1.0e-12_dp
      integer :: block_size = -1
      integer :: backend = BACKEND_AUTO
         !! One of BACKEND_*. Resolved from the deck's `backend` key.

      ! A polarizable continuum, when one was asked for. Carried whole rather
      ! than field by field: the cavity, the solvent and the charge solve travel
      ! together, and a backend either builds a continuum or does not.
      character(len=32) :: bonding_analysis = "none"
      character(len=:), allocatable :: fukui_population
      character(len=:), allocatable :: fukui_guess
         !! "neutral" (default) or "independent"; see `mqc_config_t`.
      type(deltascf_options_t) :: fukui_scf
         !! How the ion SCFs converge, already resolved against
         !! `keywords.scf`; see `scf_numerics_t`.
      character(len=:), allocatable :: charges_scheme
         !! Atomic partial charges from `properties.charges`, unallocated when
         !! none were asked for. Travels here for the same reason the bonding
         !! analysis does: partitioning a density needs the molecule and the
         !! converged density, which only the backend holds.
      logical :: bonding_energy = .false.
      logical :: bonding_no_sharing = .false.
      logical :: bonding_restrict_localization = .false.
      character(len=32) :: bonding_no_sharing_ci = "transform"
      real(dp) :: bonding_threshold = 1.0_dp
         !! kcal/mol; pairs weaker than this are counted rather than printed.
         !! A post-SCF analysis to run once the orbitals are converged, from
         !! `properties.bonding_analysis`. Travels with the SCF settings because
         !! it needs what only the backend has -- the molecule and the converged
         !! orbitals -- and because it is not a fragment property: it is a
         !! statement about the wave function this object describes.

      type(pcm_config_t) :: pcm
         !! Read by the CPU backend only; cuEST always fits.

      ! The active space, when a multiconfigurational method is what is being
      ! run. Carried whole for the same reason the continuum is: the active
      ! space, the orbital partition and the convergence thresholds only mean
      ! anything together, and every other SCF setting a CASSCF needs -- the
      ! basis, the reference SCF's iteration cap and tolerances -- is already
      ! on this type. A CASSCF *is* an SCF followed by more work, so it takes
      ! the same settings object rather than a parallel one.
      type(mcscf_config_t) :: mcscf
         !! Read by the CPU backend only; cuEST has no CI.
         !! Auxiliary (JKFIT) basis. Required: cuEST fits J and K always.
      character(len=32) :: functional = ""
         !! Exchange-correlation functional; empty means Hartree-Fock
      logical :: spherical = .true.
         !! Pure (spherical) vs Cartesian angular functions
      logical :: verbose = .false.
         !! Print the SCF iteration table
      character(len=32) :: accelerator = "diis"
         !! `keywords.scf.accelerator`: 'diis' (the default), 'adiis' or
         !! 'ediis'. The energy-based pair runs only while the error is large
         !! and hands over to DIIS, so naming one asks for a different opening,
         !! not a different endgame.
      character(len=32) :: guess = "auto"
         !! Initial guess: 'core', 'gwh', 'sac', 'sad', 'basis_set_projection',
         !! or 'auto'
      type(guess_step_t), allocatable :: guess_steps(:)
         !! The basis ladder for 'basis_set_projection', one entry per
         !! preliminary SCF in order. The target basis is `basis_set` and is not
         !! repeated, so STO-3G then 6-31G then cc-pVTZ is two steps.
         !!
         !! Carried on the shared settings type because that is where every
         !! other SCF setting lives, but only the libcint backend acts on it --
         !! the projection is built from `libcint_molecule_t` and the cuEST path
         !! refuses the guess rather than silently running a different one.
         !!
         !! 'auto' means the backend picks, because the best starting point
         !! is a property of the backend rather than of the request: the CPU
         !! path resolves it to 'sad', and cuEST to 'gwh', each having
         !! measured its own. An explicit spelling always wins over both.
      integer :: device_rank = 0
         !! Node-local MPI rank; decides which GPU this rank binds to
      logical :: unrestricted = .false.
         !! Force UHF/UKS even for a closed shell

      logical :: allow_crap_scf = .false.
         !! Keep a non-converged SCF instead of failing the fragment. Same
         !! meaning and same default as the xTB path -- one physical condition
         !! must not have two behaviours depending on which backend ran it.
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
      !! The sibling of `apply_scf_settings`, and it exists for the same
      !! reason: this block was written out by hand in each method and the
      !! three had already diverged. Kohn-Sham unpacked four of the eight and
      !! MCSCF six, so `keywords.properties.bonding_energy` on a Kohn-Sham run
      !! reached a backend that reads `settings%bonding_energy` -- the same
      !! routine Hartree-Fock uses -- and found the default. That one is worse
      !! than a dropped feature: the flag guards a *refusal*, because the
      !! bonding energy decomposition rebuilds atom energies from gas-phase
      !! operators and must not run under a continuum, so never setting it
      !! disarmed a check that exists to stop a wrong answer.
      !!
      !! `fukui_population` and `charges_scheme` are allocatable and stay
      !! guarded: an unset one must not overwrite what a backend already has.
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
      !! **The third and last copy of the shared settings.** A field reaching
      !! a backend crossed three hand-maintained layers: the options type, the
      !! `configure_*` that fills it, and this block. `scf_options_t` and
      !! `configure_scf` removed the first two; without this one a field could
      !! be declared on both methods and filled on both and still never arrive,
      !! which is exactly what `allow_crap_scf` did on the Kohn-Sham path --
      !! the deck said true, the options object said true, and the backend was
      !! never told.
      !!
      !! `backend` is deliberately not here: parsing its name can fail, and a
      !! routine that copies fields should not also be the one that reports an
      !! error. The caller keeps that line.
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
      settings%verbose = options%verbose
      settings%device_rank = options%device_rank
      settings%pcm = options%pcm
   end subroutine apply_scf_settings

   subroutine parse_backend_name(name, kind, error)
      !! A deck's backend name to one of BACKEND_*
      !!
      !! Spelled by what the thing is rather than where it runs, with `gpu` and
      !! `cpu` accepted because that is how people say it. An unknown name is
      !! refused rather than treated as `auto`: a typo that fell back to the
      !! default would be a deck asking for one backend and getting another.
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
         kind = BACKEND_LIBCINT
      case default
         call error%set(ERROR_VALIDATION, "unknown backend '"//lower//"'; expected "// &
                        "'auto', 'cuest' (or 'gpu'), or 'libcint' (or 'cpu')")
      end select
   end subroutine parse_backend_name

   pure function method_runs_on_cuest(method_type) result(offloadable)
      !! Whether cuEST has an implementation of this method
      !!
      !! An allow-list rather than a list of the refused ones, so a method added
      !! later is refused on the GPU until someone has written it there. The
      !! other way round it would inherit "offloadable" and silently be run by
      !! whatever the cuEST path does with a method it does not know.
      !!
      !! Hartree-Fock and Kohn-Sham are the two: cuEST computes the integrals
      !! and the SCF, and everything beyond it -- MP2, coupled cluster, MCSCF --
      !! is CPU-only here, as are the semi-empirical methods, which never touch
      !! this backend at all.
      use mqc_method_types, only: METHOD_TYPE_HF, METHOD_TYPE_DFT
      integer, intent(in) :: method_type
      logical :: offloadable

      offloadable = (method_type == METHOD_TYPE_HF .or. &
                     method_type == METHOD_TYPE_DFT)
   end function method_runs_on_cuest

end module mqc_cuest_iface
