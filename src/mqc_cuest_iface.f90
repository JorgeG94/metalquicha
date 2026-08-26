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
   use mqc_config_types, only: guess_step_t
   use mqc_method_config, only: pcm_config_t, mcscf_config_t
   implicit none
   private

   public :: cuest_scf_settings_t
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
      real(dp) :: level_shift = 0.0_dp
         !! Hartree added to the virtual block before each diagonalisation, and
         !! tapered to zero before the SCF exits. See `scf_config_t`.
      real(dp) :: linear_dependence = 0.0_dp
         !! Overlap eigenvalues at or below this are dropped as linearly
         !! dependent. Zero means the orthogonaliser's own cutoff. See
         !! `scf_config_t`.
      logical :: use_diis = .true.
      integer :: diis_size = 8

      integer :: radial_points = 75    !! XC grid radial points per atom
      integer :: angular_points = 302  !! XC grid Lebedev order
   end type cuest_scf_settings_t

contains

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
