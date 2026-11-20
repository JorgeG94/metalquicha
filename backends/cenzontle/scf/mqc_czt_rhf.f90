!! Closed-shell SCF on the CPU
module mqc_czt_rhf
   !! Restricted Hartree-Fock over libcint integrals, in core.
   !!
   !! Canonical orthogonalisation, DIIS on the FDS - SDF commutator, energy
   !! from the density -- the same algorithm the cuEST SCF runs, written
   !! against host arrays rather than device pointers.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm, pic_gemv, pic_syrk
   use pic_lapack_interfaces, only: pic_syev
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_convergence_report, only: convergence_header, convergence_footer
   use mqc_timing, only: timing_report_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_scf_common, only: build_orthogonalizer, report_linear_dependence, &
                             build_density_closed_shell, &
                             build_density_spin, spin_contamination, GWH_K
   use mqc_diis, only: diis_state_t, ACCEL_DIIS, ACCEL_SWITCH
   use mqc_fock_projector, only: fock_projector_t
   ! `build_fock` lives in `mqc_czt_direct` and is re-exported below: this
   ! module cannot own it any more, because the shared response product
   ! needs it and this module uses `mqc_czt_hessian`, which uses that.
   use mqc_czt_direct, only: schwarz_bounds, build_fock_direct, build_fock_direct_uhf, &
                             direct_stats_t, build_fock
   use mqc_czt_integrals, only: czt_molecule_t, build_df_tensor, build_sap_potential, &
                                df_jfit_t, build_df_jfit, df_jfit_coulomb
   use mqc_czt_xc, only: xc_context_t, xc_add_potential, xc_add_potential_uks, xc_profile_take, &
                         xc_prune_small_rho
   use mqc_czt_cosx, only: cosx_grid_t, build_cosx_grid, cosx_exchange, cosx_exchange_trace, &
                           COSX_ROLE_SCF, COSX_ROLE_FINAL, COSX_ROLE_COARSE, &
                           COSX_COARSE_UNTIL, COSX_DEFAULT_GRID_LEVEL
   use mqc_czt_pcm, only: pcm_context_t
   use mqc_program_limits, only: DF_AUX_CHUNK
   use mqc_scf_convergence, only: scf_convergence_t, parse_convergence_metric
   use mqc_scf_types, only: scf_numerics_t
   use mqc_diis, only: parse_accelerator_name
   use mqc_calculation_defaults, only: DEFAULT_SOSCF_START, DEFAULT_GDM_START, &
                                       DEFAULT_GDM_MAX_ROTATION
   use mqc_scf_gdm, only: scf_gdm_history_t, scf_gdm_preconditioner_build, &
                          scf_gdm_step_judge, scf_gdm_step_limit, scf_gdm_is_descent
   use mqc_czt_hessian, only: nuclear_response_t
   use mqc_czt_ov_hessian, only: ov_hessian_t, build_scf_ov_hessian
   use mqc_czt_soscf, only: soscf_gradient, soscf_semicanonicalize, soscf_newton_step, &
                            uhf_rotation_hessian_t, soscf_newton_step_uhf, &
                            SPIN_HESSIAN_SCALE
   use mqc_orbital_rotation, only: rotation_matrix, MAX_ROTATION, MIN_ROTATION, &
                                   ENERGY_RESOLUTION, TRUST_GROWTH, SADDLE_CURVATURE
   implicit none
   private

   integer, parameter :: DEFAULT_UHF_DIIS_START = 4
      !! First iteration the unrestricted SCF is allowed to extrapolate on.
      !!
      !! Extrapolating from the start converges an open-shell system to the
      !! wrong state: a symmetric guess gives alpha and beta the same spatial
      !! orbitals, DIIS combines Fock matrices linearly, and a history that
      !! begins symmetric stays symmetric. The restricted path needs no such
      !! delay, having no symmetry to break.

   ! Which starting point the SCF is handed. Named here rather than in the
   ! atomic-guess module because that module calls this one -- the free-atom
   ! solutions an atomic guess is built from are themselves unrestricted SCFs
   ! -- and Fortran has no circular `use`.
   integer, parameter, public :: SCF_GUESS_CORE = 0   !! F = H
   integer, parameter, public :: SCF_GUESS_GWH = 1    !! Generalized Wolfsberg-Helmholz
   integer, parameter, public :: SCF_GUESS_SAC = 2    !! Superposed atomic spin densities, unaveraged
   integer, parameter, public :: SCF_GUESS_SAD = 3    !! Superposed atomic densities
   integer, parameter, public :: SCF_GUESS_PROJ = 4   !! Projected from a smaller basis
   integer, parameter, public :: SCF_GUESS_SAP = 5    !! Superposition of atomic potentials
   integer, parameter, public :: SCF_GUESS_MINAO = 6  !! Minimal ANO atoms projected in
   character(len=*), parameter, public :: SAP_BASIS_DEFAULT = "sap_helfem_large"
      !! Which fitted atomic potential `sap` uses unless a deck names another.
      !!
      !! The fully numerical HelFEM fit rather than a GRASP one: the relativistic
      !! sets exist for heavy elements and buy nothing below them, and this is
      !! the set Psi4 defaults to, which keeps a cross-code comparison of the
      !! *guess* about the guess.

   integer, parameter :: LINE_LEN = 160
      !! Buffer length for a formatted table line handed to the logger.

   public :: rhf_result_t
   public :: run_czt_rhf
   public :: run_czt_uhf
   public :: roothaan_fock             !! Roothaan's effective Fock matrix for a restricted open shell
   public :: set_df_jfit_enabled      !! Testing switch: false forces the dense fitted tensor for pure functionals
   public :: guess_fock                !! Starting Fock for the guesses needing no atomic SCF
   public :: build_fock                !! F = H + J - K/2; with H zero it is the response operator
   public :: density_energy            !! The closed-shell energy functional of a given density
   public :: assemble_fock_uhf         !! Both spin Fock matrices and their energy, no SCF
   public :: density_pseudo_orbitals   !! Factor a guess density for the fitted exchange

   logical, protected :: df_jfit_enabled = .true.
      !! Whether a density-fitted reference with no exact exchange is built from the
      !! Coulomb-only fit rather than the dense tensor. Only `set_df_jfit_enabled` moves it,
      !! and only a test has reason to.

   type :: incremental_state_t
      !! What an incremental Fock build carries between iterations
      !!
      !! G is linear in the density, so
      !!
      !!     G(D_n) = G(D_ref) + G(D_n - D_ref)
      !!
      !! and the correction is built from the density *change*. Since
      !! `build_fock_direct` screens on the largest density element a quartet
      !! multiplies, and the change shrinks as the SCF settles, the same
      !! screening test discards more with every iteration.
      !!
      !! Only the Coulomb and exchange part. The exchange-correlation potential
      !! is not linear in the density and is rebuilt in full every iteration;
      !! `g_ref` is captured before the long-range exchange and before V_xc are
      !! added, so neither is double counted.
      !!
      !! **Both passes of a range-separated build are incremental**, and both
      !! spins of an unrestricted one, each accumulating into its own
      !! reference. All of them share `since_reset` and `active`, so a full
      !! rebuild re-syncs every reference at once.
      !!
      !! **In an unrestricted SCF `d_ref`, `g_ref` and `g_ref_lr` are the alpha
      !! channel** and the `_b` fields are beta -- the plain names are not "the
      !! total" on that path.
      real(dp), allocatable :: d_ref(:, :)   !! Density `g_ref` belongs to (alpha, if unrestricted)
      real(dp), allocatable :: g_ref(:, :)   !! Accumulated G, without H (alpha, if unrestricted)
      real(dp), allocatable :: g_ref_lr(:, :)
         !! Accumulated long-range K, allocated only for a range-separated
         !! functional. Alpha channel on the unrestricted path.
      real(dp), allocatable :: d_ref_b(:, :)  !! Beta density, unrestricted only
      real(dp), allocatable :: g_ref_b(:, :)  !! Accumulated beta G, unrestricted only
      real(dp), allocatable :: g_ref_lr_b(:, :)
         !! Accumulated beta long-range K, unrestricted and range-separated only
      integer :: since_reset = 0             !! Iterations since the last full build
      integer :: full_builds = 0
         !! Full Fock builds performed, as opposed to updates from the density
         !! change. With incremental building on this is one plus a rebuild
         !! every `INCREMENTAL_RESET` iterations until the SCF settles, and one
         !! per iteration after; with it off, one per iteration.
      integer :: updates = 0                 !! Builds from the density change
      logical :: active = .false.            !! Off until the first full build seeds it
      logical :: settled = .false.
         !! The commutator has fallen below `INCREMENTAL_SETTLED`; every build
         !! from here on is a full one
      integer :: settled_at = 0              !! The iteration that happened at; zero if never
   end type incremental_state_t

   type :: rhf_result_t
      !! What a converged closed-shell SCF leaves behind
      real(dp) :: energy = 0.0_dp              !! Total, including nuclear repulsion
      real(dp) :: electronic = 0.0_dp          !! Without it
      real(dp) :: nuclear_repulsion = 0.0_dp
      integer :: iterations = 0
      integer :: full_fock_builds = 0
         !! Fock matrices built in full, as opposed to updated from the density
         !! change. Equals `iterations` when `incremental_fock` is off. When it
         !! is on: one, a rebuild every `INCREMENTAL_RESET` iterations, and one
         !! per iteration once the commutator is below `INCREMENTAL_SETTLED`.
      integer :: incremental_updates = 0
         !! Fock matrices built from the density change. Zero when
         !! `incremental_fock` is off.
      integer :: incremental_settled_at = 0
         !! The iteration whose commutator fell below `INCREMENTAL_SETTLED`, after
         !! which every build was full. Zero if it never did, or with
         !! `incremental_fock` off.
      ! The second-order phase, when `keywords.scf.second_order` asked for one.
      ! Zero everywhere means the whole SCF was DIIS, which is the default.
      integer :: second_order_started_at = 0
         !! The DIIS iteration after which the Newton steps took over. Reported
         !! because the switch is what makes the method robust and a reader has
         !! to be able to see that it happened, and where.
      integer :: second_order_iterations = 0
         !! Newton steps taken. `iterations` is the total of both phases.
      integer :: second_order_hessian_products = 0
         !! Hessian-vector products spent inside the Newton solves, each one a
         !! Fock build over the whole basis. **Not** counted in
         !! `full_fock_builds`, which counts the builds that produced an energy:
         !! these two together are what a second-order SCF actually costs.
      real(dp) :: second_order_curvature = 0.0_dp
         !! The smallest curvature the last Newton step saw, in hartree per
         !! radian squared. An **upper bound** on the true smallest eigenvalue
         !! of the orbital Hessian, because it is a Rayleigh-Ritz value from a
         !! Krylov subspace -- negative is proof of a saddle, positive is not
         !! proof of a minimum. `keywords.scf.stability` is what answers that.
      logical :: has_second_order_curvature = .false.
      logical :: second_order_stalled = .false.
         !! No step downhill was found. The energy reported is the last one
         !! that descended, and `converged` is false.
      ! The direct-minimisation phase, when `keywords.scf.gdm` asked for one.
      integer :: gdm_started_at = 0
         !! The DIIS iteration after which direct minimisation took over; zero
         !! if it never did.
      integer :: gdm_iterations = 0
         !! Trial steps taken, rejected ones included: one Fock build each.
      integer :: gdm_rejected = 0
         !! How many of those raised the energy and were taken back.
      logical :: gdm_stalled = .false.
         !! No step downhill was found along the quasi-Newton or the
         !! steepest-descent direction. The energy reported is the last one
         !! accepted, and `converged` is false.
      integer :: gdm_off_aufbau = 0
         !! Virtual orbitals that ended below the highest occupied one. A
         !! rotation never changes which orbitals are occupied, so a
         !! minimisation started too early can end on such a state; nonzero is
         !! warned about.
      logical :: converged = .false.
      real(dp) :: commutator = huge(1.0_dp)
         !! `max|FDS - SDF|` at the last iteration, in the orthogonal basis.
         !! The density's error scales with this, so a consumer of the density
         !! reads it rather than the threshold it was asked to stop at.
      real(dp), allocatable :: orbital_energies(:)
      real(dp), allocatable :: orbitals(:, :)
      real(dp), allocatable :: density(:, :)
      integer :: n_occupied = 0
      ! Unrestricted only. The restricted path leaves these unallocated and
      ! n_occupied_beta zero, so "was this open shell?" is a question the
      ! result can answer rather than something the caller has to remember.
      real(dp), allocatable :: orbital_energies_beta(:)
      real(dp), allocatable :: orbitals_beta(:, :)
      real(dp), allocatable :: density_beta(:, :)
      integer :: n_occupied_beta = 0
      real(dp) :: spin_squared = 0.0_dp        !! <S^2>, unrestricted only
      ! Continuum solvation. Zero when the SCF ran in the gas phase; when a
      ! context was passed, `energy` and `electronic` above already include
      ! `pcm_energy` -- this is the breakdown, not an addend.
      real(dp) :: pcm_energy = 0.0_dp          !! Dielectric energy, halved already
      real(dp) :: pcm_charge = 0.0_dp          !! Total apparent surface charge
   end type rhf_result_t

   ! ---------------------------------------------------------------------------
   ! The SCF iteration, as three things rather than forty arguments
   !
   ! `run_czt_rhf` used to carry its whole iteration inline: a 123-line loop body
   ! reading some forty locals. Lifting it out is what the routine wants anyway
   ! -- and nvfortran needs it, because its optimiser segfaults on a routine this
   ! size at -O2 and above, on every version from 25.11 to 26.5.
   !
   ! Extracting it with a plain argument list would have meant a 41-argument
   ! call, which trades one unreadable thing for another. These three say what
   ! each name is FOR, which is the part an argument list cannot:
   !
   !   rhf_operators_t   fixed for the whole SCF -- read, never written
   !   rhf_controls_t    the policy the deck chose -- cheap scalars
   !   rhf_state_t       what an iteration advances, and its scratch
   !
   ! Private, and RHF's: `run_czt_uhf` has its own loop and an alpha/beta state
   ! that is genuinely a different shape, so sharing these is a later change made
   ! once UHF has said what it needs rather than guessed at now.
   ! ---------------------------------------------------------------------------

   type :: rhf_operators_t
      !! What the iteration reads and never changes
      !!
      !! The arrays are MOVED in, not copied -- `eri` is n^4 and a copy of it
      !! would dwarf everything else this routine allocates.
      real(dp), allocatable :: h(:, :)        !! Core Hamiltonian, plus `h_extra` if given
      real(dp), allocatable :: s(:, :)        !! Overlap
      real(dp), allocatable :: x(:, :)        !! Orthogonaliser, (n_ao, n_mo)
      real(dp), allocatable :: bmat(:, :)     !! Fitted B(mu nu, P); unallocated if not fitted
      real(dp), allocatable :: bmat_lr(:, :)  !! The same against erf(omega r)/r
      type(df_jfit_t) :: jfit         !! Coulomb-only fit; `jfit%raw` unallocated unless one was built
      real(dp), allocatable :: eri(:, :, :, :)  !! In-core integrals; unallocated on the direct path
      real(dp), allocatable :: bounds(:, :)   !! Schwarz bounds for screening
      type(cosx_grid_t) :: cosx
         !! The grid the exchange is integrated on during the iterations;
         !! unbuilt unless the exchange is chain-of-spheres
      type(cosx_grid_t) :: cosx_coarse
         !! A smaller grid for the iterations before the commutator falls
         !! below `COSX_COARSE_UNTIL`. Unbuilt when the final step is off.
      type(cosx_grid_t) :: cosx_final
         !! The larger grid the final energy is taken on. Unbuilt when that
         !! step is off, and then `cosx` is the larger grid itself.
      integer :: n_ao = 0
      integer :: n_mo = 0                     !! May be below n_ao where the overlap was singular
      integer :: n_occ = 0
   end type rhf_operators_t

   type :: rhf_controls_t
      !! The policy the deck asked for, fixed once before the first iteration
      logical :: incremental = .true.   !! Build each Fock from the density change
      logical :: verbose = .false.
      logical :: kohn_sham = .false.    !! Only to label the iteration table
      logical :: use_pcm = .false.
      integer :: accel = 0              !! One of ACCEL_*
      real(dp) :: shift = 0.0_dp        !! Level shift, hartree; zero is off
      real(dp) :: taper = 0.0_dp        !! Below this dD(rms) the shift is dropped
      logical :: second_order = .false.  !! Hand over to Newton once close enough
      real(dp) :: soscf_start = 0.0_dp  !! The commutator at which that happens
      logical :: gdm = .false.          !! Hand over to direct minimisation instead
      real(dp) :: gdm_start = 0.0_dp    !! The commutator at which that happens
      real(dp) :: gdm_max_rotation = 0.0_dp
         !! Largest element of one step's rotation, radians
      type(scf_convergence_t) :: conv   !! What stopping is measured against
   end type rhf_controls_t

   type :: rhf_state_t
      !! What one iteration advances, and the scratch it advances it in
      !!
      !! The scratch lives here rather than being local to the iteration so that
      !! it is allocated once for the SCF instead of once per iteration.
      real(dp), allocatable :: density(:, :)
      real(dp), allocatable :: density_old(:, :)  !! Previous iteration's, for dD(rms)
      real(dp), allocatable :: coeff(:, :)
      real(dp), allocatable :: eigenvalues(:)
      real(dp), allocatable :: fock(:, :)
      real(dp), allocatable :: fock_flat(:)       !! F as a vector, which is what DIIS stores
      real(dp), allocatable :: err(:, :)          !! FDS - SDF in the orthogonal basis
      real(dp), allocatable :: sd(:, :)           !! S D and S D S, the level shift's
      real(dp), allocatable :: sds(:, :)          !!   virtual projector; allocated on first use
      real(dp), allocatable :: v_pcm(:, :)        !! Continuum operator, allocated only with PCM
      real(dp) :: e_elec = 0.0_dp
      real(dp) :: e_old = 0.0_dp                  !! Previous iteration's, for dE
      real(dp) :: drms_prev = huge(1.0_dp)        !! Drives the level-shift taper
      type(diis_state_t) :: diis
      type(incremental_state_t) :: incr
      type(direct_stats_t) :: screening
      logical :: cosx_coarse = .true.
         !! Whether a chain-of-spheres run is still on its coarse grid. Cleared
         !! once and never set again.
   end type rhf_state_t

   type :: uhf_state_t
      !! What one unrestricted iteration advances, and the scratch it uses
      !!
      !! `rhf_state_t`'s spin-resolved twin, and for the same two reasons: the
      !! scratch is allocated once for the SCF rather than once per iteration,
      !! and grouping it is what lets the iteration be a routine of its own.
      !! The operators and controls are not duplicated -- nothing in
      !! `rhf_operators_t` or `rhf_controls_t` depends on spin, and the
      !! unrestricted driver keeps them as its own locals rather than assemble
      !! a second copy.
      real(dp), allocatable :: d_a(:, :), d_b(:, :)
      real(dp), allocatable :: d_a_old(:, :), d_b_old(:, :)  !! Previous iteration's, for dD(rms)
      real(dp), allocatable :: coeff_a(:, :), coeff_b(:, :)
      real(dp), allocatable :: eig_a(:), eig_b(:)
      real(dp), allocatable :: fock_a(:, :), fock_b(:, :)
      real(dp), allocatable :: err_a(:, :), err_b(:, :)  !! FDS - SDF per spin
      real(dp), allocatable :: fock_flat(:)  !! Both spins' F end to end, which is what DIIS stores
      real(dp), allocatable :: err_flat(:)   !! Both spins' error, likewise
      real(dp), allocatable :: dens_flat(:)  !! Both spins' D, for the energy-based schemes
      real(dp), allocatable :: sd(:, :), sds(:, :)  !! The level shift's projector; on first use
      real(dp), allocatable :: v_pcm(:, :)   !! Continuum operator, allocated only with PCM
      real(dp), allocatable :: bmat(:, :)     !! Fitted B(mu nu, P); unallocated if not fitted
      real(dp), allocatable :: bmat_lr(:, :)  !! The same against erf(omega r)/r
      real(dp) :: e_elec = 0.0_dp
      real(dp) :: e_old = 0.0_dp             !! Previous iteration's, for dE
      real(dp) :: drms_prev = huge(1.0_dp)   !! Drives the level-shift taper
      logical :: rohf = .false.
         !! Restricted open shell: both spins share one set of orbitals, the
         !! eigenvectors of Roothaan's effective Fock matrix. `fock_a` then
         !! holds that matrix from the extrapolation on, and `fock_b` a copy.
      type(diis_state_t) :: diis
      type(incremental_state_t) :: incr
   end type uhf_state_t

   integer, parameter :: INCREMENTAL_RESET = 16
      !! Iterations between full rebuilds of the accumulated G.
      !!
      !! An incremental build adds a correction per iteration, so its rounding
      !! accumulates where a full build's does not; rebuilding every sixteenth
      !! iteration keeps the drift at the level of the convergence threshold
      !! for one full build in sixteen.

   ! The screening an incremental build drops accumulates between rebuilds:
   ! about 1e-10 in the energy on CAM-B3LYP water in cc-pVDZ. A deck asking for
   ! 1e-13 then converged only when a rebuild happened to land at the right
   ! iteration -- 10 to 88 iterations run to run, and past the limit on a slow
   ! runner. Below this commutator the density change is small enough that
   ! screening already skips most of a full build, so the few that remain cost
   ! little. Psi4's `INCFOCK_CONVERGENCE` makes the same switch at 1e-5.
   real(dp), parameter :: INCREMENTAL_SETTLED = 1.0e-5_dp
      !! Commutator below which every Fock build is a full one

   ! Stage labels, named once so the per-iteration column and the summary table
   ! cannot drift apart, and so a caller can ask `clk%seconds_of(STAGE_FOCK)`.
   character(len=*), parameter :: STAGE_SETUP = "setup (1e, bounds, guess)"
   character(len=*), parameter :: STAGE_FOCK = "Fock builds"
   character(len=*), parameter :: STAGE_DIAG = "diagonalisation"
   character(len=*), parameter :: STAGE_DIIS = "DIIS"
   character(len=*), parameter :: STAGE_XC = "XC quadrature"
   character(len=*), parameter :: STAGE_PCM = "continuum operator"
   character(len=*), parameter :: STAGE_SOSCF = "second-order step"
   character(len=*), parameter :: STAGE_GDM = "direct-minimisation step"

   integer, parameter :: GDM_HISTORY = 12
      !! L-BFGS pairs the direct minimisation keeps
   real(dp), parameter :: GDM_HESSIAN_FLOOR = 0.1_dp
      !! The smallest Hessian estimate `4 (e_a - e_i)` its preconditioner
      !! divides by, hartree per radian squared: a gap of 25 millihartree.
      !! Below it a near-degenerate pair would dominate the direction and the
      !! rotation cap would shrink every other component with it.
   real(dp), parameter :: GDM_SLACK = 1.0e-10_dp
      !! The rise in energy, hartree, a trial may show and still be accepted.
      !! Every energy it compares comes from a full Fock build, so this is the
      !! rounding of one and not a noise allowance.
   real(dp), parameter :: GDM_AUFBAU_TOL = 1.0e-6_dp
      !! A virtual orbital counts as below the highest occupied one only when
      !! it is lower by more than this, hartree.

contains

   subroutine xc_breakdown_summary(verbose)
      !! Where the XC quadrature spent its time, at `verbose`
      !!
      !! Wall seconds, averaged over the threads, so the total is the threaded
      !! part of the `XC quadrature` row above. Empty, and reset, when there
      !! was no quadrature.
      logical, intent(in) :: verbose

      character(len=*), parameter :: LABEL(4) = [character(len=20) :: "basis evaluation", "density", &
                                                                               "functional (libxc)", "potential assembly"]
      real(dp) :: total
      real(dp) :: sec(4)
      character(len=LINE_LEN) :: line
      integer :: i

      call xc_profile_take(sec)
      total = sum(sec)
      if (.not. verbose .or. total <= 0.0_dp) return

      call logger%verbose("")
      call logger%verbose("  XC quadrature breakdown")
      do i = 1, 4
         write (line, "(4x,a,f10.2,a,f7.1,a)") LABEL(i), sec(i), " s", 100.0_dp*sec(i)/total, " %"
         call logger%verbose(trim(line))
      end do
      write (line, "(4x,a,f10.2,a)") "total               ", total, " s"
      call logger%verbose(trim(line))
   end subroutine xc_breakdown_summary

   subroutine screening_summary(verbose, stats)
      !! How much of the quartet space each screening test removed
      !!
      !! From the last Fock build of the SCF, which is the converged one.
      !!
      !! Split by cause: the Schwarz count is a property of the basis and the
      !! geometry and does not move between iterations, while the density count
      !! is the extra reach that weighting the bound by the density buys, and
      !! is the only figure a change to the screening can move. Counts rather
      !! than seconds, because these are exactly reproducible.
      logical, intent(in) :: verbose
      type(direct_stats_t), intent(in) :: stats

      character(len=LINE_LEN) :: line

      if (.not. verbose) return
      if (stats%quartets_total <= 0) return
      call logger%large_info("")
      call logger%large_info("  screening (last Fock build)")
      write (line, "(a,i14)") "    unique quartets      ", stats%quartets_total
      call logger%large_info(trim(line))
      write (line, "(a,i14,f9.1,a)") "    skipped, Schwarz     ", stats%screened_schwarz, &
         100.0_dp*real(stats%screened_schwarz, dp)/real(stats%quartets_total, dp), " %"
      call logger%large_info(trim(line))
      write (line, "(a,i14,f9.1,a)") "    skipped, density     ", stats%screened_density, &
         100.0_dp*real(stats%screened_density, dp)/real(stats%quartets_total, dp), " %"
      call logger%large_info(trim(line))
      write (line, "(a,i14,f9.1,a)") "    computed             ", stats%quartets_computed, &
         100.0_dp*real(stats%quartets_computed, dp)/real(stats%quartets_total, dp), " %"
      call logger%large_info(trim(line))
      ! One is perfect balance. A ratio inside one run, so a contended node
      ! largely divides out of it.
      write (line, "(a,f14.4)") "    thread imbalance     ", stats%thread_imbalance
      call logger%large_info(trim(line))
   end subroutine screening_summary

   pure function scheme_now(accel, err_norm) result(scheme)
      !! Which accelerator runs this iteration
      !!
      !! An energy-based scheme runs while the error norm is above
      !! `ACCEL_SWITCH` and hands over to DIIS below it -- on the norm rather
      !! than an iteration count, so an SCF that starts near its solution never
      !! pays for the interpolation. A caller that asked for DIIS gets DIIS at
      !! every norm.
      integer, intent(in) :: accel
      real(dp), intent(in) :: err_norm
      integer :: scheme

      scheme = ACCEL_DIIS
      if (accel /= ACCEL_DIIS .and. err_norm > ACCEL_SWITCH) scheme = accel
   end function scheme_now

   subroutine scf_table_header(verbose, kohn_sham)
      !! Column headings for the per-iteration table
      !!
      !! The frame is `mqc_convergence_report`, shared with the FMO outer loop;
      !! only the columns and their widths are the SCF's own.
      logical, intent(in) :: verbose
      logical, intent(in) :: kohn_sham
         !! Whether a functional is being evaluated. A Hartree-Fock run has no
         !! quadrature, so its `XC` column is dropped rather than printed as a
         !! column of zeros.

      if (kohn_sham) then
         call convergence_header(verbose, "SCF iterations", &
                                 "    iter                 energy          dE"// &
                                 "        diis      n       Fock         XC       rest", 93)
      else
         call convergence_header(verbose, "SCF iterations", &
                                 "    iter                 energy          dE"// &
                                 "        diis      n       Fock       rest", 82)
      end if
   end subroutine scf_table_header

   subroutine scf_table_row(verbose, iter, energy, de, gnorm, ndiis, t_fock, t_xc, &
                            t_rest, kohn_sham)
      !! One iteration's line, with the time that iteration took
      !!
      !! `STAGE_FOCK` and `STAGE_XC` are separate buckets and both are printed:
      !! on a Kohn-Sham run the quadrature is much the larger of the two, so a
      !! row carrying only the Fock time cannot be reconciled against a wall
      !! clock.
      !!
      !! Per-iteration rather than only a total, because screening tightens as
      !! the density settles: a Fock build that is flat across iterations says
      !! the density-weighted screening is not reaching.
      logical, intent(in) :: verbose
      integer, intent(in) :: iter, ndiis
      real(dp), intent(in) :: energy, de, gnorm, t_fock, t_xc, t_rest
      logical, intent(in) :: kohn_sham   !! Print the quadrature column at all

      character(len=LINE_LEN) :: line

      if (.not. verbose) return
      if (kohn_sham) then
         write (line, "(i8,f23.12,2es12.3,i7,3(f9.2,a))") &
            iter, energy, de, gnorm, ndiis, t_fock, " s", t_xc, " s", t_rest, " s"
      else
         write (line, "(i8,f23.12,2es12.3,i7,2(f9.2,a))") &
            iter, energy, de, gnorm, ndiis, t_fock, " s", t_rest, " s"
      end if
      call logger%info(trim(line))
   end subroutine scf_table_row

   subroutine scf_table_footer(verbose, converged, iterations)
      !! Close the per-iteration table
      !!
      !! The timing summary is `clk%report(...)`, which counts calls itself: the
      !! Fock stage runs `iterations + 1` times, because the converged energy is
      !! rebuilt after the loop.
      logical, intent(in) :: verbose, converged
      integer, intent(in) :: iterations

      call convergence_footer(verbose, converged, iterations, "iterations", 84)
   end subroutine scf_table_footer

   subroutine energy_components(verbose, mol, density, electronic, nuclear)
      !! The energy split into the pieces the virial theorem constrains
      !!
      !! Traces against the converged density, plus one kinetic and one
      !! nuclear-attraction build.
      !!
      !! For a variational wave function at a stationary geometry `-V/T` is
      !! exactly 2. It is not constrained to be, so drift from 2 is evidence
      !! about the calculation -- a basis too small for the cusp, or a geometry
      !! that is not a stationary point.
      logical, intent(in) :: verbose
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(in) :: electronic   !! Total electronic energy
      real(dp), intent(in) :: nuclear      !! Nuclear repulsion

      real(dp), allocatable :: kinetic(:, :), core(:, :)
      real(dp) :: t_energy, one_e, two_e, v_ne, v_ee, v_total
      character(len=128) :: line

      if (.not. verbose) return

      call mol%kinetic(kinetic)
      call mol%core_hamiltonian(core)
      t_energy = sum(density*kinetic)
      one_e = sum(density*core)
      ! The core Hamiltonian is kinetic plus nuclear attraction, so the
      ! attraction alone is the difference.
      v_ne = one_e - t_energy
      two_e = electronic - one_e
      v_ee = two_e
      v_total = v_ee + v_ne + nuclear

      call logger%large_info("")
      call logger%large_info("  energy components")
      write (line, "(a,f22.10)") "    one electron                ", one_e
      call logger%large_info(trim(line))
      write (line, "(a,f22.10)") "    two electron                ", two_e
      call logger%large_info(trim(line))
      write (line, "(a,f22.10)") "    nuclear repulsion           ", nuclear
      call logger%large_info(trim(line))
      write (line, "(a,f22.10)") "    total                       ", electronic + nuclear
      call logger%large_info(trim(line))
      write (line, "(a,f22.10)") "    electron-electron potential ", v_ee
      call logger%large_info(trim(line))
      write (line, "(a,f22.10)") "    nucleus-electron potential  ", v_ne
      call logger%large_info(trim(line))
      write (line, "(a,f22.10)") "    nucleus-nucleus potential   ", nuclear
      call logger%large_info(trim(line))
      write (line, "(a,f22.10)") "    total potential             ", v_total
      call logger%large_info(trim(line))
      write (line, "(a,f22.10)") "    total kinetic               ", t_energy
      call logger%large_info(trim(line))
      if (abs(t_energy) > 0.0_dp) then
         write (line, "(a,f22.10)") "    virial ratio (-V/T)         ", -v_total/t_energy
         call logger%large_info(trim(line))
      end if
      deallocate (kinetic, core)
   end subroutine energy_components

   subroutine set_df_jfit_enabled(enabled)
      !! Force (false) or release (true) the dense fitted tensor for pure functionals
      logical, intent(in) :: enabled

      df_jfit_enabled = enabled
   end subroutine set_df_jfit_enabled

   subroutine run_czt_rhf(mol, nelec, max_iter, energy_tol, density_tol, &
                          verbose, result, error, aux, diis_vectors, in_core, &
                          guess, guess_density, xc, h_extra, pcm, projector, &
                          level_shift, linear_dependence, b_ao_out, accelerator, grad_tol, &
                          incremental_fock, convergence, scf, second_order, soscf_start, &
                          cosx, cosx_grid_level, cosx_final_grid, gdm, gdm_start, &
                          gdm_max_rotation)
      !! Drive a closed-shell SCF to convergence
      type(czt_molecule_t), intent(in) :: mol
      integer, intent(in) :: nelec
      integer, intent(in) :: max_iter
      real(dp), intent(in) :: energy_tol, density_tol
      real(dp), intent(in), optional :: grad_tol
         !! Commutator threshold -- pyscf's `conv_tol_grad`. Zero or absent
         !! derives `sqrt(energy_tol)`, which is the right scaling *for the
         !! energy*: its error goes as `gnorm**2`.
         !!
         !! **A caller that consumes the density has to set this.** Density
         !! error goes as `gnorm`, not its square, so the derived value bounds
         !! it about three orders more loosely than the energy, and tightening
         !! `energy_tol` cannot substitute: reaching `gnorm < 1e-8` that way
         !! needs `energy_tol = 1e-16`, below what a molecular energy resolves.
      type(scf_convergence_t), intent(in), optional :: convergence
         !! The rule that decides this SCF has converged. Absent reproduces
         !! what `energy_tol` and `grad_tol` meant before this existed.
      type(scf_numerics_t), intent(in), optional :: scf
         !! The whole configuration, as one argument.
         !!
         !! **Eight of its thirteen fields are honoured here**, and passing
         !! this looks like passing everything. Read: `level_shift`,
         !! `linear_dependence`, `accelerator`, `diis_size`, `use_diis`,
         !! `incremental_fock`, `convergence_metric` and `grad_tol`. Not read:
         !!
         !! * `max_iter`, `energy_tol` and `density_tol`, which are positional
         !!   arguments here -- whichever the caller passes positionally runs.
         !! * `guess`, a deck spelling whose translation to an `SCF_GUESS_*`
         !!   kind needs `mqc_czt_atomic_guess`, which calls this routine.
         !!   The caller parses it and passes `guess=`.
         !! * `allow_crap_scf`. This routine reports `result%converged` and the
         !!   caller decides whether an unconverged SCF is acceptable.
         !!
         !! An individual optional above wins over the group's copy of the same
         !! setting, so a caller can pass the group and override one field.
      logical, intent(in), optional :: incremental_fock
         !! Build each iteration from the density change. Default true; false
         !! forces a full build every iteration, which is what to reach for when
         !! an SCF stalls or when a Fock timing needs to mean something.
      logical, intent(in) :: verbose
      real(dp), intent(in), optional :: level_shift
         !! Hartree added to the virtual block of the Fock matrix before each
         !! diagonalisation. Absent or zero is off.
      integer, intent(in), optional :: accelerator
         !! `ACCEL_DIIS`, the default, or one of the energy-based schemes to run
         !! while the error is large.
      type(rhf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: b_ao_out(:, :)
         !! The fitted `B(mu nu, P)` this SCF used, moved out rather than
         !! freed, for a correlated step that would otherwise rebuild the
         !! integrals and refit them. Only meaningful alongside `aux`.
      type(czt_molecule_t), intent(in), optional :: aux
         !! Auxiliary basis. Present means density-fitted J and K; absent means
         !! exact ERIs.
      integer, intent(in), optional :: diis_vectors
         !! Subspace size. Zero turns DIIS off.
      logical, intent(in), optional :: in_core
         !! Store every integral and contract from the tensor, instead of
         !! rebuilding the Fock matrix directly each iteration. Default is
         !! direct. This is the path validated against PySCF, and so what the
         !! direct build is checked against; it is n^4 in memory.
      integer, intent(in), optional :: guess
         !! One of SCF_GUESS_*. Defaults to Wolfsberg-Helmholz, the best guess
         !! that needs nothing but `H` and `S`. Superposition of atomic
         !! densities is better still and is what a deck gets, but
         !! `guess_density` has to be built before this routine is entered and
         !! so cannot be a default here.
      type(xc_context_t), intent(inout), optional :: xc
         !! Exchange-correlation. Present turns this into a Kohn-Sham SCF;
         !! absent leaves it Hartree-Fock.
      real(dp), intent(in), optional :: guess_density(:, :)
         !! Total starting density, required by SCF_GUESS_SAC, SCF_GUESS_SAD
         !! and SCF_GUESS_MINAO and ignored otherwise. Built by `mqc_czt_atomic_guess`, which
         !! cannot be reached from here: it runs free-atom SCFs through this
         !! very module.
      real(dp), intent(in), optional :: h_extra(:, :)
         !! An extra one-electron operator, added to H before the first Fock
         !! build and kept there; a uniform electric field `F . r` is what it is
         !! for, making this SCF the finite-field reference `..._cphf` is
         !! checked against. `result%energy` then includes the interaction with
         !! it, so the finite-field check differentiates the *dipole* and not
         !! the energy.
      type(pcm_context_t), intent(inout), optional :: pcm
         !! Continuum solvation. Present and built means every iteration solves
         !! the surface charges against its density and folds their operator
         !! into the Fock matrix; `result%energy` then includes the dielectric
         !! energy. Present-but-disabled is the same as absent, so a caller can
         !! pass its context unconditionally.
      real(dp), intent(in), optional :: linear_dependence
         !! `keywords.scf.linear_dependence_threshold`. Zero or absent leaves
         !! the orthogonaliser on its own cutoff.
      logical, intent(in), optional :: second_order
         !! `keywords.scf.second_order`. Hand over from DIIS to trust-region
         !! Newton on the orbital rotations once the commutator falls below
         !! `soscf_start`. Default false.
         !!
         !! **Not a way to start second order.** A Newton step from a poor
         !! density diverges, so the DIIS phase is not optional and the switch
         !! threshold is what makes the method robust; see `run_soscf_phase`.
      logical, intent(in), optional :: cosx
         !! `keywords.scf.exchange: "cosx"`. Build the exchange matrix by
         !! chain-of-spheres quadrature and the Coulomb matrix from the fit, so
         !! `aux` must be present. Default false.
      integer, intent(in), optional :: cosx_grid_level
         !! `keywords.scf.cosx_grid`, 1 to 3: the size of the grids that
         !! quadrature runs on. Absent takes `COSX_DEFAULT_GRID_LEVEL`. Read
         !! only with `cosx`.
      logical, intent(in), optional :: cosx_final_grid
         !! `keywords.scf.cosx_final_grid`. True, the default: iterate on the
         !! level's small grid and take the final energy on its larger one.
         !! False: iterate on the larger grid throughout, which costs more and
         !! makes the energy self-consistent on that grid. Read only with
         !! `cosx`.
      real(dp), intent(in), optional :: soscf_start
         !! The commutator at which that switch happens, as
         !! `max|FDS - SDF|` in the orthogonal basis. Absent or non-positive
         !! takes `DEFAULT_SOSCF_START`.
      logical, intent(in), optional :: gdm
         !! `keywords.scf.gdm`. Hand over from DIIS to geometric direct
         !! minimisation once the commutator falls below `gdm_start`; see
         !! `run_gdm_phase`. Default false, and refused together with
         !! `second_order`.
      real(dp), intent(in), optional :: gdm_start
         !! The commutator at which that switch happens. Absent or
         !! non-positive takes `DEFAULT_GDM_START`.
      real(dp), intent(in), optional :: gdm_max_rotation
         !! Largest element of one step's rotation, radians. Absent or
         !! non-positive takes `DEFAULT_GDM_MAX_ROTATION`.
      type(fock_projector_t), intent(in), optional :: projector
         !! A constraint on the Fock matrix, applied after every build.
         !!
         !! Unlike `h_extra` this is not an operator added to `H` but a linear
         !! map on `F` that forces it block diagonal in a basis the caller
         !! chose, so it cannot be folded into the core Hamiltonian and is
         !! reapplied each iteration. Freezing orbitals is what it is for -- see
         !! `mqc_fock_projector`.
         !!
         !! Applied before the DIIS error and extrapolation, so everything
         !! downstream sees the constrained matrix. **Not** applied to the final
         !! build: the constraint decides which determinant comes out, and the
         !! energy of that determinant is an expectation value of the real
         !! Hamiltonian.

      character(len=LINE_LEN) :: line
      integer :: diis_size, guess_kind
      logical :: use_in_core, use_cosx, use_final_grid, wants_jfit
      real(dp), allocatable :: cosx_bounds(:, :), fock_energy(:, :)
      integer :: cosx_level
      real(dp), allocatable :: bounds(:, :)
      ! TODO(mqc): `stats` is declared here and never used; the screening counts
      ! come back in `screening`, which `assemble_fock` fills.
      type(direct_stats_t) :: stats

      type(rhf_operators_t) :: ops
      type(rhf_controls_t) :: ctrl
      type(rhf_state_t) :: st
      real(dp) :: e_pcm
      real(dp) :: gtol
      logical :: accel_ok_grp, accel_second_order, accel_gdm
      integer :: metric_kind
      logical :: metric_ok
      integer :: iter
      integer :: soscf_switch   !! DIIS iteration the Newton phase took over after
      integer :: gdm_switch     !! and the one direct minimisation took over after
      type(timing_report_t) :: clk
      real(dp) :: xc_discard(4)
      real(dp) :: s_min, s_kept   !! overlap conditioning, reported before iteration 1
      real(dp) :: lindep

      if (mod(nelec, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "RHF needs an even electron count; this "// &
                        "system has an odd one and wants an unrestricted method")
         return
      end if

      diis_size = 8
      if (present(scf)) then
         diis_size = scf%diis_size
         if (.not. scf%use_diis) diis_size = 0
      end if
      if (present(diis_vectors)) diis_size = diis_vectors
      use_in_core = .false.
      ctrl%incremental = .true.
      if (present(scf)) ctrl%incremental = scf%incremental_fock
      if (present(incremental_fock)) ctrl%incremental = incremental_fock
      if (present(in_core)) use_in_core = in_core
      use_cosx = .false.
      if (present(cosx)) use_cosx = cosx
      if (use_cosx) then
         ! The Coulomb half of the method is the fitted one, so there is no
         ! chain-of-spheres run without a fit.
         if (.not. present(aux)) then
            call error%set(ERROR_VALIDATION, "chain-of-spheres exchange pairs with a "// &
                           "density-fitted Coulomb matrix and no auxiliary basis was "// &
                           "given. Set keywords.scf.density_fitting and model.aux_basis.")
            return
         end if
         ! The attenuated potential integral on the grid is not built, and
         ! falling back to the fitted long-range exchange would make the run
         ! half seminumerical without saying so.
         if (present(xc)) then
            if (xc%active .and. xc%range_separated) then
               call error%set(ERROR_VALIDATION, "chain-of-spheres exchange with a "// &
                              "range-separated functional is not implemented: the "// &
                              "long-range exchange needs the attenuated potential "// &
                              "integral on the grid. Set keywords.scf.exchange to "// &
                              "exact for this functional.")
               return
            end if
         end if
      end if
      guess_kind = SCF_GUESS_GWH
      if (present(guess)) guess_kind = guess
      ctrl%use_pcm = .false.
      if (present(pcm)) ctrl%use_pcm = pcm%enabled

      ops%n_ao = mol%nao
      ops%n_occ = nelec/2
      if (ops%n_occ < 1) then
         call error%set(ERROR_VALIDATION, "RHF: no electrons to place")
         return
      end if

      ! The count and the angular form together are what name the basis: a
      ! set read spherical where it should be Cartesian converges just as
      ! prettily, one function short. Gated on `verbose` rather than on the
      ! logger level, because the atomic guess runs one of these per element.
      if (verbose) then
         if (mol%cartesian) then
            call logger%info("  basis functions: "//to_char(ops%n_ao)//"  (Cartesian, 6d/10f)")
         else
            call logger%info("  basis functions: "//to_char(ops%n_ao)//"  (spherical, 5d/7f)")
         end if
         ! Charge derived rather than passed: `mol%charges` carries the ECP's
         ! core subtracted and `nelec` arrives reduced by the same cores, so the
         ! difference is the physical charge either way.
         call logger%info("  restricted: charge = "// &
                          to_char(nint(sum(mol%charges)) - nelec)// &
                          "  multiplicity = 1"// &
                          "  electrons = "//to_char(nelec))
      end if

      call clk%start()
      ! The quadrature accumulator is module-wide; anything an earlier caller
      ! left in it is not this SCF's time.
      call xc_profile_take(xc_discard)
      call mol%overlap(ops%s)
      call mol%core_hamiltonian(ops%h)
      if (present(h_extra)) then
         if (size(h_extra, 1) /= ops%n_ao .or. size(h_extra, 2) /= ops%n_ao) then
            call error%set(ERROR_VALIDATION, "RHF: h_extra is not n_ao square")
            return
         end if
         ops%h = ops%h + h_extra
      end if
      if (present(aux)) then
         ! A reference that takes no exchange from the fit only ever contracts
         ! B into J, so it keeps the raw integrals and the metric factor and
         ! skips the `(mn|Q) J^-1/2` transform. That is a functional with no
         ! exact exchange, and a chain-of-spheres run, whose exchange comes
         ! from the grid. Whoever asked for B back (`b_ao_out`) gets the dense
         ! tensor from the start rather than a second build at the end; a
         ! chain-of-spheres run has no use for B, so it does not wait on that.
         wants_jfit = .false.
         if (use_cosx) then
            wants_jfit = df_jfit_enabled
         else if (df_jfit_enabled .and. .not. present(b_ao_out)) then
            if (present(xc)) wants_jfit = xc%active .and. xc%exx_fraction == 0.0_dp &
                                          .and. .not. xc%range_separated
         end if
         if (wants_jfit) then
            call build_df_jfit(mol, aux, ops%jfit, error)
            if (error%has_error()) return
         else
            call build_df_tensor(mol, aux, ops%bmat, error)
         end if
         if (error%has_error()) return
         ! A second fit, against the attenuated kernel, when the functional
         ! splits its exchange by range. Built once beside the first rather than
         ! per iteration: it depends on the geometry and the auxiliary basis,
         ! neither of which moves during an SCF.
         if (present(xc)) then
            if (xc%range_separated) then
               call build_df_tensor(mol, aux, ops%bmat_lr, error, omega=xc%rs_omega)
            end if
         end if
         if (error%has_error()) return
      end if
      if (present(aux)) then
         if (use_cosx) then
            ! One set of Schwarz bounds serves every exchange grid.
            call schwarz_bounds(mol, cosx_bounds, error)
            if (error%has_error()) return
            cosx_level = COSX_DEFAULT_GRID_LEVEL
            if (present(cosx_grid_level)) cosx_level = cosx_grid_level
            use_final_grid = .true.
            if (present(cosx_final_grid)) use_final_grid = cosx_final_grid
            ! Both grids carry the overlap fit, which is what lets them be
            ! this small. The iteration grid is too coarse radially to take an
            ! energy from, so with the final step off the iterations run on
            ! the larger grid from the start.
            if (use_final_grid) then
               call build_cosx_grid(mol, ops%cosx, error, level=cosx_level, &
                                    role=COSX_ROLE_SCF, overlap=ops%s, bounds=cosx_bounds)
               if (error%has_error()) return
               call build_cosx_grid(mol, ops%cosx_final, error, level=cosx_level, &
                                    role=COSX_ROLE_FINAL, overlap=ops%s, bounds=cosx_bounds)
               if (error%has_error()) return
               call build_cosx_grid(mol, ops%cosx_coarse, error, level=cosx_level, &
                                    role=COSX_ROLE_COARSE, overlap=ops%s, bounds=cosx_bounds)
            else
               call build_cosx_grid(mol, ops%cosx, error, level=cosx_level, &
                                    role=COSX_ROLE_FINAL, overlap=ops%s, bounds=cosx_bounds)
            end if
            if (error%has_error()) return
            if (verbose) then
               call logger%info("  exchange: chain of spheres, grid level "// &
                                to_char(cosx_level)//", overlap fitted; Coulomb fitted")
               call logger%info("    iterations on "//to_char(ops%cosx%n_points)// &
                                " points")
               if (use_final_grid) then
                  call logger%info("    final energy on "// &
                                   to_char(ops%cosx_final%n_points)//" points")
               end if
            end if
         end if
      else if (use_in_core) then
         call mol%eris(ops%eri)
      else
         ! The bounds depend on the basis and the geometry, not the density, so
         ! one set serves every iteration of the SCF.
         call schwarz_bounds(mol, ops%bounds, error)
         if (error%has_error()) return
      end if

      lindep = 0.0_dp
      if (present(scf)) lindep = scf%linear_dependence
      if (present(linear_dependence)) lindep = linear_dependence
      call build_orthogonalizer(ops%s, ops%x, ops%n_mo, error, smallest_overlap=s_min, &
                                smallest_kept=s_kept, threshold=lindep)
      if (error%has_error()) return
      call report_linear_dependence(ops%n_ao, ops%n_mo, s_min, s_kept, verbose, threshold=lindep)
      if (ops%n_occ > ops%n_mo) then
         call error%set(ERROR_VALIDATION, "RHF: more occupied orbitals than the basis "// &
                        "supports after near-null modes were dropped")
         return
      end if

      allocate (st%fock(ops%n_ao, ops%n_ao), st%density(ops%n_ao, ops%n_ao), st%density_old(ops%n_ao, ops%n_ao))
      allocate (st%err(ops%n_mo, ops%n_mo), st%fock_flat(ops%n_ao*ops%n_ao))
      if (ctrl%use_pcm) allocate (st%v_pcm(ops%n_ao, ops%n_ao))

      ! The error vector lives in the orthogonal basis, where it is n_mo square
      ! rather than n_ao -- the same shape the cuEST path uses.
      ! Resolved before the subspace is built, not after it. `energy_based`
      ! decides whether the history can hold densities and energies at all, and
      ! reading `accel` here while it was still assigned some fifty lines below
      ! took that decision on an undefined value.
      ctrl%accel = ACCEL_DIIS
      accel_second_order = .false.
      accel_gdm = .false.
      if (present(scf)) call parse_accelerator_name(scf%accelerator, ctrl%accel, &
                                                    accel_ok_grp, accel_second_order, &
                                                    accel_gdm)
      if (present(accelerator)) ctrl%accel = accelerator
      call st%diis%init(diis_size, ops%n_ao*ops%n_ao, ops%n_mo*ops%n_mo, &
                        energy_based=(ctrl%accel /= ACCEL_DIIS))

      ! Every guess ends up as a starting Fock matrix, which is then diagonalised
      ! and occupied exactly as an iteration's would be. That uniformity is the
      ! point: the atomic guesses differ from the core guess only in what F is
      ! built from, so nothing below this line knows which one ran.
      select case (guess_kind)
      case (SCF_GUESS_CORE)
         st%fock = ops%h
      case (SCF_GUESS_GWH)
         call guess_fock(ops%s, ops%h, st%fock)
      case (SCF_GUESS_SAP)
         call sap_fock(mol, ops%h, st%fock, error)
         if (error%has_error()) return
      case (SCF_GUESS_SAC, SCF_GUESS_SAD, SCF_GUESS_PROJ, SCF_GUESS_MINAO)
         if (.not. present(guess_density)) then
            call error%set(ERROR_VALIDATION, "RHF: an atomic guess was asked for but no "// &
                           "guess density was supplied")
            return
         end if
         if (size(guess_density, 1) /= ops%n_ao .or. size(guess_density, 2) /= ops%n_ao) then
            call error%set(ERROR_VALIDATION, "RHF: the guess density is not the size of "// &
                           "this basis")
            return
         end if
         st%density = guess_density
         if (ops%cosx_coarse%active()) then
            call atomic_guess_fock(mol, ops%h, st%density, ops%bmat, ops%eri, ops%bounds, st%fock, &
                                   error, jfit=ops%jfit, xc=xc, cosx=ops%cosx_coarse)
         else
            call atomic_guess_fock(mol, ops%h, st%density, ops%bmat, ops%eri, ops%bounds, st%fock, &
                                   error, jfit=ops%jfit, xc=xc, cosx=ops%cosx)
         end if
         if (error%has_error()) return
      case default
         call error%set(ERROR_VALIDATION, "RHF: unknown initial guess")
         return
      end select

      call diagonalize(st%fock, ops%x, ops%n_ao, ops%n_mo, st%coeff, st%eigenvalues, error)
      if (error%has_error()) return
      call build_density_closed_shell(st%coeff, ops%n_occ, st%density)

      st%e_old = 0.0_dp
      result%converged = .false.

      ! Everything from the top of the routine to here is one-time cost: the 1e
      ! integrals, the screening bounds or the fitted/in-core tensor, the
      ! orthogonaliser and the guess.
      call clk%lap(STAGE_SETUP)
      ! `present(xc)` is not the question: the bridge passes its context
      ! unconditionally, so it is present on a Hartree-Fock run too. `active`
      ! is what says a functional was constructed.
      ctrl%kohn_sham = .false.
      if (present(xc)) ctrl%kohn_sham = xc%active
      ! Once per context, on the guess density, before the first quadrature:
      ! the atomic guess when there is one, else the density the guess Fock
      ! matrix just gave. A no-op for a cutoff of zero.
      if (ctrl%kohn_sham) then
         select case (guess_kind)
         case (SCF_GUESS_SAC, SCF_GUESS_SAD, SCF_GUESS_PROJ, SCF_GUESS_MINAO)
            call xc_prune_small_rho(xc, mol, guess_density, nelec, error)
         case default
            call xc_prune_small_rho(xc, mol, st%density, nelec, error)
         end select
         if (error%has_error()) return
      end if
      call scf_table_header(verbose, ctrl%kohn_sham)

      ! pyscf's `conv_tol_grad`: derived from the energy tolerance unless a
      ! caller states it. See the argument's own note for why a density
      ! consumer must.
      ! Assembled once. A caller that states nothing gets `CONV_METRIC_STANDARD`
      ! with the bound `sqrt(energy_tol)`.
      if (present(convergence)) then
         ctrl%conv = convergence
      else
         ctrl%conv%tolerance = energy_tol
         if (present(scf)) ctrl%conv%gradient_tolerance = scf%grad_tol
         if (present(grad_tol)) ctrl%conv%gradient_tolerance = grad_tol
         ! The metric named by the group, when one was passed. A spelling the
         ! caller has already validated; an unparseable one leaves the default
         ! rather than failing a calculation from inside the SCF.
         if (present(scf)) then
            call parse_convergence_metric(scf%convergence_metric, metric_kind, metric_ok)
            if (metric_ok) ctrl%conv%metric = metric_kind
         end if
      end if
      ! TODO(mqc): `gtol` is assigned here and never read -- `conv%is_converged`
      ! applies the commutator bound itself. Dead on both spin paths.
      gtol = ctrl%conv%commutator_bound()

      ! The group first, then any individual override on top of it.
      ctrl%shift = 0.0_dp
      if (present(scf)) ctrl%shift = scf%level_shift
      if (present(level_shift)) ctrl%shift = level_shift
      ! Refused rather than clamped: a negative shift lowers the virtuals into
      ! the occupied set, narrowing the gap the next density is built through
      ! and driving the oscillation it was asked to damp.
      if (ctrl%shift < 0.0_dp) then
         call error%set(ERROR_VALIDATION, "keywords.scf.level_shift is negative. A "// &
                        "level shift raises the virtual orbitals to widen the gap; a "// &
                        "negative one narrows it and makes convergence worse, not "// &
                        "better. Give a positive value in Hartree, or leave it out.")
         return
      end if
      ! Tapered off well before the SCF is done, because every orbital energy
      ! this routine reports has to belong to the *unshifted* operator: they are
      ! read back as MP2 and coupled-cluster denominators, as the weights of the
      ! gradient's energy-weighted density, and as `eps_occ` and the response
      ! poles of a fragment potential.
      ctrl%taper = 100.0_dp*density_tol
      st%drms_prev = huge(1.0_dp)

      ! Reported with its taper: a shift that is off by iteration three is not
      ! the shift the reader thinks they applied, and a shift dropped on the way
      ! here is otherwise indistinguishable from one that was never set.
      if (ctrl%shift > 0.0_dp) then
         write (line, "(a,f8.4,a,es9.2)") "    level shift: ", ctrl%shift, &
            " hartree, tapered off below dD ", ctrl%taper
         call logger%info(trim(line))
      end if

      ! `verbose` is a dummy of this routine, so it is the one control that is
      ! copied across rather than having been written straight into `ctrl`.
      ctrl%verbose = verbose

      ! ---- second order, and what it will not do ---------------------------
      !
      ! Refused rather than approximated. The Newton step's Hessian is the
      ! two-electron (and exchange-correlation) response of the *gas-phase*
      ! reference with no constraint on it, so neither of these would be the
      ! curvature of the energy actually being minimised, and a step taken on
      ! the wrong curvature converges to the wrong place or not at all.
      ! `accelerator: soscf` and `second_order: true` are the same request, so
      ! either one turns it on and neither can turn the other off.
      ctrl%second_order = accel_second_order
      if (present(scf)) ctrl%second_order = ctrl%second_order .or. scf%second_order
      if (present(second_order)) ctrl%second_order = second_order .or. accel_second_order
      ctrl%soscf_start = DEFAULT_SOSCF_START
      if (present(scf)) then
         if (scf%soscf_start > 0.0_dp) ctrl%soscf_start = scf%soscf_start
      end if
      if (present(soscf_start)) then
         if (soscf_start > 0.0_dp) ctrl%soscf_start = soscf_start
      end if
      if (ctrl%second_order) then
         if (ops%cosx%active()) then
            call error%set(ERROR_VALIDATION, "keywords.scf.second_order with "// &
                           "chain-of-spheres exchange: the orbital-rotation Hessian "// &
                           "is built from the exact exchange, which is not the "// &
                           "curvature of the energy this SCF minimises. Converge "// &
                           "this one by DIIS.")
            return
         end if
         if (ctrl%use_pcm) then
            call error%set(ERROR_VALIDATION, "keywords.scf.second_order with a "// &
                           "continuum solvent: the orbital-rotation Hessian carries "// &
                           "no response of the surface charges, so the Newton step "// &
                           "would be taken on the wrong curvature. Converge this one "// &
                           "by DIIS.")
            return
         end if
         if (present(projector)) then
            call error%set(ERROR_VALIDATION, "keywords.scf.second_order with a Fock "// &
                           "projector: the rotations the projector forbids are not "// &
                           "excluded from the Newton step's parameter space, so the "// &
                           "step would break the constraint. Converge this one by "// &
                           "DIIS.")
            return
         end if
      end if

      ! Geometric direct minimisation, the other way to finish what DIIS
      ! started. `accelerator: gdm` and `gdm: true` are the same request.
      ctrl%gdm = accel_gdm
      if (present(scf)) ctrl%gdm = ctrl%gdm .or. scf%gdm
      if (present(gdm)) ctrl%gdm = gdm .or. accel_gdm
      ctrl%gdm_start = DEFAULT_GDM_START
      ctrl%gdm_max_rotation = DEFAULT_GDM_MAX_ROTATION
      if (present(scf)) then
         if (scf%gdm_start > 0.0_dp) ctrl%gdm_start = scf%gdm_start
         if (scf%gdm_max_rotation > 0.0_dp) ctrl%gdm_max_rotation = scf%gdm_max_rotation
      end if
      if (present(gdm_start)) then
         if (gdm_start > 0.0_dp) ctrl%gdm_start = gdm_start
      end if
      if (present(gdm_max_rotation)) then
         if (gdm_max_rotation > 0.0_dp) ctrl%gdm_max_rotation = gdm_max_rotation
      end if
      if (ctrl%gdm) then
         if (ctrl%second_order) then
            call error%set(ERROR_VALIDATION, "keywords.scf.gdm together with "// &
                           "keywords.scf.second_order: each finishes the SCF that "// &
                           "DIIS started, and only one of them can. Choose one.")
            return
         end if
         if (ops%cosx%active()) then
            call error%set(ERROR_VALIDATION, "keywords.scf.gdm with chain-of-spheres "// &
                           "exchange: that energy is not variational in the orbitals, "// &
                           "so its Fock matrix is not the gradient of the energy a "// &
                           "step is accepted or rejected on. Converge this one by "// &
                           "DIIS.")
            return
         end if
         if (ctrl%use_pcm) then
            call error%set(ERROR_VALIDATION, "keywords.scf.gdm with a continuum "// &
                           "solvent: the minimisation rebuilds the Fock matrix "// &
                           "itself, and its builds do not carry the reaction field. "// &
                           "Converge this one by DIIS.")
            return
         end if
         if (present(projector)) then
            call error%set(ERROR_VALIDATION, "keywords.scf.gdm with a Fock projector: "// &
                           "the rotations the projector forbids are not excluded "// &
                           "from the step, which would break the constraint. "// &
                           "Converge this one by DIIS.")
            return
         end if
      end if

      ! One call per iteration, and the three views above are what makes that
      ! readable: the alternative was a forty-one argument call.
      !
      ! `exit` cannot live inside the callee, so convergence comes back on
      ! `result%converged` -- which is where it was recorded anyway -- and the
      ! error check is the caller's, exactly as it is after every other call
      ! in this routine.
      soscf_switch = 0
      gdm_switch = 0
      do iter = 1, max_iter
         call do_rhf_iteration(iter, mol, ops, ctrl, st, clk, xc, pcm, projector, &
                               result, error)
         if (error%has_error()) return
         if (result%converged) exit
         ! The handover. Tested on the commutator the iteration just measured,
         ! so the DIIS phase ends on a diagonalisation like any other and the
         ! Newton phase starts from orbitals and a density that agree.
         if (ctrl%second_order .and. iter < max_iter .and. &
             result%commutator < ctrl%soscf_start) then
            soscf_switch = iter
            exit
         end if
         if (ctrl%gdm .and. iter < max_iter .and. &
             result%commutator < ctrl%gdm_start) then
            gdm_switch = iter
            exit
         end if
      end do

      if (soscf_switch > 0) then
         call run_soscf_phase(mol, ops, ctrl, st, clk, xc, soscf_switch, &
                              max_iter - soscf_switch, result, error)
         if (error%has_error()) return
      end if
      if (gdm_switch > 0) then
         call run_gdm_phase(mol, ops, ctrl, st, clk, xc, gdm_switch, &
                            max_iter - gdm_switch, result, error)
         if (error%has_error()) return
      end if

      ! The energy that goes out belongs to the density that satisfied the
      ! test, so it is recomputed from the final Fock. Deliberately without
      ! `incr`: this one is a full build, and the number that leaves must not
      ! carry the corrections accumulated since the last reset.
      !
      ! A chain-of-spheres run takes it on the larger of its two grids when it
      ! has one: the density was converged on the smaller, and the energy's
      ! error in the density is second order while its error in the quadrature
      ! is first.
      !
      ! Only the energy is wanted from this build, so a chain-of-spheres run
      ! takes its exchange as a trace and the Fock matrix it leaves, which has
      ! no exchange in it, goes into a scratch array and not into the state.
      if (ops%cosx%active()) then
         allocate (fock_energy(ops%n_ao, ops%n_ao))
         if (ops%cosx_final%active()) then
            call assemble_fock(mol, ops%h, st%density, st%coeff, ops%n_occ, ops%bmat, ops%eri, &
                               ops%bounds, xc, fock_energy, result%electronic, error, clk=clk, &
                               bmat_lr=ops%bmat_lr, cosx=ops%cosx_final, jfit=ops%jfit, &
                               energy_only=.true.)
         else
            call assemble_fock(mol, ops%h, st%density, st%coeff, ops%n_occ, ops%bmat, ops%eri, &
                               ops%bounds, xc, fock_energy, result%electronic, error, clk=clk, &
                               bmat_lr=ops%bmat_lr, cosx=ops%cosx, jfit=ops%jfit, &
                               energy_only=.true.)
         end if
         deallocate (fock_energy)
      else
         call assemble_fock(mol, ops%h, st%density, st%coeff, ops%n_occ, ops%bmat, ops%eri, &
                            ops%bounds, xc, st%fock, result%electronic, error, clk=clk, &
                            bmat_lr=ops%bmat_lr, cosx=ops%cosx, jfit=ops%jfit)
      end if
      if (error%has_error()) return
      ! Last use of the fitted tensor: a caller that asked for it takes the
      ! allocation itself rather than a copy of several gigabytes.
      if (present(b_ao_out) .and. allocated(ops%bmat)) call move_alloc(ops%bmat, b_ao_out)
      if (ctrl%use_pcm) then
         call pcm%operator_matrix(mol, st%density, st%v_pcm, e_pcm, error)
         if (error%has_error()) return
         result%electronic = result%electronic + e_pcm
         result%pcm_energy = e_pcm
         result%pcm_charge = pcm%q_total
         ! Lapped like the iterations' operators, so the stage counts one per
         ! Fock build rather than one fewer.
         call clk%lap(STAGE_PCM)
      end if
      result%nuclear_repulsion = mol%nuclear_repulsion()
      result%energy = result%electronic + result%nuclear_repulsion
      result%n_occupied = ops%n_occ
      call move_alloc(st%eigenvalues, result%orbital_energies)
      call move_alloc(st%coeff, result%orbitals)
      call move_alloc(st%density, result%density)
      call st%diis%destroy()

      ! The final rebuild above lands in the Fock bucket, charged by
      ! `assemble_fock` itself, so there is deliberately no lap here.
      call clk%finish()
      call scf_table_footer(verbose, result%converged, result%iterations)
      ! The honest cost of a second-order run.
      !
      ! **The Hessian-vector products and nothing else**, because they are the
      ! one part of the cost nothing else reports. The timing table below counts
      ! every `assemble_fock` in its "Fock builds" row -- the DIIS iterations,
      ! the Newton phase's entry build and one per trial step -- and that row is
      ! directly comparable with a DIIS run's. The products are an integral pass
      ! each and appear in no row at all, so a reader comparing the two would
      ! otherwise miss the largest term in the second-order one.
      !
      ! `full_fock_builds` is deliberately not quoted here: it counts only what
      ! the incremental bookkeeping saw, which on the in-core and fitted paths
      ! is nothing, so it would read as a much cheaper run than the timing row
      ! says.
      if (result%gdm_started_at > 0 .and. verbose) then
         write (line, "(a,i0,a,i0,a,i0,a)") &
            "  direct minimisation: handed over after DIIS iteration ", &
            result%gdm_started_at, ", ", result%gdm_iterations, &
            " steps of one Fock build each, ", result%gdm_rejected, " rejected"
         call logger%info(trim(line))
      end if
      if (result%second_order_started_at > 0 .and. verbose) then
         write (line, "(a,i0,a,i0,a,i0,a)") &
            "  second order: handed over after DIIS iteration ", &
            result%second_order_started_at, ", ", result%second_order_iterations, &
            " Newton steps, ", result%second_order_hessian_products, &
            " Hessian-vector products (one integral pass each, counted in no"// &
            " timing row)"
         call logger%info(trim(line))
      end if
      call energy_components(verbose, mol, result%density, result%electronic, &
                             result%nuclear_repulsion)
      call screening_summary(verbose, st%screening)
      call clk%report("RHF", verbose)
      call xc_breakdown_summary(verbose)
   end subroutine run_czt_rhf

   subroutine do_rhf_iteration(iter, mol, ops, ctrl, st, clk, xc, pcm, projector, &
                               result, error)
      !! One SCF iteration: build, extrapolate, shift, diagonalise, test
      !!
      !! Lifted out of `run_czt_rhf`, whose body this was. Two reasons, and the
      !! second is not optional: a 523-line routine is hard to follow, and
      !! nvfortran's optimiser segfaults compiling one that size at -O2 and
      !! above -- every version from 25.11 to 26.5, with no `-Mno<pass>` that
      !! avoids it. Nothing in here is what upsets it; the routine was simply
      !! over some internal limit, which splitting it puts it back under.
      !!
      !! **Convergence leaves on `result%converged`, not by exiting.** A callee
      !! cannot `exit` its caller's loop, and making it return a separate flag
      !! would mean two places recording the same fact.
      integer, intent(in) :: iter
      type(czt_molecule_t), intent(in) :: mol
      type(rhf_operators_t), intent(in) :: ops
      type(rhf_controls_t), intent(in) :: ctrl
      type(rhf_state_t), intent(inout) :: st
      type(timing_report_t), intent(inout) :: clk
      type(xc_context_t), intent(inout), optional :: xc
      type(pcm_context_t), intent(inout), optional :: pcm
      type(fock_projector_t), intent(in), optional :: projector
      type(rhf_result_t), intent(inout) :: result
      type(error_t), intent(inout) :: error

      real(dp) :: de, drms, gnorm, e_pcm
      real(dp) :: shift_now
      real(dp) :: t_fock_iter, t_rest_iter, t_xc_iter
      logical :: extrapolated

      st%density_old = st%density
      ! The energy belongs to the Fock built from this density, so both come
      ! back together and before extrapolation. A DIIS-mixed Fock is a
      ! convergence device, not a state anything is the energy of.
      !
      ! Read the buckets BEFORE the build: `assemble_fock` closes the Fock
      ! stage itself, so a read afterwards would difference against itself.
      t_fock_iter = clk%seconds_of(STAGE_FOCK)
      t_xc_iter = clk%seconds_of(STAGE_XC)
      ! `incr` present is what switches incremental building on inside
      ! `assemble_fock`, so withholding it is how the deck turns it off.
      !
      ! A chain-of-spheres run starts on its coarse grid and leaves it, for
      ! good, once the previous iteration's commutator says the density is
      ! close. `result%commutator` is still the previous iteration's here.
      if (st%cosx_coarse .and. iter > 1) then
         if (result%commutator < COSX_COARSE_UNTIL) st%cosx_coarse = .false.
      end if
      if (st%cosx_coarse .and. ops%cosx_coarse%active()) then
         st%incr%full_builds = st%incr%full_builds + 1
         call assemble_fock(mol, ops%h, st%density, st%coeff, ops%n_occ, ops%bmat, ops%eri, ops%bounds, xc, &
                            st%fock, st%e_elec, error, clk=clk, screening=st%screening, &
                            bmat_lr=ops%bmat_lr, cosx=ops%cosx_coarse, jfit=ops%jfit)
      else if (ctrl%incremental) then
         call assemble_fock(mol, ops%h, st%density, st%coeff, ops%n_occ, ops%bmat, ops%eri, ops%bounds, xc, &
                            st%fock, st%e_elec, error, clk=clk, screening=st%screening, incr=st%incr, &
                            bmat_lr=ops%bmat_lr, cosx=ops%cosx, jfit=ops%jfit)
      else
         ! Counted here rather than in the assembler: without `incr` it has
         ! no state to record into, and every build here is a full one.
         st%incr%full_builds = st%incr%full_builds + 1
         call assemble_fock(mol, ops%h, st%density, st%coeff, ops%n_occ, ops%bmat, ops%eri, ops%bounds, xc, &
                            st%fock, st%e_elec, error, clk=clk, screening=st%screening, &
                            bmat_lr=ops%bmat_lr, cosx=ops%cosx, jfit=ops%jfit)
      end if
      if (error%has_error()) return
      ! The continuum, after the gas-phase pieces and before the commutator:
      ! its operator belongs to the Fock matrix DIIS extrapolates and the
      ! eigenproblem diagonalises, or the converged orbitals would not feel
      ! the solvent. The energy carries its own half already.
      !
      ! Charged to its OWN stage. `lap` also increments a stage's call count,
      ! and `assemble_fock` has already closed STAGE_FOCK, so lapping that
      ! again here would report two Fock builds per iteration.
      if (ctrl%use_pcm) then
         call pcm%operator_matrix(mol, st%density, st%v_pcm, e_pcm, error)
         if (error%has_error()) return
         st%fock = st%fock + st%v_pcm
         st%e_elec = st%e_elec + e_pcm
         call clk%lap(STAGE_PCM)
      end if
      t_fock_iter = clk%seconds_of(STAGE_FOCK) - t_fock_iter
      t_xc_iter = clk%seconds_of(STAGE_XC) - t_xc_iter

      if (present(projector)) then
         call projector%apply(st%fock, error)
         if (error%has_error()) return
      end if

      ! FDS - SDF, projected into the orthogonal basis. It vanishes exactly
      ! when F and D commute, which is what convergence means, so it is the
      ! quantity worth extrapolating against.
      call commutator(st%fock, st%density, ops%s, ops%x, st%err)
      gnorm = maxval(abs(st%err))
      if (gnorm < INCREMENTAL_SETTLED .and. .not. st%incr%settled) then
         st%incr%settled = .true.
         st%incr%settled_at = iter
      end if
      st%fock_flat = reshape(st%fock, [ops%n_ao*ops%n_ao])
      call st%diis%push(st%fock_flat, reshape(st%err, [ops%n_mo*ops%n_mo]), &
                        density=reshape(st%density, [ops%n_ao*ops%n_ao]), energy=st%e_elec)
      call st%diis%extrapolate_with(scheme_now(ctrl%accel, gnorm), st%fock_flat, extrapolated)
      if (extrapolated) st%fock = reshape(st%fock_flat, [ops%n_ao, ops%n_ao])
      t_rest_iter = clk%seconds_of(STAGE_DIIS)
      call clk%lap(STAGE_DIIS)
      t_rest_iter = clk%seconds_of(STAGE_DIIS) - t_rest_iter

      ! **After DIIS, never before it.** The error vector above is built
      ! from the unshifted Fock, and has to be: shift first and the vectors
      ! DIIS stores stop being a subspace of Fock matrices. The energy is
      ! likewise already taken, from the build above.
      !
      ! The virtual projector without the virtual orbitals: completeness
      ! gives `C_o C_o^T + C_v C_v^T = S^-1`, so
      ! `S C_v C_v^T S = S - S C_o C_o^T S = S - (1/2) S D S` for a closed
      ! shell, which is two products against the density already in hand.
      shift_now = 0.0_dp
      if (ctrl%shift > 0.0_dp .and. st%drms_prev > ctrl%taper) shift_now = ctrl%shift
      if (shift_now > 0.0_dp) then
         if (.not. allocated(st%sd)) allocate (st%sd(ops%n_ao, ops%n_ao), st%sds(ops%n_ao, ops%n_ao))
         call pic_gemm(ops%s, st%density, st%sd, beta=0.0_dp)
         call pic_gemm(st%sd, ops%s, st%sds, beta=0.0_dp)
         st%fock = st%fock + shift_now*(ops%s - 0.5_dp*st%sds)
      end if

      call diagonalize(st%fock, ops%x, ops%n_ao, ops%n_mo, st%coeff, st%eigenvalues, error)
      if (error%has_error()) return
      call build_density_closed_shell(st%coeff, ops%n_occ, st%density)
      t_rest_iter = t_rest_iter - clk%seconds_of(STAGE_DIAG)
      call clk%lap(STAGE_DIAG)
      t_rest_iter = t_rest_iter + clk%seconds_of(STAGE_DIAG)

      de = abs(st%e_elec - st%e_old)
      drms = sqrt(sum((st%density - st%density_old)**2)/real(ops%n_ao*ops%n_ao, dp))
      call scf_table_row(ctrl%verbose, iter, st%e_elec + mol%nuclear_repulsion(), de, gnorm, &
                         st%diis%count(), t_fock_iter, t_xc_iter, t_rest_iter, ctrl%kohn_sham)

      st%e_old = st%e_elec
      result%iterations = iter
      result%commutator = gnorm
      result%full_fock_builds = st%incr%full_builds
      result%incremental_settled_at = st%incr%settled_at
      result%incremental_updates = st%incr%updates
      ! `shift_now` is part of the test below rather than checked
      ! afterwards: the orbitals and eigenvalues that leave here are the ones
      ! the last diagonalisation produced, so convergence has to be declared
      ! on an unshifted iteration or every virtual energy is high by `shift`.
      st%drms_prev = drms
      ! **The energy and the commutator, and not the density.** `de` and
      ! `drms` say the iteration stopped moving; they do not say it stopped
      ! at a stationary point. `FDS - SDF` is what vanishes when F and D
      ! commute, and an SCF can hold the other two small while this one is
      ! nowhere near zero -- any scheme that interpolates rather than
      ! extrapolates will do it. `drms` is still computed, because the
      ! level-shift taper is driven from it.
      !
      ! The bound has to clear the commutator's own noise floor, which moves
      ! with thread count because the OpenMP reduction merges are unordered.
      ! `sqrt(energy_tol)` clears it by orders at the default `energy_tol`.
      !
      ! The rule lives in `scf_convergence_t`. `shift_now` stays in the
      ! caller because it is not a convergence measure: a shifted Fock matrix
      ! is a different operator, so an iterate that met the threshold under a
      ! shift has not converged the problem that was asked for.
      if (iter > 1 .and. ctrl%conv%is_converged(de, drms, gnorm) .and. &
          shift_now == 0.0_dp) then
         result%converged = .true.
      end if
   end subroutine do_rhf_iteration

   subroutine do_uhf_iteration(iter, mol, h, s, x, eri, bounds, n_ao, n_mo, n_alpha, n_beta, &
                               nsq, msq, want_incremental, use_pcm, shift, taper, start_cycle, &
                               accel, conv, verbose, kohn_sham_run, st, clk, xc, pcm, result, error)
      !! One unrestricted SCF iteration: build, extrapolate, shift, diagonalise, test
      !!
      !! `do_rhf_iteration`'s twin, lifted out of `run_czt_uhf` for the same two
      !! reasons and in the same shape. The second is not optional: nvfortran's
      !! optimiser segfaults compiling a routine that size at -O2 and above --
      !!
      !!   nvfortran-Fatal-.../tools/fort2 TERMINATED by signal 11
      !!
      !! -- on every version from 25.11 to 26.5. Nothing in here is what upsets
      !! it and no `-Mno<pass>` avoids it: measured on 26.5, only `-O1` and
      !! `-O0` compiled the unsplit routine, while `-Mnoopt` and `-Mnovect`
      !! still crashed. The unrestricted driver was simply over some internal
      !! limit that splitting it puts it back under, exactly as the restricted
      !! one was before it.
      !!
      !! **Convergence leaves on `result%converged`, not by exiting**, as on the
      !! restricted path: a callee cannot `exit` its caller's loop.
      integer, intent(in) :: iter
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: h(:, :), s(:, :), x(:, :)
      real(dp), intent(in), allocatable :: eri(:, :, :, :)
      real(dp), intent(in), allocatable :: bounds(:, :)
      integer, intent(in) :: n_ao, n_mo, n_alpha, n_beta, nsq, msq
      logical, intent(in) :: want_incremental, use_pcm
      real(dp), intent(in) :: shift, taper
      integer, intent(in) :: start_cycle, accel
      type(scf_convergence_t), intent(in) :: conv
      logical, intent(in) :: verbose, kohn_sham_run
      type(uhf_state_t), intent(inout) :: st
      type(timing_report_t), intent(inout) :: clk
      type(xc_context_t), intent(inout), optional :: xc
      type(pcm_context_t), intent(inout), optional :: pcm
      type(rhf_result_t), intent(inout) :: result
      type(error_t), intent(inout) :: error

      real(dp) :: de, drms, gnorm, e_pcm, shift_now
      real(dp) :: t_fock_iter, t_rest_iter, t_xc_iter
      logical :: extrapolated

      st%d_a_old = st%d_a
      st%d_b_old = st%d_b

      ! Read before the build, not after it: both stages are cumulative, so
      ! a reading taken once the build has charged them differences against
      ! itself.
      t_fock_iter = clk%seconds_of(STAGE_FOCK)
      t_xc_iter = clk%seconds_of(STAGE_XC)
      ! As on the restricted path: `st%incr` present is the switch, so the deck
      ! turns incremental building off by withholding it.
      if (want_incremental) then
         call assemble_fock_uhf(mol, h, st%d_a, st%d_b, eri, bounds, xc, st%fock_a, st%fock_b, &
                                st%e_elec, error, clk=clk, incr=st%incr, bmat=st%bmat, &
                                bmat_lr=st%bmat_lr)
      else
         ! No `st%incr` to record into, and every build here is a full one.
         st%incr%full_builds = st%incr%full_builds + 1
         call assemble_fock_uhf(mol, h, st%d_a, st%d_b, eri, bounds, xc, st%fock_a, st%fock_b, &
                                st%e_elec, error, clk=clk, bmat=st%bmat, bmat_lr=st%bmat_lr)
      end if
      if (error%has_error()) return
      ! The continuum sees one density -- the total -- and both spins feel
      ! one potential, as they would from any classical charge.
      if (use_pcm) then
         call pcm%operator_matrix(mol, st%d_a + st%d_b, st%v_pcm, e_pcm, error)
         if (error%has_error()) return
         st%fock_a = st%fock_a + st%v_pcm
         st%fock_b = st%fock_b + st%v_pcm
         st%e_elec = st%e_elec + e_pcm
      end if
      ! TODO(mqc): `assemble_fock_uhf` already laps STAGE_FOCK when it is
      ! given a clock, so this second lap counts two Fock builds per
      ! iteration and charges the continuum solve to the Fock bucket, where
      ! the restricted path gives it STAGE_PCM.
      call clk%lap(STAGE_FOCK)
      t_fock_iter = clk%seconds_of(STAGE_FOCK) - t_fock_iter
      t_xc_iter = clk%seconds_of(STAGE_XC) - t_xc_iter

      if (st%rohf) then
         ! One matrix and one error for both spins. The error is against the
         ! total density: in the orbital basis that density is diagonal with
         ! occupations two, one and zero, and commutes with the effective
         ! Fock matrix exactly when its three off-diagonal blocks vanish.
         call roothaan_fock(st%fock_a, st%fock_b, st%d_a, st%d_b, s)
         call commutator(st%fock_a, st%d_a + st%d_b, s, x, st%err_a)
         st%err_b = st%err_a
      else
         call commutator(st%fock_a, st%d_a, s, x, st%err_a)
         call commutator(st%fock_b, st%d_b, s, x, st%err_b)
      end if
      st%fock_flat(1:nsq) = reshape(st%fock_a, [nsq])
      st%fock_flat(nsq + 1:2*nsq) = reshape(st%fock_b, [nsq])
      st%err_flat(1:msq) = reshape(st%err_a, [msq])
      st%err_flat(msq + 1:2*msq) = reshape(st%err_b, [msq])
      gnorm = maxval(abs(st%err_flat))
      if (gnorm < INCREMENTAL_SETTLED .and. .not. st%incr%settled) then
         st%incr%settled = .true.
         st%incr%settled_at = iter
      end if
      st%dens_flat(1:nsq) = reshape(st%d_a, [nsq])
      st%dens_flat(nsq + 1:2*nsq) = reshape(st%d_b, [nsq])
      call st%diis%push(st%fock_flat, st%err_flat, density=st%dens_flat, energy=st%e_elec)
      extrapolated = .false.
      if (iter >= start_cycle) then
         call st%diis%extrapolate_with(scheme_now(accel, gnorm), st%fock_flat, extrapolated)
      end if
      if (extrapolated) then
         st%fock_a = reshape(st%fock_flat(1:nsq), [n_ao, n_ao])
         st%fock_b = reshape(st%fock_flat(nsq + 1:2*nsq), [n_ao, n_ao])
      end if
      t_rest_iter = clk%seconds_of(STAGE_DIIS)
      call clk%lap(STAGE_DIIS)
      t_rest_iter = clk%seconds_of(STAGE_DIIS) - t_rest_iter

      ! After the extrapolation and after the two commutators, for the reason
      ! `run_czt_rhf` sets out. Each spin gets its own virtual projector:
      ! `build_density_spin` gives an occupation of one, so the closed-shell
      ! factor of a half is absent and the projector is `S - S D_sigma S`.
      shift_now = 0.0_dp
      if (shift > 0.0_dp .and. st%drms_prev > taper) shift_now = shift
      if (shift_now > 0.0_dp) then
         if (.not. allocated(st%sd)) allocate (st%sd(n_ao, n_ao), st%sds(n_ao, n_ao))
         call pic_gemm(s, st%d_a, st%sd, beta=0.0_dp)
         call pic_gemm(st%sd, s, st%sds, beta=0.0_dp)
         st%fock_a = st%fock_a + shift_now*(s - st%sds)
         call pic_gemm(s, st%d_b, st%sd, beta=0.0_dp)
         call pic_gemm(st%sd, s, st%sds, beta=0.0_dp)
         st%fock_b = st%fock_b + shift_now*(s - st%sds)
      end if

      call diagonalize(st%fock_a, x, n_ao, n_mo, st%coeff_a, st%eig_a, error)
      if (error%has_error()) return
      if (st%rohf) then
         ! `fock_a` is the effective matrix, shifted against alpha's virtual
         ! projector if at all; beta takes the same orbitals.
         st%coeff_b = st%coeff_a
         st%eig_b = st%eig_a
      else
         call diagonalize(st%fock_b, x, n_ao, n_mo, st%coeff_b, st%eig_b, error)
         if (error%has_error()) return
      end if
      call build_density_spin(st%coeff_a, n_alpha, st%d_a)
      call build_density_spin(st%coeff_b, n_beta, st%d_b)
      t_rest_iter = t_rest_iter - clk%seconds_of(STAGE_DIAG)
      call clk%lap(STAGE_DIAG)
      t_rest_iter = t_rest_iter + clk%seconds_of(STAGE_DIAG)

      de = abs(st%e_elec - st%e_old)
      drms = sqrt((sum((st%d_a - st%d_a_old)**2) + sum((st%d_b - st%d_b_old)**2))/real(2*nsq, dp))
      call scf_table_row(verbose, iter, st%e_elec + mol%nuclear_repulsion(), de, gnorm, &
                         st%diis%count(), t_fock_iter, t_xc_iter, t_rest_iter, kohn_sham_run)

      st%e_old = st%e_elec
      result%iterations = iter
      result%commutator = gnorm
      result%full_fock_builds = st%incr%full_builds
      result%incremental_settled_at = st%incr%settled_at
      result%incremental_updates = st%incr%updates
      st%drms_prev = drms
      ! **The energy and the commutator, and not the density.** `de` and
      ! `drms` say the iteration stopped moving; they do not say it stopped
      ! at a stationary point. `FDS - SDF` is what vanishes when F and D
      ! commute, and an SCF can hold the other two small while this one is
      ! nowhere near zero -- any scheme that interpolates rather than
      ! extrapolates will do it. `drms` is still computed, because the
      ! level-shift taper is driven from it.
      !
      ! The bound has to clear the commutator's own noise floor, which moves
      ! with thread count because the OpenMP reduction merges are unordered.
      ! `sqrt(energy_tol)` clears it by orders at the default `energy_tol`.
      !
      ! The rule lives in `scf_convergence_t`. `shift_now` stays in the
      ! caller because it is not a convergence measure: a shifted Fock matrix
      ! is a different operator, so an iterate that met the threshold under a
      ! shift has not converged the problem that was asked for.
      if (iter > 1 .and. conv%is_converged(de, drms, gnorm) .and. &
          shift_now == 0.0_dp) then
         result%converged = .true.
      end if
   end subroutine do_uhf_iteration

   subroutine soscf_table_header(verbose)
      !! Column headings for the second-order phase
      logical, intent(in) :: verbose

      if (.not. verbose) return
      call logger%info("  second-order phase: trust-region Newton on the orbital rotations")
      call logger%info("    iter                 energy          dE     |FDS-SDF|"// &
                       "      |g|max       trust   curvature  fock")
   end subroutine soscf_table_header

   subroutine soscf_table_row(verbose, iter, energy, de, gnorm, largest, trust, &
                              lowest, builds)
      !! One line of the second-order phase
      !!
      !! `curvature` is the smallest one the Newton step saw, which is an upper
      !! bound on the true smallest eigenvalue -- see `mqc_czt_soscf`. `fock` is
      !! the Fock builds this iteration cost, energies and Hessian-vector
      !! products together, because that is the number the method has to be
      !! judged on.
      logical, intent(in) :: verbose
      integer, intent(in) :: iter, builds
      real(dp), intent(in) :: energy, de, gnorm, largest, trust, lowest

      character(len=LINE_LEN) :: line

      if (.not. verbose) return
      write (line, "(i8,f22.12,4es12.3,es12.2,i6)") iter, energy, de, gnorm, &
         largest, trust, lowest, builds
      call logger%info(trim(line))
   end subroutine soscf_table_row

   subroutine run_soscf_phase(mol, ops, ctrl, st, clk, xc, switch_iter, max_steps, &
                              result, error)
      !! Converge a closed-shell SCF by trust-region Newton, from where DIIS left it
      !!
      !! ## Why this is not where the SCF starts
      !!
      !! A Newton step is a step on a quadratic model of the energy, and that
      !! model is only the energy near the point it was built at. From an
      !! initial guess it is not, and the step it produces is long, in an
      !! arbitrary direction, and routinely uphill -- a second-order SCF
      !! started from a guess diverges, which is why every one in practice is
      !! preceded by something else. Here that something else is the DIIS SCF
      !! above, run until `max|FDS - SDF|` falls below
      !! `keywords.scf.soscf_start`; `switch_iter` is where that happened and
      !! is reported, because a reader has to be able to see it.
      !!
      !! ## What one iteration is
      !!
      !! 1. **Semicanonicalise.** The Hessian expression is only the Hessian
      !!    when the Fock matrix is diagonal within the occupied block and
      !!    within the virtual block, and a rotation destroys that. Two small
      !!    diagonalisations restore it without moving the density, so this is
      !!    a change of representation and not a step.
      !! 2. **Differentiate.** The gradient is `4 F_ai` and the Hessian is
      !!    four times `mqc_czt_ov_hessian`'s `(A+B)`, applied to a vector for
      !!    one Fock build and never written down.
      !! 3. **Step.** `soscf_newton_step` solves the level-shifted Newton
      !!    equations in a Krylov subspace, with the same `level_shifted_step`
      !!    -- the same level shift, the same saddle escape -- that the CASSCF
      !!    optimiser applies to its dense Hessian.
      !! 4. **Backtrack.** The step is scaled to the trust radius, the
      !!    orbitals are rotated by the exact matrix exponential, and the
      !!    energy is rebuilt. If it did not fall, the radius is halved and the
      !!    step retried. Nothing is accepted that does not lower the energy.
      !!
      !! ## What convergence means here
      !!
      !! The same rule the DIIS phase used -- `scf_convergence_t` on the energy
      !! change and the same commutator -- **and** no negative curvature. That
      !! second half is the entire reason to have a second-order method: a
      !! first-order SCF stops wherever the gradient vanishes, and on a
      !! symmetric guess that can be a saddle whose symmetry-breaking rotations
      !! have exactly zero gradient. The curvature test here is one sided,
      !! because the Krylov value is an upper bound: it can prove a saddle and
      !! refuse to converge, and it cannot prove a minimum. `keywords.scf.
      !! stability` is what proves a minimum.
      !!
      !! ## The honest cost
      !!
      !! Each iteration costs one Fock build per accepted or rejected trial
      !! step, plus one per Hessian-vector product inside the Newton solve. On
      !! a well-behaved closed shell that is more Fock builds than DIIS would
      !! have spent, for fewer iterations. It earns its cost where DIIS
      !! oscillates, and where the answer has to be a minimum and not merely a
      !! stationary point.
      type(czt_molecule_t), intent(in), target :: mol
      type(rhf_operators_t), intent(in), target :: ops
      type(rhf_controls_t), intent(in) :: ctrl
      type(rhf_state_t), intent(inout) :: st
      type(timing_report_t), intent(inout) :: clk
      type(xc_context_t), intent(inout), optional, target :: xc
      integer, intent(in) :: switch_iter
         !! The DIIS iteration this took over after. Reported, not used.
      integer, intent(in) :: max_steps
         !! Newton steps allowed: whatever is left of `max_iter`.
      type(rhf_result_t), intent(inout) :: result
      type(error_t), intent(inout) :: error

      ! `target`, and load-bearing: the Hessian points at the response
      ! operator, so it has to outlive every product taken through it.
      type(nuclear_response_t), target :: response
      type(ov_hessian_t) :: hessian
      real(dp), allocatable :: fock_mo(:, :), energies(:), gradient(:)
      real(dp), allocatable :: kappa(:, :), rotation(:, :), step_matrix(:, :)
      real(dp), allocatable :: fock_cur(:, :), fock_try(:, :)
      real(dp), allocatable :: coeff_try(:, :), density_try(:, :)
      real(dp) :: trust, scaling, lowest, predicted
      real(dp) :: e_cur, e_try, de, drms, gnorm, largest
      integer :: iteration, trial, products, builds, n_ao, n_mo
      logical :: accepted, saddle

      integer, parameter :: MAX_BACKTRACKS = 12

      if (error%has_error()) return
      n_ao = ops%n_ao
      n_mo = ops%n_mo
      result%second_order_started_at = switch_iter
      if (max_steps < 1) return

      allocate (fock_cur(n_ao, n_ao), fock_try(n_ao, n_ao))
      allocate (coeff_try(n_ao, n_mo), density_try(n_ao, n_ao))
      allocate (step_matrix(n_mo, n_mo))

      ! The point the Newton phase starts from, built in full. The Fock matrix
      ! left behind by the DIIS loop is extrapolated and possibly level
      ! shifted, and the energy it carries belongs to the density before the
      ! last diagonalisation -- neither is a state this can differentiate.
      call assemble_fock(mol, ops%h, st%density, st%coeff, ops%n_occ, ops%bmat, &
                         ops%eri, ops%bounds, xc, fock_cur, e_cur, error, clk=clk, &
                         bmat_lr=ops%bmat_lr, jfit=ops%jfit)
      if (error%has_error()) return
      st%incr%full_builds = st%incr%full_builds + 1

      call soscf_table_header(ctrl%verbose)
      trust = MAX_ROTATION

      do iteration = 1, max_steps
         builds = 0

         ! ---- 1. semicanonicalise, so the orbital energies mean something --
         if (allocated(fock_mo)) deallocate (fock_mo, energies)
         call soscf_semicanonicalize(fock_cur, st%coeff, ops%n_occ, fock_mo, &
                                     energies, error)
         if (error%has_error()) return

         ! ---- 2. the gradient, and the commutator the caller stops on -------
         if (allocated(gradient)) deallocate (gradient)
         call soscf_gradient(fock_mo, ops%n_occ, gradient)
         largest = maxval(abs(gradient))
         call commutator(fock_cur, st%density, ops%s, ops%x, st%err)
         gnorm = maxval(abs(st%err))
         ! Recorded here, not after the step: on the iteration that converges
         ! there is no step, and the number a caller reads has to be the one
         ! that satisfied the test.
         result%commutator = gnorm
         de = abs(e_cur - st%e_old)
         drms = sqrt(sum((st%density - st%density_old)**2)/real(n_ao*n_ao, dp))

         ! ---- 3. the step, which also reports the curvature ----------------
         !
         ! Built before the convergence test rather than after it, because the
         ! test needs the curvature: a vanished gradient at a saddle is still a
         ! vanished gradient.
         if (allocated(kappa)) deallocate (kappa)
         call build_soscf_hessian(mol, ops, ctrl, st, energies, response, hessian, &
                                  xc, error)
         if (error%has_error()) return
         call soscf_newton_step(hessian, gradient, MAX_ROTATION, kappa, lowest, &
                                predicted, products, error)
         if (error%has_error()) return
         call clk%lap(STAGE_SOSCF)
         result%second_order_hessian_products = &
            result%second_order_hessian_products + products
         result%second_order_curvature = lowest
         result%has_second_order_curvature = .true.
         builds = builds + products

         saddle = lowest < SADDLE_CURVATURE
         if (ctrl%conv%is_converged(de, drms, gnorm)) then
            if (.not. saddle) then
               result%converged = .true.
               call soscf_table_row(ctrl%verbose, iteration, &
                                    e_cur + mol%nuclear_repulsion(), de, gnorm, &
                                    largest, trust, lowest, builds)
               exit
            end if
            ! Leaving a saddle is a fresh direction, so it gets a fresh trust
            ! radius: the old one records how the previous direction behaved
            ! and is usually small by the time a run has settled.
            trust = MAX_ROTATION
            call logger%warning("    the SCF has stopped at a saddle point; "// &
                                "following the negative curvature out")
         end if

         ! ---- 4. take it, backtracking until the energy falls ---------------
         accepted = .false.
         do trial = 1, MAX_BACKTRACKS
            scaling = maxval(abs(kappa))
            if (scaling > trust) then
               step_matrix = kappa*(trust/scaling)
            else
               step_matrix = kappa
            end if

            if (allocated(rotation)) deallocate (rotation)
            call rotation_matrix(step_matrix, rotation)
            call pic_gemm(st%coeff, rotation, coeff_try, beta=0.0_dp)
            call build_density_closed_shell(coeff_try, ops%n_occ, density_try)
            ! In full, never incrementally: the trial densities of a
            ! backtracking search do not form a sequence converging on
            ! anything, so an accumulated reference would be a correction to a
            ! density nothing came back to.
            call assemble_fock(mol, ops%h, density_try, coeff_try, ops%n_occ, &
                               ops%bmat, ops%eri, ops%bounds, xc, fock_try, e_try, &
                               error, clk=clk, screening=st%screening, &
                               bmat_lr=ops%bmat_lr, jfit=ops%jfit)
            if (error%has_error()) return
            st%incr%full_builds = st%incr%full_builds + 1
            builds = builds + 1

            ! `predicted` below its own resolution is a step too small to
            ! measure, not a step that failed: near the solution a Newton step
            ! is worth `g^2/H`, which drops under what an energy difference can
            ! report while the gradient is still shrinking.
            if (e_try < e_cur .or. predicted < ENERGY_RESOLUTION) then
               st%density_old = st%density
               st%density = density_try
               st%coeff = coeff_try
               fock_cur = fock_try
               st%e_old = e_cur
               e_cur = e_try
               trust = min(MAX_ROTATION, trust*TRUST_GROWTH)
               accepted = .true.
               exit
            end if
            trust = 0.5_dp*trust
            if (trust < MIN_ROTATION) exit
         end do

         call soscf_table_row(ctrl%verbose, iteration, &
                              e_cur + mol%nuclear_repulsion(), de, gnorm, largest, &
                              trust, lowest, builds)
         result%second_order_iterations = iteration
         result%iterations = switch_iter + iteration

         if (.not. accepted) then
            ! Not an error: the energy and the density in hand are the last
            ! ones that descended, and `converged` is false, which is what the
            ! caller acts on.
            result%second_order_stalled = .true.
            call logger%warning("    the second-order SCF found no step downhill; "// &
                                "stopping. The reported energy is the last one that "// &
                                "descended.")
            exit
         end if
      end do

      ! The orbital energies that leave belong to the last Fock matrix, and are
      ! canonical only because the occupied-virtual block has vanished. On a run
      ! that stopped short of convergence they are semicanonical, which is the
      ! honest thing to report: nothing downstream can canonicalise them without
      ! a Fock matrix of its own.
      !
      ! `st%density` is deliberately left alone. This last rotation mixes
      ! occupied orbitals among themselves and virtuals among themselves, and a
      ! closed-shell density is invariant under both, so rebuilding it would
      ! produce the same matrix to rounding and lose the one the accepted step's
      ! energy was computed from.
      if (allocated(fock_mo)) deallocate (fock_mo, energies)
      call soscf_semicanonicalize(fock_cur, st%coeff, ops%n_occ, fock_mo, energies, &
                                  error)
      if (error%has_error()) return
      st%eigenvalues = energies
      result%full_fock_builds = st%incr%full_builds
      result%incremental_settled_at = st%incr%settled_at
      result%incremental_updates = st%incr%updates

      deallocate (fock_cur, fock_try, coeff_try, density_try, step_matrix)
   end subroutine run_soscf_phase

   subroutine run_soscf_phase_uhf(mol, h, s, x, eri, bounds, n_mo, n_alpha, n_beta, conv, &
                                  verbose, st, clk, xc, switch_iter, max_steps, result, &
                                  error)
      !! Converge an unrestricted SCF by trust-region Newton, from where DIIS
      !! left it
      !!
      !! `run_soscf_phase` over two spins, and that routine's notes hold here:
      !! the semicanonicalisation each iteration, the one-sided curvature
      !! test, the backtracking that accepts nothing uphill. What differs:
      !!
      !! * the rotation vector is alpha's virtual-occupied rectangle and then
      !!   beta's, with one Newton solve over both because the Hessian couples
      !!   them through the Coulomb field;
      !! * the gradient is `2 F_ai` per spin and the Hessian twice the
      !!   spin-blocked `(A+B)` of `response_product_uhf`;
      !! * the rotations searched include the ones that break spin symmetry,
      !!   so on a singlet this can leave a restricted solution for a lower
      !!   unrestricted one, which the closed-shell phase cannot.
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: h(:, :), s(:, :), x(:, :)
      real(dp), intent(in), allocatable :: eri(:, :, :, :)
      real(dp), intent(in), allocatable, target :: bounds(:, :)
      integer, intent(in) :: n_mo, n_alpha, n_beta
      type(scf_convergence_t), intent(in) :: conv
      logical, intent(in) :: verbose
      type(uhf_state_t), intent(inout), target :: st
      type(timing_report_t), intent(inout) :: clk
      type(xc_context_t), intent(inout), optional, target :: xc
      integer, intent(in) :: switch_iter
         !! The DIIS iteration this took over after. Reported, not used.
      integer, intent(in) :: max_steps
         !! Newton steps allowed: whatever is left of `max_iter`.
      type(rhf_result_t), intent(inout) :: result
      type(error_t), intent(inout) :: error

      integer, parameter :: MAX_BACKTRACKS = 12
      type(uhf_rotation_hessian_t) :: hessian
      real(dp), allocatable :: fa_cur(:, :), fb_cur(:, :), fa_try(:, :), fb_try(:, :)
      real(dp), allocatable :: ca_try(:, :), cb_try(:, :), da_try(:, :), db_try(:, :)
      real(dp), allocatable :: fmo_a(:, :), fmo_b(:, :), eps_a(:), eps_b(:)
      real(dp), allocatable :: gradient(:), step(:), kappa(:, :), rotation(:, :)
      real(dp) :: trust, scaling, factor, lowest, predicted
      real(dp) :: e_cur, e_try, de, drms, gnorm, largest
      integer :: iteration, trial, products, builds, n_ao, n_rot_a, n_rot_b, n_rot
      logical :: accepted, saddle

      if (error%has_error()) return
      n_ao = size(h, 1)
      n_rot_a = n_alpha*(n_mo - n_alpha)
      n_rot_b = n_beta*(n_mo - n_beta)
      n_rot = n_rot_a + n_rot_b
      result%second_order_started_at = switch_iter
      if (max_steps < 1 .or. n_rot < 1) return

      allocate (fa_cur(n_ao, n_ao), fb_cur(n_ao, n_ao), fa_try(n_ao, n_ao), fb_try(n_ao, n_ao))
      allocate (ca_try(n_ao, n_mo), cb_try(n_ao, n_mo), da_try(n_ao, n_ao), db_try(n_ao, n_ao))
      allocate (gradient(n_rot), kappa(n_mo, n_mo))

      ! The point the Newton phase starts from, built in full, for the reason
      ! `run_soscf_phase` gives.
      call assemble_fock_uhf(mol, h, st%d_a, st%d_b, eri, bounds, xc, fa_cur, fb_cur, &
                             e_cur, error, clk=clk)
      if (error%has_error()) return
      st%incr%full_builds = st%incr%full_builds + 1

      hessian%mol => mol
      hessian%bounds => bounds
      hessian%c_a => st%coeff_a
      hessian%c_b => st%coeff_b
      hessian%d_a => st%d_a
      hessian%d_b => st%d_b
      if (present(xc)) hessian%xc => xc
      hessian%n_occ_a = n_alpha
      hessian%n_occ_b = n_beta
      hessian%n_mo = n_mo
      allocate (hessian%zero_h(n_ao, n_ao))
      hessian%zero_h = 0.0_dp

      call soscf_table_header(verbose)
      trust = MAX_ROTATION

      do iteration = 1, max_steps
         builds = 0

         ! ---- 1. semicanonicalise each spin ---------------------------------
         call semicanonical(fa_cur, st%coeff_a, n_alpha, fmo_a, eps_a)
         if (error%has_error()) return
         call semicanonical(fb_cur, st%coeff_b, n_beta, fmo_b, eps_b)
         if (error%has_error()) return

         ! ---- 2. the gradient, and the commutator the caller stops on -------
         call spin_gradient(fmo_a, n_alpha, gradient(1:n_rot_a))
         call spin_gradient(fmo_b, n_beta, gradient(n_rot_a + 1:n_rot))
         largest = maxval(abs(gradient))
         call commutator(fa_cur, st%d_a, s, x, st%err_a)
         call commutator(fb_cur, st%d_b, s, x, st%err_b)
         gnorm = max(maxval(abs(st%err_a)), maxval(abs(st%err_b)))
         result%commutator = gnorm
         de = abs(e_cur - st%e_old)
         drms = sqrt((sum((st%d_a - st%d_a_old)**2) + sum((st%d_b - st%d_b_old)**2)) &
                     /real(2*n_ao*n_ao, dp))

         ! ---- 3. the step, which also reports the curvature ----------------
         call spin_gaps(eps_a, n_alpha, hessian%gaps_a)
         call spin_gaps(eps_b, n_beta, hessian%gaps_b)
         if (allocated(step)) deallocate (step)
         call soscf_newton_step_uhf(hessian, gradient, MAX_ROTATION, step, lowest, &
                                    predicted, products, error)
         if (error%has_error()) return
         call clk%lap(STAGE_SOSCF)
         result%second_order_hessian_products = &
            result%second_order_hessian_products + products
         result%second_order_curvature = lowest
         result%has_second_order_curvature = .true.
         builds = builds + products

         saddle = lowest < SADDLE_CURVATURE
         if (conv%is_converged(de, drms, gnorm)) then
            if (.not. saddle) then
               result%converged = .true.
               call soscf_table_row(verbose, iteration, e_cur + mol%nuclear_repulsion(), &
                                    de, gnorm, largest, trust, lowest, builds)
               exit
            end if
            trust = MAX_ROTATION
            call logger%warning("    the SCF has stopped at a saddle point; "// &
                                "following the negative curvature out")
         end if

         ! ---- 4. take it, backtracking until the energy falls ---------------
         accepted = .false.
         scaling = maxval(abs(step))
         do trial = 1, MAX_BACKTRACKS
            factor = 1.0_dp
            if (scaling > trust) factor = trust/scaling
            call rotate(st%coeff_a, n_alpha, step(1:n_rot_a), ca_try)
            call rotate(st%coeff_b, n_beta, step(n_rot_a + 1:n_rot), cb_try)
            call build_density_spin(ca_try, n_alpha, da_try)
            call build_density_spin(cb_try, n_beta, db_try)
            call assemble_fock_uhf(mol, h, da_try, db_try, eri, bounds, xc, fa_try, &
                                   fb_try, e_try, error, clk=clk)
            if (error%has_error()) return
            st%incr%full_builds = st%incr%full_builds + 1
            builds = builds + 1

            if (e_try < e_cur .or. predicted < ENERGY_RESOLUTION) then
               st%d_a_old = st%d_a
               st%d_b_old = st%d_b
               st%d_a = da_try
               st%d_b = db_try
               st%coeff_a = ca_try
               st%coeff_b = cb_try
               fa_cur = fa_try
               fb_cur = fb_try
               st%e_old = e_cur
               e_cur = e_try
               trust = min(MAX_ROTATION, trust*TRUST_GROWTH)
               accepted = .true.
               exit
            end if
            trust = 0.5_dp*trust
            if (trust < MIN_ROTATION) exit
         end do

         call soscf_table_row(verbose, iteration, e_cur + mol%nuclear_repulsion(), de, &
                              gnorm, largest, trust, lowest, builds)
         result%second_order_iterations = iteration
         result%iterations = switch_iter + iteration

         if (.not. accepted) then
            result%second_order_stalled = .true.
            call logger%warning("    the second-order SCF found no step downhill; "// &
                                "stopping. The reported energy is the last one that "// &
                                "descended.")
            exit
         end if
      end do

      ! Semicanonical orbitals and their energies for the last Fock matrices;
      ! the densities are left alone, as in the closed-shell phase.
      call semicanonical(fa_cur, st%coeff_a, n_alpha, fmo_a, eps_a)
      if (error%has_error()) return
      call semicanonical(fb_cur, st%coeff_b, n_beta, fmo_b, eps_b)
      if (error%has_error()) return
      st%eig_a = eps_a
      st%eig_b = eps_b
      st%fock_a = fa_cur
      st%fock_b = fb_cur
      st%e_elec = e_cur
      result%full_fock_builds = st%incr%full_builds
      nullify (hessian%mol, hessian%bounds, hessian%c_a, hessian%c_b, hessian%d_a, &
               hessian%d_b, hessian%xc)

   contains

      subroutine semicanonical(fock_ao, coeff, n_occ, fock_mo, energies)
         !! One spin's Fock matrix diagonalised within its occupied and
         !! virtual blocks. An empty spin is diagonalised whole.
         real(dp), intent(in) :: fock_ao(:, :)
         real(dp), allocatable, intent(inout) :: coeff(:, :)
         integer, intent(in) :: n_occ
         real(dp), allocatable, intent(inout) :: fock_mo(:, :), energies(:)

         integer :: p

         if (allocated(fock_mo)) deallocate (fock_mo)
         if (allocated(energies)) deallocate (energies)
         if (n_occ == 0) then
            call diagonalize(fock_ao, x, n_ao, n_mo, coeff, energies, error)
            if (error%has_error()) return
            allocate (fock_mo(n_mo, n_mo))
            fock_mo = 0.0_dp
            do p = 1, n_mo
               fock_mo(p, p) = energies(p)
            end do
            return
         end if
         call soscf_semicanonicalize(fock_ao, coeff, n_occ, fock_mo, energies, error)
      end subroutine semicanonical

      subroutine spin_gradient(fock_mo, n_occ, g)
         !! `2 F_ai`, virtual fastest
         real(dp), intent(in) :: fock_mo(:, :)
         integer, intent(in) :: n_occ
         real(dp), intent(out) :: g(:)

         integer :: i, a, n_vir

         n_vir = n_mo - n_occ
         do i = 1, n_occ
            do a = 1, n_vir
               g(a + (i - 1)*n_vir) = SPIN_HESSIAN_SCALE*fock_mo(n_occ + a, i)
            end do
         end do
      end subroutine spin_gradient

      subroutine spin_gaps(energies, n_occ, gaps)
         !! `eps_a - eps_i`, `(n_vir, n_occ)`
         real(dp), intent(in) :: energies(:)
         integer, intent(in) :: n_occ
         real(dp), allocatable, intent(inout) :: gaps(:, :)

         integer :: i, a, n_vir

         n_vir = n_mo - n_occ
         if (allocated(gaps)) deallocate (gaps)
         allocate (gaps(n_vir, n_occ))
         do i = 1, n_occ
            do a = 1, n_vir
               gaps(a, i) = energies(n_occ + a) - energies(i)
            end do
         end do
      end subroutine spin_gaps

      subroutine rotate(coeff, n_occ, d, coeff_new)
         !! `coeff exp(factor kappa(d))` for one spin
         real(dp), intent(in) :: coeff(:, :), d(:)
         integer, intent(in) :: n_occ
         real(dp), intent(out) :: coeff_new(:, :)

         integer :: i, a, n_vir

         n_vir = n_mo - n_occ
         if (n_occ*n_vir == 0) then
            coeff_new = coeff
            return
         end if
         kappa = 0.0_dp
         do i = 1, n_occ
            do a = 1, n_vir
               kappa(n_occ + a, i) = factor*d(a + (i - 1)*n_vir)
               kappa(i, n_occ + a) = -kappa(n_occ + a, i)
            end do
         end do
         if (allocated(rotation)) deallocate (rotation)
         call rotation_matrix(kappa, rotation)
         call pic_gemm(coeff, rotation, coeff_new, beta=0.0_dp)
      end subroutine rotate
   end subroutine run_soscf_phase_uhf

   subroutine run_gdm_phase(mol, ops, ctrl, st, clk, xc, switch_iter, max_steps, &
                            result, error)
      !! Converge a closed-shell SCF by geometric direct minimisation, from where
      !! DIIS left it
      !!
      !! Van Voorhis and Head-Gordon, Mol. Phys. 100, 1713 (2002). The energy is
      !! minimised over the occupied-virtual rotations, `C <- C exp(kappa)`,
      !! along a limited-memory BFGS direction, and a step is kept only if the
      !! energy did not rise. One Fock build per trial and none besides, which
      !! is what separates it from the second-order phase: no Hessian-vector
      !! products, so an iteration costs what a DIIS one does.
      !!
      !! One trial per iteration:
      !!
      !! 1. Rotate the accepted orbitals by `scale * d`, rebuild the density,
      !!    the Fock matrix and the energy.
      !! 2. Accept if the energy is no higher than the accepted one, to
      !!    `GDM_SLACK`. Otherwise halve the scale and try the same direction
      !!    again from the accepted orbitals. Eleven halvings in a row drop the
      !!    history and restart along preconditioned steepest descent; if that
      !!    direction breaks down too, stop.
      !! 3. On acceptance the gradient `4 F_ai` is formed in the new orbitals'
      !!    frame, the pair (step taken, gradient change) joins the history, and
      !!    the next direction is `-H g`, replaced by steepest descent if it
      !!    does not point downhill and capped at `gdm_max_rotation` always.
      !!
      !! The gradient and the history live in the frame of the point they were
      !! formed at and are compared across successive frames without transport;
      !! the frames differ by one small rotation and the history's curvature
      !! test drops a pair the approximation spoils. For the same reason the
      !! orbitals are **not** semicanonicalised between steps -- that would
      !! mix them arbitrarily within each block and make the history
      !! meaningless -- only once, at the end.
      !!
      !! Convergence is the rule the DIIS phase used, `scf_convergence_t` on
      !! the energy change, the density change and the same commutator, tested
      !! at accepted points only.
      !!
      !! **It is a minimiser and not a saddle escape.** It keeps the occupation
      !! the handover diagonalisation gave it and cannot see negative curvature
      !! along a direction whose gradient is zero. `keywords.scf.second_order`
      !! is what leaves a saddle; `keywords.scf.stability` is what finds one.
      type(czt_molecule_t), intent(in), target :: mol
      type(rhf_operators_t), intent(in), target :: ops
      type(rhf_controls_t), intent(in) :: ctrl
      type(rhf_state_t), intent(inout) :: st
      type(timing_report_t), intent(inout) :: clk
      type(xc_context_t), intent(inout), optional, target :: xc
      integer, intent(in) :: switch_iter
         !! The DIIS iteration this took over after. Reported, not used.
      integer, intent(in) :: max_steps
         !! Trial steps allowed: whatever is left of `max_iter`.
      type(rhf_result_t), intent(inout) :: result
      type(error_t), intent(inout) :: error

      type(scf_gdm_history_t) :: history
      real(dp), allocatable :: fock_acc(:, :), fock_try(:, :), fock_mo(:, :), half(:, :)
      real(dp), allocatable :: coeff_try(:, :), density_try(:, :)
      real(dp), allocatable :: kappa(:, :), rotation(:, :), energies(:)
      real(dp), allocatable :: g_acc(:), g_new(:), h0_inv(:), direction(:)
      real(dp), allocatable :: step(:), change(:), diag(:)
      real(dp) :: e_acc, e_try, de, drms, gnorm, scale, scale_used, factor
      integer :: iteration, n_ao, n_mo, n_occ, n_vir, n_rot, i, a
      logical :: accepted, broken, restarted, pushed, valid

      if (error%has_error()) return
      n_ao = ops%n_ao
      n_mo = ops%n_mo
      n_occ = ops%n_occ
      n_vir = n_mo - n_occ
      n_rot = n_occ*n_vir
      result%gdm_started_at = switch_iter
      ! Nothing to rotate: the handover's orbitals are already the answer.
      if (max_steps < 1 .or. n_rot < 1) return

      allocate (fock_acc(n_ao, n_ao), fock_try(n_ao, n_ao), fock_mo(n_mo, n_mo))
      allocate (half(n_ao, n_mo), coeff_try(n_ao, n_mo), density_try(n_ao, n_ao))
      allocate (kappa(n_mo, n_mo), diag(n_mo))
      allocate (g_acc(n_rot), g_new(n_rot), h0_inv(n_rot), direction(n_rot))
      allocate (step(n_rot), change(n_rot))
      call history%init(n_rot, GDM_HISTORY, error)
      if (error%has_error()) return

      ! The point the minimisation starts from, built in full, for the reason
      ! the second-order phase gives: the DIIS loop's Fock matrix is
      ! extrapolated and its energy belongs to the previous density.
      call assemble_fock(mol, ops%h, st%density, st%coeff, n_occ, ops%bmat, &
                         ops%eri, ops%bounds, xc, fock_acc, e_acc, error, clk=clk, &
                         bmat_lr=ops%bmat_lr, jfit=ops%jfit)
      if (error%has_error()) return
      st%incr%full_builds = st%incr%full_builds + 1
      call gdm_gradient(fock_acc, st%coeff, g_acc)
      call gdm_steepest(g_acc, direction)
      restarted = .true.
      scale = 1.0_dp
      call gdm_table_header(ctrl%verbose)

      do iteration = 1, max_steps
         ! ---- converged at the accepted point? -----------------------------
         call commutator(fock_acc, st%density, ops%s, ops%x, st%err)
         gnorm = maxval(abs(st%err))
         result%commutator = gnorm
         de = abs(e_acc - st%e_old)
         drms = sqrt(sum((st%density - st%density_old)**2)/real(n_ao*n_ao, dp))
         if (ctrl%conv%is_converged(de, drms, gnorm)) then
            result%converged = .true.
            exit
         end if

         ! ---- one trial along the current direction -------------------------
         scale_used = scale
         kappa = 0.0_dp
         do i = 1, n_occ
            do a = 1, n_vir
               kappa(n_occ + a, i) = scale_used*direction(a + (i - 1)*n_vir)
               kappa(i, n_occ + a) = -kappa(n_occ + a, i)
            end do
         end do
         if (allocated(rotation)) deallocate (rotation)
         call rotation_matrix(kappa, rotation)
         call pic_gemm(st%coeff, rotation, coeff_try, beta=0.0_dp)
         call build_density_closed_shell(coeff_try, n_occ, density_try)
         ! In full, never incrementally: a rejected trial's density is one
         ! nothing comes back to, and the accept rule compares energies to
         ! 1e-10, below what an accumulated increment carries.
         call assemble_fock(mol, ops%h, density_try, coeff_try, n_occ, ops%bmat, &
                            ops%eri, ops%bounds, xc, fock_try, e_try, error, clk=clk, &
                            screening=st%screening, bmat_lr=ops%bmat_lr, jfit=ops%jfit)
         if (error%has_error()) return
         st%incr%full_builds = st%incr%full_builds + 1
         result%gdm_iterations = iteration
         result%iterations = switch_iter + iteration

         call scf_gdm_step_judge(e_acc, e_try, GDM_SLACK, scale, accepted, broken)
         call gdm_table_row(ctrl%verbose, iteration, e_try + mol%nuclear_repulsion(), &
                            e_try - e_acc, gnorm, scale_used, history%npairs(), accepted)
         if (accepted) then
            call gdm_gradient(fock_try, coeff_try, g_new)
            step = scale_used*direction
            change = g_new - g_acc
            call history%push(step, change, pushed, error)
            if (error%has_error()) return
            st%density_old = st%density
            st%density = density_try
            st%coeff = coeff_try
            fock_acc = fock_try
            st%e_old = e_acc
            e_acc = e_try
            g_acc = g_new
            call history%direction(g_acc, h0_inv, direction, error)
            if (error%has_error()) return
            restarted = history%npairs() == 0
            if (.not. scf_gdm_is_descent(g_acc, direction)) then
               call history%reset()
               call gdm_steepest(g_acc, direction)
               restarted = .true.
            else
               call gdm_cap(direction)
            end if
         else
            result%gdm_rejected = result%gdm_rejected + 1
            if (broken) then
               if (restarted) then
                  ! Steepest descent itself found nothing downhill. Not an
                  ! error: the energy and density in hand are the last ones
                  ! accepted, and `converged` is false.
                  result%gdm_stalled = .true.
                  call logger%warning("    direct minimisation found no step downhill; "// &
                                      "stopping. The reported energy is the last one "// &
                                      "accepted.")
                  exit
               end if
               call history%reset()
               call gdm_steepest(g_acc, direction)
               restarted = .true.
            end if
         end if
         call clk%lap(STAGE_GDM)
      end do

      ! Semicanonical orbitals of the last accepted point, as the second-order
      ! phase leaves them, and for its reason: the occupied-virtual block has
      ! vanished, so these are canonical at convergence. The density is left
      ! alone -- this rotation does not change it.
      call soscf_semicanonicalize(fock_acc, st%coeff, n_occ, fock_mo, energies, error)
      if (error%has_error()) return
      st%eigenvalues = energies
      result%gdm_off_aufbau = count(energies(n_occ + 1:n_mo) &
                                    < maxval(energies(1:n_occ)) - GDM_AUFBAU_TOL)
      if (result%gdm_off_aufbau > 0) then
         call logger%warning("    direct minimisation ended with a virtual orbital "// &
                             "below the highest occupied one: it keeps the occupation "// &
                             "it was handed and may have been started too early. "// &
                             "Lower keywords.scf.gdm_start.")
      end if
      result%full_fock_builds = st%incr%full_builds
      result%incremental_settled_at = st%incr%settled_at
      result%incremental_updates = st%incr%updates
      call history%release()

   contains

      subroutine gdm_gradient(fock_ao, coeff, gradient)
         !! `4 F_ai` in the frame of `coeff`, virtual fastest, and the
         !! preconditioner from the diagonal of the same matrix
         real(dp), intent(in) :: fock_ao(:, :), coeff(:, :)
         real(dp), intent(out) :: gradient(:)

         integer :: p

         call pic_gemm(fock_ao, coeff, half, beta=0.0_dp)
         call pic_gemm(coeff, half, fock_mo, transa="T", beta=0.0_dp)
         do i = 1, n_occ
            do a = 1, n_vir
               gradient(a + (i - 1)*n_vir) = 4.0_dp*fock_mo(n_occ + a, i)
            end do
         end do
         do p = 1, n_mo
            diag(p) = fock_mo(p, p)
         end do
         call scf_gdm_preconditioner_build(diag(1:n_occ), diag(n_occ + 1:n_mo), 4.0_dp, &
                                           GDM_HESSIAN_FLOOR, h0_inv, error)
      end subroutine gdm_gradient

      subroutine gdm_steepest(gradient, d)
         !! Preconditioned steepest descent, capped
         real(dp), intent(in) :: gradient(:)
         real(dp), intent(out) :: d(:)

         d = -h0_inv*gradient
         call gdm_cap(d)
      end subroutine gdm_steepest

      subroutine gdm_cap(d)
         !! Scale `d` so no rotation exceeds `gdm_max_rotation`
         real(dp), intent(inout) :: d(:)

         call scf_gdm_step_limit(d, ctrl%gdm_max_rotation, factor, valid)
         ! A direction that is not finite is replaced, not scaled.
         if (.not. valid) then
            call history%reset()
            d = -h0_inv*g_acc
            call scf_gdm_step_limit(d, ctrl%gdm_max_rotation, factor, valid)
         end if
         d = factor*d
      end subroutine gdm_cap
   end subroutine run_gdm_phase

   subroutine run_gdm_phase_uhf(mol, h, s, x, eri, bounds, n_mo, n_alpha, n_beta, conv, &
                                verbose, max_rotation, st, clk, xc, switch_iter, &
                                max_steps, result, error)
      !! Converge an unrestricted SCF by geometric direct minimisation, from
      !! where DIIS left it
      !!
      !! `run_gdm_phase` over two spins. One vector holds both rotation spaces,
      !! alpha's occupied-virtual pairs and then beta's, virtual fastest in
      !! each, with one history and one accept rule over it because there is
      !! one energy. The gradient is `2 F_ai` per spin and the preconditioner's
      !! factor is two with it. A spin with no occupied or no virtual orbital
      !! has nothing to rotate and contributes no elements.
      !!
      !! Everything else -- the trial, the halving, the restart along steepest
      !! descent, the untransported history, the single semicanonicalisation
      !! at the end -- is as `run_gdm_phase` describes.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: h(:, :), s(:, :), x(:, :)
      real(dp), intent(in), allocatable :: eri(:, :, :, :)
      real(dp), intent(in), allocatable :: bounds(:, :)
      integer, intent(in) :: n_mo, n_alpha, n_beta
      type(scf_convergence_t), intent(in) :: conv
      logical, intent(in) :: verbose
      real(dp), intent(in) :: max_rotation   !! Largest rotation of one step, radians
      type(uhf_state_t), intent(inout) :: st
      type(timing_report_t), intent(inout) :: clk
      type(xc_context_t), intent(inout), optional :: xc
      integer, intent(in) :: switch_iter
         !! The DIIS iteration this took over after. Reported, not used.
      integer, intent(in) :: max_steps
         !! Trial steps allowed: whatever is left of `max_iter`.
      type(rhf_result_t), intent(inout) :: result
      type(error_t), intent(inout) :: error

      real(dp), parameter :: SPIN_FACTOR = 2.0_dp
         !! Between `dE/dkappa_ai` and `F_ai` for one spin
      type(scf_gdm_history_t) :: history
      real(dp), allocatable :: fa_acc(:, :), fb_acc(:, :), fa_try(:, :), fb_try(:, :)
      real(dp), allocatable :: ca_try(:, :), cb_try(:, :), da_try(:, :), db_try(:, :)
      real(dp), allocatable :: fock_mo(:, :), half(:, :), kappa(:, :), rotation(:, :)
      real(dp), allocatable :: energies(:), diag(:)
      real(dp), allocatable :: g_acc(:), g_new(:), h0_inv(:), direction(:)
      real(dp), allocatable :: step(:), change(:)
      real(dp) :: e_acc, e_try, de, drms, gnorm, scale, scale_used, factor
      integer :: iteration, n_ao, n_rot_a, n_rot_b, n_rot
      logical :: accepted, broken, restarted, pushed, valid

      if (error%has_error()) return
      n_ao = size(h, 1)
      n_rot_a = n_alpha*(n_mo - n_alpha)
      n_rot_b = n_beta*(n_mo - n_beta)
      n_rot = n_rot_a + n_rot_b
      result%gdm_started_at = switch_iter
      ! Nothing to rotate: the handover's orbitals are already the answer.
      if (max_steps < 1 .or. n_rot < 1) return

      allocate (fa_acc(n_ao, n_ao), fb_acc(n_ao, n_ao), fa_try(n_ao, n_ao), fb_try(n_ao, n_ao))
      allocate (ca_try(n_ao, n_mo), cb_try(n_ao, n_mo), da_try(n_ao, n_ao), db_try(n_ao, n_ao))
      allocate (fock_mo(n_mo, n_mo), half(n_ao, n_mo), kappa(n_mo, n_mo), diag(n_mo))
      allocate (g_acc(n_rot), g_new(n_rot), h0_inv(n_rot), direction(n_rot))
      allocate (step(n_rot), change(n_rot))
      call history%init(n_rot, GDM_HISTORY, error)
      if (error%has_error()) return

      ! The point the minimisation starts from, built in full: the DIIS
      ! loop's Fock matrices are extrapolated and its energy belongs to the
      ! previous densities.
      call assemble_fock_uhf(mol, h, st%d_a, st%d_b, eri, bounds, xc, fa_acc, fb_acc, &
                             e_acc, error, clk=clk, bmat=st%bmat, bmat_lr=st%bmat_lr)
      if (error%has_error()) return
      st%incr%full_builds = st%incr%full_builds + 1
      call gdm_gradient(fa_acc, fb_acc, st%coeff_a, st%coeff_b, g_acc)
      if (error%has_error()) return
      call gdm_steepest(g_acc, direction)
      restarted = .true.
      scale = 1.0_dp
      call gdm_table_header(verbose)

      do iteration = 1, max_steps
         ! ---- converged at the accepted point? -----------------------------
         call commutator(fa_acc, st%d_a, s, x, st%err_a)
         call commutator(fb_acc, st%d_b, s, x, st%err_b)
         gnorm = max(maxval(abs(st%err_a)), maxval(abs(st%err_b)))
         result%commutator = gnorm
         de = abs(e_acc - st%e_old)
         drms = sqrt((sum((st%d_a - st%d_a_old)**2) + sum((st%d_b - st%d_b_old)**2)) &
                     /real(2*n_ao*n_ao, dp))
         if (conv%is_converged(de, drms, gnorm)) then
            result%converged = .true.
            exit
         end if

         ! ---- one trial along the current direction -------------------------
         scale_used = scale
         call gdm_rotate(st%coeff_a, n_alpha, direction(1:n_rot_a), ca_try)
         call gdm_rotate(st%coeff_b, n_beta, direction(n_rot_a + 1:n_rot), cb_try)
         call build_density_spin(ca_try, n_alpha, da_try)
         call build_density_spin(cb_try, n_beta, db_try)
         ! In full, never incrementally, for the reason `run_gdm_phase` gives.
         call assemble_fock_uhf(mol, h, da_try, db_try, eri, bounds, xc, fa_try, fb_try, &
                                e_try, error, clk=clk, bmat=st%bmat, bmat_lr=st%bmat_lr)
         if (error%has_error()) return
         st%incr%full_builds = st%incr%full_builds + 1
         result%gdm_iterations = iteration
         result%iterations = switch_iter + iteration

         call scf_gdm_step_judge(e_acc, e_try, GDM_SLACK, scale, accepted, broken)
         call gdm_table_row(verbose, iteration, e_try + mol%nuclear_repulsion(), &
                            e_try - e_acc, gnorm, scale_used, history%npairs(), accepted)
         if (accepted) then
            call gdm_gradient(fa_try, fb_try, ca_try, cb_try, g_new)
            if (error%has_error()) return
            step = scale_used*direction
            change = g_new - g_acc
            call history%push(step, change, pushed, error)
            if (error%has_error()) return
            st%d_a_old = st%d_a
            st%d_b_old = st%d_b
            st%d_a = da_try
            st%d_b = db_try
            st%coeff_a = ca_try
            st%coeff_b = cb_try
            fa_acc = fa_try
            fb_acc = fb_try
            st%e_old = e_acc
            e_acc = e_try
            g_acc = g_new
            call history%direction(g_acc, h0_inv, direction, error)
            if (error%has_error()) return
            restarted = history%npairs() == 0
            if (.not. scf_gdm_is_descent(g_acc, direction)) then
               call history%reset()
               call gdm_steepest(g_acc, direction)
               restarted = .true.
            else
               call gdm_cap(direction)
            end if
         else
            result%gdm_rejected = result%gdm_rejected + 1
            if (broken) then
               if (restarted) then
                  ! Steepest descent itself found nothing downhill. Not an
                  ! error: the energy and densities in hand are the last ones
                  ! accepted, and `converged` is false.
                  result%gdm_stalled = .true.
                  call logger%warning("    direct minimisation found no step downhill; "// &
                                      "stopping. The reported energy is the last one "// &
                                      "accepted.")
                  exit
               end if
               call history%reset()
               call gdm_steepest(g_acc, direction)
               restarted = .true.
            end if
         end if
         call clk%lap(STAGE_GDM)
      end do

      ! Semicanonical orbitals of the last accepted point, per spin. The
      ! densities are left alone -- these rotations do not change them.
      result%gdm_off_aufbau = 0
      call gdm_finish(fa_acc, st%coeff_a, n_alpha, st%eig_a)
      if (error%has_error()) return
      call gdm_finish(fb_acc, st%coeff_b, n_beta, st%eig_b)
      if (error%has_error()) return
      if (result%gdm_off_aufbau > 0) then
         call logger%warning("    direct minimisation ended with a virtual orbital "// &
                             "below the highest occupied one of its spin: it keeps the "// &
                             "occupation it was handed and may have been started too "// &
                             "early. Lower keywords.scf.gdm_start.")
      end if
      st%fock_a = fa_acc
      st%fock_b = fb_acc
      st%e_elec = e_acc
      result%full_fock_builds = st%incr%full_builds
      call history%release()

   contains

      subroutine gdm_gradient(fock_a, fock_b, coeff_a, coeff_b, gradient)
         !! `2 F_ai` of each spin in the frame of its own orbitals, alpha then
         !! beta, and the preconditioner from the diagonals of the same
         !! matrices
         real(dp), intent(in) :: fock_a(:, :), fock_b(:, :), coeff_a(:, :), coeff_b(:, :)
         real(dp), intent(out) :: gradient(:)

         call spin_gradient(fock_a, coeff_a, n_alpha, gradient(1:n_rot_a), h0_inv(1:n_rot_a))
         if (error%has_error()) return
         call spin_gradient(fock_b, coeff_b, n_beta, gradient(n_rot_a + 1:n_rot), &
                            h0_inv(n_rot_a + 1:n_rot))
      end subroutine gdm_gradient

      subroutine spin_gradient(fock_ao, coeff, n_occ, gradient, inverse)
         !! One spin's block of the gradient and of the preconditioner
         real(dp), intent(in) :: fock_ao(:, :), coeff(:, :)
         integer, intent(in) :: n_occ
         real(dp), intent(out) :: gradient(:), inverse(:)

         integer :: i, a, p, n_vir

         n_vir = n_mo - n_occ
         if (n_occ*n_vir == 0) return
         call pic_gemm(fock_ao, coeff, half, beta=0.0_dp)
         call pic_gemm(coeff, half, fock_mo, transa="T", beta=0.0_dp)
         do i = 1, n_occ
            do a = 1, n_vir
               gradient(a + (i - 1)*n_vir) = SPIN_FACTOR*fock_mo(n_occ + a, i)
            end do
         end do
         do p = 1, n_mo
            diag(p) = fock_mo(p, p)
         end do
         call scf_gdm_preconditioner_build(diag(1:n_occ), diag(n_occ + 1:n_mo), SPIN_FACTOR, &
                                           GDM_HESSIAN_FLOOR, inverse, error)
      end subroutine spin_gradient

      subroutine gdm_rotate(coeff, n_occ, d, coeff_new)
         !! `coeff exp(scale_used kappa(d))` for one spin
         real(dp), intent(in) :: coeff(:, :), d(:)
         integer, intent(in) :: n_occ
         real(dp), intent(out) :: coeff_new(:, :)

         integer :: i, a, n_vir

         n_vir = n_mo - n_occ
         if (n_occ*n_vir == 0) then
            coeff_new = coeff
            return
         end if
         kappa = 0.0_dp
         do i = 1, n_occ
            do a = 1, n_vir
               kappa(n_occ + a, i) = scale_used*d(a + (i - 1)*n_vir)
               kappa(i, n_occ + a) = -kappa(n_occ + a, i)
            end do
         end do
         if (allocated(rotation)) deallocate (rotation)
         call rotation_matrix(kappa, rotation)
         call pic_gemm(coeff, rotation, coeff_new, beta=0.0_dp)
      end subroutine gdm_rotate

      subroutine gdm_finish(fock_ao, coeff, n_occ, eigenvalues)
         !! Semicanonicalise one spin and count its virtuals below its HOMO
         real(dp), intent(in) :: fock_ao(:, :)
         real(dp), allocatable, intent(inout) :: coeff(:, :)
         integer, intent(in) :: n_occ
         real(dp), allocatable, intent(inout) :: eigenvalues(:)

         real(dp), allocatable :: semi(:, :)

         ! An empty spin has no occupied block to separate from the virtual
         ! one; its orbitals are simply this Fock matrix's eigenvectors.
         if (n_occ == 0) then
            call diagonalize(fock_ao, x, n_ao, n_mo, coeff, eigenvalues, error)
            return
         end if
         call soscf_semicanonicalize(fock_ao, coeff, n_occ, semi, energies, error)
         if (error%has_error()) return
         eigenvalues = energies
         if (n_occ < n_mo) then
            result%gdm_off_aufbau = result%gdm_off_aufbau &
                                    + count(energies(n_occ + 1:n_mo) &
                                            < maxval(energies(1:n_occ)) - GDM_AUFBAU_TOL)
         end if
      end subroutine gdm_finish

      subroutine gdm_steepest(gradient, d)
         !! Preconditioned steepest descent, capped
         real(dp), intent(in) :: gradient(:)
         real(dp), intent(out) :: d(:)

         d = -h0_inv*gradient
         call gdm_cap(d)
      end subroutine gdm_steepest

      subroutine gdm_cap(d)
         !! Scale `d` so no rotation exceeds `max_rotation`
         real(dp), intent(inout) :: d(:)

         call scf_gdm_step_limit(d, max_rotation, factor, valid)
         ! A direction that is not finite is replaced, not scaled.
         if (.not. valid) then
            call history%reset()
            d = -h0_inv*g_acc
            call scf_gdm_step_limit(d, max_rotation, factor, valid)
         end if
         d = factor*d
      end subroutine gdm_cap
   end subroutine run_gdm_phase_uhf

   subroutine gdm_table_header(verbose)
      !! Column headings for the direct-minimisation phase
      logical, intent(in) :: verbose

      if (.not. verbose) return
      call logger%info("  direct minimisation: L-BFGS on the orbital rotations, "// &
                       "one Fock build per step")
      call logger%info("    iter                 energy          dE     |FDS-SDF|"// &
                       "       scale  pairs")
   end subroutine gdm_table_header

   subroutine gdm_table_row(verbose, iter, energy, de, gnorm, scale, pairs, accepted)
      !! One line of the direct-minimisation phase
      !!
      !! `energy` and `dE` are the trial's, against the accepted point it
      !! started from; `|FDS-SDF|` is that accepted point's. A rejected trial
      !! is marked, and the next line retries the same direction at half the
      !! scale.
      logical, intent(in) :: verbose, accepted
      integer, intent(in) :: iter, pairs
      real(dp), intent(in) :: energy, de, gnorm, scale

      character(len=LINE_LEN) :: line

      if (.not. verbose) return
      if (accepted) then
         write (line, "(i8,f22.12,3es12.3,i7)") iter, energy, de, gnorm, scale, pairs
      else
         write (line, "(i8,f22.12,3es12.3,i7,a)") iter, energy, de, gnorm, scale, pairs, &
            "  rejected"
      end if
      call logger%info(trim(line))
   end subroutine gdm_table_row

   subroutine build_soscf_hessian(mol, ops, ctrl, st, energies, response, hessian, &
                                  xc, error)
      !! Point an `ov_hessian_t` at the orbitals the Newton step is being taken from
      !!
      !! The four-way branch is on two independent facts, neither of which can
      !! be passed as an absent optional from inside a single call: whether the
      !! reference is Kohn-Sham, and whether this SCF has Schwarz bounds
      !! already. The bounds are only built on the direct path, and the Hessian
      !! needs them either way -- `build_scf_ov_hessian` builds its own when
      !! they are not handed over.
      !!
      !! **The Hessian is always built from exact integrals, direct**, whatever
      !! the reference was built from: `nuclear_apply` goes through
      !! `build_fock_direct_many`. For an in-core reference that is the same
      !! operator. For a density-fitted one it is not -- the Hessian is then the
      !! curvature of a slightly different energy surface than the one being
      !! minimised. That degrades the *rate*, not the answer: the gradient and
      !! the energy the step is accepted on are the fitted ones, so the
      !! stationary point reached is the fitted SCF's.
      type(czt_molecule_t), intent(in), target :: mol
      type(rhf_operators_t), intent(in) :: ops
      type(rhf_controls_t), intent(in) :: ctrl
      type(rhf_state_t), intent(in) :: st
      real(dp), intent(in) :: energies(:)
      type(nuclear_response_t), intent(out), target :: response
      type(ov_hessian_t), intent(out) :: hessian
      type(xc_context_t), intent(inout), optional, target :: xc
      type(error_t), intent(inout) :: error

      logical :: have_bounds

      if (error%has_error()) return
      have_bounds = allocated(ops%bounds)

      if (ctrl%kohn_sham .and. present(xc)) then
         if (have_bounds) then
            call build_scf_ov_hessian(mol, st%coeff, energies, ops%n_occ, response, &
                                      hessian, error, xc=xc, reference=st%density, &
                                      k_scale=xc%exx_fraction, rs_k_lr=xc%rs_k_lr, &
                                      rs_omega=xc%rs_omega, bounds=ops%bounds)
         else
            call build_scf_ov_hessian(mol, st%coeff, energies, ops%n_occ, response, &
                                      hessian, error, xc=xc, reference=st%density, &
                                      k_scale=xc%exx_fraction, rs_k_lr=xc%rs_k_lr, &
                                      rs_omega=xc%rs_omega)
         end if
      else
         if (have_bounds) then
            call build_scf_ov_hessian(mol, st%coeff, energies, ops%n_occ, response, &
                                      hessian, error, bounds=ops%bounds)
         else
            call build_scf_ov_hessian(mol, st%coeff, energies, ops%n_occ, response, &
                                      hessian, error)
         end if
      end if
      if (error%has_error()) then
         call error%add_context("building the orbital-rotation Hessian for the "// &
                                "second-order SCF")
      end if
   end subroutine build_soscf_hessian

   subroutine run_czt_uhf(mol, nelec, multiplicity, max_iter, energy_tol, density_tol, &
                          verbose, result, error, diis_vectors, in_core, diis_start, &
                          guess, guess_density_alpha, guess_density_beta, xc, pcm, &
                          level_shift, linear_dependence, accelerator, grad_tol, &
                          incremental_fock, convergence, scf, aux, gdm, gdm_start, &
                          gdm_max_rotation, second_order, soscf_start, restricted_open)
      !! Drive an unrestricted SCF to convergence
      !!
      !! Two Fock matrices sharing one Coulomb field: F_sigma = H + J(D_a + D_b)
      !! - K(D_sigma). They are diagonalised separately and extrapolated
      !! together, one DIIS over the pair.
      type(czt_molecule_t), intent(in) :: mol
      type(czt_molecule_t), intent(in), optional :: aux
         !! Auxiliary basis. Present fits both J and K, as on the restricted
         !! path, and excludes `in_core`. Always the dense tensor: the
         !! Coulomb-only fit a pure functional takes there is not used here.
      logical, intent(in), optional :: gdm
         !! `keywords.scf.gdm`. Hand over from DIIS to geometric direct
         !! minimisation once the commutator falls below `gdm_start`; see
         !! `run_gdm_phase_uhf`. Default false, and refused with a continuum
         !! solvent.
      real(dp), intent(in), optional :: gdm_start
         !! The commutator at which that happens. Absent or not positive
         !! leaves the default.
      real(dp), intent(in), optional :: gdm_max_rotation
         !! The largest rotation one step may make, in radians. Absent or not
         !! positive leaves the default.
      logical, intent(in), optional :: second_order
         !! `keywords.scf.second_order`. Hand over from DIIS to trust-region
         !! Newton on both spins' rotations once the commutator falls below
         !! `soscf_start`; see `run_soscf_phase_uhf`. Default false. Refused
         !! with `gdm`, a continuum solvent, or a fitted or stored-integral
         !! build.
      real(dp), intent(in), optional :: soscf_start
         !! The commutator at which that happens. Absent or not positive
         !! leaves the default.
      logical, intent(in), optional :: restricted_open
         !! Restricted open-shell Hartree-Fock or Kohn-Sham: the two spins
         !! share their orbital spaces, alpha filling `n_alpha` orbitals and
         !! beta the first `n_beta` of them. The energy is the unrestricted
         !! expression at those densities. The result is shaped as an
         !! unrestricted one -- two densities, two sets of orbitals and
         !! energies -- with each spin's orbitals diagonalising that spin's
         !! own Fock matrix within its occupied and within its virtual block.
         !! The two sets span the same spaces and are not the same vectors.
      integer, intent(in) :: nelec
      integer, intent(in) :: multiplicity   !! 2S+1
      integer, intent(in) :: max_iter
      real(dp), intent(in) :: energy_tol, density_tol
      real(dp), intent(in), optional :: grad_tol
         !! Commutator threshold -- pyscf's `conv_tol_grad`. Zero or absent
         !! derives `sqrt(energy_tol)`, which is the right scaling *for the
         !! energy*: its error goes as `gnorm**2`.
         !!
         !! **A caller that consumes the density has to set this.** Density
         !! error goes as `gnorm`, not its square, so the derived value bounds
         !! it about three orders more loosely than the energy, and tightening
         !! `energy_tol` cannot substitute: reaching `gnorm < 1e-8` that way
         !! needs `energy_tol = 1e-16`, below what a molecular energy resolves.
      type(scf_convergence_t), intent(in), optional :: convergence
         !! The rule that decides this SCF has converged. Absent reproduces
         !! what `energy_tol` and `grad_tol` meant before this existed.
      type(scf_numerics_t), intent(in), optional :: scf
         !! The whole configuration, as one argument.
         !!
         !! **Eight of its thirteen fields are honoured here**, and passing
         !! this looks like passing everything. Read: `level_shift`,
         !! `linear_dependence`, `accelerator`, `diis_size`, `use_diis`,
         !! `incremental_fock`, `convergence_metric` and `grad_tol`. Not read:
         !!
         !! * `max_iter`, `energy_tol` and `density_tol`, which are positional
         !!   arguments here -- whichever the caller passes positionally runs.
         !! * `guess`, a deck spelling whose translation to an `SCF_GUESS_*`
         !!   kind needs `mqc_czt_atomic_guess`, which calls this routine.
         !!   The caller parses it and passes `guess=`.
         !! * `allow_crap_scf`. This routine reports `result%converged` and the
         !!   caller decides whether an unconverged SCF is acceptable.
         !!
         !! An individual optional above wins over the group's copy of the same
         !! setting, so a caller can pass the group and override one field.
      logical, intent(in), optional :: incremental_fock
         !! Build each iteration from the density change. Default true; false
         !! forces a full build every iteration, which is what to reach for when
         !! an SCF stalls or when a Fock timing needs to mean something.
      logical, intent(in) :: verbose
      real(dp), intent(in), optional :: level_shift
         !! Hartree added to the virtual block of each spin's Fock matrix
         !! before its diagonalisation, against that spin's own virtual
         !! projector. Absent or zero is off; tapered off before convergence as
         !! on the restricted path.
      integer, intent(in), optional :: accelerator
         !! `ACCEL_DIIS`, the default, or one of the energy-based schemes to run
         !! while the error is large.
      type(rhf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: diis_vectors
      logical, intent(in), optional :: in_core
      integer, intent(in), optional :: diis_start
         !! First iteration allowed to extrapolate. See the default below.
      integer, intent(in), optional :: guess
         !! One of SCF_GUESS_*. Defaults to Wolfsberg-Helmholz, as the restricted
         !! path does.
      type(xc_context_t), intent(inout), optional :: xc
         !! Exchange-correlation, spin-polarised. Present makes this an
         !! unrestricted Kohn-Sham SCF; absent leaves it unrestricted
         !! Hartree-Fock. **The context must have been built with
         !! `polarized = .true.`** -- libxc fixes the spin channel at
         !! initialisation, and `xc_add_potential_uks` checks rather than
         !! assumes it.
      real(dp), intent(in), optional :: guess_density_alpha(:, :)
      real(dp), intent(in), optional :: guess_density_beta(:, :)
         !! Starting spin densities, required by SCF_GUESS_SAC, SCF_GUESS_SAD
         !! and SCF_GUESS_MINAO.
         !! The two are taken separately because whether they differ is the whole
         !! difference between those guesses: SAC hands over the free atom's own
         !! spin-polarised densities and so arrives already symmetry-broken,
         !! where SAD hands over half the spherically averaged total twice and
         !! relies on the occupations to break the symmetry.
      type(pcm_context_t), intent(inout), optional :: pcm
         !! Continuum solvation, exactly as on the restricted path. The surface
         !! charges see the total density and their operator lands in both spin
         !! Fock matrices, because the potential of a classical charge does not
         !! know spin.
      real(dp), intent(in), optional :: linear_dependence
         !! `keywords.scf.linear_dependence_threshold`. Zero or absent leaves
         !! the orthogonaliser on its own cutoff.

      integer :: diis_size, start_cycle, guess_kind
      logical :: use_in_core
      real(dp), allocatable :: bounds(:, :)
      ! TODO(mqc): `stats` is declared here and never used; the screening counts
      ! come back in `screening`, which `assemble_fock` fills.
      type(direct_stats_t) :: stats

      real(dp), allocatable :: s(:, :), h(:, :), eri(:, :, :, :)
      real(dp), allocatable :: x(:, :)
      type(uhf_state_t) :: st
      logical :: use_pcm
      real(dp) :: e_pcm
      logical :: extrapolated
      logical :: kohn_sham_run
      integer :: accel
      real(dp) :: de, drms
      real(dp) :: gnorm
      real(dp) :: gtol
      type(scf_convergence_t) :: conv
      logical :: accel_ok_grp
      integer :: metric_kind
      logical :: metric_ok
      real(dp) :: shift, shift_now, taper
      integer :: n_ao, n_mo, n_alpha, n_beta, iter, nsq, msq
      type(timing_report_t) :: clk
      real(dp) :: xc_discard(4)
      logical :: want_incremental
         !! Owned by the loop and handed only to the loop's own build. The guess
         !! and the final rebuild go without it, so the energy that leaves this
         !! routine is an exact build.
      real(dp) :: t_fock_iter, t_rest_iter, t_xc_iter
      real(dp) :: s_min, s_kept   !! overlap conditioning, reported before iteration 1
      real(dp) :: lindep
      character(len=LINE_LEN) :: line
      logical :: use_gdm, accel_gdm, use_soscf, accel_soscf
      real(dp) :: gdm_at, gdm_cap_rad, soscf_at
      integer :: gdm_switch, soscf_switch

      diis_size = 8
      if (present(scf)) then
         diis_size = scf%diis_size
         if (.not. scf%use_diis) diis_size = 0
      end if
      if (present(diis_vectors)) diis_size = diis_vectors
      start_cycle = DEFAULT_UHF_DIIS_START
      if (present(diis_start)) start_cycle = diis_start
      use_in_core = .false.
      want_incremental = .true.
      if (present(scf)) want_incremental = scf%incremental_fock
      if (present(incremental_fock)) want_incremental = incremental_fock
      if (present(in_core)) use_in_core = in_core
      guess_kind = SCF_GUESS_GWH
      if (present(guess)) guess_kind = guess
      use_pcm = .false.
      if (present(pcm)) use_pcm = pcm%enabled

      ! 2S+1 = n_alpha - n_beta + 1, so the unpaired count is multiplicity - 1.
      ! Both have to come out whole and non-negative; a deck can ask for a
      ! pairing that does not exist, four electrons in a quintet say.
      if (multiplicity < 1) then
         call error%set(ERROR_VALIDATION, "UHF: multiplicity must be at least 1")
         return
      end if
      if (mod(nelec + multiplicity - 1, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "UHF: an electron count and multiplicity "// &
                        "that cannot be paired -- their parities disagree")
         return
      end if
      n_alpha = (nelec + multiplicity - 1)/2
      n_beta = nelec - n_alpha
      if (n_beta < 0 .or. n_alpha < 0) then
         call error%set(ERROR_VALIDATION, "UHF: multiplicity asks for more unpaired "// &
                        "electrons than the system has")
         return
      end if
      if (n_alpha < 1) then
         call error%set(ERROR_VALIDATION, "UHF: no electrons to place")
         return
      end if

      n_ao = mol%nao
      if (verbose) then
         if (mol%cartesian) then
            call logger%info("  basis functions: "//to_char(n_ao)//"  (Cartesian, 6d/10f)")
         else
            call logger%info("  basis functions: "//to_char(n_ao)//"  (spherical, 5d/7f)")
         end if
         call logger%info("  "//reference_label()//": charge = "// &
                          to_char(nint(sum(mol%charges)) - nelec)// &
                          "  multiplicity = "//to_char(multiplicity)// &
                          "  electrons = "//to_char(nelec))
         call logger%info("  "//reference_label()//": n_alpha = "//to_char(n_alpha)// &
                          "  n_beta = "//to_char(n_beta))
      end if

      call clk%start()
      call xc_profile_take(xc_discard)   ! as in run_czt_rhf
      call mol%overlap(s)
      call mol%core_hamiltonian(h)
      if (present(aux)) then
         if (use_in_core) then
            call error%set(ERROR_VALIDATION, "UHF: in_core and density fitting were both "// &
                           "asked for; the fitted build never forms the four-index "// &
                           "integrals, so choose one")
            return
         end if
         call build_df_tensor(mol, aux, st%bmat, error)
         if (error%has_error()) return
         if (present(xc)) then
            if (xc%range_separated) then
               call build_df_tensor(mol, aux, st%bmat_lr, error, omega=xc%rs_omega)
               if (error%has_error()) return
            end if
         end if
      else if (use_in_core) then
         call mol%eris(eri)
      else
         call schwarz_bounds(mol, bounds, error)
         if (error%has_error()) return
      end if

      lindep = 0.0_dp
      if (present(scf)) lindep = scf%linear_dependence
      if (present(linear_dependence)) lindep = linear_dependence
      call build_orthogonalizer(s, x, n_mo, error, smallest_overlap=s_min, &
                                smallest_kept=s_kept, threshold=lindep)
      if (error%has_error()) return
      call report_linear_dependence(n_ao, n_mo, s_min, s_kept, verbose, threshold=lindep)
      if (n_alpha > n_mo) then
         call error%set(ERROR_VALIDATION, "UHF: more alpha electrons than the basis "// &
                        "supports after near-null modes were dropped")
         return
      end if

      nsq = n_ao*n_ao
      msq = n_mo*n_mo
      allocate (st%fock_a(n_ao, n_ao), st%fock_b(n_ao, n_ao))
      allocate (st%d_a(n_ao, n_ao), st%d_b(n_ao, n_ao), st%d_a_old(n_ao, n_ao), st%d_b_old(n_ao, n_ao))
      allocate (st%err_a(n_mo, n_mo), st%err_b(n_mo, n_mo))
      allocate (st%fock_flat(2*nsq), st%err_flat(2*msq), st%dens_flat(2*nsq))
      if (use_pcm) allocate (st%v_pcm(n_ao, n_ao))

      ! One subspace over both spins, so an extrapolation step moves them
      ! together. The vectors are the two Fock matrices laid end to end and the
      ! two commutators likewise.
      ! Resolved before the subspace is built, not after it. `energy_based`
      ! decides whether the history can hold densities and energies at all, and
      ! reading `accel` here while it was still assigned some fifty lines below
      ! took that decision on an undefined value.
      accel = ACCEL_DIIS
      accel_gdm = .false.
      accel_soscf = .false.
      if (present(scf)) call parse_accelerator_name(scf%accelerator, accel, accel_ok_grp, &
                                                    accel_soscf, accel_gdm)
      if (present(accelerator)) accel = accelerator
      use_soscf = accel_soscf
      if (present(scf)) use_soscf = use_soscf .or. scf%second_order
      if (present(second_order)) use_soscf = second_order .or. accel_soscf
      soscf_at = DEFAULT_SOSCF_START
      if (present(scf)) then
         if (scf%soscf_start > 0.0_dp) soscf_at = scf%soscf_start
      end if
      if (present(soscf_start)) then
         if (soscf_start > 0.0_dp) soscf_at = soscf_start
      end if
      if (use_soscf) then
         if (use_pcm) then
            call error%set(ERROR_VALIDATION, "keywords.scf.second_order with a "// &
                           "continuum solvent: the orbital Hessian carries no "// &
                           "surface-charge response. Converge this one by DIIS.")
            return
         end if
         if (allocated(st%bmat) .or. use_in_core) then
            call error%set(ERROR_VALIDATION, "keywords.scf.second_order over an "// &
                           "unrestricted reference needs the direct Fock build: its "// &
                           "Hessian-vector products are direct, and over a fitted or "// &
                           "stored-integral SCF they would belong to a different "// &
                           "energy. Run without density_fitting and in_core.")
            return
         end if
      end if
      ! `accelerator: gdm` and `gdm: true` are the same request, as on the
      ! restricted path.
      use_gdm = accel_gdm
      if (present(scf)) use_gdm = use_gdm .or. scf%gdm
      if (present(gdm)) use_gdm = gdm .or. accel_gdm
      gdm_at = DEFAULT_GDM_START
      gdm_cap_rad = DEFAULT_GDM_MAX_ROTATION
      if (present(scf)) then
         if (scf%gdm_start > 0.0_dp) gdm_at = scf%gdm_start
         if (scf%gdm_max_rotation > 0.0_dp) gdm_cap_rad = scf%gdm_max_rotation
      end if
      if (present(gdm_start)) then
         if (gdm_start > 0.0_dp) gdm_at = gdm_start
      end if
      if (present(gdm_max_rotation)) then
         if (gdm_max_rotation > 0.0_dp) gdm_cap_rad = gdm_max_rotation
      end if
      if (use_gdm .and. use_pcm) then
         call error%set(ERROR_VALIDATION, "keywords.scf.gdm with a continuum "// &
                        "solvent: the minimisation rebuilds the Fock matrix "// &
                        "itself, and its builds do not carry the reaction field. "// &
                        "Converge this one by DIIS.")
         return
      end if
      if (use_gdm .and. use_soscf) then
         call error%set(ERROR_VALIDATION, "keywords.scf.gdm together with "// &
                        "keywords.scf.second_order: each finishes the SCF that "// &
                        "DIIS started, and only one of them can. Choose one.")
         return
      end if
      st%rohf = .false.
      if (present(restricted_open)) st%rohf = restricted_open
      if (st%rohf) then
         ! Both finish an SCF over independent alpha and beta rotations, which
         ! is the freedom a restricted open shell does not have.
         if (use_gdm .or. use_soscf) then
            call error%set(ERROR_VALIDATION, "model.restricted_open with "// &
                           "keywords.scf.gdm or keywords.scf.second_order: both "// &
                           "rotate the alpha and beta orbitals independently, and a "// &
                           "restricted open shell has one set. Converge it by DIIS.")
            return
         end if
         ! The energy-based schemes interpolate Fock matrices against the
         ! energy they are the gradient of; Roothaan's effective matrix is the
         ! gradient of nothing.
         if (accel /= ACCEL_DIIS) then
            call error%set(ERROR_VALIDATION, "model.restricted_open with an "// &
                           "energy-based SCF accelerator: the matrix a restricted "// &
                           "open shell diagonalises is not the derivative of its "// &
                           "energy. Use keywords.scf.accelerator: diis.")
            return
         end if
      end if
      call st%diis%init(diis_size, 2*nsq, 2*msq, &
                        energy_based=(accel /= ACCEL_DIIS))

      ! The symmetric guesses -- core and GWH -- give alpha and beta the same
      ! orbitals, and the occupations separate them: n_alpha > n_beta puts an
      ! electron where beta has none. For n_alpha == n_beta they converge to the
      ! restricted solution, which is the right answer rather than a failure.
      ! SAC is the one that arrives already asymmetric.
      select case (guess_kind)
      case (SCF_GUESS_CORE)
         st%fock_a = h
         st%fock_b = h
      case (SCF_GUESS_GWH)
         call guess_fock(s, h, st%fock_a)
         st%fock_b = st%fock_a
      case (SCF_GUESS_SAP)
         call sap_fock(mol, h, st%fock_a, error)
         if (error%has_error()) return
         st%fock_b = st%fock_a
      case (SCF_GUESS_SAC, SCF_GUESS_SAD, SCF_GUESS_PROJ, SCF_GUESS_MINAO)
         if (.not. (present(guess_density_alpha) .and. present(guess_density_beta))) then
            call error%set(ERROR_VALIDATION, "UHF: an atomic guess was asked for but no "// &
                           "guess densities were supplied")
            return
         end if
         if (size(guess_density_alpha, 1) /= n_ao .or. size(guess_density_beta, 1) /= n_ao) then
            call error%set(ERROR_VALIDATION, "UHF: the guess densities are not the size "// &
                           "of this basis")
            return
         end if
         st%d_a = guess_density_alpha
         st%d_b = guess_density_beta
         call assemble_fock_uhf(mol, h, st%d_a, st%d_b, eri, bounds, xc, st%fock_a, st%fock_b, &
                                st%e_elec, error, bmat=st%bmat, bmat_lr=st%bmat_lr)
         if (error%has_error()) return
      case default
         call error%set(ERROR_VALIDATION, "UHF: unknown initial guess")
         return
      end select

      if (st%rohf) then
         ! An atomic guess arrives as two different Fock matrices; the
         ! symmetric ones are already one. Either way one diagonalisation.
         select case (guess_kind)
         case (SCF_GUESS_SAC, SCF_GUESS_SAD, SCF_GUESS_PROJ, SCF_GUESS_MINAO)
            call roothaan_fock(st%fock_a, st%fock_b, st%d_a, st%d_b, s)
         case default
            ! Core, GWH and SAP set both spins to one matrix above.
            continue
         end select
         call diagonalize(st%fock_a, x, n_ao, n_mo, st%coeff_a, st%eig_a, error)
         if (error%has_error()) return
         st%coeff_b = st%coeff_a
         st%eig_b = st%eig_a
      else
         call diagonalize(st%fock_a, x, n_ao, n_mo, st%coeff_a, st%eig_a, error)
         if (error%has_error()) return
         call diagonalize(st%fock_b, x, n_ao, n_mo, st%coeff_b, st%eig_b, error)
         if (error%has_error()) return
      end if
      call build_density_spin(st%coeff_a, n_alpha, st%d_a)
      call build_density_spin(st%coeff_b, n_beta, st%d_b)

      st%e_old = 0.0_dp
      result%converged = .false.

      call clk%lap(STAGE_SETUP)
      ! `present(xc)` is not the question: the bridge passes its context
      ! unconditionally, so it is present on a Hartree-Fock run too. `active`
      ! is what says a functional was constructed.
      kohn_sham_run = .false.
      if (present(xc)) kohn_sham_run = xc%active
      ! Once per context, on the total guess density; see `run_czt_rhf`.
      if (kohn_sham_run) then
         select case (guess_kind)
         case (SCF_GUESS_SAC, SCF_GUESS_SAD, SCF_GUESS_PROJ, SCF_GUESS_MINAO)
            call xc_prune_small_rho(xc, mol, guess_density_alpha + guess_density_beta, nelec, error)
         case default
            call xc_prune_small_rho(xc, mol, st%d_a + st%d_b, nelec, error)
         end select
         if (error%has_error()) return
      end if
      call scf_table_header(verbose, kohn_sham_run)

      ! pyscf's `conv_tol_grad`: derived from the energy tolerance unless a
      ! caller states it. See the argument's own note for why a density
      ! consumer must.
      ! Assembled once. A caller that states nothing gets `CONV_METRIC_STANDARD`
      ! with the bound `sqrt(energy_tol)`.
      if (present(convergence)) then
         conv = convergence
      else
         conv%tolerance = energy_tol
         if (present(scf)) conv%gradient_tolerance = scf%grad_tol
         if (present(grad_tol)) conv%gradient_tolerance = grad_tol
         ! The metric named by the group, when one was passed. A spelling the
         ! caller has already validated; an unparseable one leaves the default
         ! rather than failing a calculation from inside the SCF.
         if (present(scf)) then
            call parse_convergence_metric(scf%convergence_metric, metric_kind, metric_ok)
            if (metric_ok) conv%metric = metric_kind
         end if
      end if
      ! TODO(mqc): `gtol` is assigned here and never read -- `conv%is_converged`
      ! applies the commutator bound itself. Dead on both spin paths.
      gtol = conv%commutator_bound()

      ! The group first, then any individual override on top of it.
      shift = 0.0_dp
      if (present(scf)) shift = scf%level_shift
      if (present(level_shift)) shift = level_shift
      ! Refused rather than clamped: a negative shift lowers the virtuals into
      ! the occupied set, narrowing the gap the next density is built through
      ! and driving the oscillation it was asked to damp.
      if (shift < 0.0_dp) then
         call error%set(ERROR_VALIDATION, "keywords.scf.level_shift is negative. A "// &
                        "level shift raises the virtual orbitals to widen the gap; a "// &
                        "negative one narrows it and makes convergence worse, not "// &
                        "better. Give a positive value in Hartree, or leave it out.")
         return
      end if
      taper = 100.0_dp*density_tol
      st%drms_prev = huge(1.0_dp)
      shift_now = 0.0_dp

      ! Reported with its taper: a shift that is off by iteration three is not
      ! the shift the reader thinks they applied, and a shift dropped on the way
      ! here is otherwise indistinguishable from one that was never set.
      if (shift > 0.0_dp) then
         write (line, "(a,f8.4,a,es9.2)") "    level shift: ", shift, &
            " hartree, tapered off below dD ", taper
         call logger%info(trim(line))
      end if

      gdm_switch = 0
      soscf_switch = 0
      do iter = 1, max_iter
         call do_uhf_iteration(iter, mol, h, s, x, eri, bounds, n_ao, n_mo, n_alpha, n_beta, &
                               nsq, msq, want_incremental, use_pcm, shift, taper, start_cycle, &
                               accel, conv, verbose, kohn_sham_run, st, clk, xc, pcm, result, error)
         if (error%has_error()) return
         if (result%converged) exit
         ! The handover, on the commutator the iteration just measured, so the
         ! DIIS phase ends on a diagonalisation and the minimisation starts
         ! from orbitals and densities that agree.
         if (use_gdm .and. iter < max_iter .and. result%commutator < gdm_at) then
            gdm_switch = iter
            exit
         end if
         if (use_soscf .and. iter < max_iter .and. result%commutator < soscf_at) then
            soscf_switch = iter
            exit
         end if
      end do

      if (soscf_switch > 0) then
         call run_soscf_phase_uhf(mol, h, s, x, eri, bounds, n_mo, n_alpha, n_beta, conv, &
                                  verbose, st, clk, xc, soscf_switch, &
                                  max_iter - soscf_switch, result, error)
         if (error%has_error()) return
      end if

      if (gdm_switch > 0) then
         call run_gdm_phase_uhf(mol, h, s, x, eri, bounds, n_mo, n_alpha, n_beta, conv, &
                                verbose, gdm_cap_rad, st, clk, xc, gdm_switch, &
                                max_iter - gdm_switch, result, error)
         if (error%has_error()) return
      end if

      ! Rebuilt from the density that passed the test, as the restricted path
      ! does, so the reported energy belongs to the reported orbitals.
      call assemble_fock_uhf(mol, h, st%d_a, st%d_b, eri, bounds, xc, st%fock_a, st%fock_b, &
                             result%electronic, error, bmat=st%bmat, bmat_lr=st%bmat_lr)
      if (error%has_error()) return
      if (use_pcm) then
         call pcm%operator_matrix(mol, st%d_a + st%d_b, st%v_pcm, e_pcm, error)
         if (error%has_error()) return
         result%electronic = result%electronic + e_pcm
         result%pcm_energy = e_pcm
         result%pcm_charge = pcm%q_total
      end if
      if (st%rohf) then
         ! Each spin's orbitals are rotated to diagonalise that spin's own
         ! Fock matrix within its occupied and within its virtual block. No
         ! density moves. It is what lets everything downstream treat this as
         ! an unrestricted result: `C_occ eps C_occ^T` is then `D F D` for
         ! each spin, and their sum is the energy-weighted density of the
         ! restricted open-shell energy, because the blocks in which it
         ! differs from the Lagrangian's are the ones that vanish at
         ! convergence.
         if (use_pcm) then
            st%fock_a = st%fock_a + st%v_pcm
            st%fock_b = st%fock_b + st%v_pcm
         end if
         call spin_semicanonical(st%fock_a, n_alpha, st%coeff_a, st%eig_a)
         call spin_semicanonical(st%fock_b, n_beta, st%coeff_b, st%eig_b)
         if (error%has_error()) return
      end if
      result%nuclear_repulsion = mol%nuclear_repulsion()
      result%energy = result%electronic + result%nuclear_repulsion
      result%n_occupied = n_alpha
      result%n_occupied_beta = n_beta
      result%spin_squared = spin_contamination(st%coeff_a, st%coeff_b, s, n_alpha, n_beta)
      call clk%lap(STAGE_FOCK)
      call clk%finish()
      call scf_table_footer(verbose, result%converged, result%iterations)
      call clk%report("UHF", verbose)
      call xc_breakdown_summary(verbose)
      if (result%gdm_started_at > 0 .and. verbose) then
         write (line, "(a,i0,a,i0,a,i0,a)") "  direct minimisation took over after iteration ", &
            result%gdm_started_at, ", ", result%gdm_iterations, &
            " steps of one Fock build each, ", result%gdm_rejected, " rejected"
         call logger%info(trim(line))
      end if
      if (result%second_order_started_at > 0 .and. verbose) then
         write (line, "(a,i0,a,i0,a,i0,a,es10.2)") &
            "  second-order phase took over after iteration ", &
            result%second_order_started_at, ", ", result%second_order_iterations, &
            " Newton steps, ", result%second_order_hessian_products, &
            " Hessian products; lowest curvature seen ", result%second_order_curvature
         call logger%info(trim(line))
      end if
      if (verbose) then
         write (line, "(a,f12.8,a,f12.8)") "  <S^2> = ", result%spin_squared, &
            "   exact = ", 0.25_dp*real(n_alpha - n_beta, dp)*real(n_alpha - n_beta + 2, dp)
         call logger%info(trim(line))
      end if
      call move_alloc(st%eig_a, result%orbital_energies)
      call move_alloc(st%coeff_a, result%orbitals)
      call move_alloc(st%d_a, result%density)
      call move_alloc(st%eig_b, result%orbital_energies_beta)
      call move_alloc(st%coeff_b, result%orbitals_beta)
      call move_alloc(st%d_b, result%density_beta)
      call st%diis%destroy()

   contains

      function reference_label() result(label)
         !! What the log calls this wavefunction
         character(len=:), allocatable :: label

         label = "unrestricted"
         if (present(restricted_open)) then
            if (restricted_open) label = "restricted open shell"
         end if
      end function reference_label

      subroutine spin_semicanonical(fock_ao, n_occ, coeff, energies)
         !! One spin's Fock matrix diagonalised within its occupied and its
         !! virtual block; an empty spin is all virtual and left as it is
         real(dp), intent(in) :: fock_ao(:, :)
         integer, intent(in) :: n_occ
         real(dp), allocatable, intent(inout) :: coeff(:, :), energies(:)

         real(dp), allocatable :: fock_mo(:, :), eps(:)

         if (n_occ == 0 .or. n_occ == n_mo) return
         call soscf_semicanonicalize(fock_ao, coeff, n_occ, fock_mo, eps, error)
         if (error%has_error()) return
         energies = eps
      end subroutine spin_semicanonical
   end subroutine run_czt_uhf

   subroutine roothaan_fock(fock_a, fock_b, d_a, d_b, s)
      !! Replace both spins' Fock matrices by Roothaan's effective one
      !!
      !! With `Fc = (Fa + Fb)/2` and the projectors onto the doubly occupied,
      !! singly occupied and virtual spaces, `Pc = Db S`, `Po = (Da - Db) S`,
      !! `Pv = 1 - Da S`:
      !!
      !!              closed   open   virtual
      !!     closed     Fc      Fb      Fc
      !!     open       Fb      Fc      Fa
      !!     virtual    Fc      Fa      Fc
      !!
      !! Its off-diagonal blocks are the three conditions a restricted
      !! open-shell energy is stationary under; the diagonal blocks are a
      !! choice, this one Roothaan's, and move no density.
      real(dp), intent(inout) :: fock_a(:, :), fock_b(:, :)
         !! In: each spin's own. Out: the effective matrix, in both.
      real(dp), intent(in) :: d_a(:, :), d_b(:, :), s(:, :)

      real(dp), allocatable :: fc(:, :), pc(:, :), po(:, :), pv(:, :), eff(:, :)
      integer :: n, i

      n = size(s, 1)
      allocate (fc(n, n), pc(n, n), po(n, n), pv(n, n), eff(n, n))
      fc = 0.5_dp*(fock_a + fock_b)
      call pic_gemm(d_b, s, pc, beta=0.0_dp)
      call pic_gemm(d_a, s, pv, beta=0.0_dp)
      po = pv - pc
      pv = -pv
      do i = 1, n
         pv(i, i) = pv(i, i) + 1.0_dp
      end do
      ! Half of each diagonal block and the whole of each block below the
      ! diagonal; adding the transpose completes it.
      eff = 0.5_dp*(sandwich(pc, fc, pc) + sandwich(po, fc, po) + sandwich(pv, fc, pv)) &
            + sandwich(po, fock_b, pc) + sandwich(po, fock_a, pv) + sandwich(pv, fc, pc)
      fock_a = eff + transpose(eff)
      fock_b = fock_a

   contains

      function sandwich(left, mid, right) result(out)
         !! `left^T mid right`
         real(dp), intent(in) :: left(:, :), mid(:, :), right(:, :)
         real(dp) :: out(size(left, 2), size(right, 2))
         real(dp) :: work(size(mid, 1), size(right, 2))

         call pic_gemm(mid, right, work, beta=0.0_dp)
         call pic_gemm(left, work, out, transa="T", beta=0.0_dp)
      end function sandwich
   end subroutine roothaan_fock

   subroutine assemble_fock(mol, h, density, coeff, n_occ, bmat, eri, bounds, xc, &
                            fock, e_elec, error, clk, screening, incr, bmat_lr, jfit, cosx, &
                            energy_only)
      !! The Fock matrix for this density, and the electronic energy that belongs to it
      !!
      !! The energy comes back with the Fock matrix because for Kohn-Sham the
      !! two cannot be separated after the fact: `1/2 Tr[D (H + F)]` counts the
      !! potential at half weight, which is right for a mean field and wrong for
      !! a functional, whose energy is E_xc and not `1/2 Tr[D V_xc]`. So the
      !! energy is taken from the Fock matrix *before* V_xc is added, and E_xc
      !! is added on top.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: h(:, :), density(:, :)
      real(dp), intent(in) :: coeff(:, :)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(in) :: bmat(:, :)
      real(dp), allocatable, intent(in) :: eri(:, :, :, :)
      real(dp), allocatable, intent(in) :: bounds(:, :)
      type(df_jfit_t), intent(in), optional :: jfit
         !! The Coulomb-only fit, used when it is built (`jfit%raw` allocated)
         !! and `bmat` is not. Only a reference with no exact exchange builds one.
      real(dp), allocatable, intent(in), optional :: bmat_lr(:, :)
         !! The same fit against `erf(omega r)/r`, for a range-separated
         !! functional whose reference is density fitted. Unallocated otherwise,
         !! and only read when `xc%range_separated`.
      type(xc_context_t), intent(inout), optional :: xc
      real(dp), intent(out) :: fock(:, :)
      real(dp), intent(out) :: e_elec
      type(error_t), intent(inout) :: error
      type(timing_report_t), intent(inout), optional :: clk
         !! Present from the SCF loop, so the exchange-correlation quadrature is
         !! reported apart from the Coulomb/exchange build rather than inside it.
         !! Absent from the guess and gradient callers, which do not report.
      type(direct_stats_t), intent(out), optional :: screening
         !! How much of the quartet space each test removed. Counts rather than
         !! seconds, because they are exactly reproducible.
      type(incremental_state_t), intent(inout), optional :: incr
         !! Present from the SCF loop, which wants each build to cost only the
         !! change since the last one. Absent from the guess, the gradient
         !! callers and the final energy rebuild, all of which want an exact
         !! build.

      type(cosx_grid_t), intent(in), optional :: cosx
         !! The grid for a chain-of-spheres exchange. When it is present and
         !! built, the exchange comes from the quadrature and the Coulomb
         !! matrix from the fit: `jfit` when that was built, else `bmat`. `eri`
         !! and `bounds` are then not read.
      logical, intent(in), optional :: energy_only
         !! The caller wants `e_elec` and will not read `fock`. A
         !! chain-of-spheres build then takes the exchange energy as a trace
         !! and leaves the exchange out of `fock`, which is much cheaper.
         !! Default false, and ignored on every other path.
      type(direct_stats_t) :: stats
      real(dp), allocatable :: g_delta(:, :), d_delta(:, :), h_zero(:, :)
      real(dp), allocatable :: k_lr(:, :), k_cosx(:, :)
      logical :: use_cosx, trace_only
      real(dp) :: e_exchange, tr_dk
      logical :: use_incremental, full_build
      logical :: want_lr
      real(dp), allocatable :: v_xc(:, :)
      real(dp) :: k_scale, e_xc, n_elec
      logical :: have_lr
      logical :: kohn_sham
      logical :: use_jfit
      real(dp), allocatable :: j_fit(:, :)

      kohn_sham = .false.
      k_scale = 1.0_dp
      want_lr = .false.
      if (present(xc)) then
         if (xc%active) then
            kohn_sham = .true.
            ! Pure functionals want no Fock exchange, hybrids a fraction, and
            ! Hartree-Fock the fraction one -- hence a scale, not a branch.
            k_scale = xc%exx_fraction
            ! Settled once, because `xc` is optional and each of the three
            ! direct build sites would otherwise need its own `present` guard.
            want_lr = xc%range_separated
         end if
      end if

      ! A range-separated functional needs a second exchange matrix over the
      ! erf-attenuated kernel, and the two are combined as
      !
      !     K_eff = exx_fraction * K_full + rs_k_lr * K_lr(omega)
      !
      ! The direct build produces it as a second quartet pass with `omega` set,
      ! and the fitted build as a second fit (`bmat_lr`). The in-core tensor is
      ! refused: it is built for the full Coulomb kernel, and a second n_ao^4
      ! array is not worth having.
      if (present(xc)) then
         if (xc%range_separated .and. allocated(eri) .and. .not. allocated(bmat)) then
            call error%set(ERROR_VALIDATION, "a range-separated functional needs the "// &
                           "direct or the density-fitted Fock build: the in-core "// &
                           "integrals are built for the full Coulomb kernel, so the "// &
                           "long-range exchange would be missing. Run without in_core.")
            return
         end if
         ! `present` before `allocated`, and not merged into one expression:
         ! Fortran does not short-circuit, so `allocated(bmat_lr)` on an absent
         ! optional is a reference to something that is not there.
         if (xc%range_separated .and. allocated(bmat)) then
            have_lr = .false.
            if (present(bmat_lr)) have_lr = allocated(bmat_lr)
            if (.not. have_lr) then
               call error%set(ERROR_VALIDATION, "a range-separated functional fitted "// &
                              "with an auxiliary basis needs the attenuated tensor as "// &
                              "well, and it was not built. This is an internal "// &
                              "inconsistency rather than a deck error.")
               return
            end if
         end if
      end if

      ! Incremental building is for the direct path only: the fitted and in-core
      ! paths contract a stored tensor, whose cost does not depend on how large
      ! the density elements are.
      use_jfit = .false.
      if (present(jfit)) use_jfit = allocated(jfit%raw)
      ! Chain-of-spheres exchange goes with a fitted Coulomb matrix, from the
      ! Coulomb-only fit or, where that was not built, from the fitted tensor.
      use_cosx = .false.
      if (present(cosx)) use_cosx = cosx%active() .and. (use_jfit .or. allocated(bmat))
      e_exchange = 0.0_dp
      if (use_jfit .and. .not. use_cosx .and. (k_scale /= 0.0_dp .or. want_lr)) then
         call error%set(ERROR_VALIDATION, "a Coulomb-only density fit was built for a "// &
                        "functional that needs exact exchange, with nothing else to "// &
                        "supply it. This is an internal inconsistency rather than a "// &
                        "deck error.")
         return
      end if

      use_incremental = .false.
      if (present(incr) .and. .not. allocated(bmat) .and. .not. allocated(eri) &
          .and. .not. use_jfit) then
         use_incremental = .true.
      end if

      if (use_cosx .and. want_lr) then
         call error%set(ERROR_VALIDATION, "chain-of-spheres exchange was handed a "// &
                        "range-separated functional, which it does not build. This "// &
                        "is an internal inconsistency rather than a deck error.")
         return
      end if

      if (use_cosx) then
         ! J from the fit, K from the grid. Unfitted, the
         ! symmetrised quadrature has `Tr[D1 K(D2)] = Tr[D2 K(D1)]` and this
         ! Fock matrix is exactly the derivative of the energy taken from it
         ! below. The overlap fit breaks that by the size of its own
         ! correction, so a fitted SCF converges to a point that is stationary
         ! to within the quadrature error and no better.
         if (use_jfit) then
            allocate (j_fit(size(h, 1), size(h, 2)))
            call df_jfit_coulomb(jfit, density, j_fit)
            fock = h + j_fit
            deallocate (j_fit)
         else
            call build_fock_df(h, bmat, density, coeff, n_occ, fock, k_scale=0.0_dp)
         end if
         e_exchange = 0.0_dp
         trace_only = .false.
         if (present(energy_only)) trace_only = energy_only .and. cosx%shells%n_shell > 0
         if (k_scale /= 0.0_dp .and. trace_only) then
            call cosx_exchange_trace(mol, cosx%coords, cosx%weights, density, tr_dk, error, &
                                     cosx%shells, fit=cosx%fit, bounds=cosx%bounds, &
                                     screen_tol=cosx%screen_tol)
            if (error%has_error()) return
            ! E_x = -Tr[D K]/4, scaled by the exact-exchange fraction.
            e_exchange = -0.25_dp*k_scale*tr_dk
         else if (k_scale /= 0.0_dp) then
            ! An SCF's grids always carry the fit and the bounds.
            call cosx_exchange(mol, cosx%coords, cosx%weights, density, k_cosx, error, &
                               fit=cosx%fit, bounds=cosx%bounds, shells=cosx%shells, &
                               screen_tol=cosx%screen_tol)
            if (error%has_error()) return
            fock = fock - 0.5_dp*k_scale*k_cosx
            deallocate (k_cosx)
         end if
      else if (allocated(bmat)) then
         call build_fock_df(h, bmat, density, coeff, n_occ, fock, k_scale=k_scale)
         ! The long-range exchange, from the tensor fitted against the
         ! attenuated kernel. Same shape as the direct path's second pass: no
         ! core Hamiltonian and no Coulomb, so this returns K_lr alone.
         if (kohn_sham) then
            if (xc%range_separated) then
               block
                  real(dp), allocatable :: k_lr(:, :), h_zero(:, :)
                  allocate (k_lr(size(h, 1), size(h, 2)), h_zero(size(h, 1), size(h, 2)))
                  h_zero = 0.0_dp
                  call build_fock_df(h_zero, bmat_lr, density, coeff, n_occ, k_lr, &
                                     k_scale=xc%rs_k_lr, j_scale=0.0_dp)
                  fock = fock + k_lr
                  deallocate (k_lr, h_zero)
               end block
            end if
         end if
      else if (use_jfit) then
         allocate (j_fit(size(h, 1), size(h, 2)))
         call df_jfit_coulomb(jfit, density, j_fit)
         fock = h + j_fit
         deallocate (j_fit)
      else if (allocated(eri)) then
         call build_fock(h, eri, density, fock, k_scale=k_scale)
      else
         if (use_incremental) then
            full_build = .not. incr%active .or. incr%since_reset >= INCREMENTAL_RESET &
                         .or. incr%settled
            if (full_build) then
               incr%full_builds = incr%full_builds + 1
               call build_fock_direct(mol, h, density, bounds, fock, stats, error, &
                                      k_scale=k_scale)
               if (error%has_error()) return
               if (.not. allocated(incr%g_ref)) then
                  allocate (incr%g_ref(size(h, 1), size(h, 2)))
                  allocate (incr%d_ref(size(h, 1), size(h, 2)))
               end if
               ! G alone, so the next iteration's correction adds to a matrix that
               ! does not already contain H.
               incr%g_ref = fock - h
               incr%d_ref = density
               ! The attenuated pass, seeded from the same density in the same
               ! iteration so the two references never describe different
               ! densities. Added to `fock` here rather than by the block below,
               ! which now only runs when there is no incremental state.
               if (want_lr) then
                  if (.not. allocated(incr%g_ref_lr)) then
                     allocate (incr%g_ref_lr(size(h, 1), size(h, 2)))
                  end if
                  allocate (k_lr(size(h, 1), size(h, 2)), h_zero(size(h, 1), size(h, 2)))
                  h_zero = 0.0_dp
                  call build_fock_direct(mol, h_zero, density, bounds, k_lr, stats, &
                                         error, k_scale=xc%rs_k_lr, j_scale=0.0_dp, &
                                         omega=xc%rs_omega)
                  if (error%has_error()) return
                  incr%g_ref_lr = k_lr
                  fock = fock + k_lr
                  deallocate (k_lr, h_zero)
               end if
               incr%since_reset = 0
               incr%active = .true.
            else
               incr%updates = incr%updates + 1
               allocate (g_delta(size(h, 1), size(h, 2)), d_delta(size(h, 1), size(h, 2)))
               allocate (h_zero(size(h, 1), size(h, 2)))
               d_delta = density - incr%d_ref
               h_zero = 0.0_dp
               ! Zero for the core Hamiltonian, so this returns G(delta) and
               ! nothing has to be subtracted back out.
               !
               ! `density_screen=.true.` is the whole point of building
               ! incrementally: with it false the delta screens exactly as hard
               ! as a full build and saves nothing. The screen stays sound
               ! because the test is on the *contribution* against a fixed
               ! absolute tolerance, not on the delta relative to itself, and
               ! `INCREMENTAL_RESET` bounds how many dropped terms can
               ! accumulate before a full build re-syncs.
               !
               ! **CPHF still needs it false.** There the density is a Krylov
               ! trial vector the solver drives to zero, and the operator has to
               ! stay the same linear map from one matvec to the next.
               call build_fock_direct(mol, h_zero, d_delta, bounds, g_delta, stats, error, &
                                      k_scale=k_scale, density_screen=.true.)
               if (error%has_error()) return
               incr%g_ref = incr%g_ref + g_delta
               ! The attenuated exchange on the same delta. `erf(omega r)/r`
               ! is at most `1/r` pointwise, so the full-kernel Schwarz bounds
               ! screen it conservatively.
               if (want_lr) then
                  call build_fock_direct(mol, h_zero, d_delta, bounds, g_delta, stats, &
                                         error, k_scale=xc%rs_k_lr, j_scale=0.0_dp, &
                                         omega=xc%rs_omega, density_screen=.true.)
                  if (error%has_error()) return
                  incr%g_ref_lr = incr%g_ref_lr + g_delta
               end if
               incr%d_ref = density
               incr%since_reset = incr%since_reset + 1
               fock = h + incr%g_ref
               if (want_lr) fock = fock + incr%g_ref_lr
               deallocate (g_delta, d_delta, h_zero)
            end if
         else
            call build_fock_direct(mol, h, density, bounds, fock, stats, error, &
                                   k_scale=k_scale)
            if (error%has_error()) return
         end if

         ! The long-range exchange for the callers that build in full -- the
         ! guess, the gradient and the final energy rebuild. The incremental
         ! branch above accumulates its own into `g_ref_lr`.
         !
         ! Zero in place of the core Hamiltonian and no Coulomb term, so this
         ! pass returns the long-range exchange alone.
         if (want_lr .and. .not. use_incremental) then
            allocate (k_lr(size(h, 1), size(h, 2)), h_zero(size(h, 1), size(h, 2)))
            h_zero = 0.0_dp
            call build_fock_direct(mol, h_zero, density, bounds, k_lr, stats, &
                                   error, k_scale=xc%rs_k_lr, j_scale=0.0_dp, &
                                   omega=xc%rs_omega)
            if (error%has_error()) return
            fock = fock + k_lr
            deallocate (k_lr, h_zero)
         end if
      end if

      ! Taken here, from the Fock matrix that is still a mean field.
      e_elec = electronic_energy(h, fock, density)
      ! The exchange energy a chain-of-spheres build took as a trace, which
      ! `fock` then does not carry.
      if (use_cosx) e_elec = e_elec + e_exchange

      ! Close the Fock bucket before the quadrature opens. `lap` charges
      ! everything since the previous lap, which is in the caller, so a single
      ! `lap(STAGE_XC)` after the quadrature would bill the two-electron build
      ! to the quadrature as well and leave the caller's `lap(STAGE_FOCK)`
      ! nothing to charge.
      if (present(clk)) call clk%lap(STAGE_FOCK)

      if (kohn_sham) then
         allocate (v_xc(size(h, 1), size(h, 2)))
         call xc_add_potential(xc, mol, density, v_xc, e_xc, n_elec, error)
         if (error%has_error()) return
         fock = fock + v_xc
         e_elec = e_elec + e_xc
         if (present(clk)) call clk%lap(STAGE_XC)
         deallocate (v_xc)
      end if

      if (present(screening)) screening = stats
   end subroutine assemble_fock

   subroutine density_energy(mol, h, density, bounds, xc, energy, error)
      !! The closed-shell energy functional evaluated at `density`, no SCF
      !!
      !! `E[D; h]` in Hartree, `mol`'s nuclear repulsion included, for the
      !! functional `xc` names: the same [[assemble_fock]] the SCF's own final
      !! energy comes from, so the exact-exchange fraction, the range-separated
      !! exchange and `E_xc` are the SCF's. Direct Coulomb and exchange only.
      !! `h` is the one-electron matrix the energy is taken against and
      !! `bounds` are `mol`'s Schwarz bounds.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: h(:, :), density(:, :)
      real(dp), intent(in) :: bounds(:, :)
      type(xc_context_t), intent(inout) :: xc
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: fock(:, :), coeff(:, :), bnd(:, :)
      real(dp), allocatable :: no_bmat(:, :), no_eri(:, :, :, :)
      real(dp) :: e_elec

      energy = 0.0_dp
      allocate (fock(size(h, 1), size(h, 2)), coeff(0, 0))
      bnd = bounds
      call assemble_fock(mol, h, density, coeff, 0, no_bmat, no_eri, bnd, xc, fock, e_elec, error)
      if (error%has_error()) return
      energy = e_elec + mol%nuclear_repulsion()
   end subroutine density_energy

   subroutine assemble_fock_uhf(mol, h, d_alpha, d_beta, eri, bounds, xc, &
                                fock_a, fock_b, e_elec, error, clk, incr, bmat, bmat_lr)
      !! Both spin Fock matrices for this pair of densities, and their energy
      !!
      !! The unrestricted twin of `assemble_fock`. The energy comes back with the
      !! matrices because for Kohn-Sham it cannot be recovered afterwards:
      !! `uhf_electronic_energy` counts the potential at half weight, right for a
      !! mean field and wrong for a functional, so it is taken before V_xc is
      !! added and E_xc is added on top.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: h(:, :), d_alpha(:, :), d_beta(:, :)
      real(dp), allocatable, intent(in) :: eri(:, :, :, :)
      real(dp), allocatable, intent(in) :: bounds(:, :)
      type(xc_context_t), intent(inout), optional :: xc
      real(dp), intent(out) :: fock_a(:, :), fock_b(:, :)
      real(dp), intent(out) :: e_elec
      type(error_t), intent(inout) :: error
      type(timing_report_t), intent(inout), optional :: clk
         !! Charges the two-electron build and the quadrature to their own
         !! buckets. Optional because the guess and the final rebuild call this
         !! outside the iteration, where there is no per-iteration row to fill.
      type(incremental_state_t), intent(inout), optional :: incr
         !! Present from the SCF loop, which wants each build to cost only the
         !! change since the last one. Absent from the guess and the final energy
         !! rebuild, which want an exact build.
         !!
         !! Two spin channels and, for a range-separated functional, two kernels
         !! in each: four accumulated G in the worst case, all keyed to one
         !! `since_reset` so a rebuild re-syncs them together.
      real(dp), allocatable, intent(in), optional :: bmat(:, :)
         !! Fitted B(mu nu, P). Present and allocated builds J and K from the
         !! fit instead of the integrals; every such build is a full one.
      real(dp), allocatable, intent(in), optional :: bmat_lr(:, :)
         !! The same against the attenuated kernel, required beside `bmat`
         !! under a range-separated functional.

      type(direct_stats_t) :: stats
      logical :: fitted
      real(dp), allocatable :: v_a(:, :), v_b(:, :)
      real(dp), allocatable :: ga_delta(:, :), gb_delta(:, :)
      real(dp), allocatable :: da_delta(:, :), db_delta(:, :), h_zero(:, :)
      real(dp), allocatable :: k_lr_a(:, :), k_lr_b(:, :)
      real(dp) :: k_scale, e_xc, n_elec
      logical :: kohn_sham, want_lr, use_incremental, full_build
      integer :: n

      n = size(h, 1)
      kohn_sham = .false.
      k_scale = 1.0_dp
      want_lr = .false.
      if (present(xc)) then
         if (xc%active) then
            kohn_sham = .true.
            k_scale = xc%exx_fraction
            want_lr = xc%range_separated
         end if
      end if

      ! As in the closed-shell path: only the direct build can produce the
      ! attenuated exchange matrix, so in-core is refused rather than dropped.
      if (present(xc)) then
         if (xc%range_separated .and. allocated(eri)) then
            call error%set(ERROR_VALIDATION, "a range-separated functional needs the "// &
                           "direct Fock build: the in-core integrals are built for the "// &
                           "full Coulomb kernel, and the long-range exchange would be "// &
                           "missing. Run without in_core.")
            return
         end if
      end if

      fitted = .false.
      if (present(bmat)) fitted = allocated(bmat)
      if (fitted .and. want_lr) then
         ! Nested, because Fortran does not short-circuit and `bmat_lr` may be
         ! absent.
         full_build = .false.
         if (present(bmat_lr)) full_build = allocated(bmat_lr)
         if (.not. full_build) then
            call error%set(ERROR_VALIDATION, "a range-separated functional over a "// &
                           "fitted unrestricted build needs the tensor fitted against "// &
                           "the attenuated kernel, and none was passed")
            return
         end if
      end if

      ! Incremental building is for the direct path only, as on the restricted
      ! side.
      use_incremental = .false.
      if (present(incr) .and. .not. allocated(eri) .and. .not. fitted) use_incremental = .true.

      if (fitted) then
         if (want_lr) then
            call build_fock_df_uhf(h, bmat, d_alpha, d_beta, fock_a, fock_b, error, &
                                   k_scale=k_scale, b_lr=bmat_lr, k_lr_scale=xc%rs_k_lr)
         else
            call build_fock_df_uhf(h, bmat, d_alpha, d_beta, fock_a, fock_b, error, &
                                   k_scale=k_scale)
         end if
         if (error%has_error()) return
      else if (allocated(eri)) then
         call build_fock_uhf(h, eri, d_alpha, d_beta, fock_a, fock_b, k_scale=k_scale)
      else if (use_incremental) then
         full_build = .not. incr%active .or. incr%since_reset >= INCREMENTAL_RESET &
                      .or. incr%settled
         if (full_build) then
            call build_fock_direct_uhf(mol, h, d_alpha, d_beta, bounds, fock_a, fock_b, &
                                       stats, error, k_scale=k_scale)
            if (error%has_error()) return
            if (.not. allocated(incr%g_ref)) then
               allocate (incr%g_ref(n, n), incr%d_ref(n, n))
               allocate (incr%g_ref_b(n, n), incr%d_ref_b(n, n))
            end if
            ! G alone in each spin, so the next iteration's correction adds to a
            ! matrix that does not already contain H.
            incr%g_ref = fock_a - h
            incr%g_ref_b = fock_b - h
            incr%d_ref = d_alpha
            incr%d_ref_b = d_beta
            if (want_lr) then
               if (.not. allocated(incr%g_ref_lr)) then
                  allocate (incr%g_ref_lr(n, n), incr%g_ref_lr_b(n, n))
               end if
               ! Zero core Hamiltonian and no Coulomb, so this pass returns the
               ! long-range exchange of each spin and nothing to subtract back.
               allocate (k_lr_a(n, n), k_lr_b(n, n), h_zero(n, n))
               h_zero = 0.0_dp
               call build_fock_direct_uhf(mol, h_zero, d_alpha, d_beta, bounds, &
                                          k_lr_a, k_lr_b, stats, error, &
                                          k_scale=xc%rs_k_lr, j_scale=0.0_dp, &
                                          omega=xc%rs_omega)
               if (error%has_error()) return
               incr%g_ref_lr = k_lr_a
               incr%g_ref_lr_b = k_lr_b
               fock_a = fock_a + k_lr_a
               fock_b = fock_b + k_lr_b
               deallocate (k_lr_a, k_lr_b, h_zero)
            end if
            incr%since_reset = 0
            incr%active = .true.
         else
            allocate (ga_delta(n, n), gb_delta(n, n))
            allocate (da_delta(n, n), db_delta(n, n), h_zero(n, n))
            da_delta = d_alpha - incr%d_ref
            db_delta = d_beta - incr%d_ref_b
            h_zero = 0.0_dp
            ! Zero for the core Hamiltonian, so this returns G(delta) and
            ! nothing has to be subtracted back out.
            !
            ! `build_fock_direct_uhf` always weights its Schwarz bound by the
            ! density it multiplies, so the delta screens hard here with no flag
            ! to ask for it. The restricted path needs `density_screen=.true.`
            ! explicitly only because it also serves CPHF, which must opt out.
            call build_fock_direct_uhf(mol, h_zero, da_delta, db_delta, bounds, &
                                       ga_delta, gb_delta, stats, error, k_scale=k_scale)
            if (error%has_error()) return
            incr%g_ref = incr%g_ref + ga_delta
            incr%g_ref_b = incr%g_ref_b + gb_delta
            if (want_lr) then
               call build_fock_direct_uhf(mol, h_zero, da_delta, db_delta, bounds, &
                                          ga_delta, gb_delta, stats, error, &
                                          k_scale=xc%rs_k_lr, j_scale=0.0_dp, &
                                          omega=xc%rs_omega)
               if (error%has_error()) return
               incr%g_ref_lr = incr%g_ref_lr + ga_delta
               incr%g_ref_lr_b = incr%g_ref_lr_b + gb_delta
            end if
            incr%d_ref = d_alpha
            incr%d_ref_b = d_beta
            incr%since_reset = incr%since_reset + 1
            fock_a = h + incr%g_ref
            fock_b = h + incr%g_ref_b
            if (want_lr) then
               fock_a = fock_a + incr%g_ref_lr
               fock_b = fock_b + incr%g_ref_lr_b
            end if
            deallocate (ga_delta, gb_delta, da_delta, db_delta, h_zero)
         end if
      else
         call build_fock_direct_uhf(mol, h, d_alpha, d_beta, bounds, fock_a, fock_b, &
                                    stats, error, k_scale=k_scale)
         if (error%has_error()) return
         if (want_lr) then
            ! Zero core Hamiltonian and no Coulomb, so this pass returns the
            ! long-range exchange of each spin and nothing to subtract back.
            allocate (k_lr_a(n, n), k_lr_b(n, n), h_zero(n, n))
            h_zero = 0.0_dp
            call build_fock_direct_uhf(mol, h_zero, d_alpha, d_beta, bounds, &
                                       k_lr_a, k_lr_b, stats, error, &
                                       k_scale=xc%rs_k_lr, j_scale=0.0_dp, &
                                       omega=xc%rs_omega)
            if (error%has_error()) return
            fock_a = fock_a + k_lr_a
            fock_b = fock_b + k_lr_b
            deallocate (k_lr_a, k_lr_b, h_zero)
         end if
      end if

      ! Taken here, while the Fock matrices are still a mean field.
      e_elec = uhf_electronic_energy(h, fock_a, fock_b, d_alpha, d_beta)

      ! Close the Fock bucket before the quadrature opens, for the reason
      ! `assemble_fock` sets out.
      if (present(clk)) call clk%lap(STAGE_FOCK)

      if (kohn_sham) then
         allocate (v_a(n, n), v_b(n, n))
         call xc_add_potential_uks(xc, mol, d_alpha, d_beta, v_a, v_b, e_xc, n_elec, error)
         if (error%has_error()) return
         fock_a = fock_a + v_a
         fock_b = fock_b + v_b
         e_elec = e_elec + e_xc
         if (present(clk)) call clk%lap(STAGE_XC)
         deallocate (v_a, v_b)
      end if
   end subroutine assemble_fock_uhf

   subroutine commutator(fock, density, overlap, x, err)
      !! e = X^T (F D S - S D F) X, the DIIS error vector
      !!
      !! Projected into the orthogonal basis rather than left in the AO one,
      !! where its size would depend on the overlap's conditioning.
      real(dp), intent(in) :: fock(:, :), density(:, :), overlap(:, :), x(:, :)
      real(dp), intent(out) :: err(:, :)

      real(dp), allocatable :: fd(:, :), fds(:, :), sd(:, :), sdf(:, :), work(:, :)
      integer :: n_ao, n_mo

      n_ao = size(fock, 1)
      n_mo = size(x, 2)
      allocate (fd(n_ao, n_ao), fds(n_ao, n_ao), sd(n_ao, n_ao), sdf(n_ao, n_ao))
      allocate (work(n_ao, n_mo))

      call pic_gemm(fock, density, fd)
      call pic_gemm(fd, overlap, fds)
      call pic_gemm(overlap, density, sd)
      call pic_gemm(sd, fock, sdf)
      fds = fds - sdf

      call pic_gemm(fds, x, work)
      call pic_gemm(x, work, err, transa="T")
   end subroutine commutator

   pure subroutine guess_fock(overlap, h, fock)
      !! Generalized Wolfsberg-Helmholz starting Fock
      !!
      !!   F_ij = 1/2 K S_ij (H_ii + H_jj),   F_ii = H_ii
      !!
      !! A crude stand-in for the screening the core guess leaves out, costing
      !! nothing beyond two matrices already in hand. It needs no free-atom SCF,
      !! so it is the fallback when an atomic guess fails, and it is what cuEST
      !! starts from.
      real(dp), intent(in) :: overlap(:, :), h(:, :)
      real(dp), intent(out) :: fock(:, :)

      integer :: i, j, n

      n = size(h, 1)
      do j = 1, n
         do i = 1, n
            if (i == j) then
               fock(i, j) = h(i, j)
            else
               fock(i, j) = 0.5_dp*GWH_K*overlap(i, j)*(h(i, i) + h(j, j))
            end if
         end do
      end do
   end subroutine guess_fock

   subroutine sap_fock(mol, h, fock, error)
      !! Starting Fock from a superposition of atomic potentials
      !!
      !!   F = H + V_SAP
      !!
      !! `V_SAP` is the screening one electron feels from the electrons of the
      !! free atoms, fitted to Gaussians and contracted against three-centre
      !! integrals -- no quadrature grid and no free-atom SCF, which is what
      !! separates it from SAD. `H` already carries the nuclear attraction.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: h(:, :)
      real(dp), intent(out) :: fock(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: v(:, :)

      call build_sap_potential(mol, SAP_BASIS_DEFAULT, v, error)
      if (error%has_error()) return
      fock = h + v
   end subroutine sap_fock

   subroutine atomic_guess_fock(mol, h, density, bmat, eri, bounds, fock, error, jfit, xc, cosx)
      !! One Fock build from a guess density, through whichever path is active
      !!
      !! Brings an atomic guess to the loop in the same state as the core guess:
      !! a Fock matrix waiting to be diagonalised.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: h(:, :), density(:, :)
      real(dp), allocatable, intent(in) :: bmat(:, :)
      real(dp), allocatable, intent(in) :: eri(:, :, :, :)
      real(dp), allocatable, intent(in) :: bounds(:, :)
      real(dp), intent(out) :: fock(:, :)
      type(error_t), intent(inout) :: error
      type(df_jfit_t), intent(in), optional :: jfit
         !! Present with a Coulomb-only fit. The guess Fock is then `H + J`
         !! plus whatever stands in for the exchange the fit cannot give: the
         !! chain-of-spheres exchange when `cosx` is built, otherwise the
         !! exchange-correlation potential of the guess density. `H + J` alone
         !! is a poor starting point that costs iterations.
      type(xc_context_t), intent(inout), optional :: xc
      type(cosx_grid_t), intent(in), optional :: cosx
         !! The exchange grid of a chain-of-spheres run

      real(dp), allocatable :: pseudo(:, :), j_fit(:, :), v_xc(:, :), k_cosx(:, :)
      real(dp) :: e_xc, n_elec
      integer :: n_modes
      type(direct_stats_t) :: stats
      logical :: use_jfit, seminumerical

      use_jfit = .false.
      if (present(jfit)) use_jfit = allocated(jfit%raw)
      seminumerical = .false.
      if (present(cosx)) seminumerical = cosx%active()

      if (use_jfit) then
         allocate (j_fit(size(h, 1), size(h, 2)))
         call df_jfit_coulomb(jfit, density, j_fit)
         fock = h + j_fit
         if (seminumerical) then
            call cosx_exchange(mol, cosx%coords, cosx%weights, density, k_cosx, error, &
                               fit=cosx%fit, bounds=cosx%bounds, shells=cosx%shells, &
                               screen_tol=cosx%screen_tol)
            if (error%has_error()) return
            fock = fock - 0.5_dp*k_cosx
         else if (present(xc)) then
            allocate (v_xc(size(h, 1), size(h, 2)))
            call xc_add_potential(xc, mol, density, v_xc, e_xc, n_elec, error)
            if (error%has_error()) return
            fock = fock + v_xc
         end if
      else if (allocated(bmat)) then
         call density_pseudo_orbitals(density, pseudo, n_modes, error)
         if (error%has_error()) return
         call build_fock_df(h, bmat, density, pseudo, n_modes, fock)
      else if (allocated(eri)) then
         call build_fock(h, eri, density, fock)
      else
         call build_fock_direct(mol, h, density, bounds, fock, stats, error)
      end if
   end subroutine atomic_guess_fock

   subroutine density_pseudo_orbitals(density, coeff, n_modes, error)
      !! Columns c_i satisfying D = 2 sum_i c_i c_i^T, for a density with no orbitals
      !!
      !! `build_fock_df` gets its exchange from occupied orbitals rather than
      !! from the density. A guess density is not idempotent and so has no
      !! orbitals in the usual sense; its eigenvectors serve instead, since with
      !! D = sum_i w_i v_i v_i^T and c_i = v_i sqrt(w_i/2) the identity
      !! D = 2 sum_i c_i c_i^T holds exactly.
      !!
      !! Costs one n^3 diagonalisation, on the guess only.
      real(dp), intent(in) :: density(:, :)
      real(dp), allocatable, intent(out) :: coeff(:, :)
      integer, intent(out) :: n_modes
      type(error_t), intent(inout) :: error

      real(dp), parameter :: OCCUPATION_FLOOR = 1.0e-12_dp
         !! Occupations at or below this are dropped. A density built from
         !! converged atomic solutions is positive semi-definite in exact
         !! arithmetic, so anything here is rounding on an empty mode and
         !! contributes nothing to K; keeping one would mean taking the square
         !! root of a negative number.
      real(dp), allocatable :: vectors(:, :), values(:)
      integer :: n, i, info

      n = size(density, 1)
      allocate (vectors(n, n), values(n))
      vectors = density
      call pic_syev(vectors, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "guess: the guess density would not diagonalise")
         return
      end if

      n_modes = count(values > OCCUPATION_FLOOR)
      if (n_modes == 0) then
         call error%set(ERROR_VALIDATION, "guess: the guess density carries no occupation")
         return
      end if

      allocate (coeff(n, n_modes))
      n_modes = 0
      do i = 1, n
         if (values(i) <= OCCUPATION_FLOOR) cycle
         n_modes = n_modes + 1
         coeff(:, n_modes) = vectors(:, i)*sqrt(0.5_dp*values(i))
      end do
   end subroutine density_pseudo_orbitals

   subroutine diagonalize(fock, x, n_ao, n_mo, coeff, eigenvalues, error)
      !! F' = X^T F X, diagonalise, and bring the orbitals back: C = X C'
      real(dp), intent(in) :: fock(:, :), x(:, :)
      integer, intent(in) :: n_ao, n_mo
      real(dp), allocatable, intent(inout) :: coeff(:, :), eigenvalues(:)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: work(:, :), f_ortho(:, :)
      integer :: info

      allocate (work(n_ao, n_mo), f_ortho(n_mo, n_mo))
      call pic_gemm(fock, x, work)                       ! F X
      call pic_gemm(x, work, f_ortho, transa="T")        ! X^T F X

      if (allocated(eigenvalues)) deallocate (eigenvalues)
      allocate (eigenvalues(n_mo))
      call pic_syev(f_ortho, eigenvalues, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "RHF: Fock diagonalisation failed")
         return
      end if

      if (allocated(coeff)) deallocate (coeff)
      allocate (coeff(n_ao, n_mo))
      call pic_gemm(x, f_ortho, coeff)                   ! C = X C'
   end subroutine diagonalize

   subroutine build_fock_df(h, b, density, coeff, n_occ, fock, k_scale, j_scale)
      !! F = H + J - K/2 from the fitted tensor rather than the exact ERIs
      !!
      !! Neither term ever forms a four-index object. What is stored is B, which
      !! is n^2 by n_aux -- n^3 rather than n^4.
      !!
      !! J is two contractions with the density: c_P = sum_uv B(uv,P) D_uv, then
      !! J_uv = sum_P B(uv,P) c_P. Both are BLAS-2 against the flattened tensor,
      !! blocked over P.
      !!
      !! K goes through the occupied orbitals rather than the density. Writing
      !! D = 2 sum_i C_ui C_vi and substituting,
      !!
      !!    K_uv = sum_P sum_ls B(ul,P) D_ls B(sv,P)
      !!         = 2 sum_P sum_i W(u,i,P) W(v,i,P),   W^P = B^P C_occ
      !!
      !! which costs 2 n^2 n_occ per auxiliary function instead of the 2 n^3 of
      !! contracting the full density -- a factor of n/n_occ.
      !!
      !! **Threaded over P, not collapsed into one GEMM.** Stacking every W^P
      !! side by side would turn sum_P W^P (W^P)^T into a single large product,
      !! which is the right move against a threaded BLAS; this project links a
      !! sequential one on purpose, so one GEMM is one core however large it is.
      !! Each auxiliary function owns a private partial K, merged once at the
      !! end.
      real(dp), intent(in) :: h(:, :), density(:, :)
      real(dp), intent(in) :: b(:, :)
      contiguous :: b
         !! So that a column block `b(:, p0:p1)` reaches BLAS as a view, and
         !! `b(:, p)` reaches `df_exchange_slice` as one. Without it the compiler
         !! is entitled to copy a column per auxiliary function.
         !!
         !! On its own line because `fortitude` mis-scopes the argument list when
         !! `contiguous` shares a declaration with `intent`.
      real(dp), intent(in) :: coeff(:, :)   !! MO coefficients; only the occupied block is read
      integer, intent(in) :: n_occ
      real(dp), intent(out) :: fock(:, :)
      real(dp), intent(in), optional :: k_scale   !! Exact-exchange fraction, default one
      real(dp), intent(in), optional :: j_scale
         !! Coulomb fraction, default one. Zero is the long-range exchange pass
         !! of a range-separated functional, where the attenuated tensor must not
         !! contribute a second Coulomb term: J is already complete from the
         !! full-range pass.

      real(dp) :: kf, jf
      real(dp), allocatable :: c(:), j_flat(:), j_local(:), k(:, :), w(:, :), c_occ(:, :)
      real(dp), allocatable :: d_flat(:), k_local(:, :)
      integer :: n, naux, p, p0, p1

      n = size(h, 1)
      naux = size(b, 2)

      kf = 0.5_dp
      if (present(k_scale)) kf = 0.5_dp*k_scale
      jf = 1.0_dp
      if (present(j_scale)) jf = j_scale

      fock = h

      ! Coulomb. Skipped outright when it is scaled away, which is the
      ! long-range exchange pass of a range-separated functional.
      if (jf /= 0.0_dp) then
         allocate (c(naux), j_flat(n*n), d_flat(n*n))
         d_flat = reshape(density, [n*n])

         !$omp parallel do default(none) &
         !$omp    shared(b, d_flat, c, naux) private(p0, p1) schedule(static)
         do p0 = 1, naux, DF_AUX_CHUNK
            p1 = min(p0 + DF_AUX_CHUNK - 1, naux)
            call pic_gemv(b(:, p0:p1), d_flat, c(p0:p1), trans_a="T")
         end do
         !$omp end parallel do

         j_flat = 0.0_dp
         !$omp parallel default(none) &
         !$omp    shared(b, c, j_flat, naux, n) private(p0, p1, j_local)
         allocate (j_local(n*n))
         j_local = 0.0_dp
         !$omp do schedule(static)
         do p0 = 1, naux, DF_AUX_CHUNK
            p1 = min(p0 + DF_AUX_CHUNK - 1, naux)
            call pic_gemv(b(:, p0:p1), c(p0:p1), j_local, beta=1.0_dp)
         end do
         !$omp end do
         !$omp critical(mqc_build_fock_df_coulomb)
         j_flat = j_flat + j_local
         !$omp end critical(mqc_build_fock_df_coulomb)
         deallocate (j_local)
         !$omp end parallel

         fock = fock + jf*reshape(j_flat, [n, n])
         deallocate (c, j_flat, d_flat)
      end if

      ! Exchange. Skipped when there is none to build -- a pure functional
      ! carries no exact exchange, and this is the whole cost of it.
      if (kf /= 0.0_dp .and. n_occ > 0) then
         allocate (c_occ(n, n_occ), k(n, n))
         c_occ = coeff(:, 1:n_occ)
         k = 0.0_dp

         !$omp parallel default(none) &
         !$omp    shared(b, c_occ, k, n, n_occ, naux) private(p, w, k_local)
         allocate (w(n, n_occ), k_local(n, n))
         k_local = 0.0_dp
         !$omp do schedule(static)
         do p = 1, naux
            call df_exchange_slice(b(:, p), c_occ, w, k_local, n)
         end do
         !$omp end do
         !$omp critical(mqc_build_fock_df_exchange)
         k = k + k_local
         !$omp end critical(mqc_build_fock_df_exchange)
         deallocate (w, k_local)
         !$omp end parallel

         ! `pic_syrk` filled the lower triangle only. Mirroring it here costs
         ! one pass over n^2 where doing it inside the loop would cost n_aux.
         do p = 1, n
            do p0 = 1, p - 1
               k(p0, p) = k(p, p0)
            end do
         end do

         fock = fock - kf*k
         deallocate (c_occ, k)
      end if
   end subroutine build_fock_df

   subroutine df_exchange_slice(b_p, c_occ, w, k_local, n)
      !! One auxiliary function's contribution to K
      !!
      !! `b_p` is explicit-shape against a contiguous rank-1 actual argument, so
      !! a column of B is seen as the (n, n) matrix it already is by sequence
      !! association -- free, where a `reshape` would materialise the block.
      !! `w` is scratch that never leaves: W^P is formed and consumed here.
      integer, intent(in) :: n
      real(dp), intent(in) :: b_p(n, n)
      real(dp), intent(in) :: c_occ(:, :)
      real(dp), intent(inout) :: w(:, :)
      real(dp), intent(inout) :: k_local(:, :)

      call pic_gemm(b_p, c_occ, w, beta=0.0_dp)
      ! A rank-k update, not a general product: W W^T is symmetric. Only the
      ! lower triangle is written, and the caller mirrors it once at the end
      ! rather than n_aux times here.
      call pic_syrk(w, k_local, uplo="L", alpha=2.0_dp, beta=1.0_dp)
   end subroutine df_exchange_slice

   subroutine build_fock_df_uhf(h, b, d_alpha, d_beta, fock_a, fock_b, error, &
                                k_scale, b_lr, k_lr_scale)
      !! F_sigma = H + J(D_alpha + D_beta) - k K(D_sigma) from the fitted tensor
      !!
      !! Three calls to `build_fock_df`: the Coulomb matrix once, on the total
      !! density, then each spin's exchange alone. That routine takes its
      !! exchange from columns `c` with `D = 2 sum_i c_i c_i^T` and subtracts
      !! `k K(D)/2`, so a spin density is handed over with twice the fraction.
      !!
      !! The columns are the eigenvectors of each spin density rather than its
      !! orbitals, so this is right for any positive semi-definite density --
      !! an atomic guess, or one an accelerator has mixed -- at the cost of one
      !! diagonalisation per spin per build.
      real(dp), intent(in) :: h(:, :), d_alpha(:, :), d_beta(:, :)
      real(dp), intent(in) :: b(:, :)   !! (nao*nao, naux)
      contiguous :: b
      real(dp), intent(out) :: fock_a(:, :), fock_b(:, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: k_scale   !! Exact-exchange fraction, default one
      real(dp), intent(in), optional :: b_lr(:, :)
         !! The tensor fitted against `erf(omega r)/r`; present adds
         !! `- k_lr_scale K_lr(D_sigma)`
      contiguous :: b_lr
      real(dp), intent(in), optional :: k_lr_scale

      real(dp), allocatable :: zero(:, :), k_spin(:, :), pseudo(:, :), no_coeff(:, :)
      real(dp) :: kf
      integer :: n, n_modes

      n = size(h, 1)
      kf = 1.0_dp
      if (present(k_scale)) kf = k_scale

      ! Coulomb only; no orbitals are read when the exchange is scaled away.
      allocate (no_coeff(n, 0))
      call build_fock_df(h, b, d_alpha + d_beta, no_coeff, 0, fock_a, k_scale=0.0_dp)
      fock_b = fock_a
      if (kf == 0.0_dp .and. .not. present(b_lr)) return

      allocate (zero(n, n), k_spin(n, n))
      zero = 0.0_dp
      call spin_exchange(d_alpha, fock_a)
      if (error%has_error()) return
      call spin_exchange(d_beta, fock_b)

   contains

      subroutine spin_exchange(d_spin, fock)
         !! Subtract one spin's fitted exchange from its Fock matrix
         real(dp), intent(in) :: d_spin(:, :)
         real(dp), intent(inout) :: fock(:, :)

         ! An empty spin channel has no exchange, and no density to
         ! diagonalise.
         if (maxval(abs(d_spin)) == 0.0_dp) return
         call density_pseudo_orbitals(d_spin, pseudo, n_modes, error)
         if (error%has_error()) return
         if (kf /= 0.0_dp) then
            ! Returns `-kf K(D_sigma)`: no core Hamiltonian, no Coulomb.
            call build_fock_df(zero, b, d_spin, pseudo, n_modes, k_spin, &
                               k_scale=2.0_dp*kf, j_scale=0.0_dp)
            fock = fock + k_spin
         end if
         if (present(b_lr)) then
            call build_fock_df(zero, b_lr, d_spin, pseudo, n_modes, k_spin, &
                               k_scale=2.0_dp*k_lr_scale, j_scale=0.0_dp)
            fock = fock + k_spin
         end if
         deallocate (pseudo)
      end subroutine spin_exchange
   end subroutine build_fock_df_uhf

   pure subroutine build_fock_uhf(h, eri, d_alpha, d_beta, fock_a, fock_b, k_scale)
      !! F_sigma = H + J(D_alpha + D_beta) - K(D_sigma), straight from the ERIs
      !!
      !! The reference the direct build is checked against, and slow on purpose.
      real(dp), intent(in) :: h(:, :), eri(:, :, :, :), d_alpha(:, :), d_beta(:, :)
      real(dp), intent(out) :: fock_a(:, :), fock_b(:, :)
      real(dp), intent(in), optional :: k_scale
         !! Fraction of exact exchange. One for Hartree-Fock and the default,
         !! less for a hybrid Kohn-Sham build.

      integer :: mu, nu, la, si, n
      real(dp) :: dt, kf

      kf = 1.0_dp
      if (present(k_scale)) kf = k_scale

      n = size(h, 1)
      fock_a = h
      fock_b = h
      do nu = 1, n
         do mu = 1, n
            do si = 1, n
               do la = 1, n
                  dt = d_alpha(la, si) + d_beta(la, si)
                  ! Coulomb from both spins, exchange from the same spin only.
                  fock_a(mu, nu) = fock_a(mu, nu) + dt*eri(mu, nu, la, si) &
                                   - kf*d_alpha(la, si)*eri(mu, la, nu, si)
                  fock_b(mu, nu) = fock_b(mu, nu) + dt*eri(mu, nu, la, si) &
                                   - kf*d_beta(la, si)*eri(mu, la, nu, si)
               end do
            end do
         end do
      end do
   end subroutine build_fock_uhf

   pure function uhf_electronic_energy(h, fock_a, fock_b, d_alpha, d_beta) result(energy)
      !! E = 1/2 [ sum (D_a + D_b) H + sum D_a F_a + sum D_b F_b ]
      real(dp), intent(in) :: h(:, :), fock_a(:, :), fock_b(:, :), d_alpha(:, :), d_beta(:, :)
      real(dp) :: energy

      energy = 0.5_dp*(sum((d_alpha + d_beta)*h) + sum(d_alpha*fock_a) + sum(d_beta*fock_b))
   end function uhf_electronic_energy

   pure function electronic_energy(h, fock, density) result(energy)
      !! E = 1/2 sum_uv D_uv (H_uv + F_uv)
      real(dp), intent(in) :: h(:, :), fock(:, :), density(:, :)
      real(dp) :: energy

      energy = 0.5_dp*sum(density*(h + fock))
   end function electronic_energy

end module mqc_czt_rhf
