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
   use mqc_czt_direct, only: schwarz_bounds, build_fock_direct, build_fock_direct_uhf, &
                             direct_stats_t
   use mqc_czt_integrals, only: czt_molecule_t, build_df_tensor
   use mqc_czt_xc, only: xc_context_t, xc_add_potential, xc_add_potential_uks
   use mqc_czt_pcm, only: pcm_context_t
   use mqc_program_limits, only: DF_AUX_CHUNK
   use mqc_scf_convergence, only: scf_convergence_t, parse_convergence_metric
   use mqc_scf_types, only: scf_numerics_t
   use mqc_diis, only: parse_accelerator_name
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
   integer, parameter, public :: SCF_GUESS_SAC = 2    !! Superposed atomic coefficients
   integer, parameter, public :: SCF_GUESS_SAD = 3    !! Superposed atomic densities
   integer, parameter, public :: SCF_GUESS_PROJ = 4   !! Projected from a smaller basis

   integer, parameter :: LINE_LEN = 160
      !! Buffer length for a formatted table line handed to the logger.

   public :: rhf_result_t
   public :: run_czt_rhf
   public :: run_czt_uhf
   public :: guess_fock                !! Starting Fock for the guesses needing no atomic SCF
   public :: build_fock                !! F = H + J - K/2; with H zero it is the response operator
   public :: density_pseudo_orbitals   !! Factor a guess density for the fitted exchange

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
         !! every `INCREMENTAL_RESET` iterations; with it off, one per
         !! iteration.
         !!
         !! Counted because it is the only observable difference the setting
         !! makes: an incremental build is exact to the convergence threshold,
         !! so the energy, the iteration count and every property come out
         !! identical either way.
      integer :: updates = 0                 !! Builds from the density change
      logical :: active = .false.            !! Off until the first full build seeds it
   end type incremental_state_t

   type :: rhf_result_t
      !! What a converged closed-shell SCF leaves behind
      real(dp) :: energy = 0.0_dp              !! Total, including nuclear repulsion
      real(dp) :: electronic = 0.0_dp          !! Without it
      real(dp) :: nuclear_repulsion = 0.0_dp
      integer :: iterations = 0
      integer :: full_fock_builds = 0
         !! Fock matrices built in full, as opposed to updated from the density
         !! change. Equals `iterations` when `incremental_fock` is off, and one
         !! plus a rebuild every `INCREMENTAL_RESET` iterations when it is on.
      integer :: incremental_updates = 0
         !! Fock matrices built from the density change. Zero when
         !! `incremental_fock` is off.
      logical :: converged = .false.
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
      real(dp), allocatable :: eri(:, :, :, :)  !! In-core integrals; unallocated on the direct path
      real(dp), allocatable :: bounds(:, :)   !! Schwarz bounds for screening
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
   end type rhf_state_t

   integer, parameter :: INCREMENTAL_RESET = 16
      !! Iterations between full rebuilds of the accumulated G.
      !!
      !! An incremental build adds a correction per iteration, so its rounding
      !! accumulates where a full build's does not; rebuilding every sixteenth
      !! iteration keeps the drift at the level of the convergence threshold
      !! for one full build in sixteen.

   ! Stage labels, named once so the per-iteration column and the summary table
   ! cannot drift apart, and so a caller can ask `clk%seconds_of(STAGE_FOCK)`.
   character(len=*), parameter :: STAGE_SETUP = "setup (1e, bounds, guess)"
   character(len=*), parameter :: STAGE_FOCK = "Fock builds"
   character(len=*), parameter :: STAGE_DIAG = "diagonalisation"
   character(len=*), parameter :: STAGE_DIIS = "DIIS"
   character(len=*), parameter :: STAGE_XC = "XC quadrature"
   character(len=*), parameter :: STAGE_PCM = "continuum operator"

contains

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
      call logger%performance("")
      call logger%performance("  screening (last Fock build)")
      write (line, "(a,i14)") "    unique quartets      ", stats%quartets_total
      call logger%performance(trim(line))
      write (line, "(a,i14,f9.1,a)") "    skipped, Schwarz     ", stats%screened_schwarz, &
         100.0_dp*real(stats%screened_schwarz, dp)/real(stats%quartets_total, dp), " %"
      call logger%performance(trim(line))
      write (line, "(a,i14,f9.1,a)") "    skipped, density     ", stats%screened_density, &
         100.0_dp*real(stats%screened_density, dp)/real(stats%quartets_total, dp), " %"
      call logger%performance(trim(line))
      write (line, "(a,i14,f9.1,a)") "    computed             ", stats%quartets_computed, &
         100.0_dp*real(stats%quartets_computed, dp)/real(stats%quartets_total, dp), " %"
      call logger%performance(trim(line))
      ! One is perfect balance. A ratio inside one run, so a contended node
      ! largely divides out of it.
      write (line, "(a,f14.4)") "    thread imbalance     ", stats%thread_imbalance
      call logger%performance(trim(line))
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

      call logger%info("")
      call logger%info("  energy components")
      write (line, "(a,f22.10)") "    one electron                ", one_e
      call logger%info(trim(line))
      write (line, "(a,f22.10)") "    two electron                ", two_e
      call logger%info(trim(line))
      write (line, "(a,f22.10)") "    nuclear repulsion           ", nuclear
      call logger%info(trim(line))
      write (line, "(a,f22.10)") "    total                       ", electronic + nuclear
      call logger%info(trim(line))
      write (line, "(a,f22.10)") "    electron-electron potential ", v_ee
      call logger%info(trim(line))
      write (line, "(a,f22.10)") "    nucleus-electron potential  ", v_ne
      call logger%info(trim(line))
      write (line, "(a,f22.10)") "    nucleus-nucleus potential   ", nuclear
      call logger%info(trim(line))
      write (line, "(a,f22.10)") "    total potential             ", v_total
      call logger%info(trim(line))
      write (line, "(a,f22.10)") "    total kinetic               ", t_energy
      call logger%info(trim(line))
      if (abs(t_energy) > 0.0_dp) then
         write (line, "(a,f22.10)") "    virial ratio (-V/T)         ", -v_total/t_energy
         call logger%info(trim(line))
      end if
      deallocate (kinetic, core)
   end subroutine energy_components

   subroutine run_czt_rhf(mol, nelec, max_iter, energy_tol, density_tol, &
                          verbose, result, error, aux, diis_vectors, in_core, &
                          guess, guess_density, xc, h_extra, pcm, projector, &
                          level_shift, linear_dependence, b_ao_out, accelerator, grad_tol, &
                          incremental_fock, convergence, scf)
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
         !! Total starting density, required by SCF_GUESS_SAC and SCF_GUESS_SAD
         !! and ignored otherwise. Built by `mqc_czt_atomic_guess`, which
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
      logical :: use_in_core
      real(dp), allocatable :: bounds(:, :)
      ! TODO(mqc): `stats` is declared here and never used; the screening counts
      ! come back in `screening`, which `assemble_fock` fills.
      type(direct_stats_t) :: stats

      type(rhf_operators_t) :: ops
      type(rhf_controls_t) :: ctrl
      type(rhf_state_t) :: st
      real(dp) :: e_pcm
      real(dp) :: gtol
      logical :: accel_ok_grp
      integer :: metric_kind
      logical :: metric_ok
      integer :: iter
      type(timing_report_t) :: clk
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
      call mol%overlap(ops%s)
      call mol%core_hamiltonian(ops%h)
      if (present(h_extra)) then
         if (size(h_extra, 1) /= ops%n_ao .or. size(h_extra, 2) /= ops%n_ao) then
            call error%set(ERROR_VALIDATION, "RHF: h_extra is not ops%n_ao square")
            return
         end if
         ops%h = ops%h + h_extra
      end if
      if (present(aux)) then
         call build_df_tensor(mol, aux, ops%bmat, error)
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
      if (present(scf)) call parse_accelerator_name(scf%accelerator, ctrl%accel, accel_ok_grp)
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
      case (SCF_GUESS_SAC, SCF_GUESS_SAD, SCF_GUESS_PROJ)
         if (.not. present(guess_density)) then
            call error%set(ERROR_VALIDATION, "RHF: an atomic guess was asked for but no "// &
                           "guess st%density was supplied")
            return
         end if
         if (size(guess_density, 1) /= ops%n_ao .or. size(guess_density, 2) /= ops%n_ao) then
            call error%set(ERROR_VALIDATION, "RHF: the guess st%density is not the size of "// &
                           "this basis")
            return
         end if
         st%density = guess_density
         call atomic_guess_fock(mol, ops%h, st%density, ops%bmat, ops%eri, ops%bounds, st%fock, error)
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
                        "level ctrl%shift raises the virtual orbitals to widen the gap; a "// &
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
         write (line, "(a,f8.4,a,es9.2)") "    level ctrl%shift: ", ctrl%shift, &
            " hartree, tapered off below dD ", ctrl%taper
         call logger%info(trim(line))
      end if

      ! `verbose` is a dummy of this routine, so it is the one control that is
      ! copied across rather than having been written straight into `ctrl`.
      ctrl%verbose = verbose

      ! One call per iteration, and the three views above are what makes that
      ! readable: the alternative was a forty-one argument call.
      !
      ! `exit` cannot live inside the callee, so convergence comes back on
      ! `result%converged` -- which is where it was recorded anyway -- and the
      ! error check is the caller's, exactly as it is after every other call
      ! in this routine.
      do iter = 1, max_iter
         call do_rhf_iteration(iter, mol, ops, ctrl, st, clk, xc, pcm, projector, &
                               result, error)
         if (error%has_error()) return
         if (result%converged) exit
      end do

      ! The energy that goes out belongs to the density that satisfied the
      ! test, so it is recomputed from the final Fock. Deliberately without
      ! `incr`: this one is a full build, and the number that leaves must not
      ! carry the corrections accumulated since the last reset.
      call assemble_fock(mol, ops%h, st%density, st%coeff, ops%n_occ, ops%bmat, ops%eri, ops%bounds, xc, &
                         st%fock, result%electronic, error, clk=clk, bmat_lr=ops%bmat_lr)
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
      call energy_components(verbose, mol, result%density, result%electronic, &
                             result%nuclear_repulsion)
      call screening_summary(verbose, st%screening)
      call clk%report("RHF", verbose)
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
      if (ctrl%incremental) then
         call assemble_fock(mol, ops%h, st%density, st%coeff, ops%n_occ, ops%bmat, ops%eri, ops%bounds, xc, &
                            st%fock, st%e_elec, error, clk=clk, screening=st%screening, incr=st%incr, &
                            bmat_lr=ops%bmat_lr)
      else
         ! Counted here rather than in the assembler: without `incr` it has
         ! no state to record into, and every build here is a full one.
         st%incr%full_builds = st%incr%full_builds + 1
         call assemble_fock(mol, ops%h, st%density, st%coeff, ops%n_occ, ops%bmat, ops%eri, ops%bounds, xc, &
                            st%fock, st%e_elec, error, clk=clk, screening=st%screening, &
                            bmat_lr=ops%bmat_lr)
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
      result%full_fock_builds = st%incr%full_builds
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

   subroutine run_czt_uhf(mol, nelec, multiplicity, max_iter, energy_tol, density_tol, &
                          verbose, result, error, diis_vectors, in_core, diis_start, &
                          guess, guess_density_alpha, guess_density_beta, xc, pcm, &
                          level_shift, linear_dependence, accelerator, grad_tol, &
                          incremental_fock, convergence, scf)
      !! Drive an unrestricted SCF to convergence
      !!
      !! Two Fock matrices sharing one Coulomb field: F_sigma = H + J(D_a + D_b)
      !! - K(D_sigma). They are diagonalised separately and extrapolated
      !! together, one DIIS over the pair.
      !!
      !! **No density fitting on this path.** The fitted J and K are written for
      !! one density; asking for both is refused rather than silently run
      !! restricted.
      type(czt_molecule_t), intent(in) :: mol
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
         !! Starting spin densities, required by SCF_GUESS_SAC and SCF_GUESS_SAD.
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
      real(dp), allocatable :: x(:, :), fock_a(:, :), fock_b(:, :)
      real(dp), allocatable :: d_a(:, :), d_b(:, :), d_a_old(:, :), d_b_old(:, :)
      real(dp), allocatable :: coeff_a(:, :), coeff_b(:, :), eig_a(:), eig_b(:)
      real(dp), allocatable :: err_a(:, :), err_b(:, :), fock_flat(:), err_flat(:)
      real(dp), allocatable :: dens_flat(:)
      real(dp), allocatable :: v_pcm(:, :)
      logical :: use_pcm
      real(dp) :: e_pcm
      type(diis_state_t) :: diis
      logical :: extrapolated
      logical :: kohn_sham_run
      integer :: accel
      real(dp) :: e_elec, e_old, de, drms
      real(dp) :: gnorm
      real(dp) :: gtol
      type(scf_convergence_t) :: conv
      logical :: accel_ok_grp
      integer :: metric_kind
      logical :: metric_ok
      real(dp) :: shift, shift_now, drms_prev, taper
      real(dp), allocatable :: sd(:, :), sds(:, :)
      integer :: n_ao, n_mo, n_alpha, n_beta, iter, nsq, msq
      type(timing_report_t) :: clk
      logical :: want_incremental
      type(incremental_state_t) :: incr
         !! Owned by the loop and handed only to the loop's own build. The guess
         !! and the final rebuild go without it, so the energy that leaves this
         !! routine is an exact build.
      real(dp) :: t_fock_iter, t_rest_iter, t_xc_iter
      real(dp) :: s_min, s_kept   !! overlap conditioning, reported before iteration 1
      real(dp) :: lindep
      character(len=LINE_LEN) :: line

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
         call logger%info("  unrestricted: charge = "// &
                          to_char(nint(sum(mol%charges)) - nelec)// &
                          "  multiplicity = "//to_char(multiplicity)// &
                          "  electrons = "//to_char(nelec))
         call logger%info("  unrestricted: n_alpha = "//to_char(n_alpha)// &
                          "  n_beta = "//to_char(n_beta))
      end if

      call clk%start()
      call mol%overlap(s)
      call mol%core_hamiltonian(h)
      if (use_in_core) then
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
      allocate (fock_a(n_ao, n_ao), fock_b(n_ao, n_ao))
      allocate (d_a(n_ao, n_ao), d_b(n_ao, n_ao), d_a_old(n_ao, n_ao), d_b_old(n_ao, n_ao))
      allocate (err_a(n_mo, n_mo), err_b(n_mo, n_mo))
      allocate (fock_flat(2*nsq), err_flat(2*msq), dens_flat(2*nsq))
      if (use_pcm) allocate (v_pcm(n_ao, n_ao))

      ! One subspace over both spins, so an extrapolation step moves them
      ! together. The vectors are the two Fock matrices laid end to end and the
      ! two commutators likewise.
      ! Resolved before the subspace is built, not after it. `energy_based`
      ! decides whether the history can hold densities and energies at all, and
      ! reading `accel` here while it was still assigned some fifty lines below
      ! took that decision on an undefined value.
      accel = ACCEL_DIIS
      if (present(scf)) call parse_accelerator_name(scf%accelerator, accel, accel_ok_grp)
      if (present(accelerator)) accel = accelerator
      call diis%init(diis_size, 2*nsq, 2*msq, &
                     energy_based=(accel /= ACCEL_DIIS))

      ! The symmetric guesses -- core and GWH -- give alpha and beta the same
      ! orbitals, and the occupations separate them: n_alpha > n_beta puts an
      ! electron where beta has none. For n_alpha == n_beta they converge to the
      ! restricted solution, which is the right answer rather than a failure.
      ! SAC is the one that arrives already asymmetric.
      select case (guess_kind)
      case (SCF_GUESS_CORE)
         fock_a = h
         fock_b = h
      case (SCF_GUESS_GWH)
         call guess_fock(s, h, fock_a)
         fock_b = fock_a
      case (SCF_GUESS_SAC, SCF_GUESS_SAD, SCF_GUESS_PROJ)
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
         d_a = guess_density_alpha
         d_b = guess_density_beta
         call assemble_fock_uhf(mol, h, d_a, d_b, eri, bounds, xc, fock_a, fock_b, &
                                e_elec, error)
         if (error%has_error()) return
      case default
         call error%set(ERROR_VALIDATION, "UHF: unknown initial guess")
         return
      end select

      call diagonalize(fock_a, x, n_ao, n_mo, coeff_a, eig_a, error)
      if (error%has_error()) return
      call diagonalize(fock_b, x, n_ao, n_mo, coeff_b, eig_b, error)
      if (error%has_error()) return
      call build_density_spin(coeff_a, n_alpha, d_a)
      call build_density_spin(coeff_b, n_beta, d_b)

      e_old = 0.0_dp
      result%converged = .false.

      call clk%lap(STAGE_SETUP)
      ! `present(xc)` is not the question: the bridge passes its context
      ! unconditionally, so it is present on a Hartree-Fock run too. `active`
      ! is what says a functional was constructed.
      kohn_sham_run = .false.
      if (present(xc)) kohn_sham_run = xc%active
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
      drms_prev = huge(1.0_dp)
      shift_now = 0.0_dp

      ! Reported with its taper: a shift that is off by iteration three is not
      ! the shift the reader thinks they applied, and a shift dropped on the way
      ! here is otherwise indistinguishable from one that was never set.
      if (shift > 0.0_dp) then
         write (line, "(a,f8.4,a,es9.2)") "    level shift: ", shift, &
            " hartree, tapered off below dD ", taper
         call logger%info(trim(line))
      end if

      do iter = 1, max_iter
         d_a_old = d_a
         d_b_old = d_b

         ! Read before the build, not after it: both stages are cumulative, so
         ! a reading taken once the build has charged them differences against
         ! itself.
         t_fock_iter = clk%seconds_of(STAGE_FOCK)
         t_xc_iter = clk%seconds_of(STAGE_XC)
         ! As on the restricted path: `incr` present is the switch, so the deck
         ! turns incremental building off by withholding it.
         if (want_incremental) then
            call assemble_fock_uhf(mol, h, d_a, d_b, eri, bounds, xc, fock_a, fock_b, &
                                   e_elec, error, clk=clk, incr=incr)
         else
            ! No `incr` to record into, and every build here is a full one.
            incr%full_builds = incr%full_builds + 1
            call assemble_fock_uhf(mol, h, d_a, d_b, eri, bounds, xc, fock_a, fock_b, &
                                   e_elec, error, clk=clk)
         end if
         if (error%has_error()) return
         ! The continuum sees one density -- the total -- and both spins feel
         ! one potential, as they would from any classical charge.
         if (use_pcm) then
            call pcm%operator_matrix(mol, d_a + d_b, v_pcm, e_pcm, error)
            if (error%has_error()) return
            fock_a = fock_a + v_pcm
            fock_b = fock_b + v_pcm
            e_elec = e_elec + e_pcm
         end if
         ! TODO(mqc): `assemble_fock_uhf` already laps STAGE_FOCK when it is
         ! given a clock, so this second lap counts two Fock builds per
         ! iteration and charges the continuum solve to the Fock bucket, where
         ! the restricted path gives it STAGE_PCM.
         call clk%lap(STAGE_FOCK)
         t_fock_iter = clk%seconds_of(STAGE_FOCK) - t_fock_iter
         t_xc_iter = clk%seconds_of(STAGE_XC) - t_xc_iter

         call commutator(fock_a, d_a, s, x, err_a)
         call commutator(fock_b, d_b, s, x, err_b)
         fock_flat(1:nsq) = reshape(fock_a, [nsq])
         fock_flat(nsq + 1:2*nsq) = reshape(fock_b, [nsq])
         err_flat(1:msq) = reshape(err_a, [msq])
         err_flat(msq + 1:2*msq) = reshape(err_b, [msq])
         gnorm = maxval(abs(err_flat))
         dens_flat(1:nsq) = reshape(d_a, [nsq])
         dens_flat(nsq + 1:2*nsq) = reshape(d_b, [nsq])
         call diis%push(fock_flat, err_flat, density=dens_flat, energy=e_elec)
         extrapolated = .false.
         if (iter >= start_cycle) then
            call diis%extrapolate_with(scheme_now(accel, gnorm), fock_flat, extrapolated)
         end if
         if (extrapolated) then
            fock_a = reshape(fock_flat(1:nsq), [n_ao, n_ao])
            fock_b = reshape(fock_flat(nsq + 1:2*nsq), [n_ao, n_ao])
         end if
         t_rest_iter = clk%seconds_of(STAGE_DIIS)
         call clk%lap(STAGE_DIIS)
         t_rest_iter = clk%seconds_of(STAGE_DIIS) - t_rest_iter

         ! After the extrapolation and after the two commutators, for the reason
         ! `run_czt_rhf` sets out. Each spin gets its own virtual projector:
         ! `build_density_spin` gives an occupation of one, so the closed-shell
         ! factor of a half is absent and the projector is `S - S D_sigma S`.
         shift_now = 0.0_dp
         if (shift > 0.0_dp .and. drms_prev > taper) shift_now = shift
         if (shift_now > 0.0_dp) then
            if (.not. allocated(sd)) allocate (sd(n_ao, n_ao), sds(n_ao, n_ao))
            call pic_gemm(s, d_a, sd, beta=0.0_dp)
            call pic_gemm(sd, s, sds, beta=0.0_dp)
            fock_a = fock_a + shift_now*(s - sds)
            call pic_gemm(s, d_b, sd, beta=0.0_dp)
            call pic_gemm(sd, s, sds, beta=0.0_dp)
            fock_b = fock_b + shift_now*(s - sds)
         end if

         call diagonalize(fock_a, x, n_ao, n_mo, coeff_a, eig_a, error)
         if (error%has_error()) return
         call diagonalize(fock_b, x, n_ao, n_mo, coeff_b, eig_b, error)
         if (error%has_error()) return
         call build_density_spin(coeff_a, n_alpha, d_a)
         call build_density_spin(coeff_b, n_beta, d_b)
         t_rest_iter = t_rest_iter - clk%seconds_of(STAGE_DIAG)
         call clk%lap(STAGE_DIAG)
         t_rest_iter = t_rest_iter + clk%seconds_of(STAGE_DIAG)

         de = abs(e_elec - e_old)
         drms = sqrt((sum((d_a - d_a_old)**2) + sum((d_b - d_b_old)**2))/real(2*nsq, dp))
         call scf_table_row(verbose, iter, e_elec + mol%nuclear_repulsion(), de, gnorm, &
                            diis%count(), t_fock_iter, t_xc_iter, t_rest_iter, kohn_sham_run)

         e_old = e_elec
         result%iterations = iter
         result%full_fock_builds = incr%full_builds
         result%incremental_updates = incr%updates
         drms_prev = drms
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
            exit
         end if
      end do

      ! Rebuilt from the density that passed the test, as the restricted path
      ! does, so the reported energy belongs to the reported orbitals.
      call assemble_fock_uhf(mol, h, d_a, d_b, eri, bounds, xc, fock_a, fock_b, &
                             result%electronic, error)
      if (error%has_error()) return
      if (use_pcm) then
         call pcm%operator_matrix(mol, d_a + d_b, v_pcm, e_pcm, error)
         if (error%has_error()) return
         result%electronic = result%electronic + e_pcm
         result%pcm_energy = e_pcm
         result%pcm_charge = pcm%q_total
      end if
      result%nuclear_repulsion = mol%nuclear_repulsion()
      result%energy = result%electronic + result%nuclear_repulsion
      result%n_occupied = n_alpha
      result%n_occupied_beta = n_beta
      result%spin_squared = spin_contamination(coeff_a, coeff_b, s, n_alpha, n_beta)
      call clk%lap(STAGE_FOCK)
      call clk%finish()
      call scf_table_footer(verbose, result%converged, result%iterations)
      call clk%report("UHF", verbose)
      if (verbose) then
         write (line, "(a,f12.8,a,f12.8)") "  <S^2> = ", result%spin_squared, &
            "   exact = ", 0.25_dp*real(n_alpha - n_beta, dp)*real(n_alpha - n_beta + 2, dp)
         call logger%info(trim(line))
      end if
      call move_alloc(eig_a, result%orbital_energies)
      call move_alloc(coeff_a, result%orbitals)
      call move_alloc(d_a, result%density)
      call move_alloc(eig_b, result%orbital_energies_beta)
      call move_alloc(coeff_b, result%orbitals_beta)
      call move_alloc(d_b, result%density_beta)
      call diis%destroy()
   end subroutine run_czt_uhf

   subroutine assemble_fock(mol, h, density, coeff, n_occ, bmat, eri, bounds, xc, &
                            fock, e_elec, error, clk, screening, incr, bmat_lr)
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

      type(direct_stats_t) :: stats
      real(dp), allocatable :: g_delta(:, :), d_delta(:, :), h_zero(:, :)
      real(dp), allocatable :: k_lr(:, :)
      logical :: use_incremental, full_build
      logical :: want_lr
      real(dp), allocatable :: v_xc(:, :)
      real(dp) :: k_scale, e_xc, n_elec
      logical :: have_lr
      logical :: kohn_sham

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
      use_incremental = .false.
      if (present(incr) .and. .not. allocated(bmat) .and. .not. allocated(eri)) then
         use_incremental = .true.
      end if

      if (allocated(bmat)) then
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
      else if (allocated(eri)) then
         call build_fock(h, eri, density, fock, k_scale=k_scale)
      else
         if (use_incremental) then
            full_build = .not. incr%active .or. incr%since_reset >= INCREMENTAL_RESET
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

   subroutine assemble_fock_uhf(mol, h, d_alpha, d_beta, eri, bounds, xc, &
                                fock_a, fock_b, e_elec, error, clk, incr)
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

      type(direct_stats_t) :: stats
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

      ! Incremental building is for the direct path only, as on the restricted
      ! side.
      use_incremental = .false.
      if (present(incr) .and. .not. allocated(eri)) use_incremental = .true.

      if (allocated(eri)) then
         call build_fock_uhf(h, eri, d_alpha, d_beta, fock_a, fock_b, k_scale=k_scale)
      else if (use_incremental) then
         full_build = .not. incr%active .or. incr%since_reset >= INCREMENTAL_RESET
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

   subroutine atomic_guess_fock(mol, h, density, bmat, eri, bounds, fock, error)
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

      real(dp), allocatable :: pseudo(:, :)
      integer :: n_modes
      type(direct_stats_t) :: stats

      if (allocated(bmat)) then
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

   subroutine build_fock(h, eri, density, fock, k_scale)
      !! `F = H + J - K/2` from a stored two-electron tensor
      !!
      !! **Correct for an antisymmetric density as well as a symmetric one.**
      !! Nothing here assumes either, and the Coulomb term vanishes of its own
      !! accord for an antisymmetric density because the integral is symmetric
      !! in the contracted pair.
      !!
      !! **Blocked over the ket pair, not over the target element.** For a fixed
      !! `(c, d)` the slice `eri(:, :, c, d)` is contiguous and carries both
      !! terms: Coulomb is `J += D(c,d) * eri(:,:,c,d)`, an axpy over the block,
      !! and exchange is `K(:,c) += eri(:,:,c,d) . D(:,d)`, a matrix-vector
      !! product on the same block. The tensor is walked once, in order, and
      !! each block used twice while it is still in cache.
      real(dp), intent(in) :: h(:, :), eri(:, :, :, :), density(:, :)
      real(dp), intent(out) :: fock(:, :)
      real(dp), intent(in), optional :: k_scale   !! Exact-exchange fraction, default one

      integer :: a, b, c, d, n
      real(dp) :: kf, dcd
      real(dp), allocatable :: j_mat(:, :), k_mat(:, :), j_local(:, :)

      n = size(h, 1)
      kf = 0.5_dp
      if (present(k_scale)) kf = 0.5_dp*k_scale

      allocate (j_mat(n, n), k_mat(n, n))
      j_mat = 0.0_dp
      k_mat = 0.0_dp

      ! Threaded over `c`, which makes the exchange update safe without a
      ! reduction: `K(:, c)` belongs to exactly one thread. Coulomb accumulates
      ! over every `(c, d)` and does need one, as a single n*n array per thread
      ! merged once.
      !
      ! The `if` keeps genuinely tiny systems serial: the response solver enters
      ! this region tens of thousands of times, and for a handful of basis
      ! functions the barriers cost more than the arithmetic.
      !$omp parallel default(none) private(a, b, c, d, dcd, j_local) &
      !$omp shared(n, eri, density, j_mat, k_mat) if(n >= 16)
      allocate (j_local(n, n))
      j_local = 0.0_dp
      !$omp do schedule(static)
      do c = 1, n
         do d = 1, n
            dcd = density(c, d)
            do b = 1, n
               do a = 1, n
                  j_local(a, b) = j_local(a, b) + dcd*eri(a, b, c, d)
               end do
            end do
            ! beta = 1 to accumulate over d rather than overwrite.
            call pic_gemv(eri(:, :, c, d), density(:, d), k_mat(:, c), beta=1.0_dp)
         end do
      end do
      !$omp end do
      !$omp critical(mqc_build_fock_coulomb)
      j_mat = j_mat + j_local
      !$omp end critical(mqc_build_fock_coulomb)
      deallocate (j_local)
      !$omp end parallel

      fock = h + j_mat - kf*k_mat
      deallocate (j_mat, k_mat)
   end subroutine build_fock

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
