!! Closed-shell SCF on the CPU
module mqc_libcint_rhf
   !! Restricted Hartree-Fock over libcint integrals, in core.
   !!
   !! The same algorithm the cuEST SCF runs -- canonical orthogonalisation,
   !! DIIS on the FDS - SDF commutator, energy from the density -- written
   !! against host arrays instead of device pointers. It is not a shared
   !! implementation, and that is deliberate: the two agreeing is worth more
   !! than the duplication costs, because agreement between two independent
   !! codes is the only check the GPU path has ever had.
   !!
   !! `mqc_diis` is shared, though. That one is already backend-neutral and
   !! there is nothing to learn from writing it twice.
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
   use mqc_diis, only: diis_state_t
   use mqc_fock_projector, only: fock_projector_t
   use mqc_libcint_direct, only: schwarz_bounds, build_fock_direct, build_fock_direct_uhf, &
                                 direct_stats_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_df_tensor
   use mqc_libcint_xc, only: xc_context_t, xc_add_potential, xc_add_potential_uks
   use mqc_libcint_pcm, only: pcm_context_t
   use mqc_program_limits, only: DF_AUX_CHUNK
   implicit none
   private

   !> First iteration the unrestricted SCF is allowed to extrapolate on.
   !>
   !> Extrapolating from the start converges an open-shell system to the wrong
   !> state. The core guess gives alpha and beta the same spatial orbitals and
   !> leaves degenerate shells degenerate, and DIIS combines Fock matrices
   !> linearly -- so a history that begins symmetric stays symmetric, and the
   !> iteration never reaches a solution that has to break it. OH in def2-SVP
   !> converges tidily to -75.167538 that way, a sigma hole with its pi pair
   !> intact to six digits, where the pi-hole ground state is -75.325108.
   !>
   !> Measured rather than chosen: on that case start=1 and 2 both give the
   !> wrong state and every value from 3 to 14 gives the right one, in 14 to 16
   !> iterations against the 14 the wrong answer took. Four is three there and
   !> one spare, and the spare is free -- the whole range costs the same.
   !>
   !> This is why the restricted path has no such delay. It has no symmetry to
   !> break, so the undamped iterations would buy it nothing.
   integer, parameter :: DEFAULT_UHF_DIIS_START = 4

   !> Which starting point the SCF is handed.
   !>
   !> Named here rather than in the atomic-guess module because that module
   !> calls this one -- the free-atom solutions an atomic guess is built from
   !> are themselves unrestricted SCFs -- and Fortran has no circular `use`.
   !> cuEST splits them the same way and for the same reason.
   integer, parameter, public :: SCF_GUESS_CORE = 0   !! F = H
   integer, parameter, public :: SCF_GUESS_GWH = 1    !! Generalized Wolfsberg-Helmholz
   integer, parameter, public :: SCF_GUESS_SAC = 2    !! Superposed atomic coefficients
   integer, parameter, public :: SCF_GUESS_SAD = 3    !! Superposed atomic densities
   integer, parameter, public :: SCF_GUESS_PROJ = 4   !! Projected from a smaller basis

   !> Buffer length for a formatted table line handed to the logger.
   integer, parameter :: LINE_LEN = 160

   public :: rhf_result_t
   public :: run_libcint_rhf
   public :: run_libcint_uhf
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
      !! and the correction can be built from the density *change* rather than
      !! the density. That is the whole point: `build_fock_direct` screens on the
      !! largest density element a quartet will multiply, and once the SCF is
      !! settling the change is orders of magnitude smaller than the density. The
      !! same screening test then discards far more, and discards more of it with
      !! every iteration -- which is the compounding the plain density weighting
      !! could not give on its own.
      !!
      !! Only the Coulomb and exchange part. The exchange-correlation potential is
      !! not linear in the density and is rebuilt in full every iteration.
      !!
      !! **Kohn-Sham works and is checked, but gains little.** `g_ref` is captured
      !! before the long-range exchange and before V_xc are added, so neither is
      !! double counted, and `check_dft` agrees with PySCF to 1e-11 for a pure
      !! functional, a hybrid and two range-separated ones through this path.
      !! What it does not do is help much: for water in cc-pVDZ the quadrature is
      !! 0.44 s of a 0.48 s SCF and the Fock build is too small to measure, so
      !! this accelerates 8% of the run. The balance moves with size -- the
      !! quadrature grows about linearly in the basis where the screened Fock
      !! build tends toward quadratic -- but where they cross is not measured.
      !! For a Kohn-Sham calculation the quadrature is the thing to look at.
      !!
      !! Two gaps, both deliberate rather than overlooked. A range-separated
      !! functional needs a second exchange matrix over the attenuated kernel, and
      !! that one is still rebuilt from the full density every iteration -- it
      !! could carry its own `g_ref`, and until it does such a functional saves on
      !! one of its two passes. And `assemble_fock_uhf` takes no state, so every
      !! open-shell calculation builds in full; the mechanism carries over, with
      !! two spin densities and two accumulated G, but it is a separate change.
      real(dp), allocatable :: d_ref(:, :)   !! Density `g_ref` belongs to
      real(dp), allocatable :: g_ref(:, :)   !! Accumulated G, without H
      integer :: since_reset = 0             !! Iterations since the last full build
      logical :: active = .false.            !! Off until the first full build seeds it
   end type incremental_state_t

   type :: rhf_result_t
      !! What a converged closed-shell SCF leaves behind
      real(dp) :: energy = 0.0_dp              !! Total, including nuclear repulsion
      real(dp) :: electronic = 0.0_dp          !! Without it
      real(dp) :: nuclear_repulsion = 0.0_dp
      integer :: iterations = 0
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

   !> Iterations between full rebuilds of the accumulated G.
   !>
   !> Incremental building adds a correction per iteration, so its rounding
   !> accumulates where a full build's does not. Rebuilding periodically bounds
   !> that. Sixteen is the usual choice and costs one full build in sixteen --
   !> about 6% of the saving handed back for a drift that stays at the level of
   !> the convergence threshold rather than growing past it.
   integer, parameter :: INCREMENTAL_RESET = 16

   !> Stage labels, named once so the per-iteration column and the summary table
   !> cannot drift apart, and so a caller can ask `clk%seconds_of(STAGE_FOCK)`.
   !>
   !> The split is by what a change would move. The Fock build is the integral
   !> work and the only stage anything in `INTEGRALS_MQC.md` touches; the
   !> diagonalisation is LAPACK and n^3 regardless; setup is the one-time 1e
   !> integrals, screening bounds and guess. Reporting them apart is what makes
   !> "the integrals got 20% faster" a statement about the integrals rather than
   !> about the whole run.
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
      !! From the last Fock build of the SCF, which is the converged one and so
      !! the one whose density the screening will see for the rest of any
      !! follow-on work.
      !!
      !! Split by cause because the two answer different questions. The Schwarz
      !! count is a property of the basis and the geometry and does not move
      !! between iterations; the density count is the extra reach that weighting
      !! the bound by the density buys, and it is the only figure that a change
      !! to the screening can move. Reported as counts rather than inferred from
      !! wall time: on a shared node the run-to-run spread in seconds is larger
      !! than the effect, and these are exactly reproducible.
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
      ! One is perfect balance. This is the figure a work-ordering change has to
      ! move, and unlike wall time it is a ratio inside one run, so a contended
      ! node largely divides out of it.
      write (line, "(a,f14.4)") "    thread imbalance     ", stats%thread_imbalance
      call logger%performance(trim(line))
   end subroutine screening_summary

   subroutine scf_table_header(verbose)
      !! Column headings for the per-iteration table
      !!
      !! Through the logger like everything else the program says: the level
      !! decides whether it is seen, and a run redirected to a log file gets the
      !! table with it rather than only on a terminal that is no longer there.
      !! The frame is `mqc_convergence_report`, shared with the FMO outer loop;
      !! only the columns and their widths are the SCF's own.
      logical, intent(in) :: verbose

      call convergence_header(verbose, "SCF iterations", &
                    "    iter                 energy          dE          dD   diis       Fock         XC       rest", 93)
   end subroutine scf_table_header

   subroutine scf_table_row(verbose, iter, energy, de, drms, ndiis, t_fock, t_xc, t_rest)
      !! One iteration's line, with the time that iteration took
      !!
      !! **The quadrature has a column of its own, and needs one.** `STAGE_FOCK`
      !! and `STAGE_XC` are separate buckets -- correctly, since the commit that
      !! stopped the two-electron build being billed to the quadrature -- but the
      !! row only ever read the first of them. On a Kohn-Sham run that is the
      !! smaller half by a wide margin: 0.92 s of Fock printed beside an
      !! iteration that took thirty, because the thirty were in the quadrature
      !! and the quadrature was not on the line. The totals at the end were right
      !! the whole time, which is what made it hard to see; a table nobody can
      !! reconcile against a wall clock is worse than no table.
      !!
      !! Per-iteration rather than only a total because the first iterations of a
      !! direct SCF are not the same price as the last: screening tightens as the
      !! density settles, so a run whose Fock build is flat across iterations is
      !! telling you the density-weighted screening of `INTEGRALS_MQC.md` §6.1 is
      !! missing, which is exactly the case here.
      logical, intent(in) :: verbose
      integer, intent(in) :: iter, ndiis
      real(dp), intent(in) :: energy, de, drms, t_fock, t_xc, t_rest

      character(len=LINE_LEN) :: line

      if (.not. verbose) return
      write (line, "(i8,f23.12,2es12.3,i7,3(f9.2,a))") &
         iter, energy, de, drms, ndiis, t_fock, " s", t_xc, " s", t_rest, " s"
      call logger%info(trim(line))
   end subroutine scf_table_row

   subroutine scf_table_footer(verbose, converged, iterations)
      !! Close the per-iteration table
      !!
      !! The timing summary is `clk%report(...)`, which counts calls itself --
      !! the Fock stage runs `iterations + 1` times because the converged energy
      !! is rebuilt after the loop, and dividing by `iterations` reported a
      !! per-build cost about 12% high.
      logical, intent(in) :: verbose, converged
      integer, intent(in) :: iterations

      call convergence_footer(verbose, converged, iterations, "iterations", 84)
   end subroutine scf_table_footer

   subroutine energy_components(verbose, mol, density, electronic, nuclear)
      !! The energy split into the pieces the virial theorem constrains
      !!
      !! Every term here is already available at convergence -- these are traces
      !! against the density, not new integrals -- and printing them costs one
      !! kinetic and one nuclear-attraction build. What they buy is a check the
      !! total energy cannot give: the virial ratio.
      !!
      !! For a variational wave function at a stationary geometry `-V/T` is
      !! exactly 2. It is not constrained to be, so a value that drifts from 2
      !! is evidence about the calculation rather than a restatement of it --
      !! a basis too small to describe the cusp, or a geometry that is not a
      !! stationary point. GAMESS prints the same block for the same reason.
      logical, intent(in) :: verbose
      type(libcint_molecule_t), intent(in) :: mol
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
      ! attraction alone is the difference -- no third integral build.
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

   subroutine run_libcint_rhf(mol, nelec, max_iter, energy_tol, density_tol, &
                              verbose, result, error, aux, diis_vectors, in_core, &
                              guess, guess_density, xc, h_extra, pcm, projector, &
                              level_shift, linear_dependence, b_ao_out)
      !! Drive a closed-shell SCF to convergence
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: nelec
      integer, intent(in) :: max_iter
      real(dp), intent(in) :: energy_tol, density_tol
      logical, intent(in) :: verbose
      real(dp), intent(in), optional :: level_shift
         !! Hartree added to the virtual block of the Fock matrix before each
         !! diagonalisation. Absent or zero is off.
      type(rhf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: b_ao_out(:, :)
         !! The fitted `B(mu nu, P)` this SCF used, handed on rather than
         !! freed, for a correlated step that would otherwise rebuild the
         !! integrals and refit them. Moved, not copied -- the SCF is finished
         !! with it. Only meaningful alongside `aux`.
      type(libcint_molecule_t), intent(in), optional :: aux
         !! Auxiliary basis. Present means density-fitted J and K, which is
         !! what cuEST always does -- so a comparison against the GPU path
         !! needs this, and a comparison against an exact-ERI code needs it
         !! absent. Both are worth being able to ask for.
      integer, intent(in), optional :: diis_vectors
         !! Subspace size. Zero turns DIIS off, which is worth being able to
         !! do: the two paths should agree with and without it, and if they
         !! do not, the extrapolation is where to look.
      logical, intent(in), optional :: in_core
         !! Store every integral and contract from the tensor, instead of
         !! rebuilding the Fock matrix directly each iteration. Default is
         !! direct. This exists because the in-core path is the one validated
         !! against PySCF, so it is what the direct build is checked against --
         !! not because anything should run with it. It is n^4 in memory.
      integer, intent(in), optional :: guess
         !! One of SCF_GUESS_*. Defaults to Wolfsberg-Helmholz, the best guess
         !! that needs nothing but `H` and `S` -- superposition of atomic
         !! densities is better still, and is what a deck gets, but
         !! `guess_density` has to be built before this routine is entered and
         !! so cannot be a default here.
         !!
         !! **It used to default to the core guess**, which is `F = H`: every
         !! electron-electron term ignored, and the worst starting point
         !! available. Iterations to 1e-12, measured:
         !!
         !!     HCN   STO-3G  B2GP-PLYP   core 15   gwh  9   sad 9
         !!     H2O   cc-pVDZ B3LYP       core 11   gwh  9   sad 9
         !!     NH3   cc-pVDZ PBE         core 12   gwh 10   sad 9
         !!     CH4   cc-pVDZ HF          core 10   gwh 10   sad 8
         !!
         !! Never worse, occasionally much better, and free -- `H` and `S` are
         !! both in hand before the first Fock build.
         !!
         !! This moves nothing a deck does: the driver resolves "auto" to SAD
         !! and always passes it. What it moves is every *library* caller that
         !! did not think about the guess -- the validation harnesses, SAPT, the
         !! EFP potential, the projection code -- each of which was silently
         !! taking the worst option.
      type(xc_context_t), intent(inout), optional :: xc
         !! Exchange-correlation. Present turns this into a Kohn-Sham SCF; absent
         !! leaves it Hartree-Fock. One argument rather than six, and one loop
         !! rather than two -- a separate `run_libcint_rks` would have to be kept
         !! in step with this one for the rest of its life.
      real(dp), intent(in), optional :: guess_density(:, :)
         !! Total starting density, required by SCF_GUESS_SAC and SCF_GUESS_SAD
         !! and ignored otherwise. Built by `mqc_libcint_atomic_guess`, which
         !! cannot be reached from here: it runs free-atom SCFs through this
         !! very module.
      real(dp), intent(in), optional :: h_extra(:, :)
         !! An extra one-electron operator, added to H before the first Fock
         !! build and kept there. A uniform electric field is what it is for:
         !! `F . r` makes this SCF the finite-field reference that `..._cphf`
         !! is checked against, and a response property computed two ways from
         !! one code path is worth more than either alone.
         !!
         !! `result%energy` then includes the interaction with whatever this
         !! is, which is why the finite-field check differentiates the *dipole*
         !! and not the energy: the dipole is unambiguous, one derivative
         !! better conditioned, and needs no bookkeeping about what the total
         !! is supposed to contain.
      type(pcm_context_t), intent(inout), optional :: pcm
      real(dp), intent(in), optional :: linear_dependence
         !! `keywords.scf.linear_dependence_threshold`. Zero or absent leaves
         !! the orthogonaliser on its own cutoff.
         !! Continuum solvation. Present and built means every iteration solves
         !! the surface charges against its density and folds their operator
         !! into the Fock matrix; `result%energy` then includes the dielectric
         !! energy. Present-but-disabled is the same as absent, so a caller can
         !! pass its context unconditionally.
      type(fock_projector_t), intent(in), optional :: projector
         !! A constraint on the Fock matrix, applied after every build.
         !!
         !! Unlike `h_extra` this is not an operator added to `H`: it is a
         !! linear map on `F` that forces it block diagonal in a basis the
         !! caller chose, so it cannot be folded into the core Hamiltonian and
         !! has to be reapplied each iteration. What it is for is freezing
         !! orbitals -- see `mqc_fock_projector` -- and this routine is
         !! deliberately told only that its Fock matrix is constrained, not why.
         !!
         !! Applied before the DIIS error and extrapolation, so that everything
         !! downstream sees the constrained matrix and not two different ones.
         !! Not applied to the final build below, which is the energy's Fock:
         !! the constraint decides which determinant comes out, and the energy
         !! of that determinant is an expectation value of the real Hamiltonian.

      character(len=LINE_LEN) :: line
      integer :: diis_size, guess_kind
      logical :: use_in_core
      real(dp), allocatable :: bounds(:, :)
      type(direct_stats_t) :: stats

      real(dp), allocatable :: s(:, :), h(:, :), eri(:, :, :, :), bmat(:, :)
      real(dp), allocatable :: bmat_lr(:, :)
      real(dp), allocatable :: x(:, :), fock(:, :), density(:, :), density_old(:, :)
      real(dp), allocatable :: coeff(:, :), eigenvalues(:)
      real(dp), allocatable :: err(:, :), fock_flat(:)
      real(dp), allocatable :: v_pcm(:, :)
      logical :: use_pcm
      real(dp) :: e_pcm
      type(diis_state_t) :: diis
      logical :: extrapolated
      real(dp) :: e_elec, e_old, de, drms
      real(dp) :: shift, shift_now, drms_prev, taper
      real(dp), allocatable :: sd(:, :), sds(:, :)
      integer :: n_ao, n_mo, n_occ, iter
      type(timing_report_t) :: clk
      type(direct_stats_t) :: screening
      type(incremental_state_t) :: incr
      real(dp) :: t_fock_iter, t_rest_iter, t_xc_iter
      real(dp) :: s_min, s_kept   !! overlap conditioning, reported before iteration 1
      real(dp) :: lindep

      if (mod(nelec, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "RHF needs an even electron count; this "// &
                        "system has an odd one and wants an unrestricted method")
         return
      end if

      diis_size = 8
      if (present(diis_vectors)) diis_size = diis_vectors
      use_in_core = .false.
      if (present(in_core)) use_in_core = in_core
      guess_kind = SCF_GUESS_GWH
      if (present(guess)) guess_kind = guess
      use_pcm = .false.
      if (present(pcm)) use_pcm = pcm%enabled

      n_ao = mol%nao
      n_occ = nelec/2
      if (n_occ < 1) then
         call error%set(ERROR_VALIDATION, "RHF: no electrons to place")
         return
      end if

      ! Reported because this is what went unnoticed for as long as it did:
      ! 6-31G* read as spherical converges just as prettily as 6-31G* read
      ! Cartesian -- one basis function short and 1.4 mHartree out -- and
      ! neither the iterations nor the final energy look wrong on their own.
      ! The count and the angular form together are what name the basis.
      ! Still gated on `verbose` rather than on the logger level alone: the atomic
      ! guess runs one of these per element, and those are not runs anyone wants a
      ! table for. The level decides how loud a reported run is; `verbose` decides
      ! whether this run reports at all.
      if (verbose) then
         if (mol%cartesian) then
            call logger%info("  basis functions: "//to_char(n_ao)//"  (Cartesian, 6d/10f)")
         else
            call logger%info("  basis functions: "//to_char(n_ao)//"  (spherical, 5d/7f)")
         end if
      end if

      call clk%start()
      call mol%overlap(s)
      call mol%core_hamiltonian(h)
      if (present(h_extra)) then
         if (size(h_extra, 1) /= n_ao .or. size(h_extra, 2) /= n_ao) then
            call error%set(ERROR_VALIDATION, "RHF: h_extra is not n_ao square")
            return
         end if
         h = h + h_extra
      end if
      if (present(aux)) then
         call build_df_tensor(mol, aux, bmat, error)
         if (error%has_error()) return
         ! A second fit, against the attenuated kernel, when the functional
         ! splits its exchange by range. Built once beside the first rather than
         ! per iteration: it depends on the geometry and the auxiliary basis,
         ! neither of which moves during an SCF.
         if (present(xc)) then
            if (xc%range_separated) then
               call build_df_tensor(mol, aux, bmat_lr, error, omega=xc%rs_omega)
            end if
         end if
         if (error%has_error()) return
      else if (use_in_core) then
         call mol%eris(eri)
      else
         ! The bounds depend on the basis and the geometry, not the density, so
         ! one set serves every iteration of the SCF.
         call schwarz_bounds(mol, bounds, error)
         if (error%has_error()) return
      end if

      lindep = 0.0_dp
      if (present(linear_dependence)) lindep = linear_dependence
      call build_orthogonalizer(s, x, n_mo, error, smallest_overlap=s_min, &
                                smallest_kept=s_kept, threshold=lindep)
      if (error%has_error()) return
      call report_linear_dependence(n_ao, n_mo, s_min, s_kept, verbose, threshold=lindep)
      if (n_occ > n_mo) then
         call error%set(ERROR_VALIDATION, "RHF: more occupied orbitals than the basis "// &
                        "supports after near-null modes were dropped")
         return
      end if

      allocate (fock(n_ao, n_ao), density(n_ao, n_ao), density_old(n_ao, n_ao))
      allocate (err(n_mo, n_mo), fock_flat(n_ao*n_ao))
      if (use_pcm) allocate (v_pcm(n_ao, n_ao))

      ! The error vector lives in the orthogonal basis, where it is n_mo
      ! square rather than n_ao -- that is the same shape the cuEST path uses
      ! and the reason both converge alike.
      call diis%init(diis_size, n_ao*n_ao, n_mo*n_mo)

      ! Every guess ends up as a starting Fock matrix, which is then diagonalised
      ! and occupied exactly as an iteration's would be. That uniformity is the
      ! point: the atomic guesses differ from the core guess only in what F is
      ! built from, so nothing below this line knows which one ran.
      select case (guess_kind)
      case (SCF_GUESS_CORE)
         fock = h
      case (SCF_GUESS_GWH)
         call guess_fock(s, h, fock)
      case (SCF_GUESS_SAC, SCF_GUESS_SAD, SCF_GUESS_PROJ)
         if (.not. present(guess_density)) then
            call error%set(ERROR_VALIDATION, "RHF: an atomic guess was asked for but no "// &
                           "guess density was supplied")
            return
         end if
         if (size(guess_density, 1) /= n_ao .or. size(guess_density, 2) /= n_ao) then
            call error%set(ERROR_VALIDATION, "RHF: the guess density is not the size of "// &
                           "this basis")
            return
         end if
         density = guess_density
         call atomic_guess_fock(mol, h, density, bmat, eri, bounds, fock, error)
         if (error%has_error()) return
      case default
         call error%set(ERROR_VALIDATION, "RHF: unknown initial guess")
         return
      end select

      call diagonalize(fock, x, n_ao, n_mo, coeff, eigenvalues, error)
      if (error%has_error()) return
      call build_density_closed_shell(coeff, n_occ, density)

      e_old = 0.0_dp
      result%converged = .false.

      ! Everything from the top of the routine to here is one-time cost: the 1e
      ! integrals, the screening bounds or the fitted/in-core tensor, the
      ! orthogonaliser and the guess.
      call clk%lap(STAGE_SETUP)
      call scf_table_header(verbose)

      shift = 0.0_dp
      if (present(level_shift)) shift = level_shift
      ! Refused rather than clamped. A negative shift lowers the virtuals into
      ! the occupied set, which is the opposite of the point: it narrows the gap
      ! the next density is built through and drives the oscillation it was
      ! asked to damp. Someone typing one has misunderstood the sign, and
      ! quietly using zero would hide that.
      if (shift < 0.0_dp) then
         call error%set(ERROR_VALIDATION, "keywords.scf.level_shift is negative. A "// &
                        "level shift raises the virtual orbitals to widen the gap; a "// &
                        "negative one narrows it and makes convergence worse, not "// &
                        "better. Give a positive value in Hartree, or leave it out.")
         return
      end if
      ! Off well before the SCF is done. The shift buys nothing near the
      ! solution -- it slows the last approach -- and, far more importantly,
      ! every orbital energy this routine reports has to belong to the
      ! *unshifted* operator: they are read back as MP2 and coupled-cluster
      ! denominators, as the weights of the gradient's energy-weighted density,
      ! and as `eps_occ` and the response poles of a fragment potential. A shift
      ! left in would move all of those by an amount nothing downstream could
      ! recognise as a shift.
      taper = 100.0_dp*density_tol
      drms_prev = huge(1.0_dp)
      shift_now = 0.0_dp

      ! Say so. A shift is asked for when an SCF is misbehaving, which is
      ! exactly when the run is being read closely, and until this line there
      ! was nothing between "the deck set one" and "the deck set one and it was
      ! silently dropped on the way here" -- which is what happened to the
      ! Fukui ions for as long as `fukui_indices` took no level shift. Printing
      ! the taper alongside the value matters as much: a shift that is off by
      ! iteration three is not the shift the reader thinks they applied.
      if (shift > 0.0_dp) then
         write (line, "(a,f8.4,a,es9.2)") "    level shift: ", shift, &
            " hartree, tapered off below dD ", taper
         call logger%info(trim(line))
      end if

      do iter = 1, max_iter
         density_old = density
         ! The energy belongs to the Fock built from this density, so both come
         ! back together and before extrapolation. A DIIS-mixed Fock is a
         ! convergence device, not a state anything is the energy of.
         ! Read the bucket BEFORE the build. `assemble_fock` closes the Fock
         ! stage itself, so by the time it returns the time is already charged
         ! and a read here would difference against itself -- which is what the
         ! per-iteration Fock column was doing when it printed 0.00 s beside a
         ! Kohn-Sham run that plainly took seconds.
         t_fock_iter = clk%seconds_of(STAGE_FOCK)
         t_xc_iter = clk%seconds_of(STAGE_XC)
         call assemble_fock(mol, h, density, coeff, n_occ, bmat, eri, bounds, xc, &
                            fock, e_elec, error, clk=clk, screening=screening, incr=incr, &
                            bmat_lr=bmat_lr)
         if (error%has_error()) return
         ! The continuum, after the gas-phase pieces and before the commutator:
         ! its operator belongs to the Fock matrix DIIS extrapolates and the
         ! eigenproblem diagonalises, or the converged orbitals would not feel
         ! the solvent. The energy carries its own half already.
         !
         ! Charged to its OWN stage, not folded into the Fock bucket. `lap`
         ! adds to a stage's seconds and increments its call count, and
         ! `assemble_fock` already closes STAGE_FOCK itself -- so lapping it
         ! again here would report two Fock builds per iteration, which is the
         ! reporting bug the Fock/quadrature split was fixed for. A separate row
         ! also says what the continuum actually costs, which turns out to be
         ! almost nothing beside the quadrature.
         if (use_pcm) then
            call pcm%operator_matrix(mol, density, v_pcm, e_pcm, error)
            if (error%has_error()) return
            fock = fock + v_pcm
            e_elec = e_elec + e_pcm
            call clk%lap(STAGE_PCM)
         end if
         t_fock_iter = clk%seconds_of(STAGE_FOCK) - t_fock_iter
         t_xc_iter = clk%seconds_of(STAGE_XC) - t_xc_iter

         if (present(projector)) then
            call projector%apply(fock, error)
            if (error%has_error()) return
         end if

         ! FDS - SDF, projected into the orthogonal basis. It vanishes exactly
         ! when F and D commute, which is what convergence means, so it is the
         ! quantity worth extrapolating against.
         call commutator(fock, density, s, x, err)
         fock_flat = reshape(fock, [n_ao*n_ao])
         call diis%push(fock_flat, reshape(err, [n_mo*n_mo]))
         call diis%extrapolate(fock_flat, extrapolated)
         if (extrapolated) fock = reshape(fock_flat, [n_ao, n_ao])
         t_rest_iter = clk%seconds_of(STAGE_DIIS)
         call clk%lap(STAGE_DIIS)
         t_rest_iter = clk%seconds_of(STAGE_DIIS) - t_rest_iter

         ! **After DIIS, never before it.** The error vector above is built from
         ! the unshifted Fock, and it has to be: shift first and the vectors
         ! DIIS stores stop being a subspace of Fock matrices, so what it
         ! extrapolates is not the thing it is meant to extrapolate. It would
         ! still converge to something, which is what makes that ordering
         ! expensive to find. The energy is likewise already taken, from the
         ! build above, for the reason stated there.
         !
         ! The virtual projector without the virtual orbitals: completeness
         ! gives `C_o C_o^T + C_v C_v^T = S^-1`, so
         ! `S C_v C_v^T S = S - S C_o C_o^T S = S - (1/2) S D S` for a closed
         ! shell, which is two products against the density already in hand.
         shift_now = 0.0_dp
         if (shift > 0.0_dp .and. drms_prev > taper) shift_now = shift
         if (shift_now > 0.0_dp) then
            if (.not. allocated(sd)) allocate (sd(n_ao, n_ao), sds(n_ao, n_ao))
            call pic_gemm(s, density, sd, beta=0.0_dp)
            call pic_gemm(sd, s, sds, beta=0.0_dp)
            fock = fock + shift_now*(s - 0.5_dp*sds)
         end if

         call diagonalize(fock, x, n_ao, n_mo, coeff, eigenvalues, error)
         if (error%has_error()) return
         call build_density_closed_shell(coeff, n_occ, density)
         t_rest_iter = t_rest_iter - clk%seconds_of(STAGE_DIAG)
         call clk%lap(STAGE_DIAG)
         t_rest_iter = t_rest_iter + clk%seconds_of(STAGE_DIAG)

         de = abs(e_elec - e_old)
         drms = sqrt(sum((density - density_old)**2)/real(n_ao*n_ao, dp))
         call scf_table_row(verbose, iter, e_elec + mol%nuclear_repulsion(), de, drms, &
                            diis%count(), t_fock_iter, t_xc_iter, t_rest_iter)

         e_old = e_elec
         result%iterations = iter
         ! `shift_now` is part of the test rather than checked afterwards. The
         ! orbitals and eigenvalues that leave here are the ones the last
         ! diagonalisation produced, so convergence has to be declared on an
         ! iteration that was not shifted -- otherwise every virtual energy is
         ! high by `shift` and everything downstream quietly inherits it. The
         ! taper makes this true on its own in every ordinary case; requiring it
         ! costs an iteration in the case where it would not have been.
         drms_prev = drms
         if (iter > 1 .and. de < energy_tol .and. drms < density_tol .and. &
             shift_now == 0.0_dp) then
            result%converged = .true.
            exit
         end if
      end do

      ! The energy that goes out is the one belonging to the density that
      ! satisfied the test, so it is recomputed from the final Fock rather
      ! than carried over from the loop.
      !
      ! Deliberately without `incr`: this one is a full build. An incremental
      ! Fock carries however many corrections have accumulated since the last
      ! reset, and the number that leaves this routine should not.
      call assemble_fock(mol, h, density, coeff, n_occ, bmat, eri, bounds, xc, &
                         fock, result%electronic, error, clk=clk, bmat_lr=bmat_lr)
      if (error%has_error()) return
      ! Last use of the fitted tensor. A caller that asked for it takes the
      ! allocation itself rather than a copy of several gigabytes.
      if (present(b_ao_out) .and. allocated(bmat)) call move_alloc(bmat, b_ao_out)
      if (use_pcm) then
         call pcm%operator_matrix(mol, density, v_pcm, e_pcm, error)
         if (error%has_error()) return
         result%electronic = result%electronic + e_pcm
         result%pcm_energy = e_pcm
         result%pcm_charge = pcm%q_total
         ! Lapped like the iterations' operators, so the stage counts one per
         ! Fock build rather than one fewer. The final rebuild is a solvated
         ! build like any other -- it has to be, or the energy reported would
         ! be the gas-phase one.
         call clk%lap(STAGE_PCM)
      end if
      result%nuclear_repulsion = mol%nuclear_repulsion()
      result%energy = result%electronic + result%nuclear_repulsion
      result%n_occupied = n_occ
      call move_alloc(eigenvalues, result%orbital_energies)
      call move_alloc(coeff, result%orbitals)
      call move_alloc(density, result%density)
      call diis%destroy()

      ! The final rebuild above is a Fock build like any other and lands in the
      ! same bucket -- charged by `assemble_fock` itself, so no lap here. One
      ! used to be, and once the routine started closing its own stage that
      ! second lap only inflated the call count: 24 builds reported for a
      ! 12-iteration SCF.
      call clk%finish()
      call scf_table_footer(verbose, result%converged, result%iterations)
      call energy_components(verbose, mol, result%density, result%electronic, &
                             result%nuclear_repulsion)
      call screening_summary(verbose, screening)
      call clk%report("RHF", verbose)
   end subroutine run_libcint_rhf

   subroutine run_libcint_uhf(mol, nelec, multiplicity, max_iter, energy_tol, density_tol, &
                              verbose, result, error, diis_vectors, in_core, diis_start, &
                              guess, guess_density_alpha, guess_density_beta, xc, pcm, &
                              level_shift, linear_dependence)
      !! Drive an unrestricted SCF to convergence
      !!
      !! Two Fock matrices sharing one Coulomb field: F_sigma = H + J(D_a + D_b)
      !! - K(D_sigma). They are diagonalised separately and extrapolated
      !! together -- one DIIS over the pair, because the two spins converge to
      !! each other and letting them extrapolate independently lets them
      !! disagree about which iteration they are on.
      !!
      !! No density fitting here yet. The fitted J and K are written for one
      !! density, and splitting them is a separate piece of work rather than an
      !! argument change; asking for both is refused rather than silently run
      !! restricted.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: nelec
      integer, intent(in) :: multiplicity   !! 2S+1
      integer, intent(in) :: max_iter
      real(dp), intent(in) :: energy_tol, density_tol
      logical, intent(in) :: verbose
      real(dp), intent(in), optional :: level_shift
         !! As `run_libcint_rhf`. Applied to each spin's Fock matrix against its
         !! own virtual projector, and tapered off before convergence for the
         !! same reason.
      type(rhf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: diis_vectors
      logical, intent(in), optional :: in_core
      integer, intent(in), optional :: diis_start
         !! First iteration allowed to extrapolate. See the default below.
      integer, intent(in), optional :: guess
         !! One of SCF_GUESS_*. Defaults to Wolfsberg-Helmholz, as the restricted
         !! path does and for the same reason -- see the note there.
      type(xc_context_t), intent(inout), optional :: xc
         !! Exchange-correlation, spin-polarised. Present makes this an
         !! unrestricted Kohn-Sham SCF; absent leaves it unrestricted
         !! Hartree-Fock. The context must have been built with polarized=.true.,
         !! which `xc_add_potential_uks` checks rather than assumes -- libxc fixes
         !! the spin channel at initialisation and a restricted context would
         !! otherwise be read with the wrong stride.
      real(dp), intent(in), optional :: guess_density_alpha(:, :)
      real(dp), intent(in), optional :: guess_density_beta(:, :)
         !! Starting spin densities, required by SCF_GUESS_SAC and SCF_GUESS_SAD.
         !! Both spins are taken separately on purpose: whether they differ is
         !! the whole difference between the two atomic guesses here. SAC hands
         !! over the free atom's own spin-polarised densities and so arrives
         !! already symmetry-broken, which is what gets a radical's sigma/pi
         !! ordering right. SAD hands over half the spherically averaged total
         !! twice, and relies on the occupations to break the symmetry -- the
         !! same position the core guess is in, but from a far better density.
      type(pcm_context_t), intent(inout), optional :: pcm
      real(dp), intent(in), optional :: linear_dependence
         !! `keywords.scf.linear_dependence_threshold`. Zero or absent leaves
         !! the orthogonaliser on its own cutoff.
         !! Continuum solvation, exactly as on the restricted path. The surface
         !! charges see the total density and their operator lands in both spin
         !! Fock matrices, because the potential of a classical charge does not
         !! know spin.

      integer :: diis_size, start_cycle, guess_kind
      logical :: use_in_core
      real(dp), allocatable :: bounds(:, :)
      type(direct_stats_t) :: stats

      real(dp), allocatable :: s(:, :), h(:, :), eri(:, :, :, :)
      real(dp), allocatable :: x(:, :), fock_a(:, :), fock_b(:, :)
      real(dp), allocatable :: d_a(:, :), d_b(:, :), d_a_old(:, :), d_b_old(:, :)
      real(dp), allocatable :: coeff_a(:, :), coeff_b(:, :), eig_a(:), eig_b(:)
      real(dp), allocatable :: err_a(:, :), err_b(:, :), fock_flat(:), err_flat(:)
      real(dp), allocatable :: v_pcm(:, :)
      logical :: use_pcm
      real(dp) :: e_pcm
      type(diis_state_t) :: diis
      logical :: extrapolated
      real(dp) :: e_elec, e_old, de, drms
      real(dp) :: shift, shift_now, drms_prev, taper
      real(dp), allocatable :: sd(:, :), sds(:, :)
      integer :: n_ao, n_mo, n_alpha, n_beta, iter, nsq, msq
      type(timing_report_t) :: clk
      real(dp) :: t_fock_iter, t_rest_iter, t_xc_iter
      real(dp) :: s_min, s_kept   !! overlap conditioning, reported before iteration 1
      real(dp) :: lindep
      character(len=LINE_LEN) :: line

      diis_size = 8
      if (present(diis_vectors)) diis_size = diis_vectors
      start_cycle = DEFAULT_UHF_DIIS_START
      if (present(diis_start)) start_cycle = diis_start
      use_in_core = .false.
      if (present(in_core)) use_in_core = in_core
      guess_kind = SCF_GUESS_GWH
      if (present(guess)) guess_kind = guess
      use_pcm = .false.
      if (present(pcm)) use_pcm = pcm%enabled

      ! 2S+1 = n_alpha - n_beta + 1, so the unpaired count is multiplicity - 1.
      ! Both have to come out whole and non-negative, and a deck can ask for a
      ! pairing that does not exist -- four electrons in a quintet, say.
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
      allocate (fock_flat(2*nsq), err_flat(2*msq))
      if (use_pcm) allocate (v_pcm(n_ao, n_ao))

      ! One subspace over both spins, so an extrapolation step moves them
      ! together. The vectors are the two Fock matrices laid end to end and the
      ! two commutators likewise.
      call diis%init(diis_size, 2*nsq, 2*msq)

      ! The symmetric guesses -- core and GWH -- give alpha and beta the same
      ! orbitals, and the occupations are what separate them: n_alpha > n_beta
      ! puts an electron somewhere beta has none, and the Fock matrices diverge
      ! from there. For n_alpha == n_beta they converge to the restricted
      ! solution, which is the right answer rather than a failure. SAC is the one
      ! that arrives already asymmetric, because a free atom's own alpha and beta
      ! densities differ wherever its shell is open.
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
      call scf_table_header(verbose)

      shift = 0.0_dp
      if (present(level_shift)) shift = level_shift
      ! Refused rather than clamped. A negative shift lowers the virtuals into
      ! the occupied set, which is the opposite of the point: it narrows the gap
      ! the next density is built through and drives the oscillation it was
      ! asked to damp. Someone typing one has misunderstood the sign, and
      ! quietly using zero would hide that.
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

      ! Say so. A shift is asked for when an SCF is misbehaving, which is
      ! exactly when the run is being read closely, and until this line there
      ! was nothing between "the deck set one" and "the deck set one and it was
      ! silently dropped on the way here" -- which is what happened to the
      ! Fukui ions for as long as `fukui_indices` took no level shift. Printing
      ! the taper alongside the value matters as much: a shift that is off by
      ! iteration three is not the shift the reader thinks they applied.
      if (shift > 0.0_dp) then
         write (line, "(a,f8.4,a,es9.2)") "    level shift: ", shift, &
            " hartree, tapered off below dD ", taper
         call logger%info(trim(line))
      end if

      do iter = 1, max_iter
         d_a_old = d_a
         d_b_old = d_b

         ! Read before the build, not after it. Both stages are cumulative, so a
         ! reading taken once the build has already charged them differences
         ! against itself -- which is what the XC column was doing.
         t_fock_iter = clk%seconds_of(STAGE_FOCK)
         t_xc_iter = clk%seconds_of(STAGE_XC)
         call assemble_fock_uhf(mol, h, d_a, d_b, eri, bounds, xc, fock_a, fock_b, &
                                e_elec, error, clk=clk)
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
         call clk%lap(STAGE_FOCK)
         t_fock_iter = clk%seconds_of(STAGE_FOCK) - t_fock_iter
         t_xc_iter = clk%seconds_of(STAGE_XC) - t_xc_iter

         call commutator(fock_a, d_a, s, x, err_a)
         call commutator(fock_b, d_b, s, x, err_b)
         fock_flat(1:nsq) = reshape(fock_a, [nsq])
         fock_flat(nsq + 1:2*nsq) = reshape(fock_b, [nsq])
         err_flat(1:msq) = reshape(err_a, [msq])
         err_flat(msq + 1:2*msq) = reshape(err_b, [msq])
         call diis%push(fock_flat, err_flat)
         extrapolated = .false.
         if (iter >= start_cycle) call diis%extrapolate(fock_flat, extrapolated)
         if (extrapolated) then
            fock_a = reshape(fock_flat(1:nsq), [n_ao, n_ao])
            fock_b = reshape(fock_flat(nsq + 1:2*nsq), [n_ao, n_ao])
         end if
         t_rest_iter = clk%seconds_of(STAGE_DIIS)
         call clk%lap(STAGE_DIIS)
         t_rest_iter = clk%seconds_of(STAGE_DIIS) - t_rest_iter

         ! After the extrapolation and after the two commutators, for the reason
         ! `run_libcint_rhf` sets out. Each spin gets its own virtual projector:
         ! `build_density_spin` gives an occupation of one, so the closed-shell
         ! factor of a half is absent and it is `S - S D_sigma S`.
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
         call scf_table_row(verbose, iter, e_elec + mol%nuclear_repulsion(), de, drms, &
                            diis%count(), t_fock_iter, t_xc_iter, t_rest_iter)

         e_old = e_elec
         result%iterations = iter
         drms_prev = drms
         if (iter > 1 .and. de < energy_tol .and. drms < density_tol .and. &
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
   end subroutine run_libcint_uhf

   subroutine assemble_fock(mol, h, density, coeff, n_occ, bmat, eri, bounds, xc, &
                            fock, e_elec, error, clk, screening, incr, bmat_lr)
      !! The Fock matrix for this density, and the electronic energy that belongs to it
      !!
      !! One place, because there used to be two: the iteration built its Fock
      !! matrix through a three-way branch and the final energy rebuilt it through
      !! an identical copy, which is a duplication that only stays correct by
      !! attention. Exchange-correlation would have made it three copies.
      !!
      !! The energy is returned with the Fock matrix because for Kohn-Sham they
      !! cannot be separated after the fact. Hartree-Fock's
      !! `1/2 Tr[D (H + F)]` counts the potential at half weight, which is right
      !! for a mean field and wrong for a functional: the exchange-correlation
      !! energy is E_xc, not `1/2 Tr[D V_xc]`. So the energy is taken from the Fock
      !! matrix *before* V_xc is added, and E_xc is added on top.
      type(libcint_molecule_t), intent(in) :: mol
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
         !! How much of the quartet space each test removed. Reported rather than
         !! discarded because wall time on a shared node varies by more than a
         !! screening change does -- these counts are exactly reproducible, so
         !! they are the measurement that can actually be trusted.
      type(incremental_state_t), intent(inout), optional :: incr
         !! Present from the SCF loop, which wants each build to cost only the
         !! change since the last one. Absent from the guess, the gradient callers
         !! and the final energy rebuild, all of which want an exact build --
         !! the last of those especially, since the energy that goes out must not
         !! carry sixteen iterations of accumulated correction.

      type(direct_stats_t) :: stats
      real(dp), allocatable :: g_delta(:, :), d_delta(:, :), h_zero(:, :)
      logical :: use_incremental, full_build
      real(dp), allocatable :: v_xc(:, :)
      real(dp) :: k_scale, e_xc, n_elec
      logical :: have_lr
      logical :: kohn_sham

      kohn_sham = .false.
      k_scale = 1.0_dp
      if (present(xc)) then
         if (xc%active) then
            kohn_sham = .true.
            ! Pure functionals want no Fock exchange at all, hybrids want a
            ! fraction, and Hartree-Fock is the fraction being one -- which is why
            ! this is a scale rather than a branch.
            k_scale = xc%exx_fraction
         end if
      end if

      ! A range-separated functional needs a second exchange matrix, over the
      ! erf-attenuated kernel, and the two are combined as
      !
      !     K_eff = exx_fraction * K_full + rs_k_lr * K_lr(omega)
      !
      ! Only the direct build can produce the second one: libcint switches kernels
      ! through `env`, so an omega pass is the same quartet loop, while the in-core
      ! tensor and the fitted tensor are both built for the full Coulomb kernel and
      ! would need a second one of their own. Refused rather than approximated.
      if (present(xc)) then
         ! The in-core tensor is still refused: it is built once for the full
         ! Coulomb kernel and a second `n_ao^4` array for the attenuated one is
         ! not worth having. The *fitted* path is answerable now -- `bmat_lr`
         ! below is the same fit against `erf(omega r)/r`, which costs another
         ! `n_ao^2 n_aux` rather than another `n_ao^4`.
         if (xc%range_separated .and. allocated(eri) .and. .not. allocated(bmat)) then
            call error%set(ERROR_VALIDATION, "a range-separated functional needs the "// &
                           "direct or the density-fitted Fock build: the in-core "// &
                           "integrals are built for the full Coulomb kernel, so the "// &
                           "long-range exchange would be missing. Run without in_core.")
            return
         end if
         ! `present` before `allocated`, and not merged into one expression:
         ! Fortran does not short-circuit, so `allocated(bmat_lr)` on an absent
         ! optional is a reference to something that is not there. The same
         ! trap the two-particle density hit in the MP2 gradient.
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

      ! Incremental building is for the direct path only. The fitted and in-core
      ! paths contract a stored tensor, so their cost does not depend on how large
      ! the density elements are and there is nothing for a small delta to save.
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
               incr%since_reset = 0
               incr%active = .true.
            else
               allocate (g_delta(size(h, 1), size(h, 2)), d_delta(size(h, 1), size(h, 2)))
               allocate (h_zero(size(h, 1), size(h, 2)))
               d_delta = density - incr%d_ref
               h_zero = 0.0_dp
               ! Zero for the core Hamiltonian, so this returns G(delta) and nothing
               ! has to be subtracted back out. The same device the range-separated
               ! branch below uses to get the long-range exchange on its own.
               !
               ! `density_screen=.true.`, and it is the whole point of building
               ! incrementally. With it false this branch screened exactly as
               ! hard as a full build -- same Schwarz bound, same 15 % skipped --
               ! so the delta cost what it was meant to save: iteration twelve
               ! was as expensive as iteration one, and the incremental
               ! machinery bought nothing. With it true the skipped fraction
               ! goes to 28 % and the per-iteration cost falls 1.87 s to 1.24 s
               ! across the SCF, which is the shape GAMESS's FDIFF shows too.
               !
               ! The worry this replaces was that a screen keyed on the
               ! magnitude of what a quartet multiplies tightens as the delta
               ! shrinks, leaving a correction that is systematically incomplete
               ! and path-dependent -- the symptom being a converged energy that
               ! moves when the whole system is translated. It does not happen,
               ! because the test is on the *contribution* against a fixed
               ! absolute tolerance, not on the delta relative to itself: a term
               ! dropped here is smaller than the tolerance in absolute terms
               ! however small the delta got, and `INCREMENTAL_RESET` bounds how
               ! many such drops can accumulate before a full build re-syncs.
               !
               ! Measured rather than argued. Translating this molecule by
               ! (37, -19, 23) Angstrom moves the converged energy by 7e-13, the
               ! energy sits 4e-12 from the same SCF run in core, and all 235
               ! CPU validation cases pass against their PySCF references with
               ! the in-core path forced off so that this branch is what runs.
               !
               ! Note this is *not* the CPHF situation, which still needs false.
               ! There the density is a Krylov trial vector the solver drives to
               ! zero and the operator has to stay the same linear map from one
               ! matvec to the next; here there is no operator to keep fixed,
               ! only a sum to accumulate to a bounded accuracy.
               call build_fock_direct(mol, h_zero, d_delta, bounds, g_delta, stats, error, &
                                      k_scale=k_scale, density_screen=.true.)
               if (error%has_error()) return
               incr%g_ref = incr%g_ref + g_delta
               incr%d_ref = density
               incr%since_reset = incr%since_reset + 1
               fock = h + incr%g_ref
               deallocate (g_delta, d_delta, h_zero)
            end if
         else
            call build_fock_direct(mol, h, density, bounds, fock, stats, error, &
                                   k_scale=k_scale)
            if (error%has_error()) return
         end if

         ! The long-range exchange is rebuilt from the full density every
         ! iteration rather than incrementally. `g_ref` is captured before this is
         ! added, so the two do not overlap -- but it does mean a range-separated
         ! functional saves on only one of its two passes.
         if (kohn_sham) then
            if (xc%range_separated) then
               block
                  real(dp), allocatable :: k_lr(:, :)
                  real(dp), allocatable :: h_zero(:, :)
                  ! Zero in place of the core Hamiltonian and no Coulomb term, so
                  ! this pass returns the long-range exchange alone and nothing has
                  ! to be subtracted back out afterwards.
                  allocate (k_lr(size(h, 1), size(h, 2)), h_zero(size(h, 1), size(h, 2)))
                  h_zero = 0.0_dp
                  call build_fock_direct(mol, h_zero, density, bounds, k_lr, stats, &
                                         error, k_scale=xc%rs_k_lr, j_scale=0.0_dp, &
                                         omega=xc%rs_omega)
                  if (error%has_error()) return
                  fock = fock + k_lr
                  deallocate (k_lr, h_zero)
               end block
            end if
         end if
      end if

      ! Taken here, from the Fock matrix that is still a mean field.
      e_elec = electronic_energy(h, fock, density)

      ! CLOSE THE FOCK BUCKET BEFORE THE QUADRATURE OPENS.
      !
      ! `lap` charges everything since the previous lap to the stage it names,
      ! and the previous lap is in the caller, before this routine was entered.
      ! So a single `lap(STAGE_XC)` after the quadrature charged the two-electron
      ! build to the quadrature as well, and the caller's `lap(STAGE_FOCK)`
      ! immediately afterwards found nothing left to charge.
      !
      ! The report said so plainly and was read as a result rather than as a
      ! bug: 14 Fock builds at 0.00 s beside 81.57 s of "XC quadrature" on a
      ! 518-function PBE0 run. A Kohn-Sham profile that shows no Fock time at
      ! all is reporting on its own instrumentation.
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
                                fock_a, fock_b, e_elec, error, clk)
      !! Both spin Fock matrices for this pair of densities, and their energy
      !!
      !! The unrestricted twin of `assemble_fock`, and it exists for the same
      !! reason: the loop, the initial guess and the final rebuild each built the
      !! Fock matrices through their own copy of the same branch, so
      !! exchange-correlation would have become a fourth. The energy comes back with
      !! the matrices because for Kohn-Sham it cannot be recovered afterwards --
      !! `uhf_electronic_energy` counts the potential at half weight, which is right
      !! for a mean field and wrong for a functional, so it is taken before V_xc is
      !! added and E_xc is added on top.
      type(libcint_molecule_t), intent(in) :: mol
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

      type(direct_stats_t) :: stats
      real(dp), allocatable :: v_a(:, :), v_b(:, :)
      real(dp) :: k_scale, e_xc, n_elec
      logical :: kohn_sham
      integer :: n

      n = size(h, 1)
      kohn_sham = .false.
      k_scale = 1.0_dp
      if (present(xc)) then
         if (xc%active) then
            kohn_sham = .true.
            k_scale = xc%exx_fraction
         end if
      end if

      ! As in the closed-shell path: the attenuated exchange matrix is something
      ! only the direct build can produce, so it is refused rather than dropped.
      if (present(xc)) then
         if (xc%range_separated .and. allocated(eri)) then
            call error%set(ERROR_VALIDATION, "a range-separated functional needs the "// &
                           "direct Fock build: the in-core integrals are built for the "// &
                           "full Coulomb kernel, and the long-range exchange would be "// &
                           "missing. Run without in_core.")
            return
         end if
      end if

      if (allocated(eri)) then
         call build_fock_uhf(h, eri, d_alpha, d_beta, fock_a, fock_b, k_scale=k_scale)
      else
         call build_fock_direct_uhf(mol, h, d_alpha, d_beta, bounds, fock_a, fock_b, &
                                    stats, error, k_scale=k_scale)
         if (error%has_error()) return
         if (kohn_sham) then
            if (xc%range_separated) then
               block
                  real(dp), allocatable :: k_lr_a(:, :), k_lr_b(:, :), h_zero(:, :)
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
               end block
            end if
         end if
      end if

      ! Taken here, while the Fock matrices are still a mean field.
      e_elec = uhf_electronic_energy(h, fock_a, fock_b, d_alpha, d_beta)

      ! Close the Fock bucket before the quadrature opens, for the reason
      ! `assemble_fock` sets out at length. The unrestricted path had neither
      ! lap: it was never handed a clock, so `STAGE_XC` was never charged at
      ! all and the caller's XC column differenced a constant against itself.
      ! It printed 0.00 s next to a grid-4 PBE run, which is not a small number
      ! reported imprecisely but no measurement at all.
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
      !! Projected into the orthogonal basis rather than left in the AO one:
      !! in the AO basis the commutator's size depends on the overlap's
      !! conditioning, so two bases describing the same molecule would
      !! converge differently for no physical reason.
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
      !! nothing beyond two matrices already in hand. It is worth having even
      !! next to the atomic guesses: it needs no free-atom SCF, so it is the
      !! fallback when one fails, and it is what cuEST starts from -- which makes
      !! it the only guess under which the two backends' iteration counts are
      !! comparable.
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
      !! Exists so the atomic guesses reach the loop in the same state as the
      !! core guess: as a Fock matrix waiting to be diagonalised. The
      !! alternative -- entering the loop with the guess density directly -- was
      !! rejected because the DF path would then have no orbitals for its
      !! exchange on the one iteration that has no orbitals to give it.
      type(libcint_molecule_t), intent(in) :: mol
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
      !! from the density, which is a factor of n/n_occ cheaper and the reason
      !! the fitted path is worth having. A guess density cannot be written that
      !! way -- superposing atomic densities, or spherically averaging one, gives
      !! something that is not idempotent, so it has no orbitals in the usual
      !! sense. Its eigenvectors serve instead: with D = sum_i w_i v_i v_i^T and
      !! c_i = v_i sqrt(w_i/2), the identity D = 2 sum_i c_i c_i^T holds exactly,
      !! and K comes out exact through code that never learns the difference.
      !!
      !! Costs one n^3 diagonalisation, once, on the guess only.
      real(dp), intent(in) :: density(:, :)
      real(dp), allocatable, intent(out) :: coeff(:, :)
      integer, intent(out) :: n_modes
      type(error_t), intent(inout) :: error

      !> Occupations below this are dropped. A density built from converged
      !> atomic solutions is positive semi-definite in exact arithmetic, so
      !> anything at or below zero here is rounding on an empty mode -- an
      !> unoccupied polarisation shell, most often -- and contributes nothing to
      !> K. Keeping them would mean taking the square root of a negative number.
      real(dp), parameter :: OCCUPATION_FLOOR = 1.0e-12_dp
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
      !! Threaded over the target element, which needs no reduction: each `(mu, nu)`
      !! accumulates into its own entry. That is worth doing because the response
      !! solver applies this thousands of times -- once per iteration of an inner
      !! solve inside every iteration of an outer one -- so it, and not the integrals,
      !! is where an in-core run spends its time. Not `pure` any more for the sake of
      !! the directive.
      !!
      !! Correct for an antisymmetric density as well as a symmetric one: nothing
      !! here assumes either, and the Coulomb term vanishes of its own accord when
      !! the density is antisymmetric because the integral is symmetric in the
      !! contracted pair.
      !!
      !! **Blocked over the ket pair, not over the target element.** Writing the
      !! sum with `(mu, nu)` outermost reads the exchange integral as
      !! `eri(mu, la, nu, si)` with `la` innermost -- a stride of `n` doubles, so
      !! every one of the n^4 loads pulls a cache line to use eight bytes of it.
      !! The tensor is 380 MB at eighty-three functions and the loop was touching
      !! it as though it were eight times that.
      !!
      !! Fixing the order costs nothing in arithmetic. For a fixed ket pair
      !! `(c, d)` the slice `eri(:, :, c, d)` is contiguous, and *both* terms are
      !! expressible on it:
      !!
      !!   * Coulomb wants `J += D(c,d) * eri(:,:,c,d)` -- an axpy over the block.
      !!   * Exchange wants `K(:,c) += eri(:,:,c,d) . D(:,d)` -- a matrix-vector
      !!     product on the same block. Relabel `(mu la|nu si)` as `(a b|c d)`:
      !!     `K(a,c) = sum_bd D(b,d) (a b|c d)`, which is that gemv summed over d.
      !!
      !! So the tensor is walked once, in order, and each 55 kB block is used
      !! twice while it is still in L2. Same n^4 flops, an eighth of the traffic.
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

      ! Threaded over `c`, which is what makes the exchange update safe without a
      ! reduction: `K(:, c)` belongs to exactly one thread. Coulomb accumulates
      ! over every `(c, d)` and so does need one, but it is a single n*n array per
      ! thread merged once, not per-element atomics.
      !
      ! The `if` keeps genuinely tiny systems serial. The response solver enters
      ! this region tens of thousands of times, so for a handful of basis functions
      ! the barriers cost more than the arithmetic; by a few dozen they do not.
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
      !! Neither term ever forms a four-index object. Density fitting is not the
      !! opposite of a direct build -- it removes the four-index integrals
      !! altogether, so the direct/stored distinction does not apply to it. What
      !! is stored is B, which is n^2 by n_aux: 40 MB where the exact tensor is
      !! 800 MB, and n^3 rather than n^4.
      !!
      !! J is two contractions with the density: c_P = sum_uv B(uv,P) D_uv, then
      !! J_uv = sum_P B(uv,P) c_P. Both are BLAS-2 against the flattened tensor,
      !! blocked over P so each thread owns a slice of the output and nothing
      !! but the final merge is shared.
      !!
      !! K goes through the occupied orbitals rather than the density. Writing
      !! D = 2 sum_i C_ui C_vi and substituting,
      !!
      !!    K_uv = sum_P sum_ls B(ul,P) D_ls B(sv,P)
      !!         = 2 sum_P sum_i W(u,i,P) W(v,i,P),   W^P = B^P C_occ
      !!
      !! which costs 2 n^2 n_occ per auxiliary function instead of the 2 n^3 of
      !! contracting the full density. The saving is a factor of n/n_occ, and it
      !! grows with basis size at fixed electron count -- which is the regime
      !! density fitting is used in. On water/cc-pVDZ that is already 4.8x.
      !!
      !! **Threaded over P, not collapsed into one GEMM.** Stacking every W^P
      !! side by side turns sum_P W^P (W^P)^T into a single large product, which
      !! is the right move against a threaded BLAS and the wrong one here: this
      !! project links a sequential BLAS on purpose, so one GEMM is one core no
      !! matter how large it is. Measured, that version ran the exchange 2.2x
      !! faster on sixteen threads where the per-P loop below runs it 12x. An
      !! auxiliary function owns a private partial K instead, and the partials
      !! are merged once at the end.
      !!
      !! Every loop here is threaded. It was written as three serial loops over
      !! n_aux with a `reshape` of a column of B inside each, which is where a
      !! fitted SCF spent nearly all of its time: at 560 functions and n_aux
      !! near 2000 that is six thousand 2.5 MB array temporaries per iteration,
      !! evaluated on one core while the rest of the node idled. The reshapes
      !! are gone -- a column of B reaches the GEMM through an explicit-shape
      !! dummy, which is sequence association and copies nothing.
      real(dp), intent(in) :: h(:, :), density(:, :)
      real(dp), intent(in) :: b(:, :)
      contiguous :: b
         !! So that a column block `b(:, p0:p1)` reaches BLAS as a view, and
         !! `b(:, p)` reaches `df_exchange_slice` as one. Without it the
         !! compiler is entitled to copy 2.5 MB per auxiliary function --
         !! reinstating, per call, the very temporaries this removed.
         !! On its own line because `fortitude` mis-scopes the argument list
         !! when `contiguous` shares a declaration with `intent`.
      real(dp), intent(in) :: coeff(:, :)   !! MO coefficients; only the occupied block is read
      integer, intent(in) :: n_occ
      real(dp), intent(out) :: fock(:, :)
      real(dp), intent(in), optional :: k_scale   !! Exact-exchange fraction, default one
      real(dp), intent(in), optional :: j_scale
         !! Coulomb fraction, default one. Zero is the long-range exchange pass
         !! of a range-separated functional: the attenuated tensor must not
         !! contribute a second Coulomb term, because `erf(omega r)/r` and
         !! `1/r` differ there too and J is already complete from the full-range
         !! pass.

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
      ! long-range exchange pass of a range-separated functional: two full
      ! passes over B is not a cheap way to compute something discarded.
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
      !! Exists so that a column of B can be seen as the (n, n) matrix it
      !! already is. `b_p` is explicit-shape against a contiguous rank-1 actual
      !! argument, so this is sequence association: the reinterpretation is
      !! free, where the `reshape` it replaces materialised the whole block --
      !! 2.5 MB at 560 functions, once per auxiliary function, per iteration.
      !!
      !! Both halves in one place because `w` is scratch that never leaves:
      !! W^P is formed and consumed here, so nothing outside needs the stack.
      integer, intent(in) :: n
      real(dp), intent(in) :: b_p(n, n)
      real(dp), intent(in) :: c_occ(:, :)
      real(dp), intent(inout) :: w(:, :)
      real(dp), intent(inout) :: k_local(:, :)

      call pic_gemm(b_p, c_occ, w, beta=0.0_dp)
      ! A rank-k update, not a general product: W W^T is symmetric, so half of
      ! what the GEMM computed was the other half's mirror image. Only the
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
         !! Fraction of exact exchange, as in the closed-shell build: one for
         !! Hartree-Fock and the default, less for a hybrid Kohn-Sham build.

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

end module mqc_libcint_rhf
