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
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_timing, only: timing_report_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_diis, only: diis_state_t
   use mqc_libcint_direct, only: schwarz_bounds, build_fock_direct, build_fock_direct_uhf, &
                                 direct_stats_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_df_tensor
   use mqc_libcint_xc, only: xc_context_t, xc_add_potential, xc_add_potential_uks
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

   !> The empirical scale on the GWH off-diagonal. 1.75 is the usual value and
   !> the one cuEST uses; the two backends have to start from the same matrix or
   !> comparing their iteration counts means nothing.
   real(dp), parameter :: GWH_K = 1.75_dp

   !> Buffer length for a formatted table line handed to the logger.
   integer, parameter :: LINE_LEN = 160

   public :: rhf_result_t
   public :: run_libcint_rhf
   public :: run_libcint_uhf
   public :: guess_fock                !! Starting Fock for the guesses needing no atomic SCF
   public :: build_fock                !! F = H + J - K/2; with H zero it is the response operator
   public :: density_pseudo_orbitals   !! Factor a guess density for the fitted exchange

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
   end type rhf_result_t

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

contains

   subroutine scf_table_header(verbose)
      !! Column headings for the per-iteration table
      !!
      !! Through the logger like everything else the program says: the level
      !! decides whether it is seen, and a run redirected to a log file gets the
      !! table with it rather than only on a terminal that is no longer there.
      logical, intent(in) :: verbose

      if (.not. verbose) return
      call logger%info("")
      call logger%info("  SCF iterations")
      call logger%info("  "//repeat("-", 84))
      call logger%info("    iter                 energy          dE          dD"// &
                       "   diis       Fock       rest")
      call logger%info("  "//repeat("-", 84))
   end subroutine scf_table_header

   subroutine scf_table_row(verbose, iter, energy, de, drms, ndiis, t_fock, t_rest)
      !! One iteration's line, with the time that iteration took
      !!
      !! Per-iteration rather than only a total because the first iterations of a
      !! direct SCF are not the same price as the last: screening tightens as the
      !! density settles, so a run whose Fock build is flat across iterations is
      !! telling you the density-weighted screening of `INTEGRALS_MQC.md` §6.1 is
      !! missing, which is exactly the case here.
      logical, intent(in) :: verbose
      integer, intent(in) :: iter, ndiis
      real(dp), intent(in) :: energy, de, drms, t_fock, t_rest

      character(len=LINE_LEN) :: line

      if (.not. verbose) return
      write (line, "(i8,f23.12,2es12.3,i7,2(f9.2,a))") &
         iter, energy, de, drms, ndiis, t_fock, " s", t_rest, " s"
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

      character(len=LINE_LEN) :: line

      if (.not. verbose) return
      call logger%info("  "//repeat("-", 84))
      if (converged) then
         write (line, "(a,i0,a)") "  converged in ", iterations, " iterations"
      else
         write (line, "(a,i0,a)") "  NOT converged after ", iterations, " iterations"
      end if
      call logger%info(trim(line))
   end subroutine scf_table_footer

   subroutine run_libcint_rhf(mol, nelec, max_iter, energy_tol, density_tol, &
                              verbose, result, error, aux, diis_vectors, in_core, &
                              guess, guess_density, xc, h_extra)
      !! Drive a closed-shell SCF to convergence
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: nelec
      integer, intent(in) :: max_iter
      real(dp), intent(in) :: energy_tol, density_tol
      logical, intent(in) :: verbose
      type(rhf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
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
         !! One of SCF_GUESS_*. Defaults to the core guess, which is the only
         !! one needing nothing but H -- the caller decides policy, because
         !! `guess_density` has to be built before this routine is entered.
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

      integer :: diis_size, guess_kind
      logical :: use_in_core
      real(dp), allocatable :: bounds(:, :)
      type(direct_stats_t) :: stats

      real(dp), allocatable :: s(:, :), h(:, :), eri(:, :, :, :), bmat(:, :)
      real(dp), allocatable :: x(:, :), fock(:, :), density(:, :), density_old(:, :)
      real(dp), allocatable :: coeff(:, :), eigenvalues(:)
      real(dp), allocatable :: err(:, :), fock_flat(:)
      type(diis_state_t) :: diis
      logical :: extrapolated
      real(dp) :: e_elec, e_old, de, drms
      integer :: n_ao, n_mo, n_occ, iter
      type(timing_report_t) :: clk
      real(dp) :: t_fock_iter, t_rest_iter

      if (mod(nelec, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "RHF needs an even electron count; this "// &
                        "system has an odd one and wants an unrestricted method")
         return
      end if

      diis_size = 8
      if (present(diis_vectors)) diis_size = diis_vectors
      use_in_core = .false.
      if (present(in_core)) use_in_core = in_core
      guess_kind = SCF_GUESS_CORE
      if (present(guess)) guess_kind = guess

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
      else if (use_in_core) then
         call mol%eris(eri)
      else
         ! The bounds depend on the basis and the geometry, not the density, so
         ! one set serves every iteration of the SCF.
         call schwarz_bounds(mol, bounds, error)
         if (error%has_error()) return
      end if

      call build_orthogonalizer(s, x, n_mo, error)
      if (error%has_error()) return
      if (n_occ > n_mo) then
         call error%set(ERROR_VALIDATION, "RHF: more occupied orbitals than the basis "// &
                        "supports after near-null modes were dropped")
         return
      end if

      allocate (fock(n_ao, n_ao), density(n_ao, n_ao), density_old(n_ao, n_ao))
      allocate (err(n_mo, n_mo), fock_flat(n_ao*n_ao))

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
      case (SCF_GUESS_SAC, SCF_GUESS_SAD)
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
      call build_density(coeff, n_occ, density)

      e_old = 0.0_dp
      result%converged = .false.

      ! Everything from the top of the routine to here is one-time cost: the 1e
      ! integrals, the screening bounds or the fitted/in-core tensor, the
      ! orthogonaliser and the guess.
      call clk%lap(STAGE_SETUP)
      call scf_table_header(verbose)

      do iter = 1, max_iter
         density_old = density
         ! The energy belongs to the Fock built from this density, so both come
         ! back together and before extrapolation. A DIIS-mixed Fock is a
         ! convergence device, not a state anything is the energy of.
         call assemble_fock(mol, h, density, coeff, n_occ, bmat, eri, bounds, xc, &
                            fock, e_elec, error, clk=clk)
         if (error%has_error()) return
         t_fock_iter = clk%seconds_of(STAGE_FOCK)
         call clk%lap(STAGE_FOCK)
         t_fock_iter = clk%seconds_of(STAGE_FOCK) - t_fock_iter

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

         call diagonalize(fock, x, n_ao, n_mo, coeff, eigenvalues, error)
         if (error%has_error()) return
         call build_density(coeff, n_occ, density)
         t_rest_iter = t_rest_iter - clk%seconds_of(STAGE_DIAG)
         call clk%lap(STAGE_DIAG)
         t_rest_iter = t_rest_iter + clk%seconds_of(STAGE_DIAG)

         de = abs(e_elec - e_old)
         drms = sqrt(sum((density - density_old)**2)/real(n_ao*n_ao, dp))
         call scf_table_row(verbose, iter, e_elec + mol%nuclear_repulsion(), de, drms, &
                            diis%count(), t_fock_iter, t_rest_iter)

         e_old = e_elec
         result%iterations = iter
         if (iter > 1 .and. de < energy_tol .and. drms < density_tol) then
            result%converged = .true.
            exit
         end if
      end do

      ! The energy that goes out is the one belonging to the density that
      ! satisfied the test, so it is recomputed from the final Fock rather
      ! than carried over from the loop.
      call assemble_fock(mol, h, density, coeff, n_occ, bmat, eri, bounds, xc, &
                         fock, result%electronic, error, clk=clk)
      if (error%has_error()) return
      result%nuclear_repulsion = mol%nuclear_repulsion()
      result%energy = result%electronic + result%nuclear_repulsion
      result%n_occupied = n_occ
      call move_alloc(eigenvalues, result%orbital_energies)
      call move_alloc(coeff, result%orbitals)
      call move_alloc(density, result%density)
      call diis%destroy()

      ! The final rebuild above is a Fock build like any other, so it lands in the
      ! same bucket; `wall` closes over everything including it.
      call clk%lap(STAGE_FOCK)
      call clk%finish()
      call scf_table_footer(verbose, result%converged, result%iterations)
      call clk%report("RHF", verbose)
   end subroutine run_libcint_rhf

   subroutine run_libcint_uhf(mol, nelec, multiplicity, max_iter, energy_tol, density_tol, &
                              verbose, result, error, diis_vectors, in_core, diis_start, &
                              guess, guess_density_alpha, guess_density_beta, xc)
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
      type(rhf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: diis_vectors
      logical, intent(in), optional :: in_core
      integer, intent(in), optional :: diis_start
         !! First iteration allowed to extrapolate. See the default below.
      integer, intent(in), optional :: guess
         !! One of SCF_GUESS_*. Defaults to the core guess.
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

      integer :: diis_size, start_cycle, guess_kind
      logical :: use_in_core
      real(dp), allocatable :: bounds(:, :)
      type(direct_stats_t) :: stats

      real(dp), allocatable :: s(:, :), h(:, :), eri(:, :, :, :)
      real(dp), allocatable :: x(:, :), fock_a(:, :), fock_b(:, :)
      real(dp), allocatable :: d_a(:, :), d_b(:, :), d_a_old(:, :), d_b_old(:, :)
      real(dp), allocatable :: coeff_a(:, :), coeff_b(:, :), eig_a(:), eig_b(:)
      real(dp), allocatable :: err_a(:, :), err_b(:, :), fock_flat(:), err_flat(:)
      type(diis_state_t) :: diis
      logical :: extrapolated
      real(dp) :: e_elec, e_old, de, drms
      integer :: n_ao, n_mo, n_alpha, n_beta, iter, nsq, msq
      type(timing_report_t) :: clk
      real(dp) :: t_fock_iter, t_rest_iter
      character(len=LINE_LEN) :: line

      diis_size = 8
      if (present(diis_vectors)) diis_size = diis_vectors
      start_cycle = DEFAULT_UHF_DIIS_START
      if (present(diis_start)) start_cycle = diis_start
      use_in_core = .false.
      if (present(in_core)) use_in_core = in_core
      guess_kind = SCF_GUESS_CORE
      if (present(guess)) guess_kind = guess

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

      call build_orthogonalizer(s, x, n_mo, error)
      if (error%has_error()) return
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
      case (SCF_GUESS_SAC, SCF_GUESS_SAD)
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

      do iter = 1, max_iter
         d_a_old = d_a
         d_b_old = d_b

         call assemble_fock_uhf(mol, h, d_a, d_b, eri, bounds, xc, fock_a, fock_b, &
                                e_elec, error)
         if (error%has_error()) return
         t_fock_iter = clk%seconds_of(STAGE_FOCK)
         call clk%lap(STAGE_FOCK)
         t_fock_iter = clk%seconds_of(STAGE_FOCK) - t_fock_iter

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
                            diis%count(), t_fock_iter, t_rest_iter)

         e_old = e_elec
         result%iterations = iter
         if (iter > 1 .and. de < energy_tol .and. drms < density_tol) then
            result%converged = .true.
            exit
         end if
      end do

      ! Rebuilt from the density that passed the test, as the restricted path
      ! does, so the reported energy belongs to the reported orbitals.
      call assemble_fock_uhf(mol, h, d_a, d_b, eri, bounds, xc, fock_a, fock_b, &
                             result%electronic, error)
      if (error%has_error()) return
      result%nuclear_repulsion = mol%nuclear_repulsion()
      result%energy = result%electronic + result%nuclear_repulsion
      result%n_occupied = n_alpha
      result%n_occupied_beta = n_beta
      result%spin_squared = spin_contamination(coeff_a, coeff_b, s, n_alpha, n_beta)
      call clk%lap(STAGE_FOCK)
      call clk%finish()
      call scf_table_footer(verbose, result%converged, result%iterations)
      call clk%report("RHF", verbose)
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
                            fock, e_elec, error, clk)
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
      type(xc_context_t), intent(inout), optional :: xc
      real(dp), intent(out) :: fock(:, :)
      real(dp), intent(out) :: e_elec
      type(error_t), intent(inout) :: error
      type(timing_report_t), intent(inout), optional :: clk
         !! Present from the SCF loop, so the exchange-correlation quadrature is
         !! reported apart from the Coulomb/exchange build rather than inside it.
         !! Absent from the guess and gradient callers, which do not report.

      type(direct_stats_t) :: stats
      real(dp), allocatable :: v_xc(:, :)
      real(dp) :: k_scale, e_xc, n_elec
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
         if (xc%range_separated .and. (allocated(bmat) .or. allocated(eri))) then
            call error%set(ERROR_VALIDATION, "a range-separated functional needs the "// &
                           "direct Fock build: the in-core and density-fitted "// &
                           "integrals are built for the full Coulomb kernel, and the "// &
                           "long-range exchange would be missing. Run without "// &
                           "in_core and without an auxiliary basis for the reference.")
            return
         end if
      end if

      if (allocated(bmat)) then
         call build_fock_df(h, bmat, density, coeff, n_occ, fock, k_scale=k_scale)
      else if (allocated(eri)) then
         call build_fock(h, eri, density, fock, k_scale=k_scale)
      else
         call build_fock_direct(mol, h, density, bounds, fock, stats, error, &
                                k_scale=k_scale)
         if (error%has_error()) return
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

      if (kohn_sham) then
         allocate (v_xc(size(h, 1), size(h, 2)))
         call xc_add_potential(xc, mol, density, v_xc, e_xc, n_elec, error)
         if (error%has_error()) return
         fock = fock + v_xc
         e_elec = e_elec + e_xc
         if (present(clk)) call clk%lap(STAGE_XC)
         deallocate (v_xc)
      end if
   end subroutine assemble_fock

   subroutine assemble_fock_uhf(mol, h, d_alpha, d_beta, eri, bounds, xc, &
                                fock_a, fock_b, e_elec, error)
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

      if (kohn_sham) then
         allocate (v_a(n, n), v_b(n, n))
         call xc_add_potential_uks(xc, mol, d_alpha, d_beta, v_a, v_b, e_xc, n_elec, error)
         if (error%has_error()) return
         fock_a = fock_a + v_a
         fock_b = fock_b + v_b
         e_elec = e_elec + e_xc
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

   subroutine build_orthogonalizer(overlap, transform, n_mo, error)
      !! Canonical orthogonaliser X = U s^(-1/2), near-null modes dropped
      !!
      !! Canonical rather than symmetric so a basis with near-linear
      !! dependence loses the offending combinations instead of amplifying
      !! them, which is the same choice the cuEST path makes.
      real(dp), intent(in) :: overlap(:, :)
      real(dp), allocatable, intent(out) :: transform(:, :)
      integer, intent(out) :: n_mo
      type(error_t), intent(inout) :: error

      real(dp), parameter :: NULL_THRESHOLD = 1.0e-7_dp
      real(dp), allocatable :: vectors(:, :), values(:)
      integer :: n_ao, i, kept, info

      n_ao = size(overlap, 1)
      allocate (vectors(n_ao, n_ao), values(n_ao))
      vectors = overlap
      call pic_syev(vectors, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "RHF: overlap diagonalisation failed")
         return
      end if

      kept = count(values > NULL_THRESHOLD)
      if (kept == 0) then
         call error%set(ERROR_VALIDATION, "RHF: the overlap matrix is singular")
         return
      end if

      allocate (transform(n_ao, kept))
      n_mo = 0
      do i = 1, n_ao
         if (values(i) <= NULL_THRESHOLD) cycle
         n_mo = n_mo + 1
         transform(:, n_mo) = vectors(:, i)/sqrt(values(i))
      end do
   end subroutine build_orthogonalizer

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

   pure subroutine build_density(coeff, n_occ, density)
      !! D = 2 C_occ C_occ^T, the closed-shell density
      real(dp), intent(in) :: coeff(:, :)
      integer, intent(in) :: n_occ
      real(dp), intent(out) :: density(:, :)

      integer :: mu, nu, i

      density = 0.0_dp
      do nu = 1, size(density, 2)
         do mu = 1, size(density, 1)
            do i = 1, n_occ
               density(mu, nu) = density(mu, nu) + 2.0_dp*coeff(mu, i)*coeff(nu, i)
            end do
         end do
      end do
   end subroutine build_density

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
      real(dp), intent(in) :: h(:, :), eri(:, :, :, :), density(:, :)
      real(dp), intent(out) :: fock(:, :)
      real(dp), intent(in), optional :: k_scale   !! Exact-exchange fraction, default one

      integer :: mu, nu, la, si, n
      real(dp) :: kf, acc

      n = size(h, 1)
      kf = 0.5_dp
      if (present(k_scale)) kf = 0.5_dp*k_scale
      ! The `if` keeps genuinely tiny systems serial. The response solver enters
      ! this region tens of thousands of times, so for a handful of basis functions
      ! the barriers cost more than the arithmetic; by a few dozen they do not.
      !$omp parallel do collapse(2) default(none) private(mu, nu, la, si, acc) &
      !$omp shared(n, h, eri, density, fock, kf) schedule(static) if(n >= 16)
      do nu = 1, n
         do mu = 1, n
            acc = 0.0_dp
            do si = 1, n
               do la = 1, n
                  ! (mu nu | la si) is Coulomb, (mu la | nu si) is exchange,
                  ! and the half is the closed-shell factor -- D already
                  ! carries the two electrons per orbital.
                  acc = acc + density(la, si) &
                        *(eri(mu, nu, la, si) - kf*eri(mu, la, nu, si))
               end do
            end do
            fock(mu, nu) = h(mu, nu) + acc
         end do
      end do
      !$omp end parallel do
   end subroutine build_fock

   subroutine build_fock_df(h, b, density, coeff, n_occ, fock, k_scale)
      !! F = H + J - K/2 from the fitted tensor rather than the exact ERIs
      !!
      !! Neither term ever forms a four-index object. Density fitting is not the
      !! opposite of a direct build -- it removes the four-index integrals
      !! altogether, so the direct/stored distinction does not apply to it. What
      !! is stored is B, which is n^2 by n_aux: 40 MB where the exact tensor is
      !! 800 MB, and n^3 rather than n^4.
      !!
      !! J is two contractions with the density: c_P = sum_uv B(uv,P) D_uv, then
      !! J_uv = sum_P B(uv,P) c_P.
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
      real(dp), intent(in) :: h(:, :), b(:, :), density(:, :)
      real(dp), intent(in) :: coeff(:, :)   !! MO coefficients; only the occupied block is read
      integer, intent(in) :: n_occ
      real(dp), intent(out) :: fock(:, :)
      real(dp), intent(in), optional :: k_scale   !! Exact-exchange fraction, default one

      real(dp) :: kf
      real(dp), allocatable :: c(:), j(:, :), k(:, :), w(:, :), c_occ(:, :)
      integer :: n, naux, p

      n = size(h, 1)
      naux = size(b, 2)
      allocate (c(naux), j(n, n), k(n, n), w(n, n_occ), c_occ(n, n_occ))

      c_occ = coeff(:, 1:n_occ)

      do p = 1, naux
         c(p) = sum(reshape(b(:, p), [n, n])*density)
      end do

      j = 0.0_dp
      do p = 1, naux
         j = j + c(p)*reshape(b(:, p), [n, n])
      end do

      k = 0.0_dp
      do p = 1, naux
         ! `associate` gives a view of the column rather than a copy: b(:,p) is
         ! contiguous, so reshaping it is free as long as nothing forces a
         ! temporary, and going through gemm keeps it out of matmul's hands.
         associate (b_p => reshape(b(:, p), [n, n]))
            w = 0.0_dp
            call pic_gemm(b_p, c_occ, w)
            call pic_gemm(w, w, k, transb="T", alpha=2.0_dp, beta=1.0_dp)
         end associate
      end do

      kf = 0.5_dp
      if (present(k_scale)) kf = 0.5_dp*k_scale
      fock = h + j - kf*k
   end subroutine build_fock_df

   pure subroutine build_density_spin(coeff, n_occ, density)
      !! D_sigma = C_occ C_occ^T, one electron per occupied orbital
      !!
      !! Deliberately not `build_density` with a factor: the closed-shell one
      !! carries the two electrons per spatial orbital, and reusing it here
      !! would double every spin density. That error converges.
      real(dp), intent(in) :: coeff(:, :)
      integer, intent(in) :: n_occ
      real(dp), intent(out) :: density(:, :)

      integer :: mu, nu, i

      density = 0.0_dp
      do nu = 1, size(density, 2)
         do mu = 1, size(density, 1)
            do i = 1, n_occ
               density(mu, nu) = density(mu, nu) + coeff(mu, i)*coeff(nu, i)
            end do
         end do
      end do
   end subroutine build_density_spin

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

   function spin_contamination(coeff_a, coeff_b, overlap, n_a, n_b) result(s2)
      !! <S^2> = S_z(S_z+1) + n_beta - sum_ij |<a_i|b_j>|^2
      !!
      !! Worth reporting rather than assuming: a UHF solution that has collapsed
      !! onto the restricted one, or onto a badly broken-symmetry state, says so
      !! here and nowhere else. The exact value for a pure doublet is 0.75.
      real(dp), intent(in) :: coeff_a(:, :), coeff_b(:, :), overlap(:, :)
      integer, intent(in) :: n_a, n_b
      real(dp) :: s2

      real(dp) :: sz, overlap_sum
      integer :: i, j, n
      real(dp), allocatable :: sc(:, :), ovl(:, :)

      n = size(overlap, 1)
      sz = 0.5_dp*real(n_a - n_b, dp)
      overlap_sum = 0.0_dp
      if (n_b > 0 .and. n_a > 0) then
         ! S C_beta, then C_alpha^T (S C_beta): the alpha-beta MO overlap block.
         allocate (sc(n, n_b), ovl(n_a, n_b))
         call pic_gemm(overlap, coeff_b(:, 1:n_b), sc)
         call pic_gemm(coeff_a(:, 1:n_a), sc, ovl, transa="T")
         do j = 1, n_b
            do i = 1, n_a
               overlap_sum = overlap_sum + ovl(i, j)**2
            end do
         end do
         deallocate (sc, ovl)
      end if
      s2 = sz*(sz + 1.0_dp) + real(n_b, dp) - overlap_sum
   end function spin_contamination

   pure function electronic_energy(h, fock, density) result(energy)
      !! E = 1/2 sum_uv D_uv (H_uv + F_uv)
      real(dp), intent(in) :: h(:, :), fock(:, :), density(:, :)
      real(dp) :: energy

      energy = 0.5_dp*sum(density*(h + fock))
   end function electronic_energy

end module mqc_libcint_rhf
