!! Where a molecule reacts: Fukui indices by finite difference in the electron count
module mqc_libcint_fukui
   !! The Fukui function is the derivative of the density with respect to the
   !! number of electrons at fixed nuclei,
   !!
   !!     f(r) = [ d rho(r) / d N ]_v
   !!
   !! and it has three values rather than one, because `rho` as a function of
   !! `N` has a kink at every integer. Adding an electron and removing one are
   !! different questions, so `f+` says where the molecule accepts charge and
   !! `f-` where it gives charge up, with `f0` the average for a radical.
   !!
   !! Condensed onto atoms, which is the form anyone uses:
   !!
   !!     f_A+ = q_A(N) - q_A(N+1)      f_A- = q_A(N-1) - q_A(N)
   !!
   !! **Computed by difference rather than from frontier orbitals**, so the
   !! relaxation of every other orbital in response to the added charge is
   !! kept. It costs three SCFs instead of one.
   !!
   !! The three calculations share one `libcint_molecule_t`, so they see the
   !! same geometry and the same basis functions by construction.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_libcint_pcm, only: pcm_context_t
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_uhf, &
                              SCF_GUESS_GWH, SCF_GUESS_SAC, SCF_GUESS_SAD
   use mqc_libcint_atomic_guess, only: build_atomic_guess
   use mqc_scf_common, only: build_density_spin
   use mqc_libcint_charges, only: mulliken_charges, chelpg_charges
   use mqc_libcint_xc, only: xc_context_t, xc_context_create, xc_available
   use mqc_libcint_mp2, only: mp2_result_t, run_libcint_mp2, run_libcint_ri_mp2, &
                              run_libcint_ump2, run_libcint_uri_mp2
   implicit none
   private

   public :: fukui_result_t
   public :: fukui_indices
   public :: print_fukui_report

   real(dp), parameter :: TOL_SUM_RULE = 1.0e-6_dp
      !! How far the condensed indices may drift from summing to one. Looser
      !! than the identities elsewhere in this code because CHELPG solves a
      !! constrained least-squares problem on a grid rather than partitioning a
      !! density exactly, so its total is right to the conditioning of that
      !! solve and not to machine precision.

   type :: fukui_result_t
      !! What the three calculations leave behind
      real(dp), allocatable :: f_plus(:)    !! Nucleophilic attack: where charge arrives
      real(dp), allocatable :: f_minus(:)   !! Electrophilic attack: where charge leaves
      real(dp), allocatable :: f_zero(:)    !! Radical attack, the average of the two
         !! **A NEGATIVE ENTRY IN ANY OF THESE IS SPURIOUS.** The exact Fukui
         !! function cannot be negative; negative condensed values are an
         !! artefact of partitioning a continuous density onto atoms. Rank
         !! sites by these rather than quoting the numbers, and read a small
         !! negative as "not reactive in this channel".
      real(dp), allocatable :: dual(:)
         !! `f+ - f-`. Positive where a site prefers to accept and negative
         !! where it prefers to donate.
      real(dp), allocatable :: q_neutral(:), q_anion(:), q_cation(:)
      real(dp) :: energy_neutral = 0.0_dp
      real(dp) :: energy_anion = 0.0_dp
      real(dp) :: energy_cation = 0.0_dp
      real(dp) :: ionisation_potential = 0.0_dp   !! E(N-1) - E(N)
      real(dp) :: electron_affinity = 0.0_dp      !! E(N) - E(N+1)
      real(dp) :: chemical_potential = 0.0_dp     !! -(IP + EA)/2
      real(dp) :: hardness = 0.0_dp               !! (IP - EA)/2
      real(dp) :: electrophilicity = 0.0_dp       !! mu^2 / 2 eta
      integer :: iterations_anion = 0
      integer :: iterations_cation = 0
         !! What each ion SCF cost, the ions being the two states that
         !! misbehave.
      logical :: anion_bound = .true.
         !! Whether `E(N+1)` came out below `E(N)`. When it did not, nothing
         !! bound the extra electron and `f+` describes whatever orbital was
         !! left over. **Two different things produce this and they cannot be
         !! told apart from here**: a basis with no diffuse functions, and a
         !! molecule with no bound anion.
      character(len=16) :: scheme = "chelpg"
      character(len=:), allocatable :: functional
         !! The functional the three states were computed with, empty for
         !! Hartree-Fock. IP, EA and hardness are differences of total energies
         !! and are not comparable across methods.
   end type fukui_result_t

contains

   subroutine fukui_indices(mol, nelec, multiplicity, neutral_density, neutral_energy, &
                            scheme, max_iter, energy_tol, density_tol, res, error, &
                            functional, grid_level, pt2_fraction, neutral_orbitals, &
                            neutral_orbital_energies, aux, n_frozen, verbose, pcm, &
                            level_shift, diis_vectors, nlc_grid_level, &
                            screening_tolerance, block_size, guess_mode, incremental_fock, &
                            accelerator, grad_tol, linear_dependence, unseeded_guess, &
                            allow_unconverged)
      !! Run the ions and condense the difference onto atoms
      !!
      !! The neutral density is handed in rather than recomputed, since the
      !! caller has just converged it. The molecule must be a closed-shell
      !! neutral, so that both ions are doublets.
      !!
      !! **Kohn-Sham ions need their own context**, which is why this takes a
      !! functional NAME and not the caller's `xc_context_t`: libxc fixes the
      !! spin channel when a functional is initialised, and the caller's
      !! context was built spin-unpolarised for the neutral. One polarised
      !! context is built here and shared by both doublets. Absent or empty
      !! `functional` is Hartree-Fock.
      !!
      !! **A double hybrid's perturbative term is evaluated here for all three
      !! states**, the neutral included, even though the caller computes the
      !! neutral's separately for the energy it reports. IP and EA are
      !! differences of total energies, so the term has to be present in all
      !! three or in none, on one code path's frozen core and auxiliary basis.
      !!
      !! `aux` present fits the correlation, absent computes it exactly.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: nelec
      integer, intent(in) :: multiplicity
      real(dp), intent(in) :: neutral_density(:, :)   !! Total density, both spins
      real(dp), intent(in) :: neutral_energy
      integer, intent(in), optional :: diis_vectors
         !! The DIIS subspace, applied to both ions.
      real(dp), intent(in), optional :: level_shift
         !! Applied to both ions, which are the open-shell SCFs a shift is
         !! asked for to stabilise.
      character(len=*), intent(in) :: scheme          !! "chelpg" or "mulliken"
      integer, intent(in) :: max_iter
      real(dp), intent(in) :: energy_tol, density_tol
      type(fukui_result_t), intent(out) :: res
      type(error_t), intent(inout) :: error
      character(len=*), intent(in), optional :: functional
      integer, intent(in), optional :: grid_level
      logical, intent(in), optional :: incremental_fock
         !! Passed to both ion SCFs, so that turning incremental building off
         !! turns it off for the whole Fukui run.
      integer, intent(in), optional :: accelerator
         !! `properties.fukui.scf.accelerator`, already parsed to an `ACCEL_*`
         !! kind by the caller. Absent leaves the SCF's own default.
      real(dp), intent(in), optional :: grad_tol
         !! `properties.fukui.scf.gradient_tolerance`. Absent derives it from
         !! the energy tolerance, as the neutral's SCF does.
      real(dp), intent(in), optional :: linear_dependence
         !! `properties.fukui.scf.linear_dependence_threshold`. A diffuse basis
         !! chosen so the anion binds its extra electron is also the one whose
         !! overlap goes near-singular.
      integer, intent(in), optional :: unseeded_guess
         !! An `SCF_GUESS_*` kind for an ion that is NOT seeded from the
         !! neutral -- so it bites under `guess_mode = "independent"`, and is
         !! ignored otherwise, because a seeded ion already has a density and
         !! nothing left for a guess to decide.
         !!
         !! Absent keeps the GWH fallback below. The caller passes this only
         !! for an explicit spelling and not for `auto`, so a deck that never
         !! mentioned a guess is not quietly moved off GWH.
      logical, intent(in), optional :: allow_unconverged
         !! `properties.fukui.scf.allow_crap_scf`. Keep an ion that did not
         !! converge instead of refusing to form the index, with a warning.
         !! Every number downstream is then a result to follow up rather than
         !! to trust.
      character(len=*), intent(in), optional :: guess_mode
         !! "neutral" (the default) or "independent".
         !!
         !! "neutral" seeds each ion from the converged neutral's orbitals at
         !! that ion's own occupation, which is usually much the better start.
         !! But it fills the *neutral's* LUMO, and where that is a poor guess
         !! for the ion's own the SCF can settle on a different stationary
         !! point -- an anion stalling a few millihartree above the answer the
         !! ordinary guess reaches, which DIIS will not rescue.
         !! "independent" gives each ion the ordinary guess instead.
      integer, intent(in), optional :: nlc_grid_level
         !! VV10's own quadrature level, forwarded from
         !! `keywords.dft.nlc_grid_level`. Absent leaves `NLC_GRID_LEVEL`, which
         !! would integrate the ions' non-local term on a different grid from
         !! the neutral's.
      real(dp), intent(in), optional :: screening_tolerance
         !! AO screening threshold for the quadrature, from
         !! `keywords.dft.screening_tolerance`. The ions have no business
         !! screening differently from the neutral.
      integer, intent(in), optional :: block_size
         !! Grid points per block, from `keywords.dft.block_size`.
      real(dp), intent(in), optional :: pt2_fraction
      real(dp), intent(in), optional :: neutral_orbitals(:, :)
      real(dp), intent(in), optional :: neutral_orbital_energies(:)
      type(libcint_molecule_t), intent(in), optional :: aux
      type(pcm_context_t), intent(inout), optional :: pcm
         !! The continuum the neutral was converged in, handed on to both ions.
         !!
         !! **All three states or none.** An ion solved in vacuum and
         !! differenced against a neutral solved in solvent is not an
         !! ionisation potential of anything -- the solvation energy of a
         !! charged species is the size of the quantity being reported. The
         !! same context object is reused, the cavity being a property of a
         !! geometry that does not move between the three states.
         !!
         !! Equilibrium solvation for all three, so the indices and potentials
         !! here are the *adiabatic* ones. A vertical quantity would need the
         !! slow polarisation frozen at the neutral's value, which needs an
         !! optical dielectric this program does not carry.
      logical, intent(in), optional :: verbose
         !! Print the two ion SCFs the way the neutral is printed. Absent is
         !! quiet.
      integer, intent(in), optional :: n_frozen

      type(rhf_result_t) :: anion, cation
      real(dp), allocatable :: overlap(:, :), total(:, :)
      real(dp), allocatable :: d_guess_a(:, :), d_guess_b(:, :)
      integer :: n_ao, n_occ_neutral, guess_kind
      type(error_t) :: guess_error
      logical :: seeded, keep_unconverged
      logical :: loud
      integer :: natm
      type(xc_context_t), target :: xc_ions
      type(xc_context_t), pointer :: xc_arg
      logical :: kohn_sham, double_hybrid
      real(dp) :: pt2, e_pt2_neutral, e_pt2_anion, e_pt2_cation
      integer :: frozen

      if (error%has_error()) return

      loud = .false.
      if (present(verbose)) loud = verbose

      if (multiplicity /= 1) then
         call error%set(ERROR_VALIDATION, "Fukui indices are computed here for a "// &
                        "closed-shell neutral, whose ions are both doublets. This "// &
                        "molecule has multiplicity "//to_char(multiplicity)//", and "// &
                        "guessing what its ions should be is the caller's decision "// &
                        "rather than this routine's.")
         return
      end if
      if (nelec < 2) then
         call error%set(ERROR_VALIDATION, "a molecule with fewer than two electrons "// &
                        "has no closed shell to remove one from.")
         return
      end if

      kohn_sham = .false.
      if (present(functional)) kohn_sham = len_trim(functional) > 0
      pt2 = 0.0_dp
      if (present(pt2_fraction)) pt2 = pt2_fraction
      double_hybrid = kohn_sham .and. pt2 /= 0.0_dp
      frozen = 0
      if (present(n_frozen)) frozen = n_frozen
      e_pt2_neutral = 0.0_dp
      e_pt2_anion = 0.0_dp
      e_pt2_cation = 0.0_dp
      res%functional = ""
      xc_arg => null()

      if (kohn_sham) then
         if (.not. xc_available()) then
            call error%set(ERROR_VALIDATION, "the Fukui ions were asked for with '"// &
                           trim(functional)//"' but this build has no libxc.")
            return
         end if
         ! Spin-polarised, because both ions are doublets, and given every grid
         ! knob the neutral was given rather than `level` alone. An absent
         ! optional passed on stays absent, so a caller that says nothing still
         ! gets the defaults.
         call xc_context_create(mol, trim(functional), xc_ions, error, &
                                level=grid_level, polarized=.true., &
                                nlc_level=nlc_grid_level, &
                                screen_tol=screening_tolerance, &
                                point_block=block_size)
         if (error%has_error()) then
            call error%add_context("Fukui indices: the ions' functional")
            return
         end if
         ! The orbitals a double hybrid's perturbative term needs are not
         ! derivable from the density this routine was handed, so they are
         ! required rather than assumed.
         if (double_hybrid) then
            if (.not. (present(neutral_orbitals) .and. present(neutral_orbital_energies))) then
               call error%set(ERROR_VALIDATION, "'"//trim(functional)//"' is a double "// &
                              "hybrid, so the Fukui analysis needs the neutral's orbitals "// &
                              "to evaluate its perturbative term on the same footing as "// &
                              "the ions. The caller passed a density only.")
               return
            end if
         end if
         xc_arg => xc_ions
         res%functional = trim(functional)
      end if

      natm = mol%natm
      res%scheme = scheme
      ! The neutral's perturbative term, on a closed shell and so restricted.
      ! Computed before the ions so a failure here is reported before two
      ! unrestricted SCFs have been paid for.
      if (double_hybrid) then
         call state_pt2(mol, aux, neutral_orbitals, neutral_orbital_energies, &
                        nelec/2, frozen, e_pt2_neutral, error)
         if (error%has_error()) then
            call error%add_context("Fukui indices: the neutral's PT2 term")
            return
         end if
         e_pt2_neutral = pt2*e_pt2_neutral
      end if
      res%energy_neutral = neutral_energy + e_pt2_neutral

      call mol%overlap(overlap)
      call condense(mol, neutral_density, overlap, scheme, 0.0_dp, res%q_neutral, error)
      if (error%has_error()) return

      ! The anion. One more electron, and a doublet because the neutral was
      ! closed shell.
      !
      ! Each ion is seeded at *its own* occupation rather than from the neutral
      ! density halved: filling the neutral's orbitals Aufbau for a doublet
      ! gives the anion (n/2 + 1, n/2) and the cation (n/2, n/2 - 1), the right
      ! electron count and the right spin, where a halved density has neither.
      !
      ! `d_guess_a` and `d_guess_b` stay unallocated only when neither the
      ! neutral's orbitals nor the free atoms can supply a density, and an
      ! unallocated allocatable actual argument is *absent* at an optional
      ! non-allocatable dummy (F2018 15.5.2.13). That is what makes one call
      ! site serve every path: on that last fallback `guess_kind` becomes GWH
      ! and the two densities vanish from the call rather than erroring.
      keep_unconverged = .false.
      if (present(allow_unconverged)) keep_unconverged = allow_unconverged
      n_ao = size(neutral_density, 1)
      n_occ_neutral = nelec/2
      guess_kind = SCF_GUESS_GWH
      if (present(unseeded_guess)) guess_kind = unseeded_guess
      seeded = present(neutral_orbitals)
      ! An explicit "independent" turns the seeding off and nothing else: the
      ! orbitals are still wanted for a double hybrid's perturbative term. Only
      ! the two guess densities go unbuilt.
      if (seeded .and. present(guess_mode)) then
         seeded = .not. (trim(guess_mode) == "independent")
      end if
      if (seeded) seeded = size(neutral_orbitals, 1) == n_ao .and. &
                           size(neutral_orbitals, 2) >= n_occ_neutral + 1
      if (seeded) then
         guess_kind = SCF_GUESS_SAD
         allocate (d_guess_a(n_ao, n_ao), d_guess_b(n_ao, n_ao))
      else if (guess_kind == SCF_GUESS_SAD .or. guess_kind == SCF_GUESS_SAC) then
         ! Only when the deck asked for one. `SCF_GUESS_SAD` is the kind that
         ! says "the starting density is supplied", and in the seeded branch
         ! above it is the neutral's orbitals supplying it; asked for here with
         ! nothing to seed from, it has to come from the free atoms. The
         ! default stays GWH, so a deck that did not ask pays for no atomic
         ! solves.
         call build_atomic_guess(mol, guess_kind, d_guess_a, d_guess_b, guess_error)
         if (guess_error%has_error()) then
            ! A free atom that will not converge is a reason to start
            ! somewhere else, not to fail the ion. Loud, because an open-shell
            ! ion can land on a different state from a different start.
            call logger%warning("Fukui: atomic guess for the ions: "// &
                                guess_error%get_message()//" -- falling back to gwh")
            guess_kind = SCF_GUESS_GWH
            if (allocated(d_guess_a)) deallocate (d_guess_a)
            if (allocated(d_guess_b)) deallocate (d_guess_b)
         end if
      end if

      if (loud) then
         call logger%info("")
         call logger%info("  Fukui: anion SCF, N+1 electrons")
      end if
      if (seeded) then
         call build_density_spin(neutral_orbitals, n_occ_neutral + 1, d_guess_a)
         call build_density_spin(neutral_orbitals, n_occ_neutral, d_guess_b)
      end if
      call run_libcint_uhf(mol, nelec + 1, 2, max_iter, energy_tol, density_tol, &
                           loud, anion, error, xc=xc_arg, pcm=pcm, &
                           diis_vectors=diis_vectors, &
                           level_shift=level_shift, &
                           guess=guess_kind, guess_density_alpha=d_guess_a, &
                           guess_density_beta=d_guess_b, &
                           incremental_fock=incremental_fock, &
                           accelerator=accelerator, grad_tol=grad_tol, &
                           linear_dependence=linear_dependence)
      if (error%has_error()) then
         call error%add_context("Fukui indices: the anion")
         return
      end if
      if (.not. anion%converged) then
         if (.not. keep_unconverged) then
            call error%set(ERROR_VALIDATION, "the anion did not converge, so f+ cannot "// &
                           "be formed. Anions converge poorly in a basis with no diffuse "// &
                           "functions, which is usually the cause rather than the SCF.")
            return
         end if
         call logger%warning("  Fukui: the anion did not converge, and "// &
                             "allow_crap_scf kept it. Every index below is "// &
                             "built from it and is a result to follow up.")
      end if
      allocate (total(size(anion%density, 1), size(anion%density, 2)))
      total = anion%density + anion%density_beta
      call condense(mol, total, overlap, scheme, -1.0_dp, res%q_anion, error)
      if (error%has_error()) return
      if (double_hybrid) then
         call state_pt2_open(mol, aux, anion, frozen, e_pt2_anion, error)
         if (error%has_error()) then
            call error%add_context("Fukui indices: the anion's PT2 term")
            return
         end if
         e_pt2_anion = pt2*e_pt2_anion
      end if
      res%energy_anion = anion%energy + e_pt2_anion
      res%iterations_anion = anion%iterations
      deallocate (total)

      ! The cation.
      if (loud) then
         call logger%info("")
         call logger%info("  Fukui: cation SCF, N-1 electrons")
      end if
      if (seeded) then
         call build_density_spin(neutral_orbitals, n_occ_neutral, d_guess_a)
         call build_density_spin(neutral_orbitals, n_occ_neutral - 1, d_guess_b)
      end if
      call run_libcint_uhf(mol, nelec - 1, 2, max_iter, energy_tol, density_tol, &
                           loud, cation, error, xc=xc_arg, pcm=pcm, &
                           diis_vectors=diis_vectors, &
                           level_shift=level_shift, &
                           guess=guess_kind, guess_density_alpha=d_guess_a, &
                           guess_density_beta=d_guess_b, &
                           incremental_fock=incremental_fock, &
                           accelerator=accelerator, grad_tol=grad_tol, &
                           linear_dependence=linear_dependence)
      if (error%has_error()) then
         call error%add_context("Fukui indices: the cation")
         return
      end if
      if (.not. cation%converged) then
         if (.not. keep_unconverged) then
            call error%set(ERROR_VALIDATION, "the cation did not converge, so f- cannot "// &
                           "be formed.")
            return
         end if
         call logger%warning("  Fukui: the cation did not converge, and "// &
                             "allow_crap_scf kept it. Every index below is "// &
                             "built from it and is a result to follow up.")
      end if
      allocate (total(size(cation%density, 1), size(cation%density, 2)))
      total = cation%density + cation%density_beta
      call condense(mol, total, overlap, scheme, 1.0_dp, res%q_cation, error)
      if (error%has_error()) return
      if (double_hybrid) then
         call state_pt2_open(mol, aux, cation, frozen, e_pt2_cation, error)
         if (error%has_error()) then
            call error%add_context("Fukui indices: the cation's PT2 term")
            return
         end if
         e_pt2_cation = pt2*e_pt2_cation
      end if
      res%energy_cation = cation%energy + e_pt2_cation
      res%iterations_cation = cation%iterations
      deallocate (total, overlap)

      allocate (res%f_plus(natm), res%f_minus(natm), res%f_zero(natm), res%dual(natm))
      res%f_plus = res%q_neutral - res%q_anion
      res%f_minus = res%q_cation - res%q_neutral
      res%f_zero = 0.5_dp*(res%f_plus + res%f_minus)
      res%dual = res%f_plus - res%f_minus

      res%ionisation_potential = res%energy_cation - res%energy_neutral
      res%electron_affinity = res%energy_neutral - res%energy_anion
      res%chemical_potential = -0.5_dp*(res%ionisation_potential + res%electron_affinity)
      ! Absolute hardness is Parr and Pearson's `(IP - EA)/2`, the second
      ! derivative of the energy with respect to electron number, as the
      ! chemical potential above is the first. The half matters: the
      ! electrophilicity `omega = mu^2 / 2 eta` is written for the halved
      ! quantity.
      res%hardness = 0.5_dp*(res%ionisation_potential - res%electron_affinity)
      if (abs(res%hardness) > 1.0e-10_dp) then
         res%electrophilicity = res%chemical_potential**2/(2.0_dp*res%hardness)
      end if
      res%anion_bound = res%electron_affinity > 0.0_dp

      call check_sum_rule(res, error)
   end subroutine fukui_indices

   subroutine state_pt2(mol, aux, coeff, eps, n_occ, frozen, e_corr, error)
      !! The UNSCALED MP2 correlation for the closed-shell neutral
      !!
      !! Unscaled: the caller multiplies by the functional's fraction.
      type(libcint_molecule_t), intent(in) :: mol
      type(libcint_molecule_t), intent(in), optional :: aux
      real(dp), intent(in) :: coeff(:, :), eps(:)
      integer, intent(in) :: n_occ, frozen
      real(dp), intent(out) :: e_corr
      type(error_t), intent(inout) :: error

      type(mp2_result_t) :: m
      e_corr = 0.0_dp
      if (present(aux)) then
         call run_libcint_ri_mp2(mol, aux, coeff, eps, n_occ, 0.0_dp, m, error, &
                                 n_frozen=frozen)
      else
         call run_libcint_mp2(mol, coeff, eps, n_occ, 0.0_dp, m, error, n_frozen=frozen)
      end if
      if (error%has_error()) return
      e_corr = m%correlation
   end subroutine state_pt2

   subroutine state_pt2_open(mol, aux, scf, frozen, e_corr, error)
      !! The same for one of the doublet ions
      !!
      !! Takes the whole SCF result rather than its pieces, so that the
      !! per-spin occupied counts come from the same place the orbitals did.
      type(libcint_molecule_t), intent(in) :: mol
      type(libcint_molecule_t), intent(in), optional :: aux
      type(rhf_result_t), intent(in) :: scf
      integer, intent(in) :: frozen
      real(dp), intent(out) :: e_corr
      type(error_t), intent(inout) :: error

      type(mp2_result_t) :: m
      e_corr = 0.0_dp
      if (present(aux)) then
         call run_libcint_uri_mp2(mol, aux, scf%orbitals, scf%orbitals_beta, &
                                  scf%orbital_energies, scf%orbital_energies_beta, &
                                  scf%n_occupied, scf%n_occupied_beta, 0.0_dp, m, error, &
                                  n_frozen=frozen)
      else
         call run_libcint_ump2(mol, scf%orbitals, scf%orbitals_beta, &
                               scf%orbital_energies, scf%orbital_energies_beta, &
                               scf%n_occupied, scf%n_occupied_beta, 0.0_dp, m, error, &
                               n_frozen=frozen)
      end if
      if (error%has_error()) return
      e_corr = m%correlation
   end subroutine state_pt2_open

   subroutine condense(mol, density, overlap, scheme, total_charge, charges, error)
      !! Atomic charges for one of the three states
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(in) :: overlap(:, :)
      character(len=*), intent(in) :: scheme
      real(dp), intent(in) :: total_charge
      real(dp), allocatable, intent(out) :: charges(:)
      type(error_t), intent(inout) :: error

      if (error%has_error()) return

      select case (trim(scheme))
      case ("chelpg")
         ! The fit is told the ion's total charge, so the difference between
         ! two states is a redistribution of exactly one electron rather than
         ! of one electron plus whatever the two fits disagreed about.
         call chelpg_charges(mol, density, charges, error, total_charge=total_charge)
      case ("mulliken")
         call mulliken_charges(mol, density, overlap, charges, error)
      case default
         call error%set(ERROR_VALIDATION, "unknown population scheme '"//trim(scheme)// &
                        "' for Fukui indices; expected 'chelpg' or 'mulliken'.")
      end select
   end subroutine condense

   subroutine check_sum_rule(res, error)
      !! Both indices must sum to one over the molecule
      !!
      !! The three states differ by exactly one electron and the charges account
      !! for all of it, so `sum_A f_A+` and `sum_A f_A-` are one by
      !! construction. What this catches is a density that was not the total,
      !! an ion run at the wrong charge, or a population scheme that lost part
      !! of the molecule.
      type(fukui_result_t), intent(in) :: res
      type(error_t), intent(inout) :: error

      character(len=300) :: message
      real(dp) :: plus, minus

      if (error%has_error()) return

      plus = sum(res%f_plus)
      minus = sum(res%f_minus)
      if (abs(plus - 1.0_dp) > TOL_SUM_RULE .or. abs(minus - 1.0_dp) > TOL_SUM_RULE) then
         write (message, "(a,f12.8,a,f12.8,a)") "the condensed Fukui indices do not "// &
            "sum to one: f+ gives ", plus, " and f- gives ", minus, &
            ". One electron was added and one removed, so the difference has to "// &
            "account for exactly that much charge."
         call error%set(ERROR_VALIDATION, trim(message))
      end if
   end subroutine check_sum_rule

   subroutine print_fukui_report(res, element_symbols)
      !! The per-atom table and the global descriptors
      type(fukui_result_t), intent(in) :: res
      character(len=*), intent(in) :: element_symbols(:)

      character(len=160) :: line
      character(len=16) :: label
      integer :: natm, a
      logical :: any_negative

      natm = size(res%f_plus)
      call logger%info("")
      ! The method belongs in the header for the same reason the scheme does:
      ! IP, EA and hardness are differences of total energies, so they are not
      ! comparable across methods.
      if (allocated(res%functional)) then
         if (len_trim(res%functional) > 0) then
            call logger%info("  Fukui ions: unrestricted Kohn-Sham, "// &
                             trim(res%functional))
         else
            call logger%info("  Fukui ions: unrestricted Hartree-Fock")
         end if
      end if
      call logger%info("  ==== Fukui indices ("//trim(res%scheme)//") "// &
                       "=================================")
      call logger%info("")
      call logger%info("     atom            f+          f-          f0        dual")
      any_negative = .false.
      do a = 1, natm
         write (label, "(a,i0)") trim(adjustl(element_symbols(a)))//" ", a
         write (line, "(4x,a8,4f12.4)") label, res%f_plus(a), res%f_minus(a), &
            res%f_zero(a), res%dual(a)
         call logger%info(trim(line))
         if (res%f_plus(a) < 0.0_dp .or. res%f_minus(a) < 0.0_dp) any_negative = .true.
      end do
      write (line, "(4x,a8,4f12.4)") "sum", sum(res%f_plus), sum(res%f_minus), &
         sum(res%f_zero), sum(res%dual)
      call logger%info(trim(line))

      call logger%info("")
      call logger%info("     f+ is where the molecule accepts charge, f- where it "// &
                       "gives it up;")
      call logger%info("     dual is f+ minus f-, positive at electrophilic sites "// &
                       "and negative at")
      call logger%info("     nucleophilic ones.")

      call logger%info("")
      write (line, "(4x,a24,f14.6,a)") "ionisation potential", &
         res%ionisation_potential, " hartree"
      call logger%info(trim(line))
      write (line, "(4x,a24,f14.6,a)") "electron affinity", res%electron_affinity, &
         " hartree"
      call logger%info(trim(line))
      write (line, "(4x,a24,f14.6,a)") "chemical potential", res%chemical_potential, &
         " hartree"
      call logger%info(trim(line))
      write (line, "(4x,a24,f14.6,a)") "hardness", res%hardness, " hartree"
      call logger%info(trim(line))
      write (line, "(4x,a24,f14.6,a)") "electrophilicity", res%electrophilicity, &
         " hartree"
      call logger%info(trim(line))

      ! A negative condensed index is not a site that repels electrons; it is
      ! the population analysis failing to divide the density sensibly, which
      ! is common with Mulliken in a basis with diffuse functions.
      if (any_negative) then
         call logger%warning("")
         call logger%warning("  NEGATIVE INDICES ABOVE ARE SPURIOUS -- take them with a "// &
                             "grain of salt.")
         call logger%warning("")
         call logger%warning("  f+ and f- are derivatives of a density with respect to "// &
                             "electron count, so")
         call logger%warning("  the exact quantities cannot be negative: no site "// &
                             "releases charge when the")
         call logger%warning("  molecule gains an electron. A negative value is an "// &
                             "artefact of dividing a")
         call logger%warning("  continuous density among atoms, not a chemical "// &
                             "finding, and it appears")
         call logger%warning("  where two states' charges are fitted independently and "// &
                             "disagree by more")
         call logger%warning("  than the electron being moved.")
         call logger%warning("")
         call logger%warning("  Read the RANKING, not the number: a site with a small "// &
                             "negative index is")
         call logger%warning("  unreactive in that channel, not anti-reactive. Do not "// &
                             "quote the value,")
         call logger%warning("  and do not build a further descriptor on it -- f0 and "// &
                             "dual inherit the")
         call logger%warning("  artefact from whichever index carried it.")
         if (trim(res%scheme) == "mulliken") then
            call logger%warning("")
            call logger%warning("  Mulliken is especially prone to this, being basis-set "// &
                                "dependent; chelpg")
            call logger%warning("  is the default for that reason and is worth rerunning "// &
                                "with before")
            call logger%warning("  concluding anything from this table.")
         else
            call logger%warning("")
            call logger%warning("  A larger basis usually shrinks it. If it survives "// &
                                "that, the atom")
            call logger%warning("  simply carries little of this channel.")
         end if
         call logger%warning("")
      end if

      ! The direct test rather than a guess from the basis-set name: if the
      ! anion sits above the neutral then nothing bound the extra electron.
      if (.not. res%anion_bound) then
         call logger%warning("  the anion lies ABOVE the neutral, so nothing bound "// &
                             "the added electron and f+ describes whatever orbital "// &
                             "was left over")
         call logger%warning("  this is a basis without diffuse functions, or a "// &
                             "molecule with no bound anion, and the two look alike "// &
                             "from here")
      end if
   end subroutine print_fukui_report

end module mqc_libcint_fukui
