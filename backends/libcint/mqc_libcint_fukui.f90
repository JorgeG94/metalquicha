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
   !! **Computed by difference rather than from frontier orbitals.** The cheap
   !! route approximates `f+` by the shape of the LUMO and `f-` by the HOMO,
   !! which costs one SCF instead of three and throws away the relaxation of
   !! every other orbital in response to the added charge. That relaxation is
   !! not a correction; on a polar molecule it is most of the answer, and it is
   !! precisely what distinguishes the sites this analysis is run to rank.
   !!
   !! The three calculations share one `libcint_molecule_t`, so they see the
   !! same geometry and the same basis functions by construction rather than by
   !! agreement.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_uhf
   use mqc_libcint_charges, only: mulliken_charges, chelpg_charges
   use mqc_libcint_xc, only: xc_context_t, xc_context_create, xc_available
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
      real(dp), allocatable :: dual(:)
         !! `f+ - f-`. Positive where a site prefers to accept and negative
         !! where it prefers to donate, so one column separates electrophilic
         !! from nucleophilic sites that both look reactive in either index
         !! alone.
      real(dp), allocatable :: q_neutral(:), q_anion(:), q_cation(:)
      real(dp) :: energy_neutral = 0.0_dp
      real(dp) :: energy_anion = 0.0_dp
      real(dp) :: energy_cation = 0.0_dp
      real(dp) :: ionisation_potential = 0.0_dp   !! E(N-1) - E(N)
      real(dp) :: electron_affinity = 0.0_dp      !! E(N) - E(N+1)
      real(dp) :: chemical_potential = 0.0_dp     !! -(IP + EA)/2
      real(dp) :: hardness = 0.0_dp               !! IP - EA
      real(dp) :: electrophilicity = 0.0_dp       !! mu^2 / 2 eta
      logical :: anion_bound = .true.
         !! Whether `E(N+1)` came out below `E(N)`. When it did not, nothing
         !! bound the extra electron and `f+` describes whatever orbital was
         !! left over. **Two different things produce this and they cannot be
         !! told apart from here**: a basis with no diffuse functions, and a
         !! molecule with no bound anion. Water is the second -- adding
         !! diffuse functions moves its affinity from -0.19 to -0.04 hartree
         !! and never makes it positive, because H2O- does not exist.
      character(len=16) :: scheme = "chelpg"
      character(len=:), allocatable :: functional
         !! The functional the three states were computed with, empty for
         !! Hartree-Fock. Recorded because an IP of 0.4 hartree means nothing
         !! without it: the number is a difference of total energies and those
         !! are not comparable across methods.
   end type fukui_result_t

contains

   subroutine fukui_indices(mol, nelec, multiplicity, neutral_density, neutral_energy, &
                            scheme, max_iter, energy_tol, density_tol, res, error, &
                            functional, grid_level)
      !! Run the ions and condense the difference onto atoms
      !!
      !! The neutral density is handed in rather than recomputed, since the
      !! caller has just converged it.
      !!
      !! KOHN-SHAM IONS NEED THEIR OWN CONTEXT, which is why this takes a
      !! functional NAME and not the caller's `xc_context_t`. libxc fixes the
      !! spin channel when a functional is initialised, and the caller's
      !! context was built for the closed-shell neutral, so it is
      !! spin-unpolarised. Handing it to the unrestricted ions does not
      !! misread it -- `xc_add_potential` refuses it outright -- but the
      !! refusal would arrive as a failed analysis rather than as a design.
      !! One polarised context is built here and shared by both doublets,
      !! which see the same grid as each other by construction.
      !!
      !! Absent or empty `functional` is Hartree-Fock, which is what this
      !! routine did before it could do anything else.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: nelec
      integer, intent(in) :: multiplicity
      real(dp), intent(in) :: neutral_density(:, :)   !! Total density, both spins
      real(dp), intent(in) :: neutral_energy
      character(len=*), intent(in) :: scheme          !! "chelpg" or "mulliken"
      integer, intent(in) :: max_iter
      real(dp), intent(in) :: energy_tol, density_tol
      type(fukui_result_t), intent(out) :: res
      type(error_t), intent(inout) :: error
      character(len=*), intent(in), optional :: functional
      integer, intent(in), optional :: grid_level

      type(rhf_result_t) :: anion, cation
      real(dp), allocatable :: overlap(:, :), total(:, :)
      integer :: natm
      type(xc_context_t), target :: xc_ions
      type(xc_context_t), pointer :: xc_arg
      logical :: kohn_sham

      if (error%has_error()) return

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
      res%functional = ""
      xc_arg => null()

      if (kohn_sham) then
         if (.not. xc_available()) then
            call error%set(ERROR_VALIDATION, "the Fukui ions were asked for with '"// &
                           trim(functional)//"' but this build has no libxc.")
            return
         end if
         ! Spin-polarised, because both ions are doublets. See the note on this
         ! routine: this is the whole reason the functional arrives by name.
         call xc_context_create(mol, trim(functional), xc_ions, error, &
                                level=grid_level, polarized=.true.)
         if (error%has_error()) then
            call error%add_context("Fukui indices: the ions' functional")
            return
         end if
         ! A double hybrid's perturbative term needs an unrestricted MP2, which
         ! the CPU path does not have. Refused rather than run, because running
         ! it would not fail: the Kohn-Sham half is correct and unrestricted, so
         ! the ions would come back converged and short by the PT2 term. IP and
         ! EA are differences of total energies, so a term missing from two of
         ! the three states lands directly in every descriptor here.
         if (xc_ions%pt2_fraction /= 0.0_dp) then
            call error%set(ERROR_VALIDATION, "'"//trim(functional)//"' is a double "// &
                           "hybrid, and the Fukui ions are open shell. Its perturbative "// &
                           "term needs an unrestricted MP2 the CPU path does not have, "// &
                           "so the ions would converge with the Kohn-Sham part alone and "// &
                           "the IP and EA would be wrong by that term. Use a functional "// &
                           "without a perturbative term.")
            return
         end if
         xc_arg => xc_ions
         res%functional = trim(functional)
      end if

      natm = mol%natm
      res%scheme = scheme
      res%energy_neutral = neutral_energy

      call mol%overlap(overlap)
      call condense(mol, neutral_density, overlap, scheme, 0.0_dp, res%q_neutral, error)
      if (error%has_error()) return

      ! The anion. One more electron, and a doublet because the neutral was
      ! closed shell.
      call run_libcint_uhf(mol, nelec + 1, 2, max_iter, energy_tol, density_tol, &
                           .false., anion, error, xc=xc_arg)
      if (error%has_error()) then
         call error%add_context("Fukui indices: the anion")
         return
      end if
      if (.not. anion%converged) then
         call error%set(ERROR_VALIDATION, "the anion did not converge, so f+ cannot "// &
                        "be formed. Anions converge poorly in a basis with no diffuse "// &
                        "functions, which is usually the cause rather than the SCF.")
         return
      end if
      allocate (total(size(anion%density, 1), size(anion%density, 2)))
      total = anion%density + anion%density_beta
      call condense(mol, total, overlap, scheme, -1.0_dp, res%q_anion, error)
      if (error%has_error()) return
      res%energy_anion = anion%energy
      deallocate (total)

      ! The cation.
      call run_libcint_uhf(mol, nelec - 1, 2, max_iter, energy_tol, density_tol, &
                           .false., cation, error, xc=xc_arg)
      if (error%has_error()) then
         call error%add_context("Fukui indices: the cation")
         return
      end if
      if (.not. cation%converged) then
         call error%set(ERROR_VALIDATION, "the cation did not converge, so f- cannot "// &
                        "be formed.")
         return
      end if
      allocate (total(size(cation%density, 1), size(cation%density, 2)))
      total = cation%density + cation%density_beta
      call condense(mol, total, overlap, scheme, 1.0_dp, res%q_cation, error)
      if (error%has_error()) return
      res%energy_cation = cation%energy
      deallocate (total, overlap)

      allocate (res%f_plus(natm), res%f_minus(natm), res%f_zero(natm), res%dual(natm))
      res%f_plus = res%q_neutral - res%q_anion
      res%f_minus = res%q_cation - res%q_neutral
      res%f_zero = 0.5_dp*(res%f_plus + res%f_minus)
      res%dual = res%f_plus - res%f_minus

      res%ionisation_potential = res%energy_cation - res%energy_neutral
      res%electron_affinity = res%energy_neutral - res%energy_anion
      res%chemical_potential = -0.5_dp*(res%ionisation_potential + res%electron_affinity)
      res%hardness = res%ionisation_potential - res%electron_affinity
      if (abs(res%hardness) > 1.0e-10_dp) then
         res%electrophilicity = res%chemical_potential**2/(2.0_dp*res%hardness)
      end if
      res%anion_bound = res%electron_affinity > 0.0_dp

      call check_sum_rule(res, error)
   end subroutine fukui_indices

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
         ! The constraint matters here in a way it does not for a neutral: the
         ! fit is told the ion's total charge, so the difference between two
         ! states is a redistribution of exactly one electron rather than of one
         ! electron plus whatever the two fits disagreed about.
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
      !! for all of it, so `sum_A f_A+` and `sum_A f_A-` are one by construction.
      !! Nothing about reactivity is being tested here -- what it catches is a
      !! density that was not the total, an ion run at the wrong charge, or a
      !! population scheme that lost part of the molecule.
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
      ! comparable across methods and a table that does not say which one it
      ! used invites exactly that comparison.
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

      ! **Said rather than hidden.** A negative condensed index is not a site
      ! that repels electrons; it is the population analysis failing to divide
      ! the density sensibly, and it is common with Mulliken in a basis with
      ! diffuse functions -- a function centred on one atom reaching over
      ! another gets charged to the wrong one.
      if (any_negative) then
         call logger%warning("  some indices came out negative, which is the "// &
                             "population analysis struggling rather than a site that "// &
                             "repels charge")
         if (trim(res%scheme) == "mulliken") then
            call logger%warning("  Mulliken is especially prone to this; chelpg is "// &
                                "the default for that reason")
         end if
      end if

      ! The direct test rather than a guess from the basis-set name: if the
      ! anion sits above the neutral then nothing bound the extra electron, and
      ! f+ describes whatever orbital the basis had left over.
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
