!! One fragment into the CPU SCF and back out
module mqc_libcint_bridge
   !! Turns a `physical_fragment_t` into a libcint molecule, runs RHF over it,
   !! and fills in the result the rest of the program reads.
   !!
   !! The mirror of `mqc_cuest_bridge`, and the same shape on purpose: one
   !! entry point taking the settings the method already assembled, so which
   !! backend runs is a choice made in one place rather than a difference
   !! spread through the method layer.
   !!
   !! Caps are ordinary atoms here. A capped fragment arrives with its cap
   !! hydrogens already in `element_numbers` and `coordinates`, and nothing
   !! about the integrals cares which atoms were in the original molecule --
   !! the redistribution of forces back onto the heavy atoms happens later and
   !! elsewhere.
   use pic_logger, only: logger => global_logger
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t, SCF_CONVERGED, SCF_NOT_CONVERGED
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_number_to_symbol
   use mqc_program_limits, only: MAX_ELEMENT_SYMBOL_LEN
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf, run_libcint_uhf, &
                              SCF_GUESS_GWH, SCF_GUESS_SAC, SCF_GUESS_SAD, SCF_GUESS_PROJ
   use mqc_libcint_projection, only: climb_basis_ladder
   use mqc_libcint_atomic_guess, only: build_atomic_guess, parse_guess_name, &
                                       guess_display_name
   use mqc_libcint_xc, only: xc_context_t, xc_context_create, xc_available
   use mqc_libcint_gradient, only: libcint_scf_gradient
   use mqc_libcint_mp2_gradient, only: libcint_mp2_gradient
   use mqc_libcint_ri_mp2_gradient, only: libcint_ri_mp2_gradient
   use mqc_libcint_mp2, only: mp2_result_t, run_libcint_mp2, run_libcint_ri_mp2
   use mqc_libcint_cc, only: cc_result_t, run_libcint_ccsd
   use mqc_libcint_rcc, only: rcc_result_t, run_libcint_rccsd
   use mqc_libcint_bonding, only: run_quao_analysis, bonding_analysis_kind, &
                                  BONDING_GMS_QUAO
   use mqc_libcint_casci, only: casci_result_t, run_libcint_casci, &
                                run_libcint_ormas_ci
   use mqc_libcint_avas, only: avas_select, avas_result_t, valence_select
   use mqc_libcint_mcscf, only: casscf_result_t, run_libcint_casscf, &
                                natural_orbitals
   use mqc_program_limits, only: MAX_LINE_LENGTH
   implicit none
   private

   public :: run_libcint_hf
   public :: run_libcint_mcscf
   public :: run_libcint_fmo
   public :: run_libcint_makefp
   public :: run_libcint_charges
   public :: run_libcint_efp
   public :: run_libcint_sapt0
   public :: libcint_backend_available

   real(dp), parameter :: CI_TOLERANCE = 1.0e-11_dp
      !! Residual the CASCI Davidson is driven to on the CASCI-only path.
      !!
      !! The same number `mqc_libcint_mcscf` pins its own inner CASCI at, and
      !! deliberately not a keyword. A CASCI *is* its CI energy -- there is
      !! nothing downstream to absorb a loose solve the way a CASSCF macro-step
      !! absorbs one -- so the only defensible setting is "tight", and 1e-11 on
      !! a hundred-hartree total is inside what a threaded Fock build
      !! reproduces.

   real(dp), parameter :: ERI_CORE_BUDGET_BYTES = 2.0e9_dp
      !! How much the stored two-electron tensor may take before the SCF goes
      !! back to rebuilding it every iteration.
      !!
      !! The tensor is the full n^4, not the eightfold-unique set, so this is
      !! reached at 128 functions. Packing it would push that to 215, which is
      !! worth doing when something needs it -- `eris_packed` already exists for
      !! the correlated methods -- but the contraction in `build_fock` addresses
      !! the tensor as four indices and would have to be rewritten with it.
      !!
      !! Deliberately a fixed budget rather than a fraction of the machine.
      !! Fragmented runs put one MPI rank per fragment on a node, and sizing
      !! this from total memory would have every rank on a node each conclude it
      !! could have all of it.

contains

   pure function eri_fits_in_core(nao) result(fits)
      !! Whether n^4 stored integrals stay inside `ERI_CORE_BUDGET_BYTES`
      !!
      !! Computed in `real` rather than integer arithmetic: n^4 at four hundred
      !! functions overflows a 32-bit integer, and the symptom of that would be
      !! a large basis quietly deciding it fitted.
      integer, intent(in) :: nao
      logical :: fits

      fits = 8.0_dp*real(nao, dp)**4 <= ERI_CORE_BUDGET_BYTES
   end function eri_fits_in_core

   pure function libcint_backend_available() result(available)
      !! Whether this build can run an SCF on the CPU
      logical :: available
      available = .true.
   end function libcint_backend_available

   subroutine run_libcint_charges(atomic_numbers, element_symbols, coordinates, &
                                  basis_name, scheme, total_charge, charges, error)
      !! Atomic charges from an RHF density, by Mulliken or CHELPG
      !!
      !! Here rather than in `src/interface/` for the reason every other entry in
      !! this file is: the molecule builder, the SCF and both partition schemes
      !! live in the CPU backend. A C API reaching for them directly compiles
      !! under CMake with the backend on and nothing else -- not the stub build,
      !! and not FPM, which only sees `src/` and `app/`.
      !!
      !! Closed shell only; an odd electron count is refused rather than paired
      !! up, because a caller fragmenting a radical needs to know this cannot
      !! answer for it.
      use pic_types, only: dp
      use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
      use mqc_libcint_charges, only: mulliken_charges, chelpg_charges
      use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)     !! (3, n), Bohr
      character(len=*), intent(in) :: basis_name
      character(len=*), intent(in) :: scheme        !! "mulliken" or "chelpg"
      integer, intent(in) :: total_charge
      real(dp), allocatable, intent(out) :: charges(:)
      type(error_t), intent(inout) :: error

      integer, parameter :: SCF_MAX_ITER = 100
      real(dp), parameter :: SCF_ENERGY_TOL = 1.0e-9_dp
      real(dp), parameter :: SCF_DENSITY_TOL = 1.0e-7_dp

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: overlap(:, :)
      integer :: nelec

      nelec = sum(atomic_numbers) - total_charge
      if (mod(nelec, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "atomic charges come from a closed-shell "// &
                        "RHF and this system has an odd number of electrons")
         return
      end if

      call build_libcint_molecule(atomic_numbers, element_symbols, coordinates, &
                                  basis_name, mol, error)
      if (error%has_error()) return

      call run_libcint_rhf(mol, nelec, SCF_MAX_ITER, SCF_ENERGY_TOL, SCF_DENSITY_TOL, &
                           .false., scf, error)
      if (error%has_error()) return
      if (.not. scf%converged) then
         call error%set(ERROR_VALIDATION, "the SCF did not converge, so there is no "// &
                        "density to partition")
         return
      end if

      if (scheme == "mulliken") then
         call mol%overlap(overlap)
         call mulliken_charges(mol, scf%density, overlap, charges, error)
      else
         call chelpg_charges(mol, scf%density, charges, error, &
                             total_charge=real(total_charge, dp))
      end if
   end subroutine run_libcint_charges

   subroutine run_libcint_efp(potentials, fragment_sizes, fragment_atoms, &
                              coordinates, terms, error)
      !! The EFP2 interaction energy of a set of fragments the deck has placed
      !!
      !! Here rather than in the driver for the same reason `run_libcint_makefp` is:
      !! every module this needs lives in the CPU backend, so a driver reaching for
      !! them directly does not compile in a build without one. The stub beside this
      !! declines with the same signature.
      !!
      !! `terms` comes back as plain numbers rather than a derived type, because such
      !! a type would have to exist on both sides of that gate and there is nothing
      !! for it to hold that six reals do not.
      use pic_types, only: dp
      use mqc_program_limits, only: N_EFP_TERMS
      use mqc_efp_read, only: efp_fragment_t, read_efp_potential
      use mqc_efp_energy, only: efp_energy_t, efp_interaction_energy, place_fragment
      use mqc_efp_rotate, only: rotate_fragment
      character(len=*), intent(in) :: potentials(:)    !! One path per fragment
      integer, intent(in) :: fragment_sizes(:)
      integer, intent(in) :: fragment_atoms(:, :)      !! (max_size, n_frag), 0-based
      real(dp), intent(in) :: coordinates(:, :)        !! (3, n_atoms), Bohr
      real(dp), intent(out) :: terms(N_EFP_TERMS)
         !! electrostatics, polarization, exchange repulsion, dispersion,
         !! charge transfer, total
      type(error_t), intent(inout) :: error

      type(efp_fragment_t), allocatable :: frags(:)
      type(efp_energy_t) :: energy
      real(dp), allocatable :: shifts(:, :), own(:, :)
      real(dp) :: rot(3, 3)
      integer :: n, k, a, natom, i

      terms = 0.0_dp
      n = size(potentials)
      allocate (frags(n), shifts(3, n))

      do k = 1, n
         if (len_trim(potentials(k)) == 0) then
            call error%set(ERROR_VALIDATION, "efp: a fragment carries no potential. "// &
                           "A mixed quantum/EFP system is not supported yet -- every "// &
                           "fragment needs one.")
            return
         end if
         call read_efp_potential(trim(potentials(k)), frags(k), error)
         if (error%has_error()) return

         natom = fragment_sizes(k)
         allocate (own(3, natom))
         do i = 1, natom
            a = fragment_atoms(i, k) + 1        ! stored 0-based
            own(:, i) = coordinates(:, a)
         end do
         call place_fragment(frags(k), own, rot, shifts(:, k), error)
         deallocate (own)
         if (error%has_error()) return

         ! Turn the fragment into the deck's orientation before anything reads it.
         ! Everything a potential carries is in its own frame -- the multipoles, the
         ! polarizabilities at every rank, and the localized orbitals -- so this is
         ! not optional for any term.
         call rotate_fragment(frags(k), rot, error)
         if (error%has_error()) return
      end do

      energy = efp_interaction_energy(frags, shifts, error)
      if (error%has_error()) return

      terms = [energy%electrostatics, energy%polarization, energy%exchange_repulsion, &
               energy%dispersion, energy%charge_transfer, energy%total]

      do k = 1, n
         call frags(k)%destroy()
      end do
      deallocate (frags, shifts)
   end subroutine run_libcint_efp
   subroutine run_libcint_sapt0(z_a, sym_a, xyz_a, z_b, sym_b, xyz_b, basis_name, &
                                terms, error)
      !! SAPT0 between two monomers, as named physical terms
      !!
      !! Here rather than in the driver for the reason `run_libcint_makefp` is:
      !! every module this needs lives in the CPU backend, so a driver reaching
      !! for them directly would not compile without one. The stub beside this
      !! declines with the same signature.
      !!
      !! `terms` crosses that boundary as plain numbers rather than
      !! `sapt_terms_t`, because the type would have to exist on both sides of a
      !! gate whose whole purpose is that one side has none of this compiled.
      use pic_types, only: dp
      use mqc_program_limits, only: N_SAPT_TERMS
      use mqc_sapt, only: sapt_molecules_t, build_sapt_molecules, sapt_terms_t, &
                          run_sapt0
      integer, intent(in) :: z_a(:), z_b(:)
      character(len=*), intent(in) :: sym_a(:), sym_b(:)
      real(dp), intent(in) :: xyz_a(:, :), xyz_b(:, :)   !! (3, n), Bohr
      character(len=*), intent(in) :: basis_name
      real(dp), intent(out) :: terms(N_SAPT_TERMS)
         !! Ordered by `SAPT_TERM_NAMES`
      type(error_t), intent(inout) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_terms_t) :: t

      terms = 0.0_dp
      call build_sapt_molecules(z_a, sym_a, xyz_a, z_b, sym_b, xyz_b, &
                                basis_name, mols, error)
      if (error%has_error()) return
      call run_sapt0(mols, t, error)
      if (error%has_error()) then
         call mols%destroy()
         return
      end if

      terms = [t%elst10, t%exch10_s2, t%exch10, t%ind20_u, t%ind20_r, &
               t%exch_ind20_u, t%exch_ind20_r, t%disp20, t%exch_disp20, &
               t%delta_hf, t%e_int_hf, t%total]
      call mols%destroy()
   end subroutine run_libcint_sapt0

   subroutine run_libcint_makefp(atomic_numbers, element_symbols, coordinates, &
                                 basis_name, name, path, error, charge, verbose, &
                                 aux_basis, guess)
      !! Build an effective fragment potential and write it
      !!
      !! Here rather than in the driver so the driver needs no knowledge of whether
      !! this build has an integral backend -- the stub next to it declines with the
      !! same signature.
      use pic_types, only: dp
      use mqc_efp_potential, only: efp_potential_t, make_efp_potential, &
                                   write_efp_potential
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)     !! (3, natm), Bohr
      character(len=*), intent(in) :: basis_name
      character(len=*), intent(in) :: name          !! Fragment name for `$FRAGNAME`
      character(len=*), intent(in) :: path
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: charge   !! Net charge; a fragment may be an ion
      logical, intent(in), optional :: verbose
      character(len=*), intent(in), optional :: aux_basis
         !! Fit the response Hessian against this basis instead of building it
         !! exactly. Absent, the build is exact.
      character(len=*), intent(in), optional :: guess
         !! Initial-guess name from the deck, forwarded to the SCF. Absent leaves
         !! `make_efp_potential` on its "auto" (SAD) default.

      type(efp_potential_t) :: pot

      if (present(aux_basis)) then
         call make_efp_potential(atomic_numbers, element_symbols, coordinates, &
                                 basis_name, name, pot, error, charge=charge, &
                                 verbose=verbose, aux_basis=aux_basis, guess=guess)
      else
         call make_efp_potential(atomic_numbers, element_symbols, coordinates, &
                                 basis_name, name, pot, error, charge=charge, &
                                 verbose=verbose, guess=guess)
      end if
      if (error%has_error()) return
      call write_efp_potential(pot, path, error)
      call pot%destroy()
   end subroutine run_libcint_makefp

   subroutine run_libcint_fmo(atomic_numbers, element_symbols, coordinates, owner, &
                              basis_name, esp, expansion, far_field, resppc, &
                              level, max_outer, outer_tol, scf_max_iter, &
                              scf_energy_tol, scf_density_tol, energy, error, comm)
      !! Run FMO2 (or EE-MBE) over a partitioned system
      !!
      !! Options arrive as plain scalars rather than the backend's own options
      !! type, so the layer above never has to see a type it cannot compile
      !! without the backend. Coordinates are Bohr; `owner(i)` is atom i's
      !! fragment, numbered from one with no gaps.
      use mqc_libcint_fmo, only: fmo_options_t, fmo_result_t, run_fmo2
      use pic_types, only: dp
      use mqc_error, only: error_t
      use pic_mpi_lib, only: comm_t
      integer, intent(in) :: atomic_numbers(:)
      character(len=*), intent(in) :: element_symbols(:)
      real(dp), intent(in) :: coordinates(:, :)
      integer, intent(in) :: owner(:)
      character(len=*), intent(in) :: basis_name, esp, expansion, far_field
      real(dp), intent(in) :: resppc
      integer, intent(in) :: level
      integer, intent(in) :: max_outer
      real(dp), intent(in) :: outer_tol
      integer, intent(in) :: scf_max_iter
      real(dp), intent(in) :: scf_energy_tol, scf_density_tol
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: error
      type(comm_t), intent(in), optional :: comm
         !! Present means distribute the fragment work over this communicator.
         !! Absent means one rank does all of it.

      type(fmo_options_t) :: opts
      type(fmo_result_t) :: res
      character(len=2), allocatable :: symbols(:)
      integer :: i

      energy = 0.0_dp
      allocate (symbols(size(atomic_numbers)))
      do i = 1, size(atomic_numbers)
         symbols(i) = adjustl(element_symbols(i))
      end do

      opts%basis = basis_name
      opts%esp = esp
      opts%expansion = expansion
      opts%far_field = far_field
      opts%resppc = resppc
      opts%level = level
      opts%max_outer = max_outer
      opts%outer_tol = outer_tol
      opts%scf_max_iter = scf_max_iter
      opts%scf_energy_tol = scf_energy_tol
      opts%scf_density_tol = scf_density_tol

      call run_fmo2(atomic_numbers, symbols, coordinates, owner, opts, res, error, comm)
      if (error%has_error()) return
      energy = res%energy
   end subroutine run_libcint_fmo

   subroutine run_libcint_hf(settings, fragment, result, want_gradient)
      !! Closed-shell HF for one fragment, on the CPU
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient

      ! `aux` and `xc` are targets because the gradient takes both as optional
      ! arguments, and a disassociated pointer is how "this SCF fitted nothing"
      ! and "this SCF had no functional" are said without four spellings of the
      ! same call.
      type(libcint_molecule_t) :: mol, corr_aux
      type(libcint_molecule_t), target :: aux
      type(libcint_molecule_t), pointer :: aux_arg
      type(rhf_result_t) :: scf
      type(error_t) :: error
      character(len=MAX_ELEMENT_SYMBOL_LEN), allocatable :: symbols(:)
      integer :: iatom, diis_size, guess_kind
      logical :: unrestricted, do_gradient
      ! Left unallocated for the guesses that need no free-atom solutions. An
      ! unallocated allocatable passed to an optional dummy is an absent
      ! argument, so the SCF calls below need no branching on which guess ran.
      real(dp), allocatable :: guess_a(:, :), guess_b(:, :), guess_total(:, :)
      type(error_t) :: guess_error
      type(error_t) :: analysis_error
      type(xc_context_t), target :: xc
      type(xc_context_t), pointer :: xc_arg
      logical :: kohn_sham
      type(timer_type) :: grad_clock

      character(len=MAX_LINE_LENGTH) :: line

      do_gradient = .false.
      if (present(want_gradient)) do_gradient = want_gradient

      ! Coupled cluster stops here: its gradient needs the Lambda equations,
      ! which do not exist on this side at all. MP2's does exist, and is taken
      ! below -- but only where the reference underneath it is one the relaxed
      ! density was written for, which is the three refusals after this.
      if (do_gradient .and. settings%run_cc) then
         call result%error%set(ERROR_VALIDATION, "a coupled cluster gradient needs the "// &
                               "Lambda amplitudes, which are not implemented. Run the "// &
                               "gradient at the Hartree-Fock or MP2 level, or coupled "// &
                               "cluster as an energy.")
         result%has_error = .true.
         return
      end if
      if (do_gradient .and. settings%run_mp2) then
         ! Four combinations of what is fitted, each a different energy with a
         ! different gradient, and all four implemented. Exact reference with
         ! exact correlation is `libcint_mp2_gradient`; exact reference with
         ! fitted correlation is an `ri-mp2` deck and is
         ! `libcint_ri_mp2_gradient`. Adding `keywords.scf.density_fitting` to
         ! either fits the reference too, which moves the response operator,
         ! both potentials built from the relaxed density, and the reference's
         ! two-electron derivative term onto the auxiliary basis.
         !
         ! What that last pair does *not* do is make the conventional gradient
         ! cheap: its two-particle density stays four-index and contracts
         ! against four-centre derivatives whatever the reference fitted. It
         ! makes it consistent, which is the point -- differentiating an exact
         ! reference under a fitted SCF is the derivative of an energy nothing
         ! computed.
         if (settings%freeze_core .and. settings%n_frozen_core /= 0) then
            call result%error%set(ERROR_VALIDATION, "the MP2 gradient does not implement "// &
                                  "a frozen core: the relaxed density gains "// &
                                  "occupied-frozen and virtual-frozen blocks that are not "// &
                                  "built. Run the MP2 gradient with freeze_core off.")
            result%has_error = .true.
            return
         end if
      end if

      ! Same rule the GPU backend applies, so a deck does not change meaning
      ! when it moves between them: an odd electron count or a multiplicity
      ! above one has no restricted solution to find, whatever the keyword says.
      unrestricted = (fragment%multiplicity /= 1) .or. (mod(fragment%nelec, 2) /= 0) &
                     .or. settings%unrestricted

      ! Continuum solvation belongs to cuEST on this side: the CPU path has no
      ! cavity, no surface charges and no attenuated one-electron integrals over
      ! them. Refused rather than ignored, because ignoring it returns a
      ! gas-phase energy for a deck that asked for solvent -- of the right
      ! magnitude, converged, and with nothing in the output to say the solvent
      ! was dropped.
      if (settings%pcm%enabled) then
         call result%error%set(ERROR_VALIDATION, "continuum solvation (keywords.pcm) is "// &
                               "implemented on the cuEST backend only; the CPU path has "// &
                               "no cavity. Refused rather than run in the gas phase, "// &
                               "which is what ignoring it would silently give you.")
         result%has_error = .true.
         return
      end if

      if (unrestricted .and. settings%density_fitting) then
         call result%error%set(ERROR_VALIDATION, "unrestricted density fitting is not "// &
                               "implemented on the CPU backend: the fitted J and K are "// &
                               "written for one density. Run unrestricted without "// &
                               "density_fitting, or restricted with it.")
         result%has_error = .true.
         return
      end if

      ! Refused rather than quietly downgraded, on the same principle as the
      ! unrestricted density fitting above. Coupled cluster here is written over a
      ! restricted reference, and the fitted route does not exist yet at all -- a
      ! deck asking for either would otherwise get a Hartree-Fock energy, or a
      ! conventional CCSD reported as RI, and nothing in the output would say so.
      if (settings%run_cc .and. unrestricted) then
         call result%error%set(ERROR_VALIDATION, "coupled cluster on the CPU backend "// &
                               "needs a restricted reference: it is written in spin "// &
                               "orbitals over RHF orbitals, and an unrestricted "// &
                               "reference needs its own alpha and beta transform. Run "// &
                               "a closed-shell system, or use density functional theory.")
         result%has_error = .true.
         return
      end if

      ! And MP2 for the same reason, which it did not used to be. The transform
      ! takes one set of orbitals and an occupied count, so an unrestricted
      ! reference gave it the alpha orbitals and `nelec/2` -- which for an odd
      ! electron count truncates. OH in sto-3g came back as a converged
      ! -74.416097 computed over three doubly-occupied orbitals instead of nine
      ! electrons, with nothing in the output to say which. It is the same defect
      ! the double-hybrid guard below refuses, and it is worse here because MP2 is
      ! asked for directly.
      if (settings%run_mp2 .and. unrestricted) then
         call result%error%set(ERROR_VALIDATION, "MP2 on the CPU backend needs a "// &
                               "restricted reference: the four-index transform takes one "// &
                               "set of orbitals with an occupied count, and an "// &
                               "unrestricted reference needs separate alpha and beta "// &
                               "transforms. Run a closed-shell system, or use density "// &
                               "functional theory, which is unrestricted here.")
         result%has_error = .true.
         return
      end if
      allocate (symbols(fragment%n_atoms))
      do iatom = 1, fragment%n_atoms
         symbols(iatom) = element_number_to_symbol(fragment%element_numbers(iatom))
      end do

      ! `ghost` present and all-false is the same molecule as no ghost at all,
      ! so an ordinary fragment is unaffected -- and the AO count and ordering
      ! are identical either way, which is what makes the counterpoise
      ! difference meaningful.
      call build_libcint_molecule(fragment%element_numbers, symbols, &
                                  fragment%coordinates, trim(settings%basis_set), &
                                  mol, error, ghost=ghost_of(fragment))
      if (error%has_error()) then
         call result%error%set(ERROR_VALIDATION, error%get_message())
         result%has_error = .true.
         return
      end if

      ! The SCF banner below reports the basis but not who it is solving for:
      ! a charged or open-shell fragment reads identically to a neutral one
      ! without this, which is exactly the case a fragmented run most needs to
      ! see. Same integers the breakdown CSV carries.
      call logger%verbose("charge "//to_text(fragment%charge)// &
                          ", multiplicity "//to_text(fragment%multiplicity)// &
                          ", electrons "//to_text(fragment%nelec))

      diis_size = settings%diis_size
      if (.not. settings%use_diis) diis_size = 0

      ! ---- Kohn-Sham or Hartree-Fock? ---------------------------------------
      !
      ! A named functional is the whole difference. Nothing else about this
      ! routine changes: the SCF takes the context as one optional argument, so
      ! Hartree-Fock is the case where there is no functional to build one from.
      kohn_sham = len_trim(settings%functional) > 0
      if (kohn_sham) then
         if (.not. xc_available()) then
            call result%error%set(ERROR_VALIDATION, "a functional was requested ('"// &
                                  trim(settings%functional)//"') but this build has no "// &
                                  "libxc: configure with -DMQC_ENABLE_LIBXC=ON")
            result%has_error = .true.
            call mol%destroy()
            return
         end if
         ! Spin-polarised exactly when the SCF is. libxc fixes the spin channel
         ! when a functional is initialised, so this cannot be decided later by
         ! whoever evaluates the potential -- and the two evaluators refuse the
         ! wrong kind of context rather than reading it with the wrong stride.
         call xc_context_create(mol, trim(settings%functional), xc, error, &
                                level=settings%grid_level, polarized=unrestricted)
         if (error%has_error()) then
            call result%error%set(ERROR_VALIDATION, error%get_message())
            result%has_error = .true.
            call mol%destroy()
            return
         end if
         ! A double hybrid on an open shell would need an unrestricted MP2, and the
         ! perturbative term below is written over restricted orbitals. Refused
         ! here rather than downstream, because downstream it would not fail: it
         ! would take the alpha orbitals for a closed-shell set, pair them with
         ! themselves and return a plausible number. The Kohn-Sham half of the
         ! functional is unrestricted and correct, which is exactly what makes the
         ! total convincing.
         if (unrestricted .and. xc%pt2_fraction /= 0.0_dp) then
            call result%error%set(ERROR_VALIDATION, "'"//trim(settings%functional)// &
                                  "' is a double hybrid and its perturbative term needs "// &
                                  "an unrestricted MP2, which the CPU path does not have. "// &
                                  "The Kohn-Sham part alone would be ~65 mHartree short, "// &
                                  "so this is refused rather than reported. Run a "// &
                                  "closed-shell system, or a functional without a "// &
                                  "perturbative term.")
            result%has_error = .true.
            call xc%destroy()
            call mol%destroy()
            return
         end if
         if (settings%verbose) then
            write (line, "(a,a,a,i0)") "  functional: ", trim(settings%functional), ", grid level ", settings%grid_level
            call logger%info(trim(line))
         end if
      end if

      ! ---- which initial guess? ---------------------------------------------
      call parse_guess_name(settings%guess, guess_kind, error)
      if (error%has_error()) then
         call result%error%set(ERROR_VALIDATION, error%get_message())
         result%has_error = .true.
         call mol%destroy()
         return
      end if

      if (guess_kind == SCF_GUESS_PROJ) then
         if (.not. allocated(settings%guess_steps)) then
            call result%error%set(ERROR_VALIDATION, "guess 'basis_set_projection' needs "// &
                                  "keywords.guess.subscf.steps; the schema normally "// &
                                  "refuses a deck without them, so this settings object "// &
                                  "was built by something other than the JSON reader")
            result%has_error = .true.
            call mol%destroy()
            return
         end if
         call climb_basis_ladder(settings%guess_steps, fragment%element_numbers, symbols, &
                                 fragment%coordinates, fragment%nelec, mol, guess_total, &
                                 guess_error, verbose=settings%verbose)
         if (guess_error%has_error()) then
            ! The same reasoning the atomic guess uses: a guess that cannot be
            ! built is a reason to start somewhere else, not to fail the
            ! calculation. Loud, because the run is now doing something other
            ! than what the deck asked for.
            call logger%warning("basis projection guess: "//guess_error%get_message()// &
                                " -- falling back to sad")
            guess_kind = SCF_GUESS_SAD
            if (allocated(guess_total)) deallocate (guess_total)
         end if
      end if

      if (guess_kind == SCF_GUESS_SAC .or. guess_kind == SCF_GUESS_SAD) then
         call build_atomic_guess(mol, guess_kind, guess_a, guess_b, guess_error)
         if (guess_error%has_error()) then
            ! A free atom that will not converge is a reason to start somewhere
            ! else, not to fail a molecular calculation -- the guess is a
            ! convergence aid and cannot change what the SCF converges *to*
            ! except by changing which solution it finds. That exception is why
            ! this is a warning and not silence: an open-shell system really can
            ! land on a different state from a different starting point, and a
            ! quietly substituted guess would be a quietly different answer.
            call logger%warning("initial guess: "//guess_error%get_message()// &
                                " -- falling back to gwh")
            guess_kind = SCF_GUESS_GWH
            if (allocated(guess_a)) deallocate (guess_a)
            if (allocated(guess_b)) deallocate (guess_b)
         else if (.not. unrestricted) then
            ! The restricted path takes one total density. Summing here rather
            ! than inside the SCF keeps the spin split where it means something.
            guess_total = guess_a + guess_b
         end if
      end if

      ! Reported rather than inferred from the iteration count. Which guess ran is
      ! the first thing worth knowing about an SCF that took longer than expected,
      ! and after a fallback it is the only place the substitution is visible.
      call logger%verbose("initial guess: "//guess_display_name(guess_kind))

      ! Density fitting is asked for, not inferred. aux_basis_set carries a
      ! default, so treating its presence as the request would mean every
      ! calculation quietly fitted -- and the difference is 5e-5 Hartree, large
      ! enough to matter and small enough to look like convergence noise.
      if (settings%density_fitting) then
         ! The auxiliary set follows the orbital set: a ghost centre carries
         ! fitting functions too, or the fitted Coulomb would see a different
         ! molecule from the one being fitted.
         call build_libcint_molecule(fragment%element_numbers, symbols, &
                                     fragment%coordinates, trim(settings%aux_basis_set), &
                                     aux, error, ghost=ghost_of(fragment))
         if (error%has_error()) then
            call result%error%set(ERROR_VALIDATION, "auxiliary basis '"// &
                                  trim(settings%aux_basis_set)//"': "//error%get_message())
            result%has_error = .true.
            call mol%destroy()
            return
         end if

         ! Refused here rather than deeper down, because this is the only place
         ! that knows both basis set *names* -- and a name is what has to change
         ! to fix it. 6-31G* is Cartesian while every JKFIT set is spherical, so
         ! this is the pairing a deck lands on most easily.
         if (mol%cartesian .neqv. aux%cartesian) then
            call result%error%set(ERROR_VALIDATION, "density fitting: the orbital basis '"// &
                                  trim(settings%basis_set)//"' is "// &
                                  angular_form_name(mol%cartesian)//" and the auxiliary basis '"// &
                                  trim(settings%aux_basis_set)//"' is "// &
                                  angular_form_name(aux%cartesian)//". The fitting integrals "// &
                                  "are built in one angular form for all three centres, so "// &
                                  "the two bases have to agree.")
            result%has_error = .true.
            call mol%destroy()
            call aux%destroy()
            return
         end if

         call run_libcint_rhf(mol, fragment%nelec, settings%max_iter, settings%energy_tol, &
                              settings%density_tol, settings%verbose, scf, error, &
                              aux=aux, diis_vectors=diis_size, guess=guess_kind, &
                              guess_density=guess_total, xc=xc)
         ! Kept alive: the gradient below has to be told the same auxiliary
         ! basis this SCF fitted with. Released once past it.
      else if (unrestricted) then
         call run_libcint_uhf(mol, fragment%nelec, fragment%multiplicity, settings%max_iter, &
                              settings%energy_tol, settings%density_tol, settings%verbose, &
                              scf, error, diis_vectors=diis_size, guess=guess_kind, &
                              guess_density_alpha=guess_a, guess_density_beta=guess_b, xc=xc)
      else
         ! Store the integrals when they fit, rather than rebuilding every
         ! quartet at every iteration.
         !
         ! A direct build costs the same on iteration twelve as on iteration one
         ! -- the screening is on the basis and the density, and neither gets
         ! much smaller -- so a twelve-iteration SCF evaluates the same million
         ! shell quartets twelve times. Storing them turns every iteration after
         ! the first into a contraction that is memory-bound rather than
         ! integral-bound: 1.92 s to 0.07 s on the case this was measured on,
         ! and the SCF as a whole from 25.6 s to 3.7 s.
         !
         ! The threshold is memory, because that is the only thing that makes
         ! this the wrong choice. It is not a fidelity trade like density
         ! fitting: the stored tensor holds the same integrals the direct build
         ! would have computed, and the energy agrees to 1e-11.
         !
         ! Range separation is the one exception, and it is a refusal rather
         ! than a preference: the SCF *errors out* on an in-core tensor with an
         ! attenuated kernel, because the tensor is built for the full Coulomb
         ! one and the long-range exchange would be silently missing. Deciding
         ! that here keeps a functional like wB97X working -- it would otherwise
         ! have started failing the moment this became the default.
         call run_libcint_rhf(mol, fragment%nelec, settings%max_iter, settings%energy_tol, &
                              settings%density_tol, settings%verbose, scf, error, &
                              diis_vectors=diis_size, guess=guess_kind, &
                              guess_density=guess_total, xc=xc, &
                              in_core=eri_fits_in_core(mol%nao) .and. .not. xc%range_separated)
      end if
      if (error%has_error()) then
         call result%error%set(ERROR_VALIDATION, error%get_message())
         result%has_error = .true.
         call aux%destroy()
         call mol%destroy()
         return
      end if

      ! Same policy as every other backend: a non-converged SCF is refused
      ! unless the run asked to keep it, because its energy has the right
      ! magnitude and nothing downstream can tell.
      if (scf%converged) then
         result%scf_status = SCF_CONVERGED
      else
         result%scf_status = SCF_NOT_CONVERGED
      end if
      result%scf_iterations = scf%iterations
      if (.not. scf%converged .and. .not. settings%allow_crap_scf) then
         call result%error%set(ERROR_VALIDATION, "SCF not converged in "// &
                               to_text(scf%iterations)//" cycles")
         result%has_error = .true.
         call aux%destroy()
         call mol%destroy()
         return
      end if

      result%energy%scf = scf%energy
      result%has_energy = .true.

      ! ---- analyses asked for alongside the energy ---------------------------
      !
      ! Run here, on the converged orbitals, because this is where they exist.
      ! A failure is reported and then dropped rather than propagated: the
      ! energy above is correct and already stored, and an analysis that cannot
      ! run -- an element the minimal basis does not cover, an open shell the
      ! construction is not written for -- is not a reason to throw it away.
      if (bonding_analysis_kind(settings%bonding_analysis) == BONDING_GMS_QUAO) then
         call analysis_error%clear()
         ! Printed whatever the verbosity. Asking for a property *is* the
         ! request for its output -- a deck that names an analysis and gets a
         ! silent run has been given nothing, and the natural reading is that
         ! the analysis found nothing to say.
         call run_quao_analysis(mol, fragment%element_numbers, symbols, &
                                fragment%coordinates, scf%orbitals, fragment%nelec, &
                                analysis_error, threshold=settings%bonding_threshold, &
                                no_sharing=settings%bonding_no_sharing)
         if (analysis_error%has_error()) then
            call logger%warning("  the bonding analysis could not run: "// &
                                analysis_error%get_message())
         end if
      end if

      ! ---- dE/dR, analytically ----------------------------------------------
      !
      ! Differentiating the energy that was just converged, which is why the
      ! auxiliary basis handed over is the one the SCF fitted with rather than
      ! whatever the deck names: a fitted energy differentiated as an exact one
      ! is a gradient of a function nobody computed.
      if (do_gradient .and. .not. settings%run_mp2) then
         ! A double hybrid's energy is the Kohn-Sham part plus a scaled PT2
         ! correlation, and the call below differentiates only the first. The
         ! second is added afterwards, once the reference gradient exists --
         ! see the block after this one. What is still refused there is the
         ! combination this cannot do, rather than the whole functional.
         !
         ! Left unchecked entirely, this used to return the hybrid's gradient
         ! against the double hybrid's energy: on water/cc-pVDZ 0.011 against a
         ! true 0.016, right in shape, wrong by a third, with nothing in the
         ! output to say which of the two numbers the other belonged to.
         if (kohn_sham .and. xc%pt2_fraction /= 0.0_dp) then
            if (unrestricted) then
               call result%error%set(ERROR_VALIDATION, "a double hybrid gradient over "// &
                                     "an open shell needs an unrestricted MP2 relaxed "// &
                                     "density and a spin-resolved response, neither of "// &
                                     "which is implemented. Run it as an energy.")
               result%has_error = .true.
               call mol%destroy()
               return
            end if
            if (xc%range_separated) then
               call result%error%set(ERROR_VALIDATION, "a range-separated double "// &
                                     "hybrid gradient would need the response operator "// &
                                     "and the potential derivative at a screened "// &
                                     "omega, which are not built. Run it as an energy.")
               result%has_error = .true.
               call mol%destroy()
               return
            end if
            if (xc%any_mgga) then
               call result%error%set(ERROR_VALIDATION, "a meta-GGA double hybrid "// &
                                     "gradient needs the exchange-correlation kernel "// &
                                     "in tau, which is not implemented. Run it as an "// &
                                     "energy, or take the gradient at a GGA-based "// &
                                     "double hybrid.")
               result%has_error = .true.
               call mol%destroy()
               return
            end if
            ! A fitted reference used to be refused here. It is not any more:
            ! the perturbative term's response equations are solved with the
            ! operator the SCF actually used, and its two-electron derivative
            ! term comes off the auxiliary basis. What still holds is that the
            ! *correlation* is exact either way -- a double hybrid's PT2 term is
            ! a conventional MP2, and in a gradient run it stays one.
            if (settings%freeze_core .and. settings%n_frozen_core /= 0) then
               call result%error%set(ERROR_VALIDATION, "a double hybrid gradient is "// &
                                     "all-electron: the relaxed density gains "// &
                                     "occupied-frozen blocks that are not built. The "// &
                                     "*energy* does honour freeze_core, so this refuses "// &
                                     "rather than returning the all-electron gradient of "// &
                                     "a frozen-core energy. Run the energy with it, or "// &
                                     "the gradient with freeze_core off.")
               result%has_error = .true.
               call mol%destroy()
               return
            end if
         end if

         call logger%info("  energy converged, computing the gradient")
         call grad_clock%start()

         aux_arg => null()
         xc_arg => null()
         if (settings%density_fitting) aux_arg => aux
         if (kohn_sham) xc_arg => xc

         if (unrestricted) then
            call libcint_scf_gradient(mol, scf%density, density_beta=scf%density_beta, &
                                      orbitals=scf%orbitals, orbitals_beta=scf%orbitals_beta, &
                                      orbital_energies=scf%orbital_energies, &
                                      orbital_energies_beta=scf%orbital_energies_beta, &
                                      n_occupied=scf%n_occupied, &
                                      n_occupied_beta=scf%n_occupied_beta, &
                                      gradient=result%gradient, error=error, &
                                      aux=aux_arg, xc=xc_arg)
         else
            call libcint_scf_gradient(mol, scf%density, orbitals=scf%orbitals, &
                                      orbital_energies=scf%orbital_energies, &
                                      n_occupied=scf%n_occupied, &
                                      gradient=result%gradient, error=error, &
                                      aux=aux_arg, xc=xc_arg)
         end if
         if (error%has_error()) then
            call result%error%set(ERROR_VALIDATION, "gradient: "//error%get_message())
            result%has_error = .true.
            if (kohn_sham) call xc%destroy()
            call aux%destroy()
            call mol%destroy()
            return
         end if
         result%has_gradient = .true.
         write (line, "(a,f10.2,a)") "  gradient done in ", grad_clock%get_elapsed_time(), " s"
         call logger%info(trim(line))
         if (settings%verbose) then
            write (line, "(a,f20.12)") "  |gradient|     ", sqrt(sum(result%gradient**2))
            call logger%info(trim(line))
         end if
      end if
      call aux%destroy()

      if (settings%run_mp2) then
         block
            type(mp2_result_t) :: mp2
            integer :: frozen

            frozen = settings%n_frozen_core
            if (frozen < 0) frozen = core_orbital_count(fragment%element_numbers)
            if (.not. settings%freeze_core) frozen = 0

            if (settings%corr_density_fitting) then
               call correlation_aux_basis(settings, fragment, symbols, corr_aux, error)
               if (error%has_error()) then
                  call result%error%set(ERROR_VALIDATION, error%get_message())
                  result%has_error = .true.
                  call mol%destroy()
                  return
               end if
               call run_libcint_ri_mp2(mol, corr_aux, scf%orbitals, scf%orbital_energies, &
                                       fragment%nelec/2, scf%energy, mp2, error, &
                                       n_frozen=frozen)
               call corr_aux%destroy()
            else
               call run_libcint_mp2(mol, scf%orbitals, scf%orbital_energies, &
                                    fragment%nelec/2, scf%energy, mp2, error, &
                                    n_frozen=frozen)
            end if
            if (error%has_error()) then
               call result%error%set(ERROR_VALIDATION, "MP2: "//error%get_message())
               result%has_error = .true.
               call mol%destroy()
               return
            end if
            ! energy_t sums scf + mp2%total(), so only the components go in.
            result%energy%mp2%ss = mp2%same_spin
            result%energy%mp2%os = mp2%opposite_spin
            result%energy%mp2%ss_scale = settings%scs_ss
            result%energy%mp2%os_scale = settings%scs_os

            ! The gradient of what was just computed, and only where that is
            ! literally true: a scaled MP2 is a different energy, and its
            ! gradient is not this one scaled -- the amplitudes enter the
            ! relaxed density and the Z-vector unscaled.
            if (do_gradient) then
               if (settings%scs_ss /= 1.0_dp .or. settings%scs_os /= 1.0_dp) then
                  call result%error%set(ERROR_VALIDATION, "a spin-scaled MP2 gradient "// &
                                        "is not the unscaled one rescaled: the "// &
                                        "amplitudes enter the relaxed density and the "// &
                                        "response equations, where the two spin cases "// &
                                        "are not separable afterwards. Run SCS-MP2 as "// &
                                        "an energy.")
                  result%has_error = .true.
                  call mol%destroy()
                  return
               end if
               call logger%info("  energy converged, computing the gradient")
               call grad_clock%start()
               ! Fitted correlation gets the fitted gradient. Not a refinement:
               ! `run_libcint_ri_mp2` above computed a fitted energy, and the
               ! conventional gradient is the derivative of a different one.
               ! The auxiliary basis is rebuilt here rather than held from the
               ! energy step, so each branch owns and destroys its own and no
               ! error return in between has to remember to.
               if (settings%corr_density_fitting) then
                  call correlation_aux_basis(settings, fragment, symbols, corr_aux, error)
                  if (error%has_error()) then
                     call result%error%set(ERROR_VALIDATION, error%get_message())
                     result%has_error = .true.
                     call mol%destroy()
                     return
                  end if
                  ! One auxiliary basis for both halves when both are fitted,
                  ! which is not a simplification: `model.aux_basis` is the only
                  ! place a fitting set is named, so it is the only combination
                  ! a deck can express.
                  call libcint_ri_mp2_gradient(mol, corr_aux, scf%orbitals, &
                                               scf%orbital_energies, fragment%nelec/2, &
                                               result%gradient, error, n_frozen=frozen, &
                                               fitted_reference=settings%density_fitting)
                  call corr_aux%destroy()
               else if (settings%density_fitting) then
                  ! Exact correlation over a fitted reference. The correlation's
                  ! four-centre integrals are still built -- fitting the SCF makes
                  ! this consistent rather than cheap -- and only the reference
                  ! half moves onto the auxiliary basis.
                  call correlation_aux_basis(settings, fragment, symbols, corr_aux, error)
                  if (error%has_error()) then
                     call result%error%set(ERROR_VALIDATION, error%get_message())
                     result%has_error = .true.
                     call mol%destroy()
                     return
                  end if
                  call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, &
                                            fragment%nelec/2, result%gradient, error, &
                                            n_frozen=frozen, aux=corr_aux)
                  call corr_aux%destroy()
               else
                  call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, &
                                            fragment%nelec/2, result%gradient, error, &
                                            n_frozen=frozen)
               end if
               if (error%has_error()) then
                  call result%error%set(ERROR_VALIDATION, "MP2 gradient: "// &
                                        error%get_message())
                  result%has_error = .true.
                  call mol%destroy()
                  return
               end if
               result%has_gradient = .true.
               write (line, "(a,f10.2,a)") "  gradient done in ", grad_clock%get_elapsed_time(), " s"
               call logger%info(trim(line))
               if (settings%verbose) then
                  write (line, "(a,f20.12)") "  |gradient|     ", sqrt(sum(result%gradient**2))
                  call logger%info(trim(line))
               end if
            end if

            if (settings%verbose) then
               write (line, "(a,i0,a,i0,a,i0)") "  MP2: frozen ", mp2%n_frozen, "  occupied ", &
                  mp2%n_occupied, "  virtual ", mp2%n_virtual
               call logger%info(trim(line))
               write (line, "(a,f20.12)") "  E(OS)          ", mp2%opposite_spin
               call logger%info(trim(line))
               write (line, "(a,f20.12)") "  E(SS)          ", mp2%same_spin
               call logger%info(trim(line))
               write (line, "(a,f20.12)") "  E(corr)        ", mp2%correlation
               call logger%info(trim(line))
               if (settings%scs_ss /= 1.0_dp .or. settings%scs_os /= 1.0_dp) then
                  write (line, "(a,f6.3,a,f6.3)") "  spin scaling:  SS x", settings%scs_ss, "   OS x", settings%scs_os
                  call logger%info(trim(line))
                  write (line, "(a,f20.12)") "  E(corr) scaled ", result%energy%mp2%total()
                  call logger%info(trim(line))
               end if
               write (line, "(a,f20.12)") "  E total        ", result%energy%total()
               call logger%info(trim(line))
            end if
         end block
      end if

      if (settings%run_cc) then
         block
            type(cc_result_t) :: cc
            type(rcc_result_t) :: rcc
            integer :: frozen, cc_iterations
            logical :: cc_converged, spin_adapted
            real(dp) :: cc_mp2, cc_singles, cc_doubles, cc_triples

            ! The same frozen-core rule the MP2 block applies, deliberately
            ! duplicated rather than hoisted: the two blocks are independent, and
            ! a shared local would silently couple them if one ever wanted a
            ! different count.
            frozen = settings%n_frozen_core
            if (frozen < 0) frozen = core_orbital_count(fragment%element_numbers)
            if (.not. settings%freeze_core) frozen = 0

            ! Which formulation. Both are exact for a closed shell and agree to
            ! machine precision -- that identity is asserted by
            ! test_mqc_libcint_rcc -- so this chooses how the same number is
            ! computed, not which number. Spatial by default because it is
            ! smaller and faster; the spin-orbital path is what a doubtful
            ! result gets checked against.
            spin_adapted = settings%cc_spin_adapted

            if (settings%corr_density_fitting) then
               call correlation_aux_basis(settings, fragment, symbols, corr_aux, error)
               if (error%has_error()) then
                  call result%error%set(ERROR_VALIDATION, error%get_message())
                  result%has_error = .true.
                  call mol%destroy()
                  return
               end if
               if (spin_adapted) then
                  call run_libcint_rccsd(mol, scf%orbitals, scf%orbital_energies, &
                                         fragment%nelec/2, frozen, settings%cc_max_iter, &
                                         settings%cc_tolerance, settings%cc_triples, &
                                         settings%verbose, rcc, error, &
                                         diis_vectors=settings%cc_diis_size, aux=corr_aux)
               else
                  call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, &
                                        fragment%nelec/2, frozen, settings%cc_max_iter, &
                                        settings%cc_tolerance, settings%cc_triples, &
                                        settings%verbose, cc, error, &
                                        diis_vectors=settings%cc_diis_size, aux=corr_aux)
               end if
               call corr_aux%destroy()
            else
               if (spin_adapted) then
                  call run_libcint_rccsd(mol, scf%orbitals, scf%orbital_energies, &
                                         fragment%nelec/2, frozen, settings%cc_max_iter, &
                                         settings%cc_tolerance, settings%cc_triples, &
                                         settings%verbose, rcc, error, &
                                         diis_vectors=settings%cc_diis_size)
               else
                  call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, &
                                        fragment%nelec/2, frozen, settings%cc_max_iter, &
                                        settings%cc_tolerance, settings%cc_triples, &
                                        settings%verbose, cc, error, &
                                        diis_vectors=settings%cc_diis_size)
               end if
            end if
            if (error%has_error()) then
               call result%error%set(ERROR_VALIDATION, "CCSD: "//error%get_message())
               result%has_error = .true.
               call mol%destroy()
               return
            end if

            ! The two result types carry the same components under the same
            ! names, so everything below this line is written once.
            if (spin_adapted) then
               cc_mp2 = rcc%e_mp2
               cc_singles = rcc%e_singles
               cc_doubles = rcc%e_doubles
               cc_triples = rcc%e_triples
               cc_iterations = rcc%iterations
               cc_converged = rcc%converged
            else
               cc_mp2 = cc%e_mp2
               cc_singles = cc%e_singles
               cc_doubles = cc%e_doubles
               cc_triples = cc%e_triples
               cc_iterations = cc%iterations
               cc_converged = cc%converged
            end if

            ! energy_t sums scf + mp2%total() + cc%total(), so only the components
            ! go in. All three are filled rather than a lumped correlation energy:
            ! a total cannot be taken apart afterwards, and the singles/doubles
            ! split is what says whether the T1 amplitudes are doing anything.
            result%energy%cc%singles = cc_singles
            result%energy%cc%doubles = cc_doubles
            result%energy%cc%triples = cc_triples
            if (settings%verbose) then
               ! Which formulation ran is said out loud. The two agree to
               ! machine precision, so nothing in the numbers below would
               ! otherwise reveal it -- and it is exactly what someone
               ! comparing a timing or a memory figure needs to know.
               if (spin_adapted) then
                  call logger%info("  CCSD: spin-adapted (spatial orbitals)")
               else
                  call logger%info("  CCSD: spin orbitals")
               end if
               write (line, "(a,i0,a,l1)") "  CCSD: iterations ", cc_iterations, "  converged ", cc_converged
               call logger%info(trim(line))
               write (line, "(a,f20.12)") "  E(MP2)         ", cc_mp2
               call logger%info(trim(line))
               write (line, "(a,f20.12)") "  E(singles)     ", cc_singles
               call logger%info(trim(line))
               write (line, "(a,f20.12)") "  E(doubles)     ", cc_doubles
               call logger%info(trim(line))
               if (settings%cc_triples) then
                  write (line, "(a,f20.12)") "  E(T)           ", cc_triples
                  call logger%info(trim(line))
               end if
               write (line, "(a,f20.12)") "  E(corr)        ", result%energy%cc%total()
               call logger%info(trim(line))
               write (line, "(a,f20.12)") "  E total        ", result%energy%total()
               call logger%info(trim(line))
            end if
         end block
      end if

      ! ---- a double hybrid's perturbative term -------------------------------
      !
      ! Part of the functional, not a correction on top of it, and the reason this
      ! cannot be left to the caller: a deck asking for B2PLYP that received only
      ! the Kohn-Sham part would get a converged energy 65 mHartree wrong, which is
      ! the exact shape of failure this backend is built to refuse.
      if (kohn_sham) then
         if (xc%pt2_fraction /= 0.0_dp) then
            block
               type(mp2_result_t) :: dh_mp2
               integer :: dh_frozen

               dh_frozen = settings%n_frozen_core
               if (dh_frozen < 0) dh_frozen = core_orbital_count(fragment%element_numbers)
               if (.not. settings%freeze_core) dh_frozen = 0

               ! Fitted when the deck *named* an auxiliary basis, exact otherwise
               ! -- and `named` rather than merely present, because
               ! `scf_config_t%aux_basis_set` carries a default that cuEST needs
               ! and every deck therefore has one. Testing for its presence
               ! fitted this term with a JKFIT set nobody asked for.
               !
               ! The reference implementations of these functionals are
               ! density-fitted, so fitting is the comparable choice rather than
               ! a compromise -- but a deck that named no auxiliary basis should
               ! get an answer rather than a refusal, and the conventional
               ! transform is the more accurate of the two anyway.
               !
               ! **A gradient run makes the same choice**, which it did not used
               ! to. The gradient below is assembled by whichever routine matches:
               ! `libcint_ri_mp2_gradient` differentiates the fitted correlation,
               ! `libcint_mp2_gradient` the exact one. Before the fitted one knew
               ! about functionals, a gradient run had to drop back to exact
               ! integrals to stay consistent with its own gradient -- which left
               ! an `Energy` run and a `Gradient` run of one deck reporting
               ! different energies, and left the expensive `n^4` two-particle
               ! density in the one place fitting was supposed to remove it.
               if (settings%aux_basis_named) then
                  call correlation_aux_basis(settings, fragment, symbols, corr_aux, error)
                  if (error%has_error()) then
                     call result%error%set(ERROR_VALIDATION, error%get_message())
                     result%has_error = .true.
                     call xc%destroy()
                     call mol%destroy()
                     return
                  end if
                  call run_libcint_ri_mp2(mol, corr_aux, scf%orbitals, &
                                          scf%orbital_energies, fragment%nelec/2, &
                                          scf%energy, dh_mp2, error, n_frozen=dh_frozen)
                  call corr_aux%destroy()
               else
                  call run_libcint_mp2(mol, scf%orbitals, scf%orbital_energies, &
                                       fragment%nelec/2, scf%energy, dh_mp2, error, &
                                       n_frozen=dh_frozen)
               end if
               if (error%has_error()) then
                  call result%error%set(ERROR_VALIDATION, "double hybrid: "// &
                                        error%get_message())
                  result%has_error = .true.
                  call xc%destroy()
                  call mol%destroy()
                  return
               end if

               ! Scaled here and stored already scaled, in the field `total` adds
               ! beside the SCF rather than in `mp2` -- see `energy_t`. Putting it in
               ! `mp2` would double count the moment a deck asked for MP2 as well.
               result%energy%dh_pt2 = xc%pt2_fraction*dh_mp2%correlation
               if (settings%verbose) then
                  write (line, "(a,f6.3,a,f20.12)") "  double hybrid: PT2 x", xc%pt2_fraction, " = ", result%energy%dh_pt2
                  call logger%info(trim(line))
                  write (line, "(a,f20.12)") "  E total        ", result%energy%total()
                  call logger%info(trim(line))
               end if

               ! The perturbative term's gradient, on top of the Kohn-Sham one
               ! already in `result%gradient`. Two calls rather than one because
               ! only the first half is variational: the reference gradient
               ! needs no response, and this one is almost entirely response.
               !
               ! **The gradient follows the energy**, on both axes. Fitted
               ! correlation is differentiated by `libcint_ri_mp2_gradient` and
               ! exact correlation by `libcint_mp2_gradient`; a fitted reference
               ! is differentiated as fitted by either. Every other branch in
               ! this file exists to keep those two in step, and this is the one
               ! that used to be out of step: the fitted correlation gradient
               ! knew nothing about functionals, so a double hybrid had to take
               ! the exact transform whatever the deck asked for.
               !
               ! What that cost was not accuracy but scale. The exact assembly
               ! forms a dense `n_ao^4` two-particle density, which is the
               ! ceiling fitting exists to remove, and it was still there in the
               ! one place a double hybrid needs it gone.
               if (do_gradient) then
                  block
                     real(dp), allocatable :: dh_grad(:, :)
                     type(libcint_molecule_t), target :: dh_aux
                     type(libcint_molecule_t), pointer :: dh_aux_arg

                     ! One auxiliary basis, serving whichever halves are fitted.
                     ! Rebuilt here for the reason the MP2 branch rebuilds its
                     ! own: the SCF's copy was destroyed above, and each branch
                     ! owning one means no error return in between has to
                     ! remember it.
                     dh_aux_arg => null()
                     if (settings%aux_basis_named .or. settings%density_fitting) then
                        call correlation_aux_basis(settings, fragment, symbols, dh_aux, error)
                        if (error%has_error()) then
                           call result%error%set(ERROR_VALIDATION, error%get_message())
                           result%has_error = .true.
                           call xc%destroy()
                           call mol%destroy()
                           return
                        end if
                        dh_aux_arg => dh_aux
                     end if
                     if (settings%aux_basis_named) then
                        call libcint_ri_mp2_gradient(mol, dh_aux, scf%orbitals, &
                                                     scf%orbital_energies, &
                                                     fragment%nelec/2, dh_grad, error, &
                                                     n_frozen=dh_frozen, &
                                                     fitted_reference=settings%density_fitting, &
                                                     xc=xc, scf_density=scf%density, &
                                                     pt2_scale=xc%pt2_fraction)
                     else
                        call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, &
                                                  fragment%nelec/2, dh_grad, error, &
                                                  xc=xc, scf_density=scf%density, &
                                                  pt2_scale=xc%pt2_fraction, aux=dh_aux_arg)
                     end if
                     if (associated(dh_aux_arg)) call dh_aux%destroy()
                     if (error%has_error()) then
                        call result%error%set(ERROR_VALIDATION, "double hybrid "// &
                                              "gradient: "//error%get_message())
                        result%has_error = .true.
                        call xc%destroy()
                        call mol%destroy()
                        return
                     end if
                     result%gradient = result%gradient + dh_grad
                     if (settings%verbose) then
                        write (line, "(a,f20.12)") "  |gradient|     ", &
                           sqrt(sum(result%gradient**2))
                        call logger%info(trim(line))
                     end if
                  end block
               end if
            end block
         end if
      end if

      ! The frontier pair, from the orbital energies the SCF already produced.
      if (allocated(scf%orbital_energies)) then
         if (scf%n_occupied >= 1 .and. scf%n_occupied < size(scf%orbital_energies)) then
            result%homo = scf%orbital_energies(scf%n_occupied)
            result%lumo = scf%orbital_energies(scf%n_occupied + 1)
            result%has_orbitals = .true.
         end if
      end if

      if (kohn_sham) call xc%destroy()
      call mol%destroy()
   end subroutine run_libcint_hf

   pure function angular_form_name(cartesian) result(name)
      !! "Cartesian" or "spherical", so an error can say which basis is which
      logical, intent(in) :: cartesian
      character(len=:), allocatable :: name

      if (cartesian) then
         name = "Cartesian"
      else
         name = "spherical"
      end if
   end function angular_form_name

   subroutine resolve_active_space(settings, fragment, n_ao, space, error)
      !! Turn the deck's active space into the four integers the CI needs
      !!
      !! Everything refusable about a CASSCF request is refusable here, before a
      !! single integral is computed, and that is deliberate: an active space is
      !! four small integers that either describe a valid CI problem or do not,
      !! and finding out that they do not after an SCF has run is a waste of the
      !! SCF and a worse error message.
      !!
      !! **The spin split comes from the molecular multiplicity, not from a
      !! keyword.** `n_active_electrons` says how many electrons the CI
      !! distributes; the multiplicity says how they are distributed between
      !! the spin channels. Every inactive orbital is doubly occupied and so
      !! contributes nothing to Ms, which means the whole of the molecule's
      !! excess alpha population has to sit in the active space -- so
      !! `n_alpha - n_beta = multiplicity - 1` exactly. A separate keyword for
      !! this would be a second place to say the same thing, and the two could
      !! disagree.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      integer, intent(in) :: n_ao
      integer, intent(out) :: space(4)
         !! n_inactive, n_active, n_alpha, n_beta -- in that order
      type(error_t), intent(inout) :: error

      integer :: n_active_electrons, n_active, n_inactive, n_alpha, n_beta
      integer :: unpaired, closed_shell_electrons

      space = 0
      n_active_electrons = settings%mcscf%n_active_electrons
      n_active = settings%mcscf%n_active_orbitals

      ! An unset active space is the one mistake worth naming in full, because
      ! it is what a first CASSCF deck gets wrong and there is no sensible
      ! default to fall back on: "all the valence orbitals" is a different
      ! number for every molecule, and guessing one would produce a converged
      ! energy for a calculation nobody asked for.
      if (n_active_electrons <= 0 .or. n_active <= 0) then
         call error%set(ERROR_VALIDATION, "a multiconfigurational method needs an "// &
                        "active space, and this deck has none. Set "// &
                        "keywords.mcscf.n_active_electrons and "// &
                        "keywords.mcscf.n_active_orbitals -- for example 6 and 6 for "// &
                        "the three bonding and three antibonding orbitals of a triple "// &
                        "bond. There is no default: the right active space is a "// &
                        "property of the chemistry, not of the molecule.")
         return
      end if

      unpaired = fragment%multiplicity - 1
      if (mod(n_active_electrons - unpaired, 2) /= 0) then
         call error%set(ERROR_VALIDATION, "an active space of "// &
                        to_text(n_active_electrons)//" electrons cannot have "// &
                        "multiplicity "//to_text(fragment%multiplicity)//": the two "// &
                        "differ in parity, so there is no way to split the electrons "// &
                        "into alpha and beta counts. Change either the electron count "// &
                        "or the multiplicity by one.")
         return
      end if
      n_alpha = (n_active_electrons + unpaired)/2
      n_beta = (n_active_electrons - unpaired)/2

      if (n_beta < 0 .or. n_alpha > n_active) then
         call error%set(ERROR_VALIDATION, "an active space of "//to_text(n_active)// &
                        " orbitals cannot hold "//to_text(n_alpha)//" alpha and "// &
                        to_text(n_beta)//" beta electrons at multiplicity "// &
                        to_text(fragment%multiplicity)//". Widen the active space or "// &
                        "put fewer electrons in it.")
         return
      end if

      ! Every electron the active space does not hold is in a doubly occupied
      ! orbital, so the inactive count follows from the electron count and does
      ! not normally have to be said. It is settable because a deck may want to
      ! freeze more orbitals than the arithmetic gives -- but then the electrons
      ! have to add up, which is what the second branch checks.
      closed_shell_electrons = fragment%nelec - n_active_electrons
      if (settings%mcscf%n_inactive_orbitals < 0) then
         if (closed_shell_electrons < 0 .or. mod(closed_shell_electrons, 2) /= 0) then
            call error%set(ERROR_VALIDATION, "cannot derive the inactive orbitals: "// &
                           to_text(fragment%nelec)//" electrons less an active space "// &
                           "of "//to_text(n_active_electrons)//" leaves "// &
                           to_text(closed_shell_electrons)//", which is not a "// &
                           "non-negative even number of doubly occupied electrons. "// &
                           "Set keywords.mcscf.n_inactive_orbitals if the partition is "// &
                           "meant to be something other than the obvious one.")
            return
         end if
         n_inactive = closed_shell_electrons/2
      else
         n_inactive = settings%mcscf%n_inactive_orbitals
         if (2*n_inactive + n_active_electrons /= fragment%nelec) then
            call error%set(ERROR_VALIDATION, "the orbital partition does not account "// &
                           "for every electron: "//to_text(n_inactive)//" inactive "// &
                           "orbitals hold "//to_text(2*n_inactive)//" electrons and "// &
                           "the active space "//to_text(n_active_electrons)//", which "// &
                           "is "//to_text(2*n_inactive + n_active_electrons)//" against "// &
                           "the molecule's "//to_text(fragment%nelec)//".")
            return
         end if
      end if

      ! The virtual space may be empty -- a CAS in a minimal basis legitimately
      ! uses every orbital -- but it may not be negative.
      if (n_inactive + n_active > n_ao) then
         call error%set(ERROR_VALIDATION, to_text(n_inactive)//" inactive plus "// &
                        to_text(n_active)//" active orbitals is more than the "// &
                        to_text(n_ao)//" the basis '"//trim(settings%basis_set)// &
                        "' provides. Use a larger basis or a smaller active space.")
         return
      end if

      space = [n_inactive, n_active, n_alpha, n_beta]
   end subroutine resolve_active_space

   subroutine run_libcint_mcscf(settings, fragment, result)
      !! CASSCF, or CASCI on the reference orbitals, for one fragment
      !!
      !! Three steps, and the middle one is the reason this is not just a call
      !! into `mqc_libcint_mcscf`: a closed-shell SCF supplies the orbitals the
      !! active space is carved out of. CASSCF then moves them and CASCI does
      !! not, which is the only difference between the two -- same wavefunction,
      !! same active space, same CI solver.
      !!
      !! **The reference is restricted, and an open-shell molecule is refused
      !! rather than run unrestricted.** The active-space integrals are built
      !! from one set of orbitals with an inactive count, exactly as the MP2
      !! transform is, so alpha and beta orbital sets have nowhere to go. That
      !! is a restriction on the *reference*, not on the state: a triplet with
      !! an even electron count is perfectly reachable here, because the
      !! multiplicity enters through the CI's alpha and beta string counts and
      !! not through the SCF.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casscf_result_t) :: casscf
      type(casci_result_t) :: casci
      type(error_t) :: error
      type(error_t) :: analysis_error
      type(avas_result_t) :: avas
      logical :: use_avas, use_valence
      real(dp), allocatable :: natural(:, :), occupations(:), reference(:, :)
      character(len=MAX_ELEMENT_SYMBOL_LEN), allocatable :: symbols(:)
      character(len=MAX_LINE_LENGTH) :: line
      integer :: iatom, diis_size
      integer :: space(4)

      if (settings%pcm%enabled) then
         call result%error%set(ERROR_VALIDATION, "continuum solvation (keywords.pcm) is "// &
                               "not implemented for multiconfigurational methods on the "// &
                               "CPU backend. Refused rather than run in the gas phase.")
         result%has_error = .true.
         return
      end if

      if (mod(fragment%nelec, 2) /= 0) then
         call result%error%set(ERROR_VALIDATION, "a multiconfigurational calculation "// &
                               "here starts from a closed-shell SCF, and "// &
                               to_text(fragment%nelec)//" electrons have no restricted "// &
                               "solution. An open-shell reference would need separate "// &
                               "alpha and beta transforms into the active space, which "// &
                               "are not implemented. An open-shell *state* on an "// &
                               "even-electron molecule is fine -- set the multiplicity "// &
                               "and the CI will find it.")
         result%has_error = .true.
         return
      end if

      allocate (symbols(fragment%n_atoms))
      do iatom = 1, fragment%n_atoms
         symbols(iatom) = element_number_to_symbol(fragment%element_numbers(iatom))
      end do

      call build_libcint_molecule(fragment%element_numbers, symbols, &
                                  fragment%coordinates, trim(settings%basis_set), &
                                  mol, error, ghost=ghost_of(fragment))
      if (error%has_error()) then
         call result%error%set(ERROR_VALIDATION, error%get_message())
         result%has_error = .true.
         return
      end if

      ! Before the SCF, so a malformed active space costs nothing to discover.
      ! Skipped when AVAS is choosing: the counts do not exist yet and cannot,
      ! since the selection needs converged orbitals to project.
      use_avas = allocated(settings%mcscf%avas_orbitals)
      use_valence = settings%mcscf%full_valence
      if (.not. use_avas .and. .not. use_valence) then
         call resolve_active_space(settings, fragment, mol%nao, space, error)
         if (error%has_error()) then
            call result%error%set(ERROR_VALIDATION, error%get_message())
            result%has_error = .true.
            call mol%destroy()
            return
         end if
      end if

      if (settings%verbose .and. .not. use_avas .and. .not. use_valence) then
         write (line, "(a,i0,a,i0,a,i0,a,i0,a,i0,a)") "  active space: CAS(", &
            settings%mcscf%n_active_electrons, ",", space(2), "), ", space(1), &
            " inactive, ", space(3), " alpha and ", space(4), " beta active electrons"
         call logger%info(trim(line))
      end if

      diis_size = settings%diis_size
      if (.not. settings%use_diis) diis_size = 0

      call run_libcint_rhf(mol, fragment%nelec, settings%max_iter, settings%energy_tol, &
                           settings%density_tol, settings%verbose, scf, error, &
                           diis_vectors=diis_size)
      if (error%has_error()) then
         call result%error%set(ERROR_VALIDATION, error%get_message())
         result%has_error = .true.
         call mol%destroy()
         return
      end if
      ! The reference is a starting point for CASSCF and the answer itself for
      ! CASCI, so a reference that did not converge is refused in both cases --
      ! `allow_crap_scf` is about reporting an unconverged SCF as a result, and
      ! nothing downstream of this one would say it had happened.
      if (.not. scf%converged) then
         call result%error%set(ERROR_VALIDATION, "the reference SCF did not converge in "// &
                               to_text(settings%max_iter)//" iterations, so there are no "// &
                               "orbitals to build an active space from")
         result%has_error = .true.
         call mol%destroy()
         return
      end if

      ! ---- AVAS, if the deck named orbitals rather than counts ---------------
      !
      ! Here rather than before the SCF because the projection needs converged
      ! orbitals: it asks which *molecular* orbitals carry the atomic character
      ! requested, and there are none to ask about until the SCF has run.
      allocate (reference(size(scf%orbitals, 1), size(scf%orbitals, 2)))
      reference = scf%orbitals
      if (use_valence) then
         ! The whole valence shell, sized by counting the free-atom minimal
         ! basis. Here rather than earlier for the same reason as AVAS: the
         ! valence-virtual orbitals are combinations of converged virtuals, and
         ! there are none to combine until the SCF has run.
         call valence_select(mol, fragment%element_numbers, symbols, &
                             fragment%coordinates, scf%orbitals, fragment%nelec, &
                             avas, error, verbose=settings%verbose)
         if (error%has_error()) then
            call result%error%set(ERROR_VALIDATION, error%get_message())
            result%has_error = .true.
            call mol%destroy()
            return
         end if
         reference = avas%orbitals
         space = [avas%n_inactive, avas%n_active, avas%n_active_electrons/2, &
                  avas%n_active_electrons/2]
      end if

      if (use_avas) then
         call avas_select(mol, fragment%element_numbers, symbols, fragment%coordinates, &
                          scf%orbitals, scf%n_occupied, settings%mcscf%avas_orbitals, &
                          avas, error, threshold=settings%mcscf%avas_threshold, &
                          verbose=settings%verbose)
         if (error%has_error()) then
            call result%error%set(ERROR_VALIDATION, error%get_message())
            result%has_error = .true.
            call mol%destroy()
            return
         end if
         reference = avas%orbitals
         space = [avas%n_inactive, avas%n_active, avas%n_active_electrons/2, &
                  avas%n_active_electrons/2]
         if (modulo(avas%n_active_electrons, 2) /= 0) then
            call result%error%set(ERROR_VALIDATION, "AVAS selected an odd number of "// &
                                  "active electrons, which a closed-shell reference "// &
                                  "cannot distribute evenly over two spins.")
            result%has_error = .true.
            call mol%destroy()
            return
         end if
      end if

      if (settings%mcscf%optimize_orbitals .and. allocated(settings%mcscf%ormas_subspaces)) then
         call run_libcint_casscf(mol, reference, space(1), space(2), space(3), space(4), &
                                 casscf, error, &
                                 max_iterations=settings%mcscf%max_macro_iter, &
                                 gradient_tol=settings%mcscf%orbital_convergence, &
                                 verbose=settings%verbose, &
                                 subspaces=settings%mcscf%ormas_subspaces, &
                                 min_electrons=settings%mcscf%ormas_min_electrons, &
                                 max_electrons=settings%mcscf%ormas_max_electrons)
      else if (settings%mcscf%optimize_orbitals) then
         call run_libcint_casscf(mol, reference, space(1), space(2), space(3), space(4), &
                                 casscf, error, &
                                 max_iterations=settings%mcscf%max_macro_iter, &
                                 gradient_tol=settings%mcscf%orbital_convergence, &
                                 verbose=settings%verbose)
         if (.not. error%has_error() .and. .not. casscf%converged) then
            if (casscf%stalled) then
               ! Distinguished from running out of iterations, because the
               ! advice is the opposite. A stall means no step downhill could be
               ! found at all: the approximate Hessian has run out of
               ! resolution, and more iterations will do nothing. It happens on
               ! flat surfaces, where the energy is converged long before the
               ! gradient is small.
               call error%set(ERROR_VALIDATION, "the orbital optimisation stopped "// &
                              "improving at a gradient of "// &
                              to_real_text(casscf%gradient_norm)//", short of the "// &
                              to_real_text(settings%mcscf%orbital_convergence)// &
                              " asked for, after "//to_text(casscf%iterations)// &
                              " macro-iterations. More iterations will not help -- "// &
                              "no step downhill could be found. Ask for a looser "// &
                              "keywords.mcscf.orbital_convergence; on a flat surface "// &
                              "the energy is converged well before the gradient is.")
            else
               call error%set(ERROR_VALIDATION, "the orbital optimisation did not reach an "// &
                              "orbital gradient of "// &
                              to_real_text(settings%mcscf%orbital_convergence)//" in "// &
                              to_text(settings%mcscf%max_macro_iter)//" macro-iterations "// &
                              "(it stopped at "//to_real_text(casscf%gradient_norm)// &
                              "). Raise keywords.mcscf.max_macro_iter.")
            end if
         end if
         if (error%has_error()) then
            call result%error%set(ERROR_VALIDATION, error%get_message())
            result%has_error = .true.
            call mol%destroy()
            return
         end if
         ! A CASSCF total is a complete energy, not a correction to a reference,
         ! so it goes in the slot `energy_total` adds unconditionally. There is
         ! no separate multiconfigurational component to report beside it: the
         ! "correlation" a CAS recovers is not separable from the reference the
         ! way an MP2 or a coupled-cluster correction is, because the orbitals
         ! underneath it are no longer the Hartree-Fock ones.
         result%energy%scf = casscf%energy
         if (settings%verbose) then
            write (line, "(a,f20.12)") "  E(CASSCF)      ", casscf%energy
            call logger%info(trim(line))
         end if
      else
         ! Not `orbital_convergence`, which is a gradient threshold and means
         ! nothing to a Davidson, and not a keyword either: the CI is the whole
         ! answer on this path, so there is no reason to solve it loosely. The
         ! same threshold the CASSCF macro loop pins its own CASCI at.
         if (allocated(settings%mcscf%ormas_subspaces)) then
            call run_libcint_ormas_ci(mol, reference, space(1), space(2), space(3), &
                                      space(4), settings%mcscf%ormas_subspaces, &
                                      settings%mcscf%ormas_min_electrons, &
                                      settings%mcscf%ormas_max_electrons, casci, error, &
                                      verbose=settings%verbose, tolerance=CI_TOLERANCE)
         else
            call run_libcint_casci(mol, reference, space(1), space(2), space(3), space(4), &
                                   casci, error, verbose=settings%verbose, &
                                   tolerance=CI_TOLERANCE)
         end if
         if (error%has_error()) then
            call result%error%set(ERROR_VALIDATION, error%get_message())
            result%has_error = .true.
            call mol%destroy()
            return
         end if
         result%energy%scf = casci%energy
         if (settings%verbose) then
            if (allocated(settings%mcscf%ormas_subspaces)) then
               write (line, "(a,f20.12)") "  E(ORMAS-CI)    ", casci%energy
            else
               write (line, "(a,f20.12)") "  E(CASCI)       ", casci%energy
            end if
            call logger%info(trim(line))
         end if
      end if

      ! ---- analyses on the correlated wave function --------------------------
      !
      ! The natural orbitals rather than the optimised ones, and the occupations
      ! with them. An MCSCF has no "occupied orbitals" the way a reference
      ! determinant does -- the active ones carry fractional occupation and the
      ! optimised set is not ordered by it -- so the analysis is given the basis
      ! in which the density is diagonal and told what the diagonal is.
      if (bonding_analysis_kind(settings%bonding_analysis) == BONDING_GMS_QUAO) then
         call analysis_error%clear()
         if (settings%mcscf%optimize_orbitals) then
            call natural_orbitals(casscf%orbitals, space(1), space(2), casscf%dm1, &
                                  natural, occupations, analysis_error)
         else
            call natural_orbitals(reference, space(1), space(2), casci%dm1, &
                                  natural, occupations, analysis_error)
         end if
         if (.not. analysis_error%has_error()) then
            ! The two-particle density goes over in the orbitals it was computed
            ! in, not the natural ones the rest of the analysis uses -- it
            ! belongs to that basis and is projected out of it. Its energy is
            ! handed over too, since the decomposition's own reference assumes a
            ! single determinant and would be short by the correlation energy.
            if (settings%mcscf%optimize_orbitals) then
               call run_quao_analysis(mol, fragment%element_numbers, symbols, &
                                      fragment%coordinates, natural, fragment%nelec, &
                                      analysis_error, &
                                      threshold=settings%bonding_threshold, &
                                      occupations=occupations, &
                                      active_orbitals=casscf%orbitals(:, space(1) + 1: &
                                                                      space(1) + space(2)), &
                                      active_dm1=casscf%dm1, active_dm2=casscf%dm2, &
                                      reference_energy=casscf%energy, &
                                      no_sharing=settings%bonding_no_sharing)
            else
               call run_quao_analysis(mol, fragment%element_numbers, symbols, &
                                      fragment%coordinates, natural, fragment%nelec, &
                                      analysis_error, &
                                      threshold=settings%bonding_threshold, &
                                      occupations=occupations, &
                                      active_orbitals=reference(:, space(1) + 1: &
                                                                space(1) + space(2)), &
                                      active_dm1=casci%dm1, active_dm2=casci%dm2, &
                                      reference_energy=casci%energy, &
                                      no_sharing=settings%bonding_no_sharing)
            end if
         end if
         if (analysis_error%has_error()) then
            call logger%warning("  the bonding analysis could not run: "// &
                                analysis_error%get_message())
         end if
      end if

      result%scf_status = SCF_CONVERGED
      result%has_energy = .true.
      call mol%destroy()
   end subroutine run_libcint_mcscf

   pure function to_real_text(value) result(out)
      !! A small real in exponent form, for an error message
      real(dp), intent(in) :: value
      character(len=:), allocatable :: out
      character(len=24) :: buffer
      write (buffer, "(es9.2)") value
      out = trim(adjustl(buffer))
   end function to_real_text

   pure function to_text(value) result(out)
      integer, intent(in) :: value
      character(len=:), allocatable :: out
      character(len=16) :: buffer
      write (buffer, "(i0)") value
      out = trim(adjustl(buffer))
   end function to_text

   subroutine correlation_aux_basis(settings, fragment, symbols, aux, error)
      !! Build the auxiliary basis the correlation step will fit with
      !!
      !! Whatever the deck names is used, including sets that have no business
      !! fitting a (ia|jb) block -- a JKFIT set is fitted for the Coulomb and
      !! exchange matrices and will give a correlation energy whose error is
      !! not the RI error it is supposed to be. That is a warning rather than a
      !! refusal: someone comparing against another program's default, or
      !! probing how bad it gets, has a real reason to ask for it, and the run
      !! is not wrong so much as poorly fitted.
      !!
      !! `model.aux_basis` is the only place an auxiliary basis is named. It used
      !! to be settable under `keywords.correlation` as well, which meant two
      !! places for one thing and a silent preference between them; a basis set
      !! belongs beside the orbital basis it fits.
      !!
      !! So a run that fits both the reference and the correlation uses one set for
      !! both, which is fine in that direction: a RIFIT set fitting J and K is
      !! ordinary practice. It is worth about 1.7 mHartree on a total energy here
      !! against exact J and K -- measured, and largely cancelling in any relative
      !! quantity, which is why it is a reasonable thing to do rather than a
      !! compromise.
      !!
      !! The other direction is the one to catch, and the warning below is for it: a
      !! JKFIT set fitting a (ia|jb) block gives a correlation energy whose error is
      !! not the RI error it is supposed to be.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      character(len=*), intent(in) :: symbols(:)
      type(libcint_molecule_t), intent(out) :: aux
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: name

      name = trim(settings%aux_basis_set)
      if (len_trim(name) == 0) then
         call error%set(ERROR_VALIDATION, "density-fitted correlation needs an "// &
                        "auxiliary basis: set model.aux_basis")
         return
      end if

      if (index(name, "rifit") == 0 .and. index(name, "-ri") == 0) then
         call logger%warning("bad basis set, calculations will be poor: '"//name// &
                             "' is not a correlation-fitting (RIFIT) set")
      end if

      call build_libcint_molecule(fragment%element_numbers, symbols, &
                                  fragment%coordinates, name, aux, error, &
                                  ghost=ghost_of(fragment))
   end subroutine correlation_aux_basis

   pure function core_orbital_count(atomic_numbers) result(n_core)
      !! How many orbitals a frozen core leaves out, summed over the atoms
      !!
      !! The count per element is the number of filled shells below the valence
      !! one: none for H and He, the 1s for Li through Ne, and so on. This is
      !! the same convention PySCF and most others use by default, which is the
      !! point -- an energy computed with a different core is not comparable to
      !! a published one, and the difference is millihartrees rather than
      !! anything that looks like a bug.
      integer, intent(in) :: atomic_numbers(:)
      integer :: n_core

      integer :: i, z

      n_core = 0
      do i = 1, size(atomic_numbers)
         z = atomic_numbers(i)
         if (z <= 10) then
            if (z > 2) n_core = n_core + 1   ! 1s, and nothing at all for H and He
         else if (z <= 18) then
            n_core = n_core + 5        ! 1s 2s 2p
         else if (z <= 36) then
            n_core = n_core + 9        ! + 3s 3p
         else if (z <= 54) then
            n_core = n_core + 18       ! + 3d 4s 4p
         else
            n_core = n_core + 27       ! + 4d 5s 5p
         end if
      end do
   end function core_orbital_count

   pure function ghost_of(fragment) result(ghost)
      !! A fragment's ghost mask, or all-false when it has none
      !!
      !! `is_ghost` is unallocated on every fragment an ordinary expansion
      !! builds, and passing an unallocated array to an optional dummy is an
      !! absent argument -- which would be fine, except three call sites would
      !! each need the same conditional. This says it once.
      type(physical_fragment_t), intent(in) :: fragment
      logical :: ghost(fragment%n_atoms)

      if (allocated(fragment%is_ghost)) then
         ghost = fragment%is_ghost
      else
         ghost = .false.
      end if
   end function ghost_of

end module mqc_libcint_bridge
