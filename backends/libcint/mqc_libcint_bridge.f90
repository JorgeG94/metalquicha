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
   use mqc_program_limits, only: MAX_LINE_LENGTH
   implicit none
   private

   public :: run_libcint_hf
   public :: run_libcint_fmo
   public :: run_libcint_makefp
   public :: run_libcint_charges
   public :: run_libcint_efp
   public :: run_libcint_sapt0
   public :: libcint_backend_available

contains

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
      type(xc_context_t), target :: xc
      type(xc_context_t), pointer :: xc_arg
      logical :: kohn_sham

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
         if (settings%density_fitting) then
            call result%error%set(ERROR_VALIDATION, "the MP2 gradient differentiates an "// &
                                  "exact-ERI reference: with keywords.scf.density_fitting "// &
                                  "on it would be the derivative of an energy nothing "// &
                                  "computed. Fitting the correlation instead -- an "// &
                                  "`ri-mp2` method, which fits nothing in the SCF -- is "// &
                                  "implemented and differentiated exactly.")
            result%has_error = .true.
            return
         end if
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

      call build_libcint_molecule(fragment%element_numbers, symbols, &
                                  fragment%coordinates, trim(settings%basis_set), &
                                  mol, error)
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
         call build_libcint_molecule(fragment%element_numbers, symbols, &
                                     fragment%coordinates, trim(settings%aux_basis_set), &
                                     aux, error)
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
         call run_libcint_rhf(mol, fragment%nelec, settings%max_iter, settings%energy_tol, &
                              settings%density_tol, settings%verbose, scf, error, &
                              diis_vectors=diis_size, guess=guess_kind, &
                              guess_density=guess_total, xc=xc)
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

      ! ---- dE/dR, analytically ----------------------------------------------
      !
      ! Differentiating the energy that was just converged, which is why the
      ! auxiliary basis handed over is the one the SCF fitted with rather than
      ! whatever the deck names: a fitted energy differentiated as an exact one
      ! is a gradient of a function nobody computed.
      if (do_gradient .and. .not. settings%run_mp2) then
         ! A double hybrid's energy is the Kohn-Sham part plus a scaled PT2
         ! correlation, and only the first of those is differentiated below.
         ! Left unchecked this returns the hybrid's gradient against the double
         ! hybrid's energy -- on water/cc-pVDZ that is 0.011 against a true
         ! 0.016, right in shape, wrong by a third, and with nothing in the
         ! output to say which of the two numbers the other belongs to.
         !
         ! The PT2 piece needs a relaxed density and a response solve over a
         ! Kohn-Sham reference, which means the exchange-correlation kernel in
         ! the coupled-perturbed operator. That kernel does not exist here yet.
         if (kohn_sham) then
            if (xc%pt2_fraction /= 0.0_dp) then
               call result%error%set(ERROR_VALIDATION, "a double hybrid gradient needs "// &
                                     "the PT2 term differentiated too: its relaxed "// &
                                     "density and response solve run over a Kohn-Sham "// &
                                     "reference, whose exchange-correlation kernel is "// &
                                     "not implemented. Run the double hybrid as an "// &
                                     "energy, or take the gradient at a hybrid.")
               result%has_error = .true.
               call mol%destroy()
               return
            end if
         end if

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
         if (settings%verbose) then
            write (*, "(a,f20.12)") "  |gradient|     ", sqrt(sum(result%gradient**2))
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
                  call libcint_ri_mp2_gradient(mol, corr_aux, scf%orbitals, &
                                               scf%orbital_energies, fragment%nelec/2, &
                                               result%gradient, error, n_frozen=frozen)
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
               if (settings%verbose) then
                  write (*, "(a,f20.12)") "  |gradient|     ", &
                     sqrt(sum(result%gradient**2))
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
            integer :: frozen

            ! The same frozen-core rule the MP2 block applies, deliberately
            ! duplicated rather than hoisted: the two blocks are independent, and
            ! a shared local would silently couple them if one ever wanted a
            ! different count.
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
               call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, &
                                     fragment%nelec/2, frozen, settings%cc_max_iter, &
                                     settings%cc_tolerance, settings%cc_triples, &
                                     settings%verbose, cc, error, &
                                     diis_vectors=settings%cc_diis_size, aux=corr_aux)
               call corr_aux%destroy()
            else
               call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, &
                                     fragment%nelec/2, frozen, settings%cc_max_iter, &
                                     settings%cc_tolerance, settings%cc_triples, &
                                     settings%verbose, cc, error, &
                                     diis_vectors=settings%cc_diis_size)
            end if
            if (error%has_error()) then
               call result%error%set(ERROR_VALIDATION, "CCSD: "//error%get_message())
               result%has_error = .true.
               call mol%destroy()
               return
            end if

            ! energy_t sums scf + mp2%total() + cc%total(), so only the components
            ! go in. All three are filled rather than a lumped correlation energy:
            ! a total cannot be taken apart afterwards, and the singles/doubles
            ! split is what says whether the T1 amplitudes are doing anything.
            result%energy%cc%singles = cc%e_singles
            result%energy%cc%doubles = cc%e_doubles
            result%energy%cc%triples = cc%e_triples
            if (settings%verbose) then
               write (line, "(a,i0,a,l1)") "  CCSD: iterations ", cc%iterations, "  converged ", cc%converged
               call logger%info(trim(line))
               write (line, "(a,f20.12)") "  E(MP2)         ", cc%e_mp2
               call logger%info(trim(line))
               write (line, "(a,f20.12)") "  E(singles)     ", cc%e_singles
               call logger%info(trim(line))
               write (line, "(a,f20.12)") "  E(doubles)     ", cc%e_doubles
               call logger%info(trim(line))
               if (settings%cc_triples) then
                  write (line, "(a,f20.12)") "  E(T)           ", cc%e_triples
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

               ! Fitted when an auxiliary basis is available, exact otherwise. The
               ! reference implementations of these functionals are density-fitted,
               ! so fitting is the comparable choice rather than a compromise -- but
               ! a deck that named no auxiliary basis should get an answer, not a
               ! refusal, since the conventional transform is the more accurate of
               ! the two anyway.
               if (len_trim(settings%aux_basis_set) > 0) then
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
                                  fragment%coordinates, name, aux, error)
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

end module mqc_libcint_bridge
