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
   use mqc_string_utils, only: int_to_text
   use pic_types, only: dp
   use pic_timer, only: timer_type
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t, SCF_CONVERGED, SCF_NOT_CONVERGED, &
                               scf_not_converged_message
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_elements, only: element_number_to_symbol
   use mqc_program_limits, only: MAX_ELEMENT_SYMBOL_LEN, MAX_LINE_LENGTH, &
                                 ERI_CORE_BUDGET_CAP, ERI_CORE_BUDGET_SHARE, &
                                 ERI_CORE_BUDGET_BLIND, SAPT_CORE_BUDGET_SHARE
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_libcint_ecp, only: ecp_refuses_auto_frozen_core
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule, &
                                    angular_form_name
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf, run_libcint_uhf, &
                              SCF_GUESS_GWH, SCF_GUESS_SAC, SCF_GUESS_SAD, SCF_GUESS_PROJ
   use mqc_libcint_projection, only: climb_basis_ladder
   use mqc_libcint_atomic_guess, only: build_atomic_guess, parse_guess_name, &
                                       guess_display_name
   use mqc_libcint_xc, only: xc_context_t, xc_context_create, xc_available
   use mqc_libcint_pcm, only: pcm_context_t
   use mqc_libcint_gradient, only: libcint_scf_gradient
   use mqc_libcint_hessian, only: rhf_hessian, ks_hessian, hessian_to_matrix, &
                                  nuclear_repulsion_hessian, response_hessian
   use mqc_libcint_mp2_hessian, only: mp2_correlation_hessian
   use mqc_libcint_mp2_gradient, only: libcint_mp2_gradient
   use mqc_libcint_ri_mp2_gradient, only: libcint_ri_mp2_gradient
   use mqc_libcint_mp2, only: mp2_result_t, run_libcint_mp2, run_libcint_ri_mp2, &
                              run_libcint_ump2, run_libcint_uri_mp2
   use mqc_libcint_cc, only: cc_result_t, run_libcint_ccsd
   use mqc_libcint_rcc, only: rcc_result_t, run_libcint_rccsd
   use mqc_libcint_bonding, only: valence_wavefunction_t, &
                                  run_quao_analysis, bonding_analysis_kind, &
                                  BONDING_GMS_QUAO
   use mqc_libcint_casci, only: casci_result_t, run_libcint_casci, &
                                run_libcint_ormas_ci
   use mqc_libcint_avas, only: avas_select, avas_result_t, valence_select
   use mqc_libcint_mcscf, only: casscf_result_t, run_libcint_casscf, &
                                natural_orbitals
   use mqc_libcint_mcscf_gradient, only: libcint_mcscf_gradient
   implicit none
   private

   public :: run_libcint_hf
   public :: run_libcint_mcscf
   public :: run_libcint_fmo
   public :: run_libcint_makefp
   public :: run_libcint_charges
   public :: run_libcint_efp
   public :: run_libcint_sapt0
   public :: run_libcint_sapt2
   public :: libcint_backend_available
   public :: xc_available
      !! Re-exported so a caller that cannot see `mqc_libcint_xc` -- anything
      !! outside this backend, since the module is not compiled without it --
      !! can still ask whether a functional can be evaluated. `mqc_version`
      !! reports it, and `run_validation.py` skips the decks that need it.

   real(dp), parameter :: CI_TOLERANCE = 1.0e-11_dp
      !! Residual the CASCI Davidson is driven to on the CASCI-only path.
      !!
      !! The same number `mqc_libcint_mcscf` pins its own inner CASCI at, and
      !! deliberately not a keyword. A CASCI *is* its CI energy -- there is
      !! nothing downstream to absorb a loose solve the way a CASSCF macro-step
      !! absorbs one -- so the only defensible setting is "tight", and 1e-11 on
      !! a hundred-hartree total is inside what a threaded Fock build
      !! reproduces.

contains

   function sapt_core_bytes(nao, want_sapt2) result(bytes)
      !! Roughly what the SAPT caches will ask for at their peak
      !!
      !! The dimer `eri` is the full `nao**4` -- no eightfold folding, because
      !! every SAPT term addresses it as four indices -- and `eri_packed` is
      !! another `n_pair**2` beside it, where `n_pair = nao(nao+1)/2`. That is
      !! about `1.25 nao**4` before a term is evaluated.
      !!
      !! SAPT2 adds one more `n_pair**2`: `build_sapt2_cache` copies the packed
      !! matrix to diagonalize it, and the copy is live while the original
      !! still is. The three-index factors that come out of it are smaller and
      !! not counted, so this errs low -- it is a floor on the requirement, not
      !! an estimate of it.
      !!
      !! `real` throughout: `nao**4` at four hundred functions overflows a
      !! 32-bit integer, and the symptom would be a large basis quietly
      !! deciding it fitted.
      integer, intent(in) :: nao
      logical, intent(in) :: want_sapt2
      real(dp) :: bytes
      real(dp) :: n, n_pair

      n = real(nao, dp)
      n_pair = n*(n + 1.0_dp)/2.0_dp
      bytes = 8.0_dp*(n**4 + n_pair**2)
      if (want_sapt2) bytes = bytes + 8.0_dp*n_pair**2
   end function sapt_core_bytes

   subroutine check_sapt_fits_in_core(nao, want_sapt2, error)
      !! Refuse a SAPT run whose stored integrals cannot fit in memory
      !!
      !! Where memory cannot be read -- anything that is not Linux -- this says
      !! nothing rather than guessing. A wrong refusal here is worse than no
      !! refusal: the calculation is possible on the machine and the message
      !! would claim otherwise.
      integer, intent(in) :: nao
      logical, intent(in) :: want_sapt2
      type(error_t), intent(inout) :: error
      real(dp) :: needed, available

      available = available_memory_bytes()
      if (available <= 0.0_dp) return

      needed = sapt_core_bytes(nao, want_sapt2)
      if (needed <= SAPT_CORE_BUDGET_SHARE*available) return

      call error%set(ERROR_VALIDATION, &
                     "sapt: the dimer basis has "//int_to_text(nao)// &
                     " functions, whose stored integrals need at least "// &
                     int_to_text(nint(needed/1.0e6_dp))//" MB, and this "// &
                     "machine has "//int_to_text(nint(available/1.0e6_dp))// &
                     " MB available. Every SAPT term contracts over the full "// &
                     "dimer tensor, so there is no direct fallback to drop to "// &
                     "-- use a smaller basis.")
   end subroutine check_sapt_fits_in_core

   pure function reuse_scf_fit(settings) result(reuse)
      !! Whether the SCF should hand its fitted tensor to the correlated step
      !!
      !! A fitted SCF followed by a fitted MP2 was building the same three-
      !! centre integrals twice -- four hundred million of them at 560
      !! functions. It need not build them even once more: the MO transform
      !! acts on the pair index and the metric on the auxiliary one, so they
      !! commute and the SCF's *fitted* tensor transforms into exactly the
      !! tensor the correlation step would have built. See `build_df_mo_block`.
      !!
      !! Both halves have to be fitted for the question to arise, and they are
      !! then necessarily the same auxiliary basis: `model.aux_basis` is the
      !! only place a fitting set is named, so there is no deck that fits both
      !! halves with different ones. If that ever stops being true this has to
      !! start comparing the two basis names.
      !!
      !! Restricted only, and MP2 only, because those are the paths wired to
      !! accept it. Carrying it further is plumbing rather than a new idea.
      type(cuest_scf_settings_t), intent(in) :: settings
      logical :: reuse

      reuse = settings%density_fitting .and. settings%corr_density_fitting &
              .and. settings%run_mp2 .and. .not. settings%unrestricted
   end function reuse_scf_fit

   function eri_fits_in_core(nao) result(fits)
      !! Whether n^4 stored integrals fit in this rank's share of memory
      !!
      !! Computed in `real` rather than integer arithmetic: n^4 at four hundred
      !! functions overflows a 32-bit integer, and the symptom of that would be
      !! a large basis quietly deciding it fitted.
      !!
      !! Not `pure` any more: it reads the machine.
      integer, intent(in) :: nao
      logical :: fits
      real(dp) :: available, budget

      available = available_memory_bytes()
      if (available > 0.0_dp) then
         budget = min(ERI_CORE_BUDGET_CAP, ERI_CORE_BUDGET_SHARE*available)
      else
         budget = ERI_CORE_BUDGET_BLIND
      end if

      fits = 8.0_dp*real(nao, dp)**4 <= budget
   end function eri_fits_in_core

   function available_memory_bytes() result(bytes)
      !! MemAvailable from /proc/meminfo, or zero where that does not exist
      !!
      !! MemAvailable rather than MemFree: free memory on a warm machine is
      !! almost nothing, because the kernel has spent it on page cache it will
      !! hand back on demand. MemFree would refuse to store integrals on a
      !! machine with plenty to spare.
      real(dp) :: bytes
      integer :: unit, stat
      character(len=MAX_LINE_LENGTH) :: line
      real(dp) :: kb

      bytes = 0.0_dp
      open (newunit=unit, file="/proc/meminfo", status="old", action="read", iostat=stat)
      if (stat /= 0) return
      do
         read (unit, "(a)", iostat=stat) line
         if (stat /= 0) exit
         if (line(1:13) == "MemAvailable:") then
            read (line(14:), *, iostat=stat) kb
            if (stat == 0) bytes = kb*1024.0_dp
            exit
         end if
      end do
      close (unit)
   end function available_memory_bytes

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
                                charge_a, charge_b, terms, error)
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
      integer, intent(in) :: charge_a, charge_b
         !! The monomers' own charges. Required rather than optional: a caller
         !! that forgets them gets a neutral monomer and a wrong number, which
         !! is the one failure mode worth making impossible to write.
      real(dp), intent(out) :: terms(N_SAPT_TERMS)
         !! Ordered by `SAPT_TERM_NAMES`
      type(error_t), intent(inout) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_terms_t) :: t

      terms = 0.0_dp
      call build_sapt_molecules(z_a, sym_a, xyz_a, z_b, sym_b, xyz_b, &
                                basis_name, mols, error, &
                                charge_a=charge_a, charge_b=charge_b)
      if (error%has_error()) return
      call check_sapt_fits_in_core(mols%dimer%nao, .false., error)
      if (error%has_error()) then
         call mols%destroy()
         return
      end if
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

   subroutine run_libcint_sapt2(z_a, sym_a, xyz_a, z_b, sym_b, xyz_b, basis_name, &
                                charge_a, charge_b, terms, error)
      !! SAPT2 between two monomers: every SAPT0 term in the same slots, plus
      !! the four intramonomer-correlation corrections, the scaled
      !! exchange-induction, and its own total
      use pic_types, only: dp
      use mqc_program_limits, only: N_SAPT2_TERMS
      use mqc_sapt, only: sapt_molecules_t, build_sapt_molecules, sapt_terms_t, &
                          run_sapt2
      integer, intent(in) :: z_a(:), z_b(:)
      character(len=*), intent(in) :: sym_a(:), sym_b(:)
      real(dp), intent(in) :: xyz_a(:, :), xyz_b(:, :)   !! (3, n), Bohr
      character(len=*), intent(in) :: basis_name
      integer, intent(in) :: charge_a, charge_b   !! The monomers' own charges
      real(dp), intent(out) :: terms(N_SAPT2_TERMS)
         !! Ordered by `SAPT2_TERM_NAMES`
      type(error_t), intent(inout) :: error

      type(sapt_molecules_t) :: mols
      type(sapt_terms_t) :: t

      terms = 0.0_dp
      call build_sapt_molecules(z_a, sym_a, xyz_a, z_b, sym_b, xyz_b, &
                                basis_name, mols, error, &
                                charge_a=charge_a, charge_b=charge_b)
      if (error%has_error()) return
      call check_sapt_fits_in_core(mols%dimer%nao, .true., error)
      if (error%has_error()) then
         call mols%destroy()
         return
      end if
      call run_sapt2(mols, t, error)
      if (error%has_error()) then
         call mols%destroy()
         return
      end if

      terms = [t%elst10, t%exch10_s2, t%exch10, t%ind20_u, t%ind20_r, &
               t%exch_ind20_u, t%exch_ind20_r, t%disp20, t%exch_disp20, &
               t%delta_hf, t%e_int_hf, t%total, &
               t%elst12, t%exch11, t%exch12, t%ind22, t%exch_ind22, &
               t%total_sapt2]
      call mols%destroy()
   end subroutine run_libcint_sapt2

   subroutine run_libcint_makefp(atomic_numbers, element_symbols, coordinates, &
                                 basis_name, name, path, error, charge, verbose, &
                                 aux_basis, guess, energy_tol, density_tol, &
                                 vdwscl, dynamic_tol, dynamic_maxiter, response, &
                                 allow_crap_response, response_batch)
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
      real(dp), intent(in), optional :: energy_tol, density_tol
         !! SCF thresholds, present only when the deck named them. Absent leaves
         !! `make_efp_potential` on the tighter pair a fragment potential needs.
      real(dp), intent(in), optional :: vdwscl
      real(dp), intent(in), optional :: dynamic_tol
      integer, intent(in), optional :: dynamic_maxiter
      integer, intent(in), optional :: response
      logical, intent(in), optional :: allow_crap_response
         !! Accept an unconverged response instead of refusing. The potential
         !! is wrong; see `efp_config_t`.
      integer, intent(in), optional :: response_batch
         !! Densities per integral pass. Tuning only; the answer is unchanged.
         !! Accept an unconverged response instead of refusing. See
         !! `efp_config_t`; the potential is wrong when this is on.
         !! The `keywords.efp` group: the screening grid's van der Waals scale, and
         !! the tolerance, iteration cap and route of the dynamic response solve.
         !! Forwarded and not read here, the same as everything else in this list.

      type(efp_potential_t) :: pot

      ! One call, not one per combination of present arguments. An absent optional
      ! dummy passed on as an actual argument arrives absent at the other end, which
      ! is what `charge`, `verbose` and `guess` were already relying on; `aux_basis`
      ! was branched around for no reason it needed, and four optionals would have
      ! meant sixteen branches to keep saying the same thing.
      call make_efp_potential(atomic_numbers, element_symbols, coordinates, &
                              basis_name, name, pot, error, charge=charge, &
                              verbose=verbose, aux_basis=aux_basis, guess=guess, &
                              energy_tol=energy_tol, density_tol=density_tol, &
                              vdwscl=vdwscl, dynamic_tol=dynamic_tol, &
                            dynamic_maxiter=dynamic_maxiter, response=response, allow_crap_response=allow_crap_response, &
                              response_batch=response_batch)
      if (error%has_error()) return
      call write_efp_potential(pot, path, error)
      call pot%destroy()
   end subroutine run_libcint_makefp

   subroutine run_libcint_fmo(atomic_numbers, element_symbols, coordinates, owner, &
                              basis_name, esp, expansion, far_field, resppc, &
                              level, max_outer, outer_tol, scf_max_iter, &
                              scf_energy_tol, scf_density_tol, bond_breaking, &
                              cap_scale, energy, error, comm)
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
      character(len=*), intent(in) :: bond_breaking
      real(dp), intent(in) :: cap_scale
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
      opts%bond_breaking = bond_breaking
      opts%cap_scale = cap_scale

      call run_fmo2(atomic_numbers, symbols, coordinates, owner, opts, res, error, comm)
      if (error%has_error()) return
      energy = res%energy
   end subroutine run_libcint_fmo

   subroutine run_libcint_hf(settings, fragment, result, want_gradient, want_hessian)
      !! Closed-shell HF for one fragment, on the CPU
      use mqc_libcint_charges, only: mulliken_charges, chelpg_charges, &
                                     mulliken_spin_populations
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient
      logical, intent(in), optional :: want_hessian
         !! Ask for the analytic Hessian. **Not a promise that one comes back.**
         !! The conditions it needs are known here and not by the caller, so a
         !! request this routine cannot honour leaves `result%has_hessian`
         !! false and sets no error -- the caller is expected to fall back to
         !! finite differences, which is a correct answer rather than a failure.

      ! `aux` and `xc` are targets because the gradient takes both as optional
      ! arguments, and a disassociated pointer is how "this SCF fitted nothing"
      ! and "this SCF had no functional" are said without four spellings of the
      ! same call.
      type(libcint_molecule_t) :: mol, corr_aux
      type(libcint_molecule_t), target :: aux
      type(libcint_molecule_t), pointer :: aux_arg
      type(rhf_result_t) :: scf
      type(pcm_context_t) :: pcm_ctx
      type(error_t) :: error
      character(len=MAX_ELEMENT_SYMBOL_LEN), allocatable :: symbols(:)
      integer :: iatom, diis_size, guess_kind
      integer :: nelec        !! Valence electrons: fragment%nelec less any ECP core
      logical :: unrestricted, do_gradient, do_hessian
      ! Left unallocated for the guesses that need no free-atom solutions. An
      ! unallocated allocatable passed to an optional dummy is an absent
      ! argument, so the SCF calls below need no branching on which guess ran.
      real(dp), allocatable :: guess_a(:, :), guess_b(:, :), guess_total(:, :)
      type(error_t) :: guess_error
      type(error_t) :: analysis_error
      real(dp), allocatable :: ieda_atom(:), ieda_free(:)
      real(dp), allocatable :: ieda_pair(:, :), ieda_classical(:, :)
      real(dp) :: ieda_formation
      type(xc_context_t), target :: xc
      type(xc_context_t), pointer :: xc_arg
      logical :: kohn_sham
      type(timer_type) :: grad_clock
      real(dp), allocatable :: scf_b_ao(:, :)
         !! The SCF's fitted tensor, taken over rather than freed when a fitted
         !! correlated step follows and would otherwise rebuild the integrals
         !! and refit them. See `reuse_scf_fit` for when that is.
      logical :: keep_fit

      character(len=MAX_LINE_LENGTH) :: line

      do_gradient = .false.
      if (present(want_gradient)) do_gradient = want_gradient
      do_hessian = .false.
      if (present(want_hessian)) do_hessian = want_hessian

      ! Coupled cluster stops here: its gradient needs the Lambda equations,
      ! which do not exist on this side at all. MP2's does exist, and is taken
      ! below -- but only where the reference underneath it is one the relaxed
      ! density was written for, which is the refusals after this.
      if (do_gradient .and. settings%run_cc) then
         call result%error%set(ERROR_VALIDATION, "a coupled cluster gradient needs the "// &
                               "Lambda amplitudes, which are not implemented. Run the "// &
                               "gradient at the Hartree-Fock or MP2 level, or coupled "// &
                               "cluster as an energy.")
         result%has_error = .true.
         return
      end if
      ! An MP2 gradient is not refused here any more, frozen core or not.
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
      !
      ! A frozen core used to be refused here, fronting all four routines at
      ! once, because neither assembly built the blocks it brings. Both do
      ! now: the amplitudes and the two-particle density span the active
      ! occupied space, and the relaxed density gains an occupied-frozen
      ! block resolved directly from the Lagrangian. So a default deck --
      ! `freeze_core` is on by default -- gets the gradient of the energy it
      ! computed; the frozen count is resolved in the MP2 block below and
      ! passed to whichever routine the deck selects. (No virtual-frozen
      ! block ever existed to build, whatever the old message said:
      ! `n_frozen_core` freezes leading core orbitals and no virtuals.)

      ! Same rule the GPU backend applies, so a deck does not change meaning
      ! when it moves between them: an odd electron count or a multiplicity
      ! above one has no restricted solution to find, whatever the keyword says.
      unrestricted = (fragment%multiplicity /= 1) .or. (mod(fragment%nelec, 2) /= 0) &
                     .or. settings%unrestricted

      ! Continuum solvation on this backend is an *energy*: the surface
      ! charges enter every SCF iteration and the total. What does not exist is
      ! their derivative -- the cavity moves with the nuclei and the charges
      ! respond -- so an analytic gradient would come back converged, plausible
      ! and missing the solvent's pull on every atom. And anything that runs
      ! further calculations on top of the solvated reference -- correlation,
      ! the Fukui ions -- needs its own decision about the continuum, which
      ! none of them has yet. Each combination is refused by name rather than
      ! run in a quietly mixed phase.
      if (settings%pcm%enabled) then
         if (do_gradient) then
            call result%error%set(ERROR_VALIDATION, "the PCM nuclear gradient is not "// &
                                  "implemented on the CPU backend -- the energy is, but "// &
                                  "its cavity and charge-response derivative terms are "// &
                                  "not built. Run the energy here, or the gradient on "// &
                                  "the cuEST backend, which has them.")
            result%has_error = .true.
            return
         end if
         if (settings%run_mp2 .or. settings%run_cc) then
            call result%error%set(ERROR_VALIDATION, "correlation on top of a solvated "// &
                                  "reference is not implemented on the CPU backend: the "// &
                                  "orbitals would carry the continuum while the "// &
                                  "correlation treatment ignored it, and which scheme "// &
                                  "that amounts to should be chosen, not implied. Run "// &
                                  "the SCF with keywords.pcm, or the correlation "// &
                                  "without it.")
            result%has_error = .true.
            return
         end if
         ! The Fukui ions used to be refused here, because they were solved in
         ! the gas phase and differenced against a solvated neutral -- an IP
         ! across two different physics, and wrong by roughly the solvation
         ! energy of a charged species, which is the size of the answer. They
         ! now take the same continuum the neutral did, so the refusal is gone
         ! rather than the reasoning: it was right, and what changed is that
         ! `fukui_indices` carries `pcm` through to both of them.
         !
         ! Equilibrium solvation in all three states, so what comes out is the
         ! adiabatic quantity -- the solvent relaxes around each charge state.
         ! A vertical IP would freeze the slow polarisation at the neutral's
         ! value and let only the fast part respond, which needs an optical
         ! dielectric `pcm_config_t` does not carry.
         if (settings%bonding_energy) then
            call result%error%set(ERROR_VALIDATION, "the bonding energy decomposition "// &
                                  "rebuilds its atom energies from gas-phase operators "// &
                                  "and would drop the dielectric term without saying "// &
                                  "so. Not wired up; run it without keywords.pcm.")
            result%has_error = .true.
            return
         end if
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
      ! unrestricted density fitting above. An unrestricted reference reaches the
      ! spin-orbital path, which carries its own alpha and beta transform; the
      ! fitted route does not, because `b_vv` has no spin blocks and would return
      ! the restricted answer built from the alpha orbitals alone.
      if (settings%run_cc .and. unrestricted .and. settings%corr_density_fitting) then
         call result%error%set(ERROR_VALIDATION, "density-fitted coupled cluster needs "// &
                               "a restricted reference: the fitted three-index block has "// &
                               "no spin blocks. Run unrestricted coupled cluster without "// &
                               "an auxiliary basis, or use a closed-shell system.")
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
                                  mol, error, ghost=ghost_of(fragment), &
                                  force_cartesian=settings%cartesian, &
                                  ecp_name=trim(settings%ecp_set))
      if (error%has_error()) then
         call result%error%set(ERROR_VALIDATION, error%get_message())
         result%has_error = .true.
         return
      end if

      ! The electrons the SCF actually solves for. An ECP replaced the core,
      ! so they are not there to be solved for -- and `fragment` cannot know
      ! that, because it is built before any basis or potential is read.
      !
      ! Every use below this point is the valence count. The parity test
      ! further up is not, and does not need to be: a core is a closed shell,
      ! so `core_electrons` is always even and the choice between RHF and UHF
      ! is the same either way.
      nelec = fragment%nelec - sum(mol%core_electrons)
      if (nelec < 0) then
         call result%error%set(ERROR_VALIDATION, "the effective core potential "// &
                               "replaces more electrons than the fragment has")
         result%has_error = .true.
         call mol%destroy()
         return
      end if

      ! The SCF banner below reports the basis but not who it is solving for:
      ! a charged or open-shell fragment reads identically to a neutral one
      ! without this, which is exactly the case a fragmented run most needs to
      ! see. Same integers the breakdown CSV carries.
      call logger%verbose("charge "//int_to_text(fragment%charge)// &
                          ", multiplicity "//int_to_text(fragment%multiplicity)// &
                          ", electrons "//int_to_text(nelec))

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
                                level=settings%grid_level, polarized=unrestricted, &
                                nlc_level=settings%nlc_grid_level, &
                                screen_tol=settings%screening_tolerance, &
                                point_block=settings%block_size)
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
         if (settings%verbose) then
            write (line, "(a,a,a,i0)") "  functional: ", trim(settings%functional), ", grid level ", settings%grid_level
            call logger%info(trim(line))
         end if
         ! A double hybrid's perturbative term is an MP2 on top of the
         ! reference, and correlation over a solvated reference is refused
         ! above by name. This is the same refusal for the case where the deck
         ! spells MP2 as part of the functional.
         if (settings%pcm%enabled .and. xc%pt2_fraction /= 0.0_dp) then
            call result%error%set(ERROR_VALIDATION, "a double hybrid under continuum "// &
                                  "solvation is correlation on a solvated reference, "// &
                                  "which is not implemented on the CPU backend -- see "// &
                                  "the MP2 refusal. Use a hybrid or GGA functional "// &
                                  "with keywords.pcm.")
            result%has_error = .true.
            if (kohn_sham) call xc%destroy()
            call mol%destroy()
            return
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
                                 fragment%coordinates, nelec, mol, guess_total, &
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

      ! ---- the continuum, once per geometry ---------------------------------
      !
      ! Built here because this is the first point where the molecule exists;
      ! the SCF then solves the surface charges against each iteration's
      ! density. The banner mirrors the cuEST driver's, so a run reads the same
      ! on either backend.
      if (settings%pcm%enabled) then
         call pcm_ctx%build(mol, fragment%element_numbers, settings%pcm, error)
         if (error%has_error()) then
            call result%error%set(ERROR_VALIDATION, error%get_message())
            result%has_error = .true.
            call mol%destroy()
            return
         end if
         write (line, "(a,f10.4,a,i0,a,f6.3,a,a)") "    continuum: eps = ", &
            settings%pcm%dielectric, "   angular points ", settings%pcm%angular_points, &
            "   radii x ", settings%pcm%radii_scale, "   model ", trim(settings%pcm%method)
         call logger%info(trim(line))
         write (line, "(a,i0,a,i0,a)") "    continuum: ", pcm_ctx%surface%n_points, &
            " surface points kept, ", pcm_ctx%surface%n_discarded, " buried"
         call logger%info(trim(line))
      end if

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
                                     aux, error, ghost=ghost_of(fragment), &
                                     force_cartesian=settings%cartesian)
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

         ! Whether to carry the three-centre integrals out of the SCF. They
         ! cost a full copy of the tensor to keep, and nothing at all if the
         ! run ends at the reference -- so it is asked for only when something
         ! downstream will transform them.
         keep_fit = reuse_scf_fit(settings)
         if (keep_fit) then
            call run_libcint_rhf(mol, nelec, settings%max_iter, settings%energy_tol, &
                                 settings%density_tol, settings%verbose, scf, error, &
                                 aux=aux, diis_vectors=diis_size, guess=guess_kind, &
                                 guess_density=guess_total, xc=xc, pcm=pcm_ctx, &
                                 level_shift=settings%level_shift, &
                                 linear_dependence=settings%linear_dependence, &
                                 b_ao_out=scf_b_ao)
         else
            call run_libcint_rhf(mol, nelec, settings%max_iter, settings%energy_tol, &
                                 settings%density_tol, settings%verbose, scf, error, &
                                 aux=aux, diis_vectors=diis_size, guess=guess_kind, &
                                 guess_density=guess_total, xc=xc, pcm=pcm_ctx, &
                                 level_shift=settings%level_shift, &
                                 linear_dependence=settings%linear_dependence)
         end if
         ! Kept alive: the gradient below has to be told the same auxiliary
         ! basis this SCF fitted with. Released once past it.
      else if (unrestricted) then
         call run_libcint_uhf(mol, nelec, fragment%multiplicity, settings%max_iter, &
                              settings%energy_tol, settings%density_tol, settings%verbose, &
                              scf, error, diis_vectors=diis_size, guess=guess_kind, &
                              guess_density_alpha=guess_a, guess_density_beta=guess_b, xc=xc, &
                              pcm=pcm_ctx, level_shift=settings%level_shift, &
                              linear_dependence=settings%linear_dependence)
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
         call run_libcint_rhf(mol, nelec, settings%max_iter, settings%energy_tol, &
                              settings%density_tol, settings%verbose, scf, error, &
                              diis_vectors=diis_size, guess=guess_kind, &
                              guess_density=guess_total, xc=xc, pcm=pcm_ctx, &
                              in_core=eri_fits_in_core(mol%nao) .and. .not. xc%range_separated, &
                              level_shift=settings%level_shift, &
                              linear_dependence=settings%linear_dependence)
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
         call result%error%set(ERROR_VALIDATION, &
                               scf_not_converged_message(scf%iterations))
         result%has_error = .true.
         call aux%destroy()
         call mol%destroy()
         return
      end if

      result%energy%scf = scf%energy
      result%has_energy = .true.
      if (settings%pcm%enabled) then
         write (line, "(a,f18.10,a,f9.4)") "    continuum: E_diel = ", scf%pcm_energy, &
            "   total surface charge ", scf%pcm_charge
         call logger%info(trim(line))
      end if

      ! ---- analyses asked for alongside the energy ---------------------------
      !
      ! Run here, on the converged orbitals, because this is where they exist.
      ! A failure is reported and then dropped rather than propagated: the
      ! energy above is correct and already stored, and an analysis that cannot
      ! run -- an element the minimal basis does not cover, an open shell the
      ! construction is not written for -- is not a reason to throw it away.
      ! Where the molecule reacts, by difference in the electron count. Two
      ! further SCFs on the same `mol`, so the ions see the same geometry and
      ! the same basis functions by construction.
      if (allocated(settings%fukui_population)) then
         block
            use mqc_libcint_fukui, only: fukui_result_t, fukui_indices, print_fukui_report
            type(fukui_result_t) :: fukui
            type(error_t) :: fukui_error
            integer :: fukui_frozen
            real(dp) :: fukui_pt2
            type(libcint_molecule_t), target :: fukui_aux
            type(libcint_molecule_t), pointer :: fukui_aux_arg
            logical :: fukui_fitted

            ! Zero unless this is a Kohn-Sham run with a perturbative term.
            ! `xc` is only built when kohn_sham, so reading its fraction
            ! unconditionally would read an uninitialised context on a
            ! Hartree-Fock run.
            fukui_pt2 = 0.0_dp
            if (kohn_sham) fukui_pt2 = xc%pt2_fraction

            ! Fitted exactly when the energy's own PT2 term is fitted. The
            ! alternative -- Fukui always exact, the energy fitted when asked --
            ! would print an IP beside a total energy that used a different
            ! correlation treatment, and the difference between them would look
            ! like a property of the molecule.
            fukui_fitted = fukui_pt2 /= 0.0_dp .and. settings%aux_basis_named
            if (fukui_fitted) then
               call correlation_aux_basis(settings, fragment, symbols, fukui_aux, fukui_error)
               if (fukui_error%has_error()) then
                  call logger%warning("  the Fukui analysis could not build its "// &
                                      "auxiliary basis: "//fukui_error%get_message())
                  fukui_fitted = .false.
                  call fukui_error%clear()
               end if
            end if
            ! Absent rather than empty when unfitted. A disassociated pointer
            ! passed to an optional dummy is not present, which is the same
            ! idiom `xc_arg` uses a few hundred lines up.
            fukui_aux_arg => null()
            if (fukui_fitted) fukui_aux_arg => fukui_aux

            ! The functional by NAME, not this routine's `xc`: that context
            ! was built spin-unpolarised for the closed-shell neutral and the
            ! ions are doublets. fukui_indices builds its own polarised one.
            ! A double hybrid needs the orbitals as well as the density: its
            ! perturbative term is evaluated for all three states inside, so
            ! that the neutral's comes from the same code path the ions' do.
            ! `dh_frozen` is resolved the same way the energy's PT2 resolves
            ! it, a few hundred lines below -- a Fukui run that froze a
            ! different core than the energy would report an IP that did not
            ! match the numbers beside it.
            fukui_frozen = settings%n_frozen_core
            if (fukui_frozen < 0) then
               if (ecp_refuses_auto_frozen_core(mol%core_electrons, error)) then
                  call result%error%set(ERROR_VALIDATION, error%get_message())
                  result%has_error = .true.
                  return
               end if
               fukui_frozen = core_orbital_count(fragment%element_numbers)
            end if
            if (.not. settings%freeze_core) fukui_frozen = 0
            call fukui_indices(mol, nelec, fragment%multiplicity, scf%density, &
                               scf%energy, settings%fukui_population, settings%max_iter, &
                               settings%energy_tol, settings%density_tol, fukui, &
                               fukui_error, functional=trim(settings%functional), &
                               grid_level=settings%grid_level, &
                               pt2_fraction=fukui_pt2, &
                               neutral_orbitals=scf%orbitals, &
                               neutral_orbital_energies=scf%orbital_energies, &
                               n_frozen=fukui_frozen, aux=fukui_aux_arg, &
                               verbose=settings%verbose, pcm=pcm_ctx, &
                               level_shift=settings%level_shift, &
                               diis_vectors=diis_size)
            if (fukui_fitted) call fukui_aux%destroy()
            if (fukui_error%has_error()) then
               call logger%warning("  the Fukui analysis could not run: "// &
                                   fukui_error%get_message())
            else
               call print_fukui_report(fukui, symbols)
               result%fukui_plus = fukui%f_plus
               result%fukui_minus = fukui%f_minus
               result%fukui_dual = fukui%dual
               result%fukui_ip = fukui%ionisation_potential
               result%fukui_ea = fukui%electron_affinity
               result%fukui_hardness = fukui%hardness
               result%fukui_electrophilicity = fukui%electrophilicity
               result%fukui_anion_bound = fukui%anion_bound
               result%fukui_scheme = fukui%scheme
               result%has_fukui = .true.
            end if
         end block
      end if

      ! ---- atomic partial charges -------------------------------------------
      !
      ! Whatever reference just converged, partitioned. Nothing below knows or
      ! cares whether the density came from Hartree-Fock or a Kohn-Sham
      ! functional, restricted or unrestricted -- both schemes take a density
      ! matrix, and the converged one is already here. That is why this costs
      ! no second SCF: a deck asking for charges is asking for the charges *of
      ! the calculation it ran*, not of a Hartree-Fock stand-in for it.
      if (allocated(settings%charges_scheme)) then
         block
            real(dp), allocatable :: q(:), q_spin(:)
            real(dp), allocatable :: total_density(:, :), spin_density(:, :), s_ao(:, :)
            logical :: open_shell

            call analysis_error%clear()
            open_shell = allocated(scf%density_beta)
            if (open_shell) then
               total_density = scf%density + scf%density_beta
               spin_density = scf%density - scf%density_beta
            else
               ! Already the total: the closed-shell build is D = 2 sum_i c_i c_i^T.
               total_density = scf%density
            end if

            select case (trim(settings%charges_scheme))
            case ("mulliken")
               call mol%overlap(s_ao)
               call mulliken_charges(mol, total_density, s_ao, q, analysis_error)
               ! Free while the overlap is in hand, and the number anyone
               ! running an open shell actually wants: where the unpaired
               ! density sits.
               if (open_shell .and. .not. analysis_error%has_error()) then
                  call mulliken_spin_populations(mol, spin_density, s_ao, q_spin, &
                                                 analysis_error)
               end if
            case ("chelpg")
               ! The fit is solved under a total-charge constraint, so the
               ! fragment's charge has to be handed over. Caps carry none, but
               ! they are part of the molecule the SCF saw and so are part of
               ! what the constraint is about.
               call chelpg_charges(mol, total_density, q, analysis_error, &
                                   total_charge=real(fragment%charge, dp))
            case default
               call analysis_error%set(ERROR_VALIDATION, "unknown charge scheme '"// &
                                       trim(settings%charges_scheme)// &
                                       "'; expected 'mulliken' or 'chelpg'.")
            end select

            if (analysis_error%has_error()) then
               call logger%warning("  atomic charges could not be computed: "// &
                                   analysis_error%get_message())
            else
               call move_alloc(q, result%atomic_charges)
               if (allocated(q_spin)) call move_alloc(q_spin, result%spin_populations)
               result%charge_scheme = trim(settings%charges_scheme)
               result%has_charges = .true.
            end if
         end block
      end if

      if (bonding_analysis_kind(settings%bonding_analysis) == BONDING_GMS_QUAO) then
         call analysis_error%clear()
         ! Printed whatever the verbosity. Asking for a property *is* the
         ! request for its output -- a deck that names an analysis and gets a
         ! silent run has been given nothing, and the natural reading is that
         ! the analysis found nothing to say.
         call run_quao_analysis(mol, fragment%element_numbers, symbols, &
                                fragment%coordinates, scf%orbitals, nelec, &
                                analysis_error, threshold=settings%bonding_threshold, &
                                energy_decomposition=settings%bonding_energy, &
                                no_sharing=settings%bonding_no_sharing, &
                                no_sharing_ci=settings%bonding_no_sharing_ci, &
                                restrict_localization=settings%bonding_restrict_localization, &
                                atom_energy=ieda_atom, free_atom_energy=ieda_free, &
                                pair_energy=ieda_pair, pair_classical=ieda_classical, &
                                formation_energy=ieda_formation)
         if (analysis_error%has_error()) then
            call logger%warning("  the bonding analysis could not run: "// &
                                analysis_error%get_message())
         end if
         call store_decomposition(result, ieda_atom, ieda_free, ieda_pair, &
                                  ieda_classical, ieda_formation)
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
            !
            ! A frozen core, by contrast, is still refused -- and no longer
            ! because the blocks are missing: the MP2 gradient routines
            ! underneath build them, and a plain MP2 deck takes them. This
            ! assembly never passes the frozen count, and the occupied-frozen
            ! resolution has been validated over a Hartree-Fock reference
            ! only, not against a Kohn-Sham operator with its kernel in the
            ! response.
            if (settings%freeze_core .and. settings%n_frozen_core /= 0) then
               call result%error%set(ERROR_VALIDATION, "a double hybrid gradient is "// &
                                     "all-electron: the perturbative term's gradient "// &
                                     "has not been taught, or validated with, a frozen "// &
                                     "core over a Kohn-Sham reference -- plain MP2 has. "// &
                                     "The *energy* does honour freeze_core, so this "// &
                                     "refuses rather than returning the all-electron "// &
                                     "gradient of a frozen-core energy. Run the energy "// &
                                     "with it, or the gradient with freeze_core off.")
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

      ! ---- the analytic Hessian, where it applies ---------------------------
      !
      ! A restricted reference over exact integrals, Hartree-Fock or Kohn-Sham.
      ! Every other case here reaches this same line and every one of them must
      ! fall through: an unrestricted reference has two densities where these
      ! assume one carrying its factor of two; a fitted SCF converged a
      ! different energy from the one these derivatives belong to; and a
      ! correlated method's Hessian needs its own response entirely. None of
      ! those would fail loudly -- each would return a plausible, converged,
      ! wrong matrix -- so the guard is a list of positives rather than a list
      ! of refusals.
      !
      ! A functional carrying VV10 takes this path too: `ks_hessian` carries
      ! the non-local term's second derivative, Fock derivative and response
      ! kernel, so a `-V` functional no longer falls back to the
      ! semi-numerical path.
      !
      ! Hydrogen caps are excluded for a different reason: the shapes match and
      ! the numbers would be right, but a capped fragment's second derivatives
      ! have to be redistributed onto the heavy atoms the caps replaced, and
      ! that is checked for the finite-difference path and not for this one.
      ! Falling back leaves a fragmented Hessian exactly as it was.
      ! A solvated SCF also falls through: these second derivatives are the
      ! gas-phase operator's, and the finite-difference fallback is *correct*
      ! for the continuum -- each displaced energy rebuilds its own cavity.
      ! The double hybrid stays on this list now VV10 has come off it, and
      ! `settings%run_mp2` does not stand in for it: a double hybrid's PT2
      ! correlation is driven from `xc%pt2_fraction` further down and never
      ! sets that flag. None of the second derivatives below carries it, so
      ! `b2plyp` taken analytically returns its underlying hybrid's Hessian --
      ! a Frobenius norm of 2.146683 against a true 2.120638 on water at
      ! STO-3G, printed beside the double hybrid's own "PT2 x 0.270" energy
      ! line with nothing to say the Hessian excludes what the energy includes.
      ! That is the failure the *gradient* guard above this one already
      ! prevents, where the hybrid's number came back against the double
      ! hybrid's energy, right in shape and wrong by a third.
      if (do_hessian .and. .not. unrestricted &
          .and. .not. settings%density_fitting .and. .not. settings%run_mp2 &
          .and. .not. settings%run_cc .and. .not. settings%pcm%enabled &
          .and. fragment%n_caps == 0 &
          .and. .not. (kohn_sham .and. xc%pt2_fraction /= 0.0_dp)) then
         block
            real(dp), allocatable :: hess4(:, :, :, :)
            type(timer_type) :: hess_clock

            call logger%info("  computing the analytic Hessian")
            call hess_clock%start()
            if (kohn_sham) then
               call ks_hessian(mol, fragment%element_numbers, scf%density, scf%orbitals, &
                               scf%orbital_energies, scf%n_occupied, xc, xc%exx_fraction, &
                               hess4, error)
            else
               call rhf_hessian(mol, fragment%element_numbers, scf%density, scf%orbitals, &
                                scf%orbital_energies, scf%n_occupied, hess4, error)
            end if
            if (error%has_error()) then
               call result%error%set(ERROR_VALIDATION, "Hessian: "//error%get_message())
               result%has_error = .true.
               call aux%destroy()
               call mol%destroy()
               return
            end if

            call hessian_to_matrix(hess4, result%hessian)
            result%has_hessian = .true.
            deallocate (hess4)
            write (line, "(a,f10.2,a)") "  Hessian done in ", hess_clock%get_elapsed_time(), " s"
            call logger%info(trim(line))
         end block
      end if

      call aux%destroy()

      if (settings%run_mp2) then
         block
            type(mp2_result_t) :: mp2
            integer :: frozen

            frozen = settings%n_frozen_core
            if (frozen < 0) then
               if (ecp_refuses_auto_frozen_core(mol%core_electrons, error)) then
                  call result%error%set(ERROR_VALIDATION, error%get_message())
                  result%has_error = .true.
                  return
               end if
               frozen = core_orbital_count(fragment%element_numbers)
            end if
            if (.not. settings%freeze_core) frozen = 0

            if (settings%corr_density_fitting) then
               call correlation_aux_basis(settings, fragment, symbols, corr_aux, error)
               if (error%has_error()) then
                  call result%error%set(ERROR_VALIDATION, error%get_message())
                  result%has_error = .true.
                  call mol%destroy()
                  return
               end if
               if (allocated(scf_b_ao)) then
                  call run_libcint_ri_mp2(mol, corr_aux, scf%orbitals, scf%orbital_energies, &
                                          nelec/2, scf%energy, mp2, error, &
                                          n_frozen=frozen, b_ao_in=scf_b_ao)
               else
                  call run_libcint_ri_mp2(mol, corr_aux, scf%orbitals, scf%orbital_energies, &
                                          nelec/2, scf%energy, mp2, error, &
                                          n_frozen=frozen)
               end if
               ! Released as soon as it has been transformed. It is gigabytes,
               ! and the gradient below builds its own.
               if (allocated(scf_b_ao)) deallocate (scf_b_ao)
               call corr_aux%destroy()
            else
               call run_libcint_mp2(mol, scf%orbitals, scf%orbital_energies, &
                                    nelec/2, scf%energy, mp2, error, &
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
                                               scf%orbital_energies, nelec/2, &
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
                                            nelec/2, result%gradient, error, &
                                            n_frozen=frozen, aux=corr_aux)
                  call corr_aux%destroy()
               else
                  call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, &
                                            nelec/2, result%gradient, error, &
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

            ! ---- the analytic Hessian, where it applies ------------------
            !
            ! The same positive-list shape as the reference gate above, for
            ! the same reason: everything not on the list would return a
            ! plausible, converged, wrong matrix rather than fail loudly. A
            ! restricted reference over exact integrals, unscaled MP2,
            ! all-electron or frozen-core -- the core count flows straight
            ! into `mp2_correlation_hessian`, whose Phase 2 rotations
            ! (the Brillouin rewrite of `U^X`, the pair augmentation and the
            ! Sylvester core<->active derivative) carry it, so an ordinary
            ! `freeze_core` deck now takes this path instead of central
            ! differences.
            if (do_hessian .and. .not. unrestricted &
                .and. .not. settings%density_fitting &
                .and. .not. settings%corr_density_fitting &
                .and. .not. settings%pcm%enabled .and. fragment%n_caps == 0 &
                .and. settings%scs_ss == 1.0_dp .and. settings%scs_os == 1.0_dp) then
               block
                  real(dp), allocatable :: hcorr(:, :, :, :), href(:, :, :, :)
                  real(dp), allocatable :: hnuc(:, :, :, :), hresp(:, :, :, :)
                  type(timer_type) :: hess_clock

                  call logger%info("  computing the analytic MP2 Hessian")
                  call hess_clock%start()
                  ! The correlation block plus the reference's AO-dependent
                  ! skeleton come from one shared integral sweep; the
                  ! reference's CPHF response and the nuclear repulsion
                  ! complete the total, which is `rhf_hessian`'s own split
                  ! (the assembly test holds the identity at 2.5e-14).
                  call mp2_correlation_hessian(mol, scf%orbitals, &
                                               scf%orbital_energies, scf%density, &
                                               scf%n_occupied, frozen, hcorr, href, error)
                  if (.not. error%has_error()) then
                     call nuclear_repulsion_hessian(fragment%element_numbers, &
                                                    mol%coords, hnuc, error)
                  end if
                  if (.not. error%has_error()) then
                     call response_hessian(mol, scf%density, scf%orbitals, &
                                           scf%orbital_energies, scf%n_occupied, &
                                           hresp, error)
                  end if
                  if (error%has_error()) then
                     call result%error%set(ERROR_VALIDATION, "MP2 Hessian: "// &
                                           error%get_message())
                     result%has_error = .true.
                     call mol%destroy()
                     return
                  end if
                  hcorr = hcorr + href + hresp + hnuc
                  call hessian_to_matrix(hcorr, result%hessian)
                  result%has_hessian = .true.
                  write (line, "(a,f10.2,a)") "  Hessian done in ", &
                     hess_clock%get_elapsed_time(), " s"
                  call logger%info(trim(line))
               end block
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
            integer :: n_alpha, n_beta
            real(dp) :: cc_mp2, cc_singles, cc_doubles, cc_triples

            ! The same frozen-core rule the MP2 block applies, deliberately
            ! duplicated rather than hoisted: the two blocks are independent, and
            ! a shared local would silently couple them if one ever wanted a
            ! different count.
            frozen = settings%n_frozen_core
            if (frozen < 0) then
               if (ecp_refuses_auto_frozen_core(mol%core_electrons, error)) then
                  call result%error%set(ERROR_VALIDATION, error%get_message())
                  result%has_error = .true.
                  return
               end if
               frozen = core_orbital_count(fragment%element_numbers)
            end if
            if (.not. settings%freeze_core) frozen = 0

            ! Which formulation. Both are exact for a closed shell and agree to
            ! machine precision -- that identity is asserted by
            ! test_mqc_libcint_rcc -- so this chooses how the same number is
            ! computed, not which number. Spatial by default because it is
            ! smaller and faster; the spin-orbital path is what a doubtful
            ! result gets checked against.
            spin_adapted = settings%cc_spin_adapted
            ! Not a preference when the reference is unrestricted. The
            ! spin-adapted formulation is derived for a closed shell and has no
            ! beta orbitals to be given; the spin-orbital one takes them.
            if (unrestricted) spin_adapted = .false.
            ! The same split the unrestricted SCF made: the excess spin is
            ! alpha, so multiplicity 2S+1 puts S more electrons there.
            n_alpha = (nelec + fragment%multiplicity - 1)/2
            n_beta = nelec - n_alpha

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
                                         nelec/2, frozen, settings%cc_max_iter, &
                                         settings%cc_tolerance, settings%cc_triples, &
                                         settings%verbose, rcc, error, &
                                         diis_vectors=settings%cc_diis_size, aux=corr_aux)
               else
                  call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, &
                                        nelec/2, frozen, settings%cc_max_iter, &
                                        settings%cc_tolerance, settings%cc_triples, &
                                        settings%verbose, cc, error, &
                                        diis_vectors=settings%cc_diis_size, aux=corr_aux)
               end if
               call corr_aux%destroy()
            else
               if (spin_adapted) then
                  call run_libcint_rccsd(mol, scf%orbitals, scf%orbital_energies, &
                                         nelec/2, frozen, settings%cc_max_iter, &
                                         settings%cc_tolerance, settings%cc_triples, &
                                         settings%verbose, rcc, error, &
                                         diis_vectors=settings%cc_diis_size)
               else if (unrestricted) then
                  call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, &
                                        n_alpha, frozen, settings%cc_max_iter, &
                                        settings%cc_tolerance, settings%cc_triples, &
                                        settings%verbose, cc, error, &
                                        diis_vectors=settings%cc_diis_size, &
                                        coeff_b=scf%orbitals_beta, &
                                        orbital_energies_b=scf%orbital_energies_beta, &
                                        n_occ_b=n_beta)
               else
                  call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, &
                                        nelec/2, frozen, settings%cc_max_iter, &
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
               if (dh_frozen < 0) then
               if (ecp_refuses_auto_frozen_core(mol%core_electrons, error)) then
                  call result%error%set(ERROR_VALIDATION, error%get_message())
                  result%has_error = .true.
                  return
               end if
               dh_frozen = core_orbital_count(fragment%element_numbers)
               end if
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
               ! Four routines rather than two, because the reference decides
               ! the correlation treatment as much as the auxiliary basis does.
               ! `nelec/2` is not an occupied count for an open shell,
               ! so the unrestricted calls take the SCF's own per-spin counts
               ! rather than deriving one.
               if (settings%aux_basis_named) then
                  call correlation_aux_basis(settings, fragment, symbols, corr_aux, error)
                  if (error%has_error()) then
                     call result%error%set(ERROR_VALIDATION, error%get_message())
                     result%has_error = .true.
                     call xc%destroy()
                     call mol%destroy()
                     return
                  end if
                  if (unrestricted) then
                     call run_libcint_uri_mp2(mol, corr_aux, scf%orbitals, scf%orbitals_beta, &
                                              scf%orbital_energies, scf%orbital_energies_beta, &
                                              scf%n_occupied, scf%n_occupied_beta, &
                                              scf%energy, dh_mp2, error, n_frozen=dh_frozen)
                  else
                     call run_libcint_ri_mp2(mol, corr_aux, scf%orbitals, &
                                             scf%orbital_energies, nelec/2, &
                                             scf%energy, dh_mp2, error, n_frozen=dh_frozen)
                  end if
                  call corr_aux%destroy()
               else
                  if (unrestricted) then
                     call run_libcint_ump2(mol, scf%orbitals, scf%orbitals_beta, &
                                           scf%orbital_energies, scf%orbital_energies_beta, &
                                           scf%n_occupied, scf%n_occupied_beta, &
                                           scf%energy, dh_mp2, error, n_frozen=dh_frozen)
                  else
                     call run_libcint_mp2(mol, scf%orbitals, scf%orbital_energies, &
                                          nelec/2, scf%energy, dh_mp2, error, &
                                          n_frozen=dh_frozen)
                  end if
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
                                                     nelec/2, dh_grad, error, &
                                                     n_frozen=dh_frozen, &
                                                     fitted_reference=settings%density_fitting, &
                                                     xc=xc, scf_density=scf%density, &
                                                     pt2_scale=xc%pt2_fraction)
                     else
                        call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, &
                                                  nelec/2, dh_grad, error, &
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
                        int_to_text(n_active_electrons)//" electrons cannot have "// &
                        "multiplicity "//int_to_text(fragment%multiplicity)//": the two "// &
                        "differ in parity, so there is no way to split the electrons "// &
                        "into alpha and beta counts. Change either the electron count "// &
                        "or the multiplicity by one.")
         return
      end if
      n_alpha = (n_active_electrons + unpaired)/2
      n_beta = (n_active_electrons - unpaired)/2

      if (n_beta < 0 .or. n_alpha > n_active) then
         call error%set(ERROR_VALIDATION, "an active space of "//int_to_text(n_active)// &
                        " orbitals cannot hold "//int_to_text(n_alpha)//" alpha and "// &
                        int_to_text(n_beta)//" beta electrons at multiplicity "// &
                        int_to_text(fragment%multiplicity)//". Widen the active space or "// &
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
                           int_to_text(fragment%nelec)//" electrons less an active space "// &
                           "of "//int_to_text(n_active_electrons)//" leaves "// &
                           int_to_text(closed_shell_electrons)//", which is not a "// &
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
                           "for every electron: "//int_to_text(n_inactive)//" inactive "// &
                           "orbitals hold "//int_to_text(2*n_inactive)//" electrons and "// &
                           "the active space "//int_to_text(n_active_electrons)//", which "// &
                           "is "//int_to_text(2*n_inactive + n_active_electrons)//" against "// &
                           "the molecule's "//int_to_text(fragment%nelec)//".")
            return
         end if
      end if

      ! The virtual space may be empty -- a CAS in a minimal basis legitimately
      ! uses every orbital -- but it may not be negative.
      if (n_inactive + n_active > n_ao) then
         call error%set(ERROR_VALIDATION, int_to_text(n_inactive)//" inactive plus "// &
                        int_to_text(n_active)//" active orbitals is more than the "// &
                        int_to_text(n_ao)//" the basis '"//trim(settings%basis_set)// &
                        "' provides. Use a larger basis or a smaller active space.")
         return
      end if

      space = [n_inactive, n_active, n_alpha, n_beta]
   end subroutine resolve_active_space

   subroutine run_libcint_mcscf(settings, fragment, result, want_gradient)
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
      logical, intent(in), optional :: want_gradient
         !! The nuclear gradient of the optimised wave function. Only a CASSCF
         !! has one here: a CASCI leaves its orbitals where the SCF put them, so
         !! it is not stationary with respect to them and the terms this omits
         !! are not zero.

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casscf_result_t) :: casscf
      type(casci_result_t) :: casci
      type(valence_wavefunction_t) :: converged
         !! Offered to the bonding analysis, which takes it up only if it is
         !! over the full valence space. Filled unconditionally because filling
         !! it is four assignments and deciding whether it qualifies is not this
         !! layer's job.
      type(error_t) :: error
      type(error_t) :: analysis_error
      real(dp), allocatable :: ieda_atom(:), ieda_free(:)
      real(dp), allocatable :: ieda_pair(:, :), ieda_classical(:, :)
      real(dp) :: ieda_formation
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

      ! Refused here rather than skipped after the fact, and refused before the
      ! CASSCF rather than after it. A deck that asked for charges and got a
      ! silent run has been told nothing, and would reasonably read the absence
      ! as "no charges to report"; making it wait for a converged active space
      ! first would be worse still. The 1-RDM this path produces is in the MO
      ! basis over orbitals with fractional occupation, so the AO density both
      ! partition schemes want has to be assembled -- which is work, not an
      ! oversight.
      if (allocated(settings%charges_scheme)) then
         call result%error%set(ERROR_VALIDATION, "atomic charges (properties.charges) "// &
                               "are not implemented for a multiconfigurational wave "// &
                               "function: the partition needs an AO density this path "// &
                               "does not form. Ask for them on the reference method "// &
                               "instead, or drop properties.charges.")
         result%has_error = .true.
         return
      end if

      if (mod(fragment%nelec, 2) /= 0) then
         call result%error%set(ERROR_VALIDATION, "a multiconfigurational calculation "// &
                               "here starts from a closed-shell SCF, and "// &
                               int_to_text(fragment%nelec)//" electrons have no restricted "// &
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
                                  mol, error, ghost=ghost_of(fragment), &
                                  force_cartesian=settings%cartesian)
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
                           diis_vectors=diis_size, &
                           linear_dependence=settings%linear_dependence)
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
                               int_to_text(settings%max_iter)//" iterations, so there are no "// &
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
               ! found at all, and more iterations will do nothing. It happens
               ! on flat surfaces, where the energy is converged long before the
               ! gradient is small.
               call error%set(ERROR_VALIDATION, "the orbital optimisation stopped "// &
                              "improving at a gradient of "// &
                              to_real_text(casscf%gradient_norm)//", short of the "// &
                              to_real_text(settings%mcscf%orbital_convergence)// &
                              " asked for, after "//int_to_text(casscf%iterations)// &
                              " macro-iterations. More iterations will not help -- "// &
                              "no step downhill could be found. Ask for a looser "// &
                              "keywords.mcscf.orbital_convergence; on a flat surface "// &
                              "the energy is converged well before the gradient is.")
            else
               call error%set(ERROR_VALIDATION, "the orbital optimisation did not reach an "// &
                              "orbital gradient of "// &
                              to_real_text(settings%mcscf%orbital_convergence)//" in "// &
                              int_to_text(settings%mcscf%max_macro_iter)//" macro-iterations "// &
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
         converged%n_inactive = space(1)
         converged%n_active = space(2)
         converged%n_alpha = space(3)
         converged%n_beta = space(4)
         if (settings%mcscf%optimize_orbitals) then
            call natural_orbitals(casscf%orbitals, space(1), space(2), casscf%dm1, &
                                  natural, occupations, analysis_error)
            if (allocated(casscf%ci_vector)) converged%ci = casscf%ci_vector
            converged%orbitals = casscf%orbitals(:, space(1) + 1:space(1) + space(2))
            converged%energy = casscf%energy
         else
            call natural_orbitals(reference, space(1), space(2), casci%dm1, &
                                  natural, occupations, analysis_error)
            ! A restricted space leaves `ci_vector` unallocated and carries
            ! its coefficients in `ci_flat`, addressed by a partition that is
            ! meaningless without it -- so both go over, and the analysis writes
            ! them out over the complete space before rotating.
            if (allocated(casci%ci_vector)) converged%ci = casci%ci_vector
            if (allocated(casci%ci_flat)) then
               converged%ci_flat = casci%ci_flat
               converged%ormas = casci%ormas
            end if
            converged%orbitals = reference(:, space(1) + 1:space(1) + space(2))
            converged%energy = casci%energy
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
                                      energy_decomposition=settings%bonding_energy, &
                                      no_sharing=settings%bonding_no_sharing, &
                                      no_sharing_ci=settings%bonding_no_sharing_ci, &
                                      restrict_localization=settings%bonding_restrict_localization, &
                                      valence_wavefunction=converged, &
                                      atom_energy=ieda_atom, &
                                      free_atom_energy=ieda_free, &
                                      pair_energy=ieda_pair, &
                                      pair_classical=ieda_classical, &
                                      formation_energy=ieda_formation)
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
                                      energy_decomposition=settings%bonding_energy, &
                                      no_sharing=settings%bonding_no_sharing, &
                                      no_sharing_ci=settings%bonding_no_sharing_ci, &
                                      restrict_localization=settings%bonding_restrict_localization, &
                                      valence_wavefunction=converged, &
                                      atom_energy=ieda_atom, &
                                      free_atom_energy=ieda_free, &
                                      pair_energy=ieda_pair, &
                                      pair_classical=ieda_classical, &
                                      formation_energy=ieda_formation)
            end if
         end if
         if (analysis_error%has_error()) then
            call logger%warning("  the bonding analysis could not run: "// &
                                analysis_error%get_message())
         end if
         call store_decomposition(result, ieda_atom, ieda_free, ieda_pair, &
                                  ieda_classical, ieda_formation)
      end if

      result%scf_status = SCF_CONVERGED
      result%has_energy = .true.

      if (present(want_gradient)) then
         if (want_gradient) then
            call mcscf_gradient_into(settings, mol, casscf, space, result)
         end if
      end if

      call mol%destroy()
   end subroutine run_libcint_mcscf

   subroutine mcscf_gradient_into(settings, mol, casscf, space, result)
      !! The nuclear gradient of a converged CASSCF, onto the result
      !!
      !! Refused rather than approximated when the optimisation did not
      !! converge. The formula differentiates a *stationary* energy: away from a
      !! stationary point the orbital-response terms it omits are not zero, and
      !! what comes back is a plausible vector that is the derivative of
      !! nothing. `casscf_result_t` says the same thing about its densities,
      !! which on an unconverged run are one orbital step behind the orbitals
      !! they are reported beside.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(libcint_molecule_t), intent(in) :: mol
      type(casscf_result_t), intent(in) :: casscf
      integer, intent(in) :: space(:)
      type(calculation_result_t), intent(inout) :: result

      type(error_t) :: error

      if (.not. settings%mcscf%optimize_orbitals) then
         call result%error%set(ERROR_VALIDATION, "a CASCI gradient needs the orbital "// &
                               "response, which is not implemented: its orbitals came "// &
                               "from the SCF and were never optimised for this active "// &
                               "space, so the wave function is not stationary with "// &
                               "respect to them. Ask for 'casscf' instead.")
         result%has_error = .true.
         return
      end if
      if (.not. casscf%converged) then
         call result%error%set(ERROR_VALIDATION, "the orbital optimisation did not "// &
                               "converge, so there is no stationary point to "// &
                               "differentiate. The gradient of an unconverged MCSCF "// &
                               "is not the gradient of anything.")
         result%has_error = .true.
         return
      end if

      call libcint_mcscf_gradient(mol, casscf%orbitals, space(1), space(2), &
                                  casscf%dm1, casscf%dm2, result%gradient, error)
      if (error%has_error()) then
         call result%error%set(ERROR_VALIDATION, "gradient: "//error%get_message())
         result%has_error = .true.
         return
      end if
      result%has_gradient = .true.
   end subroutine mcscf_gradient_into

   subroutine store_decomposition(result, atom_energy, free_atom_energy, &
                                  pair_energy, pair_classical, formation_energy)
      !! Put the energy decomposition on the result, if there is one
      !!
      !! Written once and called from both analysis sites, because the two
      !! differ only in which wave function they analysed. Nothing happens when
      !! the decomposition did not run: `energy_decomposition` is off by
      !! default and the arrays come back unallocated, which is what tells
      !! "not asked for" from "came out zero".
      type(calculation_result_t), intent(inout) :: result
      real(dp), intent(in), allocatable :: atom_energy(:), free_atom_energy(:)
      real(dp), intent(in), allocatable :: pair_energy(:, :), pair_classical(:, :)
      real(dp), intent(in) :: formation_energy

      if (.not. allocated(atom_energy)) return
      result%ieda_atom = atom_energy
      if (allocated(free_atom_energy)) result%ieda_free_atom = free_atom_energy
      if (allocated(pair_energy)) result%ieda_pair = pair_energy
      if (allocated(pair_classical)) result%ieda_classical = pair_classical
      result%ieda_formation = formation_energy
      result%has_ieda = .true.
   end subroutine store_decomposition

   pure function to_real_text(value) result(out)
      !! A small real in exponent form, for an error message
      real(dp), intent(in) :: value
      character(len=:), allocatable :: out
      character(len=24) :: buffer
      write (buffer, "(es9.2)") value
      out = trim(adjustl(buffer))
   end function to_real_text

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
                                  ghost=ghost_of(fragment), &
                                  force_cartesian=settings%cartesian)
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
