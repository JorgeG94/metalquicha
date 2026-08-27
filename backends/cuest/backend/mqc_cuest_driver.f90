!! Fragment-level entry point for cuEST-backed SCF calculations
module mqc_cuest_driver
   !! One place where a `physical_fragment_t` becomes a converged SCF result.
   !!
   !! Hartree-Fock and Kohn-Sham differ only in whether a functional is named,
   !! so both methods funnel through here rather than each carrying its own
   !! copy of the basis loading, context acquisition and teardown.
   use pic_types, only: dp
   use mqc_string_utils, only: int_to_text
   use mqc_error, only: error_t, ERROR_VALIDATION, ERROR_GENERIC
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_cgto, only: molecular_basis_type
   use mqc_elements, only: element_number_to_symbol
   use mqc_basis_utils, only: find_basis_file
   use mqc_json_basis_reader, only: build_molecular_basis_json
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t, SCF_CONVERGED, SCF_NOT_CONVERGED, &
                               scf_not_converged_message, &
                               frontier_orbitals
   use mqc_cuest_context, only: cuest_context_t, get_cuest_context
   use mqc_cuest_integrals, only: cuest_system_t
   use mqc_cuest_functionals, only: functional_name_to_id
   use mqc_cuest_scf, only: scf_result_t, run_rhf_scf, run_uks_scf, spin_occupations, &
                            SCF_GUESS_CORE, SCF_GUESS_GWH, SCF_GUESS_SAC
   use mqc_cuest_atomic_guess, only: build_sac_guess, atom_ao_counts
   use mqc_population_analysis, only: ao_owner_from_counts, mulliken_atomic_charges, &
                                      mulliken_atomic_spin_populations
   use mqc_cuest_gradient, only: compute_scf_gradient
   use mqc_cuest_iface, only: cuest_scf_settings_t
   implicit none
   private

   public :: run_cuest_scf         !! Fragment -> converged result

! cuest_scf_settings_t now lives in mqc_cuest_iface (src/), so the method
   ! files can reach it without pulling the backend into the fpm build.

contains

   subroutine run_cuest_scf(settings, fragment, result, want_gradient)
      !! Run a closed-shell cuEST SCF for one fragment
      !!
      !! With `want_gradient`, the analytic nuclear gradient is evaluated from
      !! the converged density before the cuEST objects are torn down -- they
      !! describe this geometry and this basis, so reusing them is both
      !! cheaper and safer than rebuilding.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient

      type(cuest_context_t), pointer :: context
      type(cuest_system_t) :: system
      type(scf_result_t) :: scf
      type(molecular_basis_type) :: orbital_basis, auxiliary_basis
      type(error_t) :: error
      real(dp), allocatable :: gradient(:, :)
      real(dp), allocatable :: guess_alpha(:, :), guess_beta(:, :)
      type(molecular_basis_type), allocatable :: atom_bases(:), atom_aux_bases(:)
      integer :: guess_type, n_guess_columns
      character(len=8), allocatable :: element_symbols(:)
      integer :: iatom, functional_id, n_alpha, n_beta
      logical :: need_gradient, unrestricted, occupations_ok

      need_gradient = .false.
      if (present(want_gradient)) need_gradient = want_gradient

      ! ---- which functional, if any ----------------------------------------
      if (len_trim(settings%functional) == 0) then
         functional_id = -1   ! pure Hartree-Fock, no XC plan
      else
         call functional_name_to_id(settings%functional, functional_id, error)
         if (error%has_error()) then
            call record_failure(result, error)
            return
         end if
      end if

      ! ---- basis sets -------------------------------------------------------
      allocate (element_symbols(fragment%n_atoms))
      do iatom = 1, fragment%n_atoms
         element_symbols(iatom) = element_number_to_symbol(fragment%element_numbers(iatom))
      end do

      call load_basis(settings%basis_set, element_symbols, orbital_basis, error, "orbital")
      if (error%has_error()) then
         call record_failure(result, error)
         return
      end if

      call load_basis(settings%aux_basis_set, element_symbols, auxiliary_basis, error, "auxiliary")
      if (error%has_error()) then
         call record_failure(result, error)
         return
      end if

      ! ---- per-rank handle and device scratch -------------------------------
      !
      ! One handle and one set of scratch pools per rank, reused by every
      ! fragment this rank evaluates. The node-local rank spreads ranks across
      ! the GPUs on the node; it is ignored after the first call, since the
      ! context is created once and then shared.
      call get_cuest_context(settings%device_rank, context, error)
      if (error%has_error()) then
         call record_failure(result, error)
         return
      end if

      ! ---- which initial guess? ---------------------------------------------
      select case (trim(settings%guess))
      case ("core")
         guess_type = SCF_GUESS_CORE
      case ("sac")
         guess_type = SCF_GUESS_SAC
      case default
         guess_type = SCF_GUESS_GWH
      end select

      ! ---- restricted or unrestricted? --------------------------------------
      !
      ! An odd electron count forces unrestricted regardless of what the input
      ! claims, so a mis-stated multiplicity fails as a validation error rather
      ! than silently halving an odd number of electrons.
      call spin_occupations(fragment%nelec, fragment%multiplicity, n_alpha, n_beta, occupations_ok)
      if (.not. occupations_ok) then
         call error%set(ERROR_VALIDATION, "Electron count and multiplicity are inconsistent")
         call record_failure(result, error)
         return
      end if
      unrestricted = (fragment%multiplicity /= 1) .or. (mod(fragment%nelec, 2) /= 0) &
                     .or. settings%unrestricted

      ! ---- superposed atomic guess, if asked for ----------------------------
      !
      ! Built before the molecular system so its width is known: the guess
      ! carries the summed atomic occupations, usually more columns than the
      ! molecule has electrons, and the device pools must be sized for that up
      ! front rather than grown later under a live system.
      n_guess_columns = 0
      if (guess_type == SCF_GUESS_SAC) then
         call build_atom_bases(settings, element_symbols, atom_bases, atom_aux_bases, error)
         if (.not. error%has_error()) then
            call build_sac_guess(context, fragment%element_numbers, orbital_basis, &
                                 auxiliary_basis, settings%spherical, functional_id, &
                                 settings%radial_points, settings%angular_points, &
                                 atom_bases, atom_aux_bases, guess_alpha, guess_beta, error)
         end if
         if (error%has_error()) then
            call record_failure(result, error)
            return
         end if
         n_guess_columns = size(guess_alpha, 2) + size(guess_beta, 2)
      end if

      ! ---- build, solve, tear down ------------------------------------------
      if (unrestricted) then
         call system%create(context, fragment%element_numbers, fragment%coordinates, &
                            orbital_basis, auxiliary_basis, settings%spherical, &
                            n_alpha, functional_id, settings%radial_points, &
                            settings%angular_points, error, n_occ_beta=n_beta, &
                            n_guess_columns=n_guess_columns, pcm=settings%pcm)
      else
         call system%create(context, fragment%element_numbers, fragment%coordinates, &
                            orbital_basis, auxiliary_basis, settings%spherical, &
                            fragment%nelec/2, functional_id, settings%radial_points, &
                            settings%angular_points, error, &
                            n_guess_columns=n_guess_columns, pcm=settings%pcm)
      end if

      if (.not. error%has_error()) then
         if (unrestricted) then
            if (guess_type == SCF_GUESS_SAC) then
               call run_uks_scf(system, context, fragment%element_numbers, fragment%coordinates, &
                                fragment%nelec, fragment%multiplicity, settings%max_iter, &
                                settings%energy_tol, settings%density_tol, settings%use_diis, &
                                settings%diis_size, settings%verbose, scf, error, &
                                guess=guess_type, guess_alpha=guess_alpha, guess_beta=guess_beta)
            else
               call run_uks_scf(system, context, fragment%element_numbers, fragment%coordinates, &
                                fragment%nelec, fragment%multiplicity, settings%max_iter, &
                                settings%energy_tol, settings%density_tol, settings%use_diis, &
                                settings%diis_size, settings%verbose, scf, error, guess=guess_type)
            end if
         else
            if (guess_type == SCF_GUESS_SAC) then
               call run_rhf_scf(system, context, fragment%element_numbers, fragment%coordinates, &
                                fragment%nelec, settings%max_iter, settings%energy_tol, &
                                settings%density_tol, settings%use_diis, settings%diis_size, &
                                settings%verbose, scf, error, guess=guess_type, &
                                guess_alpha=guess_alpha, guess_beta=guess_beta)
            else
               call run_rhf_scf(system, context, fragment%element_numbers, fragment%coordinates, &
                                fragment%nelec, settings%max_iter, settings%energy_tol, &
                                settings%density_tol, settings%use_diis, settings%diis_size, &
                                settings%verbose, scf, error, guess=guess_type)
            end if
         end if
      end if

      if (need_gradient .and. .not. error%has_error()) then
         call compute_scf_gradient(system, fragment%element_numbers, fragment%coordinates, &
                                   scf, gradient, error)
      end if

      call system%destroy()
      call orbital_basis%destroy()
      call auxiliary_basis%destroy()

      if (error%has_error()) then
         call record_failure(result, error)
         return
      end if

      result%energy%scf = scf%total_energy
      result%has_energy = .true.

      ! The SCF has known this all along; it simply was not carried out of
      ! here. Without it a fragment that ran out of iterations contributes its
      ! last-cycle energy to the expansion, has_error stays false, and the
      ! total is wrong by however far that SCF still had to go -- with nothing
      ! anywhere saying so.
      if (scf%converged) then
         result%scf_status = SCF_CONVERGED
      else
         result%scf_status = SCF_NOT_CONVERGED
      end if
      result%scf_iterations = scf%iterations

      ! cuEST reports eigenvalues but not occupations, so they are rebuilt
      ! from the electron count: doubly occupied up to nelec/2 for a closed
      ! shell, singly up to n_alpha when unrestricted. Same routine as the
      ! xTB path picks the pair, so one definition of "frontier" serves both.
      if (allocated(scf%orbital_energies)) then
         block
            real(dp), allocatable :: occupations(:)
            integer :: n_orb, n_filled, iorb

            n_orb = size(scf%orbital_energies)
            allocate (occupations(n_orb))
            occupations = 0.0_dp
            if (unrestricted) then
               n_filled = (fragment%nelec + (fragment%multiplicity - 1))/2
               do iorb = 1, min(n_filled, n_orb)
                  occupations(iorb) = 1.0_dp
               end do
            else
               n_filled = fragment%nelec/2
               do iorb = 1, min(n_filled, n_orb)
                  occupations(iorb) = 2.0_dp
               end do
            end if
            call frontier_orbitals(scf%orbital_energies, occupations, &
                                   result%homo, result%lumo, result%has_orbitals)
         end block
      end if

      ! Same policy as xTB, deliberately. One physical condition -- an SCF that
      ! ran out of iterations -- must not stop the job on one backend and pass
      ! silently on the other, which is what it did while this was only
      ! recorded here.
      if (.not. scf%converged .and. .not. settings%allow_crap_scf) then
         call error%set(ERROR_GENERIC, scf_not_converged_message(scf%iterations))
         call record_failure(result, error)
         return
      end if

      ! The dipole is what IR intensities are built from: the Hessian path
      ! collects it at every displacement and differences it.
      if (scf%has_dipole) then
         result%dipole = scf%dipole
         result%has_dipole = .true.
      end if

      ! Atomic charges, from the density this calculation converged to.
      !
      ! Done here rather than in the SCF because `system` is still alive and
      ! owns the overlap, and done before the gradient for no reason other than
      ! that both need the same live objects.
      !
      ! The arithmetic is `mqc_population_analysis`, the same code the CPU path
      ! runs. Only the AO-to-atom map is built differently: cuEST lays its AO
      ! basis out per atom in atom order, so a count per atom is the whole of
      ! it.
      if (allocated(settings%charges_scheme)) then
         block
            real(dp), allocatable :: q(:), q_spin(:), s_ao(:, :), spin_density(:, :)
            integer, allocatable :: counts(:), owner(:)
            type(error_t) :: analysis_error
            logical :: open_shell

            call analysis_error%clear()
            open_shell = allocated(scf%density_beta)

            select case (trim(settings%charges_scheme))
            case ("mulliken")
               allocate (counts(size(fragment%element_numbers)))
               call atom_ao_counts(orbital_basis,.not. orbital_basis%is_cartesian(), counts)
               call ao_owner_from_counts(counts, owner)
               call system%compute_overlap(s_ao, analysis_error)

               if (.not. analysis_error%has_error()) then
                  ! `scf%density` is already the total on this backend -- the
                  ! unrestricted path stores density_a + density_b there, where
                  ! the CPU path stores alpha alone. So the spin density is
                  ! total - 2*beta, and not the difference of the two stored
                  ! matrices.
                  call mulliken_atomic_charges(owner, real(fragment%element_numbers, dp), &
                                               scf%density, s_ao, q, analysis_error)
               end if
               if (open_shell .and. .not. analysis_error%has_error()) then
                  spin_density = scf%density - 2.0_dp*scf%density_beta
                  call mulliken_atomic_spin_populations(owner, size(fragment%element_numbers), &
                                                        spin_density, s_ao, q_spin, &
                                                        analysis_error)
               end if
            case ("chelpg")
               ! Not a refusal that can be lifted by wiring one more call: the
               ! fit needs the electrostatic potential at points off the atoms,
               ! which this backend has no integral for. Named rather than
               ! silently answered with Mulliken, because the two schemes
               ! disagree by design.
               call analysis_error%set(ERROR_GENERIC, "CHELPG charges are not implemented "// &
                                       "on the GPU backend; it builds no electrostatic "// &
                                       "potential integrals. Use 'mulliken' here, or run "// &
                                       "CHELPG on the CPU backend.")
            case default
               call analysis_error%set(ERROR_GENERIC, "unknown charge scheme '"// &
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

      if (need_gradient) then
         call move_alloc(gradient, result%gradient)
         result%has_gradient = .true.
      end if
   end subroutine run_cuest_scf

   subroutine build_atom_bases(settings, element_symbols, atom_bases, atom_aux_bases, error)
      !! One-atom orbital and auxiliary bases, for the free-atom guess runs
      type(cuest_scf_settings_t), intent(in) :: settings
      character(len=*), intent(in) :: element_symbols(:)
      type(molecular_basis_type), allocatable, intent(out) :: atom_bases(:), atom_aux_bases(:)
      type(error_t), intent(inout) :: error

      integer :: iatom, n_atoms

      n_atoms = size(element_symbols)
      allocate (atom_bases(n_atoms), atom_aux_bases(n_atoms))
      do iatom = 1, n_atoms
         call load_basis(settings%basis_set, element_symbols(iatom:iatom), &
                         atom_bases(iatom), error)
         if (error%has_error()) return
         call load_basis(settings%aux_basis_set, element_symbols(iatom:iatom), &
                         atom_aux_bases(iatom), error)
         if (error%has_error()) return
      end do
   end subroutine build_atom_bases

   subroutine load_basis(basis_name, element_symbols, mol_basis, error, role)
      !! Locate and parse a basis set for the atoms of a fragment
      !!
      !! Basis Set Exchange JSON is the only format; `find_basis_file` walks
      !! the search path and reports where it looked if nothing turns up.
      !!
      !! A Cartesian basis is refused here, and that one is not a gap waiting to
      !! be filled the way the CPU backend's was. cuEST builds its AO shells
      !! pure, every plan and every device buffer is sized from the spherical
      !! function count, and there is no per-call convention to route on the way
      !! libcint has. Reading 6-31G* as spherical anyway is precisely the silent
      !! wrong answer this was all fixed for, so the deck is told to change the
      !! basis instead.
      character(len=*), intent(in) :: basis_name
      character(len=*), intent(in) :: element_symbols(:)
      type(molecular_basis_type), intent(out) :: mol_basis
      type(error_t), intent(out) :: error
      character(len=*), intent(in), optional :: role  !! "orbital" or "auxiliary", for diagnostics

      character(len=:), allocatable :: basis_path, what

      if (len_trim(basis_name) == 0) then
         call error%set(ERROR_VALIDATION, "No basis set specified")
         return
      end if

      call find_basis_file(basis_name, basis_path, error)
      if (error%has_error()) return

      call build_molecular_basis_json(basis_path, element_symbols, mol_basis, error)
      if (error%has_error()) return

      if (mol_basis%is_cartesian()) then
         what = "basis set"
         if (present(role)) what = trim(role)//" basis set"
         call error%set(ERROR_VALIDATION, "the "//what//" '"//trim(basis_name)// &
                        "' is Cartesian -- "//trim(basis_path)//" marks shells "// &
                        "'gto_cartesian', which is six d functions rather than five. "// &
                        "cuEST builds spherical shells only, so this basis cannot be "// &
                        "run on the GPU backend; use a spherical set such as cc-pVDZ "// &
                        "or def2-SVP, or run it on the CPU backend, which routes both.")
         return
      end if

      call check_basis_covers_atoms(basis_name, basis_path, element_symbols, mol_basis, error, role)
   end subroutine load_basis

   subroutine check_basis_covers_atoms(basis_name, basis_path, element_symbols, mol_basis, error, role)
      !! Fail unless every atom came back with at least one shell
      !!
      !! The readers already refuse an element that is absent from the file,
      !! but each of the three has its own notion of "absent" and a block that
      !! parses to nothing would slip through all of them. An atom with no
      !! basis functions is not a small error: cuEST would build a system whose
      !! AO space simply omits that centre and converge an SCF for a different
      !! molecule than the one that was asked for. Refuse it here, once, on the
      !! only path any of the readers can reach the backend by.
      character(len=*), intent(in) :: basis_name, basis_path
      character(len=*), intent(in) :: element_symbols(:)
      type(molecular_basis_type), intent(in) :: mol_basis
      type(error_t), intent(inout) :: error
      character(len=*), intent(in), optional :: role

      character(len=:), allocatable :: what
      integer :: iatom

      what = "basis set"
      if (present(role)) what = trim(role)//" basis set"

      if (mol_basis%nelements /= size(element_symbols)) then
         call error%set(ERROR_VALIDATION, "The "//what//" "//trim(basis_name)// &
                        " produced entries for a different number of atoms than the fragment has")
         return
      end if

      do iatom = 1, size(element_symbols)
         if (mol_basis%elements(iatom)%nshells <= 0) then
            call error%set(ERROR_VALIDATION, "No "//what//" defined for element "// &
                           trim(adjustl(element_symbols(iatom)))//" (atom "//int_to_text(iatom)// &
                           ") in "//trim(basis_path)//". Add the element to the basis file or "// &
                           "choose a basis that covers it.")
            return
         end if
      end do
   end subroutine check_basis_covers_atoms

   subroutine record_failure(result, error)
      !! Mark a calculation as failed, carrying the diagnostic with it
      type(calculation_result_t), intent(inout) :: result
      type(error_t), intent(in) :: error

      result%error = error
      result%has_error = .true.
      result%has_energy = .false.
   end subroutine record_failure

end module mqc_cuest_driver
