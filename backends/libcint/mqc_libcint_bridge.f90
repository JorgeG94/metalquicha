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
                              SCF_GUESS_GWH, SCF_GUESS_SAC, SCF_GUESS_SAD
   use mqc_libcint_atomic_guess, only: build_atomic_guess, parse_guess_name, &
                                       guess_display_name
   use mqc_libcint_mp2, only: mp2_result_t, run_libcint_mp2, run_libcint_ri_mp2
   use mqc_libcint_cc, only: cc_result_t, run_libcint_ccsd
   implicit none
   private

   public :: run_libcint_hf
   public :: libcint_backend_available

contains

   pure function libcint_backend_available() result(available)
      !! Whether this build can run an SCF on the CPU
      logical :: available
      available = .true.
   end function libcint_backend_available

   subroutine run_libcint_hf(settings, fragment, result, want_gradient)
      !! Closed-shell HF for one fragment, on the CPU
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient

      type(libcint_molecule_t) :: mol, aux, corr_aux
      type(rhf_result_t) :: scf
      type(error_t) :: error
      character(len=MAX_ELEMENT_SYMBOL_LEN), allocatable :: symbols(:)
      integer :: iatom, diis_size, guess_kind
      logical :: unrestricted
      ! Left unallocated for the guesses that need no free-atom solutions. An
      ! unallocated allocatable passed to an optional dummy is an absent
      ! argument, so the SCF calls below need no branching on which guess ran.
      real(dp), allocatable :: guess_a(:, :), guess_b(:, :), guess_total(:, :)
      type(error_t) :: guess_error

      if (present(want_gradient)) then
         if (want_gradient) then
            ! libcint has the derivative entry points -- ip1 and ipovlp -- so
            ! this is a gap rather than a wall, but an untested gradient is
            ! worse than an absent one.
            call result%error%set(ERROR_VALIDATION, "the CPU backend has no gradients "// &
                                  "yet; run an energy, or build with cuEST")
            result%has_error = .true.
            return
         end if
      end if

      ! Same rule the GPU backend applies, so a deck does not change meaning
      ! when it moves between them: an odd electron count or a multiplicity
      ! above one has no restricted solution to find, whatever the keyword says.
      unrestricted = (fragment%multiplicity /= 1) .or. (mod(fragment%nelec, 2) /= 0) &
                     .or. settings%unrestricted

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
                               "a closed-shell system, or use MP2.")
         result%has_error = .true.
         return
      end if
      if (settings%run_cc .and. settings%corr_density_fitting) then
         call result%error%set(ERROR_VALIDATION, "RI coupled cluster is not implemented "// &
                               "on the CPU backend. Ask for 'ccsd' or 'ccsd(t)' rather "// &
                               "than 'ri-ccsd', or set "// &
                               "keywords.correlation.density_fitting to false.")
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

      ! ---- which initial guess? ---------------------------------------------
      call parse_guess_name(settings%guess, guess_kind, error)
      if (error%has_error()) then
         call result%error%set(ERROR_VALIDATION, error%get_message())
         result%has_error = .true.
         call mol%destroy()
         return
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
                              guess_density=guess_total)
         call aux%destroy()
      else if (unrestricted) then
         call run_libcint_uhf(mol, fragment%nelec, fragment%multiplicity, settings%max_iter, &
                              settings%energy_tol, settings%density_tol, settings%verbose, &
                              scf, error, diis_vectors=diis_size, guess=guess_kind, &
                              guess_density_alpha=guess_a, guess_density_beta=guess_b)
      else
         call run_libcint_rhf(mol, fragment%nelec, settings%max_iter, settings%energy_tol, &
                              settings%density_tol, settings%verbose, scf, error, &
                              diis_vectors=diis_size, guess=guess_kind, &
                              guess_density=guess_total)
      end if
      if (error%has_error()) then
         call result%error%set(ERROR_VALIDATION, error%get_message())
         result%has_error = .true.
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
         call mol%destroy()
         return
      end if

      result%energy%scf = scf%energy
      result%has_energy = .true.

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
            if (settings%verbose) then
               write (*, "(a,i0,a,i0,a,i0)") "  MP2: frozen ", mp2%n_frozen, &
                  "  occupied ", mp2%n_occupied, "  virtual ", mp2%n_virtual
               write (*, "(a,f20.12)") "  E(OS)          ", mp2%opposite_spin
               write (*, "(a,f20.12)") "  E(SS)          ", mp2%same_spin
               write (*, "(a,f20.12)") "  E(corr)        ", mp2%correlation
               if (settings%scs_ss /= 1.0_dp .or. settings%scs_os /= 1.0_dp) then
                  write (*, "(a,f6.3,a,f6.3)") "  spin scaling:  SS x", settings%scs_ss, &
                     "   OS x", settings%scs_os
                  write (*, "(a,f20.12)") "  E(corr) scaled ", result%energy%mp2%total()
               end if
               write (*, "(a,f20.12)") "  E total        ", result%energy%total()
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

            call run_libcint_ccsd(mol, scf%orbitals, scf%orbital_energies, &
                                  fragment%nelec/2, frozen, settings%cc_max_iter, &
                                  settings%cc_tolerance, settings%cc_triples, &
                                  settings%verbose, cc, error, &
                                  diis_vectors=settings%cc_diis_size)
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
               write (*, "(a,i0,a,l1)") "  CCSD: iterations ", cc%iterations, &
                  "  converged ", cc%converged
               write (*, "(a,f20.12)") "  E(MP2)         ", cc%e_mp2
               write (*, "(a,f20.12)") "  E(singles)     ", cc%e_singles
               write (*, "(a,f20.12)") "  E(doubles)     ", cc%e_doubles
               if (settings%cc_triples) then
                  write (*, "(a,f20.12)") "  E(T)           ", cc%e_triples
               end if
               write (*, "(a,f20.12)") "  E(corr)        ", result%energy%cc%total()
               write (*, "(a,f20.12)") "  E total        ", result%energy%total()
            end if
         end block
      end if

      ! The frontier pair, from the orbital energies the SCF already produced.
      if (allocated(scf%orbital_energies)) then
         if (scf%n_occupied >= 1 .and. scf%n_occupied < size(scf%orbital_energies)) then
            result%homo = scf%orbital_energies(scf%n_occupied)
            result%lumo = scf%orbital_energies(scf%n_occupied + 1)
            result%has_orbitals = .true.
         end if
      end if

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
      !! `keywords.correlation.aux_basis` is what this reads. It falls back to
      !! `model.aux_basis` so a deck that already names one for a fitted SCF
      !! does not have to repeat it -- and that fallback is exactly the case
      !! the warning is for, since an SCF auxiliary is usually JKFIT.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      character(len=*), intent(in) :: symbols(:)
      type(libcint_molecule_t), intent(out) :: aux
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: name

      name = trim(settings%corr_aux_basis)
      if (len_trim(name) == 0) name = trim(settings%aux_basis_set)
      if (len_trim(name) == 0) then
         call error%set(ERROR_VALIDATION, "density-fitted correlation needs an "// &
                        "auxiliary basis: set keywords.correlation.aux_basis")
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
