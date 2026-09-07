!! One place where a `physical_fragment_t` becomes a terco SCF result
module mqc_terco_driver
   !! mqc builds the basis and the initial density on the CPU, hands both to
   !! terco, and terco runs every SCF iteration on the device. Nothing comes
   !! back until the SCF has converged, so this module is a marshalling layer
   !! and not a driver in the sense the cuEST one is: there is no iteration
   !! here to interleave with.
   !!
   !! **The guess is computed on the CPU on purpose.** terco's own fallback is
   !! the core Hamiltonian, which for anything past a few atoms costs more SCF
   !! iterations than the superposition-of-atomic-densities guess saves. The
   !! CPU path already has SAD, and a density is small next to the integrals,
   !! so it is cheaper to build it here and ship it than to teach terco how.
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, &
                                                                             c_null_ptr, c_null_char, c_loc
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_elements, only: element_number_to_symbol
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_ao, only: max_ao_l
   use mqc_czt_bridge, only: core_orbital_count
   use libcint_fortran, only: LIBCINT_NCTR_OF, LIBCINT_NPRIM_OF, LIBCINT_PTR_COEFF
   use mqc_czt_atomic_guess, only: parse_guess_name
   use mqc_czt_rhf, only: SCF_GUESS_CORE, SCF_GUESS_GWH, SCF_GUESS_SAC, SCF_GUESS_SAD
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_libcint, &
                               trc_set_aux_libcint, trc_set_method, trc_set_convergence, &
                               trc_set_screening, trc_set_guess, trc_set_verbose, trc_set_rimp2, &
                               trc_run_scf, trc_run_rimp2, trc_nao, trc_energy, trc_iterations, &
                               trc_rimp2_energy, trc_message, &
                               TRC_OK, TRC_ERR_NOCONV, TRC_ERR_UNSUPPORTED, &
                               TRC_GUESS_CORE, TRC_GUESS_GWH, TRC_GUESS_SAD
   implicit none
   private

   public :: run_terco_scf

   character(len=*), parameter :: TERCO_FUNCTIONALS = &
                                  " svwn pbe blyp pbe0 b3lyp b3lyp5 lda_x "
      !! What terco's XC integrator carries. Checked here, before anything is
      !! built, so a deck naming a functional terco has never heard of is told
      !! which ones exist rather than coming back as status 3 from the library.

contains

   subroutine run_terco_scf(settings, fragment, result, want_gradient)
      !! Run one fragment's SCF inside terco, on the device
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient

      type(czt_molecule_t) :: mol
      type(error_t) :: error
      character(len=8), allocatable :: element_symbols(:)
      character(len=:), allocatable :: functional_name
      character(kind=c_char), allocatable :: functional_c(:)
      integer, parameter :: MSG_LEN = 256   !! terco's message buffer
      character(kind=c_char) :: msg_c(MSG_LEN)
      type(c_ptr) :: h
      integer(c_int) :: rc, nao_c, niter_c, nbas_c, guess_c
      real(c_double) :: energy_c, e_os_c, e_ss_c
      integer(c_int), allocatable :: atm_c(:, :), bas_c(:, :)
      real(c_double), allocatable :: env_c(:)
      integer :: iatom, nspin, n_alpha, n_beta, guess_kind, l_max
      logical :: occupations_ok, unrestricted, wants_gradient

      wants_gradient = .false.
      if (present(want_gradient)) wants_gradient = want_gradient

      ! ---- what this backend cannot do, refused before anything is built ----
      !
      ! terco has energy kernels only. A gradient request must fail loudly:
      ! answering with the energy and no gradient would let a geometry
      ! optimisation run on a result that silently lacks what it asked for.
      if (wants_gradient) then
         call result%error%set(ERROR_VALIDATION, "backend 'terco' was asked for a "// &
                               "gradient, but terco has energy kernels only -- there "// &
                               "are no derivative kernels behind this entry. Ask for "// &
                               "backend 'libcint', or request an energy.")
         result%has_error = .true.
         result%has_energy = .false.
         return
      end if

      if (settings%run_cc) then
         call result%error%set(ERROR_VALIDATION, "backend 'terco' was asked for, but "// &
                               "coupled cluster has no terco implementation -- it runs "// &
                               "through the CPU backend. Ask for 'auto', or drop it.")
         result%has_error = .true.
         result%has_energy = .false.
         return
      end if
      if (settings%run_mp2 .and. .not. settings%corr_density_fitting) then
         call result%error%set(ERROR_VALIDATION, "backend 'terco' carries RI-MP2 only: "// &
                               "ask for method 'ri-mp2' with model.aux_basis, or for "// &
                               "backend 'auto' to run exact MP2 on the CPU.")
         result%has_error = .true.
         result%has_energy = .false.
         return
      end if

      functional_name = trim(adjustl(settings%functional))
      if (len(functional_name) > 0) then
         if (index(TERCO_FUNCTIONALS, " "//lowercase(functional_name)//" ") == 0) then
            call result%error%set(ERROR_VALIDATION, "backend 'terco' does not carry the "// &
                                  "functional '"//functional_name//"'. It has SVWN, PBE, "// &
                                  "BLYP, PBE0, B3LYP, B3LYP5 and LDA_X. Ask for backend "// &
                                  "'libcint', which goes through libxc.")
            result%has_error = .true.
            result%has_energy = .false.
            return
         end if
      end if

      ! ---- the molecule, in libcint's packing ------------------------------
      allocate (element_symbols(fragment%n_atoms))
      do iatom = 1, fragment%n_atoms
         element_symbols(iatom) = element_number_to_symbol(fragment%element_numbers(iatom))
      end do

      call build_czt_molecule(fragment%element_numbers, element_symbols, &
                              fragment%coordinates, settings%basis_set, mol, error, &
                              force_cartesian=settings%cartesian, &
                              ecp_name=settings%ecp_set)
      if (error%has_error()) then
         call result%error%set(ERROR_VALIDATION, "backend 'terco': "//error%get_message())
         result%has_error = .true.
         result%has_energy = .false.
         return
      end if

      ! Angular momentum, and the one case where Cartesian and spherical are
      ! the same basis. terco's four-centre kernels stop at d, and it reads
      ! every shell Cartesian. Through p the two forms span the same space with
      ! the same count -- three p functions either way -- so a spherical deck is
      ! safe there and only there. At d they differ, five functions against six,
      ! and running one as the other answers a different question in silence.
      l_max = max_ao_l(mol)
      if (l_max > 2) then
         call result%error%set(ERROR_VALIDATION, "backend 'terco' was asked for, but this "// &
                               "basis carries angular momentum above d and terco's "// &
                               "four-centre kernels stop at d. Ask for backend 'libcint', "// &
                               "or choose a smaller basis.")
         result%has_error = .true.
         result%has_energy = .false.
         return
      end if
      if (l_max >= 2 .and. .not. mol%cartesian) then
         call result%error%set(ERROR_VALIDATION, "backend 'terco' was asked for with a "// &
                               "spherical d shell. terco builds every shell Cartesian, so "// &
                               "it would answer with six d functions where the deck asked "// &
                               "for five -- a different basis, reported as though it were "// &
                               "the requested one. Set 'model.cartesian', or ask for "// &
                               "backend 'libcint'.")
         result%has_error = .true.
         result%has_energy = .false.
         return
      end if

      ! terco's basis entry reads orbital shells only; an ECP row handed to it
      ! would be taken for a basis shell with a nonsense angular momentum.
      if (mol%necpbas > 0) then
         call result%error%set(ERROR_VALIDATION, "backend 'terco' was asked for with an "// &
                               "effective core potential, which terco does not read. Ask "// &
                               "for backend 'libcint', or drop 'model.ecp'.")
         result%has_error = .true.
         result%has_energy = .false.
         return
      end if

      ! ---- occupations ------------------------------------------------------
      call spin_occupations(fragment%nelec, fragment%multiplicity, n_alpha, n_beta, &
                            occupations_ok)
      if (.not. occupations_ok) then
         call result%error%set(ERROR_VALIDATION, "backend 'terco': a multiplicity of "// &
                               "this fragment cannot be built from its electron count.")
         result%has_error = .true.
         result%has_energy = .false.
         return
      end if
      unrestricted = (n_alpha /= n_beta)
      nspin = 1
      if (unrestricted) nspin = 2

      ! ---- the initial density, built on the CPU ---------------------------
      !
      ! `guess_total` comes back unallocated for core and GWH, which terco is
      ! told by a null pointer: it then builds the core guess itself. An atomic
      ! guess that will not build has already warned and fallen back to GWH
      ! inside `build_restricted_guess`, so there is nothing to handle here.
      ! ---- which guess. terco builds SAD itself (one free-atom SCF per
      ! element, spin-averaged); core and GWH are its own; SAC is the same
      ! spherical atom under another name; a projection ladder is refused
      ! here as it was refused by the CPU guess builder.
      call parse_guess_name(settings%guess, guess_kind, error)
      if (error%has_error()) then
         call result%error%set(ERROR_VALIDATION, "backend 'terco': "//error%get_message())
         result%has_error = .true.
         result%has_energy = .false.
         return
      end if
      select case (guess_kind)
      case (SCF_GUESS_CORE)
         guess_c = TRC_GUESS_CORE
      case (SCF_GUESS_GWH)
         guess_c = TRC_GUESS_GWH
      case (SCF_GUESS_SAD, SCF_GUESS_SAC)
         guess_c = TRC_GUESS_SAD
      case default
         call result%error%set(ERROR_VALIDATION, "backend 'terco': guess '"// &
                               trim(settings%guess)//"' has no terco counterpart; use "// &
                               "'sad', 'core' or 'gwh'.")
         result%has_error = .true.
         result%has_energy = .false.
         return
      end select

      ! ---- one context, set up in steps, then run ---------------------------
      rc = trc_create(h)
      if (rc /= TRC_OK) then
         call fail("could not create a terco context (status "//status_text(rc)//")")
         return
      end if
      rc = trc_set_molecule(h, int(fragment%n_atoms, c_int), &
                            real(fragment%element_numbers, c_double), &
                            real(reshape(fragment%coordinates, [3*fragment%n_atoms]), c_double), &
                            int(fragment%charge, c_int), int(fragment%multiplicity, c_int))
      if (rc /= TRC_OK) then
         call fail("the molecule was refused by terco (status "//status_text(rc)//")")
         return
      end if
      call segment_general_contractions(mol%bas, mol%nbas, bas_c, nbas_c)
      rc = trc_set_basis_libcint(h, int(mol%atm, c_int), int(mol%natm, c_int), bas_c, nbas_c, &
                                 real(mol%env, c_double), int(size(mol%env), c_int), 1_c_int)
      if (rc /= TRC_OK) then
         call fail("the basis was refused by terco (status "//status_text(rc)//")")
         return
      end if
      rc = trc_nao(h, nao_c)
      if (rc /= TRC_OK .or. int(nao_c) /= mol%nao) then
         call fail("terco and the CPU path disagree about the size of this basis. "// &
                   "That is a bug in the adapter, not in the deck.")
         return
      end if
      functional_c = c_string(functional_name)
      rc = trc_set_method(h, functional_c, int(settings%grid_level, c_int))
      ! terco stops on the energy change and the commutator norm, not on a
      ! density change. The norm grows with the system, so mqc's per-element
      ! density tolerance is the wrong scale for it; the gate is the square
      ! root of the energy tolerance, where PySCF puts its gradient threshold.
      ! cholesterol at 1e-6 took 41 iterations, the last dozen chasing a
      ! norm that hovered just above the gate.
      if (rc == TRC_OK) rc = trc_set_convergence(h, real(settings%energy_tol, c_double), &
                                                 real(sqrt(settings%energy_tol), c_double), &
                                                 int(settings%max_iter, c_int), int(settings%diis_size, c_int))
      if (rc == TRC_OK) rc = trc_set_screening(h, real(settings%screening_tolerance, c_double))
      if (rc == TRC_OK) rc = trc_set_guess(h, guess_c, c_null_ptr, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_verbose(h, int(merge(1, 0, settings%verbose), c_int))
      if (rc /= TRC_OK) then
         call fail("a setting was refused by terco (status "//status_text(rc)//")")
         return
      end if

      rc = trc_run_scf(h)
      select case (rc)
      case (TRC_OK, TRC_ERR_NOCONV)
         niter_c = 0
         if (trc_iterations(h, niter_c) == TRC_OK) result%scf_iterations = int(niter_c)
         if (trc_energy(h, energy_c) /= TRC_OK) then
            call fail("terco ran the SCF and would not report its energy")
            return
         end if
         if (rc == TRC_OK) then
            result%energy%scf = real(energy_c, dp)
            result%has_energy = .true.
         else if (settings%allow_crap_scf) then
            result%energy%scf = real(energy_c, dp)
            result%has_energy = .true.
            call logger%warning("backend 'terco': the SCF did not converge in the "// &
                                "iterations allowed; keeping the last iterate because "// &
                                "'allow_crap_scf' is set.")
         else
            call fail("the SCF did not converge in the iterations allowed. Raise "// &
                      "'keywords.scf.max_iter', or set 'allow_crap_scf' to keep the last iterate.")
            return
         end if
      case default
         call fail("the SCF failed (status "//status_text(rc)//"): "//terco_message())
         return
      end select

      if (settings%run_mp2 .and. result%has_energy) then
         block
            type(czt_molecule_t) :: aux
            integer(c_int), allocatable :: aux_bas_c(:, :)
            integer(c_int) :: naux_bas_c
            integer :: frozen
            character(len=:), allocatable :: aux_name
            aux_name = trim(settings%aux_basis_set)
            if (len_trim(aux_name) == 0) then
               call fail("RI-MP2 needs an auxiliary basis: set model.aux_basis")
               return
            end if
            if (index(aux_name, "rifit") == 0 .and. index(aux_name, "-ri") == 0) then
               call logger%warning("bad basis set, calculations will be poor: '"//aux_name// &
                                   "' is not a correlation-fitting (RIFIT) set")
            end if
            call build_czt_molecule(fragment%element_numbers, element_symbols, &
                                    fragment%coordinates, aux_name, aux, error, &
                                    force_cartesian=settings%cartesian)
            if (error%has_error()) then
               call fail(error%get_message())
               return
            end if
            call segment_general_contractions(aux%bas, aux%nbas, aux_bas_c, naux_bas_c)
            rc = trc_set_aux_libcint(h, int(aux%atm, c_int), int(aux%natm, c_int), aux_bas_c, &
                                     naux_bas_c, real(aux%env, c_double), &
                                     int(size(aux%env), c_int), 1_c_int)
            call aux%destroy()
            if (rc /= TRC_OK) then
               call fail("the auxiliary basis was refused by terco (status "//status_text(rc)//")")
               return
            end if
            frozen = settings%n_frozen_core
            if (frozen < 0) frozen = core_orbital_count(fragment%element_numbers)
            if (.not. settings%freeze_core) frozen = 0
            rc = trc_set_rimp2(h, int(frozen, c_int), 256_c_int)
            if (rc == TRC_OK) rc = trc_run_rimp2(h)
            if (rc == TRC_OK) rc = trc_rimp2_energy(h, e_os_c, e_ss_c)
            if (rc /= TRC_OK) then
               call fail("RI-MP2 failed (status "//status_text(rc)//"): "//terco_message())
               return
            end if
            result%energy%mp2%ss = real(e_ss_c, dp)
            result%energy%mp2%os = real(e_os_c, dp)
            result%energy%mp2%ss_scale = settings%scs_ss
            result%energy%mp2%os_scale = settings%scs_os
            call result%energy%mp2%check_stability()
         end block
      end if
      rc = trc_destroy(h)

   contains

      subroutine fail(why)
         !! Refuse, and release the context on the way out.
         character(len=*), intent(in) :: why
         integer(c_int) :: rc_d
         call result%error%set(ERROR_VALIDATION, "backend 'terco': "//why)
         result%has_error = .true.
         result%has_energy = .false.
         rc_d = trc_destroy(h)
      end subroutine fail

      function terco_message() result(text)
         !! terco's last message, or an empty string.
         character(len=:), allocatable :: text
         integer :: i
         text = ""
         if (trc_message(h, msg_c, int(MSG_LEN, c_int)) /= TRC_OK) return
         do i = 1, MSG_LEN
            if (msg_c(i) == c_null_char) exit
            text = text//msg_c(i)
         end do
      end function terco_message

   end subroutine run_terco_scf

   pure subroutine segment_general_contractions(bas, nbas, bas_out, nbas_out)
      !! One `bas` row per contracted function, which is all terco can read
      !!
      !! A generally contracted shell -- several contracted functions sharing one
      !! set of primitives, which is how cc-pVDZ writes its s space -- is stored
      !! by libcint as an `(nprim, nctr)` coefficient matrix with stride nprim,
      !! `coeff[nprim*ictr + iprim]`. terco reads `NPRIM_OF` coefficients from
      !! `PTR_COEFF` and has no notion of `NCTR_OF`, so handing it such a row
      !! builds ONE function from the first column and drops the rest: a smaller
      !! basis than the deck names, reported as the deck's own.
      !!
      !! The columns are already independent functions over shared exponents, so
      !! the fix is bookkeeping rather than arithmetic. Each column becomes its
      !! own shell: same `PTR_EXP`, `NCTR_OF` of one, and `PTR_COEFF` advanced to
      !! that column. `env` is untouched -- the split rows point into the
      !! coefficient block that is already there.
      !!
      !! **AO order is preserved**, which is what makes this safe to do behind
      !! the caller's back. libcint numbers a shell's functions with the
      !! contraction index slowest -- `libcint_cgto_cart` is `nctr*ncart` -- so
      !! consecutive split shells occupy exactly the positions the one general
      !! shell did, and a density built against `mol` still lines up.
      !!
      !! What it costs is the reason general contraction exists: the primitive
      !! work that the columns used to share is now repeated per column, and the
      !! shell-pair count grows with it. Correct and slower beats refused.
      integer, intent(in) :: bas(:, :)
      integer, intent(in) :: nbas
      integer(c_int), allocatable, intent(out) :: bas_out(:, :)
      integer(c_int), intent(out) :: nbas_out

      integer :: ish, ictr, nprim, nctr, row

      nbas_out = 0
      do ish = 1, nbas
         nbas_out = nbas_out + int(max(bas(LIBCINT_NCTR_OF, ish), 1), c_int)
      end do

      allocate (bas_out(size(bas, 1), nbas_out))
      row = 0
      do ish = 1, nbas
         nprim = bas(LIBCINT_NPRIM_OF, ish)
         nctr = max(bas(LIBCINT_NCTR_OF, ish), 1)
         do ictr = 1, nctr
            row = row + 1
            bas_out(:, row) = int(bas(:, ish), c_int)
            bas_out(LIBCINT_NCTR_OF, row) = 1_c_int
            bas_out(LIBCINT_PTR_COEFF, row) = &
               int(bas(LIBCINT_PTR_COEFF, ish) + (ictr - 1)*nprim, c_int)
         end do
      end do
   end subroutine segment_general_contractions

   pure subroutine spin_occupations(n_electrons, multiplicity, n_alpha, n_beta, ok)
      !! Alpha and beta counts from an electron count and a multiplicity
      !!
      !! Written here rather than taken from the cuEST backend: this backend
      !! must build in a tree where cuEST is switched off, and a shared eight
      !! lines is a smaller price than a dependency on a GPU module.
      integer, intent(in) :: n_electrons, multiplicity
      integer, intent(out) :: n_alpha, n_beta
      logical, intent(out) :: ok

      integer :: unpaired

      n_alpha = 0
      n_beta = 0
      ok = .false.
      if (multiplicity < 1 .or. n_electrons < 0) return

      unpaired = multiplicity - 1
      if (mod(n_electrons - unpaired, 2) /= 0) return
      if (n_electrons < unpaired) return

      n_beta = (n_electrons - unpaired)/2
      n_alpha = n_beta + unpaired
      ok = .true.
   end subroutine spin_occupations

   pure function c_string(text) result(buffer)
      !! A Fortran string as a NUL-terminated C one
      character(len=*), intent(in) :: text
      character(kind=c_char), allocatable :: buffer(:)

      integer :: i

      allocate (buffer(len(text) + 1))
      do i = 1, len(text)
         buffer(i) = text(i:i)
      end do
      buffer(len(text) + 1) = c_null_char
   end function c_string

   pure function lowercase(text) result(lowered)
      !! ASCII lower case, for a case-insensitive name comparison
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lowered

      integer :: i, code

      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar("A") .and. code <= iachar("Z")) then
            lowered(i:i) = achar(code + 32)
         else
            lowered(i:i) = text(i:i)
         end if
      end do
   end function lowercase

   pure function status_text(rc) result(text)
      !! A terco status code as something a user can act on
      integer(c_int), intent(in) :: rc
      character(len=:), allocatable :: text

      select case (rc)
      case (TRC_ERR_UNSUPPORTED)
         text = "unsupported"
      case (TRC_ERR_NOCONV)
         text = "not converged"
      case default
         text = "internal"
      end select
   end function status_text

end module mqc_terco_driver
