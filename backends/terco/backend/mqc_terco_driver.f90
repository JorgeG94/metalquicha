!! One place where a `physical_fragment_t` becomes a terco calculation
module mqc_terco_driver
   !! mqc builds the molecule and basis on the CPU, in libcint's packing, and
   !! terco runs every SCF iteration, the RI-MP2 and the continuum on the
   !! device. Nothing comes back until a step has finished, so this module
   !! marshals rather than drives: there is no iteration here to interleave
   !! with.
   !!
   !! One terco context is opened per fragment and can run more than one SCF.
   !! A later state changes only the charge and multiplicity and keeps the
   !! bases, the fit, the settings and the cavity, so several electron counts
   !! at one geometry are several states of one session rather than several
   !! sessions.
   !!
   !! Densities come back in the AO order and normalisation of the
   !! `czt_molecule_t` built here, so the population analyses run on them with
   !! the code the CPU backend uses.
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_char, &
                                                                             c_null_ptr, c_null_char, c_loc, c_associated
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_pcm_radii, only: cavity_radius
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t, frontier_orbitals
   use mqc_elements, only: element_number_to_symbol
   use mqc_xc_spec, only: xc_spec_t, xc_spec_from_name
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_ao, only: max_ao_l
   use mqc_czt_bridge, only: core_orbital_count
   use libcint_fortran, only: LIBCINT_NCTR_OF, LIBCINT_NPRIM_OF, LIBCINT_PTR_COEFF
   use mqc_czt_atomic_guess, only: parse_guess_name, build_atomic_guess
   use mqc_czt_rhf, only: SCF_GUESS_CORE, SCF_GUESS_GWH, SCF_GUESS_SAC, SCF_GUESS_SAD, &
                          SCF_GUESS_SAP, SCF_GUESS_MINAO
   use mqc_czt_charges, only: mulliken_charges, mulliken_spin_populations, chelpg_charges
   use mqc_czt_bond_orders, only: mayer_bond_orders, mayer_bond_orders_open_shell, &
                                  mayer_valences, print_mayer_report
   use mqc_czt_fukui, only: fukui_result_t, fukui_from_states, print_fukui_report
   use mqc_scf_common, only: build_density_spin
   use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_libcint, &
                               trc_set_aux_libcint, trc_set_method, trc_set_convergence, &
                               trc_set_screening, trc_set_incremental, trc_set_guess, &
                               trc_set_verbose, trc_set_rimp2, trc_set_nlc_grid, trc_set_pcm, &
                               trc_run_scf, trc_run_rimp2, trc_nao, trc_energy, trc_iterations, &
                               trc_converged, trc_density, trc_mo_energies, trc_rimp2_energy, &
                               trc_mo_coefficients, &
                               trc_pcm_energy, trc_max_angular_momentum, trc_cart_order, &
                               trc_set_df, trc_set_jkfit_libcint, trc_run_gradient, trc_gradient, &
                               trc_set_charge, trc_set_unrestricted, trc_set_level_shift, &
                               trc_set_pcm_radii, trc_abi_version, trc_functional_info, &
                               trc_functional_list, trc_nmo, trc_mo_occupations, trc_message, &
                               TRC_ABI_MAJOR, TRC_ABI_MINOR, &
                               TRC_OK, TRC_ERR_NOCONV, TRC_ERR_UNSUPPORTED, &
                               TRC_GUESS_CORE, TRC_GUESS_GWH, TRC_GUESS_SAD, TRC_GUESS_GIVEN, &
                               TRC_GUESS_SAP
   implicit none
   private

   public :: run_terco_scf

   integer, parameter :: MSG_LEN = 512   !! terco's message buffer

   real(dp), parameter :: DH_WEIGHT_AGREEMENT = 1.0e-12_dp
      !! How far terco's double-hybrid PT2 weights may sit from mqc's own
      !! table. Both are typed-in constants, so anything above rounding is a
      !! functional the two tables define differently.

   integer(c_int), parameter :: FUNCTIONAL_RUNNABLE = 1
      !! `trc_functional_info` flag: terco can run this functional today
   integer(c_int), parameter :: FUNCTIONAL_DOUBLE_HYBRID = 32
      !! `trc_functional_info` flag: an RI-MP2 term is part of the functional

   type :: terco_session_t
      !! One terco context and the CPU-side molecule it was built from
      type(c_ptr) :: handle = c_null_ptr
      type(czt_molecule_t) :: mol
         !! The orbital basis, in the AO order of every density terco returns.
      character(len=8), allocatable :: symbols(:)
      integer :: nao = 0
      real(dp) :: c_os = 0.0_dp
         !! A double hybrid's weight on opposite-spin RI-MP2, zero otherwise.
      real(dp) :: c_ss = 0.0_dp
         !! A double hybrid's weight on same-spin RI-MP2, zero otherwise.
      logical :: double_hybrid = .false.
      logical :: unrestricted = .false.
         !! Every state runs unrestricted, the closed shells included.
      logical :: needs_rimp2 = .false.
         !! RI-MP2 runs after every SCF: method `ri-mp2`, or a double hybrid.
   contains
      procedure :: close => session_close
      procedure :: message => session_message
   end type terco_session_t

   type :: terco_state_t
      !! One converged (or kept) SCF, and what was run on top of it
      real(dp) :: energy_scf = 0.0_dp
         !! The SCF energy, hartree, the continuum's dielectric term included.
      real(dp) :: e_os = 0.0_dp
         !! Unscaled opposite-spin RI-MP2, zero unless the session needs it.
      real(dp) :: e_ss = 0.0_dp
         !! Unscaled same-spin RI-MP2, zero unless the session needs it.
      real(dp) :: c_os = 0.0_dp   !! The double hybrid's weight on `e_os`, zero otherwise
      real(dp) :: c_ss = 0.0_dp   !! The double hybrid's weight on `e_ss`, zero otherwise
      real(dp) :: e_diel = 0.0_dp      !! The continuum's dielectric energy
      real(dp) :: q_surface = 0.0_dp   !! Total surface charge of the continuum
      integer :: iterations = 0
      logical :: converged = .false.
      integer :: nspin = 1
      real(dp), allocatable :: density(:, :, :)
         !! (nao, nao, nspin). Restricted holds the total, unrestricted alpha
         !! then beta.
      real(dp), allocatable :: mo_energies(:, :)   !! (nao, nspin), hartree
      real(dp), allocatable :: orbitals(:, :, :)   !! (nao, nao, nspin), columns are orbitals
      real(dp), allocatable :: occupations(:, :)   !! (nmo, nspin)
   contains
      procedure :: total => state_total
   end type terco_state_t

contains

   subroutine run_terco_scf(settings, fragment, result, want_gradient)
      !! Run one fragment through terco, on the device
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(inout) :: result
      logical, intent(in), optional :: want_gradient

      type(terco_session_t) :: session
      type(terco_state_t) :: state
      type(error_t) :: error
      logical :: wants_gradient
      character(len=128) :: line

      wants_gradient = .false.
      if (present(want_gradient)) wants_gradient = want_gradient

      call open_session(settings, fragment, session, error)
      if (.not. error%has_error() .and. wants_gradient) then
         call refuse_gradient(settings, session, error)
      end if
      if (.not. error%has_error()) then
         call set_guess(session, settings%guess, fragment%charge, fragment%multiplicity, error)
      end if
      if (.not. error%has_error()) then
         call run_state(settings, session, fragment%charge, fragment%multiplicity, &
                        settings%allow_crap_scf, state, error)
      end if
      if (error%has_error()) then
         call refuse(result, "backend 'terco': "//error%get_message())
         call session%close()
         return
      end if

      result%scf_iterations = state%iterations
      result%energy%scf = state%energy_scf
      if (settings%run_mp2) then
         result%energy%mp2%ss = state%e_ss
         result%energy%mp2%os = state%e_os
         result%energy%mp2%ss_scale = settings%scs_ss
         result%energy%mp2%os_scale = settings%scs_os
         call result%energy%mp2%check_stability()
      else if (session%double_hybrid) then
         ! Scaled and kept beside the SCF, as the CPU backend stores it.
         result%energy%dh_pt2 = state%c_os*state%e_os + state%c_ss*state%e_ss
      end if
      result%has_energy = .true.

      if (.not. state%converged) then
         call logger%warning("backend 'terco': the SCF did not converge in the "// &
                             "iterations allowed; keeping the last iterate because "// &
                             "'allow_crap_scf' is set.")
      end if
      if (settings%pcm%enabled) then
         write (line, "(a,f18.10,a,f9.4)") "    continuum: E_diel = ", state%e_diel, &
            "   total surface charge ", state%q_surface
         call logger%info(trim(line))
      end if

      if (wants_gradient) then
         call run_gradient(session, result, error)
         if (error%has_error()) then
            call refuse(result, "backend 'terco': "//error%get_message())
            call session%close()
            return
         end if
      end if

      call report_frontier(state, result)
      call analyse_density(settings, fragment, session, state, result)
      if (allocated(settings%fukui_population)) then
         call run_fukui(settings, fragment, session, state, result)
      end if

      call session%close()
   end subroutine run_terco_scf

   subroutine open_session(settings, fragment, session, error)
      !! Everything that does not depend on the electron count, set up once
      !!
      !! Refuses, before anything reaches the device, every request terco
      !! would answer with something other than what the deck asked for.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(terco_session_t), intent(inout) :: session
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: functional_name
      integer(c_int) :: rc, nao_c, nbas_c, eri_lmax, df_lmax, df_lmax_aux, major, minor
      integer(c_int), allocatable :: bas_c(:, :)
      integer :: iatom, l_max

      ! A library whose ABI is not the one the declarations were compiled from
      ! would read every argument this module passes by the wrong layout.
      rc = trc_abi_version(major, minor)
      if (major /= TRC_ABI_MAJOR .or. minor < TRC_ABI_MINOR) then
         call error%set(ERROR_VALIDATION, "the linked libterco has ABI "// &
                        int_text(int(major))//"."//int_text(int(minor))//", and this "// &
                        "build was compiled against "//int_text(int(TRC_ABI_MAJOR))//"."// &
                        int_text(int(TRC_ABI_MINOR))//". Rebuild one against the other.")
         return
      end if

      call refuse_unsupported(settings, error)
      if (error%has_error()) return

      functional_name = trim(adjustl(settings%functional))
      if (len(functional_name) > 0) then
         call describe_functional(functional_name, session, error)
         if (error%has_error()) return
      end if

      ! A double hybrid's perturbative term is RI-MP2 on top of the SCF, so it
      ! needs the auxiliary basis a plain RI-MP2 does.
      session%needs_rimp2 = settings%run_mp2 .or. session%double_hybrid
      if (session%double_hybrid .and. .not. settings%aux_basis_named) then
         call error%set(ERROR_VALIDATION, "'"//functional_name//"' is a double hybrid, "// &
                        "and terco evaluates its perturbative term by RI-MP2 only. Name "// &
                        "a correlation-fitting basis in model.aux_basis, or ask for "// &
                        "backend 'libcint' to compute it exactly.")
         return
      end if

      ! ---- the molecule, in libcint's packing ------------------------------
      allocate (session%symbols(fragment%n_atoms))
      do iatom = 1, fragment%n_atoms
         session%symbols(iatom) = element_number_to_symbol(fragment%element_numbers(iatom))
      end do
      call build_czt_molecule(fragment%element_numbers, session%symbols, &
                              fragment%coordinates, settings%basis_set, session%mol, error, &
                              force_cartesian=settings%cartesian, &
                              ecp_name=settings%ecp_set)
      if (error%has_error()) return

      ! terco's basis entry reads orbital shells only; an ECP row handed to it
      ! would be taken for a basis shell with a nonsense angular momentum.
      if (session%mol%necpbas > 0) then
         call error%set(ERROR_VALIDATION, "an effective core potential was asked for, "// &
                        "which terco does not read. Ask for backend 'libcint', or drop "// &
                        "'model.ecp'.")
         return
      end if

      ! Angular momentum against the linked build's own ceiling.
      rc = trc_max_angular_momentum(eri_lmax, df_lmax, df_lmax_aux)
      l_max = max_ao_l(session%mol)
      if (l_max > int(eri_lmax)) then
         call error%set(ERROR_VALIDATION, "this basis carries angular momentum "// &
                        int_text(l_max)//", and the linked terco was built with "// &
                        "four-centre kernels up to "//int_text(int(eri_lmax))//". Rebuild "// &
                        "terco with a higher TERCO_LMAX, ask for backend 'libcint', or "// &
                        "choose a smaller basis.")
         return
      end if
      ! A Cartesian AO-indexed array crosses in libcint's within-shell order.
      ! A terco compiled for GAMESS order would hand back densities whose d
      ! blocks are silently permuted. Spherical shells are libcint's order
      ! whatever that setting.
      if (session%mol%cartesian .and. l_max >= 2 .and. trc_cart_order() /= 0) then
         call error%set(ERROR_VALIDATION, "the linked terco was compiled with a "// &
                        "non-libcint Cartesian order (TRC_CART_ORDER), so its d-shell "// &
                        "densities would not line up with this basis. Rebuild terco "// &
                        "with the default order.")
         return
      end if

      ! ---- one context, set up in steps --------------------------------------
      rc = trc_create(session%handle)
      if (rc /= TRC_OK) then
         call error%set(ERROR_VALIDATION, "could not create a terco context (status "// &
                        status_text(rc)//")")
         return
      end if
      call set_state_molecule(session, fragment%charge, fragment%multiplicity, error)
      if (error%has_error()) return
      call segment_general_contractions(session%mol%bas, session%mol%nbas, bas_c, nbas_c)
      rc = trc_set_basis_libcint(session%handle, int(session%mol%atm, c_int), &
                                 int(session%mol%natm, c_int), bas_c, nbas_c, &
                                 real(session%mol%env, c_double), &
                                 int(size(session%mol%env), c_int), cartesian_flag(session%mol))
      if (rc /= TRC_OK) then
         call set_terco_error(session, rc, "the basis was refused", error)
         return
      end if
      rc = trc_nao(session%handle, nao_c)
      if (rc /= TRC_OK .or. int(nao_c) /= session%mol%nao) then
         call error%set(ERROR_VALIDATION, "terco and the CPU path disagree about the "// &
                        "size of this basis. That is a bug in the adapter, not in the deck.")
         return
      end if
      session%nao = int(nao_c)

      ! An unknown functional comes back as TRC_ERR_UNSUPPORTED with terco's
      ! own list in the message, so there is no second list to keep here.
      rc = trc_set_method(session%handle, c_string(functional_name), &
                          int(settings%grid_level, c_int))
      if (rc /= TRC_OK) then
         call set_terco_error(session, rc, "the method was refused", error)
         return
      end if
      if (settings%nlc_grid_level >= 0) then
         rc = trc_set_nlc_grid(session%handle, int(settings%nlc_grid_level, c_int))
         if (rc /= TRC_OK) then
            call set_terco_error(session, rc, "the VV10 grid level was refused", error)
            return
         end if
      end if

      rc = trc_set_convergence(session%handle, real(settings%energy_tol, c_double), &
                               real(commutator_gate(settings), c_double), &
                               int(settings%max_iter, c_int), int(settings%diis_size, c_int))
      if (rc == TRC_OK) rc = trc_set_screening(session%handle, &
                                               real(settings%screening_tolerance, c_double))
      if (rc == TRC_OK) rc = trc_set_incremental(session%handle, &
                                                 int(merge(1, 0, settings%incremental_fock), c_int))
      if (rc == TRC_OK) rc = trc_set_verbose(session%handle, &
                                             int(merge(1, 0, settings%verbose), c_int))
      session%unrestricted = settings%unrestricted
      if (rc == TRC_OK) rc = trc_set_unrestricted(session%handle, &
                                                  int(merge(1, 0, settings%unrestricted), c_int))
      ! terco applies a level shift only until DIIS takes over. It moves the
      ! path, not the converged state.
      if (rc == TRC_OK) rc = trc_set_level_shift(session%handle, &
                                                 real(settings%level_shift, c_double))
      if (rc /= TRC_OK) then
         call set_terco_error(session, rc, "a setting was refused", error)
         return
      end if

      ! The cavity belongs to the geometry and stays set for every state this
      ! session runs, so a difference of two states is a difference of two
      ! solvated energies.
      if (settings%pcm%enabled) then
         call set_continuum(settings, fragment, session, error)
         if (error%has_error()) return
      end if

      if (settings%density_fitting) then
         call set_fitting_basis(settings, fragment, session, error)
         if (error%has_error()) return
      end if

      if (session%needs_rimp2) call set_correlation_basis(settings, fragment, session, error)
   end subroutine open_session

   subroutine set_fitting_basis(settings, fragment, session, error)
      !! Density-fitted J and K, over the deck's auxiliary basis
      !!
      !! The same `model.aux_basis` the CPU backend fits with, built the same
      !! way, in the orbital basis's Cartesian or spherical form.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(terco_session_t), intent(inout) :: session
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: aux
      integer(c_int), allocatable :: aux_bas_c(:, :)
      integer(c_int) :: naux_bas_c, rc

      ! terco fits J and K for a closed shell only; an open shell would run
      ! exact or not at all, and neither is what the deck asked for.
      if (fragment%multiplicity /= 1 .or. settings%unrestricted) then
         call error%set(ERROR_VALIDATION, "keywords.scf.density_fitting on an "// &
                        "unrestricted reference: terco fits J and K for a restricted "// &
                        "closed shell only. Drop the fitting, or ask for backend 'libcint'.")
         return
      end if
      call build_czt_molecule(fragment%element_numbers, session%symbols, &
                              fragment%coordinates, trim(settings%aux_basis_set), aux, error, &
                              force_cartesian=settings%cartesian)
      if (error%has_error()) return
      call segment_general_contractions(aux%bas, aux%nbas, aux_bas_c, naux_bas_c)
      rc = trc_set_jkfit_libcint(session%handle, int(aux%atm, c_int), int(aux%natm, c_int), &
                                 aux_bas_c, naux_bas_c, real(aux%env, c_double), &
                                 int(size(aux%env), c_int), cartesian_flag(aux))
      call aux%destroy()
      ! A fitted build is a full build every iteration, so terco refuses it
      ! beside incremental building rather than ignoring one of the two.
      if (rc == TRC_OK) rc = trc_set_df(session%handle, 1_c_int)
      if (rc == TRC_OK) rc = trc_set_incremental(session%handle, 0_c_int)
      if (rc /= TRC_OK) call set_terco_error(session, rc, "the fitting basis '"// &
                                             trim(settings%aux_basis_set)//"' was refused", error)
   end subroutine set_fitting_basis

   subroutine set_continuum(settings, fragment, session, error)
      !! The continuum, over the CPU backend's own cavity radii
      !!
      !! terco's own table is modified Bondi, with hydrogen at 1.10 Angstrom;
      !! mqc's is Bondi's, hydrogen at 1.20. Handing over the final radii makes
      !! the two backends build one cavity rather than two.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(terco_session_t), intent(inout) :: session
      type(error_t), intent(inout) :: error

      real(c_double), allocatable, target :: radii(:)
      real(dp) :: radius
      integer(c_int) :: rc
      integer :: iatom

      allocate (radii(fragment%n_atoms))
      do iatom = 1, fragment%n_atoms
         call cavity_radius(fragment%element_numbers(iatom), settings%pcm%radii_scale, &
                            radius, error)
         if (error%has_error()) return
         radii(iatom) = real(radius, c_double)
      end do
      rc = trc_set_pcm(session%handle, c_string(trim(settings%pcm%method)), &
                       real(settings%pcm%dielectric, c_double), &
                       int(settings%pcm%angular_points, c_int), &
                       real(settings%pcm%radii_scale, c_double))
      if (rc == TRC_OK) rc = trc_set_pcm_radii(session%handle, int(fragment%n_atoms, c_int), &
                                               c_loc(radii))
      if (rc /= TRC_OK) call set_terco_error(session, rc, "the continuum was refused", error)
   end subroutine set_continuum

   subroutine describe_functional(name, session, error)
      !! What terco knows about a functional, before anything is built
      !!
      !! A double hybrid's PT2 weights come from terco's table and are checked
      !! against mqc's, which the CPU backend would apply to the same deck.
      character(len=*), intent(in) :: name
      type(terco_session_t), intent(inout) :: session
      type(error_t), intent(inout) :: error

      integer(c_int) :: rc, flags, needed
      real(c_double) :: exx, c_os, c_ss
      character(kind=c_char), allocatable :: list_c(:)
      character(len=:), allocatable :: list
      type(xc_spec_t) :: spec
      type(error_t) :: spec_error
      integer :: i

      rc = trc_functional_info(c_string(name), flags, exx, c_os, c_ss)
      if (rc /= TRC_OK .or. iand(flags, FUNCTIONAL_RUNNABLE) == 0) then
         list = ""
         allocate (list_c(1))
         rc = trc_functional_list(list_c, 0_c_int, needed)
         deallocate (list_c)
         allocate (list_c(max(needed, 1_c_int)))
         if (trc_functional_list(list_c, needed, needed) == TRC_OK) then
            do i = 1, size(list_c)
               if (list_c(i) == c_null_char) exit
               list = list//list_c(i)
            end do
         end if
         call error%set(ERROR_VALIDATION, "terco cannot run the functional '"//name// &
                        "'. It runs: "//list//". Ask for backend 'libcint', which goes "// &
                        "through libxc.")
         return
      end if

      session%double_hybrid = iand(flags, FUNCTIONAL_DOUBLE_HYBRID) /= 0
      if (.not. session%double_hybrid) return
      session%c_os = real(c_os, dp)
      session%c_ss = real(c_ss, dp)
      call xc_spec_from_name(name, spec, spec_error)
      if (spec_error%has_error()) return
      if (abs(session%c_os - spec%pt2_fraction) > DH_WEIGHT_AGREEMENT .or. &
          abs(session%c_ss - spec%pt2_fraction) > DH_WEIGHT_AGREEMENT) then
         call error%set(ERROR_VALIDATION, "terco and mqc weight the perturbative term "// &
                        "of '"//name//"' differently, so the two backends would report "// &
                        "different energies for it. That is a bug in one functional "// &
                        "table, not in the deck.")
      end if
   end subroutine describe_functional

   subroutine refuse_gradient(settings, session, error)
      !! What terco's analytic gradient does not cover
      !!
      !! Refused rather than returned without the missing term: a gradient
      !! short of its correlation or solvent part is converged, plausible and
      !! wrong, and an optimiser steps on it.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(terco_session_t), intent(in) :: session
      type(error_t), intent(inout) :: error

      if (session%needs_rimp2) then
         call error%set(ERROR_VALIDATION, "a gradient was asked for with RI-MP2 or a "// &
                        "double hybrid, and terco has no correlation gradient. Ask for "// &
                        "backend 'libcint', or request an energy.")
      else if (settings%pcm%enabled) then
         call error%set(ERROR_VALIDATION, "a gradient was asked for under keywords.pcm, "// &
                        "and terco has no continuum gradient. Ask for backend 'cuest', "// &
                        "or request an energy.")
      end if
   end subroutine refuse_gradient

   subroutine run_gradient(session, result, error)
      !! terco's analytic nuclear gradient of the state just converged
      !!
      !! `(3, n_atoms)` in hartree per bohr, in the fragment's atom order, caps
      !! included, as the CPU backend returns it. One rank: terco refuses more.
      type(terco_session_t), intent(inout) :: session
      type(calculation_result_t), intent(inout) :: result
      type(error_t), intent(inout) :: error

      real(c_double), allocatable :: gradient(:, :)
      integer(c_int) :: rc

      rc = trc_run_gradient(session%handle)
      if (rc /= TRC_OK) then
         call set_terco_error(session, rc, "the gradient failed", error)
         return
      end if
      allocate (gradient(3, session%mol%natm))
      rc = trc_gradient(session%handle, gradient)
      if (rc /= TRC_OK) then
         call set_terco_error(session, rc, "the gradient could not be read", error)
         return
      end if
      result%gradient = real(gradient, dp)
      result%has_gradient = .true.
   end subroutine run_gradient

   subroutine refuse_unsupported(settings, error)
      !! What the deck asked for that terco has no counterpart for
      !!
      !! Each of these would otherwise be dropped on the way to the device, and
      !! the run would report an answer to a different question.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(error_t), intent(inout) :: error

      if (settings%run_cc) then
         call error%set(ERROR_VALIDATION, "coupled cluster has no terco implementation; "// &
                        "it runs through the CPU backend. Ask for 'auto', or drop it.")
      else if (settings%run_mp2 .and. .not. settings%corr_density_fitting) then
         call error%set(ERROR_VALIDATION, "terco carries RI-MP2 only: ask for method "// &
                        "'ri-mp2' with model.aux_basis, or for backend 'auto' to run "// &
                        "exact MP2 on the CPU.")
      else if (settings%excited%enabled) then
         call error%set(ERROR_VALIDATION, "keywords.excited_states has no terco "// &
                        "implementation; the linear-response solver lives on the CPU "// &
                        "backend. Ask for backend 'libcint', or drop the excited states.")
      else if (settings%stability) then
         call error%set(ERROR_VALIDATION, "keywords.scf.stability has no terco "// &
                        "implementation. Ask for backend 'libcint', or drop it.")
      else if (settings%second_order) then
         call error%set(ERROR_VALIDATION, "the second-order SCF has no terco "// &
                        "implementation; terco converges by DIIS. Ask for backend "// &
                        "'libcint', or drop it.")
      else if (.not. settings%use_diis .or. trim(settings%accelerator) /= "diis") then
         call error%set(ERROR_VALIDATION, "terco converges by DIIS and nothing else, so "// &
                        "keywords.scf.accelerator '"//trim(settings%accelerator)// &
                        "' has no terco counterpart. Ask for backend 'libcint'.")
      else if (settings%linear_dependence > 0.0_dp) then
         call error%set(ERROR_VALIDATION, "keywords.scf.linear_dependence_threshold "// &
                        "has no terco counterpart. Drop it, or ask for backend 'libcint'.")
      else if (settings%grid_level < 0 .and. len_trim(settings%functional) > 0) then
         call error%set(ERROR_VALIDATION, "keywords.dft.radial_points/angular_points "// &
                        "have no terco counterpart; terco takes a grid level. Give "// &
                        "keywords.dft.grid_level instead.")
      else if (allocated(settings%guess_steps)) then
         if (size(settings%guess_steps) > 0) then
            call error%set(ERROR_VALIDATION, "a basis-set projection guess runs its "// &
                           "smaller-basis SCFs on the CPU and has no terco counterpart. "// &
                           "Use 'sad', 'sap', 'minao', 'sac', 'core' or 'gwh'.")
         end if
      end if
   end subroutine refuse_unsupported

   subroutine set_correlation_basis(settings, fragment, session, error)
      !! The RI-MP2 fitting basis and frozen core, for RI-MP2 or a double hybrid
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(terco_session_t), intent(inout) :: session
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: aux
      integer(c_int), allocatable :: aux_bas_c(:, :)
      integer(c_int) :: naux_bas_c, rc
      character(len=:), allocatable :: aux_name
      integer :: frozen

      aux_name = trim(settings%aux_basis_set)
      if (len_trim(aux_name) == 0) then
         call error%set(ERROR_VALIDATION, "RI-MP2 needs an auxiliary basis: set "// &
                        "model.aux_basis")
         return
      end if
      if (index(aux_name, "rifit") == 0 .and. index(aux_name, "-ri") == 0) then
         call logger%warning("bad basis set, calculations will be poor: '"//aux_name// &
                             "' is not a correlation-fitting (RIFIT) set")
      end if
      call build_czt_molecule(fragment%element_numbers, session%symbols, &
                              fragment%coordinates, aux_name, aux, error, &
                              force_cartesian=settings%cartesian)
      if (error%has_error()) return
      call segment_general_contractions(aux%bas, aux%nbas, aux_bas_c, naux_bas_c)
      rc = trc_set_aux_libcint(session%handle, int(aux%atm, c_int), int(aux%natm, c_int), &
                               aux_bas_c, naux_bas_c, real(aux%env, c_double), &
                               int(size(aux%env), c_int), cartesian_flag(aux))
      call aux%destroy()
      if (rc /= TRC_OK) then
         call set_terco_error(session, rc, "the auxiliary basis was refused", error)
         return
      end if

      frozen = settings%n_frozen_core
      if (frozen < 0) frozen = core_orbital_count(fragment%real_element_numbers())
      if (.not. settings%freeze_core) frozen = 0
      rc = trc_set_rimp2(session%handle, int(frozen, c_int), 256_c_int)
      if (rc /= TRC_OK) then
         call set_terco_error(session, rc, "the RI-MP2 settings were refused", error)
      end if
   end subroutine set_correlation_basis

   subroutine set_guess(session, guess_name, charge, multiplicity, error)
      !! A named guess, for the next state the session runs
      !!
      !! terco builds core, GWH, SAD and SAP itself, on the device. SAC and
      !! MINAO it does not have, so mqc builds them here and hands the density
      !! over.
      type(terco_session_t), intent(inout) :: session
      character(len=*), intent(in) :: guess_name
      integer, intent(in) :: charge, multiplicity
      type(error_t), intent(inout) :: error

      integer :: guess_kind, n_alpha, n_beta
      integer(c_int) :: rc
      real(dp), allocatable :: d_alpha(:, :), d_beta(:, :)
      logical :: occupations_ok

      call parse_guess_name(guess_name, guess_kind, error)
      if (error%has_error()) return

      select case (guess_kind)
      case (SCF_GUESS_CORE)
         rc = trc_set_guess(session%handle, TRC_GUESS_CORE, c_null_ptr, 1_c_int)
      case (SCF_GUESS_GWH)
         rc = trc_set_guess(session%handle, TRC_GUESS_GWH, c_null_ptr, 1_c_int)
      case (SCF_GUESS_SAD)
         rc = trc_set_guess(session%handle, TRC_GUESS_SAD, c_null_ptr, 1_c_int)
      case (SCF_GUESS_SAP)
         rc = trc_set_guess(session%handle, TRC_GUESS_SAP, c_null_ptr, 1_c_int)
      case (SCF_GUESS_SAC, SCF_GUESS_MINAO)
         call build_atomic_guess(session%mol, guess_kind, d_alpha, d_beta, error)
         if (error%has_error()) return
         call spin_occupations(nint(sum(session%mol%charges)) - charge, multiplicity, &
                               n_alpha, n_beta, occupations_ok)
         if (.not. occupations_ok) then
            call error%set(ERROR_VALIDATION, "a multiplicity of "//int_text(multiplicity)// &
                           " cannot be built at charge "//int_text(charge)//".")
            return
         end if
         call give_density(session, d_alpha, d_beta, n_alpha /= n_beta .or. session%unrestricted, &
                           error)
         return
      case default
         call error%set(ERROR_VALIDATION, "guess '"//trim(guess_name)//"' has no "// &
                        "terco counterpart; use 'sad', 'sap', 'minao', 'sac', 'core' "// &
                        "or 'gwh'.")
         return
      end select
      if (rc /= TRC_OK) call set_terco_error(session, rc, "the guess was refused", error)
   end subroutine set_guess

   subroutine give_density(session, d_alpha, d_beta, open_shell, error)
      !! Hand terco a starting density built on the CPU
      type(terco_session_t), intent(inout) :: session
      real(dp), intent(in) :: d_alpha(:, :), d_beta(:, :)
      logical, intent(in) :: open_shell
      type(error_t), intent(inout) :: error

      real(c_double), allocatable, target :: given(:, :, :)
      integer(c_int) :: rc

      ! Restricted takes the total, unrestricted alpha then beta. terco copies
      ! it at this call, so the buffer is free on return.
      if (open_shell) then
         allocate (given(session%nao, session%nao, 2))
         given(:, :, 1) = d_alpha
         given(:, :, 2) = d_beta
         rc = trc_set_guess(session%handle, TRC_GUESS_GIVEN, c_loc(given), 2_c_int)
      else
         allocate (given(session%nao, session%nao, 1))
         given(:, :, 1) = d_alpha + d_beta
         rc = trc_set_guess(session%handle, TRC_GUESS_GIVEN, c_loc(given), 1_c_int)
      end if
      if (rc /= TRC_OK) call set_terco_error(session, rc, "the guess density was refused", error)
   end subroutine give_density

   subroutine set_state_molecule(session, charge, multiplicity, error)
      !! The atoms, once, before any basis is set
      !!
      !! Nuclear charges are the molecule's own, so a ghost atom reaches terco
      !! as a centre with basis functions and no nucleus. A later electron
      !! count goes through `set_state_charge`, which keeps what this resets.
      type(terco_session_t), intent(inout) :: session
      integer, intent(in) :: charge, multiplicity
      type(error_t), intent(inout) :: error

      integer(c_int) :: rc

      rc = trc_set_molecule(session%handle, int(session%mol%natm, c_int), &
                            real(session%mol%charges, c_double), &
                            real(reshape(session%mol%coords, [3*session%mol%natm]), c_double), &
                            int(charge, c_int), int(multiplicity, c_int))
      if (rc /= TRC_OK) call set_terco_error(session, rc, "the molecule was refused", error)
   end subroutine set_state_molecule

   subroutine set_state_charge(session, charge, multiplicity, error)
      !! Point an open session at another electron count
      !!
      !! Keeps the molecule, the bases, the fit and the cavity radii. terco
      !! drops a guess density given for the old count, so the next state's
      !! guess is set after this, never before.
      type(terco_session_t), intent(inout) :: session
      integer, intent(in) :: charge, multiplicity
      type(error_t), intent(inout) :: error

      integer(c_int) :: rc

      rc = trc_set_charge(session%handle, int(charge, c_int), int(multiplicity, c_int))
      if (rc /= TRC_OK) call set_terco_error(session, rc, "the charge was refused", error)
   end subroutine set_state_charge

   subroutine run_state(settings, session, charge, multiplicity, keep_unconverged, state, error)
      !! One electron count on an open session: SCF, then RI-MP2 if it is owed
      !!
      !! The session must already be at `charge` and `multiplicity`, and the
      !! guess is whatever it was last given. A state that did not converge is
      !! an error unless `keep_unconverged`.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(terco_session_t), intent(inout) :: session
      integer, intent(in) :: charge, multiplicity
      logical, intent(in) :: keep_unconverged
      type(terco_state_t), intent(out) :: state
      type(error_t), intent(inout) :: error

      integer(c_int) :: rc, niter_c, flag_c, npoints_c, nmo_c
      real(c_double) :: energy_c, e_os_c, e_ss_c, total_c, e_diel_c, q_c
      integer :: n_alpha, n_beta
      logical :: occupations_ok

      call spin_occupations(nint(sum(session%mol%charges)) - charge, multiplicity, &
                            n_alpha, n_beta, occupations_ok)
      if (.not. occupations_ok) then
         call error%set(ERROR_VALIDATION, "a multiplicity of "//int_text(multiplicity)// &
                        " cannot be built at charge "//int_text(charge)//".")
         return
      end if
      state%nspin = merge(2, 1, n_alpha /= n_beta .or. settings%unrestricted)

      rc = trc_run_scf(session%handle)
      if (rc /= TRC_OK .and. rc /= TRC_ERR_NOCONV) then
         call set_terco_error(session, rc, "the SCF failed", error)
         return
      end if
      niter_c = 0
      if (trc_iterations(session%handle, niter_c) == TRC_OK) state%iterations = int(niter_c)
      flag_c = 0
      if (trc_converged(session%handle, flag_c) == TRC_OK) state%converged = flag_c /= 0
      if (rc == TRC_ERR_NOCONV) state%converged = .false.
      if (.not. state%converged .and. .not. keep_unconverged) then
         call error%set(ERROR_VALIDATION, "the SCF did not converge in the iterations "// &
                        "allowed. Raise 'keywords.scf.max_iter', or set "// &
                        "'allow_crap_scf' to keep the last iterate.")
         return
      end if
      if (trc_energy(session%handle, energy_c) /= TRC_OK) then
         call error%set(ERROR_VALIDATION, "terco ran the SCF and would not report its energy")
         return
      end if
      state%energy_scf = real(energy_c, dp)

      if (settings%pcm%enabled) then
         if (trc_pcm_energy(session%handle, e_diel_c, q_c, npoints_c) == TRC_OK) then
            state%e_diel = real(e_diel_c, dp)
            state%q_surface = real(q_c, dp)
         end if
      end if

      allocate (state%density(session%nao, session%nao, state%nspin))
      allocate (state%mo_energies(session%nao, state%nspin))
      allocate (state%orbitals(session%nao, session%nao, state%nspin))
      rc = trc_density(session%handle, state%density)
      if (rc == TRC_OK) rc = trc_mo_energies(session%handle, state%mo_energies)
      if (rc == TRC_OK) rc = trc_mo_coefficients(session%handle, state%orbitals)
      if (rc == TRC_OK) rc = trc_nmo(session%handle, nmo_c)
      if (rc == TRC_OK) then
         allocate (state%occupations(int(nmo_c), state%nspin))
         rc = trc_mo_occupations(session%handle, state%occupations)
      end if
      if (rc /= TRC_OK) then
         call set_terco_error(session, rc, "the converged density could not be read", error)
         return
      end if

      if (.not. session%needs_rimp2) return

      rc = trc_run_rimp2(session%handle)
      if (rc == TRC_OK) rc = trc_rimp2_energy(session%handle, e_os_c, e_ss_c)
      if (rc /= TRC_OK) then
         call set_terco_error(session, rc, "RI-MP2 failed", error)
         return
      end if
      state%e_os = real(e_os_c, dp)
      state%e_ss = real(e_ss_c, dp)
      state%c_os = session%c_os
      state%c_ss = session%c_ss

      ! terco weights a double hybrid's total itself; it has to be the number
      ! formed here from the parts it reported.
      if (session%double_hybrid) then
         if (trc_energy(session%handle, total_c) /= TRC_OK) then
            call error%set(ERROR_VALIDATION, "terco ran RI-MP2 and would not report "// &
                           "the double hybrid's total")
            return
         end if
         if (abs(real(total_c, dp) - state%total()) > 1.0e-10_dp*max(1.0_dp, abs(state%total()))) then
            call error%set(ERROR_VALIDATION, "terco's double-hybrid total is not its "// &
                           "SCF energy plus its own weighted RI-MP2 components. That is "// &
                           "a bug in the adapter, not in the deck.")
            return
         end if
      end if
   end subroutine run_state

   subroutine report_frontier(state, result)
      !! HOMO and LUMO from the orbital energies and occupations terco returned
      type(terco_state_t), intent(in) :: state
      type(calculation_result_t), intent(inout) :: result

      integer :: nmo

      if (.not. allocated(state%occupations)) return
      nmo = size(state%occupations, 1)
      call frontier_orbitals(state%mo_energies(1:nmo, 1), state%occupations(:, 1), &
                             result%homo, result%lumo, result%has_orbitals)
   end subroutine report_frontier

   subroutine analyse_density(settings, fragment, session, state, result)
      !! Charges and bond orders of the density terco converged
      !!
      !! A failure here is a warning, as it is on every backend: the energy
      !! stands.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(terco_session_t), intent(in) :: session
      type(terco_state_t), intent(in) :: state
      type(calculation_result_t), intent(inout) :: result

      real(dp), allocatable :: total(:, :), spin(:, :), s_ao(:, :)
      real(dp), allocatable :: q(:), q_spin(:), orders(:, :), valences(:)
      type(error_t) :: analysis_error
      logical :: open_shell

      if (.not. allocated(settings%charges_scheme) .and. &
          .not. allocated(settings%bond_order_scheme)) return

      open_shell = state%nspin == 2
      if (open_shell) then
         total = state%density(:, :, 1) + state%density(:, :, 2)
         spin = state%density(:, :, 1) - state%density(:, :, 2)
      else
         total = state%density(:, :, 1)
      end if
      call session%mol%overlap(s_ao)

      if (allocated(settings%charges_scheme)) then
         call analysis_error%clear()
         select case (trim(settings%charges_scheme))
         case ("mulliken")
            call mulliken_charges(session%mol, total, s_ao, q, analysis_error)
            if (open_shell .and. .not. analysis_error%has_error()) then
               call mulliken_spin_populations(session%mol, spin, s_ao, q_spin, analysis_error)
            end if
         case ("chelpg")
            call chelpg_charges(session%mol, total, q, analysis_error, &
                                total_charge=real(fragment%charge, dp))
         case default
            call analysis_error%set(ERROR_VALIDATION, "charge scheme '"// &
                                    trim(settings%charges_scheme)//"' is not driven "// &
                                    "through terco; 'mulliken' and 'chelpg' are.")
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
      end if

      if (allocated(settings%bond_order_scheme)) then
         call analysis_error%clear()
         select case (trim(settings%bond_order_scheme))
         case ("mayer")
            if (open_shell) then
               call mayer_bond_orders_open_shell(session%mol, state%density(:, :, 1), &
                                                 state%density(:, :, 2), s_ao, orders, &
                                                 analysis_error)
            else
               call mayer_bond_orders(session%mol, total, s_ao, orders, analysis_error)
            end if
         case default
            call analysis_error%set(ERROR_VALIDATION, "unknown bond order scheme '"// &
                                    trim(settings%bond_order_scheme)//"'; expected 'mayer'.")
         end select
         if (analysis_error%has_error()) then
            call logger%warning("  bond orders could not be computed: "// &
                                analysis_error%get_message())
         else
            call mayer_valences(orders, valences)
            call print_mayer_report(orders, valences, session%symbols)
            call move_alloc(orders, result%bond_orders)
            call move_alloc(valences, result%bond_order_valences)
            result%bond_order_scheme = trim(settings%bond_order_scheme)
            result%has_bond_orders = .true.
         end if
      end if
   end subroutine analyse_density

   subroutine run_fukui(settings, fragment, session, neutral, result)
      !! Fukui indices from the neutral just converged and two ions run here
      !!
      !! Both ions run on the session the neutral ran on, so the three share
      !! the basis, the grid, the frozen core, the RI-MP2 basis and the cavity,
      !! and IP and EA are differences of like with like. A failure here is a
      !! warning, as on the CPU backend: the energy stands.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      type(terco_session_t), intent(inout) :: session
      type(terco_state_t), intent(in) :: neutral
      type(calculation_result_t), intent(inout) :: result

      type(terco_state_t) :: anion, cation
      type(fukui_result_t) :: fukui
      type(error_t) :: fukui_error
      character(len=:), allocatable :: mode
      integer :: nelec

      nelec = nint(sum(session%mol%charges)) - fragment%charge
      if (fragment%multiplicity /= 1 .or. nelec < 2) then
         call fukui_error%set(ERROR_VALIDATION, "Fukui indices are computed for a "// &
                              "closed-shell molecule, whose ions are both doublets.")
      end if
      if (.not. fukui_error%has_error()) call set_ion_numerics(settings, session, fukui_error)

      mode = "neutral"
      if (allocated(settings%fukui_guess)) mode = trim(settings%fukui_guess)

      if (settings%verbose) call logger%info("  Fukui: anion SCF, N+1 electrons")
      if (.not. fukui_error%has_error()) then
         call set_state_charge(session, fragment%charge - 1, 2, fukui_error)
      end if
      if (.not. fukui_error%has_error()) then
         call seed_ion(settings, session, neutral, mode, fragment%charge - 1, nelec/2 + 1, &
                       nelec/2, fukui_error)
      end if
      if (.not. fukui_error%has_error()) then
         call run_state(settings, session, fragment%charge - 1, 2, &
                        settings%fukui_scf%allow_crap_scf, anion, fukui_error)
         if (fukui_error%has_error()) call fukui_error%add_context("the anion")
      end if

      if (settings%verbose) call logger%info("  Fukui: cation SCF, N-1 electrons")
      if (.not. fukui_error%has_error()) then
         call set_state_charge(session, fragment%charge + 1, 2, fukui_error)
      end if
      if (.not. fukui_error%has_error()) then
         call seed_ion(settings, session, neutral, mode, fragment%charge + 1, nelec/2, &
                       nelec/2 - 1, fukui_error)
      end if
      if (.not. fukui_error%has_error()) then
         call run_state(settings, session, fragment%charge + 1, 2, &
                        settings%fukui_scf%allow_crap_scf, cation, fukui_error)
         if (fukui_error%has_error()) call fukui_error%add_context("the cation")
      end if

      if (.not. fukui_error%has_error()) then
         if (.not. (anion%converged .and. cation%converged)) then
            call logger%warning("  Fukui: an ion did not converge, and allow_crap_scf "// &
                                "kept it. Every index below is built from it and is a "// &
                                "result to follow up.")
         end if
         call fukui_from_states(session%mol, settings%fukui_population, fragment%charge, &
                                neutral%density(:, :, 1), &
                                anion%density(:, :, 1) + anion%density(:, :, 2), &
                                cation%density(:, :, 1) + cation%density(:, :, 2), &
                                neutral%total(), anion%total(), cation%total(), &
                                fukui, fukui_error)
      end if
      if (fukui_error%has_error()) then
         call logger%warning("  the Fukui analysis could not run: "// &
                             fukui_error%get_message())
         return
      end if

      fukui%functional = trim(settings%functional)
      fukui%iterations_anion = anion%iterations
      fukui%iterations_cation = cation%iterations
      call print_fukui_report(fukui, session%symbols)
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
   end subroutine run_fukui

   subroutine set_ion_numerics(settings, session, error)
      !! `properties.fukui.scf`, for the two ion SCFs
      type(cuest_scf_settings_t), intent(in) :: settings
      type(terco_session_t), intent(inout) :: session
      type(error_t), intent(inout) :: error

      integer(c_int) :: rc
      real(dp) :: gate

      ! The same refusals the neutral gets, for the ions' own copy of the
      ! numerics: each would otherwise be dropped on the way to the device.
      if (.not. settings%fukui_scf%use_diis .or. &
          trim(settings%fukui_scf%accelerator) /= "diis" .or. &
          settings%fukui_scf%second_order) then
         call error%set(ERROR_VALIDATION, "terco converges the ions by DIIS and "// &
                        "nothing else; properties.fukui.scf asks for something else. "// &
                        "Run on backend 'libcint'.")
         return
      end if
      if (settings%fukui_scf%linear_dependence > 0.0_dp) then
         call error%set(ERROR_VALIDATION, "properties.fukui.scf."// &
                        "linear_dependence_threshold has no terco counterpart. Drop "// &
                        "it, or run on backend 'libcint'.")
         return
      end if

      gate = settings%fukui_scf%grad_tol
      if (gate <= 0.0_dp) gate = sqrt(settings%fukui_scf%energy_tol)
      rc = trc_set_convergence(session%handle, real(settings%fukui_scf%energy_tol, c_double), &
                               real(gate, c_double), int(settings%fukui_scf%max_iter, c_int), &
                               int(settings%fukui_scf%diis_size, c_int))
      if (rc == TRC_OK) rc = trc_set_level_shift(session%handle, &
                                                 real(settings%fukui_scf%level_shift, c_double))
      if (rc == TRC_OK .and. .not. settings%density_fitting) then
         rc = trc_set_incremental(session%handle, &
                                  int(merge(1, 0, settings%fukui_scf%incremental_fock), c_int))
      end if
      if (rc /= TRC_OK) call set_terco_error(session, rc, "the ions' SCF settings were refused", error)
   end subroutine set_ion_numerics

   subroutine seed_ion(settings, session, neutral, mode, charge, n_alpha, n_beta, error)
      !! The starting density of one ion
      !!
      !! "neutral" fills the neutral's own orbitals at the ion's occupation,
      !! which has the ion's electron count and spin where a halved neutral
      !! density has neither. "independent" takes `properties.fukui.scf.guess`,
      !! GWH when that is `auto`, as on the CPU backend.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(terco_session_t), intent(inout) :: session
      type(terco_state_t), intent(in) :: neutral
      character(len=*), intent(in) :: mode
      integer, intent(in) :: charge, n_alpha, n_beta
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: d_alpha(:, :), d_beta(:, :)

      if (mode == "independent") then
         if (trim(settings%fukui_scf%guess) == "auto") then
            call set_guess(session, "gwh", charge, 2, error)
         else
            call set_guess(session, settings%fukui_scf%guess, charge, 2, error)
         end if
         return
      end if
      allocate (d_alpha(session%nao, session%nao), d_beta(session%nao, session%nao))
      call build_density_spin(neutral%orbitals(:, :, 1), n_alpha, d_alpha)
      call build_density_spin(neutral%orbitals(:, :, 1), n_beta, d_beta)
      call give_density(session, d_alpha, d_beta, .true., error)
   end subroutine seed_ion

   pure function state_total(this) result(energy)
      !! The state's total energy: the SCF, and a double hybrid's weighted PT2
      class(terco_state_t), intent(in) :: this
      real(dp) :: energy

      energy = this%energy_scf + this%c_os*this%e_os + this%c_ss*this%e_ss
   end function state_total

   subroutine session_close(this)
      !! Release the context and the molecule; safe to call twice
      class(terco_session_t), intent(inout) :: this
      integer(c_int) :: rc

      if (c_associated(this%handle)) then
         rc = trc_destroy(this%handle)
         this%handle = c_null_ptr
      end if
      call this%mol%destroy()
   end subroutine session_close

   function session_message(this) result(text)
      !! terco's last message for this context, or an empty string
      class(terco_session_t), intent(in) :: this
      character(len=:), allocatable :: text

      character(kind=c_char) :: msg_c(MSG_LEN)
      integer :: i

      text = ""
      if (.not. c_associated(this%handle)) return
      if (trc_message(this%handle, msg_c, int(MSG_LEN, c_int)) /= TRC_OK) return
      do i = 1, MSG_LEN
         if (msg_c(i) == c_null_char) exit
         text = text//msg_c(i)
      end do
   end function session_message

   subroutine set_terco_error(session, rc, what, error)
      !! An error carrying terco's status and its own explanation
      type(terco_session_t), intent(in) :: session
      integer(c_int), intent(in) :: rc
      character(len=*), intent(in) :: what
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: detail

      detail = session%message()
      if (len(detail) > 0) then
         call error%set(ERROR_VALIDATION, what//" (status "//status_text(rc)//"): "//detail)
      else
         call error%set(ERROR_VALIDATION, what//" (status "//status_text(rc)//")")
      end if
   end subroutine set_terco_error

   subroutine refuse(result, why)
      !! Mark a result failed with a message
      type(calculation_result_t), intent(inout) :: result
      character(len=*), intent(in) :: why

      call result%error%set(ERROR_VALIDATION, why)
      result%has_error = .true.
      result%has_energy = .false.
   end subroutine refuse

   pure function commutator_gate(settings) result(gate)
      !! The commutator-norm threshold terco stops on, beside the energy change
      type(cuest_scf_settings_t), intent(in) :: settings
      real(dp) :: gate

      ! The norm grows with the system, so mqc's per-element density tolerance
      ! is the wrong scale for it. An explicit gradient tolerance is that norm;
      ! otherwise the square root of the energy tolerance, where PySCF puts its
      ! gradient threshold. cholesterol at 1e-6 took 41 iterations, the last
      ! dozen chasing a norm that hovered just above the gate.
      if (settings%grad_tol > 0.0_dp) then
         gate = settings%grad_tol
      else
         gate = sqrt(settings%energy_tol)
      end if
   end function commutator_gate

   pure function cartesian_flag(mol) result(flag)
      !! terco's `cartesian` argument for a basis built here
      type(czt_molecule_t), intent(in) :: mol
      integer(c_int) :: flag

      flag = int(merge(1, 0, mol%cartesian), c_int)
   end function cartesian_flag

   pure subroutine segment_general_contractions(bas, nbas, bas_out, nbas_out)
      !! One `bas` row per contracted function, which is all terco can read
      !!
      !! A generally contracted shell stores an `(nprim, nctr)` coefficient
      !! matrix; terco reads `NPRIM_OF` coefficients from `PTR_COEFF` and has
      !! no notion of `NCTR_OF`. Each column becomes its own shell over the
      !! same exponents, with `PTR_COEFF` advanced to that column.
      !!
      !! AO order is preserved: libcint numbers a shell's functions with the
      !! contraction index slowest, so consecutive split shells occupy exactly
      !! the positions the one general shell did.
      integer, intent(in) :: bas(:, :)
      integer, intent(in) :: nbas
      integer(c_int), allocatable, intent(out) :: bas_out(:, :)
      integer(c_int), intent(out) :: nbas_out

      integer :: ish, ictr, nprim, nctr, row

      ! Splitting repeats the primitive work the columns used to share and
      ! grows the shell-pair count with it. Correct and slower beats refused.
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
      integer, intent(in) :: n_electrons, multiplicity
      integer, intent(out) :: n_alpha, n_beta
      logical, intent(out) :: ok

      integer :: unpaired

      ! Written here rather than taken from the cuEST backend: this backend
      ! must build in a tree where cuEST is switched off.
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

   pure function int_text(value) result(text)
      !! An integer without padding
      integer, intent(in) :: value
      character(len=:), allocatable :: text

      character(len=16) :: buffer

      write (buffer, "(i0)") value
      text = trim(buffer)
   end function int_text

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
