!! The terco backend's dispatch and the guarantees that hold either way
module test_mqc_terco_backend
   !! **Every check here is true whether or not terco is linked**, which is what
   !! makes the file worth compiling in CI: `MQC_ENABLE_TERCO` is off by default
   !! and no CI job builds the library, so a test that only meant something in a
   !! terco build would be 200 lines that never run.
   !!
   !! Two things survive that. The backend *name* is parsed by code that is
   !! always compiled, so the deck-facing spelling is checkable anywhere. And
   !! the refusals are stated as invariants over the result -- "asking terco for
   !! a correlated gradient never yields an energy" holds in a stub build
   !! because the stub refuses the whole calculation, and in a real build
   !! because the driver refuses that gradient. The reason differs; the
   !! guarantee does not.
   !!
   !! The comparisons against the CPU path run only where terco is linked.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_cuest_iface, only: parse_backend_name, cuest_scf_settings_t, &
                              BACKEND_AUTO, BACKEND_CUEST, BACKEND_CZT, BACKEND_TERCO
   use mqc_terco_bridge, only: run_terco_scf, terco_backend_available
   use mqc_czt_bridge, only: run_czt_hf
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   implicit none
   private

   public :: collect_mqc_terco_backend_tests

contains

   subroutine collect_mqc_terco_backend_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("deck_spelling", test_deck_spelling), &
                  new_unittest("terco_is_its_own_backend", test_distinct_kind), &
                  new_unittest("unknown_name_refused", test_unknown_name), &
                  new_unittest("correlated_gradient_never_yields_energy", &
                               test_correlated_gradient_refused), &
                  new_unittest("correlated_never_yields_energy", test_correlated_refused), &
                  new_unittest("ghosted_rhf_agrees_with_the_cpu_path", test_ghosted_rhf_matches_cpu), &
                  new_unittest("ghosted_b3lyp_agrees_with_the_cpu_path", &
                               test_ghosted_b3lyp_matches_cpu), &
                  new_unittest("ghosted_continuum_never_yields_energy", &
                               test_ghosted_pcm_refused), &
                  new_unittest("missing_build_is_reported", test_missing_build), &
                  new_unittest("agrees_with_the_cpu_path", test_matches_cpu), &
                  new_unittest("agrees_on_general_contraction", test_general_contraction), &
                  new_unittest("agrees_on_spherical_shells", test_spherical_matches_cpu), &
                  new_unittest("gradient_agrees_with_the_cpu_path", test_gradient_matches_cpu), &
                  new_unittest("fitted_energy_agrees_with_the_cpu_path", test_fitted_matches_cpu) &
                  ]
   end subroutine collect_mqc_terco_backend_tests

   subroutine test_deck_spelling(error)
      !! `backend: "terco"`, in whatever case the deck wrote it
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer :: kind

      call parse_backend_name("terco", kind, err)
      call check(error,.not. err%has_error(), "'terco' was refused as a backend name")
      if (allocated(error)) return
      call check(error, kind == BACKEND_TERCO, "'terco' did not parse to BACKEND_TERCO")
      if (allocated(error)) return

      call parse_backend_name("TERCO", kind, err)
      call check(error,.not. err%has_error(), "'TERCO' was refused; the parse is "// &
                 "meant to be case-insensitive")
      if (allocated(error)) return
      call check(error, kind == BACKEND_TERCO, "'TERCO' did not parse to BACKEND_TERCO")
   end subroutine test_deck_spelling

   subroutine test_distinct_kind(error)
      !! terco is not a spelling of one of the backends that already existed
      !!
      !! Worth stating: the kinds are plain integer parameters, and a duplicated
      !! value would send every terco deck to whichever branch the `select case`
      !! reached first, with no diagnostic anywhere.
      type(error_type), allocatable, intent(out) :: error

      call check(error, BACKEND_TERCO /= BACKEND_AUTO, "BACKEND_TERCO collides with BACKEND_AUTO")
      if (allocated(error)) return
      call check(error, BACKEND_TERCO /= BACKEND_CUEST, "BACKEND_TERCO collides with BACKEND_CUEST")
      if (allocated(error)) return
      call check(error, BACKEND_TERCO /= BACKEND_CZT, "BACKEND_TERCO collides with BACKEND_CZT")
   end subroutine test_distinct_kind

   subroutine test_unknown_name(error)
      !! Adding a name must not have opened the gate to any name
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer :: kind

      call parse_backend_name("tercoo", kind, err)
      call check(error, err%has_error(), "a misspelled backend name was accepted")
   end subroutine test_unknown_name

   subroutine test_correlated_gradient_refused(error)
      !! terco has no RI-MP2 gradient, and must never answer that request
      !!
      !! The failure this guards against is not an exception -- it is a geometry
      !! optimisation that receives the SCF gradient alone under an MP2 energy,
      !! converged, plausible and missing the correlation's pull on every atom.
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call water(fragment)
      settings%basis_set = "6-31g"
      settings%run_mp2 = .true.
      settings%corr_density_fitting = .true.
      settings%aux_basis_set = "cc-pvdz-rifit"
      settings%aux_basis_named = .true.

      call run_terco_scf(settings, fragment, result, want_gradient=.true.)

      call check(error, result%has_error, "asking terco for an RI-MP2 gradient did not "// &
                 "set an error")
      if (allocated(error)) return
      call check(error,.not. result%has_energy, "a refused gradient still reported an energy")
      if (allocated(error)) return
      call check(error,.not. allocated(result%gradient), &
                 "a refused gradient request still allocated a gradient")
   end subroutine test_correlated_gradient_refused

   subroutine test_correlated_refused(error)
      !! MP2 and coupled cluster have no terco implementation
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call water(fragment)
      settings%basis_set = "sto-3g"
      settings%run_mp2 = .true.

      call run_terco_scf(settings, fragment, result)

      call check(error, result%has_error, "asking terco for MP2 did not set an error")
      if (allocated(error)) return
      call check(error,.not. result%has_energy, "a refused MP2 still reported an energy")
   end subroutine test_correlated_refused

   subroutine test_ghosted_rhf_matches_cpu(error)
      !! A monomer in the dimer's basis, Hartree-Fock, against the CPU ghosts
      !!
      !! **Only runs where terco is linked.** The CPU backend's ghost path is
      !! the reference: it zeroes a ghost's charge and keeps its basis. The
      !! failure this guards against is the one terco had before it took
      !! elements apart from charges, a ghost reaching it as a full nucleus, so
      !! the monomer in the dimer's basis came back as the whole dimer.
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings

      if (.not. terco_backend_available()) return

      settings%basis_set = "cc-pvdz"
      settings%energy_tol = 1.0e-10_dp
      settings%density_tol = 1.0e-8_dp
      call check_ghosted_against_cpu(error, settings, "RHF")
   end subroutine test_ghosted_rhf_matches_cpu

   subroutine test_ghosted_b3lyp_matches_cpu(error)
      !! The same, Kohn-Sham, so the ghosts' grid points are checked too
      !!
      !! **Only runs where terco is linked.** A ghost keeps its grid points and
      !! its Becke size, which terco reads from the element; a ghost dropped
      !! from the grid would integrate the ghost basis functions' density
      !! wrongly and show up here and not in the Hartree-Fock case.
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings

      if (.not. terco_backend_available()) return

      settings%basis_set = "cc-pvdz"
      settings%functional = "b3lyp"
      settings%grid_level = 3
      settings%energy_tol = 1.0e-10_dp
      settings%density_tol = 1.0e-8_dp
      call check_ghosted_against_cpu(error, settings, "B3LYP")
   end subroutine test_ghosted_b3lyp_matches_cpu

   subroutine test_ghosted_pcm_refused(error)
      !! A ghosted fragment under a continuum is refused, not answered
      !!
      !! The CPU backend gives a ghost centre a cavity sphere and terco gives
      !! it none, so the two would return different counterpoise corrections
      !! for the same deck. True in a stub build too, which refuses everything.
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: monomer, fragment
      type(calculation_result_t) :: result

      call water(monomer)
      call ghosted_water_dimer(monomer, fragment)
      settings%basis_set = "sto-3g"
      settings%pcm%enabled = .true.

      call run_terco_scf(settings, fragment, result)

      call check(error, result%has_error, "a ghosted fragment under a continuum "// &
                 "reached terco without an error")
      if (allocated(error)) return
      call check(error,.not. result%has_energy, "a refused ghosted continuum run "// &
                 "still reported an energy")
      if (allocated(error)) return
      if (.not. terco_backend_available()) return
      call check(error, index(result%error%get_message(), "cavity") > 0, &
                 "the refusal does not say the cavity is the reason: "// &
                 result%error%get_message())
   end subroutine test_ghosted_pcm_refused

   subroutine check_ghosted_against_cpu(error, settings, label)
      !! Water in a water dimer's basis on terco and on the CPU, to 1e-8 Eh
      !!
      !! Also checks the ghosted energy sits near the bare monomer's, below it
      !! by a basis-set superposition error of millihartree. A ghost taken as a
      !! nucleus would put it some seventy hartree away.
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t), intent(in) :: settings
      character(len=*), intent(in) :: label

      type(physical_fragment_t) :: monomer, fragment
      type(calculation_result_t) :: terco_result, cpu_result, bare_result
      real(dp) :: bsse

      call water(monomer)
      call ghosted_water_dimer(monomer, fragment)

      call run_terco_scf(settings, fragment, terco_result)
      call check(error, terco_result%has_energy, "terco returned no energy for a ghosted "// &
                 "water, "//label//": "//terco_result%error%get_message())
      if (allocated(error)) return
      call run_czt_hf(settings, fragment, cpu_result)
      call check(error, cpu_result%has_energy, "the CPU path returned no energy for a "// &
                 "ghosted water, "//label)
      if (allocated(error)) return
      ! 1e-8 is the convergence both were asked for, as in `test_matches_cpu`.
      call check(error, abs(terco_result%energy%scf - cpu_result%energy%scf) < 1.0e-8_dp, &
                 "terco and the CPU path disagree on a ghosted water, "//label)
      if (allocated(error)) return

      call run_czt_hf(settings, monomer, bare_result)
      call check(error, bare_result%has_energy, "the CPU path returned no energy for "// &
                 "the bare water, "//label)
      if (allocated(error)) return
      bsse = bare_result%energy%scf - terco_result%energy%scf
      call check(error, bsse > 0.0_dp .and. bsse < 1.0e-2_dp, &
                 "the ghosts on terco should lower the monomer by a small basis-set "// &
                 "superposition error, "//label)
   end subroutine check_ghosted_against_cpu

   subroutine ghosted_water_dimer(monomer, fragment)
      !! `monomer` with a second, ghosted copy of itself 5.5 bohr along z
      !!
      !! What VMFC builds for a monomer in a dimer's basis: the real atoms
      !! first, the partner's atoms marked ghost, the electron count the real
      !! atoms' alone.
      type(physical_fragment_t), intent(in) :: monomer
      type(physical_fragment_t), intent(out) :: fragment

      fragment%n_atoms = 2*monomer%n_atoms
      fragment%charge = monomer%charge
      fragment%multiplicity = monomer%multiplicity
      fragment%nelec = monomer%nelec
      fragment%n_caps = 0
      fragment%element_numbers = [monomer%element_numbers, monomer%element_numbers]
      allocate (fragment%coordinates(3, fragment%n_atoms))
      fragment%coordinates(:, 1:monomer%n_atoms) = monomer%coordinates
      fragment%coordinates(:, monomer%n_atoms + 1:) = monomer%coordinates
      fragment%coordinates(3, monomer%n_atoms + 1:) = &
         fragment%coordinates(3, monomer%n_atoms + 1:) + 5.5_dp
      allocate (fragment%is_ghost(fragment%n_atoms))
      fragment%is_ghost = .false.
      fragment%is_ghost(monomer%n_atoms + 1:) = .true.
   end subroutine ghosted_water_dimer

   subroutine test_missing_build(error)
      !! On a build without terco, a terco deck is told so rather than ignored
      !!
      !! Skipped rather than inverted when the library IS linked: there is no
      !! missing build to report, and asserting the opposite would only be
      !! asserting that this build is the one it is.
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      if (terco_backend_available()) return

      call water(fragment)
      settings%basis_set = "sto-3g"

      call run_terco_scf(settings, fragment, result)

      call check(error, result%has_error, &
                 "a terco deck on a build without terco reported no error")
      if (allocated(error)) return
      call check(error,.not. result%has_energy, &
                 "a build without terco still reported an energy")
      if (allocated(error)) return
      call check(error, index(result%error%get_message(), "MQC_ENABLE_TERCO") > 0, &
                 "the refusal does not say how to build the backend in")
   end subroutine test_missing_build

   subroutine test_matches_cpu(error)
      !! terco and the CPU path must answer the same number for the same basis
      !!
      !! **Only runs where terco is linked**, which is nowhere in CI. It is here
      !! rather than in `validation/` because it is the one check that says the
      !! adapter marshals correctly -- every other test in this file constrains
      !! the refusals, and a basis handed over with a transposed array or a
      !! dropped contraction would pass all of them.
      !!
      !! Hartree-Fock on purpose: no quadrature, so a disagreement is the
      !! marshalling and not a difference in grid pruning between two codes that
      !! happen to share a grid implementation. 6-31G is segmented, which the
      !! driver requires, and stops at p so the Cartesian/spherical question
      !! does not arise either.
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: terco_result, cpu_result

      if (.not. terco_backend_available()) return

      call water(fragment)
      settings%basis_set = "6-31g"
      settings%energy_tol = 1.0e-10_dp
      settings%density_tol = 1.0e-8_dp

      call run_terco_scf(settings, fragment, terco_result)
      call check(error, terco_result%has_energy, "terco returned no energy for water/6-31G")
      if (allocated(error)) return

      call run_czt_hf(settings, fragment, cpu_result)
      call check(error, cpu_result%has_energy, "the CPU path returned no energy for water/6-31G")
      if (allocated(error)) return

      ! 1e-8 is the convergence both were asked for, not a fitted bound: the two
      ! solve the same equations to that threshold, so anything larger is a
      ! different Hamiltonian rather than a different iterate.
      call check(error, abs(terco_result%energy%scf - cpu_result%energy%scf) < 1.0e-8_dp, &
                 "terco and the CPU path disagree on water/6-31G Hartree-Fock")
   end subroutine test_matches_cpu

   subroutine test_general_contraction(error)
      !! cc-pVDZ, whose s space is one shell carrying three contracted functions
      !!
      !! terco has no `NCTR_OF`, so the driver splits each column into its own
      !! shell before handing the basis over. The failure that guards against is
      !! silent: terco would otherwise build one function per shell from the
      !! first column, run a perfectly convergent SCF over a smaller basis, and
      !! report the energy as the deck's own.
      !!
      !! Cartesian because terco reads every shell that way and cc-pVDZ has d
      !! functions -- six against five is a different basis, and the driver
      !! refuses that mismatch rather than papering over it.
      !!
      !! **Only runs where terco is linked.**
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: terco_result, cpu_result

      if (.not. terco_backend_available()) return

      call water(fragment)
      settings%basis_set = "cc-pvdz"
      settings%cartesian = .true.
      settings%energy_tol = 1.0e-10_dp
      settings%density_tol = 1.0e-8_dp

      call run_terco_scf(settings, fragment, terco_result)
      call check(error, terco_result%has_energy, "terco returned no energy for water/cc-pVDZ")
      if (allocated(error)) return

      call run_czt_hf(settings, fragment, cpu_result)
      call check(error, cpu_result%has_energy, "the CPU path returned no energy for water/cc-pVDZ")
      if (allocated(error)) return

      call check(error, abs(terco_result%energy%scf - cpu_result%energy%scf) < 1.0e-8_dp, &
                 "terco and the CPU path disagree on water/cc-pVDZ -- a dropped "// &
                 "contraction column would show up exactly here")
   end subroutine test_general_contraction

   subroutine test_spherical_matches_cpu(error)
      !! cc-pVDZ with five d functions, handed to terco spherical
      !!
      !! **Only runs where terco is linked.** A spherical shell read as
      !! Cartesian is a larger basis and a lower energy, converged and silent.
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: terco_result, cpu_result

      if (.not. terco_backend_available()) return

      call water(fragment)
      settings%basis_set = "cc-pvdz"
      settings%energy_tol = 1.0e-10_dp
      settings%density_tol = 1.0e-8_dp

      call run_terco_scf(settings, fragment, terco_result)
      call check(error, terco_result%has_energy, "terco returned no spherical cc-pVDZ energy")
      if (allocated(error)) return

      call run_czt_hf(settings, fragment, cpu_result)
      call check(error, cpu_result%has_energy, "the CPU path returned no spherical energy")
      if (allocated(error)) return

      call check(error, abs(terco_result%energy%scf - cpu_result%energy%scf) < 1.0e-8_dp, &
                 "terco and the CPU path disagree on spherical water/cc-pVDZ")
   end subroutine test_spherical_matches_cpu

   subroutine test_gradient_matches_cpu(error)
      !! terco's analytic Hartree-Fock gradient against the CPU path's
      !!
      !! **Only runs where terco is linked.** The two converge the same
      !! equations, so their gradients agree to the convergence asked for;
      !! measured 1e-10 on ammonia/6-31G*.
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: terco_result, cpu_result

      if (.not. terco_backend_available()) return

      call water(fragment)
      settings%basis_set = "6-31g"
      settings%energy_tol = 1.0e-10_dp
      settings%density_tol = 1.0e-8_dp

      call run_terco_scf(settings, fragment, terco_result, want_gradient=.true.)
      call check(error, terco_result%has_gradient, "terco returned no gradient for water/6-31G")
      if (allocated(error)) return

      call run_czt_hf(settings, fragment, cpu_result, want_gradient=.true.)
      call check(error, cpu_result%has_gradient, "the CPU path returned no gradient for water/6-31G")
      if (allocated(error)) return

      call check(error, maxval(abs(terco_result%gradient - cpu_result%gradient)) < 1.0e-7_dp, &
                 "terco and the CPU path disagree on the water/6-31G Hartree-Fock gradient")
   end subroutine test_gradient_matches_cpu

   subroutine test_fitted_matches_cpu(error)
      !! Density-fitted Hartree-Fock on both, over the same auxiliary basis
      !!
      !! **Only runs where terco is linked.** terco reads the fitting basis from
      !! its Basis Set Exchange file and mqc from the same file, so a fit that
      !! reached the device with a different basis shows up here.
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: terco_result, cpu_result
      character(len=96) :: line

      if (.not. terco_backend_available()) return

      ! Spherical: def2-universal-jkfit carries d functions and above, so a
      ! fitting basis that reached the device Cartesian would show up here.
      call water(fragment)
      settings%basis_set = "6-31g"
      settings%density_fitting = .true.
      settings%aux_basis_set = "def2-universal-jkfit"
      settings%aux_basis_named = .true.
      settings%energy_tol = 1.0e-10_dp
      settings%density_tol = 1.0e-8_dp

      call run_terco_scf(settings, fragment, terco_result)
      call check(error, terco_result%has_energy, "terco returned no fitted energy")
      if (allocated(error)) return

      call run_czt_hf(settings, fragment, cpu_result)
      call check(error, cpu_result%has_energy, "the CPU path returned no fitted energy")
      if (allocated(error)) return

      write (line, "(2(a,f20.12))") " terco ", terco_result%energy%scf, " cpu ", &
         cpu_result%energy%scf
      call check(error, abs(terco_result%energy%scf - cpu_result%energy%scf) < 1.0e-8_dp, &
                 "terco and the CPU path disagree on fitted water/6-31G Hartree-Fock:"// &
                 trim(line))
   end subroutine test_fitted_matches_cpu

   subroutine water(fragment)
      !! One water, near equilibrium, in Bohr
      type(physical_fragment_t), intent(out) :: fragment

      fragment%n_atoms = 3
      fragment%charge = 0
      fragment%multiplicity = 1
      fragment%nelec = 10
      fragment%n_caps = 0

      allocate (fragment%element_numbers(3))
      allocate (fragment%coordinates(3, 3))

      fragment%element_numbers = [8, 1, 1]
      fragment%coordinates(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      fragment%coordinates(:, 2) = [0.0_dp, 1.4304_dp, 1.1071_dp]
      fragment%coordinates(:, 3) = [0.0_dp, -1.4304_dp, 1.1071_dp]
   end subroutine water

end module test_mqc_terco_backend

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_terco_backend, only: collect_mqc_terco_backend_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_terco_backend", collect_mqc_terco_backend_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
