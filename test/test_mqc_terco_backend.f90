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
   !! a gradient never yields an energy" holds in a stub build because the stub
   !! refuses the whole calculation, and in a real build because the driver
   !! refuses the gradient. The reason differs; the guarantee does not.
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
                  new_unittest("gradient_never_yields_energy", test_gradient_refused), &
                  new_unittest("correlated_never_yields_energy", test_correlated_refused), &
                  new_unittest("missing_build_is_reported", test_missing_build), &
                  new_unittest("agrees_with_the_cpu_path", test_matches_cpu) &
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

   subroutine test_gradient_refused(error)
      !! terco has energy kernels only, and must never answer a gradient request
      !!
      !! The failure this guards against is not an exception -- it is a geometry
      !! optimisation that receives a result carrying an energy, no gradient and
      !! no error, and steps on whatever `result%gradient` happened to hold.
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call water(fragment)
      settings%basis_set = "sto-3g"

      call run_terco_scf(settings, fragment, result, want_gradient=.true.)

      call check(error, result%has_error, "asking terco for a gradient did not set an error")
      if (allocated(error)) return
      call check(error,.not. result%has_energy, "a refused gradient still reported an energy")
      if (allocated(error)) return
      call check(error,.not. allocated(result%gradient), &
                 "a refused gradient request still allocated a gradient")
   end subroutine test_gradient_refused

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
