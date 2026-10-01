!! SA-CASSCF gradients through the driver path (`mqc_czt_bridge`'s `run_czt_mcscf`)
module test_mqc_sa_gradient_driver
   !! Phase 6 of `SA_CASSCF_GRADIENT_PLAN.md`: the Z-vector gradient
   !! (`czt_sa_casscf_gradients`, phase 5) is correct on its own -- this
   !! exercises the bridge dispatch a Gradient driver deck actually reaches
   !! (`run_czt_mcscf`), which phase 6 wires up for the first time: the
   !! top-level gradient, `result%mcscf_state_gradients` and
   !! `result%mcscf_gradient_roots`, and the refusals that guard
   !! `keywords.mcscf.gradient_roots` and unequal weights.
   !!
   !! LiH/STO-3G SA-2-CAS(2,2), same geometry as `test_mqc_sa_casscf.f90` and
   !! `test_mqc_sa_gradient_pyscf_long.f90`.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_czt_bridge, only: run_czt_mcscf
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mcscf, only: run_czt_casscf, casscf_result_t
   use mqc_czt_sa_gradient, only: czt_sa_casscf_gradient
   implicit none
   private

   public :: collect_mqc_sa_gradient_driver_tests

   real(dp), parameter :: LIH_Z_TO_H = 3.0139241961656_dp
      !! Bohr; `1.5949` Angstrom, matching `tools/sa_casscf/lih.xyz`.
   real(dp), parameter :: WEIGHTS(2) = [0.5_dp, 0.5_dp]

contains

   subroutine collect_mqc_sa_gradient_driver_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("driver_gradient_matches_direct_call", &
                               test_driver_matches_direct), &
                  new_unittest("driver_top_level_is_weighted_sum", &
                               test_weighted_sum), &
                  new_unittest("gradient_roots_with_one_state_refused", &
                               test_roots_needs_sa), &
                  new_unittest("gradient_roots_off_a_gradient_driver_refused", &
                               test_roots_needs_gradient_driver), &
                  new_unittest("unequal_weights_gradient_refused", &
                               test_unequal_weights_refused) &
                  ]
   end subroutine collect_mqc_sa_gradient_driver_tests

   subroutine lih_fragment(fragment)
      type(physical_fragment_t), intent(out) :: fragment

      fragment%n_atoms = 2
      fragment%charge = 0
      fragment%multiplicity = 1
      fragment%nelec = 4
      fragment%n_caps = 0
      allocate (fragment%element_numbers(2))
      allocate (fragment%coordinates(3, 2))
      fragment%element_numbers = [3, 1]
      fragment%coordinates(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      fragment%coordinates(:, 2) = [0.0_dp, 0.0_dp, LIH_Z_TO_H]
   end subroutine lih_fragment

   subroutine sa2_settings(settings, want_roots)
      !! The settings every test here starts from
      type(cuest_scf_settings_t), intent(out) :: settings
      logical, intent(in) :: want_roots

      settings%basis_set = "sto-3g"
      settings%energy_tol = 1.0e-12_dp
      settings%density_tol = 1.0e-10_dp
      settings%mcscf%n_active_electrons = 2
      settings%mcscf%n_active_orbitals = 2
      settings%mcscf%max_macro_iter = 400
      settings%mcscf%orbital_convergence = 1.0e-10_dp
      settings%mcscf%n_states = 2
      settings%mcscf%state_weights = WEIGHTS
      if (want_roots) settings%mcscf%gradient_roots = [1, 2]
   end subroutine sa2_settings

   subroutine test_driver_matches_direct(error)
      !! Every root the bridge stores matches a direct
      !! `czt_sa_casscf_gradient` call on the same converged CASSCF
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casscf_result_t) :: casscf
      type(error_t) :: err
      real(dp), allocatable :: g_direct(:, :)
      real(dp) :: worst
      integer :: root

      call lih_fragment(fragment)
      call sa2_settings(settings, want_roots=.true.)

      call run_czt_mcscf(settings, fragment, result, want_gradient=.true.)
      call check(error,.not. result%has_error, err_message(result))
      if (allocated(error)) return
      call check(error, result%has_gradient, "the driver should return a gradient")
      if (allocated(error)) return
      call check(error, allocated(result%mcscf_state_gradients), &
                 "per-root gradients should be stored on the result")
      if (allocated(error)) return
      call check(error, size(result%mcscf_gradient_roots), 2)
      if (allocated(error)) return
      call check(error, result%mcscf_gradient_roots(1), 1)
      if (allocated(error)) return
      call check(error, result%mcscf_gradient_roots(2), 2)
      if (allocated(error)) return

      ! Independently reconverge the same SA-CASSCF and ask
      ! `czt_sa_casscf_gradient` for each root directly, off the same
      ! converged orbitals and CI vectors -- what the bridge's dispatch is
      ! meant to reproduce.
      call build_czt_molecule(fragment%element_numbers, ["Li", "H "], &
                              fragment%coordinates, "sto-3g", mol, err)
      call run_czt_rhf(mol, fragment%nelec, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call run_czt_casscf(mol, scf%orbitals, 1, 2, 1, 1, casscf, err, &
                          max_iterations=400, gradient_tol=1.0e-10_dp, n_states=2, &
                          weights=WEIGHTS)
      call check(error,.not. err%has_error() .and. casscf%converged, &
                 "the independent SA-CASSCF should converge")
      if (allocated(error)) return

      do root = 1, 2
         call czt_sa_casscf_gradient(mol, casscf%orbitals, 1, 2, 1, 1, &
                                     casscf%ci_vectors, casscf%energies, WEIGHTS, &
                                     root, g_direct, err)
         call check(error,.not. err%has_error(), "the direct root gradient should build")
         if (allocated(error)) return
         worst = maxval(abs(result%mcscf_state_gradients(:, :, root) - g_direct))
         call check(error, worst < 1.0e-8_dp, &
                    "the driver's stored root gradient should match the direct call")
         if (allocated(error)) return
      end do
      call mol%destroy()
   end subroutine test_driver_matches_direct

   subroutine test_weighted_sum(error)
      !! The top-level gradient is `dE_SA/dR`: with every root requested and
      !! equal weights, that is exactly `sum_I w_I g_I`
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result
      real(dp), allocatable :: weighted(:, :)
      real(dp) :: worst

      call lih_fragment(fragment)
      call sa2_settings(settings, want_roots=.true.)

      call run_czt_mcscf(settings, fragment, result, want_gradient=.true.)
      call check(error,.not. result%has_error, err_message(result))
      if (allocated(error)) return

      weighted = WEIGHTS(1)*result%mcscf_state_gradients(:, :, 1) + &
                 WEIGHTS(2)*result%mcscf_state_gradients(:, :, 2)
      worst = maxval(abs(result%gradient - weighted))
      call check(error, worst < 1.0e-9_dp, &
                 "sum_I w_I g_I should equal the top-level SA gradient")
   end subroutine test_weighted_sum

   subroutine test_roots_needs_sa(error)
      !! `gradient_roots` with `n_states = 1` is refused, not silently ignored
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call lih_fragment(fragment)
      settings%basis_set = "sto-3g"
      settings%mcscf%n_active_electrons = 2
      settings%mcscf%n_active_orbitals = 2
      settings%mcscf%gradient_roots = [1]
      ! n_states left at its default of 1.

      call run_czt_mcscf(settings, fragment, result, want_gradient=.true.)
      call check(error, result%has_error, &
                 "gradient_roots with n_states = 1 should be refused")
   end subroutine test_roots_needs_sa

   subroutine test_roots_needs_gradient_driver(error)
      !! `gradient_roots` off a Gradient driver is refused, not silently
      !! ignored
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call lih_fragment(fragment)
      call sa2_settings(settings, want_roots=.true.)

      call run_czt_mcscf(settings, fragment, result, want_gradient=.false.)
      call check(error, result%has_error, &
                 "gradient_roots without a Gradient driver should be refused")
      if (allocated(error)) return

      call result%reset()
      call run_czt_mcscf(settings, fragment, result)
      call check(error, result%has_error, &
                 "gradient_roots with want_gradient absent should be refused too")
   end subroutine test_roots_needs_gradient_driver

   subroutine test_unequal_weights_refused(error)
      !! A state-averaged CASSCF gradient needs equal weights
      type(error_type), allocatable, intent(out) :: error
      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment
      type(calculation_result_t) :: result

      call lih_fragment(fragment)
      call sa2_settings(settings, want_roots=.false.)
      settings%mcscf%state_weights = [0.7_dp, 0.3_dp]

      call run_czt_mcscf(settings, fragment, result, want_gradient=.true.)
      call check(error, result%has_error, &
                 "unequal weights and a Gradient driver should be refused")
   end subroutine test_unequal_weights_refused

   function err_message(result) result(msg)
      type(calculation_result_t), intent(in) :: result
      character(len=:), allocatable :: msg
      if (result%has_error) then
         msg = result%error%get_message()
      else
         msg = ""
      end if
   end function err_message

end module test_mqc_sa_gradient_driver

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_sa_gradient_driver, only: collect_mqc_sa_gradient_driver_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_sa_gradient_driver", &
                               collect_mqc_sa_gradient_driver_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
