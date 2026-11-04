!! The fitted correlation Hessian as a deck reaches it, through the bridge
module test_mqc_czt_ri_hessian_bridge
   !! What `test_mqc_ri_mp2_hessian` and `test_mqc_ri_mp2_hessian_dh` cannot
   !! reach: the routine is right, but whether a run takes it, declines it and
   !! says why, or hands a double hybrid's dipole derivatives back is decided in
   !! `run_czt_hf`, one level up from anything those call.
   !!
   !! **The declined metric.** `fixtures/test-dup-rifit.json` is cc-pVDZ-RIFIT
   !! for H and O with the first oxygen `s` shell listed twice. The span of the
   !! fitting space is unchanged, so the fitted energy is the same as with the
   !! real set, but the metric has an exact null mode, which is the one
   !! condition under which the analytic fitted Hessian declines. The tests run
   !! it through the method a deck's `calc_hessian` resolves to -- the only
   !! place the semi-numerical Hessian is entered from -- and require that the
   !! run succeeds, that a Hessian comes back, that it is the Hessian of the
   !! real set to within the finite-difference error, and that the one warning
   !! the run logged names the metric rather than some other reason for the
   !! fallback. Without that last check a deck that fell back for an unrelated
   !! reason would pass.
   !!
   !! A fixture file rather than a library change: nothing in `basis_sets/` can
   !! hold it, since the JSON files under `basis_sets/` are derived from the Basis Set Exchange
   !! bundle at configure time and removed when it is not in it. The test copies
   !! the fixture into `./basis_sets/` of its working directory -- the first
   !! place a basis is looked for -- and removes it again.
   !!
   !! **The dipole derivatives.** A double hybrid's dipole derivative is the
   !! reference's plus the perturbative term's share, and the latter exists in
   !! the bridge only as an addition to the former. The fitted run is compared
   !! with the conventional one, the same deck without a named auxiliary basis,
   !! and the two must agree to the fitting error and not to the last digit.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger, info_level, performance_level
   use mqc_czt_bridge, only: run_czt_hf, xc_available
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_method_hf, only: hf_method_t
   use mqc_method_dft, only: dft_method_t
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   implicit none
   private

   public :: collect_mqc_czt_ri_hessian_bridge_tests

   real(dp), parameter :: WATER_BOHR(3, 3) = reshape( &
                          [0.05_dp, 0.02_dp, -0.03_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.04_dp, 1.4308_dp, 1.1078_dp], [3, 3])
   character(len=*), parameter :: FIXTURE_AUX = "test-dup-rifit"
   character(len=*), parameter :: REAL_AUX = "cc-pvdz-rifit"
   character(len=*), parameter :: LOG_PATH = "output_ri_hessian_bridge.log"
   character(len=*), parameter :: METRIC_REASON = "the auxiliary metric drops"
   real(dp), parameter :: FD_TOL_MP2 = 1.0e-4_dp
      !! Largest element of the differenced Hessian against the analytic one of
      !! the real set. Measured 1.4e-5: central differences of an analytic
      !! gradient, which itself sits about 1e-7 from the derivative of the energy
      !! it reports. A wrong block is order 1e-2 or more.
   real(dp), parameter :: FD_TOL_DH = 1.0e-3_dp
      !! The same for a double hybrid, measured 1.8e-4: its differenced gradient
      !! moves the quadrature grid with the nuclei and its analytic Hessian
      !! holds it fixed.
   real(dp), parameter :: FITTING_LIMIT = 5.0e-4_dp
      !! Relative, fitted against conventional. Measured 2.4e-5 for the Hessian
      !! with cc-pVDZ-RIFIT on 6-31G water; a wrong block is order one.

contains

   subroutine collect_mqc_czt_ri_hessian_bridge_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("ri_mp2_dropped_mode_takes_the_semi_numerical_hessian", &
                               mp2_declined), &
                  new_unittest("b2plyp_dropped_mode_takes_the_semi_numerical_hessian", &
                               dh_declined), &
                  new_unittest("b2plyp_fitted_hessian_and_dipole_derivatives", &
                               dh_fitted_against_conventional) &
                  ]
   end subroutine collect_mqc_czt_ri_hessian_bridge_tests

   subroutine water_fragment(fragment)
      type(physical_fragment_t), intent(out) :: fragment

      fragment%n_atoms = 3
      fragment%charge = 0
      fragment%multiplicity = 1
      fragment%nelec = 10
      fragment%n_caps = 0
      allocate (fragment%element_numbers(3), fragment%coordinates(3, 3))
      fragment%element_numbers = [8, 1, 1]
      fragment%coordinates = WATER_BOHR
   end subroutine water_fragment

   subroutine install_fixture()
      call execute_command_line("mkdir -p basis_sets && cp fixtures/test-dup-rifit.json "// &
                                "basis_sets/test-dup-rifit.json")
   end subroutine install_fixture

   subroutine remove_fixture()
      call execute_command_line("rm -f basis_sets/test-dup-rifit.json")
   end subroutine remove_fixture

   subroutine log_to_file()
      call logger%configure_file_output(LOG_PATH, level=info_level)
   end subroutine log_to_file

   function log_contains(text) result(found)
      !! Whether any line of the log written since `log_to_file` has `text` in it
      character(len=*), intent(in) :: text
      logical :: found

      character(len=512) :: line
      integer :: unit, stat

      found = .false.
      call logger%close_log_file()
      open (newunit=unit, file=LOG_PATH, status="old", action="read", iostat=stat)
      if (stat /= 0) return
      do
         read (unit, "(a)", iostat=stat) line
         if (stat /= 0) exit
         if (index(line, text) > 0) then
            found = .true.
            exit
         end if
      end do
      close (unit, status="delete")
   end function log_contains

   subroutine compare_hessians(error, test, reference, tol, what)
      type(error_type), allocatable, intent(out) :: error
      real(dp), intent(in) :: test(:, :), reference(:, :)
      real(dp), intent(in) :: tol
      character(len=*), intent(in) :: what

      real(dp) :: worst

      call check(error, all(shape(test) == shape(reference)), what//": the two Hessians differ in shape")
      if (allocated(error)) return
      worst = maxval(abs(test - reference))
      call check(error, worst < tol, what//": largest element differs by "// &
                 trim(real_text(worst))//", limit "//trim(real_text(tol)))
      if (allocated(error)) return
      worst = maxval(abs(test - transpose(test)))
      call check(error, worst < 1.0e-6_dp, what//": the Hessian is not symmetric, by "// &
                 trim(real_text(worst)))
   end subroutine compare_hessians

   function real_text(x) result(text)
      real(dp), intent(in) :: x
      character(len=32) :: text

      write (text, "(es10.3)") x
   end function real_text

   subroutine mp2_options(method, aux)
      type(hf_method_t), intent(inout) :: method
      character(len=*), intent(in) :: aux

      method%options%basis_set = "6-31g"
      method%options%run_mp2 = .true.
      method%options%corr_density_fitting = .true.
      method%options%aux_basis_set = aux
      method%options%aux_basis_named = .true.
      method%options%energy_tol = 1.0e-12_dp
      method%options%density_tol = 1.0e-10_dp
      method%options%max_iter = 200
   end subroutine mp2_options

   subroutine mp2_declined(error)
      !! RI-MP2 through `hf_method_t%calc_hessian`: the metric check ahead of the
      !! MP2 Hessian
      type(error_type), allocatable, intent(out) :: error

      type(hf_method_t) :: fixture, reference
      type(calculation_result_t) :: declined, analytic
      type(physical_fragment_t) :: fragment
      integer :: threads

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call water_fragment(fragment)
      call mp2_options(reference, REAL_AUX)
      call reference%calc_hessian(fragment, analytic)

      call install_fixture()
      call log_to_file()
      call mp2_options(fixture, FIXTURE_AUX)
      call fixture%calc_hessian(fragment, declined)
      call remove_fixture()
      call omp_set_num_threads(threads)

      call check(error,.not. analytic%has_error .and. analytic%has_hessian, &
                 "control: the real auxiliary basis gave no Hessian")
      if (allocated(error)) return
      call check(error,.not. declined%has_error, "the dropped-mode run failed")
      if (allocated(error)) return
      call check(error, declined%has_hessian, "the dropped-mode run returned no Hessian")
      if (allocated(error)) return
      call check(error, log_contains(METRIC_REASON), &
                 "the run did not log the metric as its reason for the fallback")
      if (allocated(error)) return
      call compare_hessians(error, declined%hessian, analytic%hessian, FD_TOL_MP2, "RI-MP2")
   end subroutine mp2_declined

   subroutine dh_options(method, aux, named)
      type(dft_method_t), intent(inout) :: method
      character(len=*), intent(in) :: aux
      logical, intent(in) :: named

      method%options%basis_set = "6-31g"
      method%options%functional = "b2plyp"
      method%options%grid_level = 3
      method%options%freeze_core = .false.
      method%options%aux_basis_set = aux
      method%options%aux_basis_named = named
      method%options%energy_tol = 1.0e-12_dp
      method%options%density_tol = 1.0e-10_dp
      method%options%max_iter = 200
   end subroutine dh_options

   subroutine dh_declined(error)
      !! The double hybrid's own metric check, which runs before the reference
      !! Hessian is built
      type(error_type), allocatable, intent(out) :: error

      type(dft_method_t) :: fixture, reference
      type(calculation_result_t) :: declined, analytic
      type(physical_fragment_t) :: fragment
      integer :: threads

      if (.not. xc_available()) return
      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call water_fragment(fragment)
      call dh_options(reference, REAL_AUX, .true.)
      call reference%calc_hessian(fragment, analytic)

      call install_fixture()
      call log_to_file()
      call dh_options(fixture, FIXTURE_AUX, .true.)
      call fixture%calc_hessian(fragment, declined)
      call remove_fixture()
      call omp_set_num_threads(threads)

      call check(error,.not. analytic%has_error .and. analytic%has_hessian, &
                 "control: the real auxiliary basis gave no Hessian")
      if (allocated(error)) return
      call check(error,.not. declined%has_error, "the dropped-mode run failed")
      if (allocated(error)) return
      call check(error, declined%has_hessian, "the dropped-mode run returned no Hessian")
      if (allocated(error)) return
      call check(error, log_contains(METRIC_REASON), &
                 "the run did not log the metric as its reason for the fallback")
      if (allocated(error)) return
      call compare_hessians(error, declined%hessian, analytic%hessian, FD_TOL_DH, "B2PLYP")
   end subroutine dh_declined

   subroutine dh_run(aux, named, result)
      character(len=*), intent(in) :: aux
      logical, intent(in) :: named
      type(calculation_result_t), intent(out) :: result

      type(cuest_scf_settings_t) :: settings
      type(physical_fragment_t) :: fragment

      call water_fragment(fragment)
      settings%basis_set = "6-31g"
      settings%functional = "b2plyp"
      settings%grid_level = 3
      settings%freeze_core = .false.
      settings%aux_basis_set = aux
      settings%aux_basis_named = named
      settings%energy_tol = 1.0e-12_dp
      settings%density_tol = 1.0e-10_dp
      settings%grad_tol = 1.0e-9_dp
      settings%max_iter = 200
      call run_czt_hf(settings, fragment, result, want_gradient=.true., want_hessian=.true.)
   end subroutine dh_run

   subroutine dh_fitted_against_conventional(error)
      !! The Hessian and the dipole derivatives of a double hybrid with a named
      !! auxiliary basis, against the same deck with none
      type(error_type), allocatable, intent(out) :: error

      type(calculation_result_t) :: fitted, exact
      real(dp) :: rel, diff
      integer :: threads

      if (.not. xc_available()) return
      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call dh_run(REAL_AUX, .true., fitted)
      call dh_run("def2-universal-jkfit", .false., exact)
      call omp_set_num_threads(threads)

      call check(error,.not. fitted%has_error .and. .not. exact%has_error, "a run failed")
      if (allocated(error)) return
      call check(error, fitted%has_hessian .and. exact%has_hessian, &
                 "a run returned no analytic Hessian")
      if (allocated(error)) return
      call check(error, fitted%has_dipole_derivatives .and. exact%has_dipole_derivatives, &
                 "a run returned no dipole derivatives")
      if (allocated(error)) return

      diff = norm2(fitted%hessian - exact%hessian)
      rel = diff/norm2(exact%hessian)
      call check(error, rel < FITTING_LIMIT .and. diff > 1.0e-9_dp, &
                 "Hessian, fitted against conventional: relative "//trim(real_text(rel)))
      if (allocated(error)) return
      diff = norm2(fitted%dipole_derivatives - exact%dipole_derivatives)
      rel = diff/norm2(exact%dipole_derivatives)
      call check(error, rel < FITTING_LIMIT .and. diff > 1.0e-9_dp, &
                 "dipole derivatives, fitted against conventional: relative "//trim(real_text(rel)))
   end subroutine dh_fitted_against_conventional

end module test_mqc_czt_ri_hessian_bridge

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_ri_hessian_bridge, only: collect_mqc_czt_ri_hessian_bridge_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_ri_hessian_bridge", collect_mqc_czt_ri_hessian_bridge_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
