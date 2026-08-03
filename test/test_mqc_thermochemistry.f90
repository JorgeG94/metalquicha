module test_mqc_thermochemistry
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_thermochemistry, only: thermochemistry_result_t, &
                                  compute_moments_of_inertia, compute_rotational_constants, &
                                  compute_zpe, compute_translational_thermo, &
                                  compute_rotational_thermo, compute_vibrational_thermo, &
                                  compute_electronic_entropy, compute_thermochemistry, &
                                  print_thermochemistry
   use mqc_physical_constants, only: BOHR_TO_ANGSTROM, ANGSTROM_TO_BOHR, &
                                     HARTREE_TO_KCALMOL, CM1_TO_KELVIN, KB_HARTREE
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_thermochemistry_tests

   ! Tolerance for floating point comparisons
   real(dp), parameter :: tol = 1.0e-6_dp
   real(dp), parameter :: tol_loose = 1.0e-3_dp  ! For comparing with Gaussian

contains

   subroutine collect_mqc_thermochemistry_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("zpe_diatomic", test_zpe_diatomic), &
                  new_unittest("moments_of_inertia_h2o", test_moments_of_inertia_h2o), &
                  new_unittest("rotational_constants", test_rotational_constants), &
                  new_unittest("linear_molecule_detection", test_linear_molecule_detection), &
                  new_unittest("translational_thermo", test_translational_thermo), &
                  new_unittest("vibrational_thermo", test_vibrational_thermo), &
                  new_unittest("electronic_entropy", test_electronic_entropy), &
                  new_unittest("full_thermochemistry_h2o", test_full_thermochemistry_h2o) &
                  ]
   end subroutine collect_mqc_thermochemistry_tests

   subroutine test_zpe_diatomic(error)
      !! Test ZPE calculation for a simple diatomic (H2-like)
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: frequencies(1)
      real(dp) :: zpe_hartree, zpe_kcalmol
      real(dp) :: expected_zpe_hartree
      integer :: n_real

      ! H2 stretching frequency ~4400 cm^-1
      frequencies(1) = 4400.0_dp

      call compute_zpe(frequencies, 1, n_real, zpe_hartree, zpe_kcalmol)

      ! ZPE = 0.5 * h * c * nu = 0.5 * nu * CM1_TO_KELVIN * KB_HARTREE
      expected_zpe_hartree = 0.5_dp*4400.0_dp*CM1_TO_KELVIN*KB_HARTREE

      call check(error, n_real == 1, "Should count 1 real frequency")
      if (allocated(error)) return

      call check(error, abs(zpe_hartree - expected_zpe_hartree) < tol, &
                 "ZPE calculation should match expected value")
      if (allocated(error)) return

      ! Check kcal/mol conversion
      call check(error, abs(zpe_kcalmol - zpe_hartree*HARTREE_TO_KCALMOL) < tol, &
                 "ZPE kcal/mol conversion should be correct")

   end subroutine test_zpe_diatomic

   subroutine test_moments_of_inertia_h2o(error)
      !! Test moments of inertia for water molecule
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: coords(3, 3)
      integer :: atomic_numbers(3)
      real(dp) :: center_of_mass(3)
      real(dp) :: moments(3)
      real(dp) :: principal_axes(3, 3)
      logical :: is_linear
      real(dp) :: total_mass

      ! Water molecule geometry (approximately, in Angstrom, then convert to Bohr)
      ! O at origin, H atoms at typical positions
      ! O-H distance ~0.96 A, H-O-H angle ~104.5 degrees
      atomic_numbers = [8, 1, 1]  ! O, H, H

      ! Geometry in Angstrom
      coords(:, 1) = [0.0_dp, 0.0_dp, 0.1173_dp]  ! O
      coords(:, 2) = [0.7572_dp, 0.0_dp, -0.4692_dp]  ! H
      coords(:, 3) = [-0.7572_dp, 0.0_dp, -0.4692_dp]  ! H

      ! Convert to Bohr
      coords = coords*ANGSTROM_TO_BOHR

      call compute_moments_of_inertia(coords, atomic_numbers, 3, &
                                      center_of_mass, moments, principal_axes, is_linear, total_mass)

      ! Check total mass (O=15.999, H=1.008 x2 ~ 18.015)
      call check(error, abs(total_mass - 18.015_dp) < 0.01_dp, &
                 "Total mass should be ~18.015 amu for H2O")
      if (allocated(error)) return

      ! Water is NOT linear
      call check(error,.not. is_linear, "H2O should not be detected as linear")
      if (allocated(error)) return

      ! Moments should be positive and in ascending order
      call check(error, moments(1) > 0.0_dp .and. moments(2) > 0.0_dp .and. moments(3) > 0.0_dp, &
                 "All moments should be positive for nonlinear molecule")
      if (allocated(error)) return

      call check(error, moments(1) <= moments(2) .and. moments(2) <= moments(3), &
                 "Moments should be in ascending order")

   end subroutine test_moments_of_inertia_h2o

   subroutine test_rotational_constants(error)
      !! Test rotational constant calculation
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: moments(3)
      real(dp) :: rot_const(3)
      real(dp) :: expected_B
      logical :: is_linear

      ! For a simple test: I = 10 amu*A^2 should give B = 505379.07 / 10 = 50537.907 GHz
      moments = [10.0_dp, 10.0_dp, 10.0_dp]
      is_linear = .false.

      call compute_rotational_constants(moments, is_linear, rot_const)

      expected_B = 505379.07_dp/10.0_dp

      call check(error, abs(rot_const(1) - expected_B) < 0.1_dp, &
                 "Rotational constant calculation should match expected value")

   end subroutine test_rotational_constants

   subroutine test_linear_molecule_detection(error)
      !! Test linear molecule detection (CO2-like geometry)
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: coords(3, 3)
      integer :: atomic_numbers(3)
      real(dp) :: center_of_mass(3)
      real(dp) :: moments(3)
      real(dp) :: principal_axes(3, 3)
      logical :: is_linear
      real(dp) :: total_mass

      ! CO2 molecule: C at origin, O atoms on z-axis
      atomic_numbers = [6, 8, 8]  ! C, O, O

      ! Linear geometry along z-axis (in Angstrom, then convert)
      coords(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]      ! C
      coords(:, 2) = [0.0_dp, 0.0_dp, 1.16_dp]     ! O
      coords(:, 3) = [0.0_dp, 0.0_dp, -1.16_dp]    ! O

      ! Convert to Bohr
      coords = coords*ANGSTROM_TO_BOHR

      call compute_moments_of_inertia(coords, atomic_numbers, 3, &
                                      center_of_mass, moments, principal_axes, is_linear, total_mass)

      ! CO2 IS linear
      call check(error, is_linear, "CO2 should be detected as linear")
      if (allocated(error)) return

      ! For linear molecule, smallest moment should be ~0
      call check(error, moments(1) < 1.0e-6_dp, &
                 "Smallest moment should be ~0 for linear molecule")

   end subroutine test_linear_molecule_detection

   subroutine test_translational_thermo(error)
      !! Test translational thermodynamic contributions
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: E, S, Cv
      real(dp) :: total_mass, temperature, pressure
      real(dp) :: R_calmolk

      R_calmolk = 1.98720425864_dp  ! cal/(mol*K)

      total_mass = 18.015_dp  ! Water
      temperature = 298.15_dp
      pressure = 1.0_dp

      call compute_translational_thermo(total_mass, temperature, pressure, E, S, Cv)

      ! Cv_trans should be exactly 3/2 * R
      call check(error, abs(Cv - 1.5_dp*R_calmolk) < tol, &
                 "Translational Cv should be 3/2 R")
      if (allocated(error)) return

      ! E_trans should be positive
      call check(error, E > 0.0_dp, "Translational energy should be positive")
      if (allocated(error)) return

      ! S_trans should be positive and reasonable (typically 30-40 cal/mol/K for small molecules)
      call check(error, S > 20.0_dp .and. S < 50.0_dp, &
                 "Translational entropy should be in reasonable range")

   end subroutine test_translational_thermo

   subroutine test_vibrational_thermo(error)
      !! Test vibrational thermodynamic contributions
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: frequencies(3)
      real(dp) :: E, S, Cv
      real(dp) :: temperature

      ! Use H2O-like frequencies (approximate, cm^-1)
      frequencies = [1595.0_dp, 3657.0_dp, 3756.0_dp]
      temperature = 298.15_dp

      call compute_vibrational_thermo(frequencies, 3, temperature, E, S, Cv)

      ! At room temperature, high-frequency modes contribute little
      ! Vibrational entropy should be small but positive
      call check(error, S >= 0.0_dp, "Vibrational entropy should be non-negative")
      if (allocated(error)) return

      ! Vibrational Cv should be non-negative
      call check(error, Cv >= 0.0_dp, "Vibrational Cv should be non-negative")
      if (allocated(error)) return

      ! E_vib (thermal, excluding ZPE) should be small at 298K for high-freq modes
      call check(error, E >= 0.0_dp, "Vibrational thermal energy should be non-negative")

   end subroutine test_vibrational_thermo

   subroutine test_electronic_entropy(error)
      !! Test electronic entropy calculation
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: S_elec
      real(dp) :: R_calmolk

      R_calmolk = 1.98720425864_dp

      ! Singlet (multiplicity = 1): S = R * ln(1) = 0
      call compute_electronic_entropy(1, S_elec)
      call check(error, abs(S_elec) < tol, "Singlet should have zero electronic entropy")
      if (allocated(error)) return

      ! Doublet (multiplicity = 2): S = R * ln(2)
      call compute_electronic_entropy(2, S_elec)
      call check(error, abs(S_elec - R_calmolk*log(2.0_dp)) < tol, &
                 "Doublet electronic entropy should be R*ln(2)")
      if (allocated(error)) return

      ! Triplet (multiplicity = 3): S = R * ln(3)
      call compute_electronic_entropy(3, S_elec)
      call check(error, abs(S_elec - R_calmolk*log(3.0_dp)) < tol, &
                 "Triplet electronic entropy should be R*ln(3)")

   end subroutine test_electronic_entropy

   subroutine test_full_thermochemistry_h2o(error)
      !! Test full thermochemistry calculation for water
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: coords(3, 3)
      integer :: atomic_numbers(3)
      real(dp) :: frequencies(3)
      type(thermochemistry_result_t) :: result

      ! Water molecule geometry (in Bohr)
      atomic_numbers = [8, 1, 1]

      ! Geometry in Angstrom, then convert
      coords(:, 1) = [0.0_dp, 0.0_dp, 0.1173_dp]
      coords(:, 2) = [0.7572_dp, 0.0_dp, -0.4692_dp]
      coords(:, 3) = [-0.7572_dp, 0.0_dp, -0.4692_dp]
      coords = coords*ANGSTROM_TO_BOHR

      ! H2O frequencies (approximate, cm^-1)
      frequencies = [1595.0_dp, 3657.0_dp, 3756.0_dp]

      call compute_thermochemistry(coords, atomic_numbers, frequencies, 3, 3, result)

      ! Check basic sanity of results
      call check(error, result%temperature == 298.15_dp, &
                 "Default temperature should be 298.15 K")
      if (allocated(error)) return

      call check(error, result%pressure == 1.0_dp, &
                 "Default pressure should be 1.0 atm")
      if (allocated(error)) return

      call check(error, result%n_real_freqs == 3, &
                 "Should have 3 real frequencies for H2O")
      if (allocated(error)) return

      call check(error,.not. result%is_linear, &
                 "H2O should not be linear")
      if (allocated(error)) return

      ! ZPE should be positive and in reasonable range
      ! H2O ZPE ~ 13.0 kcal/mol (Gaussian reference)
      call check(error, result%zpe_kcalmol > 10.0_dp .and. result%zpe_kcalmol < 20.0_dp, &
                 "H2O ZPE should be in reasonable range (10-20 kcal/mol)")
      if (allocated(error)) return

      ! Thermal corrections should follow: Energy < Enthalpy < Gibbs (in magnitude)
      ! Actually: G = H - TS, so at 298K with positive S, G < H
      call check(error, result%thermal_correction_enthalpy > result%thermal_correction_energy, &
                 "Enthalpy correction should be > Energy correction (includes RT)")
      if (allocated(error)) return

      ! Gibbs should be less than Enthalpy (G = H - TS, S > 0)
      call check(error, result%thermal_correction_gibbs < result%thermal_correction_enthalpy, &
                 "Gibbs correction should be < Enthalpy correction")

   end subroutine test_full_thermochemistry_h2o

end module test_mqc_thermochemistry

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_thermochemistry, only: collect_mqc_thermochemistry_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_thermochemistry", collect_mqc_thermochemistry_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if

end program tester
