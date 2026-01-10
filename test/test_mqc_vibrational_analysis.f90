module test_mqc_vibrational_analysis
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_vibrational_analysis, only: compute_vibrational_frequencies, mass_weight_hessian, &
                                       project_translation_rotation, compute_vibrational_analysis, &
                                       compute_reduced_masses, compute_force_constants, &
                                       compute_cartesian_displacements, print_vibrational_analysis
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_vibrational_analysis_tests

   ! Tolerance for floating point comparisons
   real(dp), parameter :: tol = 1.0e-8_dp

contains

   subroutine collect_mqc_vibrational_analysis_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("mass_weighting_symmetry", test_mass_weighting_symmetry), &
                  new_unittest("mass_weighting_values", test_mass_weighting_values), &
                  new_unittest("diatomic_frequencies", test_diatomic_frequencies), &
                  new_unittest("frequency_ordering", test_frequency_ordering), &
                  new_unittest("imaginary_frequencies", test_imaginary_frequencies), &
                  new_unittest("trans_rot_projection", test_trans_rot_projection), &
                  new_unittest("reduced_masses", test_reduced_masses), &
                  new_unittest("force_constants", test_force_constants), &
                  new_unittest("cartesian_displacements", test_cartesian_displacements), &
                  new_unittest("full_vibrational_analysis", test_full_vibrational_analysis), &
                  new_unittest("print_vibrational_output", test_print_vibrational_output) &
                  ]
   end subroutine collect_mqc_vibrational_analysis_tests

   subroutine test_mass_weighting_symmetry(error)
      !! Test that mass-weighted Hessian preserves symmetry
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: hessian(6, 6), mw_hessian_ij, mw_hessian_ji
      real(dp), allocatable :: mw_hessian(:, :)
      integer :: element_numbers(2)
      integer :: i, j

      ! Create a simple symmetric 2-atom Hessian (H2)
      element_numbers = [1, 1]  ! Two hydrogen atoms

      ! Build a symmetric test Hessian
      hessian = 0.0_dp
      do i = 1, 6
         do j = i, 6
            hessian(i, j) = real(i + j, dp)*0.1_dp
            hessian(j, i) = hessian(i, j)
         end do
      end do

      call mass_weight_hessian(hessian, element_numbers, mw_hessian)

      ! Check symmetry of mass-weighted Hessian
      do i = 1, 6
         do j = i + 1, 6
            mw_hessian_ij = mw_hessian(i, j)
            mw_hessian_ji = mw_hessian(j, i)
            call check(error, abs(mw_hessian_ij - mw_hessian_ji) < tol, &
                       "Mass-weighted Hessian should be symmetric")
            if (allocated(error)) return
         end do
      end do

      deallocate (mw_hessian)

   end subroutine test_mass_weighting_symmetry

   subroutine test_mass_weighting_values(error)
      !! Test mass weighting calculation with known values
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: hessian(3, 3)
      real(dp), allocatable :: mw_hessian(:, :)
      integer :: element_numbers(1)
      real(dp) :: expected_value, mass_H

      ! Single hydrogen atom (3x3 Hessian)
      element_numbers = [1]  ! Hydrogen, mass ~1.008 amu
      mass_H = 1.008_dp

      ! Identity Hessian
      hessian = 0.0_dp
      hessian(1, 1) = 1.0_dp
      hessian(2, 2) = 1.0_dp
      hessian(3, 3) = 1.0_dp

      call mass_weight_hessian(hessian, element_numbers, mw_hessian)

      ! For identity matrix, H_mw(i,i) = 1.0 / mass_H
      expected_value = 1.0_dp/mass_H
      call check(error, abs(mw_hessian(1, 1) - expected_value) < tol, &
                 "Diagonal mass weighting should be 1/mass")
      if (allocated(error)) return

      call check(error, abs(mw_hessian(2, 2) - expected_value) < tol, &
                 "Diagonal mass weighting should be consistent")
      if (allocated(error)) return

      ! Off-diagonal should still be zero
      call check(error, abs(mw_hessian(1, 2)) < tol, &
                 "Off-diagonal should remain zero")

      deallocate (mw_hessian)

   end subroutine test_mass_weighting_values

   subroutine test_diatomic_frequencies(error)
      !! Test frequency calculation for a simple diatomic (H2-like)
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: hessian(6, 6)
      real(dp), allocatable :: frequencies(:)
      integer :: element_numbers(2)
      real(dp) :: force_constant
      integer :: n_positive

      ! Two hydrogen atoms
      element_numbers = [1, 1]

      ! Create a diatomic Hessian along z-axis
      ! For a diatomic with force constant k:
      ! H = [ k  0  0 -k  0  0]
      !     [ 0  k  0  0 -k  0]
      !     [ 0  0  k  0  0 -k]
      !     [-k  0  0  k  0  0]
      !     [ 0 -k  0  0  k  0]
      !     [ 0  0 -k  0  0  k]
      force_constant = 0.5_dp  ! Hartree/Bohr^2

      hessian = 0.0_dp
      ! Diagonal blocks
      hessian(1, 1) = force_constant
      hessian(2, 2) = force_constant
      hessian(3, 3) = force_constant
      hessian(4, 4) = force_constant
      hessian(5, 5) = force_constant
      hessian(6, 6) = force_constant
      ! Off-diagonal blocks
      hessian(1, 4) = -force_constant
      hessian(4, 1) = -force_constant
      hessian(2, 5) = -force_constant
      hessian(5, 2) = -force_constant
      hessian(3, 6) = -force_constant
      hessian(6, 3) = -force_constant

      call compute_vibrational_frequencies(hessian, element_numbers, frequencies)

      ! Should get 6 frequencies
      call check(error, size(frequencies) == 6, "Should have 6 frequencies for 2 atoms")
      if (allocated(error)) return

      ! Count positive (real) frequencies
      n_positive = count(frequencies > 0.0_dp)

      ! For a stable diatomic, we expect at least one positive frequency
      ! (the stretching mode). Translation and rotation modes should be near zero.
      call check(error, n_positive >= 1, "Should have at least one positive frequency")
      if (allocated(error)) return

      ! Check that frequencies are in reasonable range (not NaN or infinity)
      call check(error, all(abs(frequencies) < 1.0e10_dp), &
                 "Frequencies should be finite")

      deallocate (frequencies)

   end subroutine test_diatomic_frequencies

   subroutine test_frequency_ordering(error)
      !! Test that frequencies are returned in ascending order (eigenvalue order)
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: hessian(6, 6)
      real(dp), allocatable :: frequencies(:)
      integer :: element_numbers(2)
      integer :: i

      ! Two carbon atoms (heavier than H)
      element_numbers = [6, 6]

      ! Simple diagonal Hessian with distinct values
      hessian = 0.0_dp
      hessian(1, 1) = 0.1_dp
      hessian(2, 2) = 0.2_dp
      hessian(3, 3) = 0.3_dp
      hessian(4, 4) = 0.4_dp
      hessian(5, 5) = 0.5_dp
      hessian(6, 6) = 0.6_dp

      call compute_vibrational_frequencies(hessian, element_numbers, frequencies)

      ! LAPACK returns eigenvalues in ascending order, so frequencies should be ascending
      do i = 1, size(frequencies) - 1
         call check(error, frequencies(i) <= frequencies(i + 1), &
                    "Frequencies should be in ascending order")
         if (allocated(error)) return
      end do

      deallocate (frequencies)

   end subroutine test_frequency_ordering

   subroutine test_imaginary_frequencies(error)
      !! Test that negative eigenvalues produce negative (imaginary) frequencies
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: hessian(3, 3)
      real(dp), allocatable :: frequencies(:)
      integer :: element_numbers(1)

      ! Single atom with one negative eigenvalue (transition state-like)
      element_numbers = [1]  ! Hydrogen

      ! Hessian with one negative diagonal (saddle point)
      hessian = 0.0_dp
      hessian(1, 1) = -0.1_dp  ! Negative curvature
      hessian(2, 2) = 0.2_dp
      hessian(3, 3) = 0.3_dp

      call compute_vibrational_frequencies(hessian, element_numbers, frequencies)

      ! First eigenvalue is negative, so first frequency should be negative
      call check(error, frequencies(1) < 0.0_dp, &
                 "Negative eigenvalue should give negative (imaginary) frequency")
      if (allocated(error)) return

      ! Other frequencies should be positive
      call check(error, frequencies(2) > 0.0_dp, &
                 "Positive eigenvalue should give positive frequency")
      if (allocated(error)) return

      call check(error, frequencies(3) > 0.0_dp, &
                 "Positive eigenvalue should give positive frequency")

      deallocate (frequencies)

   end subroutine test_imaginary_frequencies

   subroutine test_trans_rot_projection(error)
      !! Test that translation/rotation projection produces 6 zero eigenvalues for a triatomic
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: hessian(9, 9)
      real(dp) :: coordinates(3, 3)
      real(dp), allocatable :: frequencies(:)
      real(dp), allocatable :: eigenvalues(:)
      integer :: element_numbers(3)
      integer :: n_near_zero
      integer :: i

      ! Water-like molecule: O at origin, two H atoms
      ! Use 3D geometry (not planar) to have all 6 trans/rot modes
      element_numbers = [8, 1, 1]  ! O, H, H

      ! Coordinates (Bohr) - 3D bent water geometry
      coordinates(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]         ! O
      coordinates(:, 2) = [1.8_dp, 0.0_dp, 0.5_dp]         ! H (out of xy plane)
      coordinates(:, 3) = [-0.6_dp, 1.7_dp, -0.3_dp]       ! H (out of xy plane)

      ! Create an identity Hessian - after projection, 6 eigenvalues should be 0
      ! and 3 should remain as 1.0
      hessian = 0.0_dp
      do i = 1, 9
         hessian(i, i) = 1.0_dp
      end do

      ! Compute frequencies WITH projection
      call compute_vibrational_frequencies(hessian, element_numbers, frequencies, &
                                           eigenvalues_out=eigenvalues, &
                                           coordinates=coordinates, &
                                           project_trans_rot=.true.)

      ! Should get 9 frequencies
      call check(error, size(frequencies) == 9, "Should have 9 frequencies for 3 atoms")
      if (allocated(error)) return

      ! Count near-zero eigenvalues (the 6 projected trans/rot modes)
      ! Using a threshold of 1e-6 for "near zero"
      n_near_zero = count(abs(eigenvalues) < 1.0e-6_dp)

      ! For a non-linear molecule, we expect exactly 6 near-zero eigenvalues
      ! (3 translations + 3 rotations)
      call check(error, n_near_zero == 6, &
                 "Should have 6 near-zero eigenvalues after projection")
      if (allocated(error)) return

      ! The remaining 3 should be positive (vibrational modes)
      call check(error, count(eigenvalues > 1.0e-6_dp) == 3, &
                 "Should have 3 positive eigenvalues (vibrational modes)")

      deallocate (frequencies, eigenvalues)

   end subroutine test_trans_rot_projection

   subroutine test_reduced_masses(error)
      !! Test reduced mass calculation
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: hessian(6, 6)
      real(dp), allocatable :: frequencies(:), eigenvalues(:), eigenvectors(:, :)
      real(dp), allocatable :: reduced_masses(:)
      integer :: element_numbers(2)
      real(dp) :: force_constant

      ! Two hydrogen atoms - diatomic
      element_numbers = [1, 1]
      force_constant = 0.5_dp

      ! Diatomic Hessian
      hessian = 0.0_dp
      hessian(1, 1) = force_constant
      hessian(2, 2) = force_constant
      hessian(3, 3) = force_constant
      hessian(4, 4) = force_constant
      hessian(5, 5) = force_constant
      hessian(6, 6) = force_constant
      hessian(1, 4) = -force_constant
      hessian(4, 1) = -force_constant
      hessian(2, 5) = -force_constant
      hessian(5, 2) = -force_constant
      hessian(3, 6) = -force_constant
      hessian(6, 3) = -force_constant

      call compute_vibrational_frequencies(hessian, element_numbers, frequencies, &
                                           eigenvalues_out=eigenvalues, &
                                           eigenvectors=eigenvectors)

      call compute_reduced_masses(eigenvectors, element_numbers, reduced_masses)

      ! Should get 6 reduced masses
      call check(error, size(reduced_masses) == 6, "Should have 6 reduced masses")
      if (allocated(error)) return

      ! All reduced masses should be positive
      call check(error, all(reduced_masses > 0.0_dp), "Reduced masses should be positive")
      if (allocated(error)) return

      ! For H2, the vibrational mode reduced mass should be approximately μ = m1*m2/(m1+m2) = 0.504 amu
      ! The stretching mode (highest frequency) should have reduced mass near this value
      call check(error, reduced_masses(6) < 2.0_dp, &
                 "Stretching mode reduced mass should be less than single atom mass")

      deallocate (frequencies, eigenvalues, eigenvectors, reduced_masses)

   end subroutine test_reduced_masses

   subroutine test_force_constants(error)
      !! Test force constant calculation
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: hessian(6, 6)
      real(dp), allocatable :: frequencies(:), eigenvalues(:), eigenvectors(:, :)
      real(dp), allocatable :: reduced_masses(:), force_constants(:), force_constants_mdyne(:)
      integer :: element_numbers(2)
      real(dp) :: input_k

      ! Two hydrogen atoms
      element_numbers = [1, 1]
      input_k = 0.5_dp  ! Hartree/Bohr^2

      ! Diatomic Hessian
      hessian = 0.0_dp
      hessian(1, 1) = input_k
      hessian(2, 2) = input_k
      hessian(3, 3) = input_k
      hessian(4, 4) = input_k
      hessian(5, 5) = input_k
      hessian(6, 6) = input_k
      hessian(1, 4) = -input_k
      hessian(4, 1) = -input_k
      hessian(2, 5) = -input_k
      hessian(5, 2) = -input_k
      hessian(3, 6) = -input_k
      hessian(6, 3) = -input_k

      call compute_vibrational_frequencies(hessian, element_numbers, frequencies, &
                                           eigenvalues_out=eigenvalues, &
                                           eigenvectors=eigenvectors)

      call compute_reduced_masses(eigenvectors, element_numbers, reduced_masses)
      call compute_force_constants(eigenvalues, reduced_masses, force_constants, &
                                   force_constants_mdyne)

      ! Should get 6 force constants
      call check(error, size(force_constants) == 6, "Should have 6 force constants")
      if (allocated(error)) return

      call check(error, size(force_constants_mdyne) == 6, "Should have 6 force constants in mdyne/A")
      if (allocated(error)) return

      ! Force constants should be related to eigenvalues by k = eigenvalue * reduced_mass
      ! The stretching modes should have positive force constants
      call check(error, force_constants(6) > 0.0_dp, &
                 "Stretching mode force constant should be positive")
      if (allocated(error)) return

      ! mdyne/A should be ~ 15.57 times atomic units
      call check(error, abs(force_constants_mdyne(6) - force_constants(6)*15.569141_dp) < 1.0e-6_dp, &
                 "mdyne/A conversion should be correct")

      deallocate (frequencies, eigenvalues, eigenvectors, reduced_masses, force_constants, &
                  force_constants_mdyne)

   end subroutine test_force_constants

   subroutine test_cartesian_displacements(error)
      !! Test Cartesian displacement calculation
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: hessian(6, 6)
      real(dp), allocatable :: frequencies(:), eigenvectors(:, :)
      real(dp), allocatable :: cart_disp(:, :)
      integer :: element_numbers(2)
      real(dp) :: force_constant, max_disp
      integer :: k

      ! Two hydrogen atoms
      element_numbers = [1, 1]
      force_constant = 0.5_dp

      ! Diatomic Hessian
      hessian = 0.0_dp
      hessian(1, 1) = force_constant
      hessian(2, 2) = force_constant
      hessian(3, 3) = force_constant
      hessian(4, 4) = force_constant
      hessian(5, 5) = force_constant
      hessian(6, 6) = force_constant
      hessian(1, 4) = -force_constant
      hessian(4, 1) = -force_constant
      hessian(2, 5) = -force_constant
      hessian(5, 2) = -force_constant
      hessian(3, 6) = -force_constant
      hessian(6, 3) = -force_constant

      call compute_vibrational_frequencies(hessian, element_numbers, frequencies, &
                                           eigenvectors=eigenvectors)

      call compute_cartesian_displacements(eigenvectors, element_numbers, cart_disp)

      ! Should have 6x6 matrix
      call check(error, size(cart_disp, 1) == 6, "Should have 6 rows in Cartesian displacements")
      if (allocated(error)) return
      call check(error, size(cart_disp, 2) == 6, "Should have 6 columns in Cartesian displacements")
      if (allocated(error)) return

      ! Each mode should be normalized so max displacement = 1 (Gaussian convention)
      do k = 1, 6
         max_disp = maxval(abs(cart_disp(:, k)))
         call check(error, abs(max_disp - 1.0_dp) < 1.0e-10_dp .or. max_disp < 1.0e-10_dp, &
                    "Cartesian displacements should be normalized to max=1")
         if (allocated(error)) return
      end do

      deallocate (frequencies, eigenvectors, cart_disp)

   end subroutine test_cartesian_displacements

   subroutine test_full_vibrational_analysis(error)
      !! Test the complete vibrational analysis wrapper
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: hessian(9, 9)
      real(dp) :: coordinates(3, 3)
      real(dp), allocatable :: frequencies(:), reduced_masses(:), force_constants(:)
      real(dp), allocatable :: cart_disp(:, :), force_constants_mdyne(:)
      integer :: element_numbers(3)
      integer :: i

      ! Water-like molecule
      element_numbers = [8, 1, 1]

      ! Coordinates (Bohr)
      coordinates(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      coordinates(:, 2) = [1.8_dp, 0.0_dp, 0.5_dp]
      coordinates(:, 3) = [-0.6_dp, 1.7_dp, -0.3_dp]

      ! Identity Hessian for testing
      hessian = 0.0_dp
      do i = 1, 9
         hessian(i, i) = 1.0_dp
      end do

      call compute_vibrational_analysis(hessian, element_numbers, frequencies, &
                                        reduced_masses, force_constants, cart_disp, &
                                        coordinates=coordinates, &
                                        project_trans_rot=.true., &
                                        force_constants_mdyne=force_constants_mdyne)

      ! Check all outputs have correct size
      call check(error, size(frequencies) == 9, "Should have 9 frequencies")
      if (allocated(error)) return

      call check(error, size(reduced_masses) == 9, "Should have 9 reduced masses")
      if (allocated(error)) return

      call check(error, size(force_constants) == 9, "Should have 9 force constants")
      if (allocated(error)) return

      call check(error, size(force_constants_mdyne) == 9, "Should have 9 force constants in mdyne/A")
      if (allocated(error)) return

      call check(error, size(cart_disp, 1) == 9, "Should have 9 rows in Cartesian displacements")
      if (allocated(error)) return

      call check(error, size(cart_disp, 2) == 9, "Should have 9 columns in Cartesian displacements")

      deallocate (frequencies, reduced_masses, force_constants, cart_disp, force_constants_mdyne)

   end subroutine test_full_vibrational_analysis

   subroutine test_print_vibrational_output(error)
      !! Test that print_vibrational_analysis runs without errors
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: hessian(6, 6)
      real(dp), allocatable :: frequencies(:), reduced_masses(:), force_constants(:)
      real(dp), allocatable :: cart_disp(:, :), force_constants_mdyne(:)
      integer :: element_numbers(2)
      real(dp) :: force_constant

      ! Two hydrogen atoms
      element_numbers = [1, 1]
      force_constant = 0.5_dp

      ! Diatomic Hessian
      hessian = 0.0_dp
      hessian(1, 1) = force_constant
      hessian(2, 2) = force_constant
      hessian(3, 3) = force_constant
      hessian(4, 4) = force_constant
      hessian(5, 5) = force_constant
      hessian(6, 6) = force_constant
      hessian(1, 4) = -force_constant
      hessian(4, 1) = -force_constant
      hessian(2, 5) = -force_constant
      hessian(5, 2) = -force_constant
      hessian(3, 6) = -force_constant
      hessian(6, 3) = -force_constant

      ! Compute full vibrational analysis
      call compute_vibrational_analysis(hessian, element_numbers, frequencies, &
                                        reduced_masses, force_constants, cart_disp, &
                                        force_constants_mdyne=force_constants_mdyne)

      ! Print should run without errors
      call print_vibrational_analysis(frequencies, reduced_masses, force_constants, &
                                      cart_disp, element_numbers, &
                                      force_constants_mdyne=force_constants_mdyne)

      ! If we got here without crashing, test passed
      call check(error, .true., "Print vibrational analysis should complete without errors")

      deallocate (frequencies, reduced_masses, force_constants, cart_disp, force_constants_mdyne)

   end subroutine test_print_vibrational_output

end module test_mqc_vibrational_analysis

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_vibrational_analysis, only: collect_mqc_vibrational_analysis_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_vibrational_analysis", collect_mqc_vibrational_analysis_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program tester
